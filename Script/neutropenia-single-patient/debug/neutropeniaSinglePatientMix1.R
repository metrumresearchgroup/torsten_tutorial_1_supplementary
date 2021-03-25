rm(list = ls())
gc()

modelName <- "neutropeniaSinglePatientMix1"
simModelName <- "neutropeniaSinglePatient1Sim"

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- scriptDir
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- projectDir
outDir <- file.path(modelDir, modelName)
## toolsDir <- file.path(scriptDir, "tools")

.libPaths("~/.R/R_libs")

library("cmdstanr")
library("bayesplot")
library("tidyverse")
library(bayesplot)
## Go back to default ggplot2 theme that was overridden by bayesplot
theme_set(theme_gray())
library(tidyverse)
library(parallel)
## source(file.path(toolsDir, "stanTools.R"))
## source(file.path(toolsDir, "functions.R"))

options(mc.cores = parallel::detectCores())

set.seed(10271998) ## not required but assures repeatable results

################################################################################################
### Simulate ME-2 plasma concentrations and ANC values

## Parameter values

ka = 2.0
CL = 10 # L/h
V1 = 35 # L
V2 = 105 # L
Q = 15
sigma = 0.1

alpha <- 3E-4

## PD parameters based on J Clin Oncol 20:4713-4721 (2002)

## drug-independent parameters
mtt <- 125 # exp(4.43)
circ0 <- 5 # exp(2.33)
gamma <- 0.17 # 0.12

sigmaNeut <- 0.1

## Observation and dosing times
doseTimes <- seq(0, 168, by = 12)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), c(xpk, 12, 18, 24) + 168)
xneut <- seq(0, 672, by = 24)
time <- sort(unique(c(xpk, xneut, doseTimes)))

## Assemble data set for simulation using Stan
obsData <- data.frame(time = time) %>% mutate(amt = 0, cmt = 1, evid = 0)

doseData <- data.frame(time = doseTimes) %>%
    mutate(amt = 80 * 1000, # mcg
           cmt = 1,
           evid = 1)

allData <- doseData %>%
    bind_rows(obsData) %>%
    arrange(time, desc(evid))

nt <- nrow(allData)

dataSim <- with(allData,
                list(nt = nt,
                     amt = amt,
                     cmt = cmt,
                     evid = evid,
                     time = time,
                     CL = CL,
                     Q = Q,
                     V1 = V1,
                     V2 = V2,
                     ka = ka,
                     circ0 = circ0,
                     mtt = mtt,
                     alpha = alpha,
                     gamma = gamma,
                     sigma = sigma,
                     sigmaNeut = sigmaNeut))

### Simulate using Stan
mod <- cmdstan_model(file.path(modelDir, paste(simModelName, ".stan", sep = "")), quiet=FALSE)
sim <- mod$sample(data = dataSim, fixed_param=TRUE, iter_sampling=1, chains=1)

################################################################################################
### Assemble data set for fitting via Stan

xdata <- allData %>%
    bind_cols(as.data.frame(sim$draws(variables = "cObs")) %>%
              gather(factor_key = TRUE) %>%
              select(cObs = value)) %>%
    bind_cols(as.data.frame(sim$draws(variables = "neutObs")) %>%
              gather(factor_key = TRUE) %>%
              select(neutObs = value))

xdata <- xdata %>%
    mutate(cObs = ifelse(time %in% xpk & time != 0 & evid == 0, cObs, NA),
           neutObs = ifelse(time %in% xneut & evid == 0, neutObs, NA))

head(xdata)

dir.create(figDir, recursive=TRUE)
dir.create(tabDir, recursive=TRUE)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

p1 <- ggplot(xdata %>% filter(!is.na(cObs)), aes(x = time, y = cObs))
p1 + geom_point() + geom_line() +
    labs(x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8))

p1 <- ggplot(xdata %>% filter(!is.na(neutObs)), aes(x = time, y = neutObs))
p1 + geom_point() + geom_line() +
    labs(x = "time (h)",
         y = "ANC") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8))

## Indices of records containing observed concentrations
iObsPK <- with(xdata, (1:nrow(xdata))[!is.na(cObs) & evid == 0])
nObsPK <- length(iObsPK)
## Indices of records containing observed neutrophil counts
iObsPD <- with(xdata, (1:nrow(xdata))[!is.na(neutObs) & evid == 0])
nObsPD <- length(iObsPD)

## Parameters for informative priors

CLPrior = 10
QPrior = 15
V1Prior = 35
V2Prior = 105
kaPrior = 2
CLPriorCV = 0.5 ## 0.10
QPriorCV = 0.5 ## 0.18
V1PriorCV = 0.5 ## 0.14
V2PriorCV = 0.5 ## 0.17
kaPriorCV = 0.5 ## 0.16

circ0Prior <- 5
circ0PriorCV <- 0.20
mttPrior <- 125
mttPriorCV <- 0.2
gammaPrior <- 0.17
gammaPriorCV <- 0.2
alphaPrior <- 3.0E-4
alphaPriorCV <- 1 ## 0.2

## create data set
data <- with(xdata,
             list(nt = nt,
                  nObsPK = nObsPK,
                  iObsPK = iObsPK,
                  nObsPD = nObsPD,
                  iObsPD = iObsPD,
                  amt = amt,
                  cmt = cmt,
                  evid = evid,
                  time = time,
                  cObs = cObs[iObsPK],
                  neutObs = neutObs[iObsPD],
                  CLPrior = CLPrior,
                  QPrior = QPrior,
                  V1Prior = V1Prior,
                  V2Prior = V2Prior,
                  kaPrior = kaPrior,
                  CLPriorCV = CLPriorCV,
                  QPriorCV = QPriorCV,
                  V1PriorCV = V1PriorCV,
                  V2PriorCV = V2PriorCV,
                  kaPriorCV = kaPriorCV,
                  circ0Prior = circ0Prior,
                  circ0PriorCV = circ0PriorCV,
                  mttPrior = mttPrior,
                  mttPriorCV = mttPriorCV,
                  gammaPrior = gammaPrior,
                  gammaPriorCV = gammaPriorCV,
                  alphaPrior = alphaPrior,
                  alphaPriorCV = alphaPriorCV
                  ))

### create initial estimates
init <- function(){
  CL = exp(rnorm(1, log(CLPrior), CLPriorCV))
  Q = exp(rnorm(1, log(QPrior), QPriorCV))
  V1 = exp(rnorm(1, log(V1Prior), V1PriorCV))
  V2 = exp(rnorm(1, log(V2Prior), V2PriorCV))
  lambda1 <- (CL / V1 + Q / V1 + Q / V2 +
                sqrt((CL / V1 + Q / V1 + Q / V2)^2 -
                       4 * CL / V1 * Q / V2)) / 2
  list(CL = CL,
       Q = Q,
       V1 = V1,
       V2 = V2,
       ka = lambda1 + exp(rnorm(1, log(kaPrior), kaPriorCV)),
       sigma = 0.2,
       alpha = exp(rnorm(1, log(alphaPrior), alphaPriorCV)),
       mtt = exp(rnorm(1, log(mttPrior), mttPriorCV)),
       circ0 = exp(rnorm(1, log(circ0Prior), circ0PriorCV)),
       gamma = exp(rnorm(1, log(gammaPrior), gammaPriorCV)),
       sigmaNeut = 0.2)
}

### Specify the variables for which you want history and density plots

parametersToPlot <- c("CL", "Q", "V1", "V2", "ka",
                      "sigma", "alpha", "mtt", "circ0",
                      "gamma", "sigmaNeut")

## Additional variables to monitor
otherRVs <- c("cObsPred", "neutObsPred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

mod.fit <- cmdstan_model(file.path(modelDir, paste(modelName, ".stan", sep = "")), quiet=FALSE)
fit <- mod.fit$sample(data = data,
                      iter_sampling = 500,
                      iter_warmup = 1000,
                      thin = 1,
                      init = init,
                      chains = 4,
                      cores = 4,
                      refresh = 10,
                      adapt_delta = 0.90, step_size = 0.01)

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

options(bayesplot.base_size = 12,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "brightblue")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 5, myTheme = myTheme)
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, byChain = TRUE,
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16,
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))

pairs(fit, pars = parametersToPlot[!grepl("rho", parametersToPlot)])

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
### posterior predictive distributions

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata) %>%
  filter(time <= max(xpk))

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(title = "individual predictions",
       x = "time (h)",
       y = "ME-2 plasma concentration (ng/mL)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))

print(p1 + geom_line(aes(x = time, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))

pred <- as.data.frame(fit, pars = "neutObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

p1 <- ggplot(pred, aes(x = time, y = neutObs))
p1 <- p1 + geom_point() +
  labs(title = "individual predictions",
       x = "time (h)",
       y = "ANC") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8))

print(p1 + geom_line(aes(x = time, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))

dev.off()
