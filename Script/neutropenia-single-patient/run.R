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

## additional packages
.libPaths("~/.R/R_libs")

library("cmdstanr")
library("bayesplot")
library("posterior")
library("tidyverse")
library(bayesplot)
## Go back to default ggplot2 theme that was overridden by bayesplot
theme_set(theme_gray())
library(tidyverse)
library(parallel)

options(mc.cores = parallel::detectCores())

bayesplot::color_scheme_set("mix-blue-green")

set.seed(10271998) ## not required but assures repeatable results

set_cmdstan_path("../../Torsten/cmdstan")

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
obsData <- data.frame(time = time) %>%
    mutate(amt = 0,
           cmt = 1,
           evid = 0)

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

mod.sim <- cmdstan_model(file.path(modelDir, paste(simModelName, ".stan", sep = "")), quiet=FALSE)
sim <- mod.sim$sample(data = dataSim, fixed_param=TRUE, iter_sampling=1, chains=1)

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

dir.create(figDir)
dir.create(tabDir)

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
                  alphaPriorCV = alphaPriorCV,
                  rtol = 1.e-6,
                  atol = 1.e-6,
                  max_num_step = 10000
                  ))

### Specify the variables for which you want history and density plots

parametersToPlot <- c("CL", "Q", "V1", "V2", "ka",
                      "sigma", "alpha", "mtt", "circ0",
                      "gamma", "sigmaNeut")

## Additional variables to monitor
otherRVs <- c("cObsPred", "neutObsPred")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

mod.fit <- cmdstan_model(file.path(modelDir, paste(modelName, ".stan", sep = "")), quiet=FALSE)
fit <- mod.fit$sample(data = data,
                      iter_sampling = 500,
                      iter_warmup = 1000,
                      thin = 1,
                      init = init,
                      chains = 4,
                      parallel_chains = 4,
                      refresh = 10,
                      adapt_delta = 0.80, step_size = 0.1)

dir.create(outDir)
fit$save_output_files()

################################################################################################
## posterior predictive checks

c.rep <- fit$draws(variables = c("cObsPred")) %>% as_draws_df() %>% as.matrix()
c.rep <- c.rep[ , data$iObsPK]

neut.rep <- fit$draws(variables = c("neutObsPred")) %>% as_draws_df() %>% as.matrix()
neut.rep <- neut.rep[ , data$iObsPD]

c.obs <- data$cObs
neut.obs <- data$neutObs
time <- data$time[-1]

## bayesplot::ppc_intervals(y = c.obs, yrep = c.rep, x = data$time[data$iObsPK])
ppc.pk <- bayesplot::ppc_ribbon(y = c.obs, yrep = c.rep, x = data$time[data$iObsPK]) +
    scale_x_continuous(name="time (h)") +
    scale_y_continuous(name="drug plasma concentration (ng/mL)") + theme(axis.title=element_text(size=20),axis.text=element_text(size=20),legend.text=element_text(size=20))
ggsave(file.path(figDir, "neutrophil_ppc_pk.pdf"))

ppc.pd <- bayesplot::ppc_ribbon(y = neut.obs, yrep = neut.rep, x = data$time[data$iObsPD]) +
    scale_x_continuous(name="time (h)") +
    scale_y_continuous(name="Neutrophil counts") + theme(axis.title=element_text(size=20),axis.text=element_text(size=20),legend.text=element_text(size=20))
ggsave(file.path(figDir, "neutrophil_ppc_pd.pdf"))
