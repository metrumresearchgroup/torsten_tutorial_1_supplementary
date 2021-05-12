
rm(list = ls())
gc()
set.seed(1954)

# adjust to your settings
##.libPaths("~/Rlib/")
##setwd("~/Code/torsten_tutorial_psp/Script/twoCptPop")
scriptDir <- getwd() ## assumes you navigated to twoCptPop directory

library(tidyverse)
library(cmdstanr)
library(rjson)
library(bayesplot)
library(posterior)
library(vpc)

set_cmdstan_path(file.path(dirname(dirname(scriptDir)), "Torsten", "cmdstan"))
bayesplot::color_scheme_set("mix-blue-green")

qnorm.trunc = function(p,mean=0,sd=1,lower=-Inf,upper=Inf)
  qnorm(p*pnorm(upper,mean,sd)+(1-p)*pnorm(lower,mean,sd),mean,sd)

rnorm.trunc = function(n,mean=0,sd=1,lower=-Inf,upper=Inf)
  qnorm.trunc(runif(n),mean,sd,lower,upper)

##########################################################################
## Read in data and create inits.

## NONMEM-style data set
data1 <- read_csv("twoCptPop.data.csv")

## Reformat for Stan
nEvent <- nrow(data1)
start <- (1:nEvent)[!duplicated(data1$id)]
end <- c(start[-1] - 1, nEvent)
n_subjects <- length(unique(data1$id))

## Indices of records containing observed  data
iObs <- with(data1, (1:nEvent)[!is.na(cObs) & evid == 0])
nObs <- length(iObs)

data <- c(as.list(data1 %>% select(-cObs, -weight)),
          list(start = start,
               end = end,
               iObs = iObs,
               nEvent = nEvent,
               nSubjects = n_subjects,
               nIIV = 5,
               nObs = nObs,
               weight = with(data1, weight[!duplicated(id)]),
               cObs = data1$cObs[iObs]))

# Draw initial conditions from the prior
# draw parameters from prior
init <- function () {
  n_subjects <- data$nSubjects
  pop_var <- c(0.25, 0.5, 0.25, 0.5, 0.25)
  
  CL_pop <- exp(rnorm(1, log(10), pop_var[1]))
  Q_pop <- exp(rnorm(1, log(15), pop_var[2]))
  VC_pop <- exp(rnorm(1, log(35), pop_var[3]))
  VP_pop <- exp(rnorm(1, log(105), pop_var[4]))
  ka_pop <- exp(rnorm.trunc(1, log(2.5), pop_var[5],
                            lower = log((CL_pop / VC_pop + Q_pop / VC_pop + Q_pop / VP_pop +
                                       sqrt((CL_pop / VC_pop + Q_pop / VC_pop + Q_pop / VP_pop)^2 -
                                              4 * CL_pop / VC_pop * Q_pop / VP_pop)) / 2)))
  omega <- abs(rnorm(5, 0, pop_var))
  
  theta_pop <- c(CL_pop, Q_pop, VC_pop, VP_pop, ka_pop)
  theta <- matrix(NA, n_subjects, length(theta_pop))
  for (j in 1:n_subjects) {
    theta[j, ] <- exp(rnorm(length(theta_pop), log(theta_pop), omega))
  }
##  eta <- matrix(NA, n_subjects, length(theta_pop))
##  for (j in 1:n_subjects) {
##    eta[j, ] <- exp(rnorm(length(theta_pop), 0, 1))
##  }
  
  list(CL_pop = CL_pop, Q_pop = Q_pop, VC_pop = VC_pop, VP_pop = VP_pop,
       ka_pop = ka_pop, omega = omega, theta = theta,
       sigma = abs(rnorm(1, 0, 1)))   
}

##########################################################################
## Compile and fit the model
model_name <- "twoCptPop"
mod <- cmdstan_model(paste0(model_name, ".stan"))

dir.create("deliv")
n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = n_chains,
                  iter_warmup = 250, iter_sampling = 250,
                  seed = 123, adapt_delta = 0.8,
                  refresh = 10,
                  output_dir = "deliv")

fit$save_object(file.path("deliv", paste0(model_name, ".fit.RDS")))

# The saved fit object can be read using the below line
# fit <- readRDS(file.path("deliv", paste0(model_name, ".fit.RDS")))

##########################################################################
## Check the inference

print(fit$time(), digits = 3)

pars = c("lp__", "CL_pop", "Q_pop", "VC_pop", "VP_pop",
         "ka_pop", "sigma")
fit$summary(c(pars, "omega"))
bayesplot::mcmc_trace(fit$draws(pars))
bayesplot::mcmc_trace(fit$draws("omega"))
bayesplot::mcmc_dens_overlay(fit$draws(pars))
bayesplot::mcmc_dens_overlay(fit$draws("omega"))

##########################################################################
## Posterior predictive checks

yrep <- as_draws_matrix(fit$draws(variables = c("concentrationObsPred")))

yobs <- data1$cObs[iObs]
time <- data1$time[iObs]
patientID <- data1$id[iObs]

# within patient predictions
bayesplot::ppc_intervals_grouped(y = yobs, yrep = yrep, x = time, patientID)

# predictions for new patient
yrepNew <- as_draws_matrix(fit$draws(variables = c("cObsNewPred")))

bayesplot::ppc_intervals_grouped(y = yobs, yrep = yrepNew, x = time, patientID)

## Combine both types of predictions
predInd <-  posterior::as_draws_df(fit$draws("concentrationObsPred")) %>%
  pivot_longer(cols = starts_with("concentrationObsPred"),
               names_transform = list(name = forcats::fct_inorder)) %>%
  group_by(name) %>%
  summarize(lbInd = quantile(value, probs = 0.05, na.rm = TRUE),
            medianInd = quantile(value, probs = 0.5, na.rm = TRUE),
            ubInd = quantile(value, probs = 0.95, na.rm = TRUE))

predPop <- posterior::as_draws_df(fit$draws("cObsNewPred")) %>%
  pivot_longer(cols = starts_with("cObsNewPred"),
               names_transform = list(name = forcats::fct_inorder)) %>%
  group_by(name) %>%
  summarize(lbPop = quantile(value, probs = 0.05, na.rm = TRUE),
            medianPop = quantile(value, probs = 0.5, na.rm = TRUE),
            ubPop = quantile(value, probs = 0.95, na.rm = TRUE))

predAll <- bind_cols(data1[iObs,], predInd, predPop)

ppc1 <- ggplot(predAll, 
                   aes(x = time, y = cObs)) +
  geom_line(aes(x = time, y = medianPop, 
                color = "population")) +
  geom_ribbon(aes(ymin = lbPop, ymax = ubPop, 
                  fill = "population"), 
              alpha = 0.25) +
  geom_line(aes(x = time, y = medianInd, 
                color = "individual")) +
  geom_ribbon(aes(ymin = lbInd, ymax = ubInd, 
                  fill = "individual"), 
              alpha = 0.25) +
  scale_color_brewer(name  ="prediction",
                     breaks=c("individual", "population"),
                     palette = "Set1") +
  scale_fill_brewer(name  ="prediction",
                    breaks=c("individual", "population"),
                    palette = "Set1")
ppc1 + geom_point(na.rm = TRUE) +
  scale_y_log10() +
  labs(x = "time (h)",
       y = "plasma drug concentration") +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 8)) +
  facet_wrap(~ id)

nPost <- fit$metadata()$iter_sampling / fit$metadata()$thin

predPop <- posterior::as_draws_df(fit$draws("cObsNewPred")) %>%
  pivot_longer(cols = starts_with("cObsNewPred"),
               names_transform = list(name = forcats::fct_inorder)) %>%
  mutate(id = rep(data1$id[iObs], n_chains * nPost),
         time = rep(data1$time[iObs], n_chains * nPost),
         cObs = rep(data1$cObs[iObs], n_chains * nPost))

vpcPlot <- vpc(sim = predPop %>% filter(!is.na(cObs)),
               obs = data1 %>% filter(!is.na(cObs)),
               obs_cols = list(id = "id", dv = "cObs", idv = "time"),
               sim_cols = list(is = "id", dv = "value", idv = "time"),
               ci = c(0.1, 0.9), pi = c(0.1, 0.9),
               show = list(obs_dv = FALSE, obs_ci = TRUE, 
                           pi = FALSE, pi_as_area = FALSE, pi_ci = TRUE, 
                           obs_median = TRUE, sim_median = FALSE, sim_median_ci = TRUE),
               vpc_theme = new_vpc_theme(update = list(
                 sim_pi_fill = "red",
                 sim_pi_alpha = 0.25,
                 obs_ci_size = 1,
                 obs_ci_color = "red",
                 obs_median_color = "blue",
                 obs_ci_linetype = 1,
                 sim_median_fill = "blue",
                 sim_median_alpha = 0.25)))

vpcPlot + geom_point(data = data1, aes(x = time, y = cObs), alpha = 0.1) +
  scale_y_log10() +
  labs(x = "time (h)",
       y = "plasma drug concentration") +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12))

