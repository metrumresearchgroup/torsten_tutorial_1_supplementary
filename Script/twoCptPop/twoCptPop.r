
rm(list = ls())
gc()
set.seed(1954)

# adjust to your settings
.libPaths("~/Rlib/")
setwd("~/Code/torsten_tutorial_psp/Script/twoCptPop")

library(cmdstanr)
library(rjson)
library(bayesplot)
library(posterior)

set_cmdstan_path("../../Torsten/cmdstan/")
bayesplot::color_scheme_set("mix-blue-green")

##########################################################################
## Read in data and create inits.

data <- fromJSON(file = "twoCptPop.data.json")

# Draw initial conditions from the prior
# draw parameters from prior
init <- function () {
  n_subjects <- data$nSubjects
  pop_var <- c(0.25, 0.5, 0.25, 0.5, 1)
  
  CL_pop <- exp(rnorm(1, log(10), pop_var[1]))
  Q_pop <- exp(rnorm(1, log(15), pop_var[2]))
  VC_pop <- exp(rnorm(1, log(35), pop_var[3]))
  VP_pop <- exp(rnorm(1, log(105), pop_var[4]))
  ka_pop <- exp(rnorm(1, log(2.5), pop_var[5]))
  omega <- abs(rnorm(5, 0, pop_var))
  
  theta_pop <- c(CL_pop, Q_pop, VC_pop, VP_pop, ka_pop)
  theta <- matrix(NA, n_subjects, length(theta_pop))
  for (j in 1:n_subjects) {
    theta[j, ] <- exp(rnorm(length(theta_pop), log(theta_pop), omega))
  }
  
  list(CL_pop = CL_pop, Q_pop = Q_pop, VC_pop = VC_pop, VP_pop = VP_pop,
       ka_pop = ka_pop, omega = omega, theta = theta,
       sigma = abs(rnorm(1, 0, 1)))   
}

##########################################################################
## Compile and fit the model
model_name <- "twoCptPop"
mod <- cmdstan_model(paste0(model_name, ".stan"))

n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = n_chains,
                  iter_warmup = 500, iter_sampling = 500,
                  seed = 123, adapt_delta = 0.95)


dir.create("deliv")
fit$save_object(paste0("deliv/", model_name, ".fit.RDS"))

# The saved fit object can be read using the below line
# fit <- readRDS(paste0("saved_fit/", model_name, ".fit.RDS"))

##########################################################################
## Check the inference

print(fit$time(), digits = 3)

pars = c("lp__", "CL_pop", "Q_pop", "VC_pop", "VP_pop",
         "ka_pop", "sigma")
fit$summary(pars)
bayesplot::mcmc_trace(fit$draws(), pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

##########################################################################
## Posterior predictive checks

yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))))
yrep <- yrep[, -((ncol(yrep) - 2):ncol(yrep))]

yobs <- data$cObs
time <- data$time[data$iObs]
patientID <- with(data, rep(1:nSubjects, each = nObs / nSubjects))

# within patient predictions
bayesplot::ppc_intervals_grouped(y = yobs, yrep = yrep, x = time, patientID)
bayesplot::ppc_ribbon_grouped(y = yobs, yrep = yrep, x = time, patientID)

# predictions for new patient
yrepNew <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("cObsNewPred"))))
yrepNew <- yrepNew[, -((ncol(yrepNew) - 2):ncol(yrepNew))]

bayesplot::ppc_intervals_grouped(y = yobs, yrep = yrepNew, x = time, patientID)
bayesplot::ppc_ribbon_grouped(y = yobs, yrep = yrepNew, x = time, patientID)
