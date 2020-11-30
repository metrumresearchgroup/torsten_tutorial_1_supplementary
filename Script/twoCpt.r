
rm(list = ls())
gc()
set.seed(1954)

# adjust to your settings
.libPaths("~/Rlib/")
setwd("~/Code/torsten_tutorial_psp/Script")

library(cmdstanr)
library(rjson)
library(bayesplot)
library(posterior)

set_cmdstan_path("../Torsten/cmdstan/")
# bayesplot::color_scheme_set("viridisC")
bayesplot::color_scheme_set("mix-blue-green")

model_name <- "twoCpt"

##########################################################################
## Read in data and create inits.

data <- fromJSON(file = "data/twoCpt.data.json")

# Draw initial conditions from the prior
init <- function() {
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(15), 0.2)),
       VC = exp(rnorm(1, log(35), 0.2)),
       VP = exp(rnorm(1, log(105), 0.2)),
       ka = exp(rnorm(1, log(2), 0.2)),
       sigma = abs(rcauchy(1, 0, 1)))
}

##########################################################################
## Compile and fit the model

mod <- cmdstan_model(paste0("model/", model_name, ".stan"))

n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains,
                  parallel_chains = n_chains,
                  iter_warmup = 500, iter_sampling = 500)


fit$save_object(paste0("saved_fit/", model_name, ".fit.RDS"))

# The saved fit object can be read using the below line
# fit <- readRDS(paste0("saved_fit/", model_name, ".fit.RDS"))

##########################################################################
## Check the inference

print(fit$time(), digits = 3)

pars = c("lp__", "CL", "Q", "VC", "VP", "ka", "sigma")
fit$summary(pars)
bayesplot::mcmc_trace(fit$draws(), pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)

##########################################################################
## Posterior predictive checks

# the predictive draws need to be turned into a matrix.
yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))
    ))[, -(52:54)]

yobs <- data$cObs
time <- data$time[-1]

# Bayesplot offers various functions we can experiment with.
bayesplot::ppc_intervals(y = yobs, yrep = yrep, x = time)
bayesplot::ppc_ribbon(y = yobs, yrep = yrep, x = time)

