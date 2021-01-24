
rm(list = ls())
gc()
set.seed(1954)

# adjust to your settings
.libPaths("~/Rlib/")
setwd("~/Code/torsten_tutorial_psp/Script/twoCpt")

library(cmdstanr)
library(rjson)
library(bayesplot)
library(posterior)
library(loo)
library(ggplot2)

# set_cmdstan_path("../../Torsten/cmdstan/")
set_cmdstan_path("~/Code//torsten_tutorial_psp/Torsten/cmdstan/")
bayesplot::color_scheme_set("mix-blue-green")

##########################################################################
## Read in data and create inits.

data <- fromJSON(file = "twoCpt.data.json")

# Draw initial conditions from the prior
init <- function() {
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(15), 0.2)),
       VC = exp(rnorm(1, log(35), 0.2)),
       VP = exp(rnorm(1, log(105), 0.2)),
       ka = exp(rnorm(1, log(2), 0.2)),
       sigma = abs(rnorm(1, 0, 1)))
}

##########################################################################
## Compile and fit the model
model_name <- "twoCpt"
model_name2 <- "oneCpt"

mod <- cmdstan_model(paste0(model_name, ".stan"))

n_chains <- 4
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = n_chains,
                  iter_warmup = 500, iter_sampling = 500)


dir.create("deliv")
fit$save_object(paste0("deliv/", model_name, ".fit.RDS"))

# The saved fit object can be read using the below line
# fit <- readRDS(paste0("saved_fit/", model_name, ".fit.RDS"))

##########################################################################
## Check the inference

print(fit$time(), digits = 3)

pars = c("CL", "Q", "VC", "VP", "ka", "sigma")
# pars = c("lp__", "CL", "VC", "ka", "sigma")
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
bayesplot::ppc_ribbon(y = yobs, yrep = yrep, x = time, y_draw = "point")

# compute PSIS-loo estimate
log_lik_draws <- fit$draws("log_lik")
loo_estimate <- loo(log_lik_draws, r_eff = relative_eff(log_lik_draws))

print(loo_estimate)

##########################################################################
## Fit the model with a one compartment model and compare elpd_loo.

mod2 <- cmdstan_model(paste0(model_name2, ".stan"))

n_chains <- 4
fit2 <- mod2$sample(data = data, chains = n_chains, init = init,
                    parallel_chains = n_chains,
                    iter_warmup = 500, iter_sampling = 500)


fit2$save_object(paste0("deliv/", model_name2, ".fit.RDS"))

log_lik_draws2 <- fit2$draws("log_lik")
loo_estimate2 <- loo(log_lik_draws2, r_eff = relative_eff(log_lik_draws2))

## ppc for second model
yrep2 <- as.matrix(
  as_draws_df(
    fit2$draws(variables = c("concentrationObsPred"))
  ))[, -(52:54)]

bayesplot::ppc_ribbon(y = yobs, yrep = yrep2, x = time)
bayesplot::ppc_ribbon(y = yobs, yrep = yrep2, x = time, y_draw = "point") +
  xlab("Time (h)") + ylab("Drug concentration (mg/L)")

print(loo_estimate2)

# save loo estimates in table
loo_results <- data.frame(rbind(loo_estimate$estimate[1, ],
                                loo_estimate2$estimate[1, ]))
loo_results$model <- c(model_name, model_name2)
loo_results

loo_results$model <- c("two compartment", "one compartment")

p <- ggplot(loo_results, aes(x = model, y = Estimate)) + theme_bw() +
  geom_pointrange(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                  colour = "dark green") + ylab("elpd loo")
p + theme(text = element_text(size = 16))
