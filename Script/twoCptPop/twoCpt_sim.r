
rm(list = ls())
gc()
set.seed(1954)

.libPaths("~/Rlib")
setwd("~/Code/torsten_tutorial_psp/Script/twoCptPop")

library(cmdstanr)
library(rjson)
set_cmdstan_path("../../Torsten/cmdstan/")

n_subjects <- 3
n_obs_persub <- 51
n_event <- n_obs_persub + 1
start <- c(1, 1:(n_subjects - 1)  * n_event + 1)

set.seed(1954)

time_after_dose <-
  c(0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)

data <- list(
  addl = rep(c(14, rep(0, n_event - 1)), n_subjects),
  amt  = rep(c(1200, rep(0, n_event - 1)), n_subjects),
  cmt = rep(c(1, rep(2, n_event - 1)), n_subjects),
  evid = rep(c(1, rep(0, n_event - 1)), n_subjects),
  rate = rep(0, n_event * n_subjects),
  ss = rep(0, n_event * n_subjects),

  # use dummy value for observations
  cObs = rep(1, n_subjects * n_obs_persub),

  time = rep(c(0, time_after_dose, 12 + time_after_dose,
               seq(from = 24, to = 156, by = 12),
               168 + time_after_dose, 180, 186, 192), n_subjects),

  start = c(1, 1:(n_subjects - 1)  * n_event + 1),
  end = 1:n_subjects * n_event,
  ii = rep(c(12, rep(0, n_event - 1)), n_subjects),
  iObs = (1:(n_event * n_subjects))[-start],
  nEvent = n_event * n_subjects,
  nSubjects = n_subjects,
  nIIV = 5,
  nObs = n_obs_persub * n_subjects
)


# draw parameters from prior
init <- function () {
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
## Simulate data with Stan
model_name <- "twoCptPop"
mod <- cmdstan_model(paste0(model_name, ".stan"))

n_chains <- 1
fit <- mod$sample(data = data, chains = n_chains, init = init,
                  parallel_chains = n_chains,
                  iter_warmup = 0, iter_sampling = 1,
                  fixed_param = TRUE, seed = 123)

yrep <- as.matrix(
  as_draws_df(
    fit$draws(variables = c("concentrationObsPred"))
  ))[1, 1:data$nObs]

plot(x = data$time[data$iObs], y = yrep)

data$cObs <- unname(yrep)

write_stan_json(data, "twoCptPop.data.json")
