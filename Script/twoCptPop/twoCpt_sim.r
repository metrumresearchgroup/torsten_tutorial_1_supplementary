
rm(list = ls())
gc()
set.seed(1954)
set.seed(11191951)

##.libPaths("~/Rlib")
##setwd("~/Code/torsten_tutorial_psp/Script/twoCptPop")
scriptDir <- getwd() ## assumes you navigated to twoCptPop directory

library(tidyverse)
library(cmdstanr)
library(posterior)

set_cmdstan_path(file.path(dirname(dirname(scriptDir)), "Torsten", "cmdstan"))

n_subjects <- 20
n_obs_persub <- 51
n_event <- n_obs_persub + 1
start <- c(1, 1:(n_subjects - 1)  * n_event + 1)

set.seed(1954)

time_after_dose <-
  c(0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)

## Create skeleton of NONMEM-style data set
data1 <- tibble(id = rep(1:n_subjects, each = n_event),
                time = rep(c(0, time_after_dose, 12 + time_after_dose,
                             seq(from = 24, to = 156, by = 12),
                             168 + time_after_dose, 180, 186, 192), n_subjects),
                addl = rep(c(14, rep(0, n_event - 1)), n_subjects),
                ii = rep(c(12, rep(0, n_event - 1)), n_subjects),
                amt  = rep(c(1200, rep(0, n_event - 1)), n_subjects),
                cmt = rep(c(1, rep(2, n_event - 1)), n_subjects),
                evid = rep(c(1, rep(0, n_event - 1)), n_subjects),
                rate = rep(0, n_event * n_subjects),
                ss = rep(0, n_event * n_subjects),
                ## use dummy value for observations
                cObs = rep(1, n_event * n_subjects))


## Reformat for Stan
start = c(1, 1:(n_subjects - 1)  * n_event + 1)
end = 1:n_subjects * n_event
iObs = (1:(n_event * n_subjects))[-start]
data <- c(as.list(data1 %>% select(-cObs)),
          list(start = start,
               end = end,
               iObs = iObs,
               nEvent = n_event * n_subjects,
               nSubjects = n_subjects,
               nIIV = 5,
               nObs = n_obs_persub * n_subjects,
               ## use dummy value for observations
               cObs = data1$cObs[iObs]))

init <- function () {
  CL_pop <- 10
  Q_pop <- 15
  VC_pop <- 35
  VP_pop <- 105
  ka_pop <- 2.5
  omega <- c(0.25, 0.5, 0.25, 0.5, 0.25)
  
  list(CL_pop = CL_pop, Q_pop = Q_pop, VC_pop = VC_pop, VP_pop = VP_pop,
       ka_pop = ka_pop, omega = omega,
       sigma = 0.1)   
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

data2 <- data1 %>%
  mutate(cObs = NA)
data2$cObs[iObs] = fit$draws("concentrationObsPred")[1, 1, ]

ggplot(data2, aes(x = time, y = cObs, group = id)) +
  geom_line()

write_csv(data2, "twoCptPop.data.csv")

