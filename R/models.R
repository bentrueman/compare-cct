
# TO DO: rerun models "part", "diss_noar", and "sim_part" with 2000 sampling iterations per chain to meet min ESS

# fit GAMs to data and simulated data (models are loaded from CSVs by default):

#------------------ setup ------------------

source("R/formulas.R")
library("bgamcar1")
library("dplyr")

options(mc.cores = parallel::detectCores())

#------------------ inputs ------------------

model_in <- readr::read_csv("data-clean/model_in.csv") %>%
  mutate(location = as.factor(location))
model_in_sim <- readr::read_csv("data-clean/simulated-data.csv") %>%
  mutate(location = as.factor(location))

knots_yday <- c(0, 1)

#------------------ models ------------------

# n.b., each takes between 30 mins and 1.5 hrs to fit on a 2017 Macbook Pro

model_diss <- fit_stan_model(
  file = "models/model_diss",
  seed = stan_seed,
  bform = form_diss,
  bdata = model_in,
  bpriors = prior_diss,
  knots = list(date_yday = knots_yday),
  backend = "cmdstanr",
  save_warmup = FALSE,
  max_treedepth = 12
)

model_part <- fit_stan_model(
  file = "models/model_part",
  seed = stan_seed,
  bform = form_part,
  bdata = model_in,
  bpriors = prior_part,
  backend = "cmdstanr",
  iter_sampling = 2000,
  save_warmup = FALSE
)

model_diss_noar <- fit_stan_model(
  file = "models/model_diss_noar",
  seed = stan_seed,
  bform = form_diss_noar,
  bdata = model_in,
  bpriors = prior_diss[-1,],
  knots = list(date_yday = knots_yday),
  backend = "cmdstanr",
  save_warmup = FALSE,
  max_treedepth = 12,
  iter_sampling = 2000
)

model_part_noar <- fit_stan_model(
  file = "models/model_part_noar",
  seed = stan_seed,
  bform = form_part_noar,
  bdata = model_in,
  bpriors = prior_part[-1,],
  backend = "cmdstanr",
  save_warmup = FALSE
)

model_sim_diss <- fit_stan_model(
  file = "models/model_sim_diss",
  seed = stan_seed,
  bform = form_diss,
  bdata = model_in_sim,
  bpriors = prior_diss,
  knots = list(date_yday = knots_yday),
  backend = "cmdstanr",
  save_warmup = FALSE,
)

model_sim_part <- fit_stan_model(
  file = "models/model_sim_part",
  seed = stan_seed,
  bform = form_part_sim,
  bdata = model_in_sim,
  bpriors = prior_part,
  save_warmup = FALSE,
  backend = "cmdstanr",
  iter_sampling = 2000
)

model_noprior_diss <- fit_stan_model(
  file = "models/model_noprior_diss",
  seed = stan_seed,
  bform = form_diss,
  bdata = model_in,
  knots = list(date_yday = knots_yday),
  backend = "cmdstanr",
  save_warmup = FALSE,
  max_treedepth = 12
)

model_noprior_part <- fit_stan_model(
  file = "models/model_noprior_part",
  seed = stan_seed,
  bform = form_part,
  bdata = model_in,
  backend = "cmdstanr",
  save_warmup = FALSE
)
