
# fit GAMs to data and simulated data (models are loaded from CSVs by default):

#------------------ setup ------------------

source("R/formulas.R")
library("rstan")
library("bgamcar1")

options(mc.cores = parallel::detectCores())

#------------------ inputs ------------------

model_in <- readr::read_csv("data-clean/model_in.csv")
model_in_sim <- readr::read_csv("data-clean/simulated-data.csv")

knots_yday <- c(0, 1)

#------------------ models ------------------

# n.b., each takes ~30 min to fit on a 2017 Macbook Pro

model_diss <- fit_stan_model(
  "models/model_diss",
  seed = stan_seed,
  form_diss, model_in, prior_diss,
  save_warmup = FALSE,
  control = list(max_treedepth = 12),
  knots = list(date_yday = knots_yday)
)

model_part <- fit_stan_model( 
  "models/model_part",
  seed = stan_seed,
  form_part, model_in, prior_part, 
  save_warmup = FALSE
)

model_diss_noar <- fit_stan_model(
  "models/model_diss_noar",
  seed = stan_seed,
  form_diss_noar, model_in, prior_diss[-1,],
  save_warmup = FALSE,
  control = list(max_treedepth = 12), 
  knots = list(date_yday = knots_yday)
)

model_part_noar <- fit_stan_model(
  "models/model_part_noar",
  seed = stan_seed,
  form_part_noar, model_in, prior_part[-1,],
  save_warmup = FALSE
)

model_sim_diss <- fit_stan_model(
  "models/model_sim_diss",
  seed = stan_seed,
  form_diss, model_in_sim, prior_diss,
  save_warmup = FALSE,
  knots = list(date_yday = knots_yday)
)

model_sim_part <- fit_stan_model(
  "models/model_sim_part",
  seed = stan_seed,
  form_part_sim, model_in_sim, prior_part,
  save_warmup = FALSE
)

model_noprior_diss <- fit_stan_model(
  "models/model_noprior_diss",
  seed = stan_seed,
  form_diss, model_in, 
  save_warmup = FALSE,
  control = list(max_treedepth = 12), 
  knots = list(date_yday = knots_yday) 
)

model_noprior_part <- fit_stan_model(
  "models/model_noprior_part",
  seed = stan_seed,
  form_part, model_in, 
  save_warmup = FALSE
)
