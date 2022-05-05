
# prior predictive checks (models are loaded from CSVs by default):

#------------------ setup ------------------

library("bgamcar1")
source("R/formulas.R")
library("rstan")
library("ggplot2")
library("tibble")
library("tidyr")
library("dplyr")

options(mc.cores = parallel::detectCores())

#------------------ inputs ------------------

model_in <- readr::read_csv("data-clean/model_in.csv")

#------------------ models ------------------

knots_yday <- c(0, 1)

# prior simulation:

model_prior_part <- fit_stan_model(
  "data-clean/prior_part",
  seed = stan_seed,
  form_part, model_in, prior_part,
  save_warmup = FALSE,
  sample_prior = "only"
)

model_prior_diss <- fit_stan_model(
  "data-clean/prior_diss",
  seed = stan_seed,
  form_diss, model_in, prior_diss,
  save_warmup = FALSE,
  sample_prior = "only",
  knots = list(date_yday = knots_yday)
)

# smooths:

conditional_smooths(model_prior_part)
conditional_smooths(model_prior_diss)

model_prior_part$prior
model_prior_diss$prior
