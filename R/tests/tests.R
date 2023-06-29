
# this script includes:
# (1) a few additional tests of the "bgamcar1" functions
# (2) model diagnostics

#------------------ tests ------------------

library("testthat")
source("R/models.R", echo = TRUE)
library("purrr")
library("dplyr")
library("posterior")
library("rstan")

#------------------ tests ------------------

# loo_cv():

loo1 <- loo(model_diss, incl_autocor = FALSE)
loo2 <- loo_cv(model_in, model_diss, car1 = FALSE)

loo3 <- loo(model_part, incl_autocor = FALSE)
loo4 <- loo_cv(model_in, model_part, car1 = FALSE)

test_that("brms::loo() yields the same results as loo_cv()", {
  expect_equal(loo1$estimate, loo2$estimate)
  expect_equal(loo3$estimate, loo4$estimate)
})

# --------------- model diagnostics ---------------

model_list <- list(
  # full models:
  "diss" = model_diss,
  "part" = model_part,
  # no CAR1 term:
  "diss_noar" = model_diss_noar,
  "part_noar" = model_part_noar,
  # simulated data
  "sim_diss" = model_sim_diss,
  "sim_part" = model_sim_part,
  # default priors:
  "diss_noprior" = model_noprior_diss,
  "part_noprior" = model_noprior_part
)

sampler_params <- model_list %>%
  map_dfr(
    # warmup draws should all be excluded b/c they aren't saved in the CSV files
    ~ get_sampler_params(.x$fit, inc_warmup = FALSE) %>%
      map(as_tibble) %>%
      bind_rows(.id = ".chain") %>%
      mutate(max_treedepth = .x$fit@stan_args[[4]]$max_depth),
    .id = "model"
  )

draw_summary <- model_list %>%
  map_dfr(
    ~ as_draws(.x) %>%
      posterior::summarise_draws(),
    .id = "model"
  )

test_that("max_treedepth never exceeded", {
  exceedances <- sampler_params %>%
    filter(treedepth__ >= max_treedepth)
  n <- nrow(exceedances)
  expect_equal(n, 0)
})

test_that("no divergent transitions", {
  divergences <- sampler_params %>%
    filter(divergent__ != 0)
  n <- nrow(divergences)
  expect_equal(n, 0)
})

test_that("no rhats greater than 1.05", {
  hi_rhat <- draw_summary %>%
    filter(rhat >= 1.05)
  n <- nrow(hi_rhat)
  expect_equal(n, 0)
})

test_that("no low ESS (bulk)", {
  lo_ess <- draw_summary %>%
    filter(ess_bulk < 400)
  n <- nrow(lo_ess)
  expect_equal(n, 0)
})

test_that("no low ESS (tail)", {
  lo_ess <- draw_summary %>%
    filter(ess_tail < 400)
  n <- nrow(lo_ess)
  expect_equal(n, 0)
})
