
# simulate data from the final model for fitting a test model:

source("R/models.R", echo = TRUE)
library("dplyr")
library("purrr")
library("withr")
library("tidyr")
library("readr")
library("assertr", include.only = "verify")

#------------------ simulate ------------------

preds_raw <- list("lead_dissolved" = model_diss, "lead_part" = model_part) %>% 
  map_dfr(
    ~ withr::with_seed(432, {
      add_pred_draws_car1(
        model_in, .x,
        draw_ids = 3386, type = "prediction"
      )
    }),
    .id = "type"
  )

preds_raw %>%
  ungroup() %>%
  distinct(type, `ar[1]`, sigma, nu)

preds <- preds_raw %>% 
  pivot_wider(
    id_cols = names(model_in),
    names_from = type, values_from = .prediction, names_prefix = ".prediction_"
  )

cens_lim_scaled <- model_in %>%
  filter(cens_lead_dissolved == "left", lead_dissolved == min(lead_dissolved)) %>%
  distinct(lead_dissolved, scaled_lead_dissolved)

simdat <- preds %>%
  ungroup() %>%
  mutate(
    # modify censoring indicators:
    across(
      starts_with(".prediction_"),
      ~ if_else(.x < cens_lim_scaled$scaled_lead_dissolved, "left", "none"), 
      .names = "cens_{str_remove(.col, '.prediction_')}"
    ),
    # recensor:
    across(
      starts_with(".prediction_"), 
      ~ pmax(.x, cens_lim_scaled$scaled_lead_dissolved),
      .names = "scaled_{str_remove(.col, '.prediction_')}"
    ),
  )

# final checks:

simdat %>%
  assertr::verify(scaled_lead_dissolved < 10 * max(model_in$scaled_lead_dissolved)) %>%
  assertr::verify(scaled_lead_part < 10 * max(model_in$scaled_lead_part)) %>%
  assertr::verify(length(scaled_lead_part) == length(model_in$scaled_lead_part)) %>%
  assertr::verify(!is.na(scaled_lead_part)) %>%
  assertr::verify(!is.na(scaled_lead_dissolved)) %>%
  # censored values are all 2e-4, uncensored greater:
  assertr::verify(
    (cens_lead_dissolved == "left" & scaled_lead_dissolved == cens_lim_scaled$scaled_lead_dissolved) |
      (cens_lead_dissolved == "none" & scaled_lead_dissolved > cens_lim_scaled$scaled_lead_dissolved)
  ) %>%
  assertr::verify(
    (cens_lead_part == "left" & scaled_lead_part == cens_lim_scaled$scaled_lead_dissolved) |
      (cens_lead_part == "none" & scaled_lead_part > cens_lim_scaled$scaled_lead_dissolved)
  ) %>%
  readr::write_csv("data-clean/simulated-data.csv")
