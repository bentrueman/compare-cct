
# generate model inputs:

library("dplyr")
library("tidyr")
library("stringr")
library("testthat")
library("assertr", include.only = "verify")

#------------------ functions ------------------

# used for scaling the upper bound of censored intervals:

scale_fun <- function(x, y = NULL, log_trans = FALSE) {
  if(log_trans) {
    x <- log(x)
  }
  if(is.null(y)) {
    y <- x
  }
  m <- mean(y)
  sdev <- sd(y)
  
  (x - m) / sdev
  
}

test_that("scale_fun() returns expected result", {
  x <- runif(10, 1, 2)
  t1 <- scale_fun(x, log_trans = TRUE)
  t2 <- scale(log(x))[,1]
  t3 <- scale_fun(x, y = log(x), log_trans = TRUE)
  t4 <- scale_fun(x)
  t5 <- scale(x)[,1]
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t4, t5)
})

#------------------ data cleaning ------------------

pdat <- readr::read_csv("data-clean/pdat.csv")

pdat_clean <- pdat %>%
  mutate(param = str_replace(param, "Lead, dissolved", "Lead Dissolved")) %>%
  filter(
    str_detect(location, "^LP[1-3][1235]$"),
    param %in% c("Lead", "Lead Dissolved"),
    !is.na(value)
  ) %>%
  transmute(
    type = str_to_lower(param) %>%
      str_replace("\\s", "_"),
    ortho_dose,
    pipe_material,
    date,
    location,
    value,
    bdl
  ) %>%
  distinct() %>%
  assertr::verify(value >= 2e-4 | (value == 2e-4 & bdl)) %>%
  assertr::verify(date > "2017-01-01") %>%
  assertr::verify(date <= lubridate::today())

pdat_part <- pdat_clean %>%
  pivot_wider(names_from = type, values_from = c(value, bdl), values_fn = median) %>% 
  # the median of a censored and uncensored value should be censored at the median:
  mutate(
    # code as censored any median-aggregated value where the median is censored 
    # or the midpoint of a censored and uncensored value:
    across(starts_with("bdl"), ~ as.logical(if_else(.x >= .5, 1, 0))),
    # calculate particulate concentration
    value_lead_part = value_lead - value_lead_dissolved,
    cens_lead_part = case_when(
      # censor any particulate concentrations < DL for dissolved:
      value_lead_part < min(value_lead_dissolved) ~ "left", 
      !bdl_lead & !bdl_lead_dissolved ~ "none",
      !bdl_lead & bdl_lead_dissolved ~ "interval",
      bdl_lead & !bdl_lead_dissolved ~ "left",
      bdl_lead & bdl_lead_dissolved ~ "left"
    ),
    value_lead_part = pmax(value_lead_part, min(value_lead_dissolved)),
    cens_lead_dissolved = if_else(bdl_lead_dissolved, "left", "none") 
  ) %>% 
  select(-c(bdl_lead, bdl_lead_dissolved))

test_that("Cu/Pb-Sn pipes alone have > 50% censoring", {
  hi_cens <- pdat_part %>% 
    group_by(ortho_dose, pipe_material) %>% 
    summarize(
      across(starts_with("cens_"), ~ 1 - mean(.x == "none"), .names = "pcent_{.col}"),
      .groups = "drop"
    ) %>% 
    filter(if_all(starts_with("pcent_"), ~ .x > .5))
  expect_equal(unique(hi_cens$pipe_material), "Cu/Pb-Sn")
  expect_equal(unique(hi_cens$ortho_dose), c("0-0.5", "1.0", "2.0-0.75"))
})

pdat_part <- pdat_part %>% 
  filter(pipe_material != "Cu/Pb-Sn")

model_in <- pdat_part %>%
  rename_all(~ str_remove(.x, "value_")) %>% 
  mutate(
    # normalized day of year:
    ndays = if_else(lubridate::leap_year(date), 366, 365),
    date_yday = lubridate::yday(date) / ndays, # fraction of year
    # convert date to numeric date:
    year = lubridate::year(date),
    year_counter = year - min(year),
    date_numeric = date_yday + year_counter,
    # scaled variables:
    across(starts_with("lead_"), ~ scale_fun(.x, log_trans = TRUE), .names = "scaled_{.col}"),
    scaled_lead = scale_fun(lead, y = log(lead_part), log_trans = TRUE)
  ) %>%
  arrange(location, date) %>%
  # determine CAR(1) exponent:
  group_by(location) %>%
  mutate(
    d_x = date - dplyr::lag(date),
    d_x = as.numeric(d_x),
    d_x = d_x / 7, # most common time difference
    d_x = replace_na(d_x, 0) # NA is due to dplyr::lag()
  ) %>%
  ungroup() %>%
  select(-c(ndays, year, year_counter)) %>% 
  assertr::verify(
    (scaled_lead > scaled_lead_part & cens_lead_part == "interval") | 
      scaled_lead >= scaled_lead_part
  ) %>% 
  assertr::verify(d_x > 0 | (d_x == 0 & date == min(date))) %>%
  assertr::verify(lead_part >= 0 & lead_part < 1e4) %>% 
  assertr::verify(lead_dissolved >= 2e-4 & lead_dissolved < 1e4)

# tests:

test_that("Current version matches previous version.", {
  
  df_test <- list(
    "pdat_clean_ref" = readr::read_csv("data-clean/pdat_clean.csv"),
    "model_in_ref" = readr::read_csv("data-clean/model_in.csv"),
    "pdat_clean" = pdat_clean,
    "model_in" = model_in
  ) %>% 
    purrr::map(
      ~ .x %>% 
        arrange(date, location) %>% 
        data.frame()
    )
  
  expect_equal(df_test$pdat_clean_ref, df_test$pdat_clean)
  expect_equal(df_test$model_in_ref, df_test$model_in)
  
})

# write:

readr::write_csv(pdat_clean, "data-clean/pdat_clean.csv")
readr::write_csv(model_in, "data-clean/model_in.csv")
