
# approximate leave-one-out cross-validation for models (runs in ~ 10 mins):

source("R/models.R", echo = TRUE)
library("dplyr")

#------------------ loo ------------------

loo_diss_car <- loo_cv(model_in, model_diss)
loo_diss_noar <- loo_cv(model_in, model_diss_noar, car1 = FALSE)

loo_part_car <- loo_cv(model_in, model_part)
loo_part_noar <- loo_cv(model_in, model_part_noar, car1 = FALSE)

# CAR(1) models fit better:
loo_compare(loo_diss_car, loo_diss_noar)
loo_compare(loo_part_car, loo_part_noar)

# write pareto-k values to file:
model_in %>%
  transmute(
    location, ortho_dose, pipe_material,
    date, date_numeric, 
    pk_vals_diss = loo::pareto_k_influence_values(loo_diss_car),
    pk_vals_part = loo::pareto_k_influence_values(loo_part_car)
  ) %>%
  readr::write_csv("data-clean/pareto-k-values.csv")
