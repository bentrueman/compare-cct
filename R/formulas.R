
# write formulas for use in "R/models.R":

library("brms")
library("stringr")

stan_seed <- 236 # random seed

# formulas:

ar <- 'ar(time = date_numeric, gr = location)'
seasonal <- 's(date_yday, bs = "cc")'

form_part_noar <- 'scaled_lead_part | cens(cens_lead_part, y2 = scaled_lead) ~
    1 + (1 | location) +
    s(date_numeric) +
    s(date_numeric, by = location, m = 1)'

form_diss_noar <- paste(
  str_replace_all(form_part_noar, "lead_part", "lead_dissolved") %>% 
    str_remove_all(", y2 = scaled_lead"), 
  seasonal, sep = " + "
)

form_part <- bf(paste(form_part_noar, ar, sep = " + "))
form_diss <- bf(paste(form_diss_noar, ar, sep = " + "))
form_part_sim <- bf(
  paste(str_remove_all(form_part_noar, ", y2 = scaled_lead"), ar, sep = " + ")
)

form_part_noar <- bf(form_part_noar)
form_diss_noar <- bf(form_diss_noar)

# prior:

prior_diss <- c(
  prior(normal(0.5, .25), class = ar, lb = 0, ub = 1),
  prior(student_t(3, 0, 2.5), class = b),
  prior(gamma(2, .5), class = nu),
  prior(normal(0, .5), class = sds, coef = `s(date_yday, bs = "cc")`)
)

prior_part <- prior_diss[-4,] # remove seasonal smooth
