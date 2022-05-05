
# generate a density overlay plot with all of the draws (runs in ~ 10 mins.)

#------------------ setup ------------------

source("R/models.R", echo = TRUE)
library("ggplot2")
library("dplyr")

theme_set(theme_bw() + theme(legend.position = "bottom"))
palette <- wesanderson::wes_palette("Zissou1")

#------------------ predict ------------------

preds <- add_pred_draws_car1(
  model_in, model_diss,
  # draw_ids = withr::with_seed(12478, {sample(1:4000, 50)}),
  type = "prediction"
)

cens_bound <- model_in %>%
  filter(cens_lead_dissolved == "left") %>%
  summarize(bound = min(scaled_lead_dissolved))

yrep <- preds %>%
  ungroup() %>%
  # recensor:
  mutate(scaled_lead_dissolved = if_else(.prediction < cens_bound$bound, cens_bound$bound, .prediction)) %>%
  filter(scaled_lead_dissolved < 10)

#------------------ plot ------------------

model_in %>%
  ggplot(aes(scaled_lead_dissolved)) +
  geom_density(
    data = yrep,
    aes(group = .draw),
    col = palette[2],
    size = .2
  ) +
  geom_density() +
  coord_cartesian(xlim = c(-2, 10)) +
  labs(x = "Observation/prediction", y = "Density")

