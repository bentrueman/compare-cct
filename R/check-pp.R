
# generate a density overlay plot with all of the draws (runs in ~ 10 mins.)

#------------------ setup ------------------

source("R/models.R", echo = TRUE)
library("ggplot2")
library("tidyr")
library("dplyr")
library("purrr")
library("bayesplot")

theme_set(theme_bw() + theme(legend.position = "bottom"))
palette <- wesanderson::wes_palette("Zissou1")

#------------------ predict ------------------

preds <- map(
  list("diss" = model_diss, "part" = model_part),
  ~ add_pred_draws_car1(
    model_in, .x,
    # draw_ids = withr::with_seed(12478, {sample(1:4000, 50)}),
    type = "prediction"
  )
)

cens_bound <- model_in %>%
  summarize(
    part = min(scaled_lead_part[cens_lead_part == "left"]),
    diss = min(scaled_lead_dissolved[cens_lead_dissolved == "left"])
  )

yrep <- map2(
  preds, list("diss", "part"),
  ~ .x %>%
    ungroup() %>%
    # recensor:
    mutate(
      recensored = if_else(
        .prediction < cens_bound[[.y]], cens_bound[[.y]], .prediction
      ),
    )
)

#------------------ density plot ------------------

density_plot <- function(xvar, fraction) {
  model_in %>%
    ggplot(aes({{xvar}})) +
    geom_density(
      data = yrep[[fraction]] %>%
        filter(recensored < 10),
      aes(recensored, group = .draw),
      col = palette[2],
      size = .2,
      inherit.aes = FALSE
    ) +
    geom_density() +
    labs(x = "Observation/prediction", y = "Density")
}

density_plot(scaled_lead_dissolved, "diss") +
  coord_cartesian(xlim = c(cens_bound$diss, 10))

density_plot(scaled_lead_part, "part") +
  coord_cartesian(xlim = c(cens_bound$part, 10))

#------------------ proportion censored ------------------

preds_mat <- preds %>%
  map(
    ~ .x %>%
      pivot_wider(
        id_cols = c(location, date),
        names_from = .draw,
        values_from = .prediction
      ) %>%
      ungroup() %>%
      select(where(is.numeric)) %>%
      as.matrix() %>%
      t()
  )

bayesplot::ppc_stat(
  1 * (model_in$cens_lead_dissolved == "left"),
  1 * (preds_mat$diss < cens_bound$diss)
)

bayesplot::ppc_stat(
  1 * (model_in$cens_lead_part == "left"),
  1 * (preds_mat$part < cens_bound$part)
)

