
# generate figures for "paper.Rmd" and "paper-si.Rmd"

#------------------ setup ------------------

source("R/models.R", echo = TRUE)
library("forcats")
library("ggplot2")
library("dplyr")
library("tidyr")
library("purrr")
library("patchwork")
library("readr")
library("posterior")
library("withr", include.only = "with_seed")
library("ggdist", include.only = c("stat_halfeye", "median_qi"))
library("glue", include.only = "glue")
library("assertr", include.only = "verify")
library("devtools", include.only = "session_info")
library("janitor", include.only = "clean_names")
library("wesanderson", include.only = "wes_palette")
library("ggh4x", include.only = "facetted_pos_scales")
library("mgcv", include.only = "gam")
library("broom", include.only = "tidy")
library("lubridate", include.only = c("year", "yday"))
library("ggtext", include.only = c("element_markdown", "geom_richtext"))
library("tibble", include.only = c("tribble", "rowid_to_column"))

theme_set(theme_bw() + theme(legend.position = "bottom"))
palette <- wesanderson::wes_palette("Zissou1", 6, "continuous")

#------------------ load ------------------

pdat <- readr::read_csv("data-clean/pdat.csv")
pdat_clean <- readr::read_csv("data-clean/pdat_clean.csv")
pk_vals <- readr::read_csv("data-clean/pareto-k-values.csv")

locations <- model_in %>%
  distinct(location, ortho_dose, pipe_material)

#------------------ figure 1 ------------------

annotations <- tibble::tribble(
  ~date, ~value, ~label,
  "2018-12-01", 0.01, "1. Non-linearity",
  "2020-10-01", 1e-4, "2. Left-censoring",
  "2020-03-01", 0.04, "3. Autocorrelation\n(clustering of similar\nvalues in time)",
  "2018-04-01", .7, "4. Irregular sampling\nfrequency",
  "2021-05-01", .5, "5. Extreme\nvalues",
) %>%
  mutate(date = as.Date(date))

fig1 <- model_in %>%
  filter(location == "LP31") %>%
  pivot_longer(c(lead_part, lead_dissolved), names_to = "type") %>%
  ggplot(aes(date, value)) +
  geom_rug(sides = "t", color = "grey", length = unit(.02, "npc")) +
  geom_label(
    data = annotations,
    aes(label = label),
    alpha = .8, label.r = unit(0, "cm"),
    label.size = 0,
    size = 2
  ) +
  geom_line(aes(col = type)) +
  scale_y_log10(
    limits = c(1e-4, 2),
    labels = function(breaks) 1e3 * breaks
  ) +
  scale_color_manual(
    values = palette[c(6, 1)],
    labels = c("<0.45 µm", ">0.45 µm")
  ) +
  labs(x = NULL, y = expression("[Pb] (µg L"^-1*")"), col = NULL)

ggsave("Rmarkdown/figures/figure-1.png", fig1, device = "png", width = 3.33, height = 2.5, dpi = 600)

#------------------ figure 2 ------------------

linecol <- "black" # points/lines

models_initial_in <- pdat %>%
  filter(
    str_detect(location, "LP11"),
    param == "Lead Dissolved"
  ) %>%
  group_by(date) %>%
  summarize(value = median(value), bdl = mean(bdl)) %>%
  ungroup()

models_initial <- models_initial_in %>%
  mutate(
    numeric_date = as.numeric(date),
    numeric_date = numeric_date - min(numeric_date) + 1
  ) %>%
  nest(data = everything()) %>%
  mutate(
    # models:
    constant = map(data, ~ lm(log(value) ~ 1, data = .x)),
    linear = map(data, ~ lm(log(value) ~ numeric_date, data = .x)),
    gam = map(data, ~ mgcv::gam(log(value) ~ s(numeric_date, k = 10, bs = "cr"), data = .x)),
    # fitted:
    fit_const = map(constant, fitted),
    fit_lin = map(linear, fitted),
    fit_gam = map(gam, fitted),
    # residuals:
    res_const = map(constant, residuals),
    res_lin = map(linear, residuals),
    res_gam = map(gam, residuals)
  )

p1 <- models_initial %>%
  unnest(c(data, starts_with("fit_"), starts_with("res_"))) %>%
  select(where(~ !is.list(.x))) %>%
  pivot_longer(
    c(starts_with("fit_"), starts_with("res_")),
    names_to = c(".value", "model"),
    names_pattern = "(.+)_(.+)"
  ) %>%
  pivot_longer(c(fit, res), values_to = "out") %>%
  mutate(
    model = fct_recode(model, "GAM" = "gam", "Constant" = "const", "Linear" = "lin"),
    model = fct_relevel(model, "GAM", after = Inf),
    name = fct_recode(name, "'['*Pb*']'~'('*mg~L^-1*')'" = "fit", "Residual" = "res")
  ) %>%
  ggplot(aes(date)) +
  facet_grid(
    rows = vars(name), cols = vars(model),
    scales = "free_y",
    labeller = label_parsed
  ) +
  geom_line(
    data = function(x) {
      x %>%
        filter(name != "Residual")
    },
    aes(y = log(value)),
    col = linecol,
    size = .4
  ) +
  geom_line(
    data = function(x) {
      x %>%
        filter(name != "Residual")
    },
    aes(y = out),
    col = palette[1], size = 1
  ) +
  geom_line(
    data = function(x) {
      x %>%
        filter(name == "Residual")
    },
    aes(y = out),
    col = linecol,
    size = .4
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      name == "'['*Pb*']'~'('*mg~L^-1*')'" ~ scale_y_continuous(
        breaks = log(c(1e-3, 1e-2)),
        labels = function(breaks) exp(breaks)
      )
    )
  ) +
  labs(x = NULL, y = NULL)

p2 <- models_initial %>%
  unnest(c(data, starts_with("res_"))) %>%
  select(where(~ !is.list(.x))) %>%
  pivot_longer(
    starts_with("res_"),
    values_to = "residual",
    names_to = "model"
  ) %>%
  arrange(model, date) %>%
  group_by(model) %>%
  mutate(lag1 = lag(residual)) %>%
  ungroup() %>%
  assertr::verify(sum(is.na(lag1)) == 3) %>%
  filter(!is.na(lag1)) %>%
  mutate(
    model = fct_recode(model, "GAM" = "res_gam", "Constant" = "res_const", "Linear" = "res_lin"),
    model = fct_relevel(model, "GAM", after = Inf)
  ) %>%
  ggplot(aes(lag1, residual)) +
  facet_wrap(vars(model)) +
  geom_point(alpha = .5, shape = 16) +
  geom_smooth(method = "lm", col = palette[1], se = FALSE) +
  geom_label(
    data = function(x) {
      x %>%
        group_by(model) %>%
        summarize(
          r = cor(residual, lag1, use = "complete"),
          txt = glue::glue("italic(r)~'='~{signif(r, 2)}"),
        )
    },
    aes(x = Inf, y = -Inf, label = txt),
    hjust = "inward", vjust = "inward",
    label.size = 0, col = palette[6],
    parse = TRUE
  ) +
  labs(
    x = expression(italic(epsilon["t-1"])),
    y = expression(italic(epsilon[t]))
  )

p3 <- models_initial %>%
  unnest(starts_with("res_")) %>%
  select(where(~ !is.list(.x))) %>%
  pivot_longer(
    starts_with("res_"),
    values_to = "residual",
    names_to = "model"
  ) %>%
  group_by(model) %>%
  summarize(
    acf_raw = list(acf(residual, plot = FALSE)),
    acf_tidy = map(acf_raw, broom::tidy),
    bound = 1.96 / sqrt(length(residual))
  ) %>%
  ungroup() %>%
  unnest(acf_tidy) %>%
  mutate(
    model = fct_recode(model, "GAM" = "res_gam", "Constant" = "res_const", "Linear" = "res_lin"),
    model = fct_relevel(model, "GAM", after = Inf)
  ) %>%
  ggplot(aes(lag, acf)) +
  facet_wrap(vars(model)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = -bound, ymax = bound), alpha = .25) +
  geom_line() +
  geom_point() +
  geom_point(
    data = function(x) {
      x %>%
        filter(lag == 1)
    },
    aes(col = "italic(r)~from~panel~bold(b)"),
    size = 2
  ) +
  scale_color_manual(values = palette[6], labels = function(breaks) parse(text = breaks)) +
  scale_x_continuous(labels = function(breaks) paste0("k = ", breaks)) +
  labs(
    x = expression(italic(epsilon["t-k"])),
    y = expression("Correlation with" ~ italic(epsilon[t])),
    col = NULL
  ) +
  theme(
    legend.position = c(.9, .8),
    legend.box.margin = margin(),
    legend.margin = margin()
  )

fig2 <- patchwork::wrap_plots(p1, p2, p3, heights = c(2, 1, 1)) +
  patchwork::plot_annotation(tag_level = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("Rmarkdown/figures/figure-2.png", fig2, width = 7, height = 6, dpi = 600)

#------------------ figure 3 ------------------

gam_basis <- models_initial %>%
  mutate(
    coefs = map(gam, ~ coefficients(.x)),
    lpmat = map(gam, ~ predict(.x, type = "lpmatrix")),
    # weight basis functions by coefficients:
    lpmat_wt = map2(lpmat, coefs, ~ apply(.x, 1, function(u) u * .y) %>% t()),
    lpmat = map(lpmat, as_tibble),
    lpmat_wt = map(lpmat_wt, as_tibble),
  )

fig3 <- gam_basis %>%
  select(data, lpmat, lpmat_wt) %>%
  pivot_longer(-data) %>%
  unnest(c(data, value)) %>%
  pivot_longer(
    matches("\\(.+\\)"),
    names_to = "basis_name",
    values_to = "basis_value"
  ) %>%
  filter(basis_name != "(Intercept)") %>%
  mutate(
    name = fct_recode(name,
      "Unweighted basis functions of x" = "lpmat",
      "Weighted basis functions of x" = "lpmat_wt",
    )
  ) %>%
  ggplot(aes(date)) +
  facet_wrap(vars(name), scales = "free_y", ncol = 1) +
  geom_line(aes(y = basis_value, col = basis_name)) +
  geom_line(
    data = function(x) {
      x %>%
        group_by(name, date) %>%
        summarize(pred = sum(basis_value)) %>%
        ungroup() %>%
        filter(!str_detect(name, "Unweighted"))
    },
    aes(y = pred, col = "sum(italic(b[j]*'('*x[t]*')'*beta[j]), j == 1, k)"),
    size = 1
  ) +
  scale_color_manual(
    values = c(wesanderson::wes_palette("Zissou1", 9, type = "continuous"), "black"),
    labels = function(breaks) {
      breaks %>%
        str_replace("s\\(numeric_date\\)\\.", "b") %>%
        str_replace("(b)(\\d)", "italic(\\1[\\2]*'('*x[t]*')')") %>%
        parse(text = .)
    }
  ) +
  theme(
    plot.margin = margin(r = .4, l = .1, b = .1, t = .1, unit = "cm"),
    legend.text.align = 0
  ) +
  labs(
    x = NULL,
    y = expression("log([Pb] (µg L"^-1 * ")"),
    col = NULL
  ) +
  guides(col = guide_legend(nrow = 4))

ggsave("Rmarkdown/figures/figure-3.png", fig3, width = 3.33, height = 4, dpi = 600)

#------------------ figure 4 ------------------

n_smooths <- 200

smooth_terms_part <- conditional_smooths(model_part)
smooth_terms_diss <- conditional_smooths(model_diss)

# repeated for supplementary figures of the models fitted to simulated data:

calc_spag_smooths <- function(model_diss, model_part) {
  expand_grid(
    model = list(model_diss, model_part),
    smooth_term = c(
      "s(date_numeric)",
      "s(date_yday,bs=\"cc\")",
      "s(date_numeric,by=location,m=1)"
    )
  ) %>%
    mutate(resp = map(model, ~ attributes(.x$fit)$model_name)) %>%
    unnest(resp) %>%
    filter(!str_detect(resp, "_part$") | !smooth_term == "s(date_yday,bs=\"cc\")") %>%
    mutate(
      type = if_else(str_detect(resp, "_diss$"), "<0.45 µm", ">0.45 µm"),
      smooth = map2(
        model, smooth_term,
        ~ posterior_smooths(
          .x,
          smooth = .y,
          draw_ids = withr::with_seed(2149, {sample(1:4000, n_smooths)})
        ) %>%
          as_tibble(rownames = ".draw") %>%
          pivot_longer(-.draw, names_to = ".row"),
        .id = "smooth"
      )
    ) %>%
    select(-model)
}

smooth_terms_spag <- calc_spag_smooths(model_diss, model_part)

plot_fig4 <- function(smooth_terms_diss, smooth_terms_part, smooth_terms_spag) {

  fig4a <- smooth_terms_diss[["mu: s(date_yday,bs=\"cc\")"]] %>%
    ggplot(aes(52 * date_yday)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = .5) +
    geom_line(
      data = smooth_terms_spag %>%
        filter(smooth_term == "s(date_yday,bs=\"cc\")", str_detect(resp, "_diss$")) %>%
        unnest(smooth) %>%
        mutate(date_yday = rep(model_in$date_yday, n_smooths)) %>%
        distinct(smooth_term, .draw, date_yday, value),
      aes(y = value, group = .draw),
      col = palette[2], size = .05
    ) +
    geom_line(aes(y = estimate__)) +
    labs(x = "Week", y = "Effect\n(transformed)")

  fig4b <- list(
    ">0.45 µm" = smooth_terms_part[["mu: s(date_numeric)"]],
    "<0.45 µm" = smooth_terms_diss[["mu: s(date_numeric)"]]
  ) %>%
    bind_rows(.id = "type") %>%
    mutate(date_numeric = date_numeric + 2017) %>%
    ggplot(aes(date_numeric)) +
    facet_wrap(vars(type), ncol = 2) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = .5) +
    geom_line(
      data = smooth_terms_spag %>%
        filter(smooth_term == "s(date_numeric)") %>%
        unnest(smooth) %>%
        mutate(date_numeric = rep(model_in$date_numeric, 2 * n_smooths) + 2017),
      aes(y = value, group = .draw),
      col = palette[2], size = .05
    ) +
    geom_line(aes(y = estimate__)) +
    labs(x = "Year", y = "Effect\n(transformed)")

  fig4c <- list(
      ">0.45 µm" = smooth_terms_part[["mu: s(date_numeric,by=location,m=1)"]],
      "<0.45 µm" = smooth_terms_diss[["mu: s(date_numeric,by=location,m=1)"]]
    ) %>%
    bind_rows(.id = "type") %>%
    mutate(date_numeric = date_numeric + 2017) %>%
    left_join(locations, by = "location") %>%
    ggplot(aes(date_numeric, col = ortho_dose)) +
    facet_grid(rows = vars(type), cols = vars(pipe_material), scales = "free_y") +
    geom_line(aes(y = estimate__)) +
    labs(
      x = "Year",
      y = "Effect\n(transformed scale)",
      col = expression("mg P L"^-1)
    ) +
    scale_y_continuous(n.breaks = 4)

  patchwork::wrap_plots(
    fig4a, fig4b, fig4c,
    design = "AB\nCC",
    heights = c(2, 3),
    widths = c(1, 2)
  ) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_level = "a") &
    theme(
      plot.tag = element_text(face = "bold"),
      legend.position = "bottom"
    ) &
    scale_color_manual(values = palette[c(6,4,1)]) &
    scale_fill_manual(values = palette[c(6,4,1)])
}

fig4 <- plot_fig4(smooth_terms_diss, smooth_terms_part, smooth_terms_spag)

ggsave("Rmarkdown/figures/figure-4.png", fig4, width = 7.5, height = 4.75, dpi = 600)

#------------------ figure 5 ------------------

drawids <- 3429

pp_part <- withr::with_seed(214,
  {add_pred_draws_car1(model_in, model_part,  draw_ids = drawids, type = "prediction")}
) %>%
  impute_censored(model_in, "scaled_lead_part", "cens_lead_part", ".prediction")

pp_diss <- withr::with_seed(
  214,
  {add_pred_draws_car1(model_in, model_diss, draw_ids = drawids, type = "prediction")}
) %>%
  impute_censored(model_in, "scaled_lead_dissolved", "cens_lead_dissolved", ".prediction")

model_part_resid <- pp_part %>%
  add_resid_draws_car1(
    model_part, scaled_lead_part,
    draw_ids = (1:4000)[-drawids]
  )

model_diss_resid <- pp_diss %>%
  add_resid_draws_car1(
    model_diss, scaled_lead_dissolved,
    draw_ids = (1:4000)[-drawids]
  )

model_part_resid_noar <- pp_part %>%
  add_resid_draws_car1(
    model_part, scaled_lead_part,
    draw_ids = (1:4000)[-drawids], car1 = FALSE
  )

model_diss_resid_noar <- pp_diss %>%
  add_resid_draws_car1(
    model_diss, scaled_lead_dissolved,
    draw_ids = (1:4000)[-drawids], car1 = FALSE
  )

acf_part_resid <- calc_acf(model_part_resid, gr_vars = c(".draw", "location"))
acf_diss_resid <- calc_acf(model_diss_resid, gr_vars = c(".draw", "location"))
acf_part_resid_noar <- calc_acf(model_part_resid_noar, gr_vars = c(".draw", "location"))
acf_diss_resid_noar <- calc_acf(model_diss_resid_noar, gr_vars = c(".draw", "location"))

p1 <- bind_rows(
  "CAR(1)_>0.45 µm" = acf_part_resid,
  "CAR(1)_<0.45 µm" = acf_diss_resid,
  "GAM_>0.45 µm" = acf_part_resid_noar,
  "GAM_<0.45 µm" = acf_diss_resid_noar,
  .id = "model"
) %>%
  separate(model, c("model", "fraction"), sep = "_") %>%
  ggplot(aes(cor, model, fill = fraction)) +
  ggdist::stat_halfeye(slab_alpha = .5, position = position_dodge(width = .1), point_size = 1, show.legend = FALSE) +
  scale_fill_manual(values = palette[c(1,4)]) +
  scale_y_discrete(expand = expansion(add = c(0, .75))) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  labs(y = "Model", x = "Rank correlation", col = NULL, fill = NULL) +
  theme(plot.margin = unit(c(5.5, 6.5, 5.5, 5.5), "points"))

# residuals:

dparams <- list("<0.45 µm" = model_diss, ">0.45 µm" = model_part) %>%
  map_dfr(
    ~ as_draws_df(.x, c("nu", "sigma")) %>%
      posterior::summarize_draws(),
    .id = "type"
  ) %>%
  pivot_wider(id_cols = type, names_from = variable, values_from = median)

resids_diss <- model_diss_resid %>%
  group_by(across(matches(names(model_in)))) %>%
  summarize_preds(retrans = FALSE, pred_var = ".residual")

resids_part <- model_part_resid %>%
  group_by(across(matches(names(model_in)))) %>%
  summarize_preds(retrans = FALSE, pred_var = ".residual")

resids <- bind_rows(
  ">0.45 µm" = resids_part,
  "<0.45 µm" = resids_diss,
  .id = "type"
)

p2 <- resids %>%
  left_join(dparams) %>%
  group_by(type) %>%
  mutate(
    density = dstudent_t(.residual, nu, sigma = sigma),
    density = density / max(density),
  ) %>%
  ungroup() %>%
  ggplot(aes(x = .residual)) +
  facet_wrap(vars(type), scales = "free_x") +
  geom_histogram(
    aes(y = ..ncount..),
    binwidth = .05
  ) +
  geom_vline(xintercept = 0, col = "white", size = .3) +
  geom_line(aes(y = density), col = palette[6], size = .3) +
  labs(x = "Median residual", y = "Normalized\ncounts/density", col = NULL)

p3 <- pk_vals %>%
  pivot_longer(starts_with("pk_vals")) %>%
  mutate(
    name = fct_recode(name, "<0.45 µm" = "pk_vals_diss", ">0.45 µm" = "pk_vals_part")
  ) %>%
  ggplot(aes(date, value, col = name, group = location)) +
  geom_hline(
    yintercept = .5,
    linetype = 3
  ) +
  scale_color_manual(values = palette[c(1, 4)]) +
  geom_line(size = .2, show.legend = FALSE) +
  theme(legend.margin = margin()) +
  labs(
    x = "Year",
    y = "Pareto k\nparameters",
    col = NULL
  )

draws_pars <- list(">0.45 µm" = model_part, "<0.45 µm" = model_diss) %>%
  map_dfr(
    ~ .x %>%
      as_draws_df(variable = "ar[1]") %>%
      as_tibble() %>%
      pivot_longer(-starts_with(".")),
    .id = "model"
  ) %>%
  mutate(
    value = if_else(model == "<0.45 µm" & str_detect(name, "sigma"), exp(value), value)
  )

p4 <- draws_pars %>%
  ggplot(aes(value, fill = model)) +
  ggdist::stat_halfeye(slab_alpha = .5) +
  scale_fill_manual(values = palette[c(1,4)]) +
  scale_y_continuous(expand = expansion(add = c(.2, 0))) +
  theme(
    strip.text = ggtext::element_markdown(),
    axis.title.x = ggtext::element_markdown()
  ) +
  labs(x = "&phi;", y = "Density", fill = NULL, col = NULL)

fig5 <- patchwork::wrap_plots(
  p2, p1, p4, p3,
  nrow = 4,
  heights = c(.8,1.2,.8,.8)
) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_level = "a") &
  theme(
    plot.tag = element_text(face = "bold"),
    legend.margin = margin()
  ) &
  guides(col = guide_legend(override.aes = list(size = .5), ncol = 2))

ggsave("Rmarkdown/figures/figure-5.png", fig5, width = 3.33, height = 6, dpi = 600)

#------------------ figure 6 ------------------

# n.b., dates reflect label positions;
# correct dates are used for the vertical dashed lines

annotations <- tibble::tribble(
  ~ortho_dose, ~x, ~pipe_material, ~labels,
  # "0-0.5", "2018-03-15", "Pb #1", "P introduced",
  # "2.0-0.75", "2018-06-01", "Pb #1", "P decreased",
  "0-0.5", "2020-01-01", "Pb-Cu", "Section<br>removed",
  "1.0", "2020-01-01", "Pb-Cu", "Section<br>removed",
  "2.0-0.75", "2020-01-01", "Pb-Cu", "Section<br>removed",
  "0-0.5", "2017-11-01", "Pb #1", "0",
  "0-0.5", "2017-11-01", "Pb #2", "0",
  "0-0.5", "2017-11-01", "Pb-Cu", "0",
  "0-0.5", "2019-03-01", "Pb #1", "0.5",
  "0-0.5", "2019-03-01", "Pb #2", "0.5",
  "0-0.5", "2019-03-01", "Pb-Cu", "0.5",
  "1.0", "2017-03-01", "Pb #1", "0",
  "1.0", "2017-03-01", "Pb #2", "0",
  "1.0", "2017-03-01", "Pb-Cu", "0",
  "1.0", "2018-04-01", "Pb #1", "1",
  "1.0", "2018-04-01", "Pb #2", "1",
  "1.0", "2018-04-01", "Pb-Cu", "1",
  "2.0-0.75", "2017-03-01", "Pb #1", "0",
  "2.0-0.75", "2017-03-01", "Pb #2", "0",
  "2.0-0.75", "2017-03-01", "Pb-Cu", "0",
  "2.0-0.75", "2018-04-01", "Pb #1", "2",
  "2.0-0.75", "2018-04-01", "Pb #2", "2",
  "2.0-0.75", "2018-04-01", "Pb-Cu", "2",
  "2.0-0.75", "2019-05-01", "Pb #1", "0.75",
  "2.0-0.75", "2019-05-01", "Pb #2", "0.75",
  "2.0-0.75", "2019-05-01", "Pb-Cu", "0.75"
) %>%
  mutate(
    y = if_else(str_detect(labels, "^Section"), .1, Inf),
    ortho_dose = paste0(ortho_dose, " mg P L<sup>-1</sup>"),
    labels = if_else(
      str_detect(labels, "^Section"),
      labels,
      paste0(labels, " mg P L<sup>-1</sup>")
    )
  )

lines <- tibble::tribble(
  ~ortho_dose, ~x, ~pipe_material,
  # section removed:
  "0-0.5", "2020-11-01", "Pb-Cu",
  "1.0", "2020-11-01", "Pb-Cu",
  "2.0-0.75", "2020-11-01", "Pb-Cu",
  # P introduced (0.5 mg/L):
  "0-0.5", "2019-01-29", "Pb #1",
  "0-0.5", "2019-01-29", "Pb #2",
  "0-0.5", "2019-01-29", "Pb-Cu",
  # P introduced (other 2 doses):
  "1.0", "2018-03-13", "Pb #1",
  "1.0", "2018-03-13", "Pb #2",
  "1.0", "2018-03-13", "Pb-Cu",
  "2.0-0.75", "2018-03-13", "Pb #1",
  "2.0-0.75", "2018-03-13", "Pb #2",
  "2.0-0.75", "2018-03-13", "Pb-Cu",
  # P decreased:
  "2.0-0.75", "2019-04-16", "Pb #1",
  "2.0-0.75", "2019-04-16", "Pb #2",
  "2.0-0.75", "2019-04-16", "Pb-Cu"
) %>%
  mutate(
    x = as.Date(x),
    ortho_dose = paste0(ortho_dose, " mg P L<sup>-1</sup>")
  )

preds_diss <- add_pred_draws_car1(model_in, model_diss)
preds_part <- add_pred_draws_car1(model_in, model_part)

preds_diss_sum <- preds_diss %>%
  select(.epred) %>%
  summarize_preds(y_var = model_in$lead_dissolved, recensor = TRUE)

preds_part_sum <- preds_part %>%
  select(.epred) %>%
  summarize_preds(y_var = model_in$lead_part, recensor = TRUE)

plot_fig6 <- function(preds_diss_sum, preds_part_sum) {

  bind_rows(
    "<0.45 µm" = preds_diss_sum %>%
      pivot_longer(c(lead_dissolved, .epred_retrans)),
    ">0.45 µm" = preds_part_sum %>%
      pivot_longer(c(lead_part, .epred_retrans)),
    .id = "type"
  ) %>%
    unite(name, c(name, type)) %>%
    mutate(
      name = fct_relevel(name,
        c(".epred_retrans_>0.45 µm", ".epred_retrans_<0.45 µm"),
        after = Inf
      ),
      ortho_dose = paste0(ortho_dose, " mg P L<sup>-1</sup>")
    ) %>%
    assertr::verify(value > 0) %>%
    assertr::verify(.lower_retrans > 0) %>%
    assertr::verify(.upper_retrans > 0) %>%
    ggplot(aes(date, col = name, fill = name)) +
    facet_grid(rows = vars(pipe_material), cols = vars(ortho_dose)) +
    geom_rect(
      data = tibble(
        xmin = as.Date(c("2019-10-01", "2020-02-25")),
        xmax = as.Date(c("2019-10-07", "2020-03-10")),
        ymin = 0, ymax = Inf
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, col = "no P", linetype = NA),
      inherit.aes = FALSE, alpha = .75, fill = palette[4],
      show.legend = c("col" = FALSE, "fill" = FALSE)
    ) +
    # geom_rect(
    #   data = function(x) {
    #     x %>%
    #       group_by(pipe_material, ortho_dose) %>%
    #       summarize(xmin = min(date)) %>%
    #       mutate(
    #         xmax = as.Date("2018-03-13"),
    #         ymin = 0, ymax = Inf
    #       )
    #   },
    #   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #   inherit.aes = FALSE, alpha = .2
    # ) +
    # geom_segment(
    #   data = lines[1:3,] %>%
    #     mutate(y = 0, yend = Inf),
    #   aes(x = x, xend = x, y = y, yend = yend, col = "zPb removed"),
    #   inherit.aes = FALSE, #linetype = 3
    # ) +
    ggtext::geom_richtext(
      data = annotations,
      aes(x = as.Date(x), y = y, label = labels),
      inherit.aes = FALSE, vjust = "inward", hjust = 0,
      label.padding = unit(0.2, "lines"), label.size = 0,
      show.legend = FALSE, size = 1.75, alpha = .75,
      label.r = unit(0, "cm")
    ) +
    geom_vline(
      # data = lines[-3:-1,],
      data = lines,
      aes(xintercept = x),
      linetype = 3, show.legend = FALSE
    ) +
    geom_ribbon(
      data = function(x) {
        x %>%
          filter(str_detect(name, "^.epred"))
      },
      aes(ymin = .lower_retrans, ymax = .upper_retrans),
      alpha = .5, col = NA,
      show.legend = FALSE
    ) +
    geom_line(aes(y = value, col = name)) +
    theme(strip.text = ggtext::element_markdown()) +
    scale_color_manual(
      values = c(palette[c(6,1)], "black", "grey", palette[4:3]),
      labels = c(
        "GAM (<0.45 µm)", "GAM (>0.45 µm)", "<0.45 µm",
        ">0.45 µm", "P interruption", "Section removed"
      )
    ) +
    scale_fill_manual(values = c(palette[c(6,1)], "black", "grey")) +
    geom_rug(
      data = function(x) x %>% filter(value > 1.9),
      sides = "t",
      aes(y = value)
    ) +
    scale_y_log10(
      limits = c(1.9e-4, 1.9),
      labels = function(breaks) 1e3 * breaks,
    ) +
    labs(
      x = NULL,
      y = expression("[Pb] (µg L"^-1 * ")"),
      col = NULL,
      fill = NULL
    )
}

fig6 <- plot_fig6(preds_diss_sum, preds_part_sum)

ggsave("Rmarkdown/figures/figure-6.png", fig6, width = 7, height = 4, dpi = 600)

#------------------ figure 7 ------------------

preds_d <- bind_rows(
  "<0.45 µm" = preds_diss %>%
    ungroup() %>%
    mutate(.epred = retrans(.epred, model_in$lead_dissolved)),
  ">0.45 µm" = preds_part %>%
    ungroup() %>%
    mutate(.epred = retrans(.epred, model_in$lead_part)),
  .id = "type"
) %>%
  pivot_wider(
    id_cols = c(.draw, type, pipe_material, date),
    names_from = ortho_dose,
    values_from = .epred
  ) %>%
  transmute(
    .draw,  type, pipe_material, date,
    `0-0.5 mg P L<sup>-1</sup> /<br> 1.0 mg P L<sup>-1</sup>` = `0-0.5` / `1.0`,
    `0-0.5 mg P L<sup>-1</sup> /<br> 2.0-0.75 mg P L<sup>-1</sup>` = `0-0.5` / `2.0-0.75`,
    `1.0 mg P L<sup>-1</sup> /<br> 2.0-0.75 mg P L<sup>-1</sup>` = `1.0` / `2.0-0.75`
  ) %>%
  pivot_longer(matches("^\\d"), names_to = "difference", values_to = ".epred") %>%
  filter(!is.na(.epred))

preds_d_summ <- preds_d %>%
  group_by(pipe_material, type, date, difference) %>%
  summarize_preds(retrans = FALSE)

preds_d_spag <- preds_d %>%
  filter(.draw %in% withr::with_seed(12450, {sample(1:4000, 200)}))

lines_d <- preds_d %>%
  distinct(difference) %>%
  crossing(filter(lines, lubridate::year(x) == 2019)) %>%
  filter(str_detect(difference, ortho_dose)) %>%
  distinct(difference, x, pipe_material) %>%
  mutate(
    label = if_else(
      x == "2019-01-29",
      "P introduced:<br>0-0.5 mg P L<sup>-1</sup>",
      "P decreased:<br>2-0.75 mg P L<sup>-1</sup>"
    ),
    xlab = if_else(
      x == "2019-01-29",
      as.Date("2018-02-01"),
      as.Date("2018-11-01")
    )
  )

ratio_plot <- function(x, ...) {

  x %>%
    filter(...) %>%
    ggplot(aes(date, .epred)) +
    facet_grid(
      rows = vars(pipe_material),
      cols = vars(type)
    ) +
    geom_hline(yintercept = 1) +
    geom_vline(
      data = lines_d %>%
        filter(...),
      aes(xintercept = x), linetype = 3
    ) +
    ggtext::geom_richtext(
      data = lines_d %>%
        filter(...) %>%
        slice(1),
      aes(x = xlab, y = 0, label = label),
      inherit.aes = FALSE, vjust = "inward", hjust = 0,
      label.padding = unit(0.2, "lines"), label.size = 0,
      show.legend = FALSE, size = 2, alpha = .75,
      label.r = unit(0, "cm")
    ) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), col = NA, alpha = .5) +
    geom_line(
      data = function(x) {
        x %>%
          distinct(type) %>%
          right_join(preds_d_spag) %>%
          filter(...)
      },
      aes(group = .draw), size = .05, col = palette[2]
    ) +
    geom_line() +
    # generates missing values by design:
    geom_line(
      data = function(x) {
        x %>%
          mutate(
            sig = sign(log(.lower)) == sign(log(.upper)),
            .epred = if_else(sig, .epred, NA_real_)
          )
      },
      col = palette[6]
    ) +
    labs(x = NULL, col = NULL, fill = NULL) +
    scale_y_log10() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
}

p1 <- preds_d_summ %>%
  ratio_plot(
    pipe_material %in% c("Pb #1", "Pb #2"),
    difference == "0-0.5 mg P L<sup>-1</sup> /<br> 1.0 mg P L<sup>-1</sup>"
  ) +
  labs(y = expression(frac("0-0.5 mg P L"^-1, "1 mg P L"^-1)))

p2 <- preds_d_summ %>%
  ratio_plot(
    pipe_material %in% c("Pb #1", "Pb #2"),
    difference == "1.0 mg P L<sup>-1</sup> /<br> 2.0-0.75 mg P L<sup>-1</sup>"
  ) +
  labs(y = expression(frac("1 mg P L"^-1, "2-0.75 mg P L"^-1)))

fig7 <- patchwork::wrap_plots(p1, p2, nrow = 2) +
  patchwork::plot_annotation(tag_level = "a") +
  patchwork::plot_layout(guides = "collect") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("Rmarkdown/figures/figure-7.png", fig7, width = 3.33, height = 5, dpi = 600)

# for text:

preds_d_summ %>%
  filter(pipe_material %in% c("Pb #1", "Pb #2"), type == "<0.45 µm") %>%
  slice_max(.epred) %>%
  readr::write_csv("data-clean/figure-7-summary.csv")

# graphical abstract

toc_art <- preds_d_summ %>%
  ratio_plot(
    pipe_material %in% c("Pb #2"),
    difference == "1.0 mg P L<sup>-1</sup> /<br> 2.0-0.75 mg P L<sup>-1</sup>"
  ) +
  theme(strip.text.y = element_blank()) +
  labs(y = expression(frac("1 mg P L"^-1, "2-0.75 mg P L"^-1)))

ggsave("Rmarkdown/figures/figure-toc.png", toc_art, width = 3.33, height = 1.8, dpi = 600)

#------------------ figure 8 ------------------

epsilon <- .001

grid_1 <- with(
  model_in,
  tibble(
    date_numeric = seq(
      min(date_numeric),
      max(date_numeric),
      by = 2 * epsilon
    )
  )
)

grid_2 <- grid_1 %>%
  mutate(date_numeric = date_numeric + epsilon)

grid_avg <- (grid_1 + grid_2) / 2

slopes <- list(
  "<0.45 µm" = model_diss,
  ">0.45 µm" = model_part
) %>%
  map_dfr(
    ~ local_slope(grid_avg, .x, "date_numeric", smooth = "s(date_numeric)"),
    .id = "fraction"
  ) %>%
  rename(Smooth = smooth, "Local slope" = slope) %>%
  pivot_longer(c(Smooth, `Local slope`))

slopes_sum <- slopes %>%
  select(-.draw) %>%
  group_by(fraction, date_numeric, name) %>%
  summarize_preds(pred_var = "value", retrans = FALSE)

slopes_sig <- slopes_sum %>%
  pivot_wider(
    id_cols = c(fraction, date_numeric),
    names_from = c(name),
    values_from = c(value, .lower, .upper)
  ) %>%
  arrange(fraction, date_numeric) %>%
  mutate(
    value = if_else(
      sign(`.lower_Local slope`) == sign(`.upper_Local slope`),
      value_Smooth,
      NA_real_
    ),
    name = "Smooth"
  )

slopes_spag <- slopes %>%
  filter(.draw %in% withr::with_seed(8765309, {sample(max(slopes$.draw), 200)}))

fig8a <- slopes_sum %>%
  ggplot(aes(date_numeric + 2017)) +
  facet_grid(
    rows = vars(name = fct_relevel(name, "Smooth", after = 0L)),
    cols = vars(fraction),
    scales = "free_y",
    labeller = labeller(name = as_labeller(c("Local slope" = "Local\nslope", "Smooth" = "Trend")))
  ) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .5) +
  geom_line(
    data = slopes_spag,
    aes(y = value, group = .draw),
    size = .1, alpha = .5,
    col = palette[2]
  ) +
  geom_line(aes(y = value), show.legend = FALSE) +
  # generates missing values by design:
  geom_line(
    data = slopes_sig,
    aes(y = value), show.legend = FALSE,
    col = palette[6]
  ) +
  # generates a missing value by design:
  geom_hline(
    data = function(x) {
      x %>%
        distinct(name) %>%
        mutate(int = c(NA, 0))
    },
    aes(yintercept = int),
    col = "grey"
  ) +
  labs(x = NULL, y = "Effect") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

fig8b_in <- pdat %>%
  filter(
    param == "Iron",
    str_detect(location, "LP\\dA$")
  ) %>%
  mutate(
    yday = lubridate::yday(date),
    yr = lubridate:::year(date),
    param = glue::glue("{str_to_sentence(param)} ({units})") %>%
      str_replace("\\/L", " L<sup>-1</sup>")
  ) %>%
  bind_rows(
    smooth_terms_diss[["mu: s(date_yday,bs=\"cc\")"]] %>%
      mutate(
        param = "GAM (seasonal component)",
        yday = 365 * date_yday
      )
  ) %>%
  mutate(
    # this is an approximation for the leap year, 2020
    m_d = as.Date(glue::glue("2022-{pmax(yday, 1)}"), "%Y-%j")
  )

fig8b <- fig8b_in %>%
  ggplot(aes(m_d)) +
  facet_wrap(
    vars(param = fct_relevel(param, "GAM (seasonal component)", after = 0L)),
    scales = "free_y", ncol = 1
  ) +
  scale_color_gradientn(colours = palette[-3]) +
  scale_x_date(date_labels = "%b") +
  ggh4x::facetted_pos_scales(
    # omits values!!
    y = list(
      param == "Iron (mg L<sup>-1</sup>)"
      ~ scale_y_continuous(limits = c(0, .5))
    )
  ) +
  # adds missing values by design:
  geom_point(
    data = function(x) filter(x, !(replace_na(bdl, FALSE))),
    aes(y = value), alpha = .7, shape = 16, col = "grey"
  ) +
  # adds missing values by design:
  geom_smooth(
    data = function(x) filter(x, !(replace_na(bdl, FALSE))),
    aes(y = value), formula = y ~ s(x, bs = "cc"),
    method = "gam", se = FALSE, col = "black"
  ) +
  # show omitted values as segments:
  # (adds missing values by design):
  geom_segment(
    data = function(x) {
      x %>%
        filter(value > .5 | is.na(value)) %>%
        mutate(start = if_else(is.na(value), NA_real_, .5))
    },
    aes(xend = m_d, y = start, yend = Inf),
    show.legend = FALSE
  ) +
  # adds missing values by design:
  geom_segment(
    data = function(x) {
      x %>%
        mutate(end = if_else(bdl, value, NA_real_))
    },
    aes(xend = m_d, y = -Inf, yend = end),
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = function(x) {
      x %>%
        filter(!str_detect(param, "Iron"))
    },
    aes(ymin = lower__, ymax = upper__), alpha = .5
  ) +
  geom_line(
    data = smooth_terms_spag %>%
      filter(smooth_term == "s(date_yday,bs=\"cc\")") %>%
      unnest(smooth) %>%
      mutate(
        yday = 365 * rep(model_in$date_yday, n_smooths),
        param = "GAM (seasonal component)",
        m_d = as.Date(glue::glue("2022-{pmax(yday, 1)}"), "%Y-%j")
      ) %>%
      distinct(smooth_term, .draw, yday, value, param, m_d),
    aes(y = value, group = .draw),
    col = palette[2], size = .05
  ) +
  geom_line(aes(y = estimate__), show.legend = FALSE) +
  theme(strip.text = ggtext::element_markdown()) +
  labs(
    x = NULL,
    y = NULL,
    col = NULL
  )

fig8c <- pdat %>%
  filter(
    param == "Water Temperature",
    str_detect(location, "LP\\d[123]$")
  ) %>%
  group_by(
    date,
    pipe_material = fct_relevel(pipe_material, "Pb #1", "Pb #2", after = 0L),
    param = "Water temperature (&deg;C)"
  ) %>%
  summarize(value = median(value)) %>%
  ungroup() %>%
  ggplot(aes(date, value, col = pipe_material)) +
  facet_wrap(vars(param)) +
  geom_line(data = function(x) filter(x, date < "2020-10-01")) +
  geom_line(data = function(x) filter(x, date >= "2020-10-01")) +
  scale_color_manual(values = palette[c(6,4,1)]) +
  theme(
    legend.margin = margin(),
    strip.text = ggtext::element_markdown()
  ) +s
  labs(
    x = NULL,
    y = NULL,
    col = NULL
  )

fig8 <- patchwork::wrap_plots(
  fig8a, fig8b, fig8c,
  nrow = 3, heights = c(.8, 1.2, .5)
) +
  patchwork::plot_annotation(tag_level = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("Rmarkdown/figures/figure-8.png", fig8, width = 3.33, height = 6, dpi = 600)

#------------------ figure 9 ------------------

annotations_short <- tibble::tribble(
            ~x, ~pipe_material,                                      ~labels, ~y, ~ortho_dose,
  "2018-01-01",        "Pb #1", "P introduced:<br>0-0.5 mg P L<sup>-1</sup>", Inf,    "0-0.5",
  "2019-03-01",        "Pb #2", "P decreased:<br>2-0.75 mg P L<sup>-1</sup>", Inf, "2.0-0.75"
  )

lines_short <- lines %>%
  distinct(ortho_dose, x) %>%
  filter(x < "2020-01-01") %>%
  mutate(ortho_dose = str_remove(ortho_dose, " mg P L<sup>-1</sup>"))

preds_fdiss <- preds_diss %>%
  ungroup() %>%
  mutate(
    .epred_retrans = retrans(.epred, model_in$lead_dissolved),
    lead_part_retrans = retrans(preds_part$.epred, model_in$lead_part),
    lead_total = .epred_retrans + lead_part_retrans,
    fdiss = lead_part_retrans / lead_total
  )

preds_fdiss_summ <- preds_fdiss %>%
  group_by(across(group_vars(preds_diss))) %>%
  summarize_preds(retrans = FALSE, pred_var = "fdiss")

fig9 <- preds_fdiss_summ %>%
  ggplot(aes(date, fdiss, col = ortho_dose, fill = ortho_dose)) +
  facet_wrap(vars(pipe_material), ncol = 1) +
  geom_rect(
    data = function(x) {
      x %>%
        group_by(pipe_material, ortho_dose) %>%
        summarize(xmin = min(date)) %>%
        mutate(
          xmax = as.Date("2018-03-13"),
          ymin = 0, ymax = Inf
        )
    },
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, alpha = .2
  ) +
  geom_vline(
    data = lines_short,
    aes(xintercept = x, col = ortho_dose),
    linetype = 3, show.legend = FALSE
  ) +
  ggtext::geom_richtext(
    data = annotations_short,
    aes(x = as.Date(x), y = y, label = labels, col = ortho_dose),
    inherit.aes = FALSE, vjust = "inward", hjust = 0,
    label.padding = unit(0.2, "lines"), label.size = 0,
    show.legend = FALSE, size = 2.5, fill = alpha("white", 0.75),
    label.r = unit(0, "cm")
  ) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .5, col = NA) +
  geom_line() +
  scale_fill_manual(values = palette[c(6,4,1)]) +
  scale_color_manual(values = palette[c(6,4,1)]) +
  labs(
    x = NULL,
    col = expression("mg P L"^-1),
    fill = expression("mg P L"^-1),
    y = expression(frac("[Pb]"[">0.45 µm"], "[Pb]"[total]))
  ) +
  theme(
    legend.margin = margin(r = 25)
  )

ggsave("Rmarkdown/figures/figure-9.png", fig9, width = 3.33, height = 4.5, dpi = 600)

#------------------ figure s2 ------------------

params <- c(
  "ar[1]" = "&phi;",
  "sds_sdate_yday_1" = "&sigma;<sub>seasonal spline</sub>",
  "nu" = "&nu;"
)

compare_priors <- crossing(
  name = names(params),
  model = c("brms default", "This paper")
) %>%
  mutate(
    prior = c(
      "flat", "N(0.5,0.25)", # ar[1]
      "&Gamma;(2,0.1)", "&Gamma;(2,0.5)", # nu
      "T(3, 0, 2.5)", "N(0,0.5)" # sds_sdate_yday_1
    )
  )

fig_s2_in <- list(
  "brms default_<0.45 µm" = as_draws_df(model_noprior_diss, names(params)),
  "This paper_<0.45 µm" = as_draws_df(model_diss, names(params)),
  "brms default_>0.45 µm" = as_draws_df(model_noprior_part, names(params)[-2]),
  "This paper_>0.45 µm" = as_draws_df(model_part, names(params)[-2])
) %>%
  map_dfr(as_tibble, .id = "model") %>%
  separate(model, c("model", "fraction"), sep = "_") %>%
  pivot_longer(-c(model, fraction, starts_with("."))) %>%
  # add priors:
  left_join(compare_priors) %>%
  filter(!is.na(value)) # removes seasonal term from part model

fig_s2 <- fig_s2_in %>%
  ggplot(aes(value, model)) +
  facet_grid(
    cols = vars(name),
    rows = vars(fraction),
    scales = "free",
    labeller = labeller(.cols = as_labeller(params))
  ) +
  ggtext::geom_richtext(
    data = function(x) {
      x %>%
        group_by(model, name, fraction, prior) %>%
        summarize(value = median(value))
    },
    aes(label = prior), nudge_y = -.2,
    label.size = 0, size = 3
  ) +
  ggdist::stat_halfeye(normalize = "panels", scale = .6) +
  theme(strip.text = ggtext::element_markdown()) +
  scale_x_continuous(expand = expansion(c(.2, .2))) +
  scale_y_discrete(expand = expansion(c(.5, 0))) +
  labs(x = NULL, y = NULL)

ggsave("Rmarkdown/figures/figure-s2.png", fig_s2, width = 7.5, height = 4, dpi = 600)

#------------------ figure s3 and figure s4 ------------------

plot_fun <- function(x, yvar = reorder(param, median), ...) {
  x %>%
    mutate(low_n = ess_bulk < 400 | ess_tail < 400) %>%
    arrange(rhat) %>%
    ggplot(aes(median, {{ yvar }})) +
    geom_vline(xintercept = 0, col = "grey") +
    geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0, size = .1, col = "grey35") +
    geom_point(aes(col = rhat, shape = low_n, size = rhat), alpha = .7) +
    scale_y_discrete(
      labels = function(breaks) if_else(str_detect(breaks, "^Y|^b?s_|^r_"), "", breaks)
    ) +
    guides(size = "none", shape = "none") +
    labs(
      x = "Median",
      y = NULL,
      shape = expression("N"[eff]),
      col = expression(italic(widehat("R")) ~ " ")
    ) +
    theme(panel.grid = element_blank())
}

panel_fun <- function(model_summ) {
  model_summ %>%
    mutate(
      location = str_extract(variable, "LP\\d{2}"),
      type = case_when(
        # GAMS:
        str_detect(variable, "^s_") ~ "GAM coefficients (penalized)",
        str_detect(variable, "^sds_.+location.+") ~ "Standard deviations of\nGAM coefficients",
        str_detect(variable, "zs_") ~ "Standardized spline coefficients",
        str_detect(variable, "^r_location.+Intercept]") ~ "Random intercepts",
        str_detect(variable, "^[rzLC](or)?_\\d\\[") ~ "internal",
        TRUE ~ "Other"
      ),
      param = case_when(
        # GAM
        variable == "ar[1]" ~ "Autoregressive parameter",
        variable == "nu" ~ "Degrees of freedom (t-distribution)",
        variable == "b_Intercept" ~ "Intercept",
        variable == "sigma" ~ "Residual SD",
        variable == "sds_sdate_yday_1" ~ "SD of seasonal trend",
        variable == "sds_sdate_numeric_1" ~ "SD of multiyear trend",
        variable == "bs_sdate_numeric_1" ~ "Unpenalized GAM coefficient",
        variable == "sd_location__Intercept" ~ "SD of random intercepts",
        variable == "cor_location__Intercept__typelead_dissolved" ~ "Random slope/intercept correlation",
        TRUE ~ variable
      )
    ) %>%
    left_join(locations) %>%
    mutate(ortho_dose = paste0(ortho_dose, " mg P L<sup>-1</sup>")) %>%
    group_split(type, .keep = TRUE)
}

patchwork_fun <- function(panels, model_summ) {
  p1 <- plot_fun(panels[[1]]) + # penalized GAM coefficients
    theme(axis.ticks.y = element_blank())
  p2 <- plot_fun(panels[[3]])
  p3 <- plot_fun(panels[[4]], ortho_dose) + # rand intercepts
    facet_wrap(vars(pipe_material)) +
    scale_x_continuous(n.breaks = 4) +
    theme(axis.text.y = ggtext::element_markdown())
  p4 <- plot_fun(panels[[5]], ortho_dose) + # gam sds
    facet_wrap(vars(pipe_material)) +
    scale_x_continuous(n.breaks = 4) +
    theme(axis.text.y = ggtext::element_markdown())

  patchwork::wrap_plots(
    p1, p2, p3, p4,
    design = "AB\nAC\nAD"
  ) +
    patchwork::plot_annotation(tag_levels = "a") +
    patchwork::plot_layout(guides = "collect") &
    scale_color_gradientn(
      colors = palette[-3],
      limits = range(model_summ$rhat, na.rm = TRUE),
      n.breaks = 3
      # breaks = c(1, 1.004, 1.008)
    ) &
    scale_size(limits = range(model_summ$rhat)) &
    scale_shape_manual(
      values = c(16, 17),
      labels = c("≥ 400", "< 400"),
      limits = c(FALSE, TRUE)
    ) &
    theme(plot.tag = element_text(face = "bold")) &
    guides(shape = guide_legend(override.aes = list(size = 2)))
}

model_diss_summ <- model_diss %>%
  as_draws() %>%
  posterior::summarise_draws() %>%
  filter(!variable %in% c("lp__", "Intercept", "lprior"))

model_part_summ <- model_part %>%
  as_draws() %>%
  posterior::summarise_draws() %>%
  filter(!variable %in% c("lp__", "Intercept", "lprior"))

fig_s3 <- panel_fun(model_diss_summ) %>%
  patchwork_fun(model_diss_summ)

fig_s4 <- panel_fun(model_part_summ) %>%
  patchwork_fun(model_part_summ)

ggsave("Rmarkdown/figures/figure-s3.png", fig_s3, width = 7.5, height = 7, dpi = 600)
ggsave("Rmarkdown/figures/figure-s4.png", fig_s4, width = 7.5, height = 7, dpi = 600)

#------------------ figure s5 ------------------

smooth_terms_part_sim <- conditional_smooths(model_sim_part)
smooth_terms_diss_sim <- conditional_smooths(model_sim_diss)

smooth_terms_sim_spag <- calc_spag_smooths(model_sim_diss, model_sim_part)

figs5 <- plot_fig4(smooth_terms_diss_sim, smooth_terms_part_sim, smooth_terms_sim_spag)

ggsave("Rmarkdown/figures/figure-s5.png", figs5, width = 7.5, height = 4.75, dpi = 600)

#------------------ figure s6 and figure s7 ------------------

model_sim_diss_summ <- model_sim_diss %>%
  as_draws() %>%
  posterior::summarise_draws() %>%
  filter(!variable %in% c("lp__", "Intercept"))

model_sim_part_summ <- model_sim_part %>%
  as_draws() %>%
  posterior::summarise_draws() %>%
  filter(!variable %in% c("lp__", "Intercept"))

fig_s6 <- panel_fun(model_sim_diss_summ) %>%
  patchwork_fun(model_sim_diss_summ)

fig_s7 <- panel_fun(model_sim_part_summ) %>%
  patchwork_fun(model_sim_part_summ)

ggsave("Rmarkdown/figures/figure-s6.png", fig_s6, width = 7.5, height = 7, dpi = 600)
ggsave("Rmarkdown/figures/figure-s7.png", fig_s7, width = 7.5, height = 7, dpi = 600)

#------------------ figure s8 ------------------

fig_s8 <- resids %>%
  mutate(ortho_dose = paste0(ortho_dose, " mg P L<sup>-1</sup>")) %>%
  ggplot(aes(date, .residual, col = type)) +
  facet_grid(
    cols = vars(ortho_dose),
    rows = vars(pipe_material)
  ) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0, size = .1) +
  geom_point(size = .2) +
  scale_color_manual(values = palette[c(1,6)]) +
  theme(
    strip.text.x = ggtext::element_markdown(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(x = NULL, y = "Simulated residual (transformed scale)", col = NULL)

ggsave("Rmarkdown/figures/figure-s8.png", fig_s8, width = 7, height = 4.5, dpi = 600)

#------------------ session info ------------------

writeLines(capture.output(devtools::session_info()), "session-info.txt")
