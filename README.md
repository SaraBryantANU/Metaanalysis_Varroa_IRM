
#Meta analysis of acaricide resistance in *Varroa destructor*"

##load libraries 
```{r setup =FALSE, warning=FALSE}
library(tidyverse)
library(readxl)
library(brms)
library(tidybayes)
library(readr)
library(dplyr)
library(loo)
library(patchwork)
library(metafor)
library(Matrix)
library(grid)
```

##Read data|
```{r}
dat <- read_xlsx("final data metaanalysis SE corrected.xlsx") %>%
  mutate( # change spelling so everything uses "Tau-fluvalinate"
    acaricide = recode(acaricide, "TauFluvalinate" = "Tau-fluvalinate"),
         rot_mix = rotation | mixture,
    n   = as.integer(n),
    k   = round(prop_dead * n)     # successes = dead mites
  )

dat$acaricide <- factor(dat$acaricide, levels = c("None", "Amitraz", "Coumaphos", "Flumethrin", "Tau-fluvalinate")) 

```

##Model 1 (no concentartion term)
```{r}
priors <- c(
  prior(normal(0, 5), class = Intercept),
  prior(exponential(1), class = sd)
)

fit <- brm(
  k | trials(n) ~ 1 + (1 | study_id) +  acaricide + rotation + mixture + years_stopped + years_treated + rotation * mixture,
  data   = dat,
  family = binomial(link = "logit"),
  prior  = priors,
  chains = 4, cores = 4, seed = 202502,
  save_pars = save_pars(all = TRUE) # necessary to save all parameters for LOO
)

summary(fit)
```
##Model 2 (concentration term)
```{r}
dat <- dat %>% mutate(scaled_log_conc = scale(log(conc + 1))) 

priors2 <- c(
  prior(normal(0, 5), class = Intercept),
  prior(exponential(1), class = sd),
  prior(normal(0, 2), class = b) # better prior for scaled concentration
)

fit2 <- brm(
  k | trials(n) ~ 1 + (1 | study_id) + acaricide*scaled_log_conc + rotation + mixture + years_stopped + years_treated + rotation * mixture,
  data = dat,
  family = binomial(link = "logit"),
  prior = priors2,
  chains = 4, cores = 4, seed = 202502,
  iter = 2000,
  save_pars = save_pars(all = TRUE)
)

summary(fit2)

```
##Compare model fit
```{r}
loo(fit2, fit)
```
```{r}
pp_check(fit2)
```
##Compute heterogeneity via ICC
```{r}
vc <- VarCorr(fit)
(ICC <- vc$study_id$sd[1]/(vc$study_id$sd[1] + pi**2/3))
```
##Publication bias 
###funnel plot
```{r}
funnel_data <- posterior_summary(fit2) %>%
  as_tibble(rownames = "parameter") %>%
  filter(str_detect(parameter, "^r_study_id\\[.*\\]$")) %>%
  rename(effect_size = Estimate, std_error = Est.Error) %>%
  mutate(study_id = str_extract(parameter, "\\d+"))

ggplot(funnel_data, aes(x = effect_size, y = std_error)) +
  geom_segment(aes(x = Q2.5, xend = Q97.5, y = std_error, yend = std_error),
               color = "skyblue", alpha = 0.7) +
  geom_point(shape = 21, size = 3, fill = "lightblue", alpha = 0.8) +   
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_reverse() +
  labs(
    title = "Funnel Plot for Publication Bias",
    subtitle = "Showing 95% credible intervals for each study-specific effect",
    x = "Study-Specific Effect (Log-Odds)",
    y = "Posterior Standard Error"
  ) +
  theme_minimal()

```
### funnel plot with study specific labels
```{r}

# 1) Build a clean lookup from your raw data
name_lu <- dat %>%
  mutate(
    # standardize to a single column `study_name` no matter what the source is called
    study_name = dplyr::coalesce(Study_name, Study_name),
    study_id   = as.character(study_id)
  ) %>%
  distinct(study_id, study_name) %>%
  group_by(study_id) %>%
  slice(1) %>%
  ungroup()

# 2) Posterior-based per-study effects + join names
funnel_data <- posterior_summary(fit2) %>%
  as_tibble(rownames = "parameter") %>%
  filter(str_detect(parameter, "^r_study_id\\[.*\\]$")) %>%
  rename(effect_size = Estimate, std_error = Est.Error) %>%
  mutate(study_id = str_extract(parameter, "\\d+")) %>%
  left_join(name_lu, by = "study_id") %>%
  mutate(
    study_name = coalesce(study_name, paste0("Study ", study_id))
  )

# Build stable breaks/labels from your data
lab_df <- funnel_data |>
  distinct(study_name) |>
  mutate(label_wrapped = str_wrap(study_name, width = 22)) |>
  arrange(study_name)

ggplot(funnel_data, aes(x = effect_size, y = std_error)) +
  geom_segment(aes(x = Q2.5, xend = Q97.5, y = std_error, yend = std_error),
               color = "skyblue", alpha = 0.7) +
  geom_point(aes(fill = study_name), shape = 21, size = 3, color = "grey30") +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_reverse() +
  scale_fill_viridis_d(
    name   = "Study",
    breaks = lab_df$study_name,
    labels = lab_df$label_wrapped,
    option = "turbo", end = 0.9
  ) +
  labs(
    title = "Funnel Plot for Publication Bias",
    subtitle = "Showing 95% credible intervals for each study-specific effect",
    x = "Study-Specific Effect (Log-Odds)",
    y = "Posterior Standard Error"
  ) +
  theme_minimal() +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 9),
    legend.key.size  = unit(4, "mm"),
    legend.spacing.y = unit(1, "mm")
  ) +
  guides(
    fill = guide_legend(
      nrow = 3,         
      byrow = TRUE
    )
  )
```
### standard funnel plot
```{r}
# data

funnel_data <- posterior_summary(fit2) %>%
  as_tibble(rownames = "parameter") %>%
  filter(str_detect(parameter, "^r_study_id\\[.*\\]$")) %>%
  rename(effect_size = Estimate, std_error = Est.Error) %>%
  mutate(study_id = str_extract(parameter, "\\d+"))

#create metaanalysis object
meta_analysis <- rma(yi = effect_size, sei = std_error, data = funnel_data, method = "REML" )

#Funnel plot

funnel(meta_analysis, main = "funnel plot", xlab = "Study-Specific Effect (Log-Odds)", ylab = "Standard Error")

#eggers test for funnel plot symmetry
regtest(meta_analysis)

# Begg’s rank correlation test (complementary, different assumptions)
ranktest(meta_analysis)

```
##estimate baseline (mite mortality without acarcide or IRM)

```{r}
# Baseline probability (probability of mite mortality (no acarcide, no rotation, no mixture, 0 years treated/stopped))

base_row <- tibble(
  years_treated = 0, years_stopped = 0,
  rotation = 0, mixture = 0,
  acaricide = factor("None", levels = levels(dat$acaricide)),
  n = 1L
)

p0 <- posterior_linpred(
  fit, newdata = base_row, re_formula = NA, transform = TRUE
) |> as.numeric() |> mean()

p0

```

## Probability plots: predicted mortality vs years treated / years stopped
```{r}

### 1) Predictions grid with other covariates held at baseline: no rotation, no mixture, acaricide = None.
newdat_treat <- tibble(
  years_treated = seq(min(dat$years_treated, na.rm = TRUE),
                      max(dat$years_treated, na.rm = TRUE), by = 1),
  years_stopped = 0,
  rotation = 0, mixture = 0,
  acaricide = factor("None", levels = levels(dat$acaricide)),
  n = 1L
)

newdat_stop <- tibble(
  years_treated = 0,
  years_stopped = seq(min(dat$years_stopped, na.rm = TRUE),
                      max(dat$years_stopped, na.rm = TRUE), by = 1),
  rotation = 0, mixture = 0,
  acaricide = factor("None", levels = levels(dat$acaricide)),
  n = 1L
)

### 2) Posterior predictions on the response scale (probabilities), summarised to mean & 95% CrI
pred_treat <- add_fitted_draws(
  fit, newdata = newdat_treat, re_formula = NA, scale = "response"
) %>%
  group_by(years_treated) %>%
  summarise(
    p = mean(.value),
    l = quantile(.value, 0.025),
    u = quantile(.value, 0.975),
    .groups = "drop"
  )

pred_stop <- add_fitted_draws(
  fit, newdata = newdat_stop, re_formula = NA, scale = "response"
) %>%
  group_by(years_stopped) %>%
  summarise(
    p = mean(.value),
    l = quantile(.value, 0.025),
    u = quantile(.value, 0.975),
    .groups = "drop"
  )

p_main1 <- ggplot(pred_treat, aes(x = years_treated, y = p)) +
  geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(breaks = pretty(pred_treat$years_treated)) +
  labs(
    x = "Years treated (cumulative)",
    y = "Predicted mortality (%)"
  ) +
  theme_minimal(base_size = 12)

p_main2 <- ggplot(pred_stop, aes(x = years_stopped, y = p)) +
  geom_ribbon(aes(ymin = l, ymax = u), alpha = 0.20) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(breaks = pretty(pred_stop$years_stopped)) +
  labs(
    x = "Years since treatment ceased",
    y = "Predicted mortality (%)"
  ) +
  theme_minimal(base_size = 12)

p_main1 | p_main2
```
##Probability table: acaricides 
```{r}
### 1) Prediction grid: one row per acaricide at baseline covariates
acar_levels <- levels(dat$acaricide)
newdat_acar <- tibble::tibble(
  years_treated   = 0,
  years_stopped   = 0,
  rotation        = 0,
  mixture         = 0,
  scaled_log_conc = 0,                           # baseline (mean of scaled predictor)
  acaricide       = factor(acar_levels, levels = acar_levels),
  n               = 1L                           # epred => predicted probabilities
)

### 2) Posterior predictions (probabilities) and tidy to long format
draws <- posterior_epred(fit2, newdata = newdat_acar, re_formula = NA)
post  <- tibble::as_tibble(draws) |>
  rlang::set_names(acar_levels) |>
  dplyr::mutate(.draw = dplyr::row_number()) |>
  tidyr::pivot_longer(-.draw, names_to = "Acaricide", values_to = "p")

### 3) Baseline per draw ("None") to compute Δ vs baseline
base_draws <- post |>
  dplyr::filter(Acaricide == "None") |>
  dplyr::select(.draw, p_base = p)

### 4) Summaries by acaricide: mean % and Δ vs baseline (%), both with 95% CrI
tab_abs <- post |>
  dplyr::group_by(Acaricide) |>
  dplyr::summarise(
    MeanPct = mean(p) * 100,
    l95     = quantile(p, 0.025) * 100,
    u95     = quantile(p, 0.975) * 100,
    .groups = "drop"
  )

tab_delta <- post |>
  dplyr::left_join(base_draws, by = ".draw") |>
  dplyr::mutate(DeltaPP = (p - p_base) * 100) |>
  dplyr::group_by(Acaricide) |>
  dplyr::summarise(
    d_mean = mean(DeltaPP),
    d_l95  = quantile(DeltaPP, 0.025),
    d_u95  = quantile(DeltaPP, 0.975),
    .groups = "drop"
  )

### 5) Combine, drop "None", order nicely, format
order_vec <- c("Amitraz", "TauFluvalinate", "Flumethrin", "Coumaphos")

tab_final <- tab_abs |>
  dplyr::inner_join(tab_delta, by = "Acaricide") |>
  dplyr::filter(Acaricide != "None") |>
  dplyr::mutate(Acaricide = factor(Acaricide, levels = order_vec[order_vec %in% Acaricide])) |>
  dplyr::arrange(Acaricide) |>
  dplyr::transmute(
    Acaricide,
    `Mean mite mortality (%) [95% CrI]`       = sprintf("%.1f (%.1f–%.1f)", MeanPct, l95, u95),
    `Difference from baseline (%) [95% CrI]` = sprintf("%+.1f (%+.1f–%+.1f)", d_mean, d_l95, d_u95)
  )

### 6) Baseline ("None") mean mortality for caption
p0_fit2 <- post |>
  dplyr::filter(Acaricide == "None") |>
  dplyr::summarise(m = mean(p)) |>
  dplyr::pull(m)

### 7) Render
knitr::kable(
  tab_final,
  caption = glue::glue(
    "Table. Posterior predicted mite mortality by acaricide at baseline (no rotation/mixture, ",
    "0 years treated/stopped, scaled_log_conc = 0). Baseline (None) ≈ {scales::percent(p0_fit2, 1)}. ",
    "Δ values are percentage points vs baseline."
  )
)
```

##Probability table: IRM + treatment history

```{r}
make_main_effects_table <- function(fit, p0) {
  # 1) Fixed-effect summaries with 95% CrI
  fx <- brms::fixef(fit, probs = c(0.025, 0.975)) |> as.data.frame()
  fx$term <- rownames(fx); rownames(fx) <- NULL

  ### 2) Identify column names for lower/upper CrI (brms may label them differently)
  lower <- intersect(c("l-95% CI","Q2.5","2.5%"), names(fx))[1]
  upper <- intersect(c("u-95% CI","Q97.5","97.5%"), names(fx))[1]

  ### 3) Keep only predictors of interest, and normalise factor names
  keep <- c("years_treated","years_stopped","rotation","mixture")

  fx |>
    dplyr::mutate(term = dplyr::recode(term,
                                       "rotationTRUE" = "rotation",
                                       "mixtureTRUE"  = "mixture")) |>
    dplyr::filter(term %in% keep) |>
    ### 4) Translate log-odds effects to % point changes at baseline p0
    dplyr::transmute(
      key  = term,
      mean = (plogis(qlogis(p0) + Estimate)       - p0) * 100,
      l95  = (plogis(qlogis(p0) + .data[[lower]]) - p0) * 100,
      u95  = (plogis(qlogis(p0) + .data[[upper]]) - p0) * 100
    ) |>
    ### 5) Friendly labels + qualitative interpretation
    dplyr::mutate(
      Predictor = dplyr::recode(key,
        "years_treated" = "Years treated (+1y)",
        "years_stopped" = "Years off (+1y)",
        "rotation"      = "Rotation",
        "mixture"       = "Mixture"
      ),
      `Effect on mortality` = dplyr::case_when(
        l95 > 0  ~ "Increase",
        u95 < 0  ~ "Decrease",
        TRUE     ~ "Uncertain"
      ),
      `Change from baseline (pp)` = sprintf("%+.1f (%+.1f to %+.1f)", mean, l95, u95)
    ) |>
    dplyr::select(Predictor, `Effect on mortality`, `Change from baseline (%)`)
}

### Build table and print
tab_simple <- make_main_effects_table(fit, p0)

knitr::kable(
  tab_simple,
  caption = glue::glue(
    "Effects at baseline (mortality ≈ {scales::percent(p0, 1)}). ",
    "Positive = higher mortality; negative = lower mortality."
  )
)
```
