# =============================================================================
# Kenya HIV Treatment Cascade Analysis
# Author: Michael Mando | Epidemiologist & Biostatistician
# Data Source: Simulated data based on Kenya HIV Estimates (NACC/NASCOP 2023)
# Description: Analysis of the HIV treatment cascade across counties in Kenya,
#              examining gaps between PLHIV, ART coverage, and viral suppression.
# =============================================================================

# --- 1. LOAD LIBRARIES -------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(scales)
library(survival)
library(broom)

set.seed(42)

# --- 2. SIMULATE DATASET (based on NACC Kenya HIV Estimates 2023) -------------

n <- 1200  # patients across 6 facilities

counties    <- c("Kisumu", "Nairobi", "Mombasa", "Kisii", "Kakamega", "Homa Bay")
facilities  <- paste0("Facility_", LETTERS[1:6])

hiv_data <- tibble(
  patient_id      = sprintf("KE-%04d", 1:n),
  county          = sample(counties, n, replace = TRUE,
                           prob = c(0.20, 0.25, 0.15, 0.15, 0.13, 0.12)),
  facility        = sample(facilities, n, replace = TRUE),
  age             = round(rnorm(n, mean = 38, sd = 11)),
  sex             = sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.57, 0.43)),
  year_enrolled   = sample(2018:2023, n, replace = TRUE,
                           prob = c(0.10, 0.12, 0.15, 0.18, 0.22, 0.23)),
  on_art          = sample(c(1, 0), n, replace = TRUE, prob = c(0.86, 0.14)),
  vl_tested       = NA_real_,
  vl_suppressed   = NA_real_,
  cd4_at_entry    = round(rnorm(n, mean = 380, sd = 180)),
  transfer_out    = sample(c(0, 1), n, replace = TRUE, prob = c(0.92, 0.08)),
  ltfu            = sample(c(0, 1), n, replace = TRUE, prob = c(0.88, 0.12)),
  time_to_art_days = round(rexp(n, rate = 1/30))  # median ~30 days
) %>%
  mutate(
    age         = pmax(15, pmin(age, 75)),
    cd4_at_entry = pmax(10, cd4_at_entry),
    vl_tested   = ifelse(on_art == 1,
                         sample(c(1, 0), n, replace = TRUE, prob = c(0.78, 0.22)),
                         0),
    vl_suppressed = ifelse(vl_tested == 1,
                           sample(c(1, 0), n, replace = TRUE, prob = c(0.90, 0.10)),
                           0),
    age_group   = cut(age,
                      breaks = c(15, 25, 35, 45, 55, 75),
                      labels = c("15-24", "25-34", "35-44", "45-54", "55+"),
                      right  = FALSE),
    late_presenter = ifelse(cd4_at_entry < 200, 1, 0)
  )

cat("Dataset created:", nrow(hiv_data), "patients across", n_distinct(hiv_data$county), "counties\n")
cat("Variables:", ncol(hiv_data), "\n\n")


# --- 3. CASCADE SUMMARY -------------------------------------------------------

cascade_overall <- tibble(
  Stage = factor(
    c("PLHIV (Estimated)", "Diagnosed", "On ART", "Virally Suppressed"),
    levels = c("PLHIV (Estimated)", "Diagnosed", "On ART", "Virally Suppressed")
  ),
  Count = c(
    n,
    round(n * 0.96),
    sum(hiv_data$on_art),
    sum(hiv_data$vl_suppressed)
  )
) %>%
  mutate(
    Percent = round(Count / Count[1] * 100, 1),
    label   = paste0(comma(Count), "\n(", Percent, "%)")
  )

cat("=== HIV Treatment Cascade ===\n")
print(cascade_overall %>% select(Stage, Count, Percent))


# --- 4. COUNTY-LEVEL CASCADE --------------------------------------------------

county_cascade <- hiv_data %>%
  group_by(county) %>%
  summarise(
    total           = n(),
    on_art_n        = sum(on_art),
    vl_tested_n     = sum(vl_tested),
    vl_suppressed_n = sum(vl_suppressed),
    art_coverage    = round(on_art_n / total * 100, 1),
    vl_suppression  = round(vl_suppressed_n / on_art_n * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(art_coverage))

cat("\n=== County-Level ART Coverage ===\n")
print(county_cascade)


# --- 5. LOGISTIC REGRESSION: PREDICTORS OF VIRAL SUPPRESSION -----------------

model_data <- hiv_data %>%
  filter(on_art == 1, vl_tested == 1) %>%
  mutate(
    sex         = relevel(factor(sex), ref = "Male"),
    age_group   = relevel(factor(age_group), ref = "25-34"),
    county      = relevel(factor(county), ref = "Nairobi")
  )

vl_model <- glm(vl_suppressed ~ age_group + sex + late_presenter + year_enrolled,
                data   = model_data,
                family = binomial(link = "logit"))

model_results <- tidy(vl_model, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    significant = ifelse(p.value < 0.05, "*", ""),
    across(c(estimate, conf.low, conf.high), ~ round(.x, 2)),
    p.value = round(p.value, 3)
  )

cat("\n=== Logistic Regression: Predictors of Viral Suppression ===\n")
cat("Outcome: Viral suppression (VL < 1000 copies/mL)\n")
cat("N =", nrow(model_data), "patients on ART with VL result\n\n")
print(model_results %>% select(term, estimate, conf.low, conf.high, p.value, significant))


# --- 6. SURVIVAL ANALYSIS: TIME TO ART INITIATION ---------------------------

surv_data <- hiv_data %>%
  mutate(
    event      = on_art,
    time       = pmax(1, time_to_art_days),
    late_presenter = factor(late_presenter, labels = c("CD4 ≥ 200", "CD4 < 200 (Late)"))
  )

surv_obj <- Surv(time = surv_data$time, event = surv_data$event)
km_fit   <- survfit(surv_obj ~ late_presenter, data = surv_data)

km_summary <- summary(km_fit, times = c(14, 30, 60, 90))
cat("\n=== KM Survival: ART Initiation by CD4 Status ===\n")
cat("Probability of NOT yet on ART by day:\n")
print(data.frame(
  time        = km_summary$time,
  group       = km_summary$strata,
  survival    = round(km_summary$surv, 3),
  n_at_risk   = km_summary$n.risk
))


# --- 7. VISUALIZATIONS --------------------------------------------------------

dir.create("plots", showWarnings = FALSE)

## Plot 1: Cascade waterfall chart
p1 <- ggplot(cascade_overall, aes(x = Stage, y = Count, fill = Stage)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = label), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("#1F3864", "#2E75B6", "#4BACC6", "#70AD47")) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "HIV Treatment Cascade — Kenya (Facility Cohort)",
    subtitle = "Patient-level cascade from PLHIV to viral suppression",
    x        = NULL,
    y        = "Number of Patients",
    caption  = "Source: Simulated facility data based on NACC/NASCOP 2023 estimates"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", color = "#1F3864"),
    plot.subtitle = element_text(color = "#666666"),
    axis.text.x   = element_text(size = 11),
    panel.grid.major.x = element_blank()
  )

ggsave("plots/01_hiv_cascade.png", p1, width = 9, height = 6, dpi = 150)
cat("\nPlot saved: plots/01_hiv_cascade.png\n")


## Plot 2: County ART coverage
p2 <- county_cascade %>%
  pivot_longer(c(art_coverage, vl_suppression),
               names_to = "Indicator", values_to = "Percent") %>%
  mutate(Indicator = recode(Indicator,
    art_coverage   = "ART Coverage",
    vl_suppression = "Viral Suppression (among ART)"
  )) %>%
  ggplot(aes(x = reorder(county, Percent), y = Percent, fill = Indicator)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text", x = 1.5, y = 96.5, label = "95% UNAIDS Target",
           color = "red", size = 3.2) +
  coord_flip() +
  scale_fill_manual(values = c("#2E75B6", "#70AD47")) +
  scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "ART Coverage & Viral Suppression by County",
    subtitle = "Benchmarked against UNAIDS 95% target",
    x        = NULL, y = "Percentage (%)", fill = NULL,
    caption  = "Source: Simulated facility cohort data"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", color = "#1F3864"),
    legend.position = "bottom"
  )

ggsave("plots/02_county_cascade.png", p2, width = 9, height = 6, dpi = 150)
cat("Plot saved: plots/02_county_cascade.png\n")


## Plot 3: Forest plot — OR for viral suppression
p3 <- model_results %>%
  mutate(term = str_replace_all(term, c(
    "age_group" = "Age: ",
    "sexFemale" = "Sex: Female",
    "late_presenter" = "Late Presenter (CD4<200)",
    "year_enrolled" = "Year Enrolled"
  ))) %>%
  ggplot(aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.25, color = "#2E75B6", linewidth = 0.8) +
  geom_point(size = 3.5, color = "#1F3864") +
  scale_x_log10() +
  labs(
    title    = "Predictors of Viral Suppression — Logistic Regression",
    subtitle = "Odds Ratios with 95% Confidence Intervals (log scale)",
    x        = "Odds Ratio (log scale)", y = NULL,
    caption  = "Reference: Male sex, Age 25-34, CD4 ≥ 200"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", color = "#1F3864"))

ggsave("plots/03_forest_plot.png", p3, width = 9, height = 6, dpi = 150)
cat("Plot saved: plots/03_forest_plot.png\n")

cat("\n✓ Analysis complete. All outputs saved to /plots\n")
