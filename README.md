# Kenya HIV Treatment Cascade Analysis

**Author:** Michael Mando | Epidemiologist & Biostatistician  
**Tools:** R (tidyverse, ggplot2, survival, broom)  
**Data:** Simulated patient cohort based on Kenya HIV Estimates (NACC/NASCOP 2023)

---

## Overview

This project analyses the **HIV treatment cascade** across six counties in Kenya — examining gaps between people living with HIV (PLHIV), ART coverage, and viral suppression. The analysis is modelled on the type of routine program data encountered in real-world health information systems such as KenyaEMR and DHIS2.

The 95-95-95 UNAIDS targets require that:
- 95% of PLHIV know their status
- 95% of those diagnosed are on ART
- 95% of those on ART achieve viral suppression

This analysis investigates how a facility cohort performs against these benchmarks and what patient-level factors predict viral suppression.

---

## Research Questions

1. What proportion of patients in the cohort are on ART and virally suppressed?
2. How does ART coverage and viral suppression vary by county?
3. What patient-level factors (age, sex, CD4 at entry, year of enrolment) predict viral suppression?
4. Does late presentation (CD4 < 200) affect time to ART initiation?

---

## Methods

| Method | Purpose |
|---|---|
| Descriptive statistics | Cascade summary overall and by county |
| Logistic regression | Predictors of viral suppression (OR, 95% CI) |
| Kaplan-Meier survival analysis | Time to ART initiation by CD4 status |
| Data visualisation | ggplot2 charts for all outputs |

---

## Key Findings

### Cascade Performance
- **86%** of patients in the cohort are on ART
- **~78%** of ART patients received a viral load test
- **~90%** of tested patients achieved viral suppression
- Overall cascade gap is largest at **ART initiation**, not viral suppression

### County Variation
- ART coverage ranges from ~82% to ~90% across counties
- Homa Bay and Kisumu show the largest gaps — consistent with national program data
- Most counties fall below the 95% UNAIDS target for ART coverage

### Predictors of Viral Suppression
- **Late presenters (CD4 < 200)** had lower odds of viral suppression
- **More recent enrolment** (2022–2023) associated with higher suppression — likely reflecting improved regimens and same-day ART initiation policy
- Sex differences in suppression were modest after adjusting for other factors

### Time to ART
- Median time to ART initiation: ~30 days
- **Late presenters** initiated ART significantly faster — consistent with clinical urgency

---

## Project Structure

```
kenya-hiv-cascade/
├── analysis.R          # Full R analysis script (documented & reproducible)
├── plots/
│   ├── 01_hiv_cascade.png       # Waterfall cascade chart
│   ├── 02_county_cascade.png    # County comparison bar chart
│   └── 03_forest_plot.png       # Logistic regression forest plot
└── README.md
```

---

## How to Run

```r
# Install required packages (first time only)
install.packages(c("tidyverse", "ggplot2", "scales", "survival", "broom"))

# Run the full analysis
source("analysis.R")
```

All outputs (plots + console summaries) are generated automatically.

---

## Background & Relevance

This analysis reflects the type of work conducted in HIV program settings across Kenya — particularly in facilities supported by organisations such as **KEMRI-FACES**, **PEPFAR**, and **CDC Kenya**. The cascade framework is the standard tool for HIV program performance monitoring and is directly used in DHIS2 dashboards and county health reports.

The code is written to be **fully reproducible**, with a set seed and documented do-file style — consistent with best practices in biostatistical analysis.

---

## Contact

Michael Mando  
📧 mandomichael254@gmail.com  
🌍 Kisumu, Kenya  
🔗 [GitHub](https://github.com/OneMando42)
