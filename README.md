# Calibration Assessment in Survival Models

## Overview

This project investigates the impact of different calibration assessment methods 
in survival analysis.

We compare three approaches to evaluate absolute risk calibration:

- Cox-based (traditional)
- Pseudo-values (PV-based)
- HARE (flexible hazard regression)

The goal is to demonstrate how the choice of calibration framework influences 
the interpretation of model performance and the perceived value of biomarkers.

---

## Methods

We used the **Rotterdam breast cancer dataset** available in the `survival` R package.

Two Cox proportional hazards models were fitted:

1. **Base model**:
   - age
   - tumor size

2. **Biomarker model**:
   - age
   - tumor size
   - number of positive lymph nodes

Predicted risks at 1 year were obtained using:

- `predictSurvProb()`

---

## Calibration Approaches

Three calibration strategies were implemented:

### 1. Cox-based calibration
- Flexible Cox model using restricted cubic splines
- Observed risks derived from `survest()`
- Represents traditional model-based calibration

### 2. Pseudo-value (PV) calibration
- Uses jackknife pseudo-values via `eventglm`
- Provides data-driven estimates of observed risk
- Less dependent on model assumptions

### 3. HARE calibration
- Hazard regression using `polspline::hare`
- Semi-parametric flexible approach
- Intermediate between Cox and PV

---

## Outputs

The script produces:

- Calibration curves (Cox vs PV vs HARE)
- Histograms of predicted risks
- Decision curve plots (optional)
- Calibration metrics:
  - ICI (Integrated Calibration Index)
  - E50, E90
  - Emax
  - RMSB

---



## Requirements

The script requires the following R packages:


survival
riskRegression
eventglm
dplyr
ggplot2
rms
gtsummary
aod
polspline




