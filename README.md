# Calibration Assessment of Survival Prediction Models

## Authors

Aurora Gaeta, Xinhui Ran, Andrea Bellavia
TIMI Study Group, Brigham and Women's Hospital, Harvard Medical School

## Overview

Calibration of clinical prediction models compares predicted risk with observed risk at a clinically relevant time horizon. For censored time-to-event outcomes, however, observed risk is not directly available for all individuals and must be estimated from the observed data.

This project provides R functions to compare alternative approaches for estimating observed risk during calibration assessment:

- Empirical approaches, using grouped Kaplan–Meier or Aalen–Johansen estimates
- Flexible survival-model-based approaches, using Cox or Fine–Gray calibration models
- Pseudo-value (PV) regression, which directly models observed absolute risk or cumulative incidence at the target time horizon

The PV framework is intended for calibration assessment of existing prediction models, rather than for prediction-model development. The same predicted risks can therefore be evaluated using alternative approaches for estimating observed risk.

## R Implementation

The <a href="calibration_functions.R">R file </a> provides functions for calibration assessment in:

- Standard right-censored survival settings
- Competing-risk settings

The functions support alternative predictor scales, including the original risk and complementary log-log scales, allowing sensitivity analyses of calibration-model specification.

## Motivation

The approach used to estimate observed risk may influence the assessment of model calibration, particularly in regions of the predicted-risk distribution with limited data. It may also affect conclusions regarding the incremental prognostic value of additional predictors, including novel biomarkers.

The accompanying manuscript systematically compares these approaches in standard survival and competing-risk settings and illustrates their application to the evaluation of established clinical prediction models and novel biomarkers.

**Manuscript currently under review.**
