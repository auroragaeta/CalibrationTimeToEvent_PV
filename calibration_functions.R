# Calibration assessment functions for censored time-to-event outcomes
# -----------------------------------------------------------------------------
# These functions compare alternative approaches for estimating the observed
# risk used in calibration assessment:
#   1. pseudo-value regression;
#   2. smooth Cox or Fine-Gray calibration models;
#   3. grouped Kaplan-Meier or Aalen-Johansen estimates with LOESS smoothing.
#
# Predicted risks are assumed to correspond to the same time horizon, t0,
# used for calibration assessment.

required_packages <- c(
  "survival",
  "riskRegression",
  "eventglm",
  "rms",
  "prodlim",
  "Hmisc"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# Restrict probabilities away from 0 and 1 before applying transformations.
bound_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

# Complementary log-log transformation.
cloglog <- function(p) {
  log(-log(1 - p))
}

# Restrict fitted risks to the probability scale.
bound_fitted_risk <- function(p) {
  pmin(pmax(as.numeric(p), 0), 1)
}

# Apply the requested predictor transformation.
transform_predicted_risk <- function(p, transform = c("risk", "cloglog")) {
  transform <- match.arg(transform)

  if (transform == "cloglog") {
    cloglog(p)
  } else {
    p
  }
}


# -----------------------------------------------------------------------------
# 1. Pseudo-value calibration: standard survival outcome
# -----------------------------------------------------------------------------

make_calibration_pv <- function(
    p_pred,
    time,
    status,
    t0,
    knots = 3,
    transform = c("risk", "cloglog"),
    eps = 1e-6
) {
  transform <- match.arg(transform)

  # Predicted risk is used as the predictor; pseudo-values represent observed
  # event risk at the target horizon.
  p_pred <- bound_probability(p_pred, eps)
  eta <- transform_predicted_risk(p_pred, transform)

  calibration_data <- data.frame(
    time = time,
    status = status,
    eta = eta
  )

  # Fit the same restricted cubic spline specification used in the
  # conventional smooth calibration model.
  calibration_fit <- eventglm::cumincglm(
    survival::Surv(time, status) ~ rms::rcs(eta, knots),
    data = calibration_data,
    time = t0,
    link = "identity"
  )

  # Subject-level fitted observed risks, used for calibration metrics.
  observed_risk <- predict(
    calibration_fit,
    newdata = data.frame(eta = eta),
    type = "response"
  )
  observed_risk <- bound_fitted_risk(observed_risk)

  # Smooth calibration curve over the observed prediction range only.
  grid_pred <- seq(
    min(p_pred, na.rm = TRUE),
    max(p_pred, na.rm = TRUE),
    length.out = 200
  )
  grid_eta <- transform_predicted_risk(grid_pred, transform)

  grid_observed <- predict(
    calibration_fit,
    newdata = data.frame(eta = grid_eta),
    type = "response"
  )
  grid_observed <- bound_fitted_risk(grid_observed)

  list(
    cal_fit = calibration_fit,
    obs_hat_i = observed_risk,
    grid_p = grid_pred,
    grid_obs = grid_observed,
    transform = transform
  )
}


# -----------------------------------------------------------------------------
# 2. Conventional smooth Cox calibration model
# -----------------------------------------------------------------------------

make_calibration_cox <- function(
    p_pred,
    time,
    status,
    t0,
    knots = 3,
    transform = c("cloglog", "risk"),
    eps = 1e-6
) {
  transform <- match.arg(transform)

  p_pred <- bound_probability(p_pred, eps)
  eta <- transform_predicted_risk(p_pred, transform)

  calibration_data <- data.frame(
    time = time,
    status = status,
    eta = eta
  )

  calibration_fit <- rms::cph(
    survival::Surv(time, status) ~ rms::rcs(eta, knots),
    data = calibration_data,
    x = TRUE,
    y = TRUE,
    surv = TRUE
  )

  # Subject-level fitted observed risks, used for calibration metrics.
  subject_survival <- rms::survest(
    calibration_fit,
    newdata = data.frame(eta = eta),
    times = t0,
    conf.int = 0
  )$surv
  observed_risk <- bound_fitted_risk(1 - subject_survival)

  # Smooth calibration curve over the observed prediction range only.
  grid_pred <- seq(
    min(p_pred, na.rm = TRUE),
    max(p_pred, na.rm = TRUE),
    length.out = 200
  )
  grid_eta <- transform_predicted_risk(grid_pred, transform)

  grid_survival <- rms::survest(
    calibration_fit,
    newdata = data.frame(eta = grid_eta),
    times = t0,
    conf.int = 0
  )$surv
  grid_observed <- bound_fitted_risk(1 - grid_survival)

  list(
    cal_fit = calibration_fit,
    obs_hat_i = observed_risk,
    grid_p = grid_pred,
    grid_obs = grid_observed,
    transform = transform
  )
}


# -----------------------------------------------------------------------------
# 3. Grouped Kaplan-Meier estimates with LOESS smoothing
# -----------------------------------------------------------------------------

make_calibration_km <- function(
    p_pred,
    time,
    status,
    t0,
    n_groups = 10,
    span = 0.5,
    eps = 1e-6
) {
  p_pred <- bound_probability(p_pred, eps)

  calibration_data <- data.frame(
    time = time,
    status = status,
    pred = p_pred
  )

  # Quantile groups are used only as an empirical reference.
  cut_points <- unique(
    stats::quantile(
      calibration_data$pred,
      probs = seq(0, 1, length.out = n_groups + 1),
      na.rm = TRUE
    )
  )

  if (length(cut_points) < 3) {
    stop("Predicted risks do not contain enough unique values to form groups.")
  }

  calibration_data$group <- cut(
    calibration_data$pred,
    breaks = cut_points,
    include.lowest = TRUE,
    labels = FALSE
  )

  grouped_estimates <- do.call(
    rbind,
    lapply(split(calibration_data, calibration_data$group), function(group_data) {
      fit <- survival::survfit(
        survival::Surv(time, status) ~ 1,
        data = group_data,
        conf.type = "log-log"
      )

      estimate <- summary(fit, times = t0, extend = TRUE)

      data.frame(
        pred_mean = mean(group_data$pred, na.rm = TRUE),
        obs = 1 - estimate$surv,
        obs_lower = 1 - estimate$upper,
        obs_upper = 1 - estimate$lower,
        n = nrow(group_data)
      )
    })
  )

  grouped_estimates <- grouped_estimates[
    stats::complete.cases(grouped_estimates[, c("pred_mean", "obs")]),
  ]
  grouped_estimates <- grouped_estimates[order(grouped_estimates$pred_mean), ]

  loess_fit <- stats::loess(
    obs ~ pred_mean,
    data = grouped_estimates,
    weights = n,
    span = span,
    degree = 1,
    control = stats::loess.control(surface = "direct")
  )

  grid_pred <- seq(
    min(p_pred, na.rm = TRUE),
    max(p_pred, na.rm = TRUE),
    length.out = 200
  )

  grid_observed <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = grid_pred)
  )
  grid_observed <- bound_fitted_risk(grid_observed)

  observed_risk <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = p_pred)
  )
  observed_risk <- bound_fitted_risk(observed_risk)

  list(
    df_group = grouped_estimates,
    cal_fit = loess_fit,
    obs_hat_i = observed_risk,
    grid_p = grid_pred,
    grid_obs = grid_observed
  )
}


# -----------------------------------------------------------------------------
# 4. Pseudo-value calibration: competing-risk outcome
# -----------------------------------------------------------------------------

make_calibration_pv_competing <- function(
    p_pred,
    time,
    status,
    cause,
    t0,
    knots = 3,
    transform = c("risk", "cloglog"),
    eps = 1e-6
) {
  transform <- match.arg(transform)

  p_pred <- bound_probability(p_pred, eps)
  eta <- transform_predicted_risk(p_pred, transform)

  calibration_data <- data.frame(
    time = time,
    status = status,
    eta = eta
  )

  calibration_fit <- eventglm::cumincglm(
    survival::Surv(time, status) ~ rms::rcs(eta, knots),
    data = calibration_data,
    time = t0,
    cause = cause,
    link = "identity"
  )

  observed_risk <- predict(
    calibration_fit,
    newdata = data.frame(eta = eta),
    type = "response"
  )
  observed_risk <- bound_fitted_risk(observed_risk)

  grid_pred <- seq(
    min(p_pred, na.rm = TRUE),
    max(p_pred, na.rm = TRUE),
    length.out = 200
  )
  grid_eta <- transform_predicted_risk(grid_pred, transform)

  grid_observed <- predict(
    calibration_fit,
    newdata = data.frame(eta = grid_eta),
    type = "response"
  )
  grid_observed <- bound_fitted_risk(grid_observed)

  list(
    cal_fit = calibration_fit,
    obs_hat_i = observed_risk,
    grid_p = grid_pred,
    grid_obs = grid_observed,
    transform = transform
  )
}


# -----------------------------------------------------------------------------
# 5. Conventional smooth Fine-Gray calibration model
# -----------------------------------------------------------------------------

make_calibration_fg <- function(
    p_pred,
    time,
    status,
    cause,
    t0,
    knots = 3,
    transform = c("cloglog", "risk"),
    eps = 1e-6
) {
  transform <- match.arg(transform)

  p_pred <- bound_probability(p_pred, eps)
  eta <- transform_predicted_risk(p_pred, transform)

  # Determine knot locations once and reuse them for the subject-level
  # predictions and plotting grid.
  knot_locations <- Hmisc::rcspline.eval(
    eta,
    nk = knots,
    knots.only = TRUE
  )

  make_spline_basis <- function(x) {
    basis <- Hmisc::rcspline.eval(
      x,
      knots = knot_locations,
      inclx = TRUE
    )
    basis <- as.data.frame(basis)
    names(basis) <- paste0("s", seq_len(ncol(basis)))
    basis
  }

  subject_basis <- make_spline_basis(eta)

  fit_data <- cbind(
    data.frame(time = time, status = status),
    subject_basis
  )

  spline_terms <- names(subject_basis)
  fg_formula <- stats::as.formula(
    paste(
      "prodlim::Hist(time, status) ~",
      paste(spline_terms, collapse = " + ")
    )
  )

  calibration_fit <- riskRegression::FGR(
    formula = fg_formula,
    data = fit_data,
    cause = cause
  )

  observed_risk <- riskRegression::predictRisk(
    calibration_fit,
    newdata = subject_basis,
    times = t0,
    cause = cause
  )
  observed_risk <- bound_fitted_risk(observed_risk)

  grid_pred <- seq(
    min(p_pred, na.rm = TRUE),
    max(p_pred, na.rm = TRUE),
    length.out = 200
  )
  grid_basis <- make_spline_basis(
    transform_predicted_risk(grid_pred, transform)
  )

  grid_observed <- riskRegression::predictRisk(
    calibration_fit,
    newdata = grid_basis,
    times = t0,
    cause = cause
  )
  grid_observed <- bound_fitted_risk(grid_observed)

  list(
    cal_fit = calibration_fit,
    obs_hat_i = observed_risk,
    grid_p = grid_pred,
    grid_obs = grid_observed,
    knot_locations = knot_locations,
    transform = transform
  )
}


# -----------------------------------------------------------------------------
# 6. Grouped Aalen-Johansen estimates with LOESS smoothing
# -----------------------------------------------------------------------------

make_calibration_aj <- function(
    p_pred,
    time,
    status,
    cause,
    t0,
    n_groups = 10,
    span = 0.5,
    eps = 1e-6
) {
  p_pred <- bound_probability(p_pred, eps)

  calibration_data <- data.frame(
    time = time,
    status = status,
    pred = p_pred
  )

  cut_points <- unique(
    stats::quantile(
      calibration_data$pred,
      probs = seq(0, 1, length.out = n_groups + 1),
      na.rm = TRUE
    )
  )

  if (length(cut_points) < 3) {
    stop("Predicted risks do not contain enough unique values to form groups.")
  }

  calibration_data$group <- cut(
    calibration_data$pred,
    breaks = cut_points,
    include.lowest = TRUE,
    labels = FALSE
  )

  grouped_estimates <- do.call(
    rbind,
    lapply(split(calibration_data, calibration_data$group), function(group_data) {
      fit <- prodlim::prodlim(
        prodlim::Hist(time, status) ~ 1,
        data = group_data
      )

      estimate <- summary(fit, times = t0)
      cause_estimate <- estimate[estimate$cause == cause, ]

      if (nrow(cause_estimate) == 0) {
        return(data.frame(
          pred_mean = mean(group_data$pred, na.rm = TRUE),
          obs = NA_real_,
          obs_lower = NA_real_,
          obs_upper = NA_real_,
          n = nrow(group_data)
        ))
      }

      data.frame(
        pred_mean = mean(group_data$pred, na.rm = TRUE),
        obs = cause_estimate$cuminc[1],
        obs_lower = cause_estimate$lower[1],
        obs_upper = cause_estimate$upper[1],
        n = nrow(group_data)
      )
    })
  )

  grouped_estimates <- grouped_estimates[
    stats::complete.cases(grouped_estimates[, c("pred_mean", "obs")]),
  ]
  grouped_estimates <- grouped_estimates[order(grouped_estimates$pred_mean), ]

  loess_fit <- stats::loess(
    obs ~ pred_mean,
    data = grouped_estimates,
    weights = n,
    span = span,
    degree = 1,
    control = stats::loess.control(surface = "direct")
  )

  grid_pred <- seq(
    min(p_pred, na.rm = TRUE),
    max(p_pred, na.rm = TRUE),
    length.out = 200
  )

  grid_observed <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = grid_pred)
  )
  grid_observed <- bound_fitted_risk(grid_observed)

  observed_risk <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = p_pred)
  )
  observed_risk <- bound_fitted_risk(observed_risk)

  list(
    df_group = grouped_estimates,
    cal_fit = loess_fit,
    obs_hat_i = observed_risk,
    grid_p = grid_pred,
    grid_obs = grid_observed
  )
}


# -----------------------------------------------------------------------------
# Example calls: pseudo-value calibration under alternative predictor scales
# -----------------------------------------------------------------------------

# Standard survival outcome, original predicted-risk scale
# pv_calibration_risk <- make_calibration_pv(
#   p_pred = predicted_risk,
#   time = validation_data$follow_up_time,
#   status = validation_data$event,
#   t0 = 3 * 365.25,
#   knots = 3,
#   transform = "risk"
# )

# Standard survival outcome, complementary log-log scale
# pv_calibration_cloglog <- make_calibration_pv(
#   p_pred = predicted_risk,
#   time = validation_data$follow_up_time,
#   status = validation_data$event,
#   t0 = 3 * 365.25,
#   knots = 3,
#   transform = "cloglog"
# )

# Competing-risk outcome, original predicted-risk scale
# pv_competing_risk <- make_calibration_pv_competing(
#   p_pred = predicted_cumulative_incidence,
#   time = validation_data$follow_up_time,
#   status = validation_data$event_code,
#   cause = 1,
#   t0 = 3 * 365.25,
#   knots = 3,
#   transform = "risk"
# )

# Competing-risk outcome, complementary log-log scale
# pv_competing_cloglog <- make_calibration_pv_competing(
#   p_pred = predicted_cumulative_incidence,
#   time = validation_data$follow_up_time,
#   status = validation_data$event_code,
#   cause = 1,
#   t0 = 3 * 365.25,
#   knots = 3,
#   transform = "cloglog"
# )
