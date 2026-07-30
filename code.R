
required_packages <- c(
  "survival","riskRegression","eventglm",
  "dplyr","ggplot2","rms","gtsummary","aod","polspline", "pec"
)

for(p in required_packages){
  if(!require(p, character.only = TRUE)){
    stop(paste("Package missing:", p))
  }
}


#Survival---------
## PSEUDOVALUES  ---------------------------------------------------
make_calibration = function(p_pred, data, t0, knots = 3, eps = 1e-6, dtime, death) {
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = p_pred
  
  dcal = data.frame(
    time = dtime,
    state   = death,
    eta        = eta
  )
  
  dd  = datadist(dcal)
  old = options(datadist = "dd")
  on.exit(options(old), add = TRUE)
  
  # Fit calibration curve at time t0; identity link returns risks on probability scale
  cal_fit = eventglm::cumincglm(
    Surv(time, state) ~ rms::rcs(eta, knots),
    data = dcal,
    time = t0,
    link = "identity"
  )
  
  # Smoothed "observed" risk for each subject
  obs_hat_i = as.numeric(predict(cal_fit, newdata = data.frame(eta = dcal$eta),
                                 type = "response"))
  obs_hat_i = pmin(pmax(obs_hat_i, 0), 1)
  
  # Calibration curve on a common grid for the plot
  grid_p   = seq(0.001, 0.99, length.out = 200)
  grid_eta = grid_p
  grid_obs = as.numeric(predict(cal_fit, newdata = data.frame(eta = grid_eta),
                                type = "response"))
  grid_obs = pmin(pmax(grid_obs, 0), 1)
  
  list(
    cal_fit   = cal_fit,   #Fitted calibration curve at time t0
    obs_hat_i = obs_hat_i, #Smoothed "OBSERVED" risk for each subject
    grid_p    = grid_p,    #Grid of x values
    grid_obs  = grid_obs   #Y given f(x) of the grid
  )
}

## COX ---------------------------------------------------
make_calibration_cox = function(p_pred, data, t0, knots = 3, eps = 1e-6, dtime, death) {
  
  
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = log(-log(1 - p_pred)) # Transform probability to the Cox linear scale (cloglog)
  
  dcal = data.frame(
    time =  dtime,
    state   = death,
    eta        = eta
  )
  
  # Fit a Cox model with a flexible spline for the score
  # This is the "Traditional" recalibration/smoothing step
  cal_fit_cox = rms::cph(
    Surv(time, state) ~ rms::rcs(eta, knots),
    data = dcal,
    x = TRUE, y = TRUE, surv = TRUE
  )
  
  # Create a common grid for the X-axis (Predicted Risk)
  grid_p   = seq(0.001, 0.99, length.out = 200)
  grid_eta = log(-log(1 - grid_p))
  
  # Use the Cox model to predict "Observed" risk at time t0
  # survest returns Survival, so 1 - Survival = Risk
  grid_surv = rms::survest(cal_fit_cox,
                           newdata = data.frame(eta = grid_eta),
                           times = t0, conf.int = 0)$surv
  
  grid_obs = 1 - as.numeric(grid_surv)
  grid_obs = pmin(pmax(grid_obs, 0), 1)
  
  list(
    cal_fit  = cal_fit_cox,
    grid_p   = grid_p,    # X-axis
    grid_obs = grid_obs   # Y-axis (Cox-judged reality)
  )
}

## KM-based Assessment (group based) ---------
make_calibration_km_loess <- function(
    p_pred,
    data,
    t0,
    n_groups = 10,
    span = 0.5,
    eps = 1e-6
) {
  
  # stabilize predictions ---
  p_pred <- pmin(pmax(p_pred, eps), 1 - eps)
  
  dcal <- data.frame(
    time   = data$days2miistr,
    status = data$miistrfu,
    pred   = p_pred
  )
  
  # quantile-based grouping ---
  dcal$group <- cut(
    dcal$pred,
    breaks = quantile(
      dcal$pred,
      probs = seq(0, 1, length.out = n_groups + 1),
      na.rm = TRUE
    ),
    include.lowest = TRUE,
    labels = FALSE
  )
  
  # KM estimate within each group ---
  df_group <- do.call(
    rbind,
    lapply(split(dcal, dcal$group), function(df) {
      
      fit <- survfit(
        Surv(time, status) ~ 1,
        data = df,
        conf.type = "log-log"
      )
      
      s <- summary(fit, times = t0)
      
      if (length(s$surv) == 0) {
        
        obs <- NA
        low <- NA
        high <- NA
        
      } else {
        
        obs <- 1 - s$surv
        low <- 1 - s$upper
        high <- 1 - s$lower
        
      }
      
      data.frame(
        pred_mean = mean(df$pred),
        obs       = obs,
        obs_lower = low,
        obs_upper = high,
        n         = nrow(df)
      )
    })
  )
  
  # Remove potential NA groups
  df_group <- df_group[complete.cases(df_group[, c(
    "pred_mean", "obs"
  )]), ]
  
  # LOESS smoothing ---
  loess_fit <- loess(
    obs ~ pred_mean,
    data = df_group,
    weights = df_group$n,
    span = span,
    degree = 1,
    control = loess.control(surface = "direct")
  )
  
  # smooth grid ---
  grid_p <- seq(
    min(p_pred),
    max(p_pred),
    length.out = 200
  )
  
  grid_obs <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = grid_p)
  )
  
  grid_obs <- pmin(pmax(grid_obs, 0), 1)
  
  # subject-level fitted observed risk ---
  obs_hat_i <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = p_pred)
  )
  
  obs_hat_i <- pmin(pmax(obs_hat_i, 0), 1)
  
  # return ---
  list(
    df_group  = df_group,
    cal_fit   = loess_fit,
    obs_hat_i = obs_hat_i,
    grid_p    = grid_p,
    grid_obs  = grid_obs
  )
}

#Competing risk setting -----
##PSEUDOVALUES
make_calibration_pv_cpr = function(p_pred, data, t0, knots = 3, eps = 1e-6) {
  
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = cloglog(p_pred)
  
  dcal = data.frame(
    etime = data$etime,
    eventnewc   = data$eventnewc,
    eta        = eta
  )
  
  dd  = datadist(dcal)
  old = options(datadist = "dd")
  on.exit(options(old), add = TRUE)
  
  # Fit calibration curve at time t0; identity link returns risks on probability scale
  cal_fit = eventglm::cumincglm(
    Surv(etime, eventnewc) ~ rms::rcs(eta, knots),
    data = dcal,
    time = t0,
    cause="Event",
    link = "identity"
  )
  
  # Smoothed "observed" risk for each subject
  obs_hat_i = as.numeric(predict(cal_fit, newdata = data.frame(eta = dcal$eta), type = "response"))
  obs_hat_i = pmin(pmax(obs_hat_i, 0), 1)
  
  # Calibration curve on a common grid
  grid_p   = seq(0.001, 0.99, length.out = 200)
  grid_eta = cloglog(grid_p)
  
  grid_obs = as.numeric(predict(cal_fit, newdata = data.frame(eta = grid_eta), type = "response"))
  grid_obs = pmin(pmax(grid_obs, 0), 1)
  
  list(
    cal_fit   = cal_fit,   #Fitted calibration curve at time t0
    obs_hat_i = obs_hat_i, #Smoothed "observed" risk for each subject
    grid_p    = grid_p,    #Grid of x values
    grid_obs  = grid_obs   #Y given f(x)
  )
}

##FG-based Assessment----------------------------------
make_calibration_fg <- function(
    p_pred,
    data,
    t0,
    knots = 3,
    eps = 1e-6
) {
  
  library(riskRegression)
  library(prodlim)
  library(Hmisc)
  
  # --------------------------------------------------
  # predicted risk -> complementary log-log
  # --------------------------------------------------
  p_pred <- pmin(pmax(p_pred, eps), 1 - eps)
  
  eta <- log(-log(1 - p_pred))
  
  dcal <- data.frame(
    time   = data$etime,
    status = data$eventnew,
    eta    = eta
  )
  
  # --------------------------------------------------
  # determine knot locations ONCE
  # --------------------------------------------------
  knot_locations <- Hmisc::rcspline.eval(
    eta,
    nk = knots,
    knots.only = TRUE
  )
  
  # --------------------------------------------------
  # spline basis for fitting
  # --------------------------------------------------
  X_fit <- Hmisc::rcspline.eval(
    eta,
    knots = knot_locations,
    inclx = TRUE
  )
  
  X_fit <- as.data.frame(X_fit)
  
  spline_names <- paste0("s", seq_len(ncol(X_fit)))
  names(X_fit) <- spline_names
  
  dcal_fit <- cbind(
    dcal[, c("time", "status")],
    X_fit
  )
  
  # --------------------------------------------------
  # Fine-Gray calibration model
  # --------------------------------------------------
  rhs <- paste(spline_names, collapse = " + ")
  
  fg_formula <- as.formula(
    paste("Hist(time,status) ~", rhs)
  )
  
  cal_fit_fg <- riskRegression::FGR(
    formula = fg_formula,
    data    = dcal_fit,
    cause   = 1
  )
  
  # --------------------------------------------------
  # observed risk for each subject
  # --------------------------------------------------
  X_subject <- Hmisc::rcspline.eval(
    eta,
    knots = knot_locations,
    inclx = TRUE
  )
  
  X_subject <- as.data.frame(X_subject)
  names(X_subject) <- spline_names
  
  obs_hat_i <- as.numeric(
    riskRegression::predictRisk(
      cal_fit_fg,
      newdata = X_subject,
      times   = t0,
      cause   = 1
    )
  )
  
  # --------------------------------------------------
  # grid for calibration plot
  # --------------------------------------------------
  grid_p <- seq(0.001, 0.99, length.out = 200)
  
  grid_eta <- log(-log(1 - grid_p))
  
  X_grid <- Hmisc::rcspline.eval(
    grid_eta,
    knots = knot_locations,
    inclx = TRUE
  )
  
  X_grid <- as.data.frame(X_grid)
  names(X_grid) <- spline_names
  
  grid_obs <- as.numeric(
    riskRegression::predictRisk(
      cal_fit_fg,
      newdata = X_grid,
      times   = t0,
      cause   = 1
    )
  )
  
  grid_obs <- pmin(pmax(grid_obs, 0), 1)
  
  # --------------------------------------------------
  # Return object
  # --------------------------------------------------
  list(
    cal_fit        = cal_fit_fg,
    obs_hat_i      = obs_hat_i,
    grid_p         = grid_p,
    grid_obs       = grid_obs,
    knot_locations = knot_locations
  )
}


##  AJ-based Assessment (group based) in competing risk setting --------

make_calibration_aj_smooth <- function(
    p_pred,
    data,
    t0,
    n_groups = 10,
    span = 0.5,
    eps = 1e-6
) {
  
  library(prodlim)
  
  # stabilize predictions ---
  p_pred <- pmin(pmax(p_pred, eps), 1 - eps)
  
  dcal <- data.frame(
    time   = data$etime,
    status = data$eventnew,
    pred   = p_pred
  )
  
  # quantile-based grouping ---
  cuts <- unique(
    quantile(
      dcal$pred,
      probs = seq(0, 1, length.out = n_groups + 1),
      na.rm = TRUE
    )
  )
  
  dcal$group <- cut(
    dcal$pred,
    breaks = cuts,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  # AJ estimate within each group ---
  df_group <- do.call(
    rbind,
    lapply(split(dcal, dcal$group), function(df) {
      
      fit <- prodlim(
        Hist(time, status) ~ 1,
        data = df
      )
      
      s <- summary(fit, times = t0)
      
      s1 <- s[s$cause == 1, ]
      
      if (nrow(s1) == 0) {
        
        data.frame(
          pred_mean = mean(df$pred),
          obs       = NA_real_,
          obs_lower = NA_real_,
          obs_upper = NA_real_,
          n         = nrow(df)
        )
        
      } else {
        
        data.frame(
          pred_mean = mean(df$pred),
          obs       = s1$cuminc,
          obs_lower = s1$lower,
          obs_upper = s1$upper,
          n         = nrow(df)
        )
      }
    })
  )
  
  # Remove potential NA groups
  df_group <- df_group[
    complete.cases(df_group[, c("pred_mean", "obs")]),
  ]
  
  df_group <- df_group[order(df_group$pred_mean), ]
  
  # LOESS smoothing ---
  loess_fit <- loess(
    obs ~ pred_mean,
    data = df_group,
    weights = df_group$n,
    span = span,
    degree = 1,
    control = loess.control(surface = "direct")
  )
  
  # smooth grid ---
  grid_p <- seq(
    min(p_pred),
    max(p_pred),
    length.out = 200
  )
  
  grid_obs <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = grid_p)
  )
  
  grid_obs <- pmin(pmax(grid_obs, 0), 1)
  
  # subject-level fitted observed risk ---
  obs_hat_i <- predict(
    loess_fit,
    newdata = data.frame(pred_mean = p_pred)
  )
  
  obs_hat_i <- pmin(pmax(obs_hat_i, 0), 1)
  
  # return ---
  list(
    df_group  = df_group,
    cal_fit   = loess_fit,
    obs_hat_i = obs_hat_i,
    grid_p    = grid_p,
    grid_obs  = grid_obs
  )
}


