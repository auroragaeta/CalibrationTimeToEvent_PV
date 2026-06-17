
# SETUP------------
required_packages <- c(
  "survival","riskRegression","eventglm",
  "dplyr","ggplot2","rms","gtsummary","aod","polspline"
)

for(p in required_packages){
  if(!require(p, character.only = TRUE)){
    stop(paste("Package missing:", p))
  }
}
    
    
library(survival)
library(riskRegression)
library(eventglm)
library(dplyr)
library(ggplot2)
library(rms)
library(gtsummary)
library(aod)
library(survival)
library(polspline)

## DATA LOADING-------
set.seed(123)
names(rotterdam)
n <- nrow(rotterdam)
id_train <- sample(seq_len(n), size = round(0.7 * n))
train <- rotterdam[id_train, ]
test  <- rotterdam[-id_train, ]



## BASE COX MODEL--------
mod_base <- coxph(
  Surv(dtime, death) ~ age + size  , data = train,
  x = TRUE, y = TRUE
)
mod_base %>% tbl_regression(exponentiate = TRUE)

## BIOMARKER MODEL--------
mod_biom <- coxph(
  Surv(dtime, death) ~ age + size + nodes ,
  data = train,
  x = TRUE, y = TRUE
)
mod_biom %>% tbl_regression()




## PREDICTED RISKS--------
t0 <- 365.25
test$prob_base <- 1 - predictSurvProb(mod_base, newdata = test, times = t0)
test$prob_biom <- 1 - predictSurvProb(mod_biom, newdata = test, times = t0)
prob_base <- as.numeric(test$prob_base)
prob_biom <- as.numeric(test$prob_biom)


# CALIBRATION FUNCTIONS--------
cloglog <- function(p) log(-log(1 - p))
##MAKE CALIBRATION PSEUDOVALUES  ---------------------------------------------------
make_calibration = function(p_pred, data, t0, knots = 3, eps = 1e-6, dtime, death) {
  # p_pred = p_base
  # data = checkt
  # t0 = 365.25*3
  # knots=5
  
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = cloglog(p_pred)
  
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
  grid_eta = cloglog(grid_p)
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

##MAKE CALIBRATION COX ---------------------------------------------------
make_calibration_cox = function(p_pred, data, t0, knots = 3, eps = 1e-6, dtime, death) {
  # 1. Transform probability to the Cox linear scale (cloglog)
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = log(-log(1 - p_pred)) # This is the cloglog/Cox link
  
  dcal = data.frame(
    time =  dtime,
    state   = death,
    eta        = eta
  )
  
  # 2. Fit a Cox model with a flexible spline for the score
  # This is the "Traditional" recalibration/smoothing step
  cal_fit_cox = rms::cph(
    Surv(time, state) ~ rms::rcs(eta, knots),
    data = dcal,
    x = TRUE, y = TRUE, surv = TRUE
  )
  
  # 3. Create a common grid for the X-axis (Predicted Risk)
  grid_p   = seq(0.001, 0.99, length.out = 200)
  grid_eta = log(-log(1 - grid_p))
  
  # 4. Use the Cox model to predict "Observed" risk at time t0
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

#MAKE CALIBRATION HARE  -----
make_calibration_hare <- function(p_pred, data, t0  , dtime, death, maxdim = 5) {
  # clamp numerico
  # p_pred=p_base
  # data=test
  # t0=365.25*3
  # maxdim = 5
  eps = 1e-6
  p_pred <- pmin(pmax(p_pred, eps), 1 - eps)
  eta <- cloglog(p_pred)
  
  dcal <- data.frame(
    dtime = dtime,
    death   = death,
    eta        = eta
  )
  cc <- complete.cases(dcal$dtime, dcal$death, dcal$eta)
  dcc <- dcal[cc, , drop = FALSE]
  cal_fit =polspline::hare(
    data  = dcc$dtime,
    delta = dcc$death,
    cov   = matrix(dcc$eta, ncol = 1),
    maxdim = maxdim)
  
  # Smoothed observed risk per soggetto: P(T <= t0 | eta_i)
  obs_hat_i_cc <- as.numeric(
    polspline::phare(
      q   = t0,
      cov = matrix(dcc$eta, ncol = 1),
      fit = cal_fit
    )
    
  )
  
  obs_hat_i_cc <- pmin(pmax(obs_hat_i_cc, 0), 1)
  obs_hat_i <- rep(NA_real_, nrow(dcal))
  obs_hat_i[cc] <- obs_hat_i_cc
  grid_p   <- seq(0.001, 0.99, length.out = 200)
  grid_eta <- cloglog(grid_p)
  grid_obs <- as.numeric(
    polspline::phare(
      q   = t0,
      cov = matrix(grid_eta, ncol = 1),
      fit = cal_fit
      
    )
    
  )
  grid_obs <- pmin(pmax(grid_obs, 0), 1)
  list(
    obs_hat_i_cc=obs_hat_i_cc,
    cal_fit   = cal_fit,
    obs_hat_i = obs_hat_i,
    grid_p    = grid_p,
    grid_obs  = grid_obs,
    cc        = cc
  )
  
}




# CALIBRATION RESULTS-----------

dd <- datadist(test)
options(datadist = "dd")

cal_base <- make_calibration(prob_base, test, t0, dtime=test$dtime, death=test$death)
cal_biomm <- make_calibration(prob_biom, test, t0,  dtime=test$dtime, death=test$death)

cal_base_cox <- make_calibration_cox(prob_base, test, t0,  dtime=test$dtime, death=test$death)
cal_biomm_cox <- make_calibration_cox(prob_biom, test, t0,  dtime=test$dtime, death=test$death)

cal_base_hare <- make_calibration_hare(prob_base, data = test, t0 = t0,  dtime=test$dtime, death=test$death)
cal_bio_hare  <- make_calibration_hare(prob_biom,  data = test, t0 = t0,  dtime=test$dtime, death=test$death)


# PLOTS
df_base <- data.frame(pred = cal_base$grid_p, 
                      obs = cal_base$grid_obs)
df_biom <- data.frame(pred = cal_biomm$grid_p, 
                      obs = cal_biomm$grid_obs)

#FIGURES------------

# Extract Traditional (Cox-based) curve for Base Model
df_base_cox = data.frame(
  pred = cal_base_cox$grid_p,
  obs  = cal_base_cox$grid_obs,
  method = "Traditional (Cox)"
)

# Extract Traditional (Cox-based) curve for Biomarker Model
df_biom_cox = data.frame(
  pred = cal_biomm_cox$grid_p,
  obs  = cal_biomm_cox$grid_obs,
  method = "Traditional (Cox)"
)


# Extract PV-based curves
df_base_pv = data.frame(
  pred = cal_base$grid_p,
  obs = cal_base$grid_obs,
  method = "New (PV)")

df_biom_pv = data.frame(
  pred = cal_biomm$grid_p,
  obs = cal_biomm$grid_obs,
  method = "New (PV)")


# Extract HARE-based curves
df_base_HARE = data.frame(
  pred = cal_base_hare$grid_p,
  obs = cal_base_hare$grid_obs,
  method = "New (PV)")

df_biom_HARE = data.frame(
  pred = cal_bio_hare$grid_p,
  obs = cal_bio_hare$grid_obs,
  method = "New (PV)")




test$prob_cox_base<-1-predictSurvProb(mod_base,
                                        newdata=test,
                                        times=365.25)

test$prob_cox_biom<-1-predictSurvProb(mod_biom,
                                        newdata=test,
                                        times=365.25)


prob_cox_base=as.numeric(test$prob_cox_base)
prob_cox_biom=as.numeric(test$prob_cox_biom)


# Create a small dataframe for the rug (individual patient predictions)
df_rug_base = data.frame(pred = prob_cox_base)

#------------

#PLOT------
## with hare -----
plot_cal_with_hist_hare_added <- function(curve_cox,
                               curve_pv,
                               curve_hare,
                               pred_cox, pred_pv, pred_hare,
                               main_title,
                               x_zoom = c(0,0.1),
                               breaks = seq(0,0.1,length.out=30)) {
  

  # FILTER PREDICTIONS

  v1 <- pred_cox[!is.na(pred_cox) & pred_cox >= x_zoom[1] & pred_cox <= x_zoom[2]]
  v2 <- pred_pv [!is.na(pred_pv)  & pred_pv  >= x_zoom[1] & pred_pv  <= x_zoom[2]]
  v3 <- pred_hare[!is.na(pred_hare) & pred_hare >= x_zoom[1] & pred_hare <= x_zoom[2]]
  
  h1 <- hist(v1, breaks = breaks, plot = FALSE)
  h2 <- hist(v2, breaks = breaks, plot = FALSE)
  h3 <- hist(v3, breaks = breaks, plot = FALSE)
  
  op <- par(no.readonly = TRUE)
  layout(matrix(c(1,2), nrow = 2), heights = c(4,1))
  

  # PANEL 1: CALIBRATION

  par(mar = c(4,4,3,1))
  
  plot(curve_cox$grid_p, curve_cox$grid_obs,
       type="l", lwd=2, col="#4D4D4D",
       xlim=x_zoom, ylim=x_zoom,
       xlab="Predicted risk (3-year)",
       ylab="Observed risk",
       main=main_title)
  
  # PV
  lines(curve_pv$grid_p, curve_pv$grid_obs,
        lwd=2, col="#D95F02")
  
  # HARE
  lines(curve_hare$grid_p, curve_hare$grid_obs,
        lwd=2, col="#1B9E77")
  
  # ideal
  abline(0,1,lty=3,col="#7570B3")
  
  legend("topleft",
         legend=c("Cox-based","PV-based","HARE","Ideal"),
         col=c("#4D4D4D","#D95F02","#1B9E77","#7570B3"),
         lty=c(1,1,1,3),
         lwd=c(2,2,2,1),
         bty="n")
  

  # PANEL 2: HISTOGRAM

  par(mar = c(4,4,1,1))
  
  ymax <- max(h1$density, h2$density, h3$density)
  
  plot(NA, NA, xlim = x_zoom, ylim = c(0, ymax),
       xlab="Predicted risk", ylab="Density")
  
  for(i in seq_along(h1$density)){
    
    rect(h1$breaks[i],0,h1$breaks[i+1],h1$density[i],
         col=adjustcolor("#4D4D4D", alpha.f=0.25), border=NA)
    
    rect(h2$breaks[i],0,h2$breaks[i+1],h2$density[i],
         col=adjustcolor("#D95F02", alpha.f=0.20), border=NA)
    
    rect(h3$breaks[i],0,h3$breaks[i+1],h3$density[i],
         col=adjustcolor("#1B9E77", alpha.f=0.20), border=NA)
  }
  
  box()
  par(op)
}




plot_cal_with_hist(
  curve_cox = cal_base_cox,
  curve_pv  = cal_base,
  curve_hare = cal_base_hare,   # <-- devi averlo
  pred_cox  = prob_base,
  pred_pv   = prob_base,
  pred_hare = prob_base,
  main_title = "Calibration (Base Model)",
  x_zoom = c(0,0.1)
)


## without hare-------
plot_cal_with_hist <- function(curve_cox, curve_pv,
                               pred_cox, pred_pv,
                               main_title,
                               x_zoom = c(0,0.1),
                               breaks = seq(0,0.1,length.out=30)) {
  
  # filtro predictions
  v1 <- pred_cox[!is.na(pred_cox) & pred_cox >= x_zoom[1] & pred_cox <= x_zoom[2]]
  v2 <- pred_pv[!is.na(pred_pv) & pred_pv >= x_zoom[1] & pred_pv <= x_zoom[2]]
  
  h1 <- hist(v1, breaks = breaks, plot = FALSE)
  h2 <- hist(v2, breaks = breaks, plot = FALSE)
  
  op <- par(no.readonly = TRUE)
  layout(matrix(c(1,2), nrow = 2), heights = c(4,1))
  

  # PANEL 1: calibration

  par(mar = c(4,4,2,1))
  
  plot(curve_cox$grid_p, curve_cox$grid_obs,
       type="l", lwd=2, col="#4D4D4D",
       xlim=x_zoom, ylim=x_zoom,
       xlab="Predicted risk (3-year)",
       ylab="Observed risk",
       main=main_title)
  
  lines(curve_pv$grid_p, curve_pv$grid_obs,
        lwd=2, col="#D95F02")
  
  abline(0,1,lty=3,col="#7570B3")
  
  legend("topleft",
         legend=c("Cox-based","PV-based","Ideal"),
         col=c("#4D4D4D","#D95F02","#7570B3"),
         lty=c(1,1,3), lwd=c(2,2,1), bty="n")
  

  # PANEL 2: histogram

  par(mar = c(4,4,1,1))
  
  ymax <- max(h1$density, h2$density)
  
  plot(NA, NA, xlim = x_zoom, ylim = c(0, ymax),
       xlab="Predicted risk", ylab="Density")
  
  for(i in seq_along(h1$density)){
    rect(h1$breaks[i],0,h1$breaks[i+1],h1$density[i],
         col=adjustcolor("#4D4D4D", alpha.f=0.25), border=NA)
    
    rect(h2$breaks[i],0,h2$breaks[i+1],h2$density[i],
         col=adjustcolor("#D95F02", alpha.f=0.25), border=NA)
  }
  
  box()
  par(op)
}

plot_cal_with_hist(
  curve_cox = cal_biomm_cox,
  curve_pv  = cal_biomm,
  pred_cox  = prob_biom,
  pred_pv   = prob_biom,
  main_title = "Calibration (Biomarker Model)",
  x_zoom = c(0,0.1)
)








## fig 1 COX VS PV ------
ggplot() +
  # 1. Identity Line (Ideal calibration)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  # 2. Calibration Curves
  geom_line(data = df_base_cox, aes(x = pred, y = obs, color = "Cox-based Assessment"), size = 1.2) +
  geom_line(data = df_base_pv,  aes(x = pred, y = obs, color = "PV-based Assessment"),  size = 1.2) +
  # 3. Rug Plot (Shows data density on the X-axis)
  geom_rug(data = df_rug_base, aes(x = pred), sides = "b", alpha = 0.1, color = "grey30") +
  # 4. Styling
  scale_color_manual(values = c("Cox-based Assessment" = "#E41A1C", "PV-based Assessment" = "#377EB8")) +
  coord_cartesian(xlim = c(0, 0.10), ylim = c(0, 0.10)) +
  labs(
    title = " Comparison of Calibration Assessment Methods",
    subtitle = "Base Model (Cox-derived) judged by Traditional vs. PV logic",
    x = "Predicted 3-Year Risk",
    y = "Observed 3-Year Risk",
    color = "Assessment Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# Data for rugs
df_rugs = data.frame(
  pred = c(prob_cox_base, prob_cox_biom),
  model = rep(c("Base Model", "Biomarker Model"), each = length(prob_cox_base))
)



## fig 2 PV--------
ggplot() +
  # 1. Identity Line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  # 2. Calibration Curves (Both judged by PV)
  geom_line(data = df_base_pv, aes(x = pred, y = obs, color = "Base Model"), size = 1.2) +
  geom_line(data = df_biom_pv, aes(x = pred, y = obs, color = "Biomarker Model"), size = 1.2) +
  # 3. Rug Plots (Color-coded to match lines)
  geom_rug(data = df_rugs, aes(x = pred, color = model), sides = "b", alpha = 0.05) +
  # 4. Styling
  scale_color_manual(values = c("Base Model" = "#4DAF4A", "Biomarker Model" = "#984EA3")) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.1)) +
  labs(
    title = "Absolute Risk Calibration (Pseudo-Value Assessment)",
    subtitle = "Improvement in calibration-in-the-large with Biomarkers",
    x = "Predicted 3-Year Risk",
    y = "Observed 3-Year Risk (PV-Targeted)",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))



## fig 2 COX --------

ggplot() +
  # 1. Identity Line
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "grey50") +
  
  # 2. Calibration Curves (Cox-based)
  geom_line(data = df_base_cox,
            aes(x = pred, y = obs, color = "Base Model"),
            size = 1.2) +
  
  geom_line(data = df_biom_cox,
            aes(x = pred, y = obs, color = "Biomarker Model"),
            size = 1.2) +
  
  # 3. Rug Plot
  geom_rug(data = df_rugs,
           aes(x = pred, color = model),
           sides = "b",
           alpha = 0.05) +
  
  # 4. Styling
  scale_color_manual(values = c(
    "Base Model" = "#4DAF4A",
    "Biomarker Model" = "#984EA3"
  )) +
  
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 0.1)) +
  
  labs(
    title = "Absolute Risk Calibration (Cox-based Assessment)",
    subtitle = "Apparent calibration similarity between models",
    x = "Predicted 3-Year Risk",
    y = "Observed 3-Year Risk (Cox-Based)",
    color = "Model Type"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )


#ICI --------
#---------------------------------------------#
# Calibration error metrics: ICI, E50, E90
#---------------------------------------------#

cal_metrics = function(p_pred, obs_hat_i) {
  err = p_pred - obs_hat_i
  abs_err = abs(err)
  
  c(
    ICI  = mean(abs_err, na.rm = TRUE),
    E25  = as.numeric(quantile(abs_err, 0.25, na.rm = TRUE)),
    E30  = as.numeric(quantile(abs_err, 0.30, na.rm = TRUE)),
    E50  = median(abs_err, na.rm = TRUE),
    E75  = as.numeric(quantile(abs_err, 0.75, na.rm = TRUE)),
    E90  = as.numeric(quantile(abs_err, 0.90, na.rm = TRUE)),
    Emax = max(abs_err, na.rm = TRUE),
    RMSB = sqrt(mean(err^2, na.rm = TRUE))
  )
}

# --- A. Prepare the "Observed" judge ---

obs_base_cox = 1 - as.numeric(rms::survest(cal_base_cox$cal_fit,
                                           newdata = data.frame(eta = log(-log(1-prob_cox_base))),
                                           times = t0, conf.int = 0)$surv)

obs_biom_cox = 1 - as.numeric(rms::survest(cal_biomm_cox$cal_fit,
                                           newdata = data.frame(eta = log(-log(1-prob_cox_biom))),
                                           times = t0, conf.int = 0)$surv)


# --- B. Calculate Metrics for the 4 scenarios ---

m_base_pv = cal_metrics(p_pred = prob_cox_base, obs_hat_i = cal_base$obs_hat_i)
m_biom_pv = cal_metrics(p_pred = prob_cox_biom, obs_hat_i = cal_biomm$obs_hat_i)
m_base_traditional = cal_metrics(p_pred = prob_cox_base, obs_hat_i = obs_base_cox)
m_biom_traditional = cal_metrics(p_pred = prob_cox_biom, obs_hat_i = obs_biom_cox)

m_base_hare=cal_metrics(p_pred = prob_cox_base, obs_hat_i = cal_base_hare$obs_hat_i_cc)
m_bio_hare= cal_metrics(p_pred = prob_cox_biom, obs_hat_i = cal_bio_hare$obs_hat_i_cc) 


# --- C. Assemble the Table ---
metrics_tbl = rbind(
  "Base (PV Judge)"        = m_base_pv,
  "Biomarker (PV Judge)"   = m_biom_pv,
  "Base (Cox Judge)"       = m_base_traditional,
  "Biomarker (Cox Judge)"  = m_biom_traditional, 
  "Base (hare Judge)"        = m_base_hare,
  "Biomarker (hare Judge)"   = m_bio_hare
)

round(metrics_tbl, 4)


# 1. Compare the two nested PV models using a Wald Test
library(aod) # For the wald.test function

# Identify which coefficients are the biomarkers (e.g., the last two)
# Fit calibration curve at time t0; identity link returns risks on probability scale
total_coefs = length(coef(cal_biomm$cal_fit))
biomarker_indices = (total_coefs - 1):total_coefs
# Perform the Wald Test
pv_comparison = aod::wald.test(
  Sigma = vcov(cal_biomm$cal_fit),
  b     = coef(cal_biomm$cal_fit),
  Terms = biomarker_indices
)
# 2. Extract thcal_biomm # 2. Extract the p-value
p_value_pv = pv_comparison$result$chi2["P"]

# 3. Calculate the Delta Pseudo-R2
r2_base = 1 - (cal_base$cal_fit$deviance / cal_base$cal_fit$null.deviance)
r2_bio  = 1 - (cal_biomm$cal_fit$deviance  / cal_biomm$cal_fit$null.deviance)
delta_r2 = r2_bio - r2_base

delta_r2

