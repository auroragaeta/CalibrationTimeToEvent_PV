#CALIBRATION----------
#PREPROCESSING -------
rm(list = ls())
library(readxl)
library(survival)
library(pec)
library(haven)
library(survminer)
library(riskRegression)
library(gtsummary)
library(eventglm)
library(CalibrationCurves)
library(dplyr)
library(stringr)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(rms)  
library(survival)
library(riskRegression)
library(survival) 
library(ggplot2)
library(readr)


check <- read_csv("~/Documents/1. A TIMI/27_02_2026/training_mace_LOQ_final.csv") 
checkt  <-read_csv("~/Documents/1. A TIMI/27_02_2026/validation_mace_LOQ_final.csv")

checkt$event_3Y <- with(checkt, as.numeric(days2miistr <= 1095.75 & miistrfu == 1))
check$event_3Y <- with(check, as.numeric(days2miistr <= 1095.75 & miistrfu == 1))

df_model_train <- check %>%
  mutate(
    age_gr = case_when(
      age >= 75             ~ "[75-]",
      age >= 65 & age <= 74 ~ "[65-74]",
      TRUE                  ~ "[-65]"
    ),
    age_gr = factor(age_gr, levels = c("[-65]", "[65-74]", "[75-]")),  # ref: <65
    
    Male = if_else(Male == 1, 1L, 0L),  # ref: female
    
    BLwaist_sd = as.numeric(scale(BLwaist)),   # per 1-SD
    BLldlc_sd  = as.numeric(scale(BLldlc)),    # per 1-SD
    
    BLSBP_cate = case_when(
      BLSBP >= 160                   ~ ">=160",
      BLSBP >= 140 & BLSBP <= 159    ~ "140-159",
      TRUE                           ~ "<140"
    ),
    BLSBP_cate = factor(BLSBP_cate, levels = c("<140", "140-159", ">=160")), # ref: <140
    
    curr_smoke = if_else(curr_smoke %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    
    # ---- Diabetes 
    BLinsulinmed = if_else(BLinsulinmed %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    BLHBA1C      = as.numeric(scale(BLHBA1C)),  # per 1-SD
    
    # ---- Kidney 
    BLGF_cate = if_else(BLGFR < 60, 1L, 0L),  # ref: >=60
    
    UACR_gr = case_when(
      UACR >= 300            ~ ">=300",
      UACR >= 30 & UACR < 300 ~ "30-299",
      TRUE                   ~ "<30"
    ),
    UACR_gr = factor(UACR_gr, levels = c("<30", "30-299", ">=300")), # ref: <30
    
    hxcad   = if_else(str_to_upper(str_trim(as.character(hxcad)))   %in% c("1","YES","Y"), 1L, 0L),
    hxmi    = if_else(str_to_upper(str_trim(as.character(hxmi)))    %in% c("1","YES","Y"), 1L, 0L),
    hxpci   = if_else(str_to_upper(str_trim(as.character(hxpci)))   %in% c("1","YES","Y"), 1L, 0L),
    hxcabg  = if_else(str_to_upper(str_trim(as.character(hxcabg)))  %in% c("1","YES","Y"), 1L, 0L),
    hxistrk = if_else(str_to_upper(str_trim(as.character(hxistrk))) %in% c("1","YES","Y"), 1L, 0L),
    hxpad   = if_else(str_to_upper(str_trim(as.character(hxpad)))   %in% c("1","YES","Y"), 1L, 0L)
  )


df_model_test<- checkt %>%
  mutate(
    age_gr = case_when(
      age >= 75             ~ "[75-]",
      age >= 65 & age <= 74 ~ "[65-74]",
      TRUE                  ~ "[-65]"
    ),
    age_gr = factor(age_gr, levels = c("[-65]", "[65-74]", "[75-]")),  # ref: <65
    
    Male = if_else(Male == 1, 1L, 0L),  # ref: female
    
    BLSBP_cate = case_when(
      BLSBP >= 160                   ~ ">=160",
      BLSBP >= 140 & BLSBP <= 159    ~ "140-159",
      TRUE                           ~ "<140"
    ),
    BLSBP_cate = factor(BLSBP_cate, levels = c("<140", "140-159", ">=160")), # ref: <140
    
    curr_smoke = if_else(curr_smoke %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    
    BLinsulinmed = if_else(BLinsulinmed %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    
    BLGF_cate = if_else(BLGFR < 60, 1L, 0L),  # ref: >=60
    
    UACR_gr = case_when(
      UACR >= 300            ~ ">=300",
      UACR >= 30 & UACR < 300 ~ "30-299",
      TRUE                   ~ "<30"
    ),
    UACR_gr = factor(UACR_gr, levels = c("<30", "30-299", ">=300")), # ref: <30
    
    hxcad   = if_else(str_to_upper(str_trim(as.character(hxcad)))   %in% c("1","YES","Y"), 1L, 0L),
    hxmi    = if_else(str_to_upper(str_trim(as.character(hxmi)))    %in% c("1","YES","Y"), 1L, 0L),
    hxpci   = if_else(str_to_upper(str_trim(as.character(hxpci)))   %in% c("1","YES","Y"), 1L, 0L),
    hxcabg  = if_else(str_to_upper(str_trim(as.character(hxcabg)))  %in% c("1","YES","Y"), 1L, 0L),
    hxistrk = if_else(str_to_upper(str_trim(as.character(hxistrk))) %in% c("1","YES","Y"), 1L, 0L),
    hxpad   = if_else(str_to_upper(str_trim(as.character(hxpad)))   %in% c("1","YES","Y"), 1L, 0L)
  )




mu_waist <- mean(df_model_train$BLwaist,  na.rm = TRUE)
sd_waist <-  sd(df_model_train$BLwaist,   na.rm = TRUE)

mu_ldlc  <- mean(df_model_train$BLldlc,   na.rm = TRUE)
sd_ldlc  <-  sd(df_model_train$BLldlc,    na.rm = TRUE)

mu_hba1c <- mean(df_model_train$BLHBA1C,  na.rm = TRUE)
sd_hba1c <-  sd(df_model_train$BLHBA1C,   na.rm = TRUE)


df_model_test <- df_model_test %>%
  mutate(
    BLwaist_sd = (BLwaist  - mu_waist)/sd_waist,
    BLldlc_sd  = (BLldlc   - mu_ldlc) /sd_ldlc,
    BLHBA1C    = (BLHBA1C  - mu_hba1c)/sd_hba1c
  )

#COX COEFFICIENT ON THE OVERALL
cox_fit <- coxph(Surv(days2miistr,miistrfu) ~
                   age_gr
                 + Male
                 + BLwaist_sd
                 + BLldlc_sd
                 + BLSBP_cate
                 + curr_smoke
                 + BLinsulinmed
                 + BLHBA1C
                 + BLGF_cate
                 + UACR_gr
                 + hxcad
                 + hxmi
                 + hxpci+hxcabg+hxistrk+hxpad,
                 data = df_model_train, 
                 ties = 'breslow', x=TRUE,y=TRUE)


cox_fit %>% tbl_regression(exponentiate=T)

check<-df_model_train[complete.cases(df_model_train[,c("age_gr",
                                                       "Male" ,
                                                       "BLwaist_sd" , 
                                                       "BLldlc_sd", 
                                                       "BLSBP_cate" ,
                                                       "curr_smoke" ,
                                                       "BLinsulinmed" ,
                                                       "BLHBA1C"  ,
                                                       "BLGF_cate",
                                                       "UACR_gr" ,
                                                       "hxcad",
                                                       "hxmi" ,
                                                       "hxpci" ,
                                                       "hxcabg",
                                                       "hxistrk",
                                                       "hxpad" ,
                                                       "log_nt",
                                                       "log_trop")]),]

checkt<-df_model_test[complete.cases(df_model_test[,c("age_gr",
                                                      "Male" ,
                                                      "BLwaist_sd" , 
                                                      "BLldlc_sd", 
                                                      "BLSBP_cate" ,
                                                      "curr_smoke" ,
                                                      "BLinsulinmed" ,
                                                      "BLHBA1C"  ,
                                                      "BLGF_cate",
                                                      "UACR_gr" ,
                                                      "hxcad",
                                                      "hxmi" ,
                                                      "hxpci" ,
                                                      "hxcabg",
                                                      "hxistrk",
                                                      "hxpad" ,
                                                      "log_nt",
                                                      "log_trop")]),]

#checkt= checkt[-which(checkt$usubjid=='11864011026'),]

#checkt_validation=checkt[1:round(nrow(checkt)*0.30),]
#checkt=checkt[-c(1:round(nrow(checkt)*0.30)),]

df_model_train <- as.data.frame(df_model_train)
checkt<- as.data.frame(checkt)
#checkt_validation<- as.data.frame(checkt_validation)

pred_vars <- c(
  "age_gr","Male","BLwaist_sd","BLldlc_sd","BLSBP_cate","curr_smoke",
  "BLinsulinmed","BLHBA1C","BLGF_cate","UACR_gr",
  "hxcad","hxmi","hxpci","hxcabg","hxistrk","hxpad" )

# aligns the levels of the factors between train and checkt
for (v in pred_vars) {
  if (is.factor(df_model_train[[v]])) {
    checkt[[v]] <- factor(checkt[[v]], levels = levels(df_model_train[[v]]))
    #checkt_validation[[v]] <- factor(checkt_validation[[v]], levels = levels(df_model_train[[v]]))
  }
}


df_model_train$ln_prod<-as.numeric(predict(cox_fit, newdata=df_model_train, type="lp"))
checkt$ln_prod<-as.numeric(predict(cox_fit, newdata=checkt, type="lp"))
#checkt_validation$ln_prod<-as.numeric(predict(cox_fit, newdata=checkt_validation, type="lp"))

checkt$ln_prod_rank<-1-predictSurvProb(cox_fit,
                                       newdata=checkt,
                                       times=1095.75)

df_model_train$ln_prod_rank<-1-predictSurvProb(cox_fit,
                                               newdata=df_model_train,
                                               times=1095.75)


df_model_train$ln_prod_rank=as.numeric(df_model_train$ln_prod_rank)
checkt$ln_prod_rank=as.numeric(checkt$ln_prod_rank)


df_model_train$ln_prod_rank=as.numeric(df_model_train$ln_prod_rank)
checkt$ln_prod_rank=as.numeric(checkt$ln_prod_rank)

mod_biom = coxph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad + log_nt + log_trop,
  data = df_model_train,
  x = TRUE,
  y = TRUE
)


mod_base = coxph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad,
  data = df_model_train,
  x = TRUE,         
  y = TRUE            
)



dd = datadist(checkt); options(datadist="dd")  
cloglog = function(p) log(-log(1 - p))



##MAKE CALIBRATION PSEUDOVALUES  ---------------------------------------------------
make_calibration = function(p_pred, data, t0, knots = 3, eps = 1e-6) {
  # p_pred = p_base
  # data = checkt
  # t0 = 365.25*3
  # knots=5
  
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = cloglog(p_pred)
  
  dcal = data.frame(
    days2miistr = data$days2miistr,
    miistrfu   = data$miistrfu,
    eta        = eta
  )
  
  dd  = datadist(dcal)
  old = options(datadist = "dd")
  on.exit(options(old), add = TRUE)
  
  # Fit calibration curve at time t0; identity link returns risks on probability scale
  cal_fit = eventglm::cumincglm(
    Surv(days2miistr, miistrfu) ~ rms::rcs(eta, knots),
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
make_calibration_cox = function(p_pred, data, t0, knots = 3, eps = 1e-6) {
  # 1. Transform probability to the Cox linear scale (cloglog)
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = log(-log(1 - p_pred)) # This is the cloglog/Cox link
  
  dcal = data.frame(
    time = data$days2miistr,
    status = data$miistrfu,
    eta = eta
  )
  
  # 2. Fit a Cox model with a flexible spline for the score
  # This is the "Traditional" recalibration/smoothing step
  cal_fit_cox = rms::cph(
    Surv(time, status) ~ rms::rcs(eta, knots),
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


## MODEL SIMPLE AND COMPLETE IN TRAINING --------------------------------------
mod_biom = cph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad + log_nt + log_trop,
  data = check,
  x = TRUE,
  y = TRUE, surv=TRUE
)

mod_base = cph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad,
  data = check,
  x = TRUE,
  y = TRUE, surv=TRUE
)


## Predicted probability from COX in validation #########################

checkt$prob_cox_base<-1-predictSurvProb(mod_base,
                                        newdata=checkt,
                                        times=1095.75)

checkt$prob_cox_biom<-1-predictSurvProb(mod_biom,
                                        newdata=checkt,
                                        times=1095.75)


prob_cox_base=as.numeric(checkt$prob_cox_base)
prob_cox_biom=as.numeric(checkt$prob_cox_biom)


#-------------------------------#
# 0) Inputs and data
#-------------------------------#

t0  = 365.25 * 3      # 3-year horizon in days
eps = 1e-6            # numerical safety constant for probabilities

train = check
test  = checkt

cal_base = make_calibration(prob_cox_base, data = test, t0 = t0, knots = 3, eps = eps)

cal_bio  = make_calibration(prob_cox_biom,  data = test, t0 = t0, knots = 3, eps = eps)


cal_base_cox = make_calibration_cox(prob_cox_base, data = test, t0 = t0, knots = 3, eps = eps)

cal_biom_cox  = make_calibration_cox(prob_cox_biom,  data = test, t0 = t0, knots = 3, eps = eps)


# Note: val.surv stores its smoothed curve in $recal

# Extract Traditional (Cox-based) curve for Base Model
df_base_cox = data.frame(
  pred = cal_base_cox$grid_p, 
  obs  = cal_base_cox$grid_obs, 
  method = "Traditional (Cox)"
)

# Extract Traditional (Cox-based) curve for Biomarker Model
df_biom_cox = data.frame(
  pred = cal_biom_cox$grid_p, 
  obs  = cal_biom_cox$grid_obs, 
  method = "Traditional (Cox)"
)


# Extract PV-based curves
df_base_pv = data.frame(
  pred = cal_base$grid_p, 
  obs = cal_base$grid_obs, 
  method = "New (PV)")

df_biom_pv = data.frame(
  pred = cal_bio$grid_p,  
  obs = cal_bio$grid_obs,  
  method = "New (PV)")

# Create a small dataframe for the rug (individual patient predictions)
df_rug_base = data.frame(pred = prob_cox_base)

##fig 1------


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
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
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



## fig 2 --------
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
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  labs(
    title = "Absolute Risk Calibration (Pseudo-Value Assessment)",
    subtitle = "Improvement in calibration-in-the-large with Biomarkers",
    x = "Predicted 3-Year Risk",
    y = "Observed 3-Year Risk (PV-Targeted)",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

#ICI --------

## Should not be calculated on the grid-----
# Calculate ICI for PV-based assessment
ici_base_pv = mean(abs(cal_base$grid_p - cal_base$grid_obs))
ici_biom_pv = mean(abs(cal_bio$grid_p  - cal_bio$grid_obs))

# Calculate ICI for Traditional (Cox) assessment
ici_base_cox = mean(abs(cal_base_cox$grid_p - cal_base_cox$grid_obs))
ici_biom_cox = mean(abs(cal_biom_cox$grid_p - cal_biom_cox$grid_obs))

# Print results
cat("--- PV-BASED ASSESSMENT by grid ---\n",
    "Base Model ICI: ", round(ici_base_pv, 4), "\n",
    "Biomarker Model ICI: ", round(ici_biom_pv, 4), "\n\n",
    "--- COX-BASED ASSESSMENT by grid---\n",
    "Base Model ICI: ", round(ici_base_cox, 4), "\n",
    "Biomarker Model ICI: ", round(ici_biom_cox, 4), "\n")


#---------------------------------------------#
# CORRECT alibration error metrics: ICI, E50, E90
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

# --- A. Prepare the "Observed" judge from the Cox engine ---
# We use the individual predictions (eta) to get individual smoothed observations
# for the Cox-based judge.
obs_base_cox = 1 - as.numeric(rms::survest(cal_base_cox$cal_fit, 
                                           newdata = data.frame(eta = log(-log(1-prob_cox_base))), 
                                           times = t0, conf.int = 0)$surv)

obs_biom_cox = 1 - as.numeric(rms::survest(cal_biom_cox$cal_fit, 
                                           newdata = data.frame(eta = log(-log(1-prob_cox_biom))), 
                                           times = t0, conf.int = 0)$surv)

# --- B. Calculate Metrics for the 4 scenarios ---

# 1. Base Model assessed via PV (The Honest Way)
m_base_pv = cal_metrics(p_pred = prob_cox_base, obs_hat_i = cal_base$obs_hat_i)

# 2. Biomarker Model assessed via PV (The Modern Way)
m_biom_pv = cal_metrics(p_pred = prob_cox_biom, obs_hat_i = cal_bio$obs_hat_i)

# 3. Base Model assessed via COX (The Traditional Way)
m_base_traditional = cal_metrics(p_pred = prob_cox_base, obs_hat_i = obs_base_cox)

# 4. Biomarker Model assessed via COX (The "Biased" Way)
m_biom_traditional = cal_metrics(p_pred = prob_cox_biom, obs_hat_i = obs_biom_cox)

# --- C. Assemble the Table ---
metrics_tbl = rbind(
  "Base (PV Judge)"        = m_base_pv,
  "Biomarker (PV Judge)"   = m_biom_pv,
  "Base (Cox Judge)"       = m_base_traditional,
  "Biomarker (Cox Judge)"  = m_biom_traditional
)

round(metrics_tbl, 4)


# 1. Compare the two nested PV models using a Wald Test
# This tests if the coefficients for log_nt and log_trop are significantly different from zero
# in the context of the 3-year absolute risk.

library(aod) # For the wald.test function

# Identify which coefficients are the biomarkers (e.g., the last two)
# Let's say log_nt is the 17th and log_trop is the 18th

# Fit calibration curve at time t0; identity link returns risks on probability scale

total_coefs = length(coef(cal_bio$cal_fit))
biomarker_indices = (total_coefs - 1):total_coefs

# Perform the Wald Test
pv_comparison = aod::wald.test(
  Sigma = vcov(cal_bio$cal_fit), 
  b     = coef(cal_bio$cal_fit), 
  Terms = biomarker_indices
)

# 2. Extract the p-value
p_value_pv = pv_comparison$result$chi2["P"]

# 3. Calculate the Delta Pseudo-R2
r2_base = 1 - (cal_base$cal_fit$deviance / cal_base$cal_fit$null.deviance)
r2_bio  = 1 - (cal_bio$cal_fit$deviance  / cal_bio$cal_fit$null.deviance)
delta_r2 = r2_bio - r2_base

delta_r2 



#RECALIBRATION------
##PREPROCESSING ---------
rm(list = ls())

library(readxl)
library(survival)
library(pec)
library(haven)
library(survminer)
library(riskRegression)
library(gtsummary)
library(eventglm)
library(CalibrationCurves)
library(dplyr)
library(stringr)
library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(rms)  
library(survival)
library(riskRegression)
library(survival) 
library(ggplot2)
library(readr)


check <- read_csv("~/Documents/1. A TIMI/27_02_2026/training_mace_LOQ_final.csv") 
checkt  <-read_csv("~/Documents/1. A TIMI/27_02_2026/validation_mace_LOQ_final.csv")

checkt$event_3Y <- with(checkt, as.numeric(days2miistr <= 1095.75 & miistrfu == 1))
check$event_3Y <- with(check, as.numeric(days2miistr <= 1095.75 & miistrfu == 1))

df_model_train <- check %>%
  mutate(
    age_gr = case_when(
      age >= 75             ~ "[75-]",
      age >= 65 & age <= 74 ~ "[65-74]",
      TRUE                  ~ "[-65]"
    ),
    age_gr = factor(age_gr, levels = c("[-65]", "[65-74]", "[75-]")),  # ref: <65
    
    Male = if_else(Male == 1, 1L, 0L),  # ref: female
    
    BLwaist_sd = as.numeric(scale(BLwaist)),   # per 1-SD
    BLldlc_sd  = as.numeric(scale(BLldlc)),    # per 1-SD
    
    BLSBP_cate = case_when(
      BLSBP >= 160                   ~ ">=160",
      BLSBP >= 140 & BLSBP <= 159    ~ "140-159",
      TRUE                           ~ "<140"
    ),
    BLSBP_cate = factor(BLSBP_cate, levels = c("<140", "140-159", ">=160")), # ref: <140
    
    curr_smoke = if_else(curr_smoke %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    
    # ---- Diabetes 
    BLinsulinmed = if_else(BLinsulinmed %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    BLHBA1C      = as.numeric(scale(BLHBA1C)),  # per 1-SD
    
    # ---- Kidney 
    BLGF_cate = if_else(BLGFR < 60, 1L, 0L),  # ref: >=60
    
    UACR_gr = case_when(
      UACR >= 300            ~ ">=300",
      UACR >= 30 & UACR < 300 ~ "30-299",
      TRUE                   ~ "<30"
    ),
    UACR_gr = factor(UACR_gr, levels = c("<30", "30-299", ">=300")), # ref: <30
    
    hxcad   = if_else(str_to_upper(str_trim(as.character(hxcad)))   %in% c("1","YES","Y"), 1L, 0L),
    hxmi    = if_else(str_to_upper(str_trim(as.character(hxmi)))    %in% c("1","YES","Y"), 1L, 0L),
    hxpci   = if_else(str_to_upper(str_trim(as.character(hxpci)))   %in% c("1","YES","Y"), 1L, 0L),
    hxcabg  = if_else(str_to_upper(str_trim(as.character(hxcabg)))  %in% c("1","YES","Y"), 1L, 0L),
    hxistrk = if_else(str_to_upper(str_trim(as.character(hxistrk))) %in% c("1","YES","Y"), 1L, 0L),
    hxpad   = if_else(str_to_upper(str_trim(as.character(hxpad)))   %in% c("1","YES","Y"), 1L, 0L)
  )


df_model_test<- checkt %>%
  mutate(
    age_gr = case_when(
      age >= 75             ~ "[75-]",
      age >= 65 & age <= 74 ~ "[65-74]",
      TRUE                  ~ "[-65]"
    ),
    age_gr = factor(age_gr, levels = c("[-65]", "[65-74]", "[75-]")),  # ref: <65
    
    Male = if_else(Male == 1, 1L, 0L),  # ref: female
    
    BLSBP_cate = case_when(
      BLSBP >= 160                   ~ ">=160",
      BLSBP >= 140 & BLSBP <= 159    ~ "140-159",
      TRUE                           ~ "<140"
    ),
    BLSBP_cate = factor(BLSBP_cate, levels = c("<140", "140-159", ">=160")), # ref: <140
    
    curr_smoke = if_else(curr_smoke %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    
    BLinsulinmed = if_else(BLinsulinmed %in% c(1, "1", "Yes", "Y", "YES", "yes"), 1L, 0L), # ref: no
    
    BLGF_cate = if_else(BLGFR < 60, 1L, 0L),  # ref: >=60
    
    UACR_gr = case_when(
      UACR >= 300            ~ ">=300",
      UACR >= 30 & UACR < 300 ~ "30-299",
      TRUE                   ~ "<30"
    ),
    UACR_gr = factor(UACR_gr, levels = c("<30", "30-299", ">=300")), # ref: <30
    
    hxcad   = if_else(str_to_upper(str_trim(as.character(hxcad)))   %in% c("1","YES","Y"), 1L, 0L),
    hxmi    = if_else(str_to_upper(str_trim(as.character(hxmi)))    %in% c("1","YES","Y"), 1L, 0L),
    hxpci   = if_else(str_to_upper(str_trim(as.character(hxpci)))   %in% c("1","YES","Y"), 1L, 0L),
    hxcabg  = if_else(str_to_upper(str_trim(as.character(hxcabg)))  %in% c("1","YES","Y"), 1L, 0L),
    hxistrk = if_else(str_to_upper(str_trim(as.character(hxistrk))) %in% c("1","YES","Y"), 1L, 0L),
    hxpad   = if_else(str_to_upper(str_trim(as.character(hxpad)))   %in% c("1","YES","Y"), 1L, 0L)
  )




mu_waist <- mean(df_model_train$BLwaist,  na.rm = TRUE)
sd_waist <-  sd(df_model_train$BLwaist,   na.rm = TRUE)

mu_ldlc  <- mean(df_model_train$BLldlc,   na.rm = TRUE)
sd_ldlc  <-  sd(df_model_train$BLldlc,    na.rm = TRUE)

mu_hba1c <- mean(df_model_train$BLHBA1C,  na.rm = TRUE)
sd_hba1c <-  sd(df_model_train$BLHBA1C,   na.rm = TRUE)


df_model_test <- df_model_test %>%
  mutate(
    BLwaist_sd = (BLwaist  - mu_waist)/sd_waist,
    BLldlc_sd  = (BLldlc   - mu_ldlc) /sd_ldlc,
    BLHBA1C    = (BLHBA1C  - mu_hba1c)/sd_hba1c
  )

#COX COEFFICIENT ON THE OVERALL
cox_fit <- coxph(Surv(days2miistr,miistrfu) ~
                   age_gr
                 + Male
                 + BLwaist_sd
                 + BLldlc_sd
                 + BLSBP_cate
                 + curr_smoke
                 + BLinsulinmed
                 + BLHBA1C
                 + BLGF_cate
                 + UACR_gr
                 + hxcad
                 + hxmi
                 + hxpci+hxcabg+hxistrk+hxpad,
                 data = df_model_train, 
                 ties = 'breslow', x=TRUE,y=TRUE)


cox_fit %>% tbl_regression(exponentiate=T)


check<-df_model_train[complete.cases(df_model_train[,c("age_gr",
                                                       "Male" ,
                                                       "BLwaist_sd" , 
                                                       "BLldlc_sd", 
                                                       "BLSBP_cate" ,
                                                       "curr_smoke" ,
                                                       "BLinsulinmed" ,
                                                       "BLHBA1C"  ,
                                                       "BLGF_cate",
                                                       "UACR_gr" ,
                                                       "hxcad",
                                                       "hxmi" ,
                                                       "hxpci" ,
                                                       "hxcabg",
                                                       "hxistrk",
                                                       "hxpad" ,
                                                       "log_nt",
                                                       "log_trop")]),]

checkt<-df_model_test[complete.cases(df_model_test[,c("age_gr",
                                                      "Male" ,
                                                      "BLwaist_sd" , 
                                                      "BLldlc_sd", 
                                                      "BLSBP_cate" ,
                                                      "curr_smoke" ,
                                                      "BLinsulinmed" ,
                                                      "BLHBA1C"  ,
                                                      "BLGF_cate",
                                                      "UACR_gr" ,
                                                      "hxcad",
                                                      "hxmi" ,
                                                      "hxpci" ,
                                                      "hxcabg",
                                                      "hxistrk",
                                                      "hxpad" ,
                                                      "log_nt",
                                                      "log_trop")]),]

#PLS EXCLUDE 11864011026
#checkt= checkt[-which(checkt$usubjid=='11864011026'),]


#checkt_validation=checkt[1:round(nrow(checkt)*0.30),]
#checkt=checkt[-c(1:round(nrow(checkt)*0.30)),]

df_model_train <- as.data.frame(df_model_train)
checkt<- as.data.frame(checkt)
#checkt_validation<- as.data.frame(checkt_validation)

pred_vars <- c(
  "age_gr","Male","BLwaist_sd","BLldlc_sd","BLSBP_cate","curr_smoke",
  "BLinsulinmed","BLHBA1C","BLGF_cate","UACR_gr",
  "hxcad","hxmi","hxpci","hxcabg","hxistrk","hxpad" )

# aligns the levels of the factors between train and checkt
for (v in pred_vars) {
  if (is.factor(df_model_train[[v]])) {
    checkt[[v]] <- factor(checkt[[v]], levels = levels(df_model_train[[v]]))
    #checkt_validation[[v]] <- factor(checkt_validation[[v]], levels = levels(df_model_train[[v]]))
  }
}


df_model_train$ln_prod<-as.numeric(predict(cox_fit, newdata=df_model_train, type="lp"))

checkt$ln_prod<-as.numeric(predict(cox_fit, newdata=checkt, type="lp"))
#checkt_validation$ln_prod<-as.numeric(predict(cox_fit, newdata=checkt_validation, type="lp"))

checkt$ln_prod_rank<-1-predictSurvProb(cox_fit,
                                       newdata=checkt,
                                       times=1095.75)

df_model_train$ln_prod_rank<-1-predictSurvProb(cox_fit,
                                               newdata=df_model_train,
                                               times=1095.75)


df_model_train$ln_prod_rank=as.numeric(df_model_train$ln_prod_rank)
checkt$ln_prod_rank=as.numeric(checkt$ln_prod_rank)


df_model_train$ln_prod_rank=as.numeric(df_model_train$ln_prod_rank)
checkt$ln_prod_rank=as.numeric(checkt$ln_prod_rank)

mod_biom = coxph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad + log_nt + log_trop,
  data = df_model_train,
  x = TRUE,
  y = TRUE
)


mod_base = coxph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad,
  data = df_model_train,
  x = TRUE,         
  y = TRUE            
)

dd = datadist(checkt); options(datadist="dd")
cloglog = function(p) log(-log(1 - p))

##CALIBRATION PSEUDOVALUES  ---------------------------------------------------
make_calibration = function(p_pred, data, t0, knots = 3, eps = 1e-6) {
  # p_pred = p_base
  # data = checkt
  # t0 = 365.25*3
  # knots=5
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = cloglog(p_pred)
  
  dcal = data.frame(
    days2miistr = data$days2miistr,
    miistrfu   = data$miistrfu,
    eta        = eta
  )
  
  dd  = datadist(dcal)
  old = options(datadist = "dd")
  on.exit(options(old), add = TRUE)
  
  # Fit calibration curve at time t0; identity link returns risks on probability scale
  cal_fit = eventglm::cumincglm(
    Surv(days2miistr, miistrfu) ~ rms::rcs(eta, knots),
    data = dcal,
    time = t0,
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


##CALIBRATION COX  ---------------------------------------------------

make_calibration_cox = function(p_pred, data, t0, knots = 3, eps = 1e-6) {
  # 1. Transform probability to the Cox linear scale (cloglog)
  p_pred = pmin(pmax(p_pred, eps), 1 - eps)
  eta = log(-log(1 - p_pred)) # This is the cloglog/Cox link
  
  dcal = data.frame(
    time = data$days2miistr,
    status = data$miistrfu,
    eta = eta
  )
  
  # 2. Fit a Cox model with a flexible spline for the score
  # This is the "Traditional" recalibration/smoothing step
  cal_fit_cox = rms::cph(
    Surv(time, status) ~ rms::rcs(eta, knots),
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



## MODEL SIMPLE AND COMPLETE IN TRAINING

mod_biom = cph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad + log_nt + log_trop,
  data = check,
  x = TRUE,
  y = TRUE, surv=TRUE
)


mod_base = cph(
  Surv(days2miistr, miistrfu) ~
    age_gr + Male + BLwaist_sd + BLldlc_sd + BLSBP_cate + curr_smoke +
    BLinsulinmed + BLHBA1C + BLGF_cate + UACR_gr +
    hxcad + hxmi + hxpci + hxcabg + hxistrk + hxpad,
  data = check,
  x = TRUE,
  y = TRUE, surv=TRUE
)


checkt$ln_prod_base<-as.numeric(predict(mod_base, newdata=checkt, type="lp"))
checkt$ln_prod_biom<-as.numeric(predict(mod_biom, newdata=checkt, type="lp"))

mod_biom = cph(
  Surv(days2miistr, miistrfu) ~ ln_prod_biom,
  data = checkt,
  x = TRUE,
  y = TRUE, surv=TRUE
)


mod_base = cph(
  Surv(days2miistr, miistrfu) ~ln_prod_base,
  data = checkt,
  x = TRUE,
  y = TRUE, surv=TRUE
)



# Predicted probability from COX in validation

checkt$prob_cox_base<-1-predictSurvProb(mod_base,
                                        newdata=checkt,
                                        times=1095.75)

checkt$prob_cox_biom<-1-predictSurvProb(mod_biom,
                                        newdata=checkt,
                                        times=1095.75)


prob_cox_base=as.numeric(checkt$prob_cox_base)
prob_cox_biom=as.numeric(checkt$prob_cox_biom)


#-------------------------------#
# 0) Inputs and data
#-------------------------------#

t0  = 365.25 * 3      # 3-year horizon in days
eps = 1e-6            # numerical safety constant for probabilities

train = check
test  = checkt

cal_base = make_calibration(prob_cox_base, data = test, t0 = t0, knots = 3, eps = eps)
cal_bio  = make_calibration(prob_cox_biom,  data = test, t0 = t0, knots = 3, eps = eps)

cal_base_cox = make_calibration_cox(prob_cox_base, data = test, t0 = t0, knots = 3, eps = eps)
cal_biom_cox  = make_calibration_cox(prob_cox_biom,  data = test, t0 = t0, knots = 3, eps = eps)


#-------------------------------#

# Note: val.surv stores its smoothed curve in $recal
# Extract Traditional (Cox-based) curve for Base Model
df_base_cox = data.frame(
  pred = cal_base_cox$grid_p, 
  obs  = cal_base_cox$grid_obs, 
  method = "Traditional (Cox)"
)

# Extract Traditional (Cox-based) curve for Biomarker Model
df_biom_cox = data.frame(
  pred = cal_biom_cox$grid_p, 
  obs  = cal_biom_cox$grid_obs, 
  method = "Traditional (Cox)"
)

# Extract PV-based curves
df_base_pv = data.frame(pred = cal_base$grid_p, obs = cal_base$grid_obs, method = "New (PV)")
df_biom_pv = data.frame(pred = cal_bio$grid_p,  obs = cal_bio$grid_obs,  method = "New (PV)")

# Create a small dataframe for the rug (individual patient predictions)
df_rug_base = data.frame(pred = prob_cox_base)

##fig 1-----
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
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
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




## fig 2 ------
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
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  labs(
    title = "Absolute Risk Calibration (Pseudo-Value Assessment)",
    subtitle = "Improvement in calibration-in-the-large with Biomarkers",
    x = "Predicted 3-Year Risk",
    y = "Observed 3-Year Risk (PV-Targeted)",
    color = "Model Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

## ICI --------
## Should not be calculated from the grid -----
# Calculate ICI for PV-based assessment
ici_base_pv = mean(abs(cal_base$grid_p - cal_base$grid_obs))
ici_biom_pv = mean(abs(cal_bio$grid_p  - cal_bio$grid_obs))

# Calculate ICI for Traditional (Cox) assessment
ici_base_cox = mean(abs(cal_base_cox$grid_p - cal_base_cox$grid_obs))
ici_biom_cox = mean(abs(cal_biom_cox$grid_p - cal_biom_cox$grid_obs))

# Print results
cat("--- PV-BASED ASSESSMENT grid  ---\n",
    "Base Model ICI: ", round(ici_base_pv, 4), "\n",
    "Biomarker Model ICI: ", round(ici_biom_pv, 4), "\n\n",
    "--- COX-BASED ASSESSMENT grid ---\n",
    "Base Model ICI: ", round(ici_base_cox, 4), "\n",
    "Biomarker Model ICI: ", round(ici_biom_cox, 4), "\n")


#---------------------------------------------#
# 4) Calibration error metrics: ICI, E50, E90
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

# --- A. Prepare the "Observed" judge from the Cox engine ---
# We use the individual predictions (eta) to get individual smoothed observations
# for the Cox-based judge.
obs_base_cox = 1 - as.numeric(rms::survest(cal_base_cox$cal_fit, 
                                           newdata = data.frame(eta = log(-log(1-prob_cox_base))), 
                                           times = t0, conf.int = 0)$surv)

obs_biom_cox = 1 - as.numeric(rms::survest(cal_biom_cox$cal_fit, 
                                           newdata = data.frame(eta = log(-log(1-prob_cox_biom))), 
                                           times = t0, conf.int = 0)$surv)

# --- B. Calculate Metrics for the 4 scenarios ---

# 1. Base Model assessed via PV (The Honest Way)
m_base_pv = cal_metrics(p_pred = prob_cox_base, obs_hat_i = cal_base$obs_hat_i)

# 2. Biomarker Model assessed via PV (The Modern Way)
m_biom_pv = cal_metrics(p_pred = prob_cox_biom, obs_hat_i = cal_bio$obs_hat_i)

# 3. Base Model assessed via COX (The Traditional Way)
m_base_traditional = cal_metrics(p_pred = prob_cox_base, obs_hat_i = obs_base_cox)

# 4. Biomarker Model assessed via COX (The "Biased" Way)
m_biom_traditional = cal_metrics(p_pred = prob_cox_biom, obs_hat_i = obs_biom_cox)

# --- C. Assemble the Table ---
metrics_tbl = rbind(
  "Base (PV Judge)"        = m_base_pv,
  "Biomarker (PV Judge)"   = m_biom_pv,
  "Base (Cox Judge)"       = m_base_traditional,
  "Biomarker (Cox Judge)"  = m_biom_traditional
)

round(metrics_tbl, 4)




# 1. Compare the two nested PV models using a Wald Test
# This tests if the coefficients for log_nt and log_trop are significantly different from zero
# in the context of the 3-year absolute risk.

library(aod) # For the wald.test function

# Identify which coefficients are the biomarkers (e.g., the last two)
# Let's say log_nt is the 17th and log_trop is the 18th

# Fit calibration curve at time t0; identity link returns risks on probability scale

total_coefs = length(coef(cal_bio$cal_fit))
biomarker_indices = (total_coefs - 1):total_coefs

# Perform the Wald Test
pv_comparison = aod::wald.test(
  Sigma = vcov(cal_bio$cal_fit), 
  b     = coef(cal_bio$cal_fit), 
  Terms = biomarker_indices
)

# 2. Extract the p-value
p_value_pv = pv_comparison$result$chi2["P"]

# 3. Calculate the Delta Pseudo-R2
r2_base = 1 - (cal_base$cal_fit$deviance / cal_base$cal_fit$null.deviance)
r2_bio  = 1 - (cal_bio$cal_fit$deviance  / cal_bio$cal_fit$null.deviance)
delta_r2 = r2_bio - r2_base

delta_r2



