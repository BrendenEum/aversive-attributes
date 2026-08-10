# Load necessary libraries
rm(list=ls())
library(lme4)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# Load data for prior beliefs about effect sizes
paths    <- yaml::yaml.load_file(here::here("config/paths.yaml"))
.datadir <- file.path(paths$root, paths$data$processed)
load(file.path(.datadir, "df_t.RData"))
df_t$choice = as.numeric(df_t$ChooseLR=="L")
df_t$DwellLT_Binary <- df_t$DwellLT < mean(df_t$DwellLT)
df_t$DwellLB_Binary <- df_t$DwellLB < mean(df_t$DwellLB)
df_t$DwellRT_Binary <- df_t$DwellRT < mean(df_t$DwellRT)
df_t$DwellRB_Binary <- df_t$DwellRB < mean(df_t$DwellRB)

# Function to simulate one dataset with subject-level randomness in predictors
simulate_dataset <- function(n_subj, n_trials, 
                             b0_m, b0_sd, bV_m, bV_sd, bInt_m, bInt_sd, 
                             g0_sd, gV_sd, gInt_sd) {
  
  # Total observations
  N <- n_subj * n_trials
  
  # Average Levels
  beta0 <- rnorm(1, b0_m, b0_sd) %>% rep(N)
  betaVL <- rnorm(1, bV_m, bV_sd) %>% rep(N)
  betaVR <- rnorm(1, -bV_m, bV_sd) %>% rep(N)
  betaIntSig <- rnorm(1, bInt_m, bInt_sd) %>% rep(N)
  betaIntInsig <- rnorm(1, 0, bInt_sd) %>% rep(N)
  
  # Subject-Level deviations from average
  gam0 <- rnorm(n_subj, 0, g0_sd) %>% rep(each=n_trials)
  gamV <- rnorm(n_subj, 0, gV_sd) %>% rep(each=n_trials)
  gamInt <- rnorm(n_subj, 0, gInt_sd) %>% rep(each=n_trials)
  
  # Data
  xVal <- sample(seq(5,45,10), N*4, replace=T) %>% matrix(nrow=N, ncol=4)
  xInt  <- as.numeric(xVal==25) * rnorm(N*20, 120, 30)
  xVal <- xVal %>% matrix(nrow=N, ncol=4)
  xInt <- xInt %>% matrix(nrow=N, ncol=20)
  
  # Generate subject IDs
  subject <- rep(1:n_subj, each = n_trials)
  
  # Choice Predictions
  linpred <- beta0 + (beta0+gam0) +
    betaVL*xVal[,1] + (betaVL+gamV)*xVal[,1] +
    betaVR*xVal[,2] + (betaVR+gamV)*xVal[,2] +
    betaIntSig*xInt[,1] + (betaIntSig+gamInt)*xInt[,1] +
    betaIntInsig*xInt[,2] + (betaIntInsig+gamInt)*xInt[,2]
  
  # Convert to probability using logistic link
  prob <- 1 / (1 + exp(-linpred))
  
  # Generate binary outcome
  y <- rbinom(N, size = 1, prob = prob)
  
  # Return dataframe
  df <- data.frame(subject=subject, choice=y, 
                   xVal1=xVal[,1], xVal2=xVal[,2], xInt1=xInt[,1], xInt2=xInt[,2])
  return(df)
}

# Function to run one simulation iteration
run_simulation <- function(n_subj, n_trials, 
                           b0_m, b0_sd, bV_m, bV_sd, bInt_m, bInt_sd, 
                           g0_sd, gV_sd, gInt_sd) {
  # Simulate dataset
  sim_data <- simulate_dataset(n_subj, n_trials, 
                               b0_m, b0_sd, bV_m, bV_sd, bInt_m, bInt_sd, 
                               g0_sd, gV_sd, gInt_sd)
  
  # Fit mixed effects logistic regression
  full_formula <- 
    "choice ~ 1 + xVal1 + xVal2 + xInt1 + xInt2 + (1 + xVal1 + xVal2 + xInt1 + xInt2 || subject)" %>% 
    as.formula()
  
  model <- try(glmer(full_formula, data = sim_data, family = binomial,
                     control = glmerControl(optimizer = "bobyqa")),
               silent = TRUE)
  
  # If the model fails to converge, return NA for p-values.
  if(inherits(model, "try-error")){
    return(rep(NA, 2))
  }
  
  # Extract fixed effects coefficient table from the model summary.
  coefs <- summary(model)$coefficients
  
  # Identify rows corresponding to interaction terms.
  int_rows <- grep("xInt", rownames(coefs))
  
  # Calculate approximate p-values using the z-statistic.
  z_vals <- coefs[int_rows, "Estimate"] / coefs[int_rows, "Std. Error"]
  p_vals <- 2 * (1 - pnorm(abs(z_vals)))
  
  return(p_vals)
}

# ------------------------------------------------------------------------------

# Calculate prior beliefs about true parameters
# full_formula <- "choice ~ 1 + VLT + VLB + VRT + VRB + VLT:factor(DwellLT_Binary) + VLB:factor(DwellLB_Binary) + VRT:factor(DwellRT_Binary) + VRB:factor(DwellRB_Binary) + (1 + VLT + VLB + VRT + VRB + VLT:factor(DwellLT_Binary) + VLB:factor(DwellLB_Binary) + VRT:factor(DwellRT_Binary) + VRB:factor(DwellRB_Binary) || SubjectNumber)" %>% as.formula()
# glmer(full_formula, data=df_t, family = binomial,
#       control = glmerControl(optimizer = "bobyqa")) %>% summary()

# Set simulation parameters
n_sim <- 200            # number of simulation iterations
n_subj <- 50         # number of subjects
n_trials <- 30         # trials per subject

# True parameters
b0_m <- 0.008    # mean intercept, variation across simulations
b0_sd <- 0.02    # sd intercept, variation across simulations
bV_m <- 0.145    # mean effect of value, variation across simulations
bV_sd <- 0.003    # sd of effect of value, variation across simulations
bInt_m <- 0.012    # mean effect of attention, variation across simulations
bInt_sd <- 0.002    # sd of effect of attention, variation across simulations
g0_sd <- 0.03   # sd of deviation of intercept across subjects
gV_sd <- 0.005    # sd of deviation of effect of value across subjects
gInt_sd <- 0.005    # sd of deviation of effect of attention across subjects

# Detect number of cores and register parallel backend
n_cores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Run simulations in parallel
set.seed(4)  # for reproducibility
p_matrix <- foreach(i = 1:n_sim, .combine = rbind, .packages = c("lme4", "dplyr")) %dopar% {
  run_simulation(n_subj, n_trials, 
                 b0_m, b0_sd, bV_m, bV_sd, bInt_m, bInt_sd, 
                 g0_sd, gV_sd, gInt_sd)
}

# Stop the parallel cluster
stopCluster(cl)

# Calculate power for each interaction term
power_estimates <- c(mean(p_matrix[,1]<=0.05), mean(p_matrix[,2]>0.05))

# Display power estimates
power_df <- data.frame(
  Interaction = paste0("xInt", 1:2),
  Power = power_estimates
)
print(power_df)

# Overall average power across interaction terms
avg_power <- mean(power_estimates)
cat("Average estimated power across interaction terms:", avg_power, "\n")
