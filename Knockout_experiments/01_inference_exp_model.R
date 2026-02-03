# R script to infer the parameters of the exponential growth model for the 
# knockout experiments.

################################################################################
################################ Preparations ##################################
################################################################################

##### Load R packages #####
library("dMod")
library("dplyr")
library("foreach")
library("doParallel")
library("parallel")

##### Folder path of this file to load data and to save results #####
# folder.path <- dirname(rstudioapi::getActiveDocumentContext()$path)

###### Load functions #####
# functions.path <- ".../utils.R"
# source(functions.path)

##### Load data and information about experiments #####
# Data of knockout experiments not provided
# data_full contains seven columns:
#   - "sample": name of leukaemia sample
#   - "gene": gene number
#   - "sgRNA": sgRNA number
#   - "mouse": mouse name or in vitro number (only mice experiments considered)
#   - "organ": organ from which measurement is taken; only bone marrow (BM) 
#              measurements are considered (NA for in vitro experiment)
#   - "time": measurement day
#   - "value": measurement value

# Information about experiments, e.g. absolute input numbers, which are depicted
# in Table 3
# info_input <- readRDS(paste0(folder.path, "/RDS/", "overview_experiments.rds"))
# info_input$K <- info_input$K*1e6
# info_input$x1 <- info_input$x1*1000


##### Generate exponential growth model for dMod #####
# ODE model
ode_equation <- eqnvec(x1 = "beta1*x1", x2 = "beta2*x2")
model_dMod <- odemodel(ode_equation, modelname = "exponential_model")
x_dMod <- Xs(model_dMod)

# Observable
observable <- eqnvec(y = "log(x1/(x1+x2))-(sigma^2)/2")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(
  error_pars, f = c(observable, ode_equation), attach.input = TRUE, compile = TRUE)

# Parameter transformations for original full model 
p_full <- P(trafo = eqnvec(
  x1 = "exp(x1)", x2 = "exp(x2)", beta1 = "beta1", beta2 = "beta2", 
  sigma = "exp(sigma)"), condition = "C1")
outerpars_full <- getParameters(p_full)

# Reparametrisations to make model structurally identifiable
equations <- getEquations(p_full)
equations[["C1"]]["beta2"] <- "theta1+beta1"
equations[["C1"]]["x2"] <- "exp(theta2)*exp(x1)"
p_repar <- P(equations[["C1"]], condition = "C1")
outerpars <- getParameters(p_repar)
theta <- c("theta1", "theta2", "sigma")


################################################################################
######################### Parameter inference ##################################
################################################################################

###### Function which is parallelised by foreach #####
estimation_exp_model <- function(experiment_info){
  
  # Extract data set, prepare it for dMod and log-transform observations
  data <- subset(data_full, sample == experiment_info$Sample & 
                   gene == experiment_info$Gene & organ == "BM")
  data <- data[, c("time", "value")]
  data <- data[order(data$time),]
  data$value <- log(data$value)
  data$name <- "y"
  data$sigma <- NA
  data <- data[, c(3, 1, 2, 4)]
  data <- datalist(C1 = data)
  
  # Center of starting values for parameter estimation
  pouter <- structure(c(3, log(2/3), 0.1, 0.1, log(0.5)), names = outerpars)
  
  # Objective function for parameter estimation
  obj <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod) + constraintL2(pouter, sigma = 5)
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                  fits = 30, samplefun = "runif", iterlim = 500)
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  obj_value <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(bestfit)$value
  
  # Profile likelihoods
  profiles <- profile(obj = obj, pars = bestfit, whichPar = theta, alpha = 0.0001,
                      method = "integrate")
  
  # Profile likelihood-based confidence intervals (for 0.95 quantile of 
  # chi-squared distribution and for adapted threshold based on Cantelli's 
  # inequality)
  ci_0.99255_confint <- confint(profiles, val.column = "data", level = 0.99255) 
  ci_0.95_confint <- confint(profiles, val.column = "data", level = 0.95)
  
  # If a confidence interval crosses confidence threshold more than twice,
  # adapt confidence interval; otherwise, just transform CIs to list
  ci_0.99255_list <- ci_extraction(ci_0.99255_confint, profiles, obj_value, 0.99255)
  ci_0.95_list <- ci_extraction(ci_0.95_confint, profiles, obj_value, 0.95)
  
  # Prepare output list
  result_list <- list(
    sample = experiment_info$Sample, gene = experiment_info$Gene, bestfit = bestfit, 
    ci_0.95_list = ci_0.95_list, ci_0.99255_list = ci_0.99255_list)
  return(result_list)
}


##### Estimations #####
set.seed(100)

# Detect the number of available cores and activate cluster
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Parallelise estimations for knockout data sets
result_list <- foreach(
  i = 1:nrow(info_input), 
  .packages = c("dMod")) %dopar% {
    estimation_exp_model(
      experiment_info = info_input[i,])
  }

# Stop cluster
stopCluster(cl)

# Name result list elements by sample and gene
for (j in 1:length(result_list)){
  names(result_list)[j] <- paste0(result_list[[j]]$sample, "_", result_list[[j]]$gene)
}

# Save results
# filename <- "results_exp_model_L2=5_seed=100.rds"
# saveRDS(result_list, file = paste0(folder.path, "/RDS/", filename))


################################################################################
############ Testing chi-squared distribution assumption #######################
################################################################################

# Load estimation results 
# filename <- "results_exp_model_L2=5_seed=100.rds"
# result_list <- readRDS(paste0(folder.path, "/RDS/", filename))


##### Simulate bootstrap samples based on estimated parameters for original ##### 
##### data sets and re-estimate the model ##### 
set.seed(424)

# Number of bootstrap samples
num_data_sets <- 500 


### Function which is parallelised to re-estimate parameters for bootstrap samples 
bootstrap_exp_model <- function(data_b, pouter){
  
  # ´Prepare bootstrap sample for dMod
  data_b <- as.data.frame(data_b)
  data_b$name <- "y"
  data_b$sigma <- NA
  data_b <- data_b[, c(3, 1, 2, 4)]
  data_b <- datalist(C1 = data_b)
  
  # Objective function for parameter estimation
  obj <- normL2(data_b, g_dMod*x_dMod*p_repar, errmodel = error_dMod) + constraintL2(pouter, sigma = 5)
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10,
                  fits = 10, samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10,
                    fits = 10, samplefun = "runif", iterlim = 1000)
  }
  
  # Extract estimated parameters of best estimation run
  bestfit_b <- as.parvec(as.parframe(fits))
  
  # Prepare results
  result_b <- data.frame(
    mllik_star = as.numeric(normL2(data_b, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(pouter)$value),
    mllik = as.numeric(normL2(data_b, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(bestfit_b)$value),
    theta1 = as.numeric(bestfit_b["theta1"]),
    theta2 = as.numeric(exp(bestfit_b["theta2"])),
    theta3 = as.numeric(exp(bestfit_b["sigma"])))
  rownames(result_b) <- NULL
  
  return(result_b)
}


### Parallelise parameter estimation of bootstrap samples for each knockout experiment 
# List for results
boot_results <- list()

# Loop through knockout experiments
for (s in 1:length(result_list)) {
  
  # Load experiment information
  sample_name <- result_list[[s]]$sample
  gene_nr <- result_list[[s]]$gene

  # Set starting values to estimates for original data set
  pouter <- result_list[[s]]$bestfit
  
  # Create bootstrap samples
  time_obs <- sort(subset(
    data_full, sample == sample_name & gene == gene_nr & organ == "BM")$time)
  num_obs <- length(time_obs)
  data_sim <- data.frame(time = time_obs)
  prediction <- (g_dMod*x_dMod*p_repar)(time_obs, pouter)
  data_sim <- merge(data_sim, prediction$C1, by = "time", all.x = TRUE)
  data_array <- array(NA, dim = c(num_obs, 2, num_data_sets), 
                      dimnames = list(NULL, c("time", "value"), NULL))
  data_array[,1,] <- data_sim[,1]
  data_array[,2,] <- data_sim[,2]
  data_array[,2,] <- apply(data_array[,2,], 2, function(obs) {
    random_vector <- rnorm(num_obs, 0, exp(pouter["sigma"]))
    obs + random_vector
  })
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1) 
  registerDoParallel(cl) 
  
  # Parallelise parameter estimation for all bootstrap samples
  results <- foreach(
    i = 1:num_data_sets, 
    .packages = c("dMod")) %dopar% {
      bootstrap_exp_model(
        data_b = data_array[,,i], 
        pouter = pouter)
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Store bootstrap results for this knockout experiment
  results_df <- do.call(rbind, results)
  boot_results <- c(boot_results, list(list(results_df = results_df, data_array = data_array)))
  names(boot_results)[length(boot_results)] <- paste0(sample_name, "_", gene_nr)
}

# Save all bootstrap samples and all estimated models
# filename <- "results_exp_model_L2=5_chi_boot_seed=424.rds"
# saveRDS(boot_results, file = paste0(folder.path, "/RDS/", filename))


##### Empirical likelihood ratios (ELR), empirical cumulative distribution #####
#####  function (ECDF) and theoretical cumulative distribution function ##### 
##### (CDF) (chi_1^2-distribution) ##### 
set.seed(6)

# Load bootstrap results
# filename <- "results_exp_model_L2=5_chi_boot_seed=424.rds"
# boot_results <- readRDS(paste0(folder.path, "/RDS/", filename))


### Function which is parallelised to estimate model parameters for bootstrap 
### samples if parameters are consecutively fixed to estimated values of 
### original data set
ecdf_exp_model <- function(data_b, pouter){
  
  # Prepare data for dMod
  data_b <- as.data.frame(data_b)
  data_b$name <- "y"
  data_b$sigma <- NA
  data_b <- data_b[, c(3, 1, 2, 4)]
  data_b <- datalist(C1=data_b)
  
  # Objective function for parameter estimation
  obj <- normL2(data_b, g_dMod*x_dMod*p_repar, errmodel = error_dMod) + constraintL2(pouter, sigma = 5)
  
  
  ### theta1 fixed to "true" value
  fixed <- pouter["theta1"]
  # Multiple parameter estimation with randomly chosen starting values and theta1 fixed
  fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, 
                  rinit = 0.1, rmax = 10, studyname = "multiple_fits", fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, 
                    rinit = 0.1, rmax = 10, studyname = "multiple_fits", fits = 10, 
                    samplefun = "runif", iterlim = 1000)
  }
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  # Objective value for best estimation run without L2-contribution
  mllik_simple_theta1 <- normL2(
    data_b, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(c(bestfit, fixed))$value
  
  
  ### theta2 fixed to "true" value
  fixed <- pouter["theta2"]
  # Multiple parameter estimation with randomly chosen starting values and theta2 fixed
  fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, 
                  rinit = 0.1, rmax = 10, studyname = "multiple_fits", fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, 
                    rinit = 0.1, rmax = 10, studyname = "multiple_fits", fits = 10, 
                    samplefun = "runif", iterlim = 1000)
  }
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  # Objective value for best estimation run without L2-contribution
  mllik_simple_theta2 <- normL2(
    data_b, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(c(bestfit, fixed))$value
  
  
  ### theta3 fixed to "true" value
  fixed <- pouter["sigma"]
  # Multiple parameter estimation with randomly chosen starting values and theta2 fixed
  fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, 
                  rinit = 0.1, rmax = 10, studyname = "multiple_fits", fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, 
                    rinit = 0.1, rmax = 10, studyname = "multiple_fits", fits = 10, 
                    samplefun = "runif", iterlim = 1000)
  }
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  # Objective value for best estimation run without L2-contribution
  mllik_simple_theta3 <- normL2(
    data_b, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(c(bestfit, fixed))$value
  
  # Collect log-likelihoods 
  mllik_simple_b <- data.frame(
    mllik_theta1 = mllik_simple_theta1, 
    mllik_theta2 = mllik_simple_theta2,
    mllik_theta3 = mllik_simple_theta3)
  
  return(mllik_simple_b=mllik_simple_b)
}


### Parallelise estimations for each knockout experiment
# List for ECDFs
ecdf_results <- list()

# Loop through bootstrap results of the knockout experiments 
for (s in 1:length(boot_results)){ 
  
  # Extract bootstrap samples
  exp_name <- names(boot_results)[s]
  data_array <- boot_results[[s]]$data_array
  
  # Set starting values to estimated parameters of original data set
  pouter <- result_list[[exp_name]]$bestfit
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl) 
  
  # Parallelise estimations for the bootstrap samples of a knockout experiment
  mllik_simple_results <- foreach(
    i = 1:num_data_sets, 
    .packages = c("dMod")) %dopar% {
      ecdf_exp_model(
        data_b = data_array[,,i], 
        pouter = pouter)
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Combine obtained log-likelihoods of the samples for this knockout experiment
  mllik_simple_df <- do.call(rbind, mllik_simple_results)
  
  # Calculate ELRs
  elr_theta1 <- mllik_simple_df$mllik_theta1 - boot_results[[s]]$results_df$mllik
  elr_theta2 <- mllik_simple_df$mllik_theta2 - boot_results[[s]]$results_df$mllik
  elr_theta3 <- mllik_simple_df$mllik_theta3 - boot_results[[s]]$results_df$mllik
  
  elr_theta1_sorted <- sort(elr_theta1)
  elr_theta2_sorted <- sort(elr_theta2)
  elr_theta3_sorted <- sort(elr_theta3)
  
  # Calculate ECDFs
  ecdf_df <- data.frame(
    ecdf = ecdf_func(elr_theta1_sorted),
    theta1_theo = pchisq(elr_theta1_sorted, 1),
    theta2_theo = pchisq(elr_theta2_sorted, 1),
    theta3_theo = pchisq(elr_theta3_sorted, 1))
  
  ecdf_results <- c(ecdf_results, list(ecdf_df))
  names(ecdf_results)[length(ecdf_results)] <- exp_name
}

# Save ECDFs
# filename <- "results_exp_model_L2=5_ecdf_seed=6.rds"
# saveRDS(ecdf_results, file = paste0(folder.path, "/RDS/", filename))


##### Derive perfect consensus #####
# Generate samples from chi-squared distributio with one degree of freedom
set.seed(444)
num_samples <- 1000
num_real <- 500
chi_data <- matrix(rchisq(num_samples*num_real, 1), nrow = num_real, ncol = num_samples)
chi_data_sorted <- apply(chi_data, 2, sort)

# Calculate empirical and theoretical CDF
chi_data_ecdf <- apply(chi_data_sorted, 2, ecdf_func)
chi_data_ecdf_theoretical <- apply(chi_data_sorted, 2, pchisq, df = 1)

# Calculate deviation from theoretical chi-squared distribution with one degree of freedom
# and fit a polynomial with degree 2 to deviation
diff_emp_theo <- chi_data_ecdf - chi_data_ecdf_theoretical
diff_emp_theo_abs <- abs(diff_emp_theo)
diff_emp_theo_numeric <- as.numeric(diff_emp_theo_abs)
axis_x <- rep(seq(0, 1, length.out= (num_real+1))[-1], num_samples)
data_diff <- data.frame(diff = diff_emp_theo_numeric, q = axis_x)
poly_fit <- lm(diff ~ poly(q, 2), data_diff)

# Perfect consensus as the fitted polynomial of degree 2
pred_x <- seq(0,1, length.out = num_real + 1)
diff_pred <- predict(poly_fit, data.frame(q = pred_x))


##### Classification of pp-plot graphs #####
# Load bootstrap and ECDF results
# filename <- "results_exp_model_L2=5_chi_boot_seed=424.rds"
# boot_results <- readRDS(paste0(folder.path,"/RDS/",filename))
# num_data_sets <- dim(boot_results[[1]]$data_array)[3]
# filename <- "results_exp_model_L2=5_ecdf_seed=6.rds"
# ecdf_results <- readRDS(paste0(folder.path,"/RDS/",filename))


### Classification of pp-plots of all parameters and classification at 0.95 quantile
### 1: perfect consensus, 2:conservative, 3: anti-conservative, 4:alternating

# Create tolerance region for perfect consensus region and allow some outliers
tol <- 3*diff_pred[-1]
outlier <- 50
p <- 0.95
p_idx <- 476 # index for which: which(pred_x==0.95) 
tol_p <- tol[p_idx]

# Data frame for number of points below consensus region (i.e. in anti-conservative region)
no_points_below_consensus_df <- data.frame(matrix(NA, nrow = length(ecdf_results), ncol = 4))
colnames(no_points_below_consensus_df) <- c("experiment", "theta1", "theta2", "theta3")

# Data frame for pp-plot classification
classification_pp_df <- data.frame(matrix(NA, nrow = length(ecdf_results), ncol = 4))
colnames(classification_pp_df) <- c("experiment", "theta1", "theta2", "theta3")

# Data frame for classification at 0.95 confidence quantile
classification_pp_thresh_df <- data.frame(matrix(NA, nrow = length(ecdf_results), ncol = 4))
colnames(classification_pp_thresh_df) <- c("experiment", "theta1", "theta2", "theta3")

# Loop through the ECDFs of all knockout experiments
for (i in 1:length(ecdf_results)) {
  
  # Insert name of experiment into data frames
  no_points_below_consensus_df[i, "experiment"] <- names(ecdf_results)[i]
  classification_pp_df[i, "experiment"] <- names(ecdf_results)[i]
  classification_pp_thresh_df[i, "experiment"] <- names(ecdf_results)[i]
  
  # theta1: count number of points below consensus region and conduct pp-plot classifications
  theta1_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "theta1_theo"]
  no_points_below_consensus_df[i, "theta1"] <- sum(theta1_diff_emp_theo <= -tol)
  classification_pp_df[i, "theta1"] <- classification_pp_plot(theta1_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "theta1"] <- 
    classification_pp_plot_thresh(p, theta1_diff_emp_theo, ecdf_results[[i]][, "theta1_theo"], tol_p)
  
  # theta2: count number of points below consensus region and conduct pp-plot classifications
  theta2_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "theta2_theo"]
  no_points_below_consensus_df[i, "theta2"] <- sum(theta2_diff_emp_theo <= -tol)
  classification_pp_df[i, "theta2"] <- classification_pp_plot(theta2_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "theta2"] <- 
    classification_pp_plot_thresh(p, theta2_diff_emp_theo, ecdf_results[[i]][, "theta2_theo"], tol_p)
  
  # theta3: count number of points below consensus region and conduct pp-plot classifications
  theta3_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "theta3_theo"]
  no_points_below_consensus_df[i, "theta3"] <- sum(theta3_diff_emp_theo <= -tol)
  classification_pp_df[i, "theta3"] <- classification_pp_plot(theta3_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "theta3"] <- 
    classification_pp_plot_thresh(p, theta3_diff_emp_theo, ecdf_results[[i]][, "theta3_theo"], tol_p)
}


################################################################################
###### Parameter estimation for knockout experiments and each sgRNA used #######
################################################################################

##### Function which is parallelised by foreach #####
estimation_exp_model <- function(experiment_info){
  
  # Extract data set 
  data_raw <- subset(
    data_full, sample == experiment_info$sample & gene == experiment_info$gene & 
      sgRNA == experiment_info$sgRNA & organ == "BM")
  
  # If number observations is less than 3, we do not consider this combination
  if (nrow(data_raw) < 3){
    result_list <- NULL
  } else {
    
    # Prepare data set for dMod and log-transform observations
    data <- data_raw[, c("time", "value")]
    data <- data[order(data$time),]
    data$value <- log(data$value)
    data$name <- "y"
    data$sigma <- NA
    data <- data[, c(3, 1, 2, 4)]
    data <- datalist(C1 = data)
    
    # Center of starting values for parameter estimation
    pouter <- structure(c(3, log(2/3), 0.1, 0.1, log(0.5)), names = outerpars)
    
    # Objective function for parameter estimation
    obj <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod) + constraintL2(pouter, sigma = 5)
    
    # Multiple parameter estimation with randomly chosen starting values
    fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                    fits = 30, samplefun = "runif", iterlim = 500)
    
    # Extract estimated parameters of best estimation run
    bestfit <- as.parvec(as.parframe(fits))
    
    # Objective value for best estimation run without L2-contribution
    obj_value <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(bestfit)$value
    
    # Profile likelihoods
    profiles <- profile(obj = obj, pars = bestfit, whichPar = theta, alpha = 0.0001,
                        method = "integrate")
    
    # Profile likelihood-based confidence intervals (for 0.95 quantile of 
    # chi-squared distribution and for adapted threshold based on Cantelli's 
    # inequality)
    ci_0.99255_confint <- confint(profiles, val.column = "data", level = 0.99255) 
    ci_0.95_confint <- confint(profiles, val.column = "data", level = 0.95)
    
    # If a confidence interval crosses confidence threshold more than twice,
    # adapt confidence interval; otherwise, just transform CIs to list
    ci_0.99255_list <- ci_extraction(ci_0.99255_confint, profiles, obj_value, 0.99255)
    ci_0.95_list <- ci_extraction(ci_0.95_confint, profiles, obj_value, 0.95)
    
    # Conduct t-tests (if applicable)
    data0 <- subset(data_raw, time == 0)[, c("sgRNA", "mouse", "value")]
    data14 <- subset(data_raw, time > 0 & time < 20)[, c("sgRNA", "mouse", "value")]
    if (nrow(data14) > 1){
      data14 <- inner_join(data0, data14, by = c("sgRNA", "mouse"), suffix = c("_0", "_14"))
      t_14 <- t.test(data14$value_0, data14$value_14, paired = TRUE, 
                     alternative = "two.sided")$p.value
    } else {
      t_14 <- NA
    }
    
    if (max(data_raw$time) > 17){
      data_end <- subset(data_raw, time == max(data_raw$time))[, c("sgRNA", "mouse", "value")]
      if (nrow(data_end) > 1){
        data_end <- inner_join(data0, data_end, by = c("sgRNA", "mouse"), suffix = c("_0", "_end"))
        t_end <- t.test(data_end$value_0, data_end$value_end, paired = TRUE,
                        alternative = "two.sided")$p.value
      } else {
        t_end <- NA
      }
    } else {
      t_end <- NA
    }
    
    # Prepare output list
    result_list <- list(
      sample = experiment_info$sample, gene = experiment_info$gene, 
      sgRNA = experiment_info$sgRNA, bestfit = bestfit, ci_0.95_list = ci_0.95_list, 
      ci_0.99255_list = ci_0.99255_list, t_14 = t_14, t_end = t_end)
  }
  return(result_list)
}


##### Run parameter estimations in parallel #####
set.seed(101)

# Create overview of experiments separated by sgRNAs
experiments_sgRNA <- unique(
  subset(data_full, organ == "BM")[, c("sample", "gene", "sgRNA")])

# Detect the number of available cores and activate cluster
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Parallelise parameter estimations for each sgRNA
exp_model_results <- foreach(
  i = 1:nrow(experiments_sgRNA), 
  .packages = c("dMod", "dplyr")) %dopar% {
    estimation_exp_model(
      experiment_info = experiments_sgRNA[i,])
  }

# Stop cluster
stopCluster(cl)

# Name result list elements by sample, gene and sgRNA
for (j in 1:length(exp_model_results)){
  names(exp_model_results)[j] <- paste0(
    exp_model_results[[j]]$sample, "_", exp_model_results[[j]]$gene, "_",
    exp_model_results[[j]]$sgRNA)
}

# Save results
# filename <- "results_exp_model_sgRNA_L2=5_seed=101.rds"
# saveRDS(exp_model_results, file = paste0(folder.path,"/RDS/",filename))


################################################################################
###### Parameter estimation for knockout experiments and each mouse used #######
################################################################################

##### Function which is parallelised by foreach #####
estimation_exp_model <- function(experiment_info){
  
  # Extract data set 
  data_raw <- subset(
    data_full, sample == experiment_info$sample & gene == experiment_info$gene & 
      sgRNA == experiment_info$sgRNA & mouse == experiment_info$mouse & organ == "BM")
  
  # Prepare data for dMod and log-transform observations
  data <- data_raw[, c("time", "value")]
  data <- data[order(data$time),]
  data$value <- log(data$value)
  data$name <- "y"
  data$sigma <- NA
  data <- data[, c(3, 1, 2, 4)]
  data <- datalist(C1 = data)
  
  # Center of starting values for parameter estimation
  pouter <- structure(c(3, log(2/3), 0.1, 0.1, log(0.5)), names = outerpars)
  
  # Objective function for parameter estimation
  obj <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod) + constraintL2(pouter, sigma = 1)
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                  fits = 50, samplefun = "runif", iterlim = 500)
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  obj_value <- normL2(data, g_dMod*x_dMod*p_repar, errmodel = error_dMod)(bestfit)$value
  
  # Profile likelihoods
  profiles <- profile(obj = obj, pars = bestfit, whichPar = theta, alpha = 0.0001,
                      method = "integrate")
  
  # Profile likelihood-based confidence intervals (for 0.95 quantile of 
  # chi-squared distribution and for adapted threshold based on Cantelli's 
  # inequality)
  ci_0.99255_confint <- confint(profiles, val.column = "data", level = 0.99255) 
  ci_0.95_confint <- confint(profiles, val.column = "data", level = 0.95)
  
  # If a confidence interval crosses confidence threshold more than twice,
  # adapt confidence interval; otherwise, just transform CIs to list
  ci_0.99255_list <- ci_extraction(ci_0.99255_confint, profiles, obj_value, 0.99255)
  ci_0.95_list <- ci_extraction(ci_0.95_confint, profiles, obj_value, 0.95)
  
  # Prepare output list
  result_list <- list(
    sample = experiment_info$sample, gene = experiment_info$gene, 
    sgRNA = experiment_info$sgRNA, mouse = experiment_info$mouse, 
    t_end = max(data_raw$time), bestfit = bestfit,  ci_0.95_list = ci_0.95_list, 
    ci_0.99255_list = ci_0.99255_list)
  
  return(result_list)
}


##### Run parameter estimations in parallel #####
set.seed(102)

# create overview of experiments separated by mice
experiments_mouse <- unique(
  subset(data_full, organ == "BM")[, c("sample", "gene", "sgRNA", "mouse")])

# Detect the number of available cores and activate cluster
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl) 

# Parallelise parameter estimations for each mouse
exp_model_results <- foreach(
  i = 1:nrow(experiments_mouse), 
  .packages = c("dMod","dplyr")) %dopar% {
    estimation_exp_model(
      experiment_info = experiments_mouse[i,])
  }

# Stop cluster
stopCluster(cl)

# Name result list elements by sample, gene, sgRNA and mouse
for (j in 1:length(exp_model_results)){
  names(exp_model_results)[j] <- paste0(
    exp_model_results[[j]]$sample, "_", exp_model_results[[j]]$gene, "_",
    exp_model_results[[j]]$sgRNA, "_", exp_model_results[[j]]$mouse)
}

# Save results
# filename <- "results_exp_model_mouse_L2=1_seed=102.rds"
# saveRDS(exp_model_results, file = paste0(folder.path, "/RDS/", filename))

