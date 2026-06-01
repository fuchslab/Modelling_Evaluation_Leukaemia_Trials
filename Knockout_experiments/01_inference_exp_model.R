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
# KO data
# data_full <- readRDS(paste0(folder.path, "/RDS/", "KO_data.rds"))

# Information about experiments, e.g. absolute input numbers and maximum number 
# of cells ever measured in laboratory for specific leukaemia samples
# info_table <- readRDS(paste0(folder.path, "/RDS/", "overview_experiments.rds"))


##### Generate exponential growth model for dMod #####
# ODE model
ode_equation <- eqnvec(x1 = "beta1*x1", x2 = "beta2*x2")
model <- odemodel(ode_equation, modelname = "exponential_model") 
x_dMod <- Xs(model) 

# Observable
observable <- eqnvec(y = "x1/(x1+x2)")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(error_pars, f = c(observable, ode_equation), attach.input = T, compile = T)

# Parameter transformations for original full model 
p_full <- P(
  eqnvec(
    beta1 = "beta1", beta2 = "beta2", x1 = "exp(x1)", x2 = "exp(x2)", 
    sigma = "exp(sigma)"), condition = "full") 

# Reparametrisations to make model structurally identifiable
equations <- getEquations(p_full)
equations[["full"]]["beta2"] <- "theta1+beta1"
equations[["full"]]["x2"] <- "exp(theta2)*exp(x1)"
p_global <- P(equations[["full"]], condition = "global")


################################################################################
######################### Parameter inference ##################################
################################################################################

###### Function which is parallelised by foreach #####
estimation_exp_model <- function(experiment_info){
  
  # Prepare mouse-specific data sets and parameter transformations for dMod
  data <- subset(
    data_full, sample == experiment_info$Sample & gene == experiment_info$Gene 
    & organ == "BM")
  data$name <- "y"
  unique_mice <- unique(data$mouse)
  data_dMod <- NULL
  p_dMod <- NULL
  
  for (j in 1:length(unique_mice)){
    data_m <- data[data$mouse == unique_mice[j], ]
    data_m <- data_m[, c("name", "time", "value")]
    data_m <- data_m[order(data_m$time), ]
    data_m$sigma <- NA
    data_m <- datalist(data_m)
    names(data_m) <- unique_mice[j]
    data_dMod <- data_dMod + data_m
    
    # add initial parameters for each mouse individually to the parameter transformations
    trafo <- getEquations(p_global, conditions = "global")
    trafo["x2"] <- paste0("exp(x1)/exp(theta2_", unique_mice[j], ")")
    p_dMod <- p_dMod + P(trafo, condition = unique_mice[j])
  }
  
  # Center of starting values for parameter estimation
  outerpars <- getParameters(p_dMod)
  outerpars_theta2 <- outerpars[grepl("theta2_", outerpars)]
  pouter <- structure(
    c(0.1, 0.2, log(experiment_info$x1), rep(log(1), length(outerpars_theta2)), log(0.05)), 
    names = c("theta1", "beta1", "x1", outerpars_theta2, "sigma"))
  
  # Objective function for parameter estimation with L2-constraint
  obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = 3)
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, 
                  rmax = 10, fits = 30, samplefun = "runif", iterlim = 500)
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  obj_value <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(
    bestfit)$value
  
  # Profile likelihoods
  profiles <- profile(
    obj = obj, pars = bestfit, whichPar = c("theta1", outerpars_theta2, "sigma"),
    alpha = 0.001, method = "integrate")
  
  # Profile likelihood-based confidence intervals (for 0.95 quantile of 
  # chi-squared distribution and for adapted threshold based on Cantelli's 
  # inequality)
  ci_0.99255_confint <- confint(profiles, val.column = "data", level = 0.99255)
  ci_0.95_confint <- confint(profiles, val.column = "data", level = 0.95)
  
  # If a confidence interval crosses confidence threshold more than twice,
  # adapt confidence interval; otherwise, just transform CIs to list
  ci_0.99255_list <- ci_extraction(ci_0.99255_confint, profiles, obj_value, 0.99255)
  ci_0.95_list <- ci_extraction(ci_0.95_confint, profiles, obj_value, 0.95)
  
  # Prepare output
  result_list <- list(
    sample = experiment_info$Sample, gene = experiment_info$Gene, 
    mllik = obj_value, bestfit = bestfit, profiles = profiles, 
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
  i = 1:nrow(info_table),
  .packages = c("dMod")) %dopar% {
    estimation_exp_model(
      experiment_info = info_table[i,])
  }

# Stop cluster
stopCluster(cl)

# Name result list elements by sample and gene
for (i in 1:length(result_list)){
  names(result_list)[i] <- paste(
    result_list[[i]]$sample, result_list[[i]]$gene, sep = "_")
}

# Save results
# filename <- "results_exp_model.rds"
# saveRDS(result_list, file = paste0(folder.path, "/RDS/", filename))


################################################################################
############ Testing chi-squared distribution assumption #######################
################################################################################

# Load estimation results 
# filename <- "results_exp_model.rds"
# result_list <- readRDS(paste0(folder.path, "/RDS/", filename))


##### Simulate bootstrap samples based on estimated parameters for original ##### 
##### data sets and re-estimate the model ##### 
set.seed(424)

# Number of bootstrap samples
num_data_sets <- 500 


### Function which is parallelised to re-estimate parameters for bootstrap samples 
bootstrap_exp_model <- function(data_sim, bestfit, p_dMod){
  
  # Generate bootstrap sample by adding noise to simulated data and prepare 
  # mouse-specific data sets for dMod 
  data_sim$value <- data_sim$value + rnorm(
    nrow(data_sim), mean = 0, sd = as.numeric(exp(bestfit["sigma"])))
  unique_mice <- unique(data_sim$mouse)
  data_dMod <- NULL
  
  for (j in 1:length(unique_mice)){
    data_m <- data_sim[data_sim$mouse == unique_mice[j], ]
    data_m <- data_m[, c("name", "time", "value")]
    data_m <- data_m[order(data_m$time), ]
    data_m$sigma <- NA
    data_m <- datalist(data_m)
    names(data_m) <- unique_mice[j]
    data_dMod <- data_dMod + data_m
  }
  
  # Objective function for parameter estimation with L2-constraint
  obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    bestfit, sigma = 3)
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, bestfit, studyname = "multiple_fits", rinit = 0.1, 
                  rmax = 10, fits = 30, samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, bestfit, studyname = "multiple_fits", rinit = 0.1, 
                    rmax = 10, fits = 10, samplefun = "runif", iterlim = 1000)
  }
  
  # Extract estimated parameters of best estimation run
  bestfit_b <- as.parvec(as.parframe(fits))
  
  # Calculate objective values for bootstrap estimates and original estimates
  mllik <- as.numeric(
    normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(bestfit_b)$value)
  mllik_star <- as.numeric(
    normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(
      bestfit)$value)
  
  # Prepare output
  result_b <- list(mllik = mllik, mllik_star = mllik_star, data_dMod = data_dMod)
  
  return(result_b)
}


### Parallelise parameter estimation of bootstrap samples for each knockout experiment 
# List for results
boot_results <- list()

# Loop through knockout experiments
for (s in 1:length(result_list)) {
  # Load information about knockout data set
  sample_name <- result_list[[s]]$sample
  gene_name <- result_list[[s]]$gene
  exp_name <- paste0(sample_name, "_", gene_name)
  bestfit <- result_list[[s]]$bestfit
  
  # Extract data set and prepare mouse-specific parameter transformations for dMod
  data <- subset(data_full, sample == sample_name & gene == gene_name)
  unique_mice <- unique(data$mouse)
  p_dMod <- NULL
  
  for (j in 1:length(unique_mice)){
    trafo <- getEquations(p_global, conditions = "global")
    trafo["x2"] <- paste0("exp(x1)/exp(theta2_", unique_mice[j], ")")
    p_dMod <- p_dMod + P(trafo, condition = unique_mice[j])
  }
  
  # Prepare bootstrap samples via predicting observable values for time points
  # of original data set and estimated parameters
  time_obs <- sort(data$time)
  prediction_list <- (g_dMod*x_dMod*p_dMod)(time_obs, bestfit)
  data_prediction <- data.frame(matrix(NA, nrow = 0, ncol = 3))
  for (j in 1:length(prediction_list)){
    df <- as.data.frame(prediction_list[[j]])
    df$mouse <- names(prediction_list)[j]
    data_prediction <- rbind(data_prediction, df)
  }
  data_sim <- left_join(data, data_prediction, by = c("time", "mouse"))
  data_sim$name <- "y"
  data_sim <- data_sim[, c("mouse", "name", "time", "y")]
  colnames(data_sim)[4] <- "value"
  data_sim <- data_sim[order(data_sim$time, data_sim$mouse), ]
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1) 
  registerDoParallel(cl) 
  
  # Parallelise parameter estimation for all bootstrap samples
  results <- foreach(
    i = 1:num_data_sets,
    .packages = c("dMod")) %dopar% {
      bootstrap_exp_model(
        data_sim = data_sim, 
        bestfit = bestfit, 
        p_dMod = p_dMod)
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Store bootstrap results for this knockout experiment
  boot_results <- c(boot_results, list(results))
  names(boot_results)[length(boot_results)] <- exp_name
}

# Save bootstrap samples and estimated models
# filename <- "results_exp_model_chi_boot.rds"
# saveRDS(boot_results, file = paste0(folder.path, "/RDS/", filename))


##### Empirical likelihood ratios (ELR), empirical cumulative distribution #####
#####  function (ECDF) and theoretical cumulative distribution function ##### 
##### (CDF) (chi_1^2-distribution) ##### 
set.seed(6)

# Load bootstrap results
# filename <- "results_exp_model_chi_boot.rds"
# boot_results <- readRDS(paste0(folder.path, "/RDS/", filename))


### Function which is parallelised to estimate model parameters for bootstrap 
### samples if parameters are consecutively fixed to estimated values of 
### original data set
ecdf_exp_model <- function(boot_result, pouter, p_dMod){
  
  # Extract bootstrap sample
  data_dMod <- boot_result$data_dMod
  
  # Objective function with L2-constraint
  obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = 3)
  
  # Parameter names for which empirical likelihood ratios should be calculated 
  # and output data frame
  par_names <- c("theta1", "sigma", paste0("theta2_", getConditions(p_dMod)))
  mllik_simple_df <- data.frame(matrix(NA, nrow = 1, ncol = (length(par_names) + 1)))
  colnames(mllik_simple_df) <- c("mllik_theta_all", par_names)
  mllik_simple_df[1, "mllik_theta_all"] <- boot_result$mllik
  
  # Loop through all parameters for which empirical likelihood ratios should be
  # calculated and estimate nested models
  for (j in 1:length(par_names)){
    # Fix parameter to "true" value
    fixed <- pouter[par_names[j]]
    
    # Estimate remaining parameters
    fits <- mstrust(
      obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, rinit = 0.1, 
      rmax = 10, studyname = "multiple_fits", fits = 10, samplefun = "runif", 
      iterlim = 500)
    
    # If no estimation has converged yet, run multiple estimations until at least
    # one estimation converges
    while (inherits(try(as.parframe(fits)), "try-error")) {
      fits <- mstrust(
        obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed, rinit = 0.1, 
        rmax = 10, studyname = "multiple_fits", fits = 10, samplefun = "runif", 
        iterlim = 1000)
    }
    
    # Extract estimated parameters of best estimation run
    fit_pars <- as.parvec(as.parframe(fits))
    
    # Store results
    mllik_simple <- normL2(
      data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(c(fit_pars, fixed))$value
    mllik_simple_df[1, par_names[j]] <- mllik_simple
  }
  
  return(mllik_simple_df)
}


### Parallelise estimations for each knockout experiment
# List for ECDFs
ecdf_results <- list()

# Loop through bootstrap results of the knockout experiments 
for (s in 1:length(boot_results)){
  
  # Extract parameter estimates from original data set and create center of
  # starting values for parameter estimation
  exp_name <- names(boot_results)[s]
  pouter <- result_list[[exp_name]]$bestfit
  
  # Prepare mouse-specific parameter transformations for dMod 
  unique_mice <- names(boot_results[[s]][[1]]$data_dMod)
  p_dMod <- NULL
  for (j in 1:length(unique_mice)){
    trafo <- getEquations(p_global, conditions = "global")
    trafo["x2"] <- paste0("exp(x1)/exp(theta2_", unique_mice[j], ")")
    p_dMod <- p_dMod + P(trafo, condition = unique_mice[j])
  }

  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl) 
  
  # Parallelise estimations for the bootstrap samples of a knockout experiment
  mllik_simple_results <- foreach(
    i = 1:num_data_sets,
    .packages=c("dMod")) %dopar% {
      ecdf_exp_model(
        boot_result = boot_results[[s]][[i]], 
        pouter = pouter,
        p_dMod = p_dMod)
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Combine results for this knockout experiment
  mllik_simple_results_df <- do.call(rbind, mllik_simple_results)
  
  # ECDF
  ecdf_df <- data.frame(matrix(
    NA, nrow = nrow(mllik_simple_results_df), ncol = ncol(mllik_simple_results_df)))
  colnames(ecdf_df) <- c("ecdf", colnames(mllik_simple_results_df)[-1])
  
  for (j in 2:ncol(mllik_simple_results_df)){
    elr <- mllik_simple_results_df[, j] - mllik_simple_results_df[, "mllik_theta_all"]
    elr <- sort(elr)
    if (j == 2) ecdf_df[, 1] <- ecdf_func(elr)
    ecdf_df[, j] <- pchisq(elr, 1)
  }
  
  ecdf_results <- c(ecdf_results, list(ecdf_df))
  names(ecdf_results)[length(ecdf_results)] <- exp_name
}

# Save ECDFs
# filename <- "results_exp_model_ecdf.rds"
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
# filename <- "results_exp_model_chi_boot.rds"
# boot_results <- readRDS(paste0(folder.path, "/RDS/", filename))
# num_data_sets <- length(boot_results[[1]])
# filename <- "results_exp_model_ecdf.rds"
# ecdf_results <- readRDS(paste0(folder.path, "/RDS/", filename))


### Classification of pp-plots of all parameters and classification at 0.95 quantile
### 1: perfect consensus, 2:conservative, 3: anti-conservative, 4:alternating

# Create tolerance region for perfect consensus region and allow some outliers
tol <- 3*diff_pred[-1]
outlier <- 50
p <- 0.95
p_idx <- 476 # index for which: which(pred_x==0.95) 
tol_p <- tol[p_idx]

# Data frame for number of points below consensus region (i.e. in anti-conservative region)
no_points_below_consensus_df <- data.frame(matrix(
  NA, nrow = length(ecdf_results), ncol = (3 + length(unique(data_full$mouse)))))
colnames(no_points_below_consensus_df) <- c(
  "experiment", "theta1", "theta3", paste0("theta2_", unique(data_full$mouse)))

# Data frame for pp-plot classification
classification_pp_df <- data.frame(matrix(
  NA, nrow = length(ecdf_results), ncol = (3 + length(unique(data_full$mouse)))))
colnames(classification_pp_df) <- c(
  "experiment", "theta1", "theta3", paste0("theta2_", unique(data_full$mouse)))

# Data frame for classification at 0.95 confidence quantile
classification_pp_thresh_df <- data.frame(matrix(
  NA, nrow = length(ecdf_results), ncol = (3 + length(unique(data_full$mouse)))))
colnames(classification_pp_thresh_df) <- c(
  "experiment", "theta1", "theta3", paste0("theta2_", unique(data_full$mouse)))

# Loop through the ECDFs of all knockout experiments
for (i in 1:length(ecdf_results)) {
  
  # Insert name of experiment into data frames
  no_points_below_consensus_df[i, "experiment"] <- names(ecdf_results)[i]
  classification_pp_df[i, "experiment"] <- names(ecdf_results)[i]
  classification_pp_thresh_df[i, "experiment"] <- names(ecdf_results)[i]
  
  # theta1: count number of points below consensus region and conduct pp-plot classifications
  theta1_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "theta1"]
  no_points_below_consensus_df[i, "theta1"] <- sum(theta1_diff_emp_theo <= -tol)
  classification_pp_df[i, "theta1"] <- classification_pp_plot(theta1_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "theta1"] <- 
    classification_pp_plot_thresh(p, theta1_diff_emp_theo, ecdf_results[[i]][, "theta1"], tol_p)
  
  # theta3: count number of points below consensus region and conduct pp-plot classifications
  theta3_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "sigma"]
  no_points_below_consensus_df[i, "theta3"] <- sum(theta3_diff_emp_theo <= -tol)
  classification_pp_df[i, "theta3"] <- classification_pp_plot(theta3_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "theta3"] <- 
    classification_pp_plot_thresh(p, theta3_diff_emp_theo, ecdf_results[[i]][, "sigma"], tol_p)
  
  # theta2 versions: count number of points below consensus region and conduct pp-plot classifications
  ecdf_theta2 <- ecdf_results[[i]][, startsWith(names(ecdf_results[[i]]), "theta2_")]
  for (j in 1:ncol(ecdf_theta2)){
    theta2_name <- colnames(ecdf_theta2)[j]
    theta2_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_theta2[, theta2_name]
    no_points_below_consensus_df[i, theta2_name] <- sum(theta2_diff_emp_theo <= -tol)
    classification_pp_df[i, theta2_name] <- classification_pp_plot(theta2_diff_emp_theo, tol, outlier)
    classification_pp_thresh_df[i, theta2_name] <-
      classification_pp_plot_thresh(p, theta2_diff_emp_theo, ecdf_theta2[, theta2_name], tol_p)
  }
}


################################################################################
###### Parameter estimation for knockout experiments and each sgRNA used #######
################################################################################

##### Function which is parallelised by foreach #####
estimation_exp_model <- function(experiment_info, input_info){
  
  # Extract data set
  data <- subset(data_full, sample == experiment_info$sample 
                 & gene == experiment_info$gene  & sgRNA==experiment_info$sgRNA)
  
  # If more than one mouse was used for this sgRNA, continue
  if (nrow(data) >= 4) {
    
    # Prepare mouse-specific data sets and parameter transformations for dMod 
    data$name <- "y"
    unique_mice <- unique(data$mouse)
    data_dMod <- NULL
    p_dMod <- NULL
    
    # create datalist separate for each mouse
    for (j in 1:length(unique_mice)){
      data_m <- data[data$mouse == unique_mice[j], ]
      data_m <- data_m[, c("name", "time", "value")]
      data_m <- data_m[order(data_m$time), ]
      data_m$sigma <- NA
      data_m <- datalist(data_m)
      names(data_m) <- unique_mice[j]
      data_dMod <- data_dMod + data_m
      
      trafo <- getEquations(p_global, conditions = "global")
      trafo["x2"] <- paste0("exp(x1)/exp(theta2_", unique_mice[j], ")")
      p_dMod <- p_dMod + P(trafo, condition = unique_mice[j])
    }
    
    # Center of starting values for parameter estimation
    outerpars <- getParameters(p_dMod)
    outerpars_theta2 <- outerpars[grepl("theta2_", outerpars)]
    pouter <- structure(
      c(0.1, 0.2, log(input_info$x1), rep(log(1), length(outerpars_theta2)), log(0.05)), 
      names = c("theta1", "beta1", "x1", outerpars_theta2, "sigma"))
    
    # Objective function for parameter estimation with L2-constraint
    obj <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
      pouter, sigma = 3)
    
    # Multiple parameter estimation with randomly chosen starting values
    fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, 
                    rmax = 10, fits = 30, samplefun = "runif", iterlim = 500)
    
    # Extract estimated parameters of best estimation run
    bestfit <- as.parvec(as.parframe(fits))
    
    # Objective value for best estimation run without L2-contribution
    obj_value <- normL2(data_dMod, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(
      bestfit)$value
    
    # Profile likelihoods
    profiles <- profile(
      obj = obj, pars = bestfit, whichPar = c("theta1", outerpars_theta2, "sigma"),
      alpha = 0.001, method = "integrate")
    
    # Profile likelihood-based confidence intervals (for 0.95 quantile of 
    # chi-squared distribution and for adapted threshold based on Cantelli's 
    # inequality)
    ci_0.99255_confint <- confint(profiles, val.column = "data", level = 0.99255)
    ci_0.95_confint <- confint(profiles, val.column = "data", level = 0.95)
    
    # If a confidence interval crosses confidence threshold more than twice,
    # adapt confidence interval; otherwise, just transform CIs to list
    ci_0.99255_list <- ci_extraction(ci_0.99255_confint, profiles, obj_value, 0.99255)
    ci_0.95_list <- ci_extraction(ci_0.95_confint, profiles, obj_value, 0.95)
    
    # t-tests
    data0 <- subset(data, time == 0)[, c("sgRNA", "mouse", "value")]
    data14 <- subset(data, time > 0 & time < 20)[, c("sgRNA", "mouse", "value")]
    if (nrow(data14) > 1){
      data14 <- inner_join(data0, data14, by = c("sgRNA", "mouse"), suffix = c("_0", "_14"))
      t_14 <- t.test(data14$value_0, data14$value_14, paired = TRUE, 
                     alternative = "two.sided")$p.value
    } else {
      t_14 <- NA
    }
    
    if (max(data$time) > 17){
      data_end <- subset(data, time == max(data$time))[, c("sgRNA", "mouse", "value")]
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
    
    # prepare output
    result_list <- list(
      sample = experiment_info$sample, gene = experiment_info$gene, 
      sgRNA = experiment_info$sgRNA, mllik = obj_value, bestfit = bestfit, 
      profiles = profiles,  ci_0.95_list = ci_0.95_list, 
      ci_0.99255_list = ci_0.99255_list, t_14 = t_14, t_end = t_end)
    
  } else {
    result_list <- NULL
  }
  
  return(result_list)
}


##### Run parameter estimations in parallel #####
set.seed(101)

# Create overview of experiments separated by sgRNAs
experiments_sgRNA <- unique(data_full[, c("sample","gene","sgRNA")])

# Detect the number of available cores and activate cluster
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Parallelise parameter estimations for each sgRNA
exp_model_results <- foreach(
  i = 1:nrow(experiments_sgRNA), 
  .packages = c("dMod", "dplyr")) %dopar% { 
    estimation_exp_model(
      experiment_info = experiments_sgRNA[i,],
      input_info = subset(info_table, Sample == experiments_sgRNA[i, "sample"] & Gene == experiments_sgRNA[i, "gene"]))
  }

# Stop cluster
stopCluster(cl)

# Remove experiments for which no estimation has been conducted
exp_model_results <- exp_model_results[!sapply(exp_model_results, is.null)]

# Name result list elements by sample, gene and sgRNA
for (j in 1:length(exp_model_results)){
  names(exp_model_results)[j] <- paste(
    exp_model_results[[j]]$sample, exp_model_results[[j]]$gene, 
    exp_model_results[[j]]$sgRNA, sep = "_")
}

# Save results
# filename <- "results_exp_model_sgRNA.rds"
# saveRDS(exp_model_results, file = paste0(folder.path, "/RDS/", filename))
