# R script to infer the parameters of the logistic growth model for the 
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
folder.path <- dirname(rstudioapi::getActiveDocumentContext()$path)

###### Load functions #####
functions.path <- "C:/Users/jwaesche/sciebo/01_Promotion/01_MyPhD/01_Leukaemieprojekt_R/1_Writing/Paper/Paper_Modelling/Creating Plots/Code_upload/utils.R"
source(functions.path)

##### Load data and information about input #####
# Data of knockout experiments not provided
data_full <- readRDS(paste0(folder.path,"/RDS/ko_data_genes_anonymised.rds"))
# data_full contains seven columns:
#   - "sample": name of leukaemia sample
#   - "gene": gene number
#   - "sgRNA": sgRNA number
#   - "mouse": mouse name or in vitro number (only mice experiments considered)
#   - "organ": organ from which measurement is taken; only bone marrow (BM) 
#              measurements are considered (NA for in vitro experiment)
#   - "time": measurement day
#   - "value": measurement value

# Absolute input cell numbers for all leukaemia samples
info_input <- readRDS(paste0(folder.path,"/RDS/","overview_experiments.rds"))
info_input$K <- info_input$K*1e6
info_input$x1 <- info_input$x1*1000


##### Generate logistic growth model for dMod #####
# ODE model
ode_equation <- eqnvec(x1 = "lambda1*x1*(1-(x1+x2)/K)",
                       x2 = "lambda2*x2*(1-(x1+x2)/K)")
model_dMod <- odemodel(ode_equation, modelname = "logistic_model") 
x_dMod <- Xs(model_dMod) 

# Observable
observable <- eqnvec(y = "log(x1/(x1+x2))-(sigma^2)/2")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(
  error_pars, f = c(observable, ode_equation), attach.input = TRUE, compile = TRUE)

# Parameter transformations
p_dMod <- P(trafo = eqnvec(
  lambda1 = "lambda1", lambda2 = "lambda2", K = "exp(K)", x1 = "exp(x1)", 
  x2 = "exp(x2)", sigma = "exp(sigma)"), condition = "C1")
outerpars <- getParameters(p_dMod)


################################################################################
######################## Parameter estimation ##################################
################################################################################

###### Function which is parallelised by foreach #####
estimation_log_model <- function(experiment_info){
  
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
  
  # Specify fixed parameters
  fixed_par <- c("K", "x1")
  
  # Center of starting values for parameter estimation
  pouter <- structure(
    c(0, 0.1, log(experiment_info$K), log(experiment_info$x1),
      log(experiment_info$x1), log(0.5)), names = outerpars)
  
  # Objective function for parameter estimation
  obj <- normL2(data, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = c(3, 3, 5, 5, 5, 3)) 
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 30, 
                  samplefun = "runif", iterlim = 500) 
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective function for parameter estimation
  obj_value <- normL2(data, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(c(bestfit, pouter[fixed_par]))$value
  
  # Profile likelihoods
  profiles <- profile(obj = obj, pars = bestfit, fixed = pouter[fixed_par], 
                      whichPar = c("lambda1", "lambda2", "x2", "sigma"), 
                      alpha = 0.001, method = "integrate")
  
  # Profile likelihood-based confidence intervals
  ci_confint <- confint(profiles, val.column = "data", level = 0.95)
  
  # If a confidence interval crosses confidence threshold more than twice,
  # adapt confidence interval; otherwise, just transform CIs to list
  ci_0.95_list <- ci_extraction(ci_confint, profiles, obj_value, 0.95)
  
  # Prepare output
  bestfit_df <- data.frame(
    sample = experiment_info$Sample, gene = experiment_info$Gene, mllik = obj_value, 
    lambda1 = bestfit["lambda1"], lambda2 = bestfit["lambda2"], K = exp(pouter["K"]), 
    x1 = exp(pouter["x1"]), x2 = exp(bestfit["x2"]), sigma = exp(bestfit["sigma"]))
  
  # Both bestfit_df and bestfit contain estimated (reparametrised) parameter values. 
  # bestfit_df is a data frame with one row for further processing and bestfit
  # is a parvec (a dMod class).
  result_list <- list(bestfit_df=bestfit_df, bestfit=bestfit,
                      ci_0.95_list=ci_0.95_list)
  return(result_list)
}


##### Estimations #####
set.seed(200)

# Detect the number of available cores and activate cluster
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Parallelise estimations for knockout data sets
log_model_results <- foreach(
  i = 1:nrow(info_input),
  .packages = c("dMod")) %dopar% {
    estimation_log_model(
      experiment_info = info_input[i,])
  }

# Stop cluster
stopCluster(cl)

# Name result list elements by sample and gene
for (j in 1:length(log_model_results)){
  names(log_model_results)[j] <- paste(log_model_results[[j]]$bestfit_df[1, c("sample", "gene")], collapse = "_")
}

# Save results
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_seed=200.rds"
# saveRDS(log_model_results, file = paste0(folder.path, "/RDS/", filename))


################################################################################
############ Testing chi-squared distribution assumption #######################
################################################################################

# Load estimation results
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_seed=200.rds"
# result_list <- readRDS(paste0(folder.path, "/RDS/", filename))


##### Simulate bootstrap samples based on estimated parameters for original ##### 
##### data sets and re-estimate the model ##### 

set.seed(525)

# Number of bootstrap samples
num_data_sets <- 500 


### Function which is parallelised to re-estimate parameters for bootstrap samples 
bootstrap_log_model <- function(data_b, pouter, fixed_par){
  
  # ´Prepare bootstrap sample for dMod
  data_b <- as.data.frame(data_b)
  data_b$name <- "y"
  data_b$sigma <- NA
  data_b <- data_b[, c(3, 1, 2, 4)]
  data_b <- datalist(C1 = data_b)
  
  # Objective function for parameter estimation
  obj <- normL2(data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = c(3, 3, 5, 5, 5, 3)) 
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter[!names(pouter)%in%fixed_par], fixed = pouter[fixed_par],
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter)%in%fixed_par], fixed = pouter[fixed_par],
                    studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, iterlim = 1000)
  }
  
  # Extract estimated parameters of best estimation run
  bestfit_b <- as.parvec(as.parframe(fits))
  
  # Prepare results
  result_b <- data.frame(
    mllik_star = as.numeric(normL2(data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(pouter)$value),
    mllik = as.numeric(normL2(data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(c(bestfit_b, pouter[fixed_par]))$value),
    lambda1 = as.numeric(bestfit_b["lambda1"]),
    lambda2 = as.numeric(bestfit_b["lambda2"]),
    x2 = as.numeric(exp(bestfit_b["x2"])),
    sigma = as.numeric(exp(bestfit_b["sigma"])))
  rownames(result_b) <- NULL
  
  return(result_b)
}


### Parallelise parameter estimation of bootstrap samples for each knockout experiment 
# List for results
boot_results <- list()

# Loop through knockout experiments
for (s in 1:length(result_list)) {
  
  # Load estimation results of original knockout data set
  bestfit_df <- result_list[[s]]$bestfit_df
  exp_name <- paste0(bestfit_df[1, "sample"], "_", bestfit_df[1, "gene"])
  
  # Specify fixed parameters and set starting values to estimates for original data set
  fixed_par <- c("K", "x1")
  pouter <- structure(
    c(as.numeric(result_list[[s]]$bestfit), as.numeric(log(bestfit_df[1, fixed_par]))),
    names = c(names(result_list[[s]]$bestfit), fixed_par))
  
  # Create bootstrap samples
  time_obs <- sort(subset(
    data_full, sample == bestfit_df[1, "sample"] & 
      gene == bestfit_df[1, "gene"] & organ == "BM")$time)
  num_obs <- length(time_obs)
  data_sim <- data.frame(time = time_obs)
  prediction <- (g_dMod*x_dMod*p_dMod)(time_obs, pouter)
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
      bootstrap_log_model(
        data_b = data_array[,,i], 
        pouter = pouter, 
        fixed_par = fixed_par)
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Save bootstrap results for this knockout experiment
  results_df <- do.call(rbind, results)
  boot_results <- c(boot_results, list(list(results_df = results_df, data_array = data_array)))
  names(boot_results)[length(boot_results)] <- exp_name
}

# Save all bootstrap samples and all estimated models
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_chi_boot_seed=525.rds"
# saveRDS(boot_results, file = paste0(folder.path, "/RDS/", filename))


##### Empirical likelihood ratios (ELR), empirical cumulative distribution #####
#####  function (ECDF) and theoretical cumulative distribution function ##### 
##### (CDF) (chi_1^2-distribution) ##### 
set.seed(8)

# Load bootstrap results
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_chi_boot_seed=525.rds"
# boot_results <- readRDS(paste0(folder.path, "/RDS/", filename))


### Function which is parallelised to estimate model parameters for bootstrap 
### samples if parameters are consecutively fixed to estimated values of 
### original data set
ecdf_log_model <- function(data_b, pouter){
  
  # Prepare data for dMod
  data_b <- as.data.frame(data_b)
  data_b$name <- "y"
  data_b$sigma <- NA
  data_b <- data_b[, c(3, 1, 2, 4)]
  data_b <- datalist(C1 = data_b)
  
  # Objective function for parameter estimation
  obj <- normL2(data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = c(3, 3, 5, 5, 5, 3)) 
  
  
  ### lambda1 fixed to "true" value alongside K and x1
  fixed <- pouter[c("lambda1", "K", "x1")]
  
  # Multiple parameter estimation with randomly chosen starting values and additionally lambda1 fixed
  fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                    studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, iterlim = 1000)
  }
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  mllik_simple_lambda1 <- normL2(
    data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(c(bestfit, fixed))$value
  
  
  ### lambda2 fixed to "true" value alongside K and x1
  fixed <- pouter[c("lambda2", "K", "x1")]
  
  # Multiple parameter estimation with randomly chosen starting values and additionally lambda2 fixed
  fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                    studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, iterlim = 1000)
  }
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  mllik_simple_lambda2 <- normL2(
    data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(c(bestfit, fixed))$value
  
  
  ### x2 fixed to "true" value alongside K and x1
  fixed <- pouter[c("K", "x1", "x2")]
  
  # Multiple parameter estimation with randomly chosen starting values and additionally x2 fixed
  fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                    studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, iterlim = 1000)
  }
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  mllik_simple_x2 <- normL2(
    data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(c(bestfit, fixed))$value
  
  
  ### sigma fixed to "true" value alongside K and x1
  fixed <- pouter[c("K", "x1", "sigma")]
  
  # Multiple parameter estimation with randomly chosen starting values and additionally sigma fixed
  fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                  studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, 
                  samplefun = "runif", iterlim = 500)
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% names(fixed)], fixed = fixed,
                    studyname = "multiple_fits", rinit = 0.1, rmax = 10, fits = 10, iterlim = 1000)
  }
  
  # Extract estimated parameters of best estimation run
  bestfit <- as.parvec(as.parframe(fits))
  
  # Objective value for best estimation run without L2-contribution
  mllik_simple_sigma <- normL2(
    data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod)(c(bestfit, fixed))$value
  
  # Collect log-likelihoods 
  mllik_simple_b <- data.frame(
    mllik_lambda1 = mllik_simple_lambda1, 
    mllik_lambda2 = mllik_simple_lambda2,
    mllik_x2 = mllik_simple_x2,
    mllik_sigma = mllik_simple_sigma
  )
  
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
  pouter <- structure(
    c(as.numeric(result_list[[exp_name]]$bestfit), 
      log(result_list[[exp_name]]$bestfit_df[1, "K"]),
      log(result_list[[exp_name]]$bestfit_df[1, "x1"])),
    names = c(names(result_list[[exp_name]]$bestfit), "K", "x1"))
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl) 
  
  # Parallelise estimations for the bootstrap samples of a knockout experiment
  mllik_simple_results <- foreach(
    i = 1:num_data_sets,
    .packages = c("dMod")) %dopar% {
      ecdf_log_model(
        data_b = data_array[,,i], 
        pouter = pouter)
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Combine obtained log-likelihoods of the samples for this knockout experiment
  mllik_simple_df <- do.call(rbind, mllik_simple_results)
  
  # Calculate ELRs
  elr_lambda1 <- mllik_simple_df$mllik_lambda1 - boot_results[[s]]$results_df$mllik
  elr_lambda2 <- mllik_simple_df$mllik_lambda2 - boot_results[[s]]$results_df$mllik
  elr_x2 <- mllik_simple_df$mllik_x2 - boot_results[[s]]$results_df$mllik
  elr_sigma <- mllik_simple_df$mllik_sigma - boot_results[[s]]$results_df$mllik
  
  elr_lambda1_sorted <- sort(elr_lambda1)
  elr_lambda2_sorted <- sort(elr_lambda2)
  elr_x2_sorted <- sort(elr_x2)
  elr_sigma_sorted <- sort(elr_sigma)
  
  # Calculate ECDFs
  ecdf_df <- data.frame(
    ecdf = ecdf_func(elr_lambda1_sorted),
    lambda1_theo = pchisq(elr_lambda1_sorted, 1),
    lambda2_theo = pchisq(elr_lambda2_sorted, 1),
    x2_theo = pchisq(elr_x2_sorted, 1),
    sigma_theo = pchisq(elr_sigma_sorted, 1))
  
  ecdf_results <- c(ecdf_results, list(ecdf_df))
  names(ecdf_results)[length(ecdf_results)] <- exp_name
}

# Save ECDFs
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_ecdf_seed=8_genes.rds"
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
axis_x <- rep(seq(0, 1, length.out = (num_real + 1))[-1], num_samples)
data_diff <- data.frame(diff = diff_emp_theo_numeric, q = axis_x)
poly_fit <- lm(diff ~ poly(q, 2), data_diff)

# Perfect consensus as the fitted polynomial of degree 2
pred_x <- seq(0, 1, length.out = num_real + 1)
diff_pred <- predict(poly_fit, data.frame(q = pred_x))


##### Classification of pp-plot graphs #####
# Load bootstrap and ECDF results
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_chi_boot_seed=525.rds"
# boot_results <- readRDS(paste0(folder.path, "/RDS/", filename))
# num_data_sets <- dim(boot_results[[1]]$data_array)[3]
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_ecdf_seed=8_genes.rds"
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
no_points_below_consensus_df <- data.frame(matrix(NA, nrow = length(ecdf_results), ncol = 5))
colnames(no_points_below_consensus_df) <- c(
  "experiment", "lambda1", "lambda2", "x2", "sigma")

# Data frame for pp-plot classification
classification_pp_df <- data.frame(matrix(NA, nrow = length(ecdf_results), ncol = 5))
colnames(classification_pp_df) <- c(
  "experiment", "lambda1", "lambda2", "x2", "sigma")

# Data frame for classification at confidence 0.95 quantile
classification_pp_thresh_df <- data.frame(matrix(NA, nrow = length(ecdf_results), ncol = 5))
colnames(classification_pp_thresh_df) <- c(
  "experiment", "lambda1", "lambda2", "x2", "sigma")

# Loop through the ECDFs of all knockout experiments
for (i in 1:length(ecdf_results)) {
  
  # Insert name of experiment into data frames
  no_points_below_consensus_df[i, "experiment"] <- names(ecdf_results)[i]
  classification_pp_df[i, "experiment"] <- names(ecdf_results)[i]
  classification_pp_thresh_df[i, "experiment"] <- names(ecdf_results)[i]
  
  # lambda1: count number of points below consensus region and conduct pp-plot classifications
  lambda1_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "lambda1_theo"]
  no_points_below_consensus_df[i, "lambda1"] <- sum(lambda1_diff_emp_theo <= -tol)
  classification_pp_df[i, "lambda1"] <- classification_pp_plot(lambda1_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "lambda1"] <- 
    classification_pp_plot_thresh(p, lambda1_diff_emp_theo, ecdf_results[[i]][, "lambda1_theo"], tol_p)
  
  # lambda2: count number of points below consensus region and conduct pp-plot classifications
  lambda2_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "lambda2_theo"]
  no_points_below_consensus_df[i, "lambda2"] <- sum(lambda2_diff_emp_theo <= -tol)
  classification_pp_df[i, "lambda2"] <- classification_pp_plot(lambda2_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "lambda2"] <- 
    classification_pp_plot_thresh(p, lambda2_diff_emp_theo, ecdf_results[[i]][, "lambda2_theo"], tol_p)
  
  # x2: count number of points below consensus region and conduct pp-plot classifications
  x2_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "x2_theo"]
  no_points_below_consensus_df[i, "x2"] <- sum(x2_diff_emp_theo <= -tol)
  classification_pp_df[i, "x2"] <- classification_pp_plot(x2_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "x2"] <- 
    classification_pp_plot_thresh(p, x2_diff_emp_theo, ecdf_results[[i]][, "x2_theo"], tol_p)
  
  # sigma: count number of points below consensus region and conduct pp-plot classifications
  sigma_diff_emp_theo <- ecdf_results[[i]][, "ecdf"] - ecdf_results[[i]][, "sigma_theo"]
  no_points_below_consensus_df[i, "sigma"] <- sum(sigma_diff_emp_theo <= -tol)
  classification_pp_df[i, "sigma"] <- classification_pp_plot(sigma_diff_emp_theo, tol, outlier)
  classification_pp_thresh_df[i, "sigma"] <- 
    classification_pp_plot_thresh(p, sigma_diff_emp_theo, ecdf_results[[i]][, "sigma_theo"], tol_p)
}


################################################################################
################ Bootstrapping for lambda2 - lambda1 ###########################
################################################################################

##### Load estimation results #####
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_seed=200.rds"
# result_list <- readRDS(paste0(folder.path, "/RDS/", filename))


##### Function which is parallelised by foreach #####
KO_bootstrap_log_model <- function(
    data, exp_name, groups, sample_sizes, pouter, fixed_par){
  
  # Bootstrap sampling of original design
  data_b_index <- unlist(lapply(unique(groups), function(gr) {
    group_indices <- which(groups == gr)
    n_samples <- sample_sizes[gr]
    sample(group_indices, n_samples, replace = TRUE)
  }))
  
  # Prepare bootstrap sample for dMod
  data_b <- data[data_b_index,]
  rownames(data_b) <- NULL
  data_b <- datalist(C1 = data_b[order(data_b$time),])
  
  # Objective function for parameter estimation
  obj <- normL2(data_b, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(
    pouter, sigma = c(3, 3, 5, 5, 5, 3))
  
  # Multiple parameter estimation with randomly chosen starting values
  fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
                  rinit = 0.1, rmax = 10, fits = 20, samplefun = "runif", iterlim = 500,
                  studyname = "multiple_fits")
  
  # If no estimation has converged yet, run multiple estimations until at least
  # one estimation converges
  while (inherits(try(as.parframe(fits)), "try-error")) {
    fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par], 
                    rinit = 0.1, rmax = 10, fits = 10, samplefun = "runif", iterlim = 1000,
                    studyname = "multiple_fits")
  }
  
  # Extract estimated parameters of best estimation run
  bestfit_b <- as.parvec(as.parframe(fits))
  
  # Compute difference of lambda2 and lambda1
  diff_lambda <- as.numeric(bestfit_b["lambda2"] - bestfit_b["lambda1"])
  
  return(diff_lambda)
}


##### Parallelise bootstrapping for each knockout experiment #####
set.seed(201)

# Number of bootstrap samples
B <- 999
boot_results <- list()

# Loop through all knockout experiments
for (s in 1:length(result_list)){ 
  
  # Load estimation results of original knockout data set
  bestfit_df <- result_list[[s]]$bestfit_df
  exp_name <- paste0(bestfit_df[1, "sample"], "_", bestfit_df[1, "gene"])
  
  # Extract data set, log-transform observations and prepare it for dMod
  data <- subset(data_full, sample == bestfit_df[1, "sample"] & 
                   gene == bestfit_df[1, "gene"] & organ == "BM")
  data <- data[, c("time", "value")]
  data <- data[order(data$time),]
  data$value <- log(data$value)
  data$name <- "y"
  data$sigma <- NA
  data <- data[, c(3, 1, 2, 4)]
  
  # Design framework for bootstrap samples based on original data set
  groups <- as.character(data[, "time"])
  sample_sizes <- sapply(unique(groups), function(x) 
    return(sum(data[, "time"] == as.numeric(x))))
  
  # Specify fixed parameters
  fixed_par <- c("K", "x1")
  
  # Center of starting values for parameter estimation
  pouter <- structure(
    c(as.numeric(result_list[[s]]$bestfit), as.numeric(log(bestfit_df[1, fixed_par]))),
    names = c(names(result_list[[s]]$bestfit), fixed_par))
  
  # Detect the number of available cores and activate cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl) 
  
  # Parallelise estimation for bootstrap samples of the knockout experiment
  boot_results_exp <- foreach(
    i = 1:B,
    .packages = c("dMod")) %dopar% {
      KO_bootstrap_log_model(
        data = data,
        exp_name = exp_name,
        groups = groups,
        sample_sizes = sample_sizes,
        pouter = pouter,
        fixed_par = fixed_par)
    }
  
  # Stop cluster
  stopCluster(cl)
  
  # Save bootstrap results for this knockout experiment
  boot_results <- c(boot_results, list(unlist(boot_results_exp)))
  names(boot_results)[length(boot_results)] <- exp_name
}

# Save bootstrap results
# filename <- "results_log_model_fixed=(K,x1)_L2=(3,3,5,5,5,3)_boot_seed=201.rds"
# saveRDS(boot_results, file = paste0(folder.path, "/RDS/", filename))

