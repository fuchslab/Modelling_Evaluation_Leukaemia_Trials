# R script to investigate structural identifiability of the exponential and the
# logistic growth model. 

################################################################################
################################# Preparation ##################################
################################################################################

##### Load R packages #####
library("dMod")
library("ggplot2")
library("latex2exp")


################################################################################
###################### Exponential growth model ################################
################################################################################

##### Generate exponential growth model for dMod #####
# ODE model
ode_equation <- eqnvec(x1 = "beta1*x1", x2 = "beta2*x2")
model_dMod <- odemodel(ode_equation, modelname = "exponential_model")
x_dMod <- Xs(model_dMod)

# Observable model (for both estimation and simulation)
observable <- eqnvec(y = "log(x1/(x1+x2))-(sigma^2)/2")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)
observable_sim <- eqnvec(y = "x1/(x1+x2)")
g_sim <- Y(observable_sim, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(
  error_pars, f = c(observable, ode_equation), attach.input = TRUE, compile = TRUE)

# Parameter transformations for original full model 
p_full <- P(trafo=eqnvec(
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


##### Simulate data for structural identifiability analysis #####
set.seed(90)

# Parameter values
pars <- c(x1 = 10, x2 = 10, beta1 = 0.1, beta2 = 0.2)
sigma <- 0.2

# Time range
times <- seq(0, 50, 1)

# Simulation with dMod 
out <- (g_sim*x_dMod)(times, pars)
data_sim <- do.call(rbind, replicate(20, out[[1]], simplify = FALSE))
data_sim <- data_sim[order(data_sim[,1]),]

# Multiplying simulated observable with log-normally distributed measurement error
data_sim[,2] <- data_sim[,2]*rlnorm(nrow(data_sim), (-(sigma^2)/2), sigma)
data <- data.frame(name = rep("y", nrow(data_sim)), time = data_sim[,1], 
                   value = log(data_sim[,2]), sigma = rep(NA,nrow(data_sim)))
data <- datalist(C1 = data)


##### Identifiability of full model #####
### Parameter estimation
set.seed(99)

# Center of starting values for parameter estimation
pouter <- structure(
  c(log(pars["x1"]), log(pars["x2"]), pars["beta1"], pars["beta2"],
    log(sigma)), names = outerpars_full)

# Objective function
obj <- normL2(
  data, g_dMod*x_dMod*p_full, errmodel = error_dMod) + constraintL2(pouter, sigma = 5)

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                cores = 4, fits = 50, samplefun = "rnorm", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(obj = obj, pars = bestfit, whichPar = names(bestfit), 
                    alpha = 0.05, cores = 4)
plotProfile(profiles, mode == "data")


##### Identifiability of reparametrised model #####
### Parameter estimation
set.seed(30)

# Center of starting values for parameter estimation
pouter <- structure(c(log(pars["x1"]), log(pars["x2"]/pars["x1"]), pars["beta1"],
                      pars["beta2"]-pars["beta1"], log(sigma)), names = outerpars)

# Objective function
obj <- normL2(
  data, g_dMod*x_dMod*p_repar, errmodel = error_dMod) + constraintL2(pouter, sigma = 3)

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                cores = 4, fits = 50, samplefun = "runif", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(obj = obj, pars = bestfit, whichPar = theta, alpha = 0.05,
                    cores = 4)
plotProfile(profiles, mode == "data")


################################################################################
######################### Logistic growth model ################################
################################################################################

##### Generate logistic growth model for dMod #####
# ODE model
ode_equation <- eqnvec(
  x1 = "lambda1*x1*(1-(x1+x2)/K)", x2 = "lambda2*x2*(1-(x1+x2)/K)")
model_dMod <- odemodel(ode_equation, modelname = "logistic_model") 
x_dMod <- Xs(model_dMod) 

# Observable (for both estimation and simulation)
observable <- eqnvec(y = "log(x1/(x1+x2))-(sigma^2)/2")
g_dMod <- Y(observable, x_dMod, compile = T, attach.input = F)
observable_sim <- eqnvec(y = "x1/(x1+x2)")
g_sim <- Y(observable_sim, x_dMod, compile = T, attach.input = F)

# Error model
error_pars <- eqnvec(y = "sigma")
error_dMod <- Y(
  error_pars, f = c(observable, ode_equation), attach.input = TRUE,  compile = TRUE)

# Parameter transformations
p_dMod <- P(trafo=eqnvec(
  lambda1 = "lambda1", lambda2 = "lambda2", K = "exp(K)", x1 = "exp(x1)", 
  x2 = "exp(x2)", sigma = "exp(sigma)"), condition = "C1")
outerpars <- getParameters(p_dMod)


##### Simulate data for structural identifiability analysis #####
set.seed(59)

# Parameter values
pars <- c(x1 = 10, x2 = 10, lambda1 = 0.1, lambda2 = 0.2, K = 1e3) 
sigma <- 0.2

# Time range
times <- seq(0, 50, 1)

# Simulation with dMod 
out <- (g_sim*x_dMod)(times, pars)
data_sim <- do.call(rbind, replicate(20, out[[1]], simplify = FALSE))
data_sim <- data_sim[order(data_sim[,1]),]

# Multiplying simulated observable with log-normally distributed measurement error
data_sim[,2] <- data_sim[,2]*rlnorm(nrow(data_sim), (-(sigma^2)/2), sigma)
data <- data.frame(name = rep("y", nrow(data_sim)), time = data_sim[,1], 
                   value = log(data_sim[,2]), sigma = rep(NA,nrow(data_sim)))
data <- datalist(C1 = data)


##### Identifiability of full model #####
### Parameter estimation
set.seed(100)

# Center of starting values for parameter estimation
pouter <- structure(
  c(pars["lambda1"], pars["lambda2"], log(pars["K"]), log(pars["x1"]), 
    log(pars["x2"]), log(sigma)), names = outerpars)

# Objective function
obj <- normL2(
  data, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(pouter, sigma = 10)

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(obj, pouter, studyname = "multiple_fits", rinit = 0.1, rmax = 10, 
                cores = 4, fits = 50, samplefun = "rnorm", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run 
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(obj = obj, pars = bestfit, whichPar = names(bestfit), 
                    alpha = 0.05, cores = 4)
plotProfile(profiles, mode == "data")


##### Identifiability with fixed x1(0) #####
### Parameter estimation
set.seed(100)

# Specify fixed parameter x1
fixed_par <- c("x1") 

# Center of starting values for parameter estimation
pouter <- structure(
  c(pars["lambda1"], pars["lambda2"], log(pars["K"]), log(pars["x1"]), 
    log(pars["x2"]), log(sigma)), names = outerpars)

# Objective function
obj <- normL2(
  data, g_dMod*x_dMod*p_dMod, errmodel = error_dMod) + constraintL2(pouter, sigma = 10) 

# Multiple parameter estimations with randomly chosen starting values
fits <- mstrust(obj, pouter[!names(pouter) %in% fixed_par], fixed = pouter[fixed_par],  
                studyname = "multiple_fits", rinit = 0.1, rmax = 10, cores = 4,
                fits = 50, samplefun = "rnorm", iterlim = 500)
summary(fits)

# Extract estimated parameters of best estimation run 
bestfit <- as.parvec(as.parframe(fits))

### Profile likelihoods
profiles <- profile(obj = obj, pars = bestfit, fixed = pouter[fixed_par],
                    whichPar = names(bestfit), alpha = 0.05, cores = 4)
plotProfile(profiles, mode == "data")

