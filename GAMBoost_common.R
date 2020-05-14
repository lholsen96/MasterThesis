##### File that contains common functions for both 
##### GAMBoost_splines and GAMBoost_stumps
### create_data_combined
### create_total_design_matrix
### GAMBoost_plot_predictor_contribution

create_data_combined = function(n, seed_number = NA, cn = NA, stnr = NA, 
                                type = c("smooth", "step", 'advancedstep'),
                                family = gaussian()) {
  # Function that creates the data in accordance to 
  # the setup in chapter 4 in thesis.
  
  # n: number of observations
  # seed_number: the seed for the random sampling generator
  # cn: constant in eta that alters the signal to noise ratio. Set to 3 in paper.
  # stnr: desired signal to noise ratio. Either cn is provided or this.
  # type: If the functions should be smooth or stepwise
  # family: The distribution of the response; gaussian, binomial or poisson.
  
  # Check valid type input
  type = match.arg(type)
  
  if ((is.na(cn) & is.na(stnr)) | (!is.na(cn) & !is.na(stnr))) {
    stop("Either cn or stnr needs to be provided, and not both.")
  }
  
  # See if cn was provided
  cn_recieved = ifelse(is.na(cn), FALSE, TRUE)
  
  # Change from characters to function 
  if (is.character(family)) {
    family = switch(family, gaussian = gaussian(), binomial = binomial(), poisson = poisson())
  }
  
  # Set the seed for reproducibility
  if (!is.na(seed_number)){
    set.seed(seed_number)
  }
  
  # Create needed auxillary functions
  if (type == 'step') {
    # Create the stepfunction f(x) which is -1 if x<= 0 and 1 if x>0.
    sf = stepfun(c(0), c(-1, 1), right=TRUE)
  } else if (type == 'advancedstep') {
    # Create the piecewise constant functions
    step1 = approxfun(x=seq(-1,1,0.4), y = c(-2,-1,0,1,2,2), method = "constant")
    step2 = approxfun(x=seq(-1, 1, length.out = 7), y = c(1,-1,2,-2,3,-3,-3), method = "constant")
    step3 = approxfun(x=seq(-1,1,0.4), y = c(2,-2,0,-2,2,2), method = "constant")
    step4 = approxfun(x=seq(-1,1, length.out = 7), y = c(1,0,-1, -1, 0, 1,1), method = "constant")
    step5 = approxfun(x=seq(-1,1,2), y = c(0,0), method = "constant")
  }
  
  # Create the observations.
  X = matrix(runif(n*5, -1, 1), nrow = n, ncol = 5)
  
  # Find the cn which yields the desired signal to noise ratio
  if (!is.na(stnr)) {
    if (type == 'smooth') {
      optim_func_smooth = function(cc) {
        eta_temp = cc*(-0.7 + X[,1] + 2*X[,3]^2 + sin(5*X[,5])) 
        mu_temp = family$linkinv(eta_temp)
        (stnr - sum((mu_temp - mean(mu_temp))^2) / sum(family$variance(mu_temp)))^2
      }
      cn = optimize(optim_func_smooth, c(0,100))$minimum
    } else if (type == 'step') {
      optim_func_step = function(cc) {
        # Create the stepfunction f(x) which is -1 if x<= 0 and 1 if x>0.
        #sf = stepfun(c(0), c(-1, 1), right=TRUE)
        # Create the linear predictor eta.
        eta_temp = cc*(0.5*sf(X[,1]) + 0.25*sf(X[,3]) + sf(X[,5])) 
        mu_temp = family$linkinv(eta_temp)
        (stnr - sum((mu_temp - mean(mu_temp))^2) / sum(family$variance(mu_temp)))^2
      }
      cn = optimize(optim_func_step, c(0,100))$minimum
    } else if (type == 'advancedstep') {
      optim_func_advstep = function(cc) {
        eta_temp = cc*(step1(X[,1]) + step2(X[,2]) + step3(X[,3]) + step4(X[,4]) + step5(X[,5]))
        mu_temp = family$linkinv(eta_temp)
        (stnr - sum((mu_temp - mean(mu_temp))^2) / sum(family$variance(mu_temp)))^2
      }
      cn = optimize(optim_func_advstep, c(0,100))$minimum
    }
  }
  
  if (type == 'smooth') {
    # Create the linear predictor eta.
    eta = cn*(-0.7 + X[,1] + 2*X[,3]^2 + sin(5*X[,5]))    
  } else if (type == 'step') {
    # Create the stepfunction f(x) which is -1 if x<= 0 and 1 if x>0.
    # sf = stepfun(c(0), c(-1, 1), right=TRUE)
    # Create the linear predictor eta.
    eta = cn*(0.5*sf(X[,1]) + 0.25*sf(X[,3]) + sf(X[,5])) 
  } else if (type == 'advancedstep') {
    # step1 = approxfun(x=seq(-1,1,0.4), y = c(-2,-1,0,1,2,2), method = "constant")
    # step2 = approxfun(x=seq(-1, 1, length.out = 7), y = c(1,-1,2,-2,3,-3,-3), method = "constant")
    # step3 = approxfun(x=seq(-1,1,0.4), y = c(2,-2,0,-2,2,2), method = "constant")
    # step4 = approxfun(x=seq(-1,1, length.out = 7), y = c(1,0.25,-0.5, -0.5, 0.25, 1,1), method = "constant")
    # step5 = approxfun(x=seq(-1,1,2), y = c(0,0), method = "constant")
    
    # Create the linear predictor eta.
    eta = cn*(step1(X[,1]) + step2(X[,2]) + step3(X[,3]) + step4(X[,4]) + step5(X[,5]))
  }
  
  # Calculate the expected response.
  mu = family$linkinv(eta)
  
  # If cn recieved, then compute SNTR.
  if (cn_recieved) {
    stnr = sum((mu - mean(mu))^2) / sum(family$variance(mu))
  }
  
  # Randomly generate the response based on family and mu.
  if (family$family == "gaussian") {
    y = rnorm(n, mean=mu, sd=1)
  }
  if (family$family == "binomial") {
    y = rbinom(n, size = 1, prob = mu)
  }
  if (family$family == "poisson") {
    y = rpois(n, lambda = mu)
  }
  return(list("X" = X, "eta" = eta, "mu" = mu, "y" = y, "type" = type, "family" = family,
              "cn" = cn, 'stnr' = stnr, "seed_number" = seed_number))
}

create_total_design_matrix = function(X) {
  # X is a matrix containing the traing data
  # Here we create the design matrix assuming 
  # the stump/tree paradigm. I.e. create a matrix
  # where we assume that there exists a split at all 
  # observed values in the training data
  
  # Get the number of observations and number of parameters
  n = nrow(X)
  p = ncol(X)
  
  # Create a matrix on the right form
  total_design_matrix = matrix(NA, nrow = n, ncol = 2*n*p)#, sparse = TRUE)
  
  # Iterate over the p parameters
  for (j in 1:p){
    
    # Extract all the x_ij, for i = 1,2,...,n.
    xx_param_j = X[,j]
    
    # Calculate the offset of where to put the new values into the 
    # total_design_matrix
    offset = (j-1)*(2*n)
    
    # Iterate over the n observations if the j'th parameter
    for (i in 1:n) {
      
      # Find the observations that are smaller than x_ij.
      smaller_than_x_ij = xx_param_j <= xx_param_j[i]
      
      # Apply the indiacator functions
      total_design_matrix[, offset + 2*i - 1] = smaller_than_x_ij
      total_design_matrix[, offset + 2*i] = !smaller_than_x_ij
    }
  }
  return(total_design_matrix*1) # Multiply with one to get 0-1 values.
}

