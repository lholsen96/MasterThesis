FisherScoring_normal = function(X,
                                y,
                                beta_hat,
                                sigma_squared,
                                lambda,
                                num_iter,
                                print=FALSE) {
  
  # X: matrix of data
  # y: vector of the response
  # beta_hat: The value where we want to start the fisher
  # sigma_squared: the variance of var(y_i). 
  # lambda: penalty/smoothing parameter
  # num_iter: number of iterations
  # print: whether or not we print the results
  
  # Get the number of observations
  n = dim(X)[1]
  # Get the number of covariates (intercept is not included)
  p = dim(X)[2]
  
  # Create the design matrix with intercept. Add columns with 1s
  Z = matrix(c(rep(1, n), X), nrow = n)
  
  # Create the block diagonal penatlty matrix, with 1x1 block 0
  # and a pxp block with 1s one diagonal. Do not penalize intercept
  P = diag(p+1)
  P[1,1] = 0
  
  # Calculate eta_hat based on beta_hat
  eta_hat = Z %*% beta_hat
  
  # Since we are working with normal data (h(x) = x) we get that
  # the current estimate for the mean mu is
  mu_hat = eta_hat
  
  # Matrix versions
  beta_hat_matrix = matrix(nrow = (p+1), ncol = num_iter+1)
  beta_hat_matrix[,1] = beta_hat
  mu_hat_matrix = matrix(nrow = n, ncol = num_iter+1)
  mu_hat_matrix[,1] = mu_hat
  
  
  # Calculate the penalized log likelihood
  # Should increase (or stay put) for each iteration since
  # we update beta with respsect to maximize penalized logl. 
  loglike = -n/2 * log(2*pi*sigma_squared) - sum((y-mu_hat)^2)/(2*sigma_squared)
  loglike_penalized = loglike - lambda/2 * t(beta_hat) %*% P %*% beta_hat
  
  # Print the values if we are asked to
  if (print) {
    cat(sprintf("Iteration: %3d", 0))
    cat(sprintf("  logl_pen:"), c(round(loglike_penalized, 4)))
    cat(sprintf("  Beta:"), c(round(beta_hat, 3), sprintf("\n")))
  }
  
  ### Have to create the matrices W, D and SIGMA
  ### These should be in the for loop but since they
  ### don't depend on beta slash eta for normal data
  ### I moved them outside to remove redundant calculations.
  # SIGMA is diag(var(y_1), var(y_2), ..., var(y_n))
  # and we know that var(y_i) = sigma_squared = 1 
  # could have estimated it if it was uknown.
  SIGMA = diag(n)*sigma_squared
  
  # D is a diagonal matrix with dh(eta_i)/(d*eta_i) on the diagonal
  # for the normal distirbution we have that h(x) = x, identity
  # function, so the derivative of it is 1 on the diagonal.
  # Could also think that we calcultate, for i = 1,2,...,n
  # dh(eta_i)/(d*eta) = [dh(eta_i)/(d*eta_1), ..., dh(eta_i)/(d*eta_n)]
  # and only take the ith component of the ith row. Get the same.
  D = diag(n)
  
  # W is defined to be the product of D * SIGMA^-1 * D
  W = D %*% solve(SIGMA) %*% D
  
  for (iter in 1:num_iter){
    # The score is given by Z*W*D^-1*(y-mu)
    score = t(Z) %*% W %*% solve(D) %*% (y-mu_hat)
    score_pen = score - lambda * P %*% beta_hat
    
    # Calculate the Fisher matrix
    # They write X in the paper, but they mean Z
    # I.e. they want the design matrix with the intercept
    fisher = t(Z) %*% W %*% Z
    
    # Add the penalization to the fisher matrix
    fisher_pen = fisher + lambda*P
    
    # Find the update change in beta
    update_term = solve(fisher_pen) %*% score_pen
    
    # Calculate the new estimates for beta an mu
    beta_hat = beta_hat + update_term
    mu_hat = Z %*% beta_hat
    
    # Matrix versions
    beta_hat_matrix[, iter+1] = beta_hat
    mu_hat_matrix[, iter+1] = mu_hat
    
    # Since we are maximizing the likelihood, we expect to see
    # that the likelihood increase for each iteration
    loglike = -n/2 * log(2*pi*sigma_squared) - sum((y-mu_hat)^2)/(2*sigma_squared)
    loglike_penalized = loglike - lambda/2 * t(beta_hat) %*% P %*% beta_hat
    
    # We print every 10% percent iteration
    if (print && !(iter %% (num_iter %/% 10))) {
      cat(sprintf("Iteration: %3d", iter))
      cat(sprintf("  logl_pen:"), c(round(loglike_penalized,4)))
      cat(sprintf("  Beta:"), c(round(beta_hat, 3), sprintf("\n")))
    }
  }
  
  # Return the estimates for beta and mu
  return(list(beta_hat = beta_hat, mu_hat = mu_hat, beta_hat_matrix=beta_hat_matrix,
              mu_hat_matrix=mu_hat_matrix))
}
