GenPartBoostR_normal = function(X,
                                y,
                                beta_hat,
                                sigma_squared,
                                lambda,
                                num_iter,
                                num_iter_init,
                                V_compulsory,
                                print=FALSE) {
  # X: matrix of data
  # y: vector of the response
  # beta_hat: The value where we want to start the fisher
  # sigma_squared: the variance of var(y_i). 
  # lambda: penalty/smoothing parameter
  # num_iter: number of iterations. Assume num_iter > 1.
  # num_iter_init: number of iterations in the initialization of beta_0_hat
  # V_compulsory: list of covariates that always gets updated. Domain: {0,1,..,p}, where 0 is intercept. 
  # print: whether or not we print the results
  
  ##### Step 0: Setup
  n = dim(X)[1]
  p = dim(X)[2] 
  
  # Create Z that is equal to X, but with the first column beeing 1, intercept  
  Z = matrix(c(rep(1, n), X), nrow = n)
  
  # Check for invalid compulsory variables 
  if(length(V_compulsory) != 0 && (range(V_compulsory)[1] < 0 ||  range(V_compulsory)[2] > p)) {
    stop("Invalid compulsory variable indices.")
  }
  
  # Ensure no duplicates in V_compulsory
  V_compulsory = sort(unique(c(V_compulsory)))
  
  # Check if we are in the extreme case of GenBoostR
  GenBoostR_algo = length(V_compulsory) == p+1
  
  # Calculate the number of parameters to be updated in each iteration
  num_param_update = min(length(V_compulsory) + 1, p+1)
  
  # List over indices not in the compulsory list
  if (length(V_compulsory) == 0) { # There are no compulsory variables
    V_not_compulsory = 0:p                  
  } else if (GenBoostR_algo) { # All the variables are compulsory -> GenBoostR
    V_not_compulsory = c()                  
  } else { # There are some compulsory variables
    V_not_compulsory = (0:p)[-(V_compulsory+1)] 
  }
  
  # Create the possible subsets
  subsets = list()
  if (GenBoostR_algo) { # The only subset is the whole set
    subsets[[1]] = V_compulsory
  } else {
    for (j in 1:length(V_not_compulsory)) { # Combine compulsory and potential
      subsets[[j]] = c(V_compulsory, V_not_compulsory[j])
    }
  }
  
  # Create matrices to store the updates 
  beta_hat_matrix = matrix(nrow = (p+1), ncol = (num_iter + 1))
  mu_hat_matrix = matrix(nrow = n, ncol = (num_iter + 1))
  
  # Variable to store the total acumulated Beta values
  beta_hat_tot = NULL
  
  # matrix to store which variable we update
  update_trace = matrix(rep(NA, num_iter*num_param_update), nrow = num_iter)
  
  ##### Step 1: Initialization
  # Chose a random value for beta_hat_intercept
  beta_hat_intercept = 0
  init = mean(y) # Only holds for normal distrbution
  
  # We know create the beta_0_hat vector consisting of the fitted
  # intercept and the rest of the parameters are set to zero
  beta_hat_matrix[,1] = c(init, rep(0, p))
  
  # Can now calculate the linear product, eta
  eta_hat = Z %*% beta_hat_matrix[,1]
  
  # Since we are deling with normal data
  mu_hat_matrix[,1] = Z %*% beta_hat_matrix[,1]
  
  ##### Step 2: Iterations
  for (m in 1:num_iter) {
    #print(c("Working on iteration:", m))
    
    # Array to store the l2 loss for each of the candidate sets
    L2_loss_array = c()
    
    # iterate over all subsets we need to evaluate
    for (j in 1:length(subsets)) {
      # Create the partial design matrix
      X_Vm_j = as.matrix(Z[ ,subsets[[j]]+1])
    
      # Since we can update the intercept, we use Z instead of X.
      # Note that the paper use X since they have defined it to cointain
      # the intercept, while we have not. 
      # Since we include the intercept and R counts from 1 we have to add
      # 1 to each component. I.e. to get covariate 2 we need to take column
      # number 2+1 = 3.

      # Check if it is beta_0 (intercept) we are updating.
      beta_0_included = 0 %in% subsets[[j]] 
      

      # Compute the penalized Fisher scoring update
      # Get the number of observations and covariates
      nn = dim(X_Vm_j)[1]
      pp = dim(X_Vm_j)[2]
        
      # Create the penalization matrix (Can be 1x1 if no compulsory variables)
      P = diag(pp)
      if (beta_0_included) {
        # Do not penalize the intercept
        P[1,1] = 0
      }
        
      # Get the matrices. See Fisher_algorithm_normal_data() for explanation.
      # ONLY FOR NORMAL DATA!!!
      SIGMA = diag(nn)*sigma_squared
      D = diag(nn)
      W = D %*% solve(SIGMA) %*% D
        
      # Calculate the score. Penalty disappear since -lambda*P*beta = 0
      # by defavult since we assume we start at beta == 0 in this update.
      score_p_vm_j = t(X_Vm_j) %*% W %*% solve(D) %*% (y-mu_hat_matrix[,m])
      
      # Calculate the fisher matrix witht the penalty
      fisher_p_vm_j = t(X_Vm_j) %*% W %*% X_Vm_j
      fisher_p_vm_j_penalty = fisher_p_vm_j + lambda * P
        
      # Calculate the update for beta through one step of Fisher scoring
      beta_hat_temp_update = solve(fisher_p_vm_j_penalty) %*% score_p_vm_j
  
      # Get the new temp possible beta parameter array
      beta_hat_temp = beta_hat_matrix[,m]
      beta_hat_temp[subsets[[j]]+1] = beta_hat_temp[subsets[[j]]+1] + beta_hat_temp_update
      
      # Compute the L2 loss
      eta_hat_m_j = Z %*% beta_hat_temp
      mu_hat_m_j = eta_hat_m_j
      
      # Compute the L2 loss since we are working with normal data. Otherwise we use deviance
      L2_loss_temp = sum((mu_hat_m_j - y)^2)
      
      # Add this possible loss to the loss array
      L2_loss_array = c(L2_loss_array, L2_loss_temp)
    }
    
    # Need to find the update that yields the best fit. min L2 loss
    min_idx = which.min(L2_loss_array)
    
    # Partial matrix for best update
    update_idx = subsets[[min_idx]]  # Find the variables with best update
    update_trace[m, ] = update_idx   # Record which variables that got updated
    Z_m = as.matrix(Z[, update_idx+1]) # Get the best partial design matrix
    
    # Check if it is beta_0 (intercept) we are updating.
    beta_0_included = 0 %in% update_idx
    
    # Conduct the penalized Fisher scoring update
    # Get the number of observations and covariates
    nn = dim(Z_m)[1]
    pp = dim(Z_m)[2]
    
    # Create the penalization matrix (Can be 1x1 if no compulsory variables)
    P = diag(pp)
    if (beta_0_included) {
      # Do not penalize the intercept
      P[1,1] = 0
    }
    
    # Get the matrices. See Fisher_algorithm_normal_data() for explanation.
    # ONLY FOR NORMAL DATA!!!
    SIGMA = diag(nn)*sigma_squared
    D = diag(nn)
    W = D %*% solve(SIGMA) %*% D
    
    # Calculate the score. Penalty disappear since -lambda*P*beta = 0
    # by defavult since we assume we start at beta == 0 in this update.
    score_p_vm_j = t(Z_m) %*% W %*% solve(D) %*% (y-mu_hat_matrix[,m])
    
    # Calculate the fisher matrix witht the penalty
    fisher_p_vm_j = t(Z_m) %*% W %*% Z_m
    fisher_p_vm_j_penalty = fisher_p_vm_j + lambda * P
    
    # Calculate the update for beta through one step of Fisher scoring
    best_update = solve(fisher_p_vm_j_penalty) %*% score_p_vm_j
    
    # New estimates
    beta_hat_matrix[,m+1] = beta_hat_matrix[,m]
    beta_hat_matrix[update_idx+1,m+1] = beta_hat_matrix[update_idx+1, m+1] + best_update
    mu_hat_matrix[,m+1] = Z %*% beta_hat_matrix[,m+1]
    
    # We print every 10% percent iteration
    if (print && !(m %% (num_iter %/% 10))) {
      cat(sprintf("Iteration: %3d", m))
      cat(sprintf("  MSE:"), mean((y - round(mu_hat_matrix[,m+1],3))^2))
      cat(sprintf("  Beta:"), c(round(beta_hat_matrix[,m+1], 4), sprintf("\n")))
    }
  }
  return(list(beta_hat_matrix = beta_hat_matrix, mu_hat_matrix = mu_hat_matrix, 
              update_trace = update_trace))
}



