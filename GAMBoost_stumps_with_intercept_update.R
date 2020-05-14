Fisher_update_intercept = function(y,
                                   current_intercept,
                                   current_eta_hat,
                                   family) {
  
  # y: The observed response
  # lambda: The amount of penalization
  # beta_hat: The current estimate for the parameter vector
  # eta_hat: The current estimate for eta, linear combination of X and beta
  # mu_hat: The current estimate for mean mu, mu_hat = h(eta_hat)
  # sigma_squared: The variance of the response, var(y_i) = sigma_squared
  # beta_hat_includes_beta_0_first: Bool to show if beta_0 is recieved.
  #                                 Need to know this since it is not
  #                                 supposed to be penlized.
  
  
  # Note the since we are in the special case of standard normal data we
  # have that mu_hat = eta_hat and that W and D is always the identity matrix.
  # So this function needs to be generalized to take into account other
  # distributions.
  
  
  
  # Get the number of observations and covariates
  n = length(current_eta_hat)
  p = 1
  X_mat = matrix(1, nrow = n)
  
  # Get current_mu_hat
  current_mu_hat = family$linkinv(current_eta_hat)
  
  # SIGMA is a diagonal matrix with the variance of y_i on the i'th diagonal entry.
  SIGMA = family$variance(current_mu_hat)
  SIGMA_inv = (SIGMA)^(-1)
  
  # D is a diagonal matrix with dh(eta_i)/(d*eta_i) on the i'th entry.
  D = family$mu.eta(current_eta_hat)
  
  # Calculate W = D * SIGMA^(-1) * D
  W = D * SIGMA_inv * D # Simpler W = D
  
  # Calculate the score.
  score = t(X_mat) %*% diag(W) %*% solve(diag(D)) %*% (y-current_mu_hat)
  score = sum(y-current_mu_hat)
  
  # Calculate the fisher matrix
  fisher = t(X_mat) %*% diag(W) %*% X_mat
  fisher = sum(W)
  
  # Calculate the update for beta through one step of Fisher scoring
  intercept_update = solve(fisher) %*% score
  intercept_update = score / fisher
  
  return(intercept_update)
}




GAMBoost_stumps_with_intercept_update = function(X, y, num_iter, lambda, family = gaussian(),
                                                 print_msg = FALSE, tiny_return = FALSE){
  # Conduct the GAMBoost procedure with P-stumps.
  # Returns the parameter build-up, Hat matrix and so on
  # More memory efficient. Does not save the whole hat matrices,
  # but only the "variance" diag(Q * Q^T). Return only the final hat matrix
  
  # X: The matrix of the observed values
  # y: The responses asscosiated with the explanatory varaibles
  # num_iter: the number of iterations we perform. Bias-variance trade-off.
  # num_iter_init: number of iterations of step 1. Intercept adjustment/fitting.
  # lambda: the amount of penalization we apply in Fisher scoring
  # sigma_squared: the varaince of the reponse y, i.e. sigma_squared = var(y)
  #                assumed to be one.
  # family: The distribution of the response. Code only supports guassian.
  #           Need to rewrite funtion to make it work in non-gaussian case.
  # print_msg; If we want information during the fitting.
  # tiny_return: If we only want the essential return. Memory efficient.
  
  # # Check valid family input
  # if (is.character(family)) {
  #   family = switch(family, gaussian = gaussian(), binomial = binomial(), poisson = poisson())
  # } else {
  #   if (is.function(family))
  #     family = family()
  # }
  # family = match.arg(family)
  # canonical.link = switch(family$family, gaussian = "identity",
  #                          binomial = "logit", poisson = "log")
  # if (is.null(canonical.link) || family$link != canonical.link) {
  #   warning(paste("GAMBoost expects the canonical link for family '",
  #                 family$family, "' and does not honor dispersion parameters.\n",
  #                 sep = ""))
  # }
  
  ### Step 0: Setup. Create vectors, matrices, et cetera to store quantities.
  {
    # Dimensions
    n = nrow(X)
    p = ncol(X)
    
    # Weights used in calculating deviance. Assumed to be identical.
    weights = rep(1, n)
    
    # The penalization term used in penalized Fisher scoring
    pen_term = lambda * matrix(c(1,-1,-1,1), ncol = 2)
    
    # Create the total transformed design matrix
    Z_matrix = create_total_design_matrix(X)
    
    # Matrices to store the model parameters, predictor eta, response mu = h(eta) = g^(-1)(eta).
    alpha_hat_mat = matrix(0, ncol = 2*n*p, nrow = num_iter + 1)
    eta_hat_mat   = matrix(0, ncol = n, nrow = num_iter + 1)
    mu_hat_mat    = matrix(0, ncol = n, nrow = num_iter + 1)
    
    # Matix to store the aprroximate response obtained the hat mat.
    mu_hat_mat_approx = matrix(0, ncol = n, nrow = num_iter + 1)
    
    # Matrix to store which split variable and observation
    split_points = matrix(NA, nrow = num_iter, ncol = 2, dimnames = list(NULL, c("Var", "Obs")))
    
    # Arrays to store deviance, trace/DF, AIC and BIC (approximated from hat mat)
    deviance = rep(NA, num_iter+1)
    trace    = rep(NA, num_iter+1)
    AIC      = rep(NA, num_iter+1)
    BIC      = rep(NA, num_iter+1)
    
    # Array to store the updated intercept
    intercept_array = rep(NA, num_iter+1)
    
    # Array to store the deviance obtained for the aprroximate response from the hat mat.
    deviance_approx = rep(NA, num_iter+1)
    
    # Array to store the difference between the predicted response from algorithm
    # and the approximate predicted response based on the hat matrix.
    deviance_difference = rep(NA, num_iter+1)
    
    # List to store the iterative approximate hat matrices
    hat_mat = list()
    
    # Variables to store the matrices in the calculation for the hat matrix
    M_mat_new = NA
    M_mat_old = 1/n * matrix(1, nrow = n, ncol = n)
    
    # Used to iteratively caclulate the new term in the hat matrix efficiently
    cumulative_product = diag(n)
    
    # Matrix to store the sum-to-zero constraints
    stz = matrix(NA, ncol = p, nrow = num_iter+1)
    
    # List to store lists containing the Q_MJ matrices. p lists with n matrices
    # of size nxn. Used to derive the standard error bands. Not memory effiecient.
    Q_mj_list = list()
    for (i in 1:p) {
      Q_mj_list[[i]] = list()
      names(Q_mj_list)[i] = paste("var", i, sep = "")
    }
    
    # List to store p matrices of dim nxnum_iter which will store the diagonal
    # elements of the Q_mj matrices. These values are the approximate variance
    # for each of the n observations in each iteration.
    obsvar = list()
    for (i in 1:p) {
      obsvar[[i]] = matrix(NA, nrow = n, ncol = num_iter+1)
      names(obsvar)[i] = paste("var", i, sep = "")
    }
  }
  
  
  ### Step 1: Initialization
  {
    # Fit the intercept model
    eta_0 = rep(family$linkfun(mean(y)), n)
    intercept_array[1] = eta_0[1]
    
    # Insert the initial predictor and reponse.
    eta_hat_mat[1, ] = eta_0
    mu_hat_mat[1, ]  = family$linkinv(eta_hat_mat[1, ])
    
    # Calculate the first hat matrix, since mu = \bar{y}
    hat_mat[[1]] = 1/n * matrix(1, nrow = n, ncol = n)
    
    # Calculate the initial quantites
    deviance[1] = sum(family$dev.resids(y, mu_hat_mat[1, ], weights))
    trace[1]    = sum(diag(hat_mat[[1]]))
    AIC[1]      = deviance[1] + 2*(trace[1])      # Binder adds 1 to trace if Guassian. Not sure why.
    BIC[1]      = deviance[1] + log(n)*(trace[1]) # Binder adds 1 to trace if Guassian. Not sure why.
    # Could calculate AIC corrected for Gaussian, but have not seen Binders version before.
    
    # Compute the deviance between the response and approiximate reponse.
    deviance_difference[1] = sum(family$dev.resids(mu_hat_mat[1, ], hat_mat[[1]] %*% y, weights))
    
    # Compute the approximate responses.
    mu_hat_mat_approx[1, ] = hat_mat[[1]] %*% y
    deviance_approx[1]     = sum(family$dev.resids(y, mu_hat_mat_approx[1, ], weights))
    
    # Do not update the f_j functions, hence they are zero.
    stz[1, ] = rep(0, p)
    
  }
  # Small printout to the user
  if (print_msg) {
    cat(sprintf('Iter: %-3d Var: %-2d  Dev: %8.3f  DF: %5.2f  AIC: %7.2f  BIC: %7.2f  DevDif: %10.3g  stz = %g\n',
                0, NA, deviance[1], trace[1], AIC[1], BIC[1], deviance_difference[1], NA))
  }
  
  
  ### Step 2 (Model fit):
  # Iterate over the number of iterations
  for (l in 1:num_iter){
    
    # Matrix to store the L2 loss of the p*n potential updates.
    deviance_loss_matrix = matrix(NA, nrow = p, ncol = n)
    
    ### In Fisher scoring we need the matrices; D, SIGMA and W. These are
    ### only iteration specific and do not depent on split variable and value.
    ### Thus, we save computational time by computing them here.
    {
      # SIGMA is a diagonal matrix with the variance of y_i on the i'th diagonal entry.
      SIGMA = family$variance(mu_hat_mat[l,])
      SIGMA_inv = (SIGMA)^(-1)
      
      # D is a diagonal matrix with dh(eta_i)/(d*eta_i) on the i'th entry.
      D = family$mu.eta(eta_hat_mat[l, ])
      
      # Calculate W = D * SIGMA^(-1) * D
      W = D * SIGMA_inv * D # Simpler W = D
    }
    
    # Iterate over the number of parameters to find split variable.
    for (s in 1:p) {
      
      # Iterate over the number of observations to find split value.
      for (d in 1:n) {
        
        ### Calculate the potential update
        {
          # Extract the partial transformed design matrix Z_sd of dimension nx2.
          offset = 2*n*(s-1) + 2*d - 1
          Z_sd = Z_matrix[,offset:(offset+1)]
          
          # Calculate param update based on one iteration of penalized fisher scoring starting at zero.
          F_penalized = t(Z_sd * W) %*% Z_sd + pen_term
          score = t(Z_sd * (W / D)) %*% (y-mu_hat_mat[l,]) # Simpler score = t(Z_sd) %*% (y-mu_hat_mat[l,])
          alpha_hat_update = solve(F_penalized) %*% score
          
          # Update the parameter vector
          alpha_hat_temp = alpha_hat_mat[l, ]
          alpha_hat_temp[offset:(offset+1)] = alpha_hat_temp[offset:(offset+1)] + alpha_hat_update
        }
        
        ### Evaluate the potential update based on deviance
        {
          # Update the predictor and include the intercept
          eta_hat_temp = Z_matrix %*% alpha_hat_temp + rep(intercept_array[l], n)
          
          # Calculate the predicted response mu = h(eta) = g^(-1)(eta)
          mu_hat_temp = family$linkinv(eta_hat_temp)
          
          # Calculate the deviance. L2 loss in the case of Gaussian
          deviance_loss_matrix[s,d] = sum(family$dev.resids(y, mu_hat_temp, weights))
        }
      }
    }
    
    ### Find best split variable and value
    {
      best_indices = arrayInd(which.min(deviance_loss_matrix), dim(deviance_loss_matrix))
      s_new = best_indices[1]
      d_new = best_indices[2]
      
      # Add that we use this split point and split variable
      split_points[l, ] = c(s_new, d_new)
    }
    
    ### Calculate the best update
    {
      # Get partial transformed design matrix
      offset = 2*n*(s_new-1) + 2*d_new - 1
      Z_sd = Z_matrix[ ,offset:(offset+1)]
      
      # Calculate update based on one iteration of penalized fisher scoring starting at zero
      F_penalized = t(Z_sd * W) %*% Z_sd + pen_term        # Can exchange W with D or SIGMA
      score = t(Z_sd * (W / D)) %*% (y-mu_hat_mat[l,])     # Simpler score = t(Z_sd) %*% (y-mu_hat_mat[l,])
      alpha_hat_update = solve(F_penalized) %*% score
      
      # Create the new parameter vector
      alpha_hat_mat[l+1, ] = alpha_hat_mat[l, ]
      alpha_hat_mat[l+1, offset:(offset+1)] = alpha_hat_mat[l+1, offset:(offset+1)] + alpha_hat_update
    }
    
    # Compute the intercept update. Will be 0 for Gaussian response
    {
      # Compute the new eta with the previous intercept
      current_eta_hat = Z_matrix %*% alpha_hat_mat[l+1, ] + rep(intercept_array[l], n)
      
      # calculate the response mu = h(eta) = g^(-1)(eta)
      current_mu_hat = family$linkinv(current_eta_hat)
      
      # create the design matrix for the intercept. Column of ones
      X_mat = matrix(1, nrow = n)
      
      # Compute the new SIGMA, D and W matrices.
      SIGMA2 = array(family$variance(current_mu_hat))
      SIGMA_inv2 = (SIGMA2)^(-1)
      
      # D is a diagonal matrix with dh(eta_i)/(d*eta_i) on the i'th entry.
      D2 = array(family$mu.eta(current_eta_hat))
      
      # Calculate W = D * SIGMA^(-1) * D
      W2 = array(D2 * SIGMA_inv2 * D2) # Simpler W = D
      
      # Calculate the score. # Due to the simiplicity, a lot cancels
      score = t(X_mat) %*% diag(W2) %*% solve(diag(D2)) %*% (y-current_mu_hat)
      #score = sum(y-current_mu_hat)
      
      # Calculate the fisher matrix
      fisher = t(X_mat) %*% diag(W2) %*% X_mat
      #fisher = sum(W2)
      
      # Calculate the update for beta through one step of Fisher scoring
      intercept_update = solve(fisher) %*% score
      #intercept_update = score / fisher
      
      # Add the new intercept to the array
      intercept_array[l+1] = intercept_array[l] + intercept_update
    }
    
    # Compute the new eta with the previous intercept
    eta_hat_mat[l+1, ] = Z_matrix %*% alpha_hat_mat[l+1, ] + rep(intercept_array[l+1], n)
    
    # calculate the response mu = h(eta) = g^(-1)(eta)
    mu_hat_mat[l+1, ] = family$linkinv(eta_hat_mat[l+1, ])
    sum(family$dev.resids(y, mu_hat_mat[l+1, ], weights))
    
    
    # Do not need to calcuate the hat matrices and so on.
    if (!tiny_return) {
      ### Calculate the new hat matrix
      {
        # Compute the M matrix, a component in the approximate hat matrix
        H_bar_mat = sqrt(W) * Z_sd %*% solve(t(Z_sd * W) %*% Z_sd + pen_term) %*% t(Z_sd * sqrt(W))
        M_mat_new = diag(sqrt(SIGMA)) %*% H_bar_mat %*% sqrt(diag(SIGMA_inv))
        
        # update the cumulative product
        cumulative_product = (diag(n) - M_mat_old) %*% cumulative_product
        
        # calcualte the new hat matrix
        hat_mat[[l+1]] = hat_mat[[l]] + M_mat_new %*% cumulative_product
        
        # Add the intecept update term to the hat matrix
        hat_mat[[l+1]] = hat_mat[[l+1]] + (
          (family$linkinv(intercept_array[l+1]) - family$linkinv(intercept_array[l]))
          * 1/apply(hat_mat[[l+1]] %*% y, 1, mean) * hat_mat[[l+1]])
        
        # Set the old M matrix to be the new one
        M_mat_old = M_mat_new
      }
      
      ### Calculate the covariate specific eta contribution matrices Q_mj
      {
        # Calculate an auxillary matrix  PS: CAN REMOVE W * SIGMA_inv AS IT BECOMES 1.
        R_mat_s =  Z_sd %*% solve(t(Z_sd * W) %*% Z_sd + pen_term) %*% t(Z_sd * W * SIGMA_inv)
        
        # Get the previous elements
        if (length(Q_mj_list[[s_new]]) == 0) {
          prev_sum = matrix(0, nrow = n, ncol = n)
        } else {
          prev_sum = Q_mj_list[[s_new]][[length(Q_mj_list[[s_new]])]]
        }
        
        # Update the Q_mj list. PS. YES, it is l and not l+1 here. Use the lth iter to predict the l+1th.
        Q_mj_list[[s_new]][[l]] = prev_sum + R_mat_s %*% (diag(n) - hat_mat[[l]])
      }
      
      
      ### calculate the new obsvar (observation (co)variance) which are used to calculate SE bands
      {
        # Copy over the old approximate variance for covariates that were not updated
        for (j_idx in 1:p) {
          if (j_idx == s_new){
            # Eq 11 in GAMBoost paper
            obsvar[[j_idx]][ ,l+1] = diag(Q_mj_list[[s_new]][[l]] %*%
                                            diag(family$variance(mu_hat_mat[l+1, ]))
                                          %*% t(Q_mj_list[[s_new]][[l]]))
          } else {
            obsvar[[j_idx]][ ,l+1] = obsvar[[j_idx]][ ,l]
          }
        }
      }
      
      ### Add new quantites to the different "tracking" vectors
      {
        deviance[l+1] = sum(family$dev.resids(y, mu_hat_mat[l+1, ], weights))
        trace[l+1]    = sum(diag(hat_mat[[l+1]]))
        AIC[l+1]      = deviance[l+1] + 2*(trace[l+1])      # Binder adds 1 to trace if Guassian. Not sure why.
        BIC[l+1]      = deviance[l+1] + log(n)*(trace[l+1]) # Binder adds 1 to trace if Guassian. Not sure why.
        
        # Calculate the predicted response from the hat matrix
        pred_hat = hat_mat[[l+1]] %*% y
        if (family$family == 'binomial') {
          # Ensure values between 0 and 1. Add small epsilon value to no get 'inf'
          pred_hat[pred_hat < 0] = 10e-12
          pred_hat[pred_hat > 1] = 1 - 10e-12
        }
        if (family$family == 'poisson') {
          # Ensure postitive values. Add small epsilon value to no get 'inf'.
          pred_hat[pred_hat < 0] = 10e-12
        }
        deviance_difference[l+1]   = sum(family$dev.resids(mu_hat_mat[l+1, ], pred_hat, weights))
        deviance_approx[l+1]       = sum(family$dev.resids(y, pred_hat, weights))
        mu_hat_mat_approx[l+1, ]   = hat_mat[[l]] %*% y
        
        
        # The sum-to-zero values for the lth iteration
        stz[l+1, ] = stz[l, ]
        stz[l+1, s_new] = sum(Z_matrix[, (2*n*(s_new-1)+1):(2*n*s_new)] %*%
                                alpha_hat_mat[l+1, (2*n*(s_new-1)+1):(2*n*s_new)])
        
        # Compute the approximate responses
        mu_hat_mat_approx[1, ] = hat_mat[[1]] %*% y
        deviance_approx[1]     = sum(family$dev.resids(y, mu_hat_mat_approx[1, ], weights))
      }
    }
    
    # Printout to the user
    if (print_msg) {
      cat(sprintf('Iter: %-3d Var: %-2d  Dev: %8.3f  DF: %5.2f  AIC: %7.2f  BIC: %7.2f  DevDif: %10.3g  stz = %10.4g\n',
                  l, s_new, deviance[l+1], trace[l+1], AIC[l+1], BIC[l+1], deviance_difference[l+1], stz[l+1, s_new]))
    }
  }
  
  
  if (!tiny_return){
    # Create matrix of the diagonal vaules in the hat matrices. Used to calcualte CI predicted response
    list_diag_hat_mat = lapply(1:(num_iter + 1), function(x) diag(hat_mat[[x]]))
    matrix_diag_hat_mat = matrix(unlist(list_diag_hat_mat), ncol = n, byrow = TRUE)
    
  }
  
  
  return_list = list(
    "model_name" = 'GAMBoost_stumps',
    "family" = family,
    "penalty" = lambda,
    "n" = n,
    "p" = p,
    "X" = X,
    "y" = y,
    "num_iter" = num_iter,
    "intercept" = eta_0[1],
    "intercept_array" = intercept_array,
    "alpha" = alpha_hat_mat)
  
  if (!tiny_return) {
    return_list_add = list(
      "eta_hat" = eta_hat_mat,
      "mu_hat" = mu_hat_mat,
      "mu_hat_approx" = mu_hat_mat_approx,
      "deviance" = deviance,
      "deviance_difference" = deviance_difference,
      "deviance_approx" = deviance_approx,
      "AIC" = AIC,
      "BIC" = BIC,
      "trace" = trace,
      "split_points" = split_points,
      "hat_mat_last" = hat_mat[[l+1]],
      "matrix_diag_hat_mat" = matrix_diag_hat_mat,
      #"hat_mat" = hat_mat,
      #"QMJ" = Q_mj_list,
      "obsvar" = obsvar,
      "sum_to_zero" = stz)
    
    # Combine the two return lists
    return_list = c(return_list, return_list_add)
  }
  
  # Return the fitted model
  return(return_list)
}


GAMBoost_stumps_inter_update_predict = function(model, X_new = NA, after_iter = NA,
                                   version = c('response', 'link', 'terms')) {
  
  # model: A fitted GAMBoost_stumps model.
  # X_new: The dataset to be predicted. Either a matrix or a vector.
  #        If not provided we use the training data
  # after_iter: For which iteration of model we should use. We use the
  #             last iteration if no additional information is provided.
  # version: type of prediction to be returned: 
  #          "link" gives prediction at the level of the predictor,
  #          "response" at the response level. 
  #          "terms" returns individual contributions of the components to the predictor.
  #
  
  # Verify the input.
  version = match.arg(version)
  
  # If no iteration is provided we use the last iteration
  if (is.na(after_iter)) {
    after_iter = model$num_iter
  } else{
    # Check if at_step is a legit value
    stopifnot((after_iter >= 0) && (after_iter <= model$num_iter))
  }
  
  # If no new data is provided we predict the training data
  if (!is.matrix(X_new)){
    if (is.na(X_new)) {
      X_new = model$X
    } 
  }
  
  # Get the number of observations and the number of parameters
  n = model$n
  p = model$p
  
  # Get the model parameter vector at after_iter.
  # a vector of size 2np
  alpha_vec = model$alpha[after_iter+1, ]
  
  # Check if we are predicting one or many values. Transform to matrix
  if (is.matrix(X_new)){
    # Check that it has the right number of columns
    stopifnot(ncol(X_new) == p)
    n_new = nrow(X_new)
  } else {
    # Only a signle new observation.
    stopifnot(length(X_new) == p)
    X_new = matrix(X_new, ncol = p)
    n_new = 1
  }
  
  # Create a matrix to store the predicted value for each explanatory variable
  # for each of the observations.
  values_mat = matrix(NA, nrow = n_new, ncol = p)
  
  # Iterate over all the observations
  for (i in 1:n_new) {
    
    # Extract the current new value. Call it x_i
    x_new_value = X_new[i,]
    
    # For each of the p predictors we check if x_ij is smaller than or equal
    # to the observations of the j'th variable in the training data.
    smallerequal = matrix(x_new_value <= c(t(model$X)), byrow=T, ncol = p)
    
    # For each of the p predictors we check if x_ij is larger than 
    # to the observations of the j'th variable in the training data.
    # NB. SAME AS !smallerequal
    larger = matrix(x_new_value  > c(t(model$X)), byrow=T, ncol = p)
    
    # Combine the two matrices.
    # such that we get rows: smallerequal[1, ], larger[1,], smallerequal[2, ], larger[2,], ....
    # Becomes a matrix of shape 2n times p.
    design_test = rbind(smallerequal, larger)[order(rep(1:nrow(smallerequal), 2)), ]
    
    # Calculat the new values for each of the predictors
    # We split alpha_vec into p parts of length 2n and then multiply
    # these with the transformed design matrix for each of the p
    # coavriates and get the eta contribution from each covariate.
    values_mat[i, ] = apply(alpha_vec * design_test, 2, sum)
  }
  
  # version == 'terms'. Return f_1, f_2, ..., f_p values.
  return_object = values_mat
  
  if (version == 'link' | version == 'response') {
    # Additive model so we can sum the effect of the explanatory variables together
    eta_sum = apply(values_mat, 1, sum)
    
    # Add the intercept. We have not the total eta
    return_object = eta_sum + rep(model$intercept_array[after_iter + 1], length(eta_sum))
    
    # If we want to return the response mu.
    if (version == 'response'){
      # apply mu = h(eta) = g^(-1)(eta)
      return_object = model$family$linkinv( return_object)
    }
  } 
  return(return_object)
}
