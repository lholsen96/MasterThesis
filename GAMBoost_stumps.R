##### File containint the GAMBoost_stumps algorithm 
##### and additionally helping functions
### GAMBoost_stumps
###
###

require(graphics)
nice_colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)


GAMBoost_stumps = function(X, y, num_iter, lambda, family = gaussian(),
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
  # }
  # else {
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
    Z_matrix = GAMBoost_stumps_create_total_design_matrix(X)
    
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
          eta_hat_temp = Z_matrix %*% alpha_hat_temp + eta_0
          
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
    
    # Compute the new eta with ther intercept
    eta_hat_mat[l+1, ] = Z_matrix %*% alpha_hat_mat[l+1, ] + eta_0
    
    # calculate the response mu = h(eta) = g^(-1)(eta)
    mu_hat_mat[l+1, ] = family$linkinv(eta_hat_mat[l+1, ])
    
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
      "hat_mat" = hat_mat,
      #"QMJ" = Q_mj_list,
      "obsvar" = obsvar,
      "sum_to_zero" = stz)
    
    # Combine the two return lists
    return_list = c(return_list, return_list_add)
  }
  
  # Return the fitted model
  return(return_list)
}

GAMBoost_stumps_create_total_design_matrix = function(X) {
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

GAMBoost_stumps_predict = function(model, X_new = NA, after_iter = NA,
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
    return_object = eta_sum + rep(model$intercept, length(eta_sum))
    
    # If we want to return the response mu.
    if (version == 'response'){
      # apply mu = h(eta) = g^(-1)(eta)
      return_object = model$family$linkinv( return_object)
    }
  } 
  return(return_object)
}

GAMBoost_stumps_boot_cofidence_bands = function(models, X_original, after_iter) {
  # models: List of GAMBoost_stumps models fitted to bootstrap samples. Should be 
  #         at least 200 for somewhat decent estimates.
  # X_original: The full observed training data which the boot samples are based on. 
  # after_iter: Specify the iteration of the bootstrap models we want to use. 
  
  # For each model we use the predict function to obtain the predicted eta contr
  # for each of the original observed values for each covariate.
  # Then we have B values for each of the n values for each of the p covariate 
  # functions. We then extract the empirical 95% CI for each of these quantities.
  
  # Assume all bootmodels are fitted to the same number of iterations. 
  # Check for invalid after_iter argument.
  if (after_iter > models[[1]]$num_iter) {
    stop("To high after_iter. Bootstrap models are not fitted for this iteration.")
  }
  
  # Get number of bootstraped models
  B = length(models)
  
  # Get the dimensions
  n = nrow(X_original)
  p = ncol(X_original)
  
  # Find the order of the observed values and sort them.
  # Not needed, but makes thing more logical.
  # X_original_order = apply(X_original, 2, order)
  X_original_sort  = apply(X_original, 2, sort)
  
  # Create a list containing p matrices of dimension n times B.
  # The ith row of the jth matrix will contain the B different 
  # predicted values for the ith value of the jth covariate function.   
  eta_contribution_list = list()
  for (j in 1:p) {
    eta_contribution_list[[j]] = matrix(NA, nrow = n, ncol = B)
  }
  
  # Iterate through all the bootstrap models.
  for (b in 1:B) {
    # Get the current bootstrapped model.
    current_model = models[[b]]
    
    # Predict the eta contribution for each f_j at the n original values.
    eta_contribution_values = GAMBoost_stumps_predict(current_model,
                                                      X_new = X_original_sort,
                                                      after_iter = after_iter,
                                                      version = 'terms')
    
    # Add the predicted contributions to the list of matrices
    for (j in 1:p) {
      eta_contribution_list[[j]][ ,b] = eta_contribution_values[,j] 
    }
  }
  
  # Create a new list to store the p matrices of dimension 2xn,
  # where the jth matrix contain the 95% empirical CI for the n
  # observations at the original observed values fo the jth covariate
  boot_confidence_bands_list = list()
  for (j in 1:p) {
    boot_confidence_bands_list[[j]] = apply(eta_contribution_list[[j]],
                                            1, quantile, c(0.025, 0.975))
  }
  
  #plot(curr_model, select = 2)
  #lines(X_original_sort[,2], boot_se_bands_list_splines[[2]][1,], col = 'orange', lty = 2)
  #lines(X_original_sort[,2], boot_se_bands_list_splines[[2]][2,], col = 'orange', lty = 2)
  
  return(boot_confidence_bands_list)
}

GAMBoost_stumps_coverage_covariate = function(model, boot_models, cn,
                                              after_iter = NA, type = c('smooth', 'step')) {
  # model: A GAMBoost_stumps model
  # boot_models: A list of GAMBoost_stumps models fitted to bootstrap 
  #              samples from the training data used on fitting 'model'.
  # cn: The scalar used in generating the linear predictor eta in the training data.
  # type: If the true lines are the smooth or the step functions.
  # cn: The signal to noise ratio used in creating the training data.
  
  if (is.na(after_iter)) {
    after_iter = model$num_iter
  } else {
    stopifnot(after_iter >= 0 & after_iter <= model$num_iter)
  }
  
  # Find the x_ij values we want to evaluate
  x_values_sort = apply(model$X, 2, sort)
  
  # Find the true values of the true f_j functions.
  true_f_values = matrix(0, nrow = model$n, ncol = model$p)
  if (type == 'smooth') {
    true_f_values[,1] = cn*x_values_sort[,1]
    true_f_values[,3] = cn*(2*x_values_sort[,3]^2) - 2*cn/3 # add term such that integrates to 0
    true_f_values[,5] = cn*sin(5*x_values_sort[,5])
  } 
  if (type == 'step') {
    sf = stepfun(c(0), c(-1, 1), right=TRUE)
    true_f_values[,1] = cn*(0.5*sf(x_values_sort[,1]))
    true_f_values[,3] = cn*(0.25*sf(x_values_sort[,3]))
    true_f_values[,5] = cn*(1*sf(x_values_sort[,5]))
  }
  
  # Predict the eta contribution for each f_j
  eta_contribution_values = GAMBoost_stumps_predict(model, 
                                                    X_new = x_values_sort,
                                                    after_iter = after_iter,
                                                    version = 'terms')
  # Compute the empirical confidence bands
  boot_bands = GAMBoost_stumps_boot_cofidence_bands(boot_models, model$X, after_iter)
  
  
  # End coverage 
  coverage = list('Approximate' = rep(NA, model$p), 'Empirical' = rep(NA, model$p))
  
  par(mfrow = c(1, model$p), mar = c(3.5,3.5,2,1))
  for (s in 1:model$p) {
    ### Calculate the approximate pointwise confidence bands.
    {
      stand_err = sqrt(model$obsvar[[s]][ ,after_iter+1])
      stand_err = stand_err[order(model$X[,s])]
      # Check if the variable has not yet been fitted
      if (is.na(stand_err[1])) {
        plot(0, type = "n", ylim = c(0,1), xlim = c(0,model$n), xlab = "", ylab = "")
        title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
        next
      }
      lower = eta_contribution_values[,s] - 2*stand_err
      upper = eta_contribution_values[,s] + 2*stand_err
      cofidence_band_approx = rbind(lower, upper)
      
      # Find their coverage
      inside_cofidence_band_approx = (cofidence_band_approx[1,] <= true_f_values[,s]) &
        (true_f_values[,s] <= cofidence_band_approx[2,])
      
      inside_cofidence_band_approx_frequency = cumsum(inside_cofidence_band_approx) / (1:model$n)
    }
    
    ### Look at the empirical confidence bands
    {
      # Find their coverage
      inside_cofidence_band_empirical = (boot_bands[[s]][1,] <= true_f_values[,s]) &
        (true_f_values[,s] <= boot_bands[[s]][2,])
      
      inside_cofidence_band_empirical_frequency = cumsum(inside_cofidence_band_empirical) / (1:model$n)
      
    }
    
    # Plot the coverage development.
    plot(x_values_sort[,s], inside_cofidence_band_approx_frequency, type = 'l', ylim = c(0,1), col = nice_colors[1],
         xlab = "", ylab = "")
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    lines(x_values_sort[,s], inside_cofidence_band_empirical_frequency, lty = 1, type = 'l', col = nice_colors[4])
    abline(h = 0.95, lty = 4)
    
    # Add title to the middle most plot
    if (s == floor((model$p+1)/2)) {
      title(main = paste("Iteration:", after_iter))
    }
    
    coverage[['Approximate']][s] = inside_cofidence_band_approx_frequency[model$n]
    coverage[['Empirical']][s]   = inside_cofidence_band_empirical_frequency[model$n]
  }
  
  return(coverage)
}

GAMBoost_stumps_plot_predictor_contribution = function(model,
                                                       type = c('smooth', 'step', 'advancedstep'), 
                                                       boot_models = NULL,
                                                       cn = NULL,
                                                       after_iter = NULL,
                                                       save_filename = NULL,
                                                       xmin = -1,
                                                       xmax = 1,
                                                       ymin = NULL,
                                                       ymax = NULL,
                                                       lwd = 1,
                                                       colors = TRUE) {
  
  # model: A fitted GAMBoost_stumps model
  # boot_models: bootstrap models based on the data in model.
  # type: Either the smooth seting or the stepwise setting.
  # cn: The scalar used in generating the linear predictor eta in the training data.
  # after_iter: For which iteration we want see the functions. 
  #             Must be higher or equal to 1.
  
  # Check for valid type
  type = match.arg(type)
  
  # If we are not provided with bootstrap models we cannot compute 
  # bootstrap confidence bands.
  if (is.null(boot_models)) {
    boot_bands = NULL
  }
  
  # Get the dimesniosn of the dataset
  X_train = model$X
  n = nrow(X_train)
  p = ncol(X_train)
  
  # Get the number of iterations that was perfomed in the boosting fitting
  num_iter = model$num_iter
  
  if (is.null(after_iter)) {
    # If at_step is not chosen we use the last iterations parameters
    after_iter = num_iter
  }
  
  # Check if after_iter are legit values
  stopifnot(sum((after_iter > 0)) == length(after_iter) &&
              sum((after_iter <= num_iter)) == length(after_iter))
  
  # Create the values that are going to be evaluated
  x_values = seq(xmin, xmax, length.out = 500)
  
  # Find the values of the true f_j functions.
  true_y_values = matrix(NA, nrow = length(x_values), ncol = p)
  if (type == 'smooth') {
    true_y_values[,1] = cn*x_values
    true_y_values[,2] = cn*0*x_values
    true_y_values[,3] = cn*(2*x_values^2) - 2*cn/3 # add term such that integrates to 0
    true_y_values[,4] = cn*0*x_values
    true_y_values[,5] = cn*sin(5*x_values)
  } 
  if (type == 'step') {
    sf = stepfun(c(0), c(-1, 1), right=TRUE)
    true_y_values[,1] = cn*(0.5*sf(x_values))
    true_y_values[,2] = cn*0*x_values
    true_y_values[,3] = cn*(0.25*sf(x_values))
    true_y_values[,4] = cn*0*x_values
    true_y_values[,5] = cn*(1*sf(x_values))
  }
  if (type == 'advancedstep') {
    step1 = approxfun(x=seq(-1,1,0.4), y = c(-2,-1,0,1,2,2), method = "constant")
    step2 = approxfun(x=seq(-1, 1, length.out = 7), y = c(1,-1,2,-2,3,-3,-3), method = "constant")
    step3 = approxfun(x=seq(-1,1,0.4), y = c(2,-2,0,-2,2,2), method = "constant")
    step4 = approxfun(x=seq(-1,1, length.out = 7), y = c(1,0,-1, -1, 0, 1,1), method = "constant")
    step5 = approxfun(x=seq(-1,1,2), y = c(0,0), method = "constant")
    
    true_y_values[,1] = cn*step1(x_values)
    true_y_values[,2] = cn*step2(x_values)
    true_y_values[,3] = cn*step3(x_values)
    true_y_values[,4] = cn*step4(x_values)
    true_y_values[,5] = cn*step5(x_values)
  }
  
  if (is.null(ymin) | is.null(ymax)) {
    # Find the range of y values
    range_y = range(true_y_values)
    max_y = max(abs(true_y_values))
    range_length = range_y[2] - range_y[1]
    ymin = -max_y - 0.2*range_length
    ymax = max_y + 0.2*range_length
  }
  
  # 
  #X_new = model$X[apply(model$X, 2, order), ]
  # X_original_order = apply(X_original, 2, order)
  X_new = apply(model$X, 2, sort)
  
  # matrix(rep(x_values, each = p), ncol=p, byrow = T)
  
  # Check if we are supposed to save it or just display it
  if (!is.null(save_filename)){
    png(paste(save_filename, ".png", sep = ""), width = 3500, height = 1500*length(after_iter)^(2/3), res = 350)
  }
  
  #png(paste("temp10.png", sep = ""), width = 3500, height = 1500*length(after_iter)^(2/3), res = 350)
  # Set some plot parameters
  par(mfrow = c(length(after_iter), p), mar = c(3.5,3.5,2,1))
  
  # Iterate over the different boosting iteration numbers 
  for (ai in after_iter) {
    
    # Predict the eta contribution for each f_j
    eta_contribution_values = GAMBoost_stumps_predict(model, 
                                                      X_new = X_new,
                                                      after_iter = ai,
                                                      version = 'terms')
    
    # Create the empirical pointwise bootstrap confidence bands.
    if (!is.null(boot_models)) {
      boot_bands = GAMBoost_stumps_boot_cofidence_bands(boot_models, X_train, ai)
    } 
    
    
    # Plot each of the covariate specific functions.
    for (s in 1:p) {
      # Calculate the approximate pointwise confidence bands.
      #ci_list = GAMBoost_calc_confidence_bands(model, s, ai, 1)
      stand_err = sqrt(model$obsvar[[s]][ ,ai+1])
      stand_err = stand_err[order(model$X[,s])]
      lower = eta_contribution_values[,s] - 2*stand_err
      upper = eta_contribution_values[,s] + 2*stand_err
      
      # Crete a canvas. NEED TO FIX YLIM SUCH THAT IT IS MORE GENERIC!!!!
      plot(NA, type = 'n', ylim = c(ymin, ymax), xlim = c(xmin, xmax), xlab = '', ylab = '')
      title( xlab = bquote("x"[.(s)]), ylab = expression(paste(eta, " contrib.")), line = 2.25)
      
      # Plot the true curve
      lines(x_values, true_y_values[,s], lty = 1, col='darkgray', lwd = lwd+0.5, type = 's')
      
      
      # Plot the predicted contribution curves
      lines(X_new[,s], eta_contribution_values[,s], lwd = lwd, type = 's')
      
      # Approximate confidence bands
      if (colors) {
        lines(X_new[,s], lower, lty = 2, col = nice_colors[1], lwd = lwd, type = 's')
        lines(X_new[,s], upper, lty = 2, col = nice_colors[1], lwd = lwd, type = 's')
      } else {
        lines(X_new[,s], lower, lty = 3, col = "black", lwd = lwd, type = 's')
        lines(X_new[,s], upper, lty = 3, col = "black", lwd = lwd, type = 's')
      }
      
      # Add points to see where the training data is located
      points(jitter(model$X[,s]), rep(ymin, length(model$X[,s])), cex =.5, pch ="|", col = "darkgrey")
      
      # Add the empirical pointwise bootstrap confidence bands if provided
      if (!is.null(boot_models)) {
        if (colors) {
          lines(sort(model$X[,s]), boot_bands[[s]][1,], lty = 2, col = nice_colors[4], lwd = lwd, type = 's')
          lines(sort(model$X[,s]), boot_bands[[s]][2,], lty = 2, col = nice_colors[4], lwd = lwd, type = 's')
        } else {
          lines(sort(model$X[,s]), boot_bands[[s]][1,], lty = 2, col = "black", lwd = lwd, type = 's')
          lines(sort(model$X[,s]), boot_bands[[s]][2,], lty = 2, col = "black", lwd = lwd, type = 's')
        }
      }  
      
      # Add title to the middle most plot
      if (s == floor((p+1)/2)) {
        title(main = paste("Iteration:", ai))
      }
    }
  }
  
  if (!is.null(save_filename)){
    dev.off()
  }
}


GAMBoost_stumps_inside_conf_bands = function(model,
                                             cn,
                                             training_data,
                                             boot_models = NA,
                                             after_iter = NA,
                                             type = c('smooth', 'step', 'advancedstep'),
                                             colors = TRUE,
                                             ymin = NULL,
                                             ymax = NULL,
                                             figure_type = 1){
  # model: A GAMBoost model fitted to the best iteration
  # boot_models: A list of GAMBoost models fitted to bootstrap 
  #              samples from the training data used on fitting 'model'.
  # cn: The scalar used in generating the linear predictor eta in the training data.
  # training data: The full generated training data list which the model is fitted to
  # After_iter: the optimal iteration
  # type: If the true lines are the smooth or the step functions.
  
  if (is.na(after_iter)) {
    after_iter = model$num_iter
  } else {
    stopifnot(after_iter >= 0 & after_iter <= model$num_iter)
  }
  
  # Get the dimensions
  n = model$n
  p = model$p
  
  # Find the x_ij values we want to evaluate
  x_values_sort = apply(model$X, 2, sort)
  
  # Find the true values of the true f_j functions.
  true_f_values = matrix(0, nrow = n, ncol = p)
  if (type == 'smooth') {
    true_f_values[,1] = cn*x_values_sort[,1]
    true_f_values[,3] = cn*(2*x_values_sort[,3]^2) - 2*cn/3 # add term such that integrates to 0
    true_f_values[,5] = cn*sin(5*x_values_sort[,5])
  } else if (type == 'step') {
    sf = stepfun(c(0), c(-1, 1), right=TRUE)
    true_f_values[,1] = cn*(0.5*sf(x_values_sort[,1]))
    true_f_values[,3] = cn*(0.25*sf(x_values_sort[,3]))
    true_f_values[,5] = cn*(1*sf(x_values_sort[,5]))
  } else if (type == 'advancedstep') {
    step1 = approxfun(x=seq(-1,1,0.4), y = c(-2,-1,0,1,2,2), method = "constant")
    step2 = approxfun(x=seq(-1, 1, length.out = 7), y = c(1,-1,2,-2,3,-3,-3), method = "constant")
    step3 = approxfun(x=seq(-1,1,0.4), y = c(2,-2,0,-2,2,2), method = "constant")
    step4 = approxfun(x=seq(-1,1, length.out = 7), y = c(1,0,-1, -1, 0, 1,1), method = "constant")
    step5 = approxfun(x=seq(-1,1,2), y = c(0,0), method = "constant")
    
    true_f_values[,1] = cn*step1(x_values_sort[,1])
    true_f_values[,2] = cn*step2(x_values_sort[,2])
    true_f_values[,3] = cn*step3(x_values_sort[,3])
    true_f_values[,4] = cn*step4(x_values_sort[,4])
    true_f_values[,5] = cn*step5(x_values_sort[,5])
  }
  
  # Predict the eta contribution for each f_j
  eta_contribution_values = GAMBoost_stumps_predict(model, X_new = x_values_sort,
                                                    after_iter = after_iter,
                                                    version = 'terms')
  
  
  # Compute the empirical confidence bands
  if (is.list(boot_models)) {
    boot_bands = GAMBoost_stumps_boot_cofidence_bands(boot_models, model$X, after_iter)
    
    # End coverage 
    coverage = list('Approximate'    = matrix(NA, ncol = p, nrow = n),
                    'Empirical'      = matrix(NA, ncol = p, nrow = n),
                    'Approximate_mu' = rep(NA, n),
                    'Empirical_mu'   = rep(NA, n))
  } else {
    coverage = list('Approximate'    = matrix(NA, ncol = p, nrow = n),
                    'Approximate_mu' = rep(NA, n))
  }
  
  
  for (s in 1:model$p) {
    ### Calculate the approximate pointwise confidence bands.
    {
      # Calculate the approximate pointwise confidence bands.
      stand_err = sqrt(model$obsvar[[s]][ ,after_iter+1])
      stand_err = stand_err[order(model$X[,s])]
      lower = eta_contribution_values[,s] - 2*stand_err
      upper = eta_contribution_values[,s] + 2*stand_err
      cofidence_band_approx = rbind(lower, upper)
      
      # Find their coverage
      inside_cofidence_band_approx = (cofidence_band_approx[1,] <= true_f_values[,s]) &
        (true_f_values[,s] <= cofidence_band_approx[2,])
      
      coverage[['Approximate']][,s] = inside_cofidence_band_approx
    }
    
    ### Look at the empirical confidence bands
    if (is.list(boot_models)) {
      # Find their coverage
      inside_cofidence_band_empirical = (boot_bands[[s]][1,] <= true_f_values[,s]) &
        (true_f_values[,s] <= boot_bands[[s]][2,])
      
      coverage[['Empirical']][,s] = inside_cofidence_band_empirical
    }
  }
  
  # Predict the mu
  predicted_mu = model$family$linkinv(model$eta_hat[after_iter+1, ])
  stumps_pointwise_var = sqrt(diag(model$hat_mat_last %*%
                                     diag(model$family$variance(predicted_mu)) %*%
                                     t(model$hat_mat_last)))
  
  stumps_lower = predicted_mu - 2*stumps_pointwise_var
  stumps_upper = predicted_mu + 2*stumps_pointwise_var
  coverage[['Approximate_mu']] = stumps_lower <= training_data$mu & training_data$mu <= stumps_upper
  
  # Compute the bootsrap predictions
  if (is.list(boot_models)) {
    boot_predictions = matrix(NA, ncol = n, nrow = length(boot_models))
    {
      for (m in 1:length(boot_models)) {
        curr_mod = boot_models[[m]]
        boot_predictions[m,] =  GAMBoost_stumps_predict(curr_mod,
                                                        X_new = model$X,
                                                        after_iter = after_iter,
                                                        version = 'respons')
      }
    }
    boot_ci = apply(boot_predictions, 2, quantile, c(0.025, 0.975))
    coverage[['Empirical_mu']] = boot_ci[1,] <= training_data$mu & training_data$mu <= boot_ci[2,]
  }
  
  
  ### Plot 
  # Create the title
  if (training_data$type == 'smooth') {
    if (training_data$family$family == 'gaussian') {
      main_tit = expression(paste("Gaussian Smooth Stumps: Empirical and Approximate CI for ", mu[i]))
    } else if (training_data$family$family == 'binomial') {
      main_tit = expression(paste("Binomial Smooth Stumps: Empirical and Approximate CI for ", mu[i]))      
    } else {
      main_tit = expression(paste("Poisson Smooth Stumps: Empirical and Approximate CI for ", mu[i]))     
    }
  } else if (training_data$type == 'step') {
    if (training_data$family$family == 'gaussian') {
      main_tit = expression(paste("Gaussian Stepwise Stumps: Empirical and Approximate CI for ", mu[i]))
    } else if (training_data$family$family == 'binomial') {
      main_tit = expression(paste("Binomial Stepwise Stumps: Empirical and Approximate CI for ", mu[i]))      
    } else {
      main_tit = expression(paste("Poisson Stepwise Stumps: Empirical and Approximate CI for ", mu[i]))     
    }
  } else if (training_data$type == 'advancedstep') {
    if (training_data$family$family == 'gaussian') {
      main_tit = expression(paste("Gaussian Advanced Stepwise Stumps: Empirical and Approximate CI for ", mu[i]))
    } else if (training_data$family$family == 'binomial') {
      main_tit = expression(paste("Binomial Advanced Stepwise Stumps: Empirical and Approximate CI for ", mu[i]))      
    } else {
      main_tit = expression(paste("Poisson Advanced Stepwise Stumps: Empirical and Approximate CI for ", mu[i]))     
    }
  }
  
  if (is.null(ymin) | is.null(ymax)) {
    ymin = min(stumps_lower)
    ymax = max(stumps_upper)
  }
  stumpsorder = order(training_data$mu)
  if (figure_type == 1) {
    #par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
    plot(training_data$mu[stumpsorder], type = 'l', col = 'darkgray', lwd = 1.5,
         ylim = c(ymin, ymax), xlab = '', ylab = '',
         main = main_tit, cex.main = 1.13)
    points(training_data$mu[stumpsorder], col = 'darkgray', pch = 19, cex = 0.5)
    title(xlab = bquote("Observations in Increasing Order W.R.T. True"~mu[i]),
          ylab = bquote(mu[i]), line = 2)
    lines(predicted_mu[stumpsorder], lty = 1, col = 'black')
    points(predicted_mu[stumpsorder], col = 'black', pch = 19, cex = 0.4)
    
    if (colors) {
      lines(stumps_upper[stumpsorder], lty = 2, col = nice_colors[1])
      lines(stumps_lower[stumpsorder], lty = 2, col = nice_colors[1])    
    } else {
      lines(stumps_upper[stumpsorder], lty = 3, col = 'black')
      lines(stumps_lower[stumpsorder], lty = 3, col = 'black')
    }
    
    if (is.list(boot_models)) {
      if (colors) {
        lines(boot_ci[1,stumpsorder], lty = 2, col = nice_colors[4])
        lines(boot_ci[2,stumpsorder], lty = 2, col = nice_colors[4])
        legend("topleft", lty = c(1,1,2,2), 
               col = c('darkgray', 'black', nice_colors[4], nice_colors[1]),
               legend = c(expression(paste("True ", mu[i])),
                          expression(paste("Predicted ", mu[i])),
                          'Empirical', 'Approximate'), cex = 0.8)      
      } else {
        lines(boot_ci[1,stumpsorder], lty = 2, col = "#969696")
        lines(boot_ci[2,stumpsorder], lty = 2, col = "#969696")
        legend("topleft", lty = c(1,1,2,3),
               col = c('darkgray', 'black', "#969696", "black"),
               legend = c(expression(paste("True ", mu[i])),
                          expression(paste("Predicted ", mu[i])),
                          'Empirical', 'Approximate'), cex = 0.8)
      }
      
    } else {
      if (colors) {
        legend("topleft", lty = c(1,1,2), col = c('darkgray', 'black', nice_colors[1]),
               legend = c(expression(paste("True ", mu[i])),
                          expression(paste("Predicted ", mu[i])),
                          'Approximate'),
               cex = 0.8)
      } else {
        legend("topleft", lty = c(1,1,3), col = c('darkgray', 'black', 'black'),
               legend = c(expression(paste("True ", mu[i])),
                          expression(paste("Predicted ", mu[i])),
                          'Approximate'),
               cex = 0.8) 
      }
    }
  } else {
    plot(1:n, training_data$mu[stumpsorder], col = 'darkgray',
         ylim=c(ymin, ymax), pch=19, xlab="", ylab="",
         main = main_tit, cex.main = 1.13)
    title(xlab = bquote("Observations in Increasing Order W.R.T. True"~mu[i]),
          ylab = bquote(mu[i]), line = 2)
    points(1:n, predicted_mu[stumpsorder], col = 'black', pch = 17)
    # hack: we draw arrows but with very special "arrowheads"
    arrows(1:n - 0.0, stumps_lower[stumpsorder],
           1:n - 0.0, stumps_upper[stumpsorder],
           length=0.02, angle=90, code=3, col = 'black', lwd = 2.5)
    arrows(1:n + 0.0, boot_ci[1,stumpsorder],
           1:n + 0.0, boot_ci[2,stumpsorder],
           length=0.02, angle=90, code=3, col = "#969696", lty = 1, lwd = 1)
    legend("topleft", lty = c(1,1,2,1),
           pch = NULL,
           col = c('darkgray', 'black', "#969696", 'black'),
           legend = c(expression(paste("True ", mu[i])),
                      expression(paste("Predicted ", mu[i])),
                      'Empirical',
                      'Approximate'), cex = 0.75)
  }
  
  return(coverage)
}



