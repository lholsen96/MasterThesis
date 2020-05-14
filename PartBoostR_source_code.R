# File contains three funcitons
# PartBoostR: Fit a PartBoostR model
# PartBoostR_predict: takes a PartBoostR model and predicts the outcomes
# plot_coeff_build_up: plot the coff buildup of a PartBoostR model

PartBoostR = function (X, y, num_iter, lambda, V_compulsory = c(), 
                       return_hatmat = FALSE, print = FALSE) {
  # Function that fits a linear model, without intercept, based on the 
  # PartBoostR algorithm. Does not supporst blocks, only componentwise.
  # X: the n times p design matrix, assumed normalized.
  # y: the n-dimensional response vector.
  # num_iter: number of boosting iterations.
  # lambda: penalization factor used in ridge regression.
  # V_compusory: Set of indices of mandatory variables in each boosting iteration.
  
  # Find the number of observations
  n = nrow(X)
  p = ncol(X)
  
  # Check for invalid compulsory variables 
  if(length(V_compulsory) != 0 && (range(V_compulsory)[1] < 1 ||  range(V_compulsory)[2] > p)) {
    stop("Invalid compulsory variable indices.")
  }
  
  # Ensure no duplicates in V_compulsory
  V_compulsory = sort(unique(c(V_compulsory)))
  
  # Check if we are in the extreme case of BoostR
  BoostR_algo = length(V_compulsory) == p
  
  # Calculate the number of parameters to be updated in each iteration
  num_param_update = min(length(V_compulsory) + 1, p)
  
  # List over indices not in the compulsory list
  if (length(V_compulsory) == 0) { # There are no compulsory variables
    V_not_compulsory = 1:p                  
  } else if (BoostR_algo) { # All the variables are compulsory -> BoostR
    V_not_compulsory = c()                  
  } else { # There are some compulsory variables
    V_not_compulsory = (1:p)[-V_compulsory] 
  }
  
  # Create the possible subsets
  subsets = list()
  if (BoostR_algo) { # The only subset is the whole set
    subsets[[1]] = V_compulsory
  } else {
    for (j in 1:length(V_not_compulsory)) { # Combine compulsory and potential
      subsets[[j]] = c(V_compulsory, V_not_compulsory[j])
    }
  }
  
  # Variables to store the acumulated values of quantities of interest
  beta_hat = matrix(rep(NA, (num_iter+1)*p), nrow = num_iter+1)
  mu_hat = matrix(rep(NA, (num_iter+1)*n), nrow = num_iter+1)
  hat_mat = list()
  df_hat = rep(NA, num_iter+1)
  mse_training = rep(NA, num_iter+1)
  update_trace = matrix(rep(NA, num_iter*num_param_update), nrow = num_iter)
  cumulative_product = diag(n) # Will store iterative products for the hat matrix
  
  ### Step 1: Initialization
  B0 = solve(t(X) %*% X + lambda * diag(p)) %*% t(X)
  beta_hat[1, ] = B0 %*% y
  mu_hat[1, ] = X %*% beta_hat[1, ]
  mse_training[1] = mean((mu_hat[1, ] - y)^2) 
  
  # Compute the hat matrix
  S_old = X %*% B0
  hat_mat[[1]] = S_old
  df_hat[1] = sum(diag(hat_mat[[1]]))
  
  ### Step 2: Iteration
  if (num_iter >= 1) {
    for (m in 1:num_iter) {
      # Array to store all the L2 losses to see which update is optimal
      L2_loss_array = c()
      
      # iterate over all subsets we need to evaluate
      for (j in 1:length(subsets)) {
        # Create the partial design matrix
        X_Vm_j = as.matrix(X[ ,subsets[[j]]])
        
        # Calculate the ridge solution of the resiudal model
        beta_Vm_j = solve(t(X_Vm_j) %*% X_Vm_j + lambda*diag(num_param_update)) %*% t(X_Vm_j) %*% (y-mu_hat[m, ])
        
        # The new estimates
        mu_Vm_j = mu_hat[m, ] + X_Vm_j %*% beta_Vm_j
        
        # Calculate the L2-loss
        L2_loss_temp = sum((mu_Vm_j - y)^2)
        L2_loss_array = c(L2_loss_array, L2_loss_temp)
        
      }
      # Need to find the update that yields the best fit. min L2 loss
      min_idx = which.min(L2_loss_array)
      
      # Partial matrix for best update
      update_idx = subsets[[min_idx]]  # Find the variables with best update
      update_trace[m, ] = update_idx   # Record which variables that got updated
      X_m = as.matrix(X[, update_idx]) # Get the best partial design matrix
      
      # Obtain the parameter update
      B_m = solve(t(X_m) %*% X_m + lambda*diag(num_param_update)) %*% t(X_m) 
      beta_hat[m+1, ] = beta_hat[m, ]
      beta_hat[m+1, update_idx] = beta_hat[m+1, update_idx] + B_m %*% (y-mu_hat[m, ])
      
      # The new estimates
      mu_hat[m+1, ] = X %*% beta_hat[m+1, ] 
      
      # Calculate the training MSE
      mse_training[m+1] = mean((mu_hat[m+1, ] - y)^2) 
      
      ### Calculate the new hat matrix
      # update the cumulative product
      cumulative_product = (diag(n) - S_old) %*% cumulative_product
      S_new = X_m %*% B_m
      
      # calcualte the new hat matrix and effective degrees of freedom
      hat_mat[[m+1]] = hat_mat[[m]] + S_new %*% cumulative_product
      df_hat[m+1] = sum(diag(hat_mat[[m+1]]))
      
      # Verify the results # THIS CAN BE REMOVED
      if (!all.equal(c(hat_mat[[m+1]] %*% y), mu_hat[m+1, ])) {
        warning("Results obtained from hat matrix differ from the
                other results.")
      }
      
      # Set the old S matrix to be the new one
      S_old = S_new 
      
      # We print every 10% percent iteration
      if (print && !(m %% (num_iter %/% 10))) {
        cat(sprintf("Iteration: %3d", m))
        cat(sprintf("  MSE: "), mse_training[m+1])
        cat(sprintf("  Beta: "), c(beta_hat[m+1, ] ), sprintf("\n"))
      }
    }
  }
  
  # Creteat list to return (Maybe use structure?`)
  return_list = list(mu_hat = mu_hat, beta_hat = beta_hat, num_iter = num_iter, 
                     lambda = lambda, V_compulsory = V_compulsory,
                     df_hat = df_hat, update_trace = update_trace, 
                     mse_training = mse_training, n = n, p = p, X = X, y = y)
  
  if (return_hatmat) {
    return_list[["hat_mat"]] = hat_mat
  }
  
  return(return_list)
}



PartBoostR_predict = function(model, X_test, y_test = NULL, after_iter = NULL) {
  # Function that takes a PartBoostModel and predicts the outcome for the last iter
  # or at the specified iteration after_iter. If it is 
  # given a y_test set is also computes the mean squared error, either for iteration
  # after_iter or for all iterations if no iteration is provided.
  
  # model: A PartBoostR model
  # X_test: The test data. Assumed normalized with respect to the normalization
  #         of the training data.
  # y_test: The true y values for X_test
  # after_iter: which version of the model we want to use to predict
  
  if (is.null(after_iter)) {
    mse_loss = rep(NA, model$num_iter+1) # Variable to store the mse losses
    
    # Iterate over all the models
    for (m in 1:(model$num_iter+1)) {
      # We are only interested in a specific iteration
      beta_values = model$beta_hat[m, ]
      pred_values = X_test %*% beta_values
      
      # Check if we can calculate the mse
      if (!is.null(y_test)) {
        mse_loss[m] = mean((y_test - pred_values)^2)
      } 
    }
    # We return the predictions from the m'th (last) iteration
    after_iter = m
    
  } else {
    # We are only interested in a specific iteration
    beta_values = model$beta_hat[after_iter+1, ]
    pred_values = X_test %*% beta_values
    
    if (!is.null(y_test)) {
      mse_loss = mean((y_test - pred_values)^2)
    } 
  }
  return(list(pred_values = pred_values, after_iter = after_iter, mse_loss = mse_loss))
}



# Function which plots the coefficient build ups MOVE TO OTHER FILE
plot_coeff_build_up = function(models, ls_solutions = NULL, title = NULL, max_iter = NULL,
                               ylim_val = NULL, save_plot_name = NULL) {
  # Function to plot the build ups
  # Quite general, but makes assumptions that if models containt two
  # models, then the first is BoostR and the second PartBoostR with
  # same penalization.
  # If three or more models, then it assume that it is PartBoostR
  # with different penalities
  
  # models: list of PartBoostR with same number of parameters p.
  # ls_solutions: array of length p with the least squares solution
  # title: name of the plot
  # max_iter: if we only want to plot up to max_iter iteration
  # ylim_val: limits in the plot of ht ecoeff build up
  # save_plot_name: Save the figure to file with this name if provided
  
  num_models = length(models) 
  num_param = models[[1]]$p
  p = models[[1]]$p
  
  # Check that we got correct number of parameters if given least squares solutions.
  if (!is.null(ls_solutions) && num_param != length(ls_solutions)) {
    stop("Inconsistency between the number of parameters in models and LS solutions.")
  }
  
  # If we have several models we want to check that we have the 
  # same amount of parameters
  if (num_models > 1) {
    # Assert that the models have the same number of parameters
    if (length(unique(sapply(models, function(x) x$p) )) != 1) {
      stop("The models have different number of parameters.")
    } 
    # 
    num_iters = sapply(models, function(x) x$num_iter)
    if (length(unique(num_iters)) != 1) {
      warning("The models are run for different number of iterations.")
    }
  }
  
  # Find the penalization factors
  lambdas = sapply(models, function(x) x$lambda)
  if (length(unique(num_iters)) == 1) {
    lambdas == lambdas[1]
  }
  
  # Use the largest number of iterations
  num_iter = max(num_iters)
  
  # Find the ylim for our plot
  ranges = sapply(models, function(x) apply(x$beta_hat, 2, range))
  min_beta = min(ranges)
  max_beta = max(ranges)
  range_beta = max_beta - min_beta 
  min_beta = min_beta - 0.05 * range_beta
  max_beta = max_beta + 0.35 * range_beta
  
  ### Plot the first model by itself and then loop over the remaining
  # Check if we are going to save the plot
  if (!is.null(save_plot_name)) {
    png(paste(save_plot_name, ".png", sep = ""), width = 2500, height = 1500, res = 300)
  }
  
  par(mar = c(3.5,3.5,2,1))
  
  # Allow for some user inputs
  num_iter_plot = ifelse(!is.null(max_iter), max_iter, num_iter) 
  if (!is.null(ylim_val)) {
    ylim_val = ylim_val
  } else {
    ylim_val =  c(min_beta, max_beta)
  }
  
  # Create the plot
  plot(NA, xlim = c(0, num_iter_plot), ylim = ylim_val,
       ylab = "", xlab = "", main = "", cex.main = 1)
  title(ylab=bquote("Standardized" ~ beta['j']), xlab = "Iterations", line = 2)
  
  colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
             "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
             "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)
  colors = colors[1:p]
  
  if (is.null(title)) {
    if (num_models == 1) {
      #colors = brewer.pal(2*p, "Paired")
      #colors = c(colors[seq.int(2L, length(colors), 2L)])
      title(main = bquote("PartBoostR:" ~ lambda ~ '=' ~ .(lambda)))
    }
    if (num_models == 2) {
      #colors = brewer.pal(2*p, "Paired")
      #colors = c(colors[seq.int(2L, length(colors), 2L)], colors[seq.int(1L, length(colors), 2L)])
      title(main = "BoostR (Solid) and PartBoostR (Dashed)")
    }
    
    if (num_models > 2) {
      #colors = brewer.pal(2*p, "Paired")
      #colors = c(colors[seq.int(2L, length(colors), 2L)])
      title(main = "PartBoostR With Different Regularization")
    }
  } else {
    title(main = title)
  }

  
  legend_names = c()  # Array to store the legend names
  legend_names2 = c()
  for (j in 1:p) {
    legend_names = c(legend_names, as.expression(bquote(beta[.(j)] ~ " ")))
  }  
  for (m in 1:num_models) {
    legend_names2 = c(legend_names2, as.expression(bquote(lambda[.(m)] ~ "=" ~ .(lambdas[m]) ~ " ")))
  }
  
  # Add the first legend for the parameters
  legend("topright",lty = rep(1,p), col = colors[1:p],
         legend = legend_names, cex = 0.75)  
  
  legend("topleft",lty = 1:num_models, col = rep('black', num_models),
         legend = legend_names2, cex = 0.75)
  
  # Add least squares solutions if provided
  if (!is.null(ls_solutions)) {
    points(rep(num_iter, p), ls_solutions)
  }
  
  # Plot the curves
  # Get the betas for the first model
  curr_beta = models[[1]]$beta_hat
  
  # Pad with NA if fitted with different number of iterations
  curr_beta = rbind(curr_beta, matrix(rep(NA,  p * (num_iter + 1 - nrow(curr_beta))), ncol =  p))
  matlines(curr_beta, lty = 1, col = colors[1:p], type = 's')
  
  # Repeat the colors for the number of models
  colors = rep(colors, num_models)
  
  for (m in 2:num_models) {
    # Get the betas for the mth model
    curr_beta = models[[m]]$beta_hat
    
    # Pad with NA if fitted with different number of iterations
    curr_beta = rbind(curr_beta, matrix(rep(NA,  p * (num_iter + 1 - nrow(curr_beta))), ncol =  p))
    
    # Plot the coeff build up
    matlines(curr_beta, lty = m, col = colors[((m-1)*p +1):(m*p)], type = 's')
  }
  if (!is.null(save_plot_name)) {
    dev.off()
  }
}
