##### File containing the aid functions for GAMBoost_splines

library(GAMBoost)

GAMBoost_splines_boot_cofidence_bands = function(models, X_original, after_iter){
  
  # models: List of GAMBoost models fitted to bootstrap samples. Should be 
  #         at least 200 for somewhat decent estimates.
  # X_original: The full observed training data 
  # after_iter: Specify the iteration of the bootstrap models we want to use. 
  
  # For each model we use the predict function to obtain the predicted eta contr
  # for each of the original observed values for each covariate.
  # Then we have B values for each of the n values for each of the p covariate 
  # functions. We then extract the empirical 95% CI for each of these quantities.
  
  # Assume all bootmodels are fitted to the same number of iterations. 
  # Check for invalid after_iter argument.
  if (after_iter > models[[1]]$stepno) {
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
    eta_contribution_values = predict(current_model, 
                                      newdata = X_original_sort,
                                      type="terms",
                                      at.step = after_iter)
    
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

GAMBoost_splines_coverage_covariate = function(model, boot_models, cn,
                                               after_iter = NA, type = c('smooth', 'step')){
  # model: A GAMBoost model
  # boot_models: A list of GAMBoost models fitted to bootstrap 
  #              samples from the training data used on fitting 'model'.
  # cn: The scalar used in generating the linear predictor eta in the training data.
  # type: If the true lines are the smooth or the step functions.
  
  if (is.na(after_iter)) {
    after_iter = model$stepno
  } else {
    stopifnot(after_iter >= 0 & after_iter <= model$stepno)
  }
  
  # Get the dimensions
  n = model$n
  p = ncol(model$x)
  
  # Find the x_ij values we want to evaluate
  x_values_sort = apply(model$x, 2, sort)
  
  # Find the true values of the true f_j functions.
  true_f_values = matrix(0, nrow = n, ncol = p)
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
  eta_contribution_values = predict(model, newdata = x_values_sort, 
                                    at.step = after_iter, type = 'terms')
  
  # Compute the empirical confidence bands
  boot_bands = GAMBoost_splines_boot_cofidence_bands(boot_models, model$x, after_iter)
  
  # End coverage 
  coverage = list('Approximate' = rep(NA, p), 'Empirical' = rep(NA, p))
  
  par(mfrow = c(1, p), mar = c(3.5,3.5,2,1))
  for (s in 1:p) {
    ### Calculate the approximate pointwise confidence bands.
    {
      # Calculate the approximate pointwise confidence bands.
      ci_list = GAMBoost_calc_confidence_bands(model, s, after_iter, 1)
      
      # Check if the variable has not yet been fitted
      if (!is.list(ci_list)) {
        #if (is.na(ci_list)) {
        plot(0, type = "n", ylim = c(0,1), xlim = c(0,model$n), xlab = "", ylab = "")
        title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
        next
        #}
      }
      
      
      # Get the approximate confidence bands
      cofidence_band_approx = rbind(ci_list$lower, ci_list$upper)
      
      # Find their coverage
      inside_cofidence_band_approx = (cofidence_band_approx[1,] <= true_f_values[,s]) &
        (true_f_values[,s] <= cofidence_band_approx[2,])
      
      inside_cofidence_band_approx_frequency = cumsum(inside_cofidence_band_approx) / (1:n)
    }
    
    ### Look at the empirical confidence bands
    {
      # Find their coverage
      inside_cofidence_band_empirical = (boot_bands[[s]][1,] <= true_f_values[,s]) &
        (true_f_values[,s] <= boot_bands[[s]][2,])
      
      inside_cofidence_band_empirical_frequency = cumsum(inside_cofidence_band_empirical) / (1:model$n)
      
    }
    
    # Plot the coverage development.
    plot(ci_list$x, inside_cofidence_band_approx_frequency, type = 'l', ylim = c(0,1), col = nice_colors[1],
         xlab = "", ylab = "")
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    lines(ci_list$x, inside_cofidence_band_empirical_frequency, lty = 1, type = 'l', col = nice_colors[4])
    abline(h = 0.95, lty = 4)
    
    # Add title to the middle most plot
    if (s == floor((p+1)/2)) {
      title(main = paste("Iteration:", after_iter))
    }
    
    coverage[['Approximate']][s] = inside_cofidence_band_approx_frequency[model$n]
    coverage[['Empirical']][s]   = inside_cofidence_band_empirical_frequency[model$n]
  }
  
  return(coverage)
}

GAMBoost_calc_confidence_bands = function(object, covno, at.step, phi) 
{
  at.step <- at.step + 1
  result <- NULL
  if (length(object$obsvar[[covno + 1]]) > 1) {
    updated <- (1:at.step)[!is.na(object$obsvar[[covno + 
                                                   1]][1, 1:at.step])]
    if (length(updated) > 0) {
      actual.obsvar <- object$obsvar[[covno + 1]][, max(updated)]
      ori.x <- object$x[, covno]
      ori.eta <- predict(object, type = "terms", at.step = at.step - 
                           1)[, covno]
      unique.x <- sort(unique(ori.x))
      mean.eta <- rep(0, length(unique.x))
      upper.eta <- rep(0, length(unique.x))
      lower.eta <- rep(0, length(unique.x))
      for (j in 1:length(unique.x)) {
        this.sd <- sqrt(mean(actual.obsvar[ori.x == 
                                             unique.x[j]], na.rm = TRUE))
        mean.eta[j] <- mean(ori.eta[ori.x == unique.x[j]])
        upper.eta[j] <- mean.eta[j] + 2 * this.sd * 
          sqrt(phi)
        lower.eta[j] <- mean.eta[j] - 2 * this.sd * 
          sqrt(phi)
      }
      result <- list(x = unique.x, eta = mean.eta, upper = upper.eta, 
                     lower = lower.eta)
    }
  }
  return(result)
}

GAMBoost_plot_predictor_contribution = function(model,
                                                cn,
                                                type = c('smooth', 'step'), 
                                                boot_models = NULL,
                                                after_iter = NULL,
                                                save_filename = NULL,
                                                xmin = -1,
                                                xmax = 1,
                                                ymin = NULL,
                                                ymax = NULL,
                                                lwd = 1,
                                                color = TRUE) {
  
  # model: A fitted GAMBoost model
  # cn: the constant used in creating the training data.
  # type: Either the smooth seting or the stepwise setting.
  # boot_models: bootstrap models based on the data in model.
  # after_iter: For which iteration we want see the functions. 
  #             Must be higher or equal to 1.
  
  
  # Check for valid type
  type = match.arg(type)
  
  # If we are not provided with bootstrap models we cannot compute 
  # bootstrap confidence bands.
  if (is.null(boot_models)) {
    boot_bands = NULL
  }
  
  
  # Get the dimension of the training data.
  X_train = model$x
  n = nrow(X_train)
  p = ncol(X_train)
  
  # # Verify the size of bootstrap empirical confidence bands if provided
  # if(!is.null(boot_bands)) {
  #   stopifnot(length(boot_bands) == p)
  #   stopifnot(ncol(boot_bands[[1]]) == n)
  #   stopifnot(nrow(boot_bands[[1]]) == 2)
  # }
  
  # Get the number of iterations that was perfomed in the boosting fitting
  num_iter = model$stepno
  
  if (is.null(after_iter)) {
    # If at_step is not chosen we use the last iterations parameters
    after_iter = num_iter
  }
  
  # Check if after_iter are legit values
  stopifnot(sum((after_iter > 0)) == length(after_iter) &&
              sum((after_iter <= num_iter)) == length(after_iter))
  
  # Create the values that are going to be evaluated
  x_values = seq(xmin, xmax, length.out = 1000)
  X_new = matrix(rep(x_values, each = p), ncol=p, byrow = T)
  
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
  
  if (is.null(ymin) | is.null(ymax)) {
    # Find the range of y values
    range_y = range(true_y_values)
    max_y = max(abs(true_y_values))
    range_length = range_y[2] - range_y[1]
    ymin = -max_y - 0.2*range_length
    ymax = max_y + 0.2*range_length
  }
  
  
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
    eta_contribution_values = predict(model, 
                                      newdata = X_new,
                                      type="terms",
                                      at.step = ai)
    
    # Create the empirical pointwise bootstrap confidence bands.
    if (!is.null(boot_models)) {
      boot_bands = GAMBoost_splines_boot_cofidence_bands(boot_models, X_train, ai)
    } 
    
    # Plot each of the covariate specific functions.
    for (s in 1:p) {
      # Calculate the approximate pointwise confidence bands.
      ci_list = GAMBoost_calc_confidence_bands(model, s, ai, 1)
      
      # Crete a canvas. NEED TO FIX YLIM SUCH THAT IT IS MORE GENERIC!!!!
      plot(NA, type = 'n', ylim = c(ymin, ymax), xlim = c(xmin, xmax), xlab = '', ylab = '')
      title( xlab = bquote("x"[.(s)]), ylab = expression(paste(eta, " contrib.")), line = 2.25)
      
      # Plot the true curve
      lines(x_values, true_y_values[,s], lty = 1, col='darkgray', lwd = lwd+0.5)
      
      # Plot the predicted contribution curves
      lines(x_values, eta_contribution_values[,s], lwd = lwd)
      
      # Approximate confidence bands
      if (color) {
        lines(ci_list$x, ci_list$lower, lty = 2, col = nice_colors[1], lwd = lwd)
        lines(ci_list$x, ci_list$upper, lty = 2, col = nice_colors[1], lwd = lwd)        
      } else {
        lines(ci_list$x, ci_list$lower, lty = 3, col = 'black', lwd = lwd)
        lines(ci_list$x, ci_list$upper, lty = 3, col = 'black', lwd = lwd)
      }

      # Add points to see where the training data is located
      points(jitter(model$x[,s]), rep(ymin, n), cex =.5, pch ="|", col = "darkgrey")
      
      # Add the empirical pointwise bootstrap confidence bands if provided
      if (!is.null(boot_models)) {
        if (color) {
          lines(sort(model$x[,s]), boot_bands[[s]][1,], lty = 2, col = nice_colors[4], lwd = lwd)
          lines(sort(model$x[,s]), boot_bands[[s]][2,], lty = 2, col = nice_colors[4], lwd = lwd)
        } else {
          lines(sort(model$x[,s]), boot_bands[[s]][1,], lty = 2, col = "black", lwd = lwd)
          lines(sort(model$x[,s]), boot_bands[[s]][2,], lty = 2, col = "black", lwd = lwd)
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


GAMBoost_splines_inside_conf_bands = function(model, boot_models, cn, training_data,
                                              after_iter = NA, type = c('smooth', 'step'),
                                              colors = TRUE, ymin = NULL, ymax = NULL,
                                              figure_type = 1){
  # model: A GAMBoost model fitted to the best iteration
  # boot_models: A list of GAMBoost models fitted to bootstrap 
  #              samples from the training data used on fitting 'model'.
  # cn: The scalar used in generating the linear predictor eta in the training data.
  # training data: The full generated training data list which the model is fitted to
  # After_iter: the optimal iteration
  # type: If the true lines are the smooth or the step functions.
  
  if (is.na(after_iter)) {
    after_iter = model$stepno
  } else {
    stopifnot(after_iter >= 0 & after_iter <= model$stepno)
  }
  
  # Get the dimensions
  n = model$n
  p = ncol(model$x)
  
  # Find the x_ij values we want to evaluate
  x_values_sort = apply(model$x, 2, sort)
  
  # Find the true values of the true f_j functions.
  true_f_values = matrix(0, nrow = n, ncol = p)
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
  eta_contribution_values = predict(model, newdata = x_values_sort, 
                                    at.step = after_iter, type = 'terms')
  
  # Compute the empirical confidence bands
  boot_bands = GAMBoost_splines_boot_cofidence_bands(boot_models, model$x, after_iter)
  
  # End coverage 
  coverage = list('Approximate'    = matrix(NA, ncol = p, nrow = n),
                  'Empirical'      = matrix(NA, ncol = p, nrow = n),
                  'Approximate_mu' = rep(NA, n),
                  'Empirical_mu'   = rep(NA, n))
  
  for (s in 1:p) {
    ### Calculate the approximate pointwise confidence bands.
    {
      # Calculate the approximate pointwise confidence bands.
      ci_list = GAMBoost_calc_confidence_bands(model, s, after_iter, 1)
      
      if (is.null(ci_list)){
        inside_cofidence_band_approx  = rep(NA, n)
      } else {
        # Get the approximate confidence bands
        cofidence_band_approx = rbind(ci_list$lower, ci_list$upper)
        
        # Find their coverage
        inside_cofidence_band_approx = (cofidence_band_approx[1,] <= true_f_values[,s]) &
          (true_f_values[,s] <= cofidence_band_approx[2,]) 
      }
    }
    
    ### Look at the empirical confidence bands
    {
      # Find their coverage
      inside_cofidence_band_empirical = (boot_bands[[s]][1,] <= true_f_values[,s]) &
        (true_f_values[,s] <= boot_bands[[s]][2,])
    }
    
    
    coverage[['Approximate']][,s] = inside_cofidence_band_approx
    coverage[['Empirical']][,s]   = inside_cofidence_band_empirical
  }
  
  
  # Predict the mu
  predicted_mu = model$family$linkinv(model$eta[ , after_iter+1])
  splines_pointwise_var = sqrt(diag(model$hatmatrix %*%
                                      diag(model$family$variance(predicted_mu)) %*%
                                      t(model$hatmatrix)))
  splines_lower = predicted_mu - 2*splines_pointwise_var
  splines_upper = predicted_mu + 2*splines_pointwise_var
  coverage[['Approximate_mu']] = splines_lower <= training_data$mu & training_data$mu <= splines_upper
  
  
  # splines_order = order(predicted_mu)
  # plot(splines_lower[splines_order], type = 'l', lty = 2, ylim = c(min(splines_lower), max(splines_upper)))
  # lines(splines_upper[splines_order], lty = 2)
  # lines(predicted_mu[splines_order], lty = 1)
  # points(training_data$mu[splines_order], col = 'steelblue')
  
  
  #
  boot_predictions = matrix(NA, ncol = n, nrow = length(boot_models))
  {
    for (m in 1:length(boot_models)) {
      curr_mod = boot_models[[m]]
      boot_predictions[m,] = predict(curr_mod, newdata = model$x, at.step = after_iter, type = 'response')
    }
  }
  boot_ci = apply(boot_predictions, 2, quantile, c(0.025, 0.975))
  coverage[['Empirical_mu']] = boot_ci[1,] <= training_data$mu & training_data$mu <= boot_ci[2,]
  
  # splines_order = order(splines_pred_mu)
  # lines()
  # lines(boot_ci[1,order(apply(boot_predictions, 2, mean))], lty = 2, col = 'blue')
  # lines(boot_ci[2,order(apply(boot_predictions, 2, mean))], lty = 2, col = 'blue')
  
  splinesorder = order(training_data$mu)
  boot_ci[1,splinesorder] <= training_data$mu[splinesorder] & training_data$mu[splinesorder] <= boot_ci[2,splinesorder]
  
  
  # Create the title
  if (training_data$type == 'smooth') {
    if (training_data$family$family == 'gaussian') {
      main_tit = expression(paste("Gaussian Smooth Splines: Empirical and Approximate CI for ", mu[i]))
    } else if (training_data$family$family == 'binomial') {
      main_tit = expression(paste("Binomial Smooth Splines: Empirical and Approximate CI for ", mu[i]))      
    } else {
      main_tit = expression(paste("Poisson Smooth Splines: Empirical and Approximate CI for ", mu[i]))     
    }
  } else {
    if (training_data$family$family == 'gaussian') {
      main_tit = expression(paste("Gaussian Stepwise Splines: Empirical and Approximate CI for ", mu[i]))
    } else if (training_data$family$family == 'binomial') {
      main_tit = expression(paste("Binomial Stepwise Splines: Empirical and Approximate CI for ", mu[i]))      
    } else {
      main_tit = expression(paste("Poisson Stepwise Splines: Empirical and Approximate CI for ", mu[i]))     
    }
  }
  if (is.null(ymin) | is.null(ymax)) {
    ymin = min(splines_lower)
    ymax = max(splines_upper)
  }
  
  if (figure_type == 1) {
    #par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
    plot(training_data$mu[splinesorder], type = 'l', col = 'darkgray', lwd = 1.5,
         ylim = c(ymin, ymax), xlab = '', ylab = '',
         main = main_tit, cex.main = 1.13)
    points(training_data$mu[splinesorder], col = 'darkgray', pch = 19, cex = 0.5)
    title(xlab = bquote("Observations in Increasing Order W.R.T. True"~mu[i]),
          ylab = bquote(mu[i]), line = 2)
    
    lines(predicted_mu[splinesorder], lty = 1, col = 'black')
    points(predicted_mu[splinesorder], col = 'black', pch = 19, cex = 0.4)
    if (colors) {
      legend("topleft", lty = c(1,1,2,2), col = c('darkgray', 'black', nice_colors[4], nice_colors[1]),
             legend = c(expression(paste("True ", mu[i])),
                        expression(paste("Predicted ", mu[i])),
                        'Empirical',
                        'Approximate'), cex = 0.75)
      lines(splines_upper[splinesorder], lty = 2, col = nice_colors[1])
      lines(splines_lower[splinesorder], lty = 2, col = nice_colors[1])
      lines(boot_ci[1,splinesorder], lty = 2, col = nice_colors[4])
      lines(boot_ci[2,splinesorder], lty = 2, col = nice_colors[4])    
    } else {
      legend("topleft", lty = c(1,1,2,3), col = c('darkgray', 'black', "#969696", 'black'),
             legend = c(expression(paste("True ", mu[i])),
                        expression(paste("Predicted ", mu[i])),
                        'Empirical',
                        'Approximate'), cex = 0.75)
      lines(splines_upper[splinesorder], lty = 3, col = 'black')
      lines(splines_lower[splinesorder], lty = 3, col = 'black')
      lines(boot_ci[1,splinesorder], lty = 2, col = "#969696")
      lines(boot_ci[2,splinesorder], lty = 2, col = "#969696")
    }
  } else {
    plot(1:n, training_data$mu[splinesorder], col = 'darkgray',
         ylim=c(ymin, ymax), pch=19, xlab="", ylab="",
         main = main_tit, cex.main = 1.13)
    title(xlab = bquote("Observations in Increasing Order W.R.T. True"~mu[i]),
          ylab = bquote(mu[i]), line = 2)
    points(1:n, predicted_mu[splinesorder], col = 'black', pch = 17)
    # hack: we draw arrows but with very special "arrowheads"
    arrows(1:n - 0.0, splines_lower[splinesorder],
           1:n - 0.0, splines_upper[splinesorder],
           length=0.02, angle=90, code=3, col = 'black', lwd = 2.5)
    arrows(1:n + 0.0, boot_ci[1,splinesorder],
           1:n + 0.0, boot_ci[2,splinesorder],
           length=0.02, angle=90, code=3, col = "#969696", lty = 1, lwd = 1)
    legend("topleft", lty = c(1,1,2,1),
           pch = NULL,
           col = c('darkgray', 'black', "#969696", 'black'),
           legend = c(expression(paste("True ", mu[i])),
                      expression(paste("Predicted ", mu[i])),
                      'Empirical',
                      'Approximate'), cex = 0.75)
  }

  
  # splinesorder = 1:10
  # plot(training_data$mu[splinesorder], type = 'p', col = 'gray', lwd = 5,
  #      ylim = c(min(splines_lower), max(splines_upper)), xlab = '', ylab = '',
  #      main = "Empirical and Approximate CI For Mean",  pch = 19)
  # title(xlab = "The then first observations", ylab = "Mean", line = 2)
  # legend("topleft", lty = 1, col = c(nice_colors[1], nice_colors[4]),
  #        legend = c('Approximate', 'Empirical'), cex = 0.8)
  # points(splines_upper[splinesorder], lty = 2, col = nice_colors[1], pch = 18)
  # points(splines_lower[splinesorder], lty = 2, col = nice_colors[1], pch = 18)
  # points(boot_ci[1,splinesorder], lty = 2, col = nice_colors[4], pch = 17)
  # points(boot_ci[2,splinesorder], lty = 2, col = nice_colors[4], pch = 17)
  # points(predicted_mu[splinesorder], lty = 1, col = 'black', pch = 19, lwd = 1)
  
  return(coverage)
}
