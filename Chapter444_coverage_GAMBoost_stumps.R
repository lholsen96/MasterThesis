###
setwd("C:/Users/lars9/OneDrive/Master/Finalized R codes to put on GitHub")

source('GAMBoost_stumps.R', local = TRUE)
source('GAMBoost_common.R', local = TRUE)

# rm(list=ls())
# Include libraries
library(parallel)
library(foreach)
library(snow)
library(doSNOW)

#
nice_colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)

colors = FALSE


B_outer = 200
B_inner = 200
stepno = 250

n_train = 100
n_test = 1000
p = 5

### Define possible distributions and function types
# types 
types = c('smooth', 'step', 'advancedstep')

# Set family
families = c('gaussian', 'binomial', 'poisson')

# Signal to noise ratios
stnrs = c(1,3,10)

# set penatly lambda
lambda = 4

for (family in families) {
  # Change from string to fucntion
  family = switch(family, gaussian = gaussian(), binomial = binomial(), poisson = poisson())
  
  for (type in types) {
    for (STNR in stnrs) {
      cat(sprintf("Family: %s \t Type: %s \t STNR: %d \n", family$family, type, STNR))
      # Lists of p matrix, Where the j'th matrix of dimension B_outer x n.
      # Each outer bootstrap is a row containing zeros and ones.
      # Zero means that for this iteration the true line/value is outside
      # the confidence band, while one represents that the value is inside.
      # Approximate based on the approximate CI from hat mat while empirical
      # comes from bootstrap. 
      inside_band_approximate = list()
      inside_band_empirical   = list()
      for (j in 1:p) {
        inside_band_approximate[[j]] = matrix(NA, nrow = B_outer, ncol = n_train)
        inside_band_empirical[[j]]   = matrix(NA, nrow = B_outer, ncol = n_train)
      }
      
      mu_inside_interval_approximate = matrix(NA, nrow = B_outer, ncol = n_train)
      mu_inside_interval_empirical   = matrix(NA, nrow = B_outer, ncol = n_train)
      
      
      
      # Vector to store the best iterations. 
      best_iterations = rep(NA, B_outer)
      best_iterations_dev = rep(NA, B_outer)
      
      # matrix to store times used
      times = matrix(NA, ncol = 3, nrow = B_outer)
      loop_times = matrix(NA, ncol = 3, nrow = B_outer)
      
      ### VERSION 2
      set.seed(16180340)
      training_data = create_data_combined(n_train, stnr = STNR, type = type, family = family)
      testing_data  = create_data_combined(n_test, cn = training_data$cn, type = type, family = family)
      for (bo in 1:B_outer){
        loop_timer = system.time({
          cat(sprintf("Currently in outer bootstrap number %d of %d. ", bo, B_outer))
          ### Setup: generate data
          # Randomly generate the response based on family and mu.
          if (family$family == "gaussian") {
            training_data$y = rnorm(n_train, mean=training_data$mu, sd=1)
          }
          if (family$family == "binomial") {
            training_data$y = rbinom(n_train, size = 1, prob = training_data$mu)
          }
          if (family$family == "poisson") {
            training_data$y = rpois(n_train, lambda = training_data$mu)
          }
          
          # Fit the model
          stumps_model = GAMBoost_stumps(X         = training_data$X, 
                                         y         = training_data$y,
                                         num_iter  = stepno,
                                         lambda    = lambda,
                                         family    = family,
                                         print_msg = FALSE)
          
          # Find best iter based on lowest evaluation error
          stumps_model_deviance_test = rep(NA, stepno+1)
          for (iter in 1:(stepno+1)) {
            # Get the predicted results
            predicted_response = GAMBoost_stumps_predict(stumps_model, X_new = testing_data$X,
                                                         after_iter = (iter-1), version = 'response')
            
            # Compute the deviance
            stumps_model_deviance_test[iter] = sum(family$dev.resids(testing_data$y,
                                                                     predicted_response,
                                                                     rep(1, n_test)))
            
          }
          # Find the optimal test deviance 
          best_iterations[bo] = which.min(stumps_model_deviance_test) - 1 # Subtract 1 since intercept is iter 0.
          best_iterations_dev[bo] = min(stumps_model_deviance_test)
          
          # Fit a new GAMBoost model to the best iteration to obtain the hat matrix
          # stumps_model = GAMBoost_stumps(X         = training_data$X, 
          #                                y         = training_data$y,
          #                                num_iter  = best_iterations[bo],
          #                                lambda    = lambda,
          #                                family    = family)
          
          stumps_model$hat_mat_last = stumps_model$hat_mat[[best_iterations[bo]+1]]
          
          ##### Want to create B_inner Bootstrap samples
          ### Parallel bootstraps
          #Find the number of available cores
          cores=detectCores() 
          
          # Set up the clusters
          cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
          registerDoSNOW(cl)
          
          # Iterate over the B_inner bootstrap samples 
          parallel_timer = system.time({
            stumps_boot_models = foreach(b = 1:B_inner) %dopar% {
              # The indices for this bootstrap sample
              boot_idx = sample(1:n_train, n_train, replace = TRUE)
              
              # Get the relevant data for b
              boot_X_train = training_data$X[boot_idx, ]
              boot_y_train = training_data$y[boot_idx]
              
              # Fit the model (250 iterations = 2MB)
              boot_pstump  = GAMBoost_stumps(X = boot_X_train, 
                                             y = boot_y_train,
                                             num_iter  = best_iterations[bo],
                                             lambda    = lambda,
                                             family    = family, 
                                             print_msg = FALSE,
                                             tiny_return = TRUE)
              }
          })
          
          # Stop the cluster
          parallel::stopCluster(cl)
          
          # record the timer
          times[bo,] = parallel_timer[1:3]
          
          if (bo <= 5) {
            png(paste('stumps_', type, '_', family$family, '_stnr_', STNR,
                      '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_',
                      B_outer, '_fitted_func_', bo, '.png', sep = ''),
                width = 3500, height = 1250, res = 350)
          }
          # Just to look at the fitted model
          GAMBoost_stumps_plot_predictor_contribution(stumps_model,
                                                      type = type, 
                                                      boot_models = stumps_boot_models,
                                                      cn = training_data$cn,
                                                      after_iter = best_iterations[bo],
                                                      save_filename = NULL,
                                                      xmin = -1,
                                                      xmax = 1,
                                                      ymin = NULL,
                                                      ymax = NULL,
                                                      lwd = 1,
                                                      color = FALSE) 
          
          if (bo <= 5) {
            dev.off()
            png(paste('stumps_', type, '_', family$family, '_stnr_', STNR,
                      '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_',
                      B_outer, '_fitted_mu_', bo, '.png', sep = ''),
                width = 3500, height = 1250, res = 350)
          }
          par(mfrow=c(1,1), mar = c(3.5, 3.5, 2, 1))
          # check if the true line is inside the bands
          #par(mfrow=c(2,1), mar = c(3.5,3.5,2,1))
          coverage_temp = GAMBoost_stumps_inside_conf_bands(model = stumps_model,
                                                            boot_models = stumps_boot_models,
                                                            cn = training_data$cn, 
                                                            training_data = training_data,
                                                            after_iter = best_iterations[bo],
                                                            type = type,
                                                            colors = FALSE)
          if (bo <= 5) {
            dev.off()
          }
          
          for (j in 1:p){
            inside_band_approximate[[j]][bo, ] = coverage_temp$Approximate[,j]
            inside_band_empirical[[j]][bo, ]   = coverage_temp$Empirical[,j]
          }
          
          # Add the mu hat coverage values
          mu_inside_interval_approximate[bo, ] = coverage_temp$Approximate_mu
          mu_inside_interval_empirical[bo, ]   = coverage_temp$Empirical_mu
        })
        
        loop_times[bo,] = loop_timer[1:3]
        cat(sprintf("Best iter: %3d  Time: %5.3g   Loop: %5.3g\n",
                    best_iterations[bo], times[bo,3], loop_times[bo,3]))
      }
      
      png(paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                '_coverage_final.png',  sep = ""),
          width = 3500, height = 1250, res = 350)
      
      par(mfrow=c(1,5), mar = c(3.5,3.5,2,1))
      for (s in 1:p) {
        # Find the coverage for the n_train values
        coverage_approx_temp = apply(inside_band_approximate[[s]], 2, sum, na.rm = TRUE) / 
          (B_outer - sum(is.na(inside_band_approximate[[s]][,1])))
        
        coverage_empir_temp = apply(inside_band_empirical[[s]], 2, sum, na.rm = TRUE) / 
          (B_outer - sum(is.na(inside_band_empirical[[s]][,1])))
        
        
        plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
        title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
        title(main = bquote('Coverage f'[.(s)]), line = 1)
        points(jitter(training_data$X[,s]), rep(0, n_train), cex =.65, pch ="|", col = "#D3D3D3")
        abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
        if (colors) {
          lines(sort(training_data$X[,s]), coverage_approx_temp, col = nice_colors[1])
          lines(sort(training_data$X[,s]), coverage_empir_temp,  col = nice_colors[4])          
        } else {
          lines(sort(training_data$X[,s]), coverage_approx_temp, lty = 1, lwd = 1.5, col = 'black')
          lines(sort(training_data$X[,s]), coverage_empir_temp,  lty = 2, lwd = 1.5, col = "#676767")
        }
        
        print(c(median(coverage_approx_temp), median(coverage_empir_temp)))
      }
      dev.off()
      
      png(paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                '_coverage_mu.png',  sep = ""),
          width = 3500, height = 1250, res = 350)
      {
        par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
        coverage_approx_mu = apply(mu_inside_interval_approximate, 2, sum, na.rm = TRUE) / 
          (B_outer - sum(is.na(mu_inside_interval_approximate[,1])))
        
        coverage_empircal_mu = apply(mu_inside_interval_empirical, 2, sum, na.rm = TRUE) / 
          (B_outer - sum(is.na(mu_inside_interval_empirical[,1])))
        
        plot(NA, type = 'n', xlim = c(1, n_train), ylim = c(0,1), xlab = '', ylab = '')
        title(xlab = bquote("Observation"), ylab =  bquote("Coverage"), line = 2.25)
        title(main = bquote("Coverage" ~ mu), line = 1)
        abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
        
        if (colors) {
          lines(1:n_train, coverage_approx_mu,    col = nice_colors[1])
          lines(1:n_train, coverage_empircal_mu,  col = nice_colors[4])
          points(1:n_train, coverage_approx_mu,   col = nice_colors[1], pch = 19)
          points(1:n_train, coverage_empircal_mu, col = nice_colors[4], pch = 19) 
        } else {
          lines(1:n_train,  coverage_approx_mu,   col = 'black',   lty = 1,  lwd = 1.5)
          points(1:n_train, coverage_approx_mu,   col = 'black',   pch = 19, cex = 0.5)
          lines(1:n_train,  coverage_empircal_mu, col = "#676767", lty = 2,  lwd = 1.5)
          points(1:n_train, coverage_empircal_mu, col = "#676767", pch = 19, cex = 0.5)
        }
        
        print(c(median(coverage_approx_mu), median(coverage_empircal_mu)))
      }
      dev.off()
      
      
      saveRDS(list('times' = times, 'training' = training_data, 'testing' = testing_data,
                   'inside_band_app' = inside_band_approximate, 'inside_band_emp'= inside_band_empirical,
                   'best_iters' = best_iterations, 'family' = family, 'type' = type, 'lambda' = lambda,
                   'stnr' = STNR, 'B_out' = B_outer, 'B_in' = B_inner, 'stepno' = stepno,
                   'best_iters_dev' = best_iterations_dev, 'loop_times' = loop_times, 
                   'mu_inside_interval_app' = mu_inside_interval_approximate,
                   'mu_inside_interval_emp' = mu_inside_interval_empirical)
              ,
              paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                    '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                    "_final2", sep='')
      )
    }
  }
}




#######################################
##### CODE TO SKIP THE BOOTSTRAPS #####
#######################################
B_outer = 200
stepno = 800

n_train = 100
n_test = 1000
p = 5

### Define possible distributions and function types
# types 
types = c('smooth', 'step', 'advancedstep')

# Signal to noise ratios
stnrs = c(1,3,10)

# Set family
families = c('gaussian', 'binomial', 'poisson')

# set penatly lambda
lambda = 100

for (family in families) {
  # Change from string to fucntion
  family = switch(family, gaussian = gaussian(), binomial = binomial(), poisson = poisson())
  
  for (type in types) {
    for (STNR in stnrs) {
      cat(sprintf("Family: %s \t Type: %s \t STNR: %d \n", family$family, type, STNR))
  
      # Data structures to store the results
      inside_band_approximate = list()
      for (j in 1:p) {
        inside_band_approximate[[j]] = matrix(NA, nrow = B_outer, ncol = n_train)
      }
      mu_inside_interval_approximate = matrix(NA, nrow = B_outer, ncol = n_train)

      # Vector to store the best iterations. 
      best_iterations = rep(NA, B_outer)
      best_iterations_dev = rep(NA, B_outer)
      
      ### VERSION 2
      set.seed(16180340)
      training_data = create_data_combined(n_train, stnr = STNR, type = type, family = family)
      testing_data  = create_data_combined(n_test, cn = training_data$cn, type = type, family = family)
      
      ##### Want to create B_inner Bootstrap samples
      ##### Parallel bootstraps
      # Find the number of available cores
      cores=detectCores()

      # Set up the clusters
      cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
      registerDoSNOW(cl)
      
      # Create a progress bar
      pb = txtProgressBar(max = B_outer, style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      opts = list(progress = progress)
      
      parallel_timer = system.time({
        list_of_results = foreach(bo = 1:B_outer, .options.snow = opts) %dopar% {
          #cat(sprintf("Currently in outer bootstrap number %d of %d. ", bo, B_outer))
          ### Setup: generate data
          # Randomly generate the response based on family and mu.
          if (family$family == "gaussian") {
            training_data$y = rnorm(n_train, mean=training_data$mu, sd=1)
          }
          if (family$family == "binomial") {
            training_data$y = rbinom(n_train, size = 1, prob = training_data$mu)
          }
          if (family$family == "poisson") {
            training_data$y = rpois(n_train, lambda = training_data$mu)
          }
          
          # Fit the model
          stumps_model = GAMBoost_stumps(X         = training_data$X, 
                                         y         = training_data$y,
                                         num_iter  = stepno,
                                         lambda    = lambda,
                                         family    = family)
          
          # Find best iter based on lowest evaluation error
          stumps_model_deviance_test = rep(NA, stepno+1)
          for (iter in 1:(stepno+1)) {
            # Get the predicted results
            predicted_response = GAMBoost_stumps_predict(stumps_model, X_new = testing_data$X,
                                                         after_iter = (iter-1), version = 'response')
            
            # Compute the deviance
            stumps_model_deviance_test[iter] = sum(family$dev.resids(testing_data$y,
                                                                     predicted_response,
                                                                     rep(1, n_test)))
            
          }
          # Find the optimal test deviance 
          best_iteration = which.min(stumps_model_deviance_test) - 1 # Subtract 1 since intercept is iter 0.
          best_iteration_dev = min(stumps_model_deviance_test)
          
          # Tacky way of making GAMBoost_stumps_inside_conf_bands
          # work, since it uses hat_mat_last to compute the CI.
          stumps_model$hat_mat_last = stumps_model$hat_mat[[best_iteration+1]]
          
          # check if the true line is inside the bands
          coverage_temp = GAMBoost_stumps_inside_conf_bands(model = stumps_model,
                                                            boot_models = NA,
                                                            cn = training_data$cn, 
                                                            training_data = training_data,
                                                            after_iter = best_iteration,
                                                            type = type)

          list('best_iteration' = best_iteration, 
               'best_iteration_dev' = best_iteration_dev,
               'mu_inside' = coverage_temp$Approximate_mu,
               'fs_inside' = coverage_temp$Approximate)
        }
      })
      
      # Stop the progress bar
      close(pb)
      # Stop the cluster
      parallel::stopCluster(cl)
      
      ### Convert the list of results to the desired objects
      for (bo in 1:B_outer) {
        curr_result = list_of_results[[bo]]
        best_iterations[bo] = curr_result$best_iteration
        best_iterations_dev[bo] = curr_result$best_iteration_dev
        mu_inside_interval_approximate[bo, ] = curr_result$mu_inside
        for (j in 1:p){
          inside_band_approximate[[j]][bo, ] = curr_result$fs_inside[,j]
        }
      }
      
      png(paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                '_lambda_', lambda, '_bootouter_', B_outer, 
                '_coverage_without_empirical.png',  sep = ""),
          width = 3500, height = 1250, res = 350)
      
      par(mfrow=c(1,5), mar = c(3.5,3.5,2,1))
      for (s in 1:p) {
        # Find the coverage for the n_train values
        coverage_approx_temp = apply(inside_band_approximate[[s]], 2, sum, na.rm = TRUE) / 
          (B_outer - sum(is.na(inside_band_approximate[[s]][,1])))
        
        
        plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
        title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
        lines(sort(training_data$X[,s]), coverage_approx_temp, col = nice_colors[1])
        abline(h = 0.95, lty = 4)
        
        print(c(mean(coverage_approx_temp), median(coverage_approx_temp)))
      }
      dev.off()
      
      png(paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                '_lambda_', lambda, '_bootouter_', B_outer, 
                '_coverage_mu_without_empirical.png',  sep = ""),
          width = 3500, height = 1250, res = 350)
      {
        par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
        coverage_approx_mu = apply(mu_inside_interval_approximate, 2, sum, na.rm = TRUE) / 
          (B_outer - sum(is.na(mu_inside_interval_approximate[,1])))
        
        plot(NA, type = 'n', xlim = c(1, n_train), ylim = c(0,1), xlab = '', ylab = '')
        title(xlab = bquote("Observation"), ylab =  bquote("Coverage"), line = 2.25)
        title(main = bquote("Coverage of the predicted" ~ mu), line = 1)
        lines(1:n_train, coverage_approx_mu, col = nice_colors[1])
        points(1:n_train, coverage_approx_mu, col = nice_colors[1], pch = 19)

        abline(h = 0.95, lty = 4)
        print(c(mean(coverage_approx_mu), median(coverage_approx_mu)))
      }
      dev.off()
      
      
      saveRDS(list('training' = training_data, 'testing' = testing_data,
                   'inside_band_app' = inside_band_approximate, 
                   'best_iters' = best_iterations, 'family' = family,
                   'type' = type, 'lambda' = lambda,
                   'stnr' = STNR, 'B_out' = B_outer, 'stepno' = stepno,
                   'best_iters_dev' = best_iterations_dev, 'time' = parallel_timer, 
                   'mu_inside_interval_app' = mu_inside_interval_approximate)
              ,
              paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                    '_lambda_', lambda, '_bootouter_', B_outer, 
                    "_final_without_empirical", sep='')
      )
      
      # 
      # hist(best_iterations_dev)
      # hist(best_iterations)
      
    }
  }
}






#########################################################
##### Code to open saved lists and generate figures #####
#########################################################
setwd("C:/Users/lars9/OneDrive/Master/R_codes/GAMBoost_final/Resultater/Saves")

# Some parameters needed to read the files
# All possible cobinations of types, families, SNTRs and lambdas
types = c('smooth', 'step')
stnrs = c(1,3,10)
families = c('gaussian', 'binomial', 'poisson')
lambdas = c(2, 4, 10, 25)
B_outer = 200
B_inner = 200

# Create a list to store all the results and extract them fromt the files
# Around xxxMb
stumps = list()
for (family in families) {
  stumps[[family]] = list()
  for (type in types) {
    stumps[[family]][[type]] = list()
    for (STNR in stnrs) {
      stumps[[family]][[type]][[paste('STNR',STNR,sep = '')]] = list()
      for (lambda in lambdas) {
        stumps[[family]][[type]][[paste('STNR',STNR,sep = '')]][[paste('lam',lambda,sep = '')]] =
          readRDS(paste("stumps_", type, "_", family, "_stnr_", STNR,
                        '_lambda_', lambda, '_bootinner_', B_inner,
                        '_bootouter_', B_outer, "_final2", sep=''))
      }
    }
  }
}





