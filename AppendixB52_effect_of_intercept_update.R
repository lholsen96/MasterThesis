# rm(list=ls())

library(parallel)
#library(doParallel)
library(foreach)
library(snow)
library(doSNOW)


setwd("C:/Users/lars9/OneDrive/Master/Finalized R codes to put on GitHub")
source('GAMBoost_stumps.R', local = TRUE)
source('GAMBoost_splines.R', local = TRUE)
source('GAMBoost_common.R', local = TRUE)
source('GAMBoost_stumps_with_intercept_update.R', local = TRUE)
setwd("C:/Users/lars9/OneDrive/Master/R_codes/GAMBoost_final")



GAMBoost_intercept_update_effect = function(family = poisson(), type = 'smooth', stnr = 1,
                                            n_train = 100, n_test  = 1000,
                                            num_iter = 500, lambda  = 4,
                                            seed_train = 2020, seed_test = 2030,
                                            repetitions = 20) {
  ### Function that fit GAMBoost with and without intercept update to the same
  ### data and compute the test deviance.
  
  ## Parallel fitting of model with different lambdas
  # Find the number of available cores
  cores=detectCores() 
  
  # Set up the clusters
  cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
  registerDoSNOW(cl)
  
  # Create a progress bar
  pb = txtProgressBar(max = repetitions, style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  parallel_timer = system.time({
    # Iterate over the penalty factors
    parallel_list = foreach(i = 1:repetitions, 
                            .export = c('create_data_combined', 'GAMBoost_stumps',
                                        'GAMBoost_stumps_with_intercept_update',
                                        'GAMBoost_stumps_predict', 
                                        'GAMBoost_stumps_inter_update_predict',
                                        'GAMBoost_stumps_create_total_design_matrix',
                                        'create_total_design_matrix'),
                            .options.snow = opts) %dopar% {
      
      ### Setup: generate data
      {
        training_data = create_data_combined(n_train, seed_number = seed_train + 7*i, stnr = stnr,
                                             type = type, family = family)
        
        testing_data  = create_data_combined(n_test, seed_number = seed_test + 7*i, stnr = stnr,
                                             type = type, family = family)
      }
      
      ### fit the models
      stumps_model = GAMBoost_stumps(X           = training_data$X, 
                                     y           = training_data$y,
                                     num_iter    = num_iter,
                                     lambda      = lambda,
                                     family      = family,
                                     print_msg   = FALSE,
                                     tiny_return = TRUE)
      
      stumps_model_intercept = 
        GAMBoost_stumps_with_intercept_update(X           = training_data$X,
                                              y           = training_data$y,
                                              num_iter    = num_iter,
                                              lambda      = lambda,
                                              family      = family,
                                              print_msg   = TRUE,
                                              tiny_return = TRUE)
      
      
      ### Get the test deviance
      stumps_model_deviance_test = rep(NA, num_iter+1)
      stumps_model_deviance_test_inter = rep(NA, num_iter+1)
      for (iter in 1:(num_iter+1)) {
        ## Without the intercept
        # Get the predicted results
        predicted_response = GAMBoost_stumps_predict(model = stumps_model,
                                                     X_new = testing_data$X,
                                                     after_iter = (iter-1),
                                                     version = 'response')
        
        # Compute the deviance
        stumps_model_deviance_test[iter] = sum(family$dev.resids(testing_data$y,
                                                                 predicted_response,
                                                                 rep(1, n_test)))
        
        ## With the intercept
        # Get the predicted results
        predicted_response = GAMBoost_stumps_inter_update_predict(model = stumps_model_intercept,
                                                                  X_new = testing_data$X,
                                                                  after_iter = (iter-1),
                                                                  version = 'response')
        
        # Compute the deviance
        stumps_model_deviance_test_inter[iter] = sum(family$dev.resids(testing_data$y,
                                                                       predicted_response,
                                                                       rep(1, n_test)))
        
      }
      
      ### Save model and data in list and return it
      list('GAMBoost_stumps' = stumps_model,
           'GAMBoost_stumps_intercept' = stumps_model_intercept,
           'stumps_deviance_test' = stumps_model_deviance_test,
           'stumps_intercept_deviance_test' = stumps_model_deviance_test_inter,
           'training_data' = training_data,
           'test_data' = testing_data,
           'repetition' = i
      )
    }
  })
  
  # Stop the progress bar
  close(pb)
  # Stop the cluster
  parallel::stopCluster(cl)
  
  save_list = list('models' = parallel_list,
                   'lambda' = lambda,
                   'num_iter' = num_iter,
                   'type' = type,
                   'stnr' = stnr,
                   'family' = family$family,
                   'repetitions' = repetitions,
                   'n_train' = n_train,
                   'n_test' = n_test,
                   'seed_train' = seed_train,
                   'seed_test' = seed_test,
                   'parallel_timer' = parallel_timer)
  
  
  saveRDS(save_list,
          file = paste('GAMBoost_stumps_with_and_without_intercept_', family$family, '_', type,
                       '_stnr_', stnr, '_lambda_', lambda, '_numiter_', num_iter,
                       '_repetitions_', repetitions, sep = ''))
  
  return(save_list)
}




# Parameters
n_train = 100
n_test  = 1000
num_iter = 500
lambda  = 4
type = 'smooth'
seed_train = 2020
seed_test = 2030
repetitions = 25
stnr = 1

######################################
##### FIT GAUSSIAN WITH STNR = 1 #####
######################################
gaussian1 = GAMBoost_intercept_update_effect(family = gaussian(), stnr = stnr,
                                            type = type, n_train = n_train, n_test = n_test,
                                            num_iter = num_iter, lambda  = lambda,
                                            seed_train = seed_train, seed_test = seed_test,
                                            repetitions = repetitions)

gaussian1 = readRDS(paste('GAMBoost_stumps_with_and_without_intercept_', 'gaussian', '_', type,
                          '_stnr_', stnr, '_lambda_', lambda, '_numiter_', num_iter,
                          '_repetitions_', repetitions, sep = ''))

#####################################
##### FIT POISSON WITH STNR = 1 #####
#####################################
poisson1 = GAMBoost_intercept_update_effect(family = poisson(), stnr = stnr,
                                            type = type, n_train = n_train, n_test = n_test,
                                            num_iter = num_iter, lambda  = lambda,
                                            seed_train = seed_train, seed_test = seed_test,
                                            repetitions = repetitions)

poisson1 = readRDS(paste('GAMBoost_stumps_with_and_without_intercept_', 'poisson', '_', type,
                          '_stnr_', stnr, '_lambda_', lambda, '_numiter_', num_iter,
                          '_repetitions_', repetitions, sep = ''))

######################################
##### FIT BINOMIAL WITH STNR = 1 #####
######################################
binomial1 = GAMBoost_intercept_update_effect(family = binomial(), stnr = stnr,
                                            type = type, n_train = n_train, n_test = n_test,
                                            num_iter = num_iter, lambda  = lambda,
                                            seed_train = seed_train, seed_test = seed_test,
                                            repetitions = repetitions)
binomial1 = readRDS(paste('GAMBoost_stumps_with_and_without_intercept_', 'binomial', '_', type,
                          '_stnr_', stnr, '_lambda_', lambda, '_numiter_', num_iter,
                          '_repetitions_', repetitions, sep = ''))

##### Plot the curves
### Create the curves
{
  gaussian1_mean_test_dev = 
    apply(sapply(1:repetitions, function(x) gaussian1$models[[x]]$stumps_deviance_test),
          1,
          mean)
  
  gaussian1_mean_test_dev_inter = 
    apply(sapply(1:repetitions, function(x) gaussian1$models[[x]]$stumps_intercept_deviance_test),
          1,
          mean)
  
  poisson1_mean_test_dev = 
    apply(sapply(1:repetitions, function(x) poisson1$models[[x]]$stumps_deviance_test),
          1,
          mean)
  
  poisson1_mean_test_dev_inter = 
    apply(sapply(1:repetitions, function(x) poisson1$models[[x]]$stumps_intercept_deviance_test),
          1,
          mean)
  
  binomial1_mean_test_dev = 
    apply(sapply(1:repetitions, function(x) binomial1$models[[x]]$stumps_deviance_test),
          1,
          mean)
  
  binomial1_mean_test_dev_inter = 
    apply(sapply(1:repetitions, function(x) binomial1$models[[x]]$stumps_intercept_deviance_test),
          1,
          mean)
}


png(paste('GAMBoost_stumps_with_and_without_intercept_', family$family, '_', type,
                '_stnr_', stnr, '_lambda_', lambda, '_numiter_', num_iter,
                '_repetitions_', repetitions, '.png', sep = ''), width = 3500, height = 1125, res = 350)
par(mfrow=c(1,3), mar = c(3.5,3.5,2,1))
max_iter = 100
{
  ### Gaussian Plot
  plot(NA, xlim = c(0, max_iter), 
       ylim = c(min(gaussian1_mean_test_dev[1:(max_iter+1)], gaussian1_mean_test_dev_inter[1:(max_iter+1)]),
                max(gaussian1_mean_test_dev[1:(max_iter+1)], gaussian1_mean_test_dev_inter[1:(max_iter+1)])),
       ylab = '', xlab = '')
  title(ylab = 'Test Deviance', xlab = 'Boosting Iteration', line = 2, cex.lab = 1.35)
  title(main = paste("Gaussian Smooth:", repetitions, " Repetitions", sep = ''), cex.main = 1.5)
  lines(0:max_iter, gaussian1_mean_test_dev[1:(max_iter+1)])
  lines(0:max_iter, gaussian1_mean_test_dev_inter[1:(max_iter+1)], lty = 2)
  legend('topright', legend = c('With', 'Without'), title = 'Intercept Update', lty = c(2,1))
  
  
  ### Binomial plot
  plot(NA, xlim = c(0, max_iter), 
       ylim = c(min(binomial1_mean_test_dev[1:(max_iter+1)], binomial1_mean_test_dev_inter[1:(max_iter+1)]),
                max(binomial1_mean_test_dev[1:(max_iter+1)], binomial1_mean_test_dev_inter[1:(max_iter+1)])),
       ylab = '', xlab = '')
  title(ylab = 'Test Deviance', xlab = 'Boosting Iteration', line = 2, cex.lab = 1.35)
  title(main = paste("Binomial Smooth:", repetitions, " Repetitions", sep = ''), cex.main = 1.5)
  lines(0:max_iter, binomial1_mean_test_dev[1:(max_iter+1)])
  lines(0:max_iter, binomial1_mean_test_dev_inter[1:(max_iter+1)], lty = 2)
  legend('topright', legend = c('With', 'Without'), title = 'Intercept Update', lty = c(2,1))
  
  
  ### Poisson Plot
  plot(NA, xlim = c(0, max_iter), 
       ylim = c(min(poisson1_mean_test_dev[1:(max_iter+1)], poisson1_mean_test_dev_inter[1:(max_iter+1)]),
                max(poisson1_mean_test_dev[1:(max_iter+1)], poisson1_mean_test_dev_inter[1:(max_iter+1)])),
       ylab = '', xlab = '')
  title(ylab = 'Test Deviance', xlab = 'Boosting Iteration', line = 2, cex.lab = 1.35)
  title(main = paste("Poisson Smooth:", repetitions, " Repetitions", sep = ''), cex.main = 1.5)
  lines(0:max_iter, poisson1_mean_test_dev[1:(max_iter+1)])
  lines(0:max_iter, poisson1_mean_test_dev_inter[1:(max_iter+1)], lty = 2)
  legend('topright', legend = c('With', 'Without'), title = 'Intercept Update', lty = c(2,1))
}
dev.off()











### Create the title
if (type == 'smooth') {
  if (family$family == 'gaussian') {
    main_tit = 'Gaussian Smooth'
  } else if (family$family == 'binomial') {
    main_tit = 'Binomial Smooth'
  } else if (family$family == 'poisson') {
    main_tit = 'Poisson Smooth'
  }
} else if (type == 'step') {
  if (family$family == 'gaussian') {
    main_tit = 'Gaussian Stepwise'
  } else if (family$family == 'binomial') {
    main_tit = 'Binomial Stepwise'
  } else if (family$family == 'poisson') {
    main_tit = 'Poisson Stepwise'
  }
}