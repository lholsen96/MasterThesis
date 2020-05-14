### File containing code to demonstrate that both m_stop and lambda 
### are important tuning parameters.

library(parallel)
#library(doParallel)
library(foreach)
library(snow)
library(doSNOW)


setwd("C:/Users/lars9/OneDrive/Master/Finalized R codes to put on GitHub")
source('GAMBoost_stumps.R', local = TRUE)
source('GAMBoost_splines.R', local = TRUE)
source('GAMBoost_common.R', local = TRUE)
setwd("C:/Users/lars9/OneDrive/Master/R_codes/GAMBoost_final")

family = gaussian()

# Parameters
n_train = 100
n_test  = 1000
cn = 3
num_iter = list('splines' = 250, 'stumps' = 250)
penalty  = list('splines' = 300, 'stumps' = 4)
penalty  = list('smooth' = list('splines' = 300, 'stumps' = 4),
                'step'   = list('splines' = 300, 'stumps' = 4))

# Number of bootstrap samples for each model
B = list('splines' = 1000, 'stumps' = 1000)

seed_train = 2020
seed_test = 2030
### Setup: generate data
{
  training_data_smooth = create_data_combined(n_train, seed_number = seed_train, cn = cn,
                                              type = 'smooth', family = family)
  training_data_step = create_data_combined(n_train, seed_number = seed_train, cn = cn,
                                            type = 'step', family = family)
  
  testing_data_smooth  = create_data_combined(n_test, seed_number = seed_test, cn = cn,
                                              type = 'smooth', family = family)
  testing_data_step  = create_data_combined(n_test, seed_number = seed_test, cn = cn,
                                            type = 'step', family = family)  
}

#####################################################
##### Look at best penalty for gaussian splines #####
#####################################################
splines_lambdas = c(0.1, 1, 5, 10, 20, 30, 40, seq(50, 1000, 50))
num_iter = 1000
seed_train = 2020
seed_test = 2030

outer_repetitions = 10
splines_results = list()

for (o in 1:outer_repetitions) {
  cat(sprintf('Outer iteration: %4d\n', o))
  
  splines_results[[o]] = list()
  #splines_model_deviance_test_smooth_mat = matrix(NA, ncol = num_iter+1, nrow = length(splines_lambdas))
  #splines_model_deviance_test_step_mat   = matrix(NA, ncol = num_iter+1, nrow = length(splines_lambdas))
  
  # Generate the new data
  {
    training_data_smooth = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, cn = cn,
                                                type = 'smooth', family = family)
    training_data_step = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, cn = cn,
                                              type = 'step', family = family)
    
    testing_data_smooth  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = cn,
                                                type = 'smooth', family = family)
    testing_data_step  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = cn,
                                              type = 'step', family = family)  
  }
  
  ## Parallel fitting of model with different lambdas
  # Find the number of available cores
  cores=detectCores() 
  
  # Set up the clusters
  cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
  registerDoSNOW(cl)
  
  # Create a progress bar
  pb = txtProgressBar(max = lenght(splines_lambdas), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  parallel_timer = system.time({
    # Iterate over the penalty factors
    parallel_list = foreach(i = 1:length(splines_lambdas),
                            .export = 'GAMBoost', .options.snow = opts) %dopar% {
      # Get teh current lambda 
      cur_lambda = splines_lambdas[i]
      #cat(sprintf('Lambda: %5.1f   ', cur_lambda))
      
      ### Fit the splines
      splines_model_smooth = GAMBoost(x        = training_data_smooth$X, 
                                      y        = training_data_smooth$y,
                                      stepno   = num_iter,
                                      penalty  = cur_lambda,
                                      family   = family,
                                      calc.hat = TRUE,
                                      calc.se  = TRUE)
      
      splines_model_step = GAMBoost(x        = training_data_step$X, 
                                    y        = training_data_step$y,
                                    stepno   = num_iter,
                                    penalty  = cur_lambda,
                                    family   = family,
                                    calc.hat = TRUE,
                                    calc.se  = TRUE)
      
      # Arrays to store data
      splines_model_deviance_test_smooth_arr = rep(NA, num_iter+1)
      splines_model_deviance_test_step_arr   = rep(NA, num_iter+1)
      
      
      # Iterate over all boosting iterations and find the best iteration based on a test set.
      for (iter in 1:(num_iter+1)) {
        ## Smooth
        # Compute the deviance
        splines_model_deviance_test_smooth_arr[iter] = sum(family$dev.resids(
          testing_data_smooth$y, 
          predict(splines_model_smooth, newdata = testing_data_smooth$X,
                  at.step = (iter-1), type = "response"),
          rep(1, n_test)))
        
        ## Step
        # Compute the deviance
        splines_model_deviance_test_step_arr[iter] = sum(family$dev.resids(
          testing_data_step$y, 
          predict(splines_model_step, newdata = testing_data_step$X,
                  at.step = (iter-1), type = "response"),
          rep(1, n_test)))
      }
      
      # List of objects to be returned
      ret_list = list('lambda' = cur_lambda,
                      'deviance_smooth' = splines_model_deviance_test_smooth_arr,
                      'deviance_step' = splines_model_deviance_test_step_arr)
      
      # splines_results[[o]]$splines_model_deviance_test_smooth_mat = splines_model_deviance_test_smooth_mat
      # splines_results[[o]]$splines_model_deviance_test_step_mat = splines_model_deviance_test_step_mat
      # 
      # cat(sprintf('Smooth: Dev = %9.3f Iter = %-4d     Step: Dev = %9.3f Iter = %-4d\n',
      #             min(splines_model_deviance_test_smooth_mat[i, ]),
      #             which.min(splines_model_deviance_test_smooth_mat[i, ]) - 1,
      #             min(splines_model_deviance_test_step_mat[i, ]),
      #             which.min(splines_model_deviance_test_step_mat[i, ]) - 1))
      }
    })
  
  # Stop the progress bar
  close(pb)
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Small printout to the user
  cat(sprintf('Parallel time = %7.2f\n', parallel_timer[3]))
  
  # Save the results to the splines_results list
  splines_results[[o]]$splines_model_deviance_test_smooth_mat = 
    t(sapply(1:length(splines_lambdas), function(x) parallel_list[[x]]$deviance_smooth))
  
  splines_results[[o]]$splines_model_deviance_test_step_mat = 
    t(sapply(1:length(splines_lambdas), function(x) parallel_list[[x]]$deviance_step))
}

# saveRDS(list('splines_results' = splines_results,
#              'splines_lambdas' = splines_lambdas),
#         file = paste('Splines_outer_', outer_repetitions, '_numiter_',
#                      num_iter, '_smooth_step_results', sep = ''))

temp_list = readRDS(paste('Splines_outer_', outer_repetitions, '_numiter_',
                          num_iter, '_smooth_step_results', sep = ''))

splines_results = temp_list$splines_results

splines_best_dev_smooth = sapply(1:outer_repetitions, function(x) apply(splines_results[[x]]$splines_model_deviance_test_smooth_mat, 1, min))
splines_best_dev_step   = sapply(1:outer_repetitions, function(x) apply(splines_results[[x]]$splines_model_deviance_test_step_mat, 1, min))

splines_best_iters_smooth = sapply(1:outer_repetitions, function(x) apply(splines_results[[x]]$splines_model_deviance_test_smooth_mat, 1, which.min)) - 1 
splines_best_iters_step   = sapply(1:outer_repetitions, function(x) apply(splines_results[[x]]$splines_model_deviance_test_step_mat, 1, which.min)) - 1

splines_avg_dev_smooth = apply(splines_best_dev_smooth, 1, mean)
splines_avg_dev_step   = apply(splines_best_dev_step, 1, mean)

splines_avg_iter_smooth = apply(splines_best_iters_smooth, 1, mean)
splines_avg_iter_step   = apply(splines_best_iters_step, 1, mean)


png(paste('Splines_outer_', outer_repetitions, '_numiter_',
          num_iter, '_smooth_step_results.png', sep = ""), width = 3500, height = 2250, res = 350)
par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
{
  ### Smooth
  # Plot 1
  plot(splines_lambdas, splines_avg_dev_smooth, type = 'l',
       xlab = '', ylab = '', ylim = c(min(splines_avg_dev_smooth), 1700))
  title(main = bquote('Smooth Splines'),  cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2.1, cex.lab = 1.35)
  text(500, 1695, bquote(lambda~'='~.(splines_lambdas[which.min(splines_avg_dev_smooth)])
                         ~' min ='~.(round(min(splines_avg_dev_smooth)))), cex = 1.35)
  
  # Plot 2
  plot(splines_lambdas,
       splines_avg_iter_smooth,  cex.main = 1.25,
       ylim = c(0, 1000), type = 'l', xlab = '', ylab = '')
  title(main = bquote('Smooth Splines'),  cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2.1, cex.lab = 1.35)
  
  ### Step
  # Plot 3
  plot(splines_lambdas, splines_avg_dev_step, type = 'l',
       xlab = '', ylab = '', ylim = c(min(splines_avg_dev_step), 2100))
  title(main = bquote('Stepwise Splines'),  cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2.1, cex.lab = 1.35)
  text(500, 2095, bquote(lambda~'='~.(splines_lambdas[which.min(splines_avg_dev_step)])
                         ~' min ='~.(round(min(splines_avg_dev_step)))), cex = 1.35)
  
  # Plot 4
  plot(splines_lambdas,
       splines_avg_iter_step,
       ylim = c(0, 1000), type = 'l', xlab = '', ylab = '')
  title(main = bquote('Stepwise Splines'),  cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2.1, cex.lab = 1.35)
}
dev.off()

# 
# 
# png(paste("stumps_min_test_deviance_different_penalty.png", sep = ""), width = 3500, height = 2250, res = 350)
# par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
# {
#   splines_best_iters_test_smooth = apply(splines_model_deviance_test_smooth_mat, 1, which.min)
#   splines_best_iters_test_step   = apply(splines_model_deviance_test_step_mat, 1, which.min)
#   splines_test_dev_smooth = splines_model_deviance_test_smooth_mat[cbind(1:length(splines_lambdas), splines_best_iters_test_smooth)]
#   splines_test_dev_step = splines_model_deviance_test_step_mat[cbind(1:length(splines_lambdas), splines_best_iters_test_step)]
#   
#   splines_lambda_smooth_min = splines_lambdas[which.min(splines_model_deviance_test_smooth_mat[cbind(1:length(splines_lambdas), splines_best_iters_test_smooth)])]
#   splines_lambda_step_min = splines_lambdas[which.min(splines_model_deviance_test_step_mat[cbind(1:length(splines_lambdas), splines_best_iters_test_step)])]
#   splines_smooth_min = round(min(splines_model_deviance_test_smooth_mat[cbind(1:length(splines_lambdas), splines_best_iters_test_smooth)]))
#   splines_step_min = round(min(splines_model_deviance_test_step_mat[cbind(1:length(splines_lambdas), splines_best_iters_test_step)]))
#   
#   plot(splines_lambdas,
#        splines_test_dev_smooth,
#        type = 'l', main = bquote('Smooth Splines'),
#        xlab = '', ylab = '', xlim = c(0,1000))
#   title(xlab = 'Penalty', ylab = 'Test Deviance', line = 2)
#   #abline(v = lambdas[min(which(best_iters_test_smooth >= 1500))], lty = 2)
#   text(500, 1300, bquote(lambda~'='~.(splines_lambda_smooth_min)~' min ='~.(splines_smooth_min)))
#   
#   plot(splines_lambdas,
#        splines_best_iters_test_smooth, xlim = c(0,2000),
#        ylim = c(0, 1500), type = 'l', main = bquote('Smooth Splines'), xlab = '', ylab = '')
#   title(xlab = 'Penalty', ylab = 'Best Iteration', line = 2)
#   #abline(v = lambdas[min(which(best_iters_test_smooth >= 1500))], lty = 2)
#   
#   
#   
#   plot(splines_lambdas, splines_test_dev_step,
#        main = bquote('Stepwise Splines'),
#        xlab = '', ylab = '', pch = 19, type = 'l')
#   title(xlab = 'Penalty', ylab = 'Test Deviance', line = 2)
#   #abline(v = 1000, lty = 2)
#   text(950,1865,bquote(lambda~'='~.(splines_lambda_step_min)~' min ='~.(splines_step_min)))
#   
#   plot(splines_lambdas, splines_best_iters_test_step, ylim = c(0, 1500),
#        type = 'l', main = bquote('Stepwise Splines'), xlab = '', ylab = '', pch = 19)
#   title(xlab = 'Penalty', ylab = 'Best Iteration', line = 2)
#   #abline(v = 1000, lty = 2)
# }
# dev.off()




####################################################
##### Look at best penalty for gaussian stumps #####
####################################################
stumps_lambdas = c(0.1, 1, 2, 5, 10, 15, 20, 25, 30, 40, seq(50, 400, 25))
num_iter = 1000
num_iter = list('smooth' = 1500, 'step' = 500)
seed_train = 2020
seed_test = 2030

outer_repetitions = 10
stumps_results = list()

for (o in 1:outer_repetitions) {
  cat(sprintf('Outer iteration: %4d\n', o))
  
  stumps_results[[o]] = list()
  #stumps_model_deviance_test_smooth_mat = matrix(NA, ncol = num_iter$smooth+1, nrow = length(stumps_lambdas))
  #stumps_model_deviance_test_step_mat   = matrix(NA, ncol = num_iter$step+1, nrow = length(stumps_lambdas))
  
  # Generate the new data
  {
    training_data_smooth = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, cn = cn,
                                                type = 'smooth', family = family)
    training_data_step = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, cn = cn,
                                              type = 'step', family = family)
    
    testing_data_smooth  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = cn,
                                                type = 'smooth', family = family)
    testing_data_step  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = cn,
                                              type = 'step', family = family)  
  }
  
  ## Parallel fitting of model with different lambdas
  # Find the number of available cores
  cores=detectCores() 
  
  # Set up the clusters
  cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
  registerDoSNOW(cl)
  
  # Create a progress bar
  pb = txtProgressBar(max = lenght(splines_lambdas), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  parallel_timer = system.time({
    # Iterate over the penalty factors
    parallel_list = foreach(i = 1:length(stumps_lambdas), .export = 'GAMBoost', .options.snow = opts) %dopar% {
      # Get teh current lambda 
      cur_lambda = stumps_lambdas[i]
      #cat(sprintf('Lambda: %5.1f   ', cur_lambda))
      
      ### Fit the stumps
      stumps_model_smooth = GAMBoost_stumps(X           = training_data_smooth$X, 
                                            y           = training_data_smooth$y,
                                            num_iter    = num_iter$smooth,
                                            lambda      = cur_lambda,
                                            family      = family,
                                            print_msg   = FALSE,
                                            tiny_return = TRUE)
      
      stumps_model_step = GAMBoost_stumps(X           = training_data_step$X, 
                                          y           = training_data_step$y,
                                          num_iter    = num_iter$step,
                                          lambda      = cur_lambda,
                                          family      = family,
                                          print_msg   = FALSE,
                                          tiny_return = TRUE)
      
      
      # Arrays to store data
      stumps_model_deviance_test_smooth_arr = rep(NA, num_iter$smooth+1)
      stumps_model_deviance_test_step_arr   = rep(NA, num_iter$step+1)
      
      
      # Iterate over all boosting iterations and find the best iteration based on a test set.
      for (iter in 1:(num_iter$smooth+1)) {
        ## Smooth
        # Get the predicted results
        predicted_response = GAMBoost_stumps_predict(stumps_model_smooth,
                                                     X_new = testing_data_smooth$X,
                                                     after_iter = (iter-1),
                                                     version = 'response')
        
        # Compute the deviance
        stumps_model_deviance_test_smooth_arr[iter] = 
          sum(family$dev.resids(testing_data_smooth$y,
                                predicted_response,
                                rep(1, n_test)))
      }
      
      for (iter in 1:(num_iter$step+1)) {
        ## Step
        # Get the predicted results
        predicted_response = GAMBoost_stumps_predict(stumps_model_step,
                                                     X_new = testing_data_step$X,
                                                     after_iter = (iter-1),
                                                     version = 'response')
        
        # Compute the deviance
        stumps_model_deviance_test_step_arr[iter] =
          sum(family$dev.resids(testing_data_step$y,
                                predicted_response,
                                rep(1, n_test)))
      }
      
      # List of objects to be returned
      ret_list = list('lambda' = cur_lambda,
                      'deviance_smooth' = stumps_model_deviance_test_smooth_arr,
                      'deviance_step' = stumps_model_deviance_test_step_arr)
    }                          
  })
  
  # Stop the progress bar
  close(pb)
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Small printout to the user
  cat(sprintf('Parallel time = %7.2f\n', parallel_timer[3]))
  
  # Save the results to the stumps_results list
  stumps_results[[o]]$stumps_model_deviance_test_smooth_mat = 
    t(sapply(1:length(stumps_lambdas), function(x) parallel_list[[x]]$deviance_smooth))
  
  stumps_results[[o]]$stumps_model_deviance_test_step_mat = 
    t(sapply(1:length(stumps_lambdas), function(x) parallel_list[[x]]$deviance_step))
}

saveRDS(list('stumps_results' = stumps_results,
             'stumps_lambdas' = stumps_lambdas),
        file = paste('stumps_outer_', outer_repetitions, '_numiter_',
                     num_iter$smooth, '_', num_iter$step, '_smooth_step_results', sep = ''))

temp_list = readRDS(paste('stumps_outer_', outer_repetitions, '_numiter_',
                          num_iter$smooth, '_', num_iter$step, '_smooth_step_results', sep = ''))
stumps_results = temp_list$stumps_results


stumps_best_dev_smooth = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_smooth_mat, 1, min))
stumps_best_dev_step   = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, min))

stumps_best_iters_smooth = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_smooth_mat, 1, which.min)) - 1 
stumps_best_iters_step   = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, which.min)) - 1

stumps_avg_dev_smooth = apply(stumps_best_dev_smooth, 1, mean)
stumps_avg_dev_step   = apply(stumps_best_dev_step, 1, mean)

stumps_avg_iter_smooth = apply(stumps_best_iters_smooth, 1, mean)
stumps_avg_iter_step   = apply(stumps_best_iters_step, 1, mean)


# ### Check that we got the same values for the coincing penaltites. YES it worked
# temp = readRDS(paste('stumps_outer_', outer_repetitions, '_numiter_',
#                      num_iter$smooth, '_', num_iter$step, '_smooth_step_results', sep = ''))
# #stumps_results = temp$stumps_results
# 
# stumps_best_dev_step1   = sapply(1:outer_repetitions, function(x) apply(temp$stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, min))
# stumps_avg_dev_step1   = apply(stumps_best_dev_step1, 1, mean)
# matrix(c(temp$stumps_lambdas, stumps_avg_dev_step1), nrow = 2, byrow = TRUE)
# 
# temp2 = readRDS('stumps_only_step_outer_10_numiter_60_smooth_step_results')
# stumps_results = temp2$stumps_results
# stumps_best_dev_step2   = sapply(1:outer_repetitions, function(x) apply(temp2$stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, min))
# stumps_avg_dev_step2  = apply(stumps_best_dev_step2, 1, mean)
# matrix(c(temp2$stumps_lambdas, stumps_avg_dev_step2), nrow = 2, byrow = TRUE)
# 


png(paste('stumps_outer_', outer_repetitions, '_numiter_',
          num_iter$smooth, '_', num_iter$step, '_smooth_step_results.png', sep = ""), width = 3500, height = 2250, res = 350)
par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
{
  ### Smooth
  # Plot 1
  plot(stumps_lambdas, stumps_avg_dev_smooth, type = 'l',
       xlab = '', ylab = '', ylim = c(min(stumps_avg_dev_smooth), 2300))
  title(main = bquote('Smooth Stumps'), cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2.1, cex.lab = 1.35)
  text(200, 2294, bquote(lambda~'='~.(stumps_lambdas[which.min(stumps_avg_dev_smooth)])
                         ~' min ='~.(round(min(stumps_avg_dev_smooth)))), cex = 1.35)
  
  # Plot 2
  plot(stumps_lambdas,
       stumps_avg_iter_smooth,
       ylim = c(0, num_iter$smooth), type = 'l', xlab = '', ylab = '')
  title(main = bquote('Smooth Stumps'), cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2.1, cex.lab = 1.35)
  
  ### Step
  # Plot 3
  plot(stumps_lambdas, stumps_avg_dev_step, type = 'l',
       xlab = '', ylab = '', ylim = c(min(stumps_avg_dev_step), 1460))
  title(main = bquote('Stepwise Stumps'), cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2.1, cex.lab = 1.35)
  text(200, 1458, bquote(lambda~'='~.(stumps_lambdas[which.min(stumps_avg_dev_step)])
                         ~' min ='~.(round(min(stumps_avg_dev_step)))), cex = 1.35)
  
  # Plot 4
  plot(stumps_lambdas,
       stumps_avg_iter_step,
       ylim = c(0, num_iter$step), type = 'l', xlab = '', ylab = '')
  title(main = bquote('Stepwise Stumps'), cex.main = 1.65, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2.1, cex.lab = 1.35)
}
dev.off()










####################################################################################
##### Look at best penalty for gaussian stumps on step data with small lambdas #####
####################################################################################
stumps_lambdas = seq(0.1, 10, 0.1)
num_iter = 60
num_iter = list('smooth' = 1500, 'step' = 60)
seed_train = 2020
seed_test = 2030

outer_repetitions = 10
stumps_results = list()

for (o in 1:outer_repetitions) {
  cat(sprintf('Outer iteration: %4d\n', o))
  
  stumps_results[[o]] = list()
  #stumps_model_deviance_test_smooth_mat = matrix(NA, ncol = num_iter$smooth+1, nrow = length(stumps_lambdas))
  #stumps_model_deviance_test_step_mat   = matrix(NA, ncol = num_iter$step+1, nrow = length(stumps_lambdas))
  
  # Generate the new data
  {
    # training_data_smooth = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, cn = cn,
    #                                             type = 'smooth', family = family)
    training_data_step = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, cn = cn,
                                              type = 'step', family = family)
    
    # testing_data_smooth  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = cn,
    #                                             type = 'smooth', family = family)
    testing_data_step  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = cn,
                                              type = 'step', family = family)  
  }
  
  ## Parallel fitting of model with different lambdas
  # Find the number of available cores
  cores=detectCores() 
  
  # Set up the clusters
  cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
  registerDoSNOW(cl)
  
  # Create a progress bar
  pb = txtProgressBar(max = length(stumps_lambdas), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  parallel_timer = system.time({
    # Iterate over the penalty factors
    parallel_list = foreach(i = 1:length(stumps_lambdas), .export = 'GAMBoost', .options.snow = opts) %dopar% {
      # Get teh current lambda 
      cur_lambda = stumps_lambdas[i]
      #cat(sprintf('Lambda: %5.1f   ', cur_lambda))
      
      ### Fit the stumps
      # stumps_model_smooth = GAMBoost_stumps(X           = training_data_smooth$X, 
      #                                       y           = training_data_smooth$y,
      #                                       num_iter    = num_iter$smooth,
      #                                       lambda      = cur_lambda,
      #                                       family      = family,
      #                                       print_msg   = FALSE,
      #                                       tiny_return = TRUE)
      
      stumps_model_step = GAMBoost_stumps(X           = training_data_step$X, 
                                          y           = training_data_step$y,
                                          num_iter    = num_iter$step,
                                          lambda      = cur_lambda,
                                          family      = family,
                                          print_msg   = FALSE,
                                          tiny_return = TRUE)
      
      
      # Arrays to store data
      # stumps_model_deviance_test_smooth_arr = rep(NA, num_iter$smooth+1)
      stumps_model_deviance_test_step_arr   = rep(NA, num_iter$step+1)
      
      
      # Iterate over all boosting iterations and find the best iteration based on a test set.
      # for (iter in 1:(num_iter$smooth+1)) {
      #   ## Smooth
      #   # Get the predicted results
      #   predicted_response = GAMBoost_stumps_predict(stumps_model_smooth,
      #                                                X_new = testing_data_smooth$X,
      #                                                after_iter = (iter-1),
      #                                                version = 'response')
      #   
      #   # Compute the deviance
      #   stumps_model_deviance_test_smooth_arr[iter] = 
      #     sum(family$dev.resids(testing_data_smooth$y,
      #                           predicted_response,
      #                           rep(1, n_test)))
      # }
      
      for (iter in 1:(num_iter$step+1)) {
        ## Step
        # Get the predicted results
        predicted_response = GAMBoost_stumps_predict(stumps_model_step,
                                                     X_new = testing_data_step$X,
                                                     after_iter = (iter-1),
                                                     version = 'response')
        
        # Compute the deviance
        stumps_model_deviance_test_step_arr[iter] =
          sum(family$dev.resids(testing_data_step$y,
                                predicted_response,
                                rep(1, n_test)))
      }
      
      # List of objects to be returned
      ret_list = list('lambda' = cur_lambda,
                      'deviance_step' = stumps_model_deviance_test_step_arr)
      #'deviance_smooth' = stumps_model_deviance_test_smooth_arr,
    }                          
  })
  
  # Stop the progress bar
  close(pb)
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Small printout to the user
  cat(sprintf('Parallel time = %7.2f\n', parallel_timer[3]))
  
  # # Save the results to the stumps_results list
  # stumps_results[[o]]$stumps_model_deviance_test_smooth_mat = 
  #   t(sapply(1:length(stumps_lambdas), function(x) parallel_list[[x]]$deviance_smooth))
  
  stumps_results[[o]]$stumps_model_deviance_test_step_mat = 
    t(sapply(1:length(stumps_lambdas), function(x) parallel_list[[x]]$deviance_step))
}

saveRDS(list('stumps_results' = stumps_results,
             'stumps_lambdas' = stumps_lambdas),
        file = paste('stumps_only_step_outer_', outer_repetitions, '_numiter_',
                     num_iter$step, '_smooth_step_results', sep = ''))

stumps_results = readRDS(paste('stumps_only_step_outer_', outer_repetitions, '_numiter_',
                               num_iter$step, '_smooth_step_results', sep = ''))$stumps_results

#stumps_best_dev_smooth = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_smooth_mat, 1, min))
stumps_best_dev_step   = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, min))

#stumps_best_iters_smooth = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_smooth_mat, 1, which.min)) - 1 
stumps_best_iters_step   = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, which.min)) - 1

#stumps_avg_dev_smooth = apply(stumps_best_dev_smooth, 1, mean)
stumps_avg_dev_step   = apply(stumps_best_dev_step, 1, mean)

#stumps_avg_iter_smooth = apply(stumps_best_iters_smooth, 1, mean)
stumps_avg_iter_step   = apply(stumps_best_iters_step, 1, mean)


png(paste('stumps_only_step_outer_', outer_repetitions, '_numiter_',
          num_iter$step, '_smooth_step_results2.png', sep = ""), width = 3500, height = 1250, res = 350)
par(mfrow=c(1,2), mar = c(3.5,3.5,2,1))
{
  # ### Smooth
  # # Plot 1
  # plot(stumps_lambdas, stumps_avg_dev_smooth, type = 'l', main = bquote('Smooth stumps'),
  #      xlab = '', ylab = '', ylim = c(min(stumps_avg_dev_smooth), 2300))
  # title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2)
  # text(200, 2300, bquote(lambda~'='~.(stumps_lambdas[which.min(stumps_avg_dev_smooth)])
  #                        ~' min ='~.(round(min(stumps_avg_dev_smooth)))))
  # 
  # # Plot 2
  # plot(stumps_lambdas,
  #      stumps_avg_iter_smooth,
  #      ylim = c(0, num_iter$smooth), type = 'l', main = bquote('Smooth stumps'), xlab = '', ylab = '')
  # title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2)
  
  ### Step
  # Plot 3
  plot(stumps_lambdas, stumps_avg_dev_step, type = 'l',
       xlab = '', ylab = '')#, ylim = c(min(stumps_avg_dev_step), 1400))
  title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2.1, cex.lab = 1.15)
  title( main = bquote('Stepwise Stumps'), cex.main = 1.45, line = 0.5)
  text(5, 1474, bquote(lambda~'='~.(stumps_lambdas[which.min(stumps_avg_dev_step)])
                       ~' min ='~.(round(min(stumps_avg_dev_step)))), cex = 1.15)
  
  # Plot 4
  plot(stumps_lambdas,
       stumps_avg_iter_step,
       ylim = c(0, 10), type = 'l', xlab = '', ylab = '')
  title( main = bquote('Stepwise Stumps'), cex.main = 1.45, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2.1, cex.lab = 1.15)
}
dev.off()







#################################################################################################
##### Look at best penalty for gaussian stumps on advanced stepwise data with small lambdas #####
#################################################################################################
stumps_lambdas = seq(0.1, 10, 0.1)
stumps_lambdas = c(0.1, 1, 2, 5, 10, 15, 20, 25, 30, 40, seq(50, 400, 25))
num_iter = 500
num_iter = list('smooth' = 1500, 'step' = 500)
seed_train = 2020
seed_test = 2030
stnr = 3

outer_repetitions = 10
stumps_results = list()

for (o in 1:outer_repetitions) {
  cat(sprintf('Outer iteration: %4d\n', o))
  
  stumps_results[[o]] = list()
  #stumps_model_deviance_test_smooth_mat = matrix(NA, ncol = num_iter$smooth+1, nrow = length(stumps_lambdas))
  #stumps_model_deviance_test_step_mat   = matrix(NA, ncol = num_iter$step+1, nrow = length(stumps_lambdas))
  
  # Generate the new data
  {
    # training_data_smooth = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, cn = cn,
    #                                             type = 'smooth', family = family)
    training_data_step = create_data_combined(n_train, seed_number = seed_train + (o-1)*7, stnr = stnr,
                                              type = 'advancedstep', family = family)
    
    # testing_data_smooth  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = cn,
    #                                             type = 'smooth', family = family)
    testing_data_step  = create_data_combined(n_test, seed_number = seed_test + (o-1)*7, cn = training_data_step$cn,
                                              type = 'advancedstep', family = family)  
  }
  
  ## Parallel fitting of model with different lambdas
  # Find the number of available cores
  cores=detectCores() 
  
  # Set up the clusters
  cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
  registerDoSNOW(cl)
  
  # Create a progress bar
  pb = txtProgressBar(max = length(stumps_lambdas), style = 3)
  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  
  parallel_timer = system.time({
    # Iterate over the penalty factors
    parallel_list = foreach(i = 1:length(stumps_lambdas), .export = 'GAMBoost', .options.snow = opts) %dopar% {
      # Get teh current lambda 
      cur_lambda = stumps_lambdas[i]
      #cat(sprintf('Lambda: %5.1f   ', cur_lambda))
      
      ### Fit the stumps
      # stumps_model_smooth = GAMBoost_stumps(X           = training_data_smooth$X, 
      #                                       y           = training_data_smooth$y,
      #                                       num_iter    = num_iter$smooth,
      #                                       lambda      = cur_lambda,
      #                                       family      = family,
      #                                       print_msg   = FALSE,
      #                                       tiny_return = TRUE)
      
      stumps_model_step = GAMBoost_stumps(X           = training_data_step$X, 
                                          y           = training_data_step$y,
                                          num_iter    = num_iter$step,
                                          lambda      = cur_lambda,
                                          family      = family,
                                          print_msg   = FALSE,
                                          tiny_return = TRUE)
      
      
      # Arrays to store data
      # stumps_model_deviance_test_smooth_arr = rep(NA, num_iter$smooth+1)
      stumps_model_deviance_test_step_arr   = rep(NA, num_iter$step+1)
      
      
      # Iterate over all boosting iterations and find the best iteration based on a test set.
      # for (iter in 1:(num_iter$smooth+1)) {
      #   ## Smooth
      #   # Get the predicted results
      #   predicted_response = GAMBoost_stumps_predict(stumps_model_smooth,
      #                                                X_new = testing_data_smooth$X,
      #                                                after_iter = (iter-1),
      #                                                version = 'response')
      #   
      #   # Compute the deviance
      #   stumps_model_deviance_test_smooth_arr[iter] = 
      #     sum(family$dev.resids(testing_data_smooth$y,
      #                           predicted_response,
      #                           rep(1, n_test)))
      # }
      
      for (iter in 1:(num_iter$step+1)) {
        ## Step
        # Get the predicted results
        predicted_response = GAMBoost_stumps_predict(stumps_model_step,
                                                     X_new = testing_data_step$X,
                                                     after_iter = (iter-1),
                                                     version = 'response')
        
        # Compute the deviance
        stumps_model_deviance_test_step_arr[iter] =
          sum(family$dev.resids(testing_data_step$y,
                                predicted_response,
                                rep(1, n_test)))
      }
      
      # List of objects to be returned
      ret_list = list('lambda' = cur_lambda,
                      'deviance_step' = stumps_model_deviance_test_step_arr)
      #'deviance_smooth' = stumps_model_deviance_test_smooth_arr,
    }                          
  })
  
  # Stop the progress bar
  close(pb)
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Small printout to the user
  cat(sprintf('Parallel time = %7.2f\n', parallel_timer[3]))
  
  # # Save the results to the stumps_results list
  # stumps_results[[o]]$stumps_model_deviance_test_smooth_mat = 
  #   t(sapply(1:length(stumps_lambdas), function(x) parallel_list[[x]]$deviance_smooth))
  
  stumps_results[[o]]$stumps_model_deviance_test_step_mat = 
    t(sapply(1:length(stumps_lambdas), function(x) parallel_list[[x]]$deviance_step))
}

saveRDS(list('stumps_results' = stumps_results,
             'stumps_lambdas' = stumps_lambdas),
        file = paste('stumps_only_step_outer_', outer_repetitions, '_numiter_',
                     num_iter$step, '_advanced_step_results', sep = ''))

stumps_results = readRDS(paste('stumps_only_step_outer_', outer_repetitions, '_numiter_',
                               num_iter$step, '_advanced_step_results', sep = ''))$stumps_results

#stumps_best_dev_smooth = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_smooth_mat, 1, min))
stumps_best_dev_step   = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, min))

#stumps_best_iters_smooth = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_smooth_mat, 1, which.min)) - 1 
stumps_best_iters_step   = sapply(1:outer_repetitions, function(x) apply(stumps_results[[x]]$stumps_model_deviance_test_step_mat, 1, which.min)) - 1

#stumps_avg_dev_smooth = apply(stumps_best_dev_smooth, 1, mean)
stumps_avg_dev_step   = apply(stumps_best_dev_step, 1, mean)

#stumps_avg_iter_smooth = apply(stumps_best_iters_smooth, 1, mean)
stumps_avg_iter_step   = apply(stumps_best_iters_step, 1, mean)


png(paste('stumps_only_step_outer_', outer_repetitions, '_numiter_',
          num_iter$step, '_advanced_step_results2.png', sep = ""), width = 3500, height = 1250, res = 350)
par(mfrow=c(1,2), mar = c(3.5,3.5,2,1))
{
  # ### Smooth
  # # Plot 1
  # plot(stumps_lambdas, stumps_avg_dev_smooth, type = 'l', main = bquote('Smooth stumps'),
  #      xlab = '', ylab = '', ylim = c(min(stumps_avg_dev_smooth), 2300))
  # title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2)
  # text(200, 2300, bquote(lambda~'='~.(stumps_lambdas[which.min(stumps_avg_dev_smooth)])
  #                        ~' min ='~.(round(min(stumps_avg_dev_smooth)))))
  # 
  # # Plot 2
  # plot(stumps_lambdas,
  #      stumps_avg_iter_smooth,
  #      ylim = c(0, num_iter$smooth), type = 'l', main = bquote('Smooth stumps'), xlab = '', ylab = '')
  # title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2)
  
  ### Step
  # Plot 3
  plot(stumps_lambdas, stumps_avg_dev_step, type = 'l',
       xlab = '', ylab = '')#, ylim = c(min(stumps_avg_dev_step), 1400))
  title(xlab = 'Penalty', ylab = 'Average Test Deviance', line = 2.1, cex.lab = 1.15)
  title( main = bquote('Advanced Stepwise Stumps'), cex.main = 1.45, line = 0.5)
  text(5, 1474, bquote(lambda~'='~.(stumps_lambdas[which.min(stumps_avg_dev_step)])
                       ~' min ='~.(round(min(stumps_avg_dev_step)))), cex = 1.15)
  
  # Plot 4
  plot(stumps_lambdas,
       stumps_avg_iter_step,
       ylim = c(0, 10), type = 'l', xlab = '', ylab = '')
  title( main = bquote('Advanced Stepwise Stumps'), cex.main = 1.45, line = 0.5)
  title(xlab = 'Penalty', ylab = 'Average Best Iteration', line = 2.1, cex.lab = 1.15)
}
dev.off()













# 
# 
# # RUn the first lines in this file!
# lambdas = c(c(0.1,0.5,1,2,3,4,5,6,7,8,9), seq(10,100, length.out = 19), seq(120, 1000, length.out = 45))
# 
# lambdas = c(0.5,1,2,3,4,5,6,7,8,9,10)
# num_iter = 1500
# stumps_model_deviance_test_smooth_mat = matrix(NA, ncol = num_iter+1, nrow = length(lambdas))
# stumps_model_deviance_test_step_mat   = matrix(NA, ncol = num_iter+1, nrow = length(lambdas))
# stumps_model_smooth_aic = matrix(NA, ncol = num_iter+1, nrow = length(lambdas))
# stumps_model_smooth_bic = matrix(NA, ncol = num_iter+1, nrow = length(lambdas))
# stumps_model_step_aic = matrix(NA, ncol = num_iter+1, nrow = length(lambdas))
# stumps_model_step_bic = matrix(NA, ncol = num_iter+1, nrow = length(lambdas))
# 
# 
# for (i in 1:length(lambdas)) {
#   lambda = lambdas[i]
#   
#   ### Fit the stumps
#   stumps_model_smooth = GAMBoost_stumps(X         = training_data_smooth$X, 
#                                         y         = training_data_smooth$y,
#                                         num_iter  = num_iter,
#                                         lambda    = lambda,
#                                         family    = family, 
#                                         print_msg = TRUE)
#   
#   stumps_model_step = GAMBoost_stumps(X         = training_data_step$X, 
#                                       y         = training_data_step$y,
#                                       num_iter  = num_iter,
#                                       lambda    = lambda,
#                                       family    = family, 
#                                       print_msg = TRUE)
#   
#   for (iter in 1:(num_iter+1)) {
#     cat(sprintf("Iteration %d of %d\n", iter, num_iter+1))
    # ## Smooth
    # # Get the predicted results
    # predicted_response = GAMBoost_stumps_predict(stumps_model_smooth, X_new = testing_data_smooth$X,
    #                                              after_iter = (iter-1), version = 'response')
    # 
    # # Compute the deviance
    # stumps_model_deviance_test_smooth_mat[i, iter] = sum(family$dev.resids(testing_data_smooth$y,
    #                                                                        predicted_response,
    #                                                                        rep(1, n_test)))
    # 
    # ## Step
    # # Get the predicted results
    # predicted_response = GAMBoost_stumps_predict(stumps_model_step, X_new = testing_data_step$X,
    #                                              after_iter = (iter-1), version = 'response')
    # 
    # # Compute the deviance
    # stumps_model_deviance_test_step_mat[i, iter] = sum(family$dev.resids(testing_data_step$y,
    #                                                                      predicted_response,
    #                                                                      rep(1, n_test)))
#   }
#   stumps_model_smooth_aic[i, ] = stumps_model_smooth$AIC
#   stumps_model_smooth_bic[i, ] = stumps_model_smooth$BIC
#   stumps_model_step_aic[i, ]   = stumps_model_step$AIC
#   stumps_model_step_bic[i, ]   = stumps_model_step$BIC
# }
# 
# saveRDS(list('stumps_model_step_aic' = stumps_model_step_aic,
#              'stumps_model_step_bic' = stumps_model_step_bic,
#              'stumps_model_smooth_aic' = stumps_model_smooth_aic,
#              'stumps_model_smooth_bic' = stumps_model_smooth_aic,
#              'stumps_model_deviance_test_step_mat' = stumps_model_deviance_test_step_mat,
#              'stumps_model_deviance_test_smooth_mat' = 'stumps_model_deviance_test_smooth_mat',
#              'lambda' = lambda),
#         file = 'Best_Lambda_1500_iters_smoth_step_results')
# 
# 
# lambda_smooth_min = lambdas[which.min(stumps_model_deviance_test_smooth_mat[cbind(1:length(lambdas), best_iters_test_smooth)])]
# lambda_step_min = lambdas[which.min(stumps_model_deviance_test_step_mat[cbind(1:length(lambdas), best_iters_test_step)])]
# smooth_min = round(min(stumps_model_deviance_test_smooth_mat[cbind(1:length(lambdas), best_iters_test_smooth)]))
# step_min = round(min(stumps_model_deviance_test_step_mat[cbind(1:length(lambdas), best_iters_test_step)]))
# 
# png(paste("stumps_min_test_deviance_different_penalty.png", sep = ""), width = 3500, height = 2250, res = 350)
# par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
# {
#   best_iters_test_smooth = apply(stumps_model_deviance_test_smooth_mat, 1, which.min)
#   best_iters_test_step   = apply(stumps_model_deviance_test_step_mat, 1, which.min)
#   test_dev_smooth = stumps_model_deviance_test_smooth_mat[cbind(1:length(lambdas), best_iters_test_smooth)]
#   test_dev_step = stumps_model_deviance_test_step_mat[cbind(1:length(lambdas), best_iters_test_step)]
#   
#   plot(lambdas[1:min(which(best_iters_test_smooth >= 1500))],
#        test_dev_smooth[1:min(which(best_iters_test_smooth >= 1500))],
#        type = 'l', main = bquote('Smooth'),
#        xlab = '', ylab = '', xlim = c(0,1000))
#   lines(lambdas[min(which(best_iters_test_smooth >= 1500)):length(lambdas)],
#         test_dev_smooth[min(which(best_iters_test_smooth >= 1500)):length(lambdas)],
#         col = 'darkgray', lty = 2)
#   title(xlab = 'Penalty', ylab = 'Test Deviance', line = 2)
#   #abline(v = lambdas[min(which(best_iters_test_smooth >= 1500))], lty = 2)
#   text(500, 2345, bquote(lambda~'='~.(lambda_smooth_min)~' min ='~.(smooth_min)))
#   
#   plot(lambdas[1:min(which(best_iters_test_smooth >= 1500))],
#        best_iters_test_smooth[1:min(which(best_iters_test_smooth >= 1500))], xlim = c(0,1000),
#        ylim = c(0, 1500), type = 'l', main = bquote('Smooth'), xlab = '', ylab = '')
#   lines(lambdas[min(which(best_iters_test_smooth >= 1500)):length(lambdas)],
#         best_iters_test_smooth[min(which(best_iters_test_smooth >= 1500)):length(lambdas)],
#         col = 'darkgray', lty = 2)
#   title(xlab = 'Penalty', ylab = 'Best Iteration', line = 2)
#   #abline(v = lambdas[min(which(best_iters_test_smooth >= 1500))], lty = 2)
#   
#   
#   
#   plot(lambdas, test_dev_step,
#        main = bquote('Stepwise'),
#        xlab = '', ylab = '', pch = 19, type = 'l')
#   title(xlab = 'Penalty', ylab = 'Test Deviance', line = 2)
#   #abline(v = 1000, lty = 2)
#   text(500,1406,bquote(lambda~'='~.(lambda_step_min)~' min ='~.(step_min)))
#   
#   plot(lambdas, best_iters_test_step, ylim = c(0, 1500),
#        type = 'l', main = bquote('Stepwise'), xlab = '', ylab = '', pch = 19)
#   title(xlab = 'Penalty', ylab = 'Best Iteration', line = 2)
#   #abline(v = 1000, lty = 2)
# }
# dev.off()
# 
# 
# which.min(stumps_model_deviance_test_smooth_mat[cbind(1:length(lambdas), best_iters_test_smooth)])
# which.min(stumps_model_deviance_test_step_mat[cbind(1:length(lambdas), best_iters_test_step)])
# 
# 
# plot(lambdas, best_iters_test_step,
#      main = 'Step', xlab = '', ylab = '', pch = 19, type = 'l')
# 
# lambdas[21]
# lambdas[6]
# 
# 
# min(which(best_iters_test_smooth >= 240))
