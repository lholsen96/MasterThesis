##### INSTURCTIONS TO LARS. TO DO LIST
# 0.    Create the data.
# 1.    Fit model by running it for ish 300 iterations
# 2.    Determine the best iteration based on test deviance
# 2.5.  Can also look at best iteration based on AIC, BIC
#       (CV if we implement it, but not needed!)
# 3.    Genereate 1000 bootstap samples
# 3.5.  Fit a model to each of these samples. Run for 
#       same number of iterations (~300). 
# 4.    Plot the development of the fitted functions
#       iter 3, 5, 10, 25, 50, 100, 250 etc.
# 5.    For the best iteration we want to do a coverage
#       analysis of the approximate pointwise confidence band
#       from the hat matrix and the empirical pointwise 
#       confidence bands from the bootstrap samples. 
#       We do this by counting the number of times the true
#       function line is inside the confidence bands.
#       We evalute this at the n_train original observations
#       in the training data.
# 6.    Repeat 1-5 for Gaussian, Binomial and Poisson
# 7.    Repeat 1-6 for smooth f and step f. 
# 8.    Repeat 1-7 for splines and stumps
# 9.    Maybe alter the cn for different signal to noise rati
# 10.   Maybe look at different penalizations. now lambda stumps = 2 and splines = 30

# rm(list=ls())
setwd("C:/Users/lars9/OneDrive/Master/R_codes/GAMBoost_final")
setwd("~/Master_R_files")
setwd("D:/Rmaster")

#
nice_colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)


# Include libraries
library(GAMBoost)
library(parallel)
library(doParallel)
library(foreach)
library(doSNOW)

# Define possible distributions and function types
families = c('gaussian', 'binomial', 'poisson')
types = c('smooth', 'step')

# Parameters
n_train = 100
n_test  = 1000

# NUmber of bootstrap samples for each model
B = list('splines' = 999, 'stumps' = 999)

num_iter = list('splines' = 250, 'stumps' = 250)
penalty  = list('splines' = 300, 'stumps' = 5)

# Set the signal to noise ratio
stnr = c(1, 3, 10)
# DONE stnr = 3
# DONE
stnr = 10
# DONE stnr = 1

# Controls signal to noise. Higher = more signal.
print_msg = FALSE

set.seed(2020)
list_of_param = list()
for (type in types) {
  for (family in families) {
    temp_seeds = sample(100000, 2)
    list_of_param[[paste(type, family, sep = '')]] = list('seed_train' = temp_seeds[1],
                                                          'seed_test'  = temp_seeds[2])
  }
}
list_of_param


for (type in types) {
  for (family in families) {
    
    # Change from string to fucntion
    family = switch(family, gaussian = gaussian(), binomial = binomial(), poisson = poisson())
    
    # Get the paramaters
    typefamily = paste(type, family$family, sep = '')
    params = list_of_param[[typefamily]]
    print(typefamily)
    
    ### Setup: generate data
    training_data = create_data_combined(n_train, seed_number = params$seed_train, stnr = stnr,
                                         type = type, family = family)
    
    testing_data  = create_data_combined(n_test, seed_number = params$seed_test, cn = training_data$cn,
                                         type = type, family = family)
    
    list_of_param[[typefamily]][['training_data']] = training_data
    list_of_param[[typefamily]][['testing_data']]  = testing_data
    
    # training_data$cn
    # sum((training_data$mu - mean(training_data$mu))^2) / sum(family$variance(training_data$mu))
    # sum((testing_data$mu - mean(testing_data$mu))^2) / sum(family$variance(testing_data$mu))
    # 
    
    #################################
    ##### GAMBoost with splines #####
    #################################
    # temp = optimGAMBoostPenalty(x        = training_data$X,
    #                      y        = training_data$y,
    #                      minstepno=50, maxstepno=200,
    #                      penalty = 10, iter.max = 10,
    #                      start.penalty=100,
    #                      family   = family)
    # temp$penalty 
    # STOPS when penalty yields
    # is larger or equal to 50. With a smaller number of steps boosting may
    # become too “greedy” and show sub-optimal performance. 
    single_timer = system.time({
      splines_model = GAMBoost(x        = training_data$X,
                               y        = training_data$y,
                               stepno   = num_iter[['splines']],
                               penalty  = penalty[['splines']],
                               family   = family,
                               calc.hat = TRUE,
                               calc.se  = TRUE)
      #,AIC.type = 'classical')
    })
    
    list_of_param[[typefamily]][['splines_time_single']]  = single_timer
    
    # Save the model
    saveRDS(splines_model, 
            file = paste("splines_", type, "_", family$family, "_stnr_", stnr,
                         '_lambda_', penalty[['splines']], sep = ""))
    
    
    
    
    ##### Test Deviance Development #####
    splines_model_deviance_test = rep(NA, num_iter[['splines']]+1)
    
    for (iter in 1:(num_iter[['splines']]+1)) {
      splines_model_deviance_test[iter] = sum(family$dev.resids(
        testing_data$y, 
        predict(splines_model, newdata = testing_data$X, at.step = (iter-1), type = "response"),
        rep(1, n_test)))
    }
    # Find the optimal test deviance, AIC and BIC
    splines_model_optimal_test = which.min(splines_model_deviance_test) - 1
    splines_model_optimal_AIC  = which.min(splines_model$AIC) - 1
    splines_model_optimal_BIC  = which.min(splines_model$BIC) - 1
    
    list_of_param[[typefamily]][['splines_test']] = splines_model_optimal_test
    list_of_param[[typefamily]][['splines_AIC']]  = splines_model_optimal_AIC
    list_of_param[[typefamily]][['splines_BIC']]  = splines_model_optimal_BIC
    
    
    ### Crate a figure of the test deviance error alongside the scaled AIC and BIC
    # Scaled AIC and BIC such that the initial value is the same as test deviance.
    scaled_AIC = splines_model$AIC * splines_model_deviance_test[1]/splines_model$AIC[1]
    scaled_BIC = splines_model$BIC * splines_model_deviance_test[1]/splines_model$BIC[1]
    
    miny = min(splines_model_deviance_test, scaled_AIC, scaled_BIC)
    maxy = max(splines_model_deviance_test, scaled_AIC, scaled_BIC)
    range = maxy - miny
    
    save_filename = paste("splines_", type, "_", family$family, "_stnr_", stnr,
                          '_lambda_', penalty[['splines']], "_test_deviance_error", sep = "")
    
    
    png(paste(save_filename, ".png", sep = ""), width = 3500, height = 1500, res = 350)
    {
      par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
      plot(0:num_iter[['splines']], splines_model_deviance_test, ylab = "", xlab = "",
           type = "l", col = nice_colors[1], ylim = c(miny-0.1*range, maxy+0.1*range))
      title(xlab = 'Iteration', ylab = 'Test devaince / scaled AIC and BIC', line = 2.25)
      title(main = paste('Splines:', family$family, type))
      
      lines(0:num_iter[['splines']], scaled_AIC, lty = 1, col = nice_colors[2])
      lines(0:num_iter[['splines']], scaled_BIC, lty = 1, col = nice_colors[3])
      
      abline(v = splines_model_optimal_test, lty = 2, col = nice_colors[1])
      abline(v = splines_model_optimal_AIC,  lty = 3, col = nice_colors[2])
      abline(v = splines_model_optimal_BIC,  lty = 4, col = nice_colors[3])
      
      legend('topright', legend = c(paste('Test(', splines_model_optimal_test, ')', sep = ''),
                                    paste('AIC(', splines_model_optimal_AIC, ')', sep = ''),
                                    paste('BIC(', splines_model_optimal_BIC, ')', sep = '')),
             lty = 1, col = nice_colors, cex = 0.75)
    }    
    dev.off()
    
    
    
    ##### Create the booststrap empirical confidence intervals #####
    # # List to store all the models
    # splines_boot_models = list()
    # 
    # # Iterate over the B bootstrap samples
    # tic("seq")
    # for (b in 1:B[['splines']]) {
    #   cat(sprintf('Splines %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
    #   
    #   # The indices for this bootstrap sample
    #   boot_idx = sample(1:n_train, n_train, replace = TRUE)
    #   
    #   # Get the relevant data for b
    #   boot_X_train = training_data$X[boot_idx, ]
    #   boot_y_train = training_data$y[boot_idx]
    #   
    #   # Fit the model (100 iterations = 1.5 Mb)
    #   boot_pspline = GAMBoost(x = boot_X_train, 
    #                           y = boot_y_train,
    #                           stepno   = num_iter[['splines']],
    #                           penalty  = penalty[['splines']],
    #                           family   = family,
    #                           calc.hat = FALSE, # Do not need these as we are only interested in the 
    #                           calc.se  = FALSE) # the predicted curves.
    #   
    #   splines_boot_models[[b]] = boot_pspline
    # }
    # seq1 = toc()
    
    ### Parallel bootstraps
    #Find the number of available cores
    cores=detectCores() 
    list_of_param[[typefamily]][['splines_cores']] = cores[1] - 1
    
    # Set up the clusters
    cl = parallel::makeCluster(cores[1]-1) #not to overload your computer
    registerDoSNOW(cl)
    
    # Iterate over the B[['splines']] bootstrap samples 
    parallel_timer = system.time({
      splines_boot_models = foreach(b = 1:B[['splines']], .export = 'GAMBoost') %dopar% {
        #cat(sprintf('Splines %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
        
        # The indices for this bootstrap sample
        boot_idx = sample(1:n_train, n_train, replace = TRUE)
        
        # Get the relevant data for b
        boot_X_train = training_data$X[boot_idx, ]
        boot_y_train = training_data$y[boot_idx]
        
        # Fit the model (100 iterations = 1.5 Mb)
        boot_pspline = GAMBoost(x = boot_X_train, 
                                y = boot_y_train,
                                stepno   = num_iter[['splines']],
                                penalty  = penalty[['splines']],
                                family   = family,
                                calc.hat = FALSE, # Do not need these as we are only interested in the 
                                calc.se  = FALSE) # the predicted curves.
      }})
    
    # Stop the cluster
    parallel::stopCluster(cl)
    
    # Stop the timer
    list_of_param[[typefamily]][['splines_time']]  = parallel_timer
    
    # Save bootstrap models
    saveRDS(splines_boot_models, 
            file = paste("splines_", type, "_", family$family, "_stnr_", stnr, 
                         '_lambda_', penalty[['splines']], "_bootstrap_", B[['splines']], sep = ""))
    
    
    
    ### Plot the fitted cureves alongside with approximate and empirical confidence bands
    after_iter = round(seq(1, num_iter[['splines']], length.out = 5))
    file_name = paste(paste("splines_", type, "_", family$family, "_stnr_", stnr,
                            '_lambda_', penalty[['splines']], '_boot_', B[['splines']],
                            '_iter_', sep = ""),
                      paste(after_iter, collapse = '_'), sep = '')
    
    GAMBoost_plot_predictor_contribution(splines_model,
                                         type = type, 
                                         boot_models = splines_boot_models,
                                         cn = training_data$cn,
                                         after_iter = after_iter,
                                         save_filename = file_name,
                                         xmin = -1,
                                         xmax = 1,
                                         ymin = NULL,
                                         ymax = NULL,
                                         lwd = 1) 
    
    # Only look at the best test iteration
    after_iter = splines_model_optimal_test
    file_name = paste(paste("splines_", type, "_", family$family, "_stnr_", stnr,
                            '_lambda_', penalty[['splines']], '_boot_', B[['splines']],
                            '_iter_', sep = ""),
                      paste(after_iter, collapse = '_'), sep = '')
    
    GAMBoost_plot_predictor_contribution(splines_model,
                                         type = type, 
                                         boot_models = splines_boot_models,
                                         cn = training_data$cn,
                                         after_iter = after_iter,
                                         save_filename = file_name,
                                         xmin = -1,
                                         xmax = 1,
                                         ymin = NULL,
                                         ymax = NULL,
                                         lwd = 1) 
    
    
    ### Look at the coverage for iteration splines_model_optimal_test
    png(paste("splines_", type, "_", family$family, "_stnr_", stnr,
              '_lambda_', penalty[['splines']], '_boot_', B[['splines']],
              '_iter_', splines_model_optimal_test, '_coverage.png',  sep = ""),
        width = 3500, height = 1500, res = 350)
    
    splines_coverage = GAMBoost_splines_coverage_covariate(splines_model, splines_boot_models,
                                                           after_iter = splines_model_optimal_test, 
                                                           type = type, cn = training_data$cn)
    
    dev.off()
    
    list_of_param[[typefamily]][['splines_coverage_best_test_iter']] = splines_coverage
    
    
    
    ### Look at development of coverage
    number = 50
    p = ncol(training_data$X)
    splines_approx_cov = matrix(NA, ncol = p, nrow = number)
    splines_empirical_cov = matrix(NA, ncol = p, nrow = number)
    iters = round(seq(1, num_iter[['splines']], length.out = number))
    for (ai in 1:number) {
      splines_coverage = GAMBoost_splines_coverage_covariate(splines_model, splines_boot_models,
                                                             after_iter = iters[ai], 
                                                             type = type, cn = training_data$cn)
      splines_approx_cov[ai, ] = splines_coverage$Approximate
      splines_empirical_cov[ai, ] = splines_coverage$Empirical
    }
    
    png(paste("splines_", type, "_", family$family, "_stnr_", stnr,
              '_lambda_', penalty[['splines']], '_boot_', B[['splines']], 
              '_coverage_development.png',  sep = ""),
        width = 3500, height = 1500, res = 350)
    
    par(mfrow=c(1,5), mar = c(3.5,3.5,2,1))
    for (s in 1:p) {
      plot(iters, splines_approx_cov[,s], col = nice_colors[1], type = 'l', ylim = c(0,1), xlab = '', ylab = '')
      title(xlab = 'Bosting iter.', ylab = 'Coverage', line = 2.25)
      
      lines(iters, splines_empirical_cov[,s], col = nice_colors[4])
      
      abline(h = 0.95, lty = 4)
    }
    dev.off()
    
    
    # ##### Coverage of the pointwise confidence intervals for the predicted y_i's
    # splines_model_optimal = GAMBoost(x        = training_data$X,
    #                                  y        = training_data$y,
    #                                  stepno   = splines_model_optimal_test,
    #                                  penalty  = penalty[['splines']],
    #                                  family   = family,
    #                                  calc.hat = TRUE,
    #                                  calc.se  = TRUE)
    #                                  #,AIC.type = 'classical')
    # 
    # # Save the model
    # saveRDS(splines_model_optimal, 
    #         file = paste("splines_", type, "_", family$family, "_stnr_", stnr,
    #                      '_lambda_', penalty[['splines']], '_iter_', splines_model_optimal_test,
    #                      sep = ""))
    # 
    # 
    # # create the CI
    # splines_pred_mu = family$linkinv(splines_model_optimal$eta[ ,splines_model_optimal_test+1])
    # splines_pointwise_var = diag(splines_model_optimal$hatmatrix %*%
    #                                diag(family$variance(splines_pred_mu)) %*%
    #                                t(splines_model_optimal$hatmatrix))
    # 
    # splines_lower = splines_pred_mu - 2*splines_pointwise_var
    # splines_upper = splines_pred_mu + 2*splines_pointwise_var
    # 
    # 
    # sum(splines_lower <= training_data$mu & training_data$mu <= splines_upper)
    # sum(splines_lower <= training_data$y & training_data$y <= splines_upper)
    # splines_order = order(splines_pred_mu)
    # 
    # plot(splines_lower[splines_order], type = 'l', lty = 2)
    # lines(splines_upper[splines_order], lty = 2)
    # lines(splines_pred_mu[splines_order], lty = 1)
    # points(training_data$mu[splines_order], col = 'steelblue')
    # points(training_data$y[splines_order], col = 'orange')
    # 
    # #plot(splines_pred_mu, training_data$y)
    # 
    # 
    # pred_temp = matrix(NA, ncol = 20, nrow = 100)
    # for(i in 1:20) {
    #   pred_temp[,i] = predict(splines_boot_models[[i]], newdata = training_data$X, type = 'response')
    # }
    # temp1 = apply(pred_temp, 1, quantile, c(0.025, 0.975))
    # 
    # splines_lower1 = temp1[1,]
    # splines_upper1 = temp1[2,]
    # sum(splines_lower1 <= training_data$mu & training_data$mu <= splines_upper1)
    # sum(splines_lower1 <= training_data$y  & training_data$y  <= splines_upper1)
    # 
    # plot(splines_lower[splines_order], type = 'l', lty = 2)
    # lines(splines_upper[splines_order], lty = 2)
    # lines(splines_lower1[splines_order], lty = 2, col = 'red')
    # lines(splines_upper1[splines_order], lty = 2, col = 'red')
    # lines(splines_pred_mu[splines_order], lty = 1)
    # points(training_data$mu[splines_order], col = 'steelblue')
    # points(training_data$y[splines_order], col = 'orange')
    
    
    
    
    
    
    ################################
    ##### GAMBoost with stumps #####
    ################################
    single_timer = system.time({
      stumps_model = GAMBoost_stumps(X         = training_data$X, 
                                     y         = training_data$y,
                                     num_iter  = num_iter[['stumps']],
                                     lambda    = penalty[['stumps']],
                                     family    = family, 
                                     print_msg = print_msg)
    })
    
    # Stop the timer
    list_of_param[[typefamily]][['stumps_time_single']]  = single_timer
    
    # Save the model
    saveRDS(stumps_model, 
            file = paste("stumps_", type, "_", family$family, "_stnr_", stnr,
                         '_lambda_', penalty[['stumps']], sep = ""))
    
    
    ##### Test Deviance Development #####
    stumps_model_deviance_test = rep(NA, num_iter[['stumps']]+1)
    
    for (iter in 1:(num_iter[['stumps']]+1)) {
      # Get the predicted results
      predicted_response = GAMBoost_stumps_predict(stumps_model, X_new = testing_data$X,
                                                   after_iter = (iter-1), version = 'response')
      
      # Compute the deviance
      stumps_model_deviance_test[iter] = sum(family$dev.resids(testing_data$y,
                                                               predicted_response,
                                                               rep(1, n_test)))
    }
    
    # Find the optimal number of iterations
    stumps_model_optimal_test = which.min(stumps_model_deviance_test) - 1
    stumps_model_optimal_AIC = which.min(stumps_model$AIC) - 1
    stumps_model_optimal_BIC = which.min(stumps_model$BIC) - 1
    
    list_of_param[[typefamily]][['stumps_test']] = stumps_model_optimal_test
    list_of_param[[typefamily]][['stumps_AIC']]  = stumps_model_optimal_AIC
    list_of_param[[typefamily]][['stumps_BIC']]  = stumps_model_optimal_BIC
    
    
    
    ### Crate a figure of the test deviance error
    # SCaled AIC and BIC such that the initial value is the same as test deviance.
    scaled_AIC = stumps_model$AIC * stumps_model_deviance_test[1]/stumps_model$AIC[1]
    scaled_BIC = stumps_model$BIC * stumps_model_deviance_test[1]/stumps_model$BIC[1]
    
    miny = min(stumps_model_deviance_test, scaled_AIC, scaled_BIC)
    maxy = max(stumps_model_deviance_test, scaled_AIC, scaled_BIC)
    range = maxy - miny
    
    save_filename = paste("stumps_", type, "_", family$family, "_stnr_", stnr,
                          '_lambda_', penalty[['stumps']], "_test_deviance_error", sep = "")
    
    png(paste(save_filename, ".png", sep = ""), width = 3500, height = 1500, res = 350)
    {
      par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
      plot(0:num_iter[['stumps']], stumps_model_deviance_test, ylab = "", xlab = "",
           type = "l", col = nice_colors[1], ylim = c(miny-0.1*range, maxy+0.1*range))
      title(xlab = 'Iteration', ylab = 'Test devaince / scaled AIC and BIC', line = 2.25)
      title(main = paste('Stumps:', family$family, type))
      
      lines(0:num_iter[['stumps']], scaled_AIC, lty = 1, col = nice_colors[2])
      lines(0:num_iter[['stumps']], scaled_BIC, lty = 1, col = nice_colors[3])
      
      abline(v = stumps_model_optimal_test, lty = 2, col = nice_colors[1])
      abline(v = stumps_model_optimal_AIC,  lty = 3, col = nice_colors[2])
      abline(v = stumps_model_optimal_BIC,  lty = 4, col = nice_colors[3])
      
      legend('topright', legend = c(paste('Test(', stumps_model_optimal_test, ')', sep = ''),
                                    paste('AIC(', stumps_model_optimal_AIC, ')', sep = ''),
                                    paste('BIC(', stumps_model_optimal_BIC, ')', sep = '')),
             lty = 1, col = nice_colors, cex = 0.75)
    }    
    dev.off()
    
    
    
    ##### Create the booststrap empirical confidence intervals #####
    # # List to store all the models
    # stumps_boot_models = list()
    # start_time_stumps_boot_models = Sys.time()
    # 
    # # Iterate over the B bootstrap samples
    # for (b in 1:B[['stumps']]) {
    #   cat(sprintf('stumps %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
    #   
    #   # The indices for this bootstrap sample
    #   boot_idx = sample(1:n_train, n_train, replace = TRUE)
    #   
    #   # Get the relevant data for b
    #   boot_X_train = training_data$X[boot_idx, ]
    #   boot_y_train = training_data$y[boot_idx]
    #   
    #   # Fit the model (250 iterations = 2MB)
    #   boot_pstump  = GAMBoost_stumps(X = boot_X_train, 
    #                                  y = boot_y_train,
    #                                  num_iter  = num_iter[['stumps']],
    #                                  lambda    = penalty[['stumps']],
    #                                  family    = family, 
    #                                  print_msg = FALSE,
    #                                  tiny_return = TRUE)
    #   
    #   stumps_boot_models[[b]] = boot_pstump
    #   
    #   end_time_stumps_boot_models_temp = Sys.time()
    #   
    #   (end_time_stumps_boot_models_temp - start_time_stumps_boot_models) / b
    # }
    # 
    # end_time_stumps_boot_models = Sys.time()
    # 
    # end_time_stumps_boot_models - start_time_stumps_boot_models
    # # Time difference of 8.134289 hours (Gaussian smooth)
    # 
    # # Save bootstrap models
    # saveRDS(stumps_boot_models, file = 
    #        paste("stumps_", type, "_", family$family, "_stnr_", stnr,
    #              '_lambda_', penalty[['stumps']], "_bootstrap_", B[['stumps']], sep = ""))
    # # stumps_boot_models = readRDS(paste("stumps_", type, "_", family$family, "_stnr_", stnr,
    # #   '_lambda_', penalty[['stumps']], "_bootstrap_", B[['stumps']], sep = ""))
    
    
    
    ### TRY to paralellize it
    ### Parallel bootstraps
    #Find the number of available cores
    cores=detectCores() 
    list_of_param[[typefamily]][['stumps_cores']] = cores[1] - 1
    
    # Set up the clusters
    cl = parallel::makeCluster(cores[1]-1) #not to overload your computer
    registerDoSNOW(cl)
    
    
    # Iterate over the B[['stumps']] bootstrap samples 
    parallel_timer_stumps = system.time({
      stumps_boot_models = foreach(b = 1:B[['stumps']]) %dopar% {
        #cat(sprintf('stumps %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
        
        # The indices for this bootstrap sample
        boot_idx = sample(1:n_train, n_train, replace = TRUE)
        
        # Get the relevant data for b
        boot_X_train = training_data$X[boot_idx, ]
        boot_y_train = training_data$y[boot_idx]
        
        # Fit the model (250 iterations = 2MB)
        boot_pstump  = GAMBoost_stumps(X = boot_X_train, 
                                       y = boot_y_train,
                                       num_iter  = num_iter[['stumps']],
                                       lambda    = penalty[['stumps']],
                                       family    = family, 
                                       print_msg = FALSE,
                                       tiny_return = TRUE)
      }})
    # Stop the cluster
    parallel::stopCluster(cl)
    
    # Stop the timer
    list_of_param[[typefamily]][['stumps_time']] = parallel_timer_stumps
    
    saveRDS(stumps_boot_models, file = 
                      paste("stumps_", type, "_", family$family, "_stnr_", stnr,
                            '_lambda_', penalty[['stumps']], "_bootstrap_", B[['stumps']], sep = ""))
              
              
              
    ### Plot the fitted cureves alongside with approximate and empirical confidence bands
    after_iter = round(seq(1, num_iter[['splines']], length.out = 5))
    file_name = paste(paste("stumps_", type, "_", family$family, "_stnr_", stnr,
                                '_lambda_', penalty[['stumps']], '_boot_', B[['stumps']],
                                    '_iter_', sep = ""),
                              paste(after_iter, collapse = '_'), sep = '')
            
    GAMBoost_stumps_plot_predictor_contribution(stumps_model,
                                                        type = type, 
                                                        boot_models = stumps_boot_models,
                                                        cn = training_data$cn,
                                                        after_iter = after_iter,
                                                        save_filename = file_name,
                                                        xmin = -1,
                                                        xmax = 1,
                                                        ymin = NULL,
                                                        ymax = NULL,
                                                        lwd = 1) 
            
            # Only look at the best iteration with respect to the test deviance.
            after_iter = stumps_model_optimal_test
            file_name = paste(paste("stumps_", type, "_", family$family, "_stnr_", stnr,
                                    '_lambda_', penalty[['stumps']], '_boot_', B[['stumps']],
                                    '_iter_', sep = ""),
                              paste(after_iter, collapse = '_'), sep = '')
            
            
            GAMBoost_stumps_plot_predictor_contribution(stumps_model,
                                                        type = type, 
                                                        boot_models = stumps_boot_models,
                                                        cn = training_data$cn,
                                                        after_iter = after_iter,
                                                        save_filename = file_name,
                                                        xmin = -1,
                                                        xmax = 1,
                                                        ymin = NULL,
                                                        ymax = NULL,
                                                        lwd = 1) 
            
            
            
            
            ### Look at the coverage for iteration stumps_model_optimal_test
            png(paste("stumps_", type, "_", family$family, "_stnr_", stnr,
                      '_lambda_', penalty[['stumps']], '_boot_', B[['stumps']],
                      '_iter_', stumps_model_optimal_test, '_coverage.png',  sep = ""),
                width = 3500, height = 1500, res = 350)
            
            stumps_coverage = GAMBoost_stumps_coverage_covariate(stumps_model, stumps_boot_models,
                                                                 after_iter = stumps_model_optimal_test, 
                                                                 type = type, cn = training_data$cn)
            
            dev.off()
            
            list_of_param[[typefamily]][['stumps_coverage_best_test_iter']] = stumps_coverage
            
            
            ### Look at development of coverage
            {
              number = 50
              stumps_approx_cov = matrix(NA, ncol = p, nrow = number)
              stumps_empirical_cov = matrix(NA, ncol = p, nrow = number)
              iters = round(seq(1, num_iter[['stumps']], length.out = number))
              for (ai in 1:number) {
                stumps_coverage = GAMBoost_stumps_coverage_covariate(stumps_model, stumps_boot_models,
                                                                     after_iter = iters[ai], 
                                                                     type = type, cn = training_data$cn)
                stumps_approx_cov[ai, ] = stumps_coverage$Approximate
                stumps_empirical_cov[ai, ] = stumps_coverage$Empirical
              }
              
              png(paste("stumps_", type, "_", family$family, "_stnr_", stnr,
                        '_lambda_', penalty[['stumps']], '_boot_', B[['stumps']], 
                        '_coverage_development.png',  sep = ""),
                  width = 3500, height = 1500, res = 350)
              
              par(mfrow=c(1,5), mar = c(3.5,3.5,2,1))
              for (s in 1:p) {
                plot(iters, stumps_approx_cov[,s], col = nice_colors[1], type = 'l', ylim = c(0,1), xlab = '', ylab = '')
                title(xlab = 'Bosting iter.', ylab = 'Coverage', line = 2.25)
                
                lines(iters, stumps_empirical_cov[,s], col = nice_colors[4])
                
                abline(h = 0.95, lty = 4)
              }
              dev.off()
            }
  }
}



# Save the model
saveRDS(list_of_param, 
        file = paste("list_of_param_stnr_", stnr, '_lambdastumps_', penalty[['stumps']],
                     '_lambdasplines_', penalty[['splines']], sep = ""))



list_of_param$smoothgaussian$stumps_time
list_of_param[[1]]$stumps_time
list_of_param$smoothgaussian$splines_time

  sapply(1:6, function(x) list_of_param[[x]]$stumps_time)










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



create_data_combined = function(n, seed_number = NA, cn = NA, stnr = NA, type = c("smooth", "step"),
                                family = gaussian()) {
  # Function that creates the data in accordance to 
  # the setup in 2.2.3 in GAM paper.
  
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
  
  # Create the observations.
  X = matrix(runif(n*5, -1, 1), nrow = n, ncol = 5)
  
  # Find the cn which yields the desired signal to noise ratio
  if (!is.na(stnr)) {
    if (type == 'smooth') {
      optim_func_smooth = function(c) {
        eta_temp = c*(-0.7 + X[,1] + 2*X[,3]^2 + sin(5*X[,5])) 
        mu_temp = family$linkinv(eta_temp)
        (stnr - sum((mu_temp - mean(mu_temp))^2) / sum(family$variance(mu_temp)))^2
      }
      cn = optimize(optim_func_smooth, c(0,100))$minimum
      
    } else {
      optim_func_step = function(c) {
        # Create the stepfunction f(x) which is -1 if x<= 0 and 1 if x>0.
        sf = stepfun(c(0), c(-1, 1), right=TRUE)
        # Create the linear predictor eta.
        eta_temp = c*(0.5*sf(X[,1]) + 0.25*sf(X[,3]) + sf(X[,5])) 
        mu_temp = family$linkinv(eta_temp)
        (stnr - sum((mu_temp - mean(mu_temp))^2) / sum(family$variance(mu_temp)))^2
      }
      cn = optimize(optim_func_step, c(0,100))$minimum
    }
  }
  
  if (type == 'smooth') {
    # Create the linear predictor eta.
    eta = cn*(-0.7 + X[,1] + 2*X[,3]^2 + sin(5*X[,5]))    
  }
  if (type == 'step') {
    # Create the stepfunction f(x) which is -1 if x<= 0 and 1 if x>0.
    sf = stepfun(c(0), c(-1, 1), right=TRUE)
    # Create the linear predictor eta.
    eta = cn*(0.5*sf(X[,1]) + 0.25*sf(X[,3]) + sf(X[,5])) 
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


GAMBoost_stumps_plot_predictor_contribution = function(model,
                                                       type = c('smooth', 'step'), 
                                                       boot_models = NULL,
                                                       cn = NULL,
                                                       after_iter = NULL,
                                                       save_filename = NULL,
                                                       xmin = -1,
                                                       xmax = 1,
                                                       ymin = NULL,
                                                       ymax = NULL,
                                                       lwd = 1) {
  
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
      
      # Plot the predicted contribution curves
      lines(X_new[,s], eta_contribution_values[,s], lwd = lwd, type = 's')
      
      # Approximate confidence bands
      lines(X_new[,s], lower, lty = 2, col = nice_colors[1], lwd = lwd, type = 's')
      lines(X_new[,s], upper, lty = 2, col = nice_colors[1], lwd = lwd, type = 's')
      
      # Add points to see where the training data is located
      points(jitter(model$X[,s]), rep(ymin, length(model$X[,s])), cex =.5, pch ="|", col = "darkgrey")
      
      # Add the empirical pointwise bootstrap confidence bands if provided
      if (!is.null(boot_models)) {
        lines(sort(model$X[,s]), boot_bands[[s]][1,], lty = 2, col = nice_colors[4], lwd = lwd, type = 's')
        lines(sort(model$X[,s]), boot_bands[[s]][2,], lty = 2, col = nice_colors[4], lwd = lwd, type = 's')
      }
      
      # Plot the true curve
      lines(x_values, true_y_values[,s], lty = 1, col='darkgray', lwd = lwd, type = 's')
      
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
                                                lwd = 1) {
  
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
      
      # Plot the predicted contribution curves
      lines(x_values, eta_contribution_values[,s], lwd = lwd)
      
      # Approximate confidence bands
      lines(ci_list$x, ci_list$lower, lty = 2, col = nice_colors[1], lwd = lwd)
      lines(ci_list$x, ci_list$upper, lty = 2, col = nice_colors[1], lwd = lwd)
      
      # Add points to see where the training data is located
      points(jitter(model$x[,s]), rep(ymin, n), cex =.5, pch ="|", col = "darkgrey")
      
      # Add the empirical pointwise bootstrap confidence bands if provided
      if (!is.null(boot_models)) {
        lines(sort(model$x[,s]), boot_bands[[s]][1,], lty = 2, col = nice_colors[4], lwd = lwd)
        lines(sort(model$x[,s]), boot_bands[[s]][2,], lty = 2, col = nice_colors[4], lwd = lwd)
      }
      
      # Plot the true curve
      lines(x_values, true_y_values[,s], lty = 1, col='darkgray', lwd = lwd)
      
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
