### File to create figures of fitting procedure of GAMBoost for
### a fixed penalty, cn = 3 and gaussian distribution
library(GAMBoost)

setwd("C:/Users/lars9/OneDrive/Master/Finalized R codes to put on GitHub")
source('GAMBoost_stumps.R', local = TRUE)
source('GAMBoost_splines.R', local = TRUE)
source('GAMBoost_common.R', local = TRUE)



library(parallel)
library(foreach)
library(snow)
library(doSNOW)

# rm(list=ls())
family = gaussian()

#
nice_colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)

colors = FALSE

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
data_list = list('training_data_smooth' = training_data_smooth,
                 'training_data_step'   = training_data_step,
                 'testing_data_smooth'  = testing_data_smooth,
                 'testing_data_step'    = testing_data_step)

saveRDS(data_list, 
        file = paste("splines_stumps_smooth_step_", family$family, "_cn_", cn, 
                     '_lambda_', penalty[['smooth']][['splines']], "_",
                     penalty[['smooth']][['stumps']], "_bootstrap_", B[['splines']],
                     "_",  B[['stumps']], sep = ""))



##### Fit the models to both sets of training data #####
{
  ### Fit the splines
  splines_model_smooth = GAMBoost(x        = training_data_smooth$X,
                                  y        = training_data_smooth$y,
                                  stepno   = num_iter[['splines']],
                                  penalty  = penalty[['smooth']][['splines']],
                                  family   = family,
                                  calc.hat = TRUE,
                                  calc.se  = TRUE,
                                  AIC.type = 'classical')
  
  splines_model_step = GAMBoost(x        = training_data_step$X,
                                y        = training_data_step$y,
                                stepno   = num_iter[['splines']],
                                penalty  = penalty[['step']][['splines']],
                                family   = family,
                                calc.hat = TRUE,
                                calc.se  = TRUE,
                                AIC.type = 'classical')
  
  ### Fit the stumps
  stumps_model_smooth = GAMBoost_stumps(X         = training_data_smooth$X, 
                                        y         = training_data_smooth$y,
                                        num_iter  = num_iter[['stumps']],
                                        lambda    = penalty[['smooth']][['stumps']],
                                        family    = family, 
                                        print_msg = TRUE)
  
  stumps_model_step = GAMBoost_stumps(X         = training_data_step$X, 
                                      y         = training_data_step$y,
                                      num_iter  = num_iter[['stumps']],
                                      lambda    = penalty[['step']][['stumps']],
                                      family    = family, 
                                      print_msg = TRUE)
}

#####################################
##### Test Deviance Development #####
#####################################
{
  ### Splines
  splines_model_deviance_test_smooth = rep(NA, num_iter[['splines']]+1)
  splines_model_deviance_test_step   = rep(NA, num_iter[['splines']]+1)
  
  for (iter in 1:(num_iter[['splines']]+1)) {
    # Smooth
    splines_model_deviance_test_smooth[iter] = sum(family$dev.resids(
      testing_data_smooth$y, 
      predict(splines_model_smooth, newdata = testing_data_smooth$X, at.step = (iter-1), type = "response"),
      rep(1, n_test)))
    
    # Step
    splines_model_deviance_test_step[iter] = sum(family$dev.resids(
      testing_data_step$y, 
      predict(splines_model_step, newdata = testing_data_step$X, at.step = (iter-1), type = "response"),
      rep(1, n_test)))
  }
  
  # Find the optimal test deviance, AIC and BIC
  splines_model_optimal_test_smooth = which.min(splines_model_deviance_test_smooth) - 1
  splines_model_optimal_AIC_smooth  = which.min(splines_model_smooth$AIC) - 1
  splines_model_optimal_BIC_smooth  = which.min(splines_model_smooth$BIC) - 1
  
  splines_model_optimal_test_step = which.min(splines_model_deviance_test_step) - 1
  splines_model_optimal_AIC_step  = which.min(splines_model_step$AIC) - 1
  splines_model_optimal_BIC_step  = which.min(splines_model_step$BIC) - 1
  
  
  ### Stumps
  stumps_model_deviance_test_smooth = rep(NA, num_iter[['stumps']]+1)
  stumps_model_deviance_test_step   = rep(NA, num_iter[['stumps']]+1)
  
  for (iter in 1:(num_iter[['stumps']]+1)) {
    cat(sprintf("Iteration %d of %d\n", iter, num_iter[['stumps']]+1))
    ## Smooth
    # Get the predicted results
    predicted_response = GAMBoost_stumps_predict(stumps_model_smooth, X_new = testing_data_smooth$X,
                                                 after_iter = (iter-1), version = 'response')

    # Compute the deviance
    stumps_model_deviance_test_smooth[iter] = sum(family$dev.resids(testing_data_smooth$y,
                                                                    predicted_response,
                                                                    rep(1, n_test)))

    ## Step
    # Get the predicted results
    predicted_response = GAMBoost_stumps_predict(stumps_model_step, X_new = testing_data_step$X,
                                                 after_iter = (iter-1), version = 'response')
    
    # Compute the deviance
    stumps_model_deviance_test_step[iter] = sum(family$dev.resids(testing_data_step$y,
                                                                  predicted_response,
                                                                  rep(1, n_test)))
  }
  # Find the optimal number of iterations
  stumps_model_optimal_test_smooth = which.min(stumps_model_deviance_test_smooth) - 1
  stumps_model_optimal_AIC_smooth = which.min(stumps_model_smooth$AIC) - 1
  stumps_model_optimal_BIC_smooth = which.min(stumps_model_smooth$BIC) - 1
  
  stumps_model_optimal_test_step = which.min(stumps_model_deviance_test_step) - 1
  stumps_model_optimal_AIC_step = which.min(stumps_model_step$AIC) - 1
  stumps_model_optimal_BIC_step = which.min(stumps_model_step$BIC) - 1
}

#####################################################
##### Crate a figure of the test deviance error #####
#####################################################
### Crate a figure of the test deviance error alongside the scaled AIC and BIC
{
  # Scaled AIC and BIC such that the initial value is the same as test deviance.
  splines_scaled_AIC_smooth = splines_model_smooth$AIC *
    splines_model_deviance_test_smooth[1]/splines_model_smooth$AIC[1]
  splines_scaled_BIC_smooth = splines_model_smooth$BIC *
    splines_model_deviance_test_smooth[1]/splines_model_smooth$BIC[1]
  splines_scaled_AIC_step = splines_model_step$AIC *
    splines_model_deviance_test_step[1]/splines_model_step$AIC[1]
  splines_scaled_BIC_step = splines_model_step$BIC *
    splines_model_deviance_test_step[1]/splines_model_step$BIC[1]
  
  splines_miny_smooth = min(splines_model_deviance_test_smooth, splines_scaled_AIC_smooth, splines_scaled_BIC_smooth)
  splines_maxy_smooth = max(splines_model_deviance_test_smooth, splines_scaled_AIC_smooth, splines_scaled_BIC_smooth)
  splines_miny_step   = min(splines_model_deviance_test_step, splines_scaled_AIC_step, splines_scaled_BIC_step)
  splines_maxy_step   = max(splines_model_deviance_test_step, splines_scaled_AIC_step, splines_scaled_BIC_step)
  
  # SCaled AIC and BIC such that the initial value is the same as test deviance.
  stumps_scaled_AIC_smooth = stumps_model_smooth$AIC * stumps_model_deviance_test_smooth[1]/stumps_model_smooth$AIC[1]
  stumps_scaled_BIC_smooth = stumps_model_smooth$BIC * stumps_model_deviance_test_smooth[1]/stumps_model_smooth$BIC[1]
  stumps_scaled_AIC_step   = stumps_model_step$AIC * stumps_model_deviance_test_step[1]/stumps_model_step$AIC[1]
  stumps_scaled_BIC_step   = stumps_model_step$BIC * stumps_model_deviance_test_step[1]/stumps_model_step$BIC[1]
  
  stumps_miny_smooth = min(stumps_model_deviance_test_smooth, stumps_scaled_AIC_smooth, stumps_scaled_BIC_smooth)
  stumps_maxy_smooth = max(stumps_model_deviance_test_smooth, stumps_scaled_AIC_smooth, stumps_scaled_BIC_smooth)
  stumps_miny_step   = min(stumps_model_deviance_test_step, stumps_scaled_AIC_step, stumps_scaled_BIC_step)
  stumps_maxy_step   = max(stumps_model_deviance_test_step, stumps_scaled_AIC_step, stumps_scaled_BIC_step)
  
  miny_smooth  = min(stumps_miny_smooth, splines_miny_smooth)
  maxy_smooth  = min(stumps_maxy_smooth, splines_maxy_smooth)
  range_smooth = maxy_smooth - miny_smooth
  
  miny_step  = min(stumps_miny_step, splines_miny_step)
  maxy_step  = min(stumps_maxy_step, splines_maxy_step)
  range_step = maxy_step - miny_step
}

png(paste("GAMBoost_cn3_lambda_4_300_Gaussian_test_deviance.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
  
  # Plot 1
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_smooth-0.1*range_smooth), maxy_smooth+0.1*range_smooth))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Smooth Splines'))
  
  lines(0:num_iter[['splines']], splines_model_deviance_test_smooth, lty = 1, col = nice_colors[1], lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_AIC_smooth, lty = 1, col = nice_colors[2], lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_BIC_smooth, lty = 1, col = nice_colors[3], lwd = 1)

  abline(v = splines_model_optimal_test_smooth, lty = 2, col = nice_colors[1])
  abline(v = splines_model_optimal_AIC_smooth,  lty = 3, col = nice_colors[2])
  abline(v = splines_model_optimal_BIC_smooth,  lty = 4, col = nice_colors[3])
  
  legend('topright', legend = c(paste('Test(', splines_model_optimal_test_smooth, ')', sep = ''),
                                paste('AIC(', splines_model_optimal_AIC_smooth, ')', sep = ''),
                                paste('BIC(', splines_model_optimal_BIC_smooth, ')', sep = '')),
         title = as.expression(bquote('Min Test ='~.(round(min(splines_model_deviance_test_smooth),0)))),
         lty = 1, col = nice_colors, cex = 1)
  
  
  # Plot 2
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_smooth-0.1*range_smooth), maxy_smooth+0.1*range_smooth))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Smooth Stumps'))
  
  lines(0:num_iter[['stumps']], stumps_model_deviance_test_smooth, lty = 1, col = nice_colors[1], lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_AIC_smooth, lty = 1, col = nice_colors[2], lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_BIC_smooth, lty = 1, col = nice_colors[3], lwd = 1)
  
  abline(v = stumps_model_optimal_test_smooth, lty = 2, col = nice_colors[1])
  abline(v = stumps_model_optimal_AIC_smooth,  lty = 3, col = nice_colors[2])
  abline(v = stumps_model_optimal_BIC_smooth,  lty = 4, col = nice_colors[3])
  
  legend('topright', legend = c(paste('Test(', stumps_model_optimal_test_smooth, ')', sep = ''),
                                paste('AIC(', stumps_model_optimal_AIC_smooth, ')', sep = ''),
                                paste('BIC(', stumps_model_optimal_BIC_smooth, ')', sep = '')),
         title = as.expression(bquote('Min Test ='~.(round(min(stumps_model_deviance_test_smooth),0)))),
         lty = 1, col = nice_colors, cex = 1)
  
  
  # Plot 3
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_step-0.1*range_step), maxy_step+0.1*range_step))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Stepwise Splines'))
  
  lines(0:num_iter[['splines']], splines_model_deviance_test_step, lty = 1, col = nice_colors[1], lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_AIC_step, lty = 1, col = nice_colors[2], lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_BIC_step, lty = 1, col = nice_colors[3], lwd = 1)
  
  abline(v = splines_model_optimal_test_step, lty = 2, col = nice_colors[1])
  abline(v = splines_model_optimal_AIC_step,  lty = 3, col = nice_colors[2])
  abline(v = splines_model_optimal_BIC_step,  lty = 4, col = nice_colors[3])
  
  legend('topright', legend = c(paste('Test(', splines_model_optimal_test_step, ')', sep = ''),
                                paste('AIC(', splines_model_optimal_AIC_step, ')', sep = ''),
                                paste('BIC(', splines_model_optimal_BIC_step, ')', sep = '')),
         title = as.expression(bquote('Min Test ='~.(round(min(splines_model_deviance_test_step),0)))),
         lty = 1, col = nice_colors, cex = 1)
  
  
  # Plot 4
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_step-0.1*range_step), maxy_step+0.1*range_step))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Stepwise Stumps'))
  
  lines(0:num_iter[['stumps']], stumps_model_deviance_test_step, lty = 1, col = nice_colors[1], lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_AIC_step, lty = 1, col = nice_colors[2], lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_BIC_step, lty = 1, col = nice_colors[3], lwd = 1)
  
  abline(v = stumps_model_optimal_test_step, lty = 2, col = nice_colors[1])
  abline(v = stumps_model_optimal_AIC_step,  lty = 3, col = nice_colors[2])
  abline(v = stumps_model_optimal_BIC_step,  lty = 4, col = nice_colors[3])
  
  legend('topright', legend = c(paste('Test(', stumps_model_optimal_test_step, ')', sep = ''),
                                paste('AIC(', stumps_model_optimal_AIC_step, ')', sep = ''),
                                paste('BIC(', stumps_model_optimal_BIC_step, ')  ', sep = '')),
         title = as.expression(bquote('Min Test ='~.(round(min(stumps_model_deviance_test_step),0)))),
         lty = 1, col = nice_colors, cex = 1)
} 
dev.off()


png(paste("GAMBoost_cn3_lambda_4_300_Gaussian_test_deviance_colorless.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
  
  # Plot 1
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_smooth-0.1*range_smooth), maxy_smooth+0.1*range_smooth))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Smooth Splines'))
  
  lines(0:num_iter[['splines']], splines_model_deviance_test_smooth, lty = 1, col = 'black', lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_AIC_smooth, lty = 2, col = 'black', lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_BIC_smooth, lty = 3, col = 'black', lwd = 1)
  
  abline(v = splines_model_optimal_test_smooth, lty = 1, col = 'black')
  abline(v = splines_model_optimal_AIC_smooth,  lty = 2, col = 'black')
  abline(v = splines_model_optimal_BIC_smooth,  lty = 3, col = 'black')
  
  l = legend('topright', legend = c(paste('Test(', splines_model_optimal_test_smooth, ')', sep = ''),
                                paste('AIC(', splines_model_optimal_AIC_smooth, ')', sep = ''),
                                paste('BIC(', splines_model_optimal_BIC_smooth, ')', sep = '')),
             title = as.expression(bquote('Min Test ='~.(round(min(splines_model_deviance_test_smooth),0)))),
             lty = c(1,2,3), col = 'black', cex = 1)
  
  # temp = legend(-0100, -1000, legend = NA,  title = as.expression(bquote('Min Test ='~.(round(min(splines_model_deviance_test_smooth),0)))))
  # legend(l$rect$left + (l$rect$w - temp$rect$w), l$rect$top - l$rect$h, cex = 1,
  #        legend = NA,  title = as.expression(bquote('Min Test ='~.(round(min(splines_model_deviance_test_smooth),0)))))

  
  
  
  # Plot 2
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_smooth-0.1*range_smooth), maxy_smooth+0.1*range_smooth))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Smooth Stumps'))
  
  lines(0:num_iter[['stumps']], stumps_model_deviance_test_smooth, lty = 1, col = 'black', lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_AIC_smooth, lty = 2, col = 'black', lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_BIC_smooth, lty = 3, col = 'black', lwd = 1)
  
  abline(v = stumps_model_optimal_test_smooth, lty = 1, col = 'black')
  abline(v = stumps_model_optimal_AIC_smooth,  lty = 2, col = 'black')
  abline(v = stumps_model_optimal_BIC_smooth,  lty = 3, col = 'black')
  
  legend('topright', legend = c(paste('Test(', stumps_model_optimal_test_smooth, ')', sep = ''),
                                paste('AIC(', stumps_model_optimal_AIC_smooth, ')', sep = ''),
                                paste('BIC(', stumps_model_optimal_BIC_smooth, ')', sep = '')),
         title = as.expression(bquote('Min Test ='~.(round(min(stumps_model_deviance_test_smooth),0)))),
         lty = c(1,2,3), col = 'black', cex = 1)
  
  
  # Plot 3
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_step-0.1*range_step), maxy_step+0.1*range_step))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Stepwise Splines'))
  
  lines(0:num_iter[['splines']], splines_model_deviance_test_step, lty = 1, col = 'black', lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_AIC_step, lty = 2, col = 'black', lwd = 1)
  lines(0:num_iter[['splines']], splines_scaled_BIC_step, lty = 3, col = 'black', lwd = 1)
  
  abline(v = splines_model_optimal_test_step, lty = 1, col = 'black')
  abline(v = splines_model_optimal_AIC_step,  lty = 2, col = 'black')
  abline(v = splines_model_optimal_BIC_step,  lty = 3, col = 'black')
  
  legend('topright', legend = c(paste('Test(', splines_model_optimal_test_step, ')', sep = ''),
                                paste('AIC(', splines_model_optimal_AIC_step, ')', sep = ''),
                                paste('BIC(', splines_model_optimal_BIC_step, ')', sep = '')),
         title = as.expression(bquote('Min Test ='~.(round(min(splines_model_deviance_test_step),0)))),
         lty = 1:3, col = 'black', cex = 1)
  
  
  # Plot 4
  plot(0:num_iter[['splines']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_step-0.1*range_step), maxy_step+0.1*range_step))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2.25)
  title(main = paste('Gaussian Stepwise Stumps'))
  
  lines(0:num_iter[['stumps']], stumps_model_deviance_test_step, lty = 1, col = 'black', lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_AIC_step, lty = 2, col = 'black', lwd = 1)
  lines(0:num_iter[['stumps']], stumps_scaled_BIC_step, lty = 3, col = 'black', lwd = 1)
  
  abline(v = stumps_model_optimal_test_step, lty = 1, col = 'black')
  abline(v = stumps_model_optimal_AIC_step,  lty = 2, col = 'black')
  abline(v = stumps_model_optimal_BIC_step,  lty = 3, col = 'black')
  
  legend('topright', legend = c(paste('Test(', stumps_model_optimal_test_step, ')', sep = ''),
                                paste('AIC(', stumps_model_optimal_AIC_step, ')', sep = ''),
                                paste('BIC(', stumps_model_optimal_BIC_step, ')  ', sep = '')),
         title = as.expression(bquote('Min Test ='~.(round(min(stumps_model_deviance_test_step),0)))),
         lty = c(1,2,3), col = 'black', cex = 1)
} 
dev.off()


#######################################################
##### Create the Bootstrap models for all 4 cases #####
#######################################################
# Find the number of available cores
cores=detectCores() 

### First we do the splines
{
  ## Smooth
  cat(sprintf("Currently fitting bootstraped versions of Splines Smooth\n"))
  # Set up the clusters
  cl = parallel::makeCluster(cores[1]-1) #not to overload your computer
  registerDoSNOW(cl)
  
  # Iterate over the B[['splines']] bootstrap samples 
  splines_parallel_timer_smooth = system.time({
    splines_boot_models_smooth = foreach(b = 1:B[['splines']], .export = 'GAMBoost') %dopar% {
      #cat(sprintf('Splines %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
      
      # The indices for this bootstrap sample
      boot_idx = sample(1:n_train, n_train, replace = TRUE)
      
      # Get the relevant data for b
      boot_X_train = training_data_smooth$X[boot_idx, ]
      boot_y_train = training_data_smooth$y[boot_idx]
      
      # Fit the model (100 iterations = 1.5 Mb)
      boot_pspline = GAMBoost(x = boot_X_train, 
                              y = boot_y_train,
                              stepno   = num_iter[['splines']],
                              penalty  = penalty[['smooth']][['splines']],
                              family   = family,
                              calc.hat = FALSE, # Do not need these as we are only interested in the 
                              calc.se  = FALSE) # the predicted curves.
    }})
  
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Save bootstrap models
  saveRDS(splines_boot_models_smooth, 
          file = paste("splines_smooth_", family$family, "_cn_", cn, 
                       '_lambda_', penalty[['smooth']][['splines']],
                       "_bootstrap_", B[['splines']], sep = ""))
  
  ## Step
  cat(sprintf("Currently fitting bootstraped versions of Splines Step\,"))
  # Set up the clusters
  cl = parallel::makeCluster(cores[1]-1) #not to overload your computer
  registerDoSNOW(cl)
  
  # Iterate over the B[['splines']] bootstrap samples 
  splines_parallel_timer_step = system.time({
    splines_boot_models_step = foreach(b = 1:B[['splines']], .export = 'GAMBoost') %dopar% {
      #cat(sprintf('Splines %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
      
      # The indices for this bootstrap sample
      boot_idx = sample(1:n_train, n_train, replace = TRUE)
      
      # Get the relevant data for b
      boot_X_train = training_data_step$X[boot_idx, ]
      boot_y_train = training_data_step$y[boot_idx]
      
      # Fit the model (100 iterations = 1.5 Mb)
      boot_pspline = GAMBoost(x = boot_X_train, 
                              y = boot_y_train,
                              stepno   = num_iter[['splines']],
                              penalty  = penalty[['step']][['splines']],
                              family   = family,
                              calc.hat = FALSE, # Do not need these as we are only interested in the 
                              calc.se  = FALSE) # the predicted curves.
    }})
  
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Save bootstrap models
  saveRDS(splines_boot_models_step, 
          file = paste("splines_step_", family$family, "_cn_", cn, 
                       '_lambda_', penalty[['step']][['splines']],
                       "_bootstrap_", B[['splines']], sep = ""))
  
}

### Then we do the stumps
{
  ## Smooth
  cat(sprintf("Currently fitting bootstraped versions of Stumps Smooth\n"))
  # Set up the clusters
  cl = parallel::makeCluster(cores[1]-1) #not to overload your computer
  registerDoSNOW(cl)
  
  # Iterate over the B[['stumps']] bootstrap samples 
  stumps_parallel_timer_smooth = system.time({
    stumps_boot_models_smooth = foreach(b = 1:B[['stumps']]) %dopar% {
      #cat(sprintf('stumps %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
      
      # The indices for this bootstrap sample
      boot_idx = sample(1:n_train, n_train, replace = TRUE)
      
      # Get the relevant data for b
      boot_X_train = training_data_smooth$X[boot_idx, ]
      boot_y_train = training_data_smooth$y[boot_idx]
      
      # Fit the model (250 iterations = 2MB)
      boot_pstump  = GAMBoost_stumps(X = boot_X_train, 
                                     y = boot_y_train,
                                     num_iter  = num_iter[['stumps']],
                                     lambda    = penalty[['smooth']][['stumps']],
                                     family    = family, 
                                     print_msg = FALSE,
                                     tiny_return = TRUE)
    }})
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Save bootstrap models
  saveRDS(stumps_boot_models_smooth, 
          file = paste("stumps_smooth_", family$family, "_cn_", cn, 
                       '_lambda_', penalty[['smooth']][['stumps']],
                       "_bootstrap_", B[['stumps']], sep = ""))
  
  
  ## Step
  cat(sprintf("Currently fitting bootstraped versions of Stumps Step\n"))
  # Set up the clusters
  cl = parallel::makeCluster(cores[1]-1) #not to overload your computer
  registerDoSNOW(cl)
  
  # Iterate over the B[['stumps']] bootstrap samples 
  stumps_parallel_timer_steps = system.time({
    stumps_boot_models_step = foreach(b = 1:B[['stumps']]) %dopar% {
      #cat(sprintf('stumps %10s %10s bootstrap sample: %-4d\n', type, family$family, b))
      
      # The indices for this bootstrap sample
      boot_idx = sample(1:n_train, n_train, replace = TRUE)
      
      # Get the relevant data for b
      boot_X_train = training_data_step$X[boot_idx, ]
      boot_y_train = training_data_step$y[boot_idx]
      
      # Fit the model (250 iterations = 2MB)
      boot_pstump  = GAMBoost_stumps(X = boot_X_train, 
                                     y = boot_y_train,
                                     num_iter  = num_iter[['stumps']],
                                     lambda    = penalty[['step']][['stumps']],
                                     family    = family, 
                                     print_msg = FALSE,
                                     tiny_return = TRUE)
    }})
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Save bootstrap models
  saveRDS(stumps_boot_models_step, 
          file = paste("stumps_step_", family$family, "_cn_", cn, 
                       '_lambda_', penalty[['step']][['stumps']],
                       "_bootstrap_", B[['stumps']], sep = ""))
}

data_list = list('training_data_smooth' = training_data_smooth,
                 'training_data_step'   = training_data_step,
                 'testing_data_smooth'  = testing_data_smooth,
                 'testing_data_step'    = testing_data_step,
                 'stumps_parallel_timer_smooth' = stumps_parallel_timer_smooth,
                 'stumps_parallel_timer_steps' = stumps_parallel_timer_steps,
                 'splines_parallel_timer_smooth' = splines_parallel_timer_smooth,
                 'splines_parallel_timer_step' = splines_parallel_timer_step)

saveRDS(data_list, 
        file = paste("splines_stumps_smooth_step_", family$family, "_cn_", cn, 
                     '_lambda_Smooth_', penalty[['smooth']][['splines']], "_", penalty[['smooth']][['stumps']],
                     '_lambda_step_', penalty[['step']][['splines']], "_", penalty[['step']][['stumps']],
                     "_bootstrap_", B[['splines']], "_",  B[['stumps']], sep = ""))


# If the we have saveded the bootstrap lists, we can load them
# sourcepath = "C:/Users/lars9/OneDrive/Master/R_codes/GAMBoost_final/Saves/"
# 
# splines_boot_models_smooth = readRDS(paste(sourcepath,"splines_smooth_", family$family, "_cn_", cn,
#                                                  '_lambda_', penalty[['smooth']][['splines']],
#                                                  "_bootstrap_", B[['splines']], sep = ""))
# 
# splines_boot_models_step =readRDS(paste(sourcepath,"splines_step_", family$family, "_cn_", cn,
#                                           '_lambda_', penalty[['step']][['splines']],
#                                           "_bootstrap_", B[['splines']], sep = ""))
# 
# stumps_boot_models_smooth = readRDS(paste(sourcepath, "stumps_smooth_", family$family, "_cn_", cn,
#                                           '_lambda_', penalty[['smooth']][['stumps']],
#                                           "_bootstrap_", B[['stumps']], sep = ""))
# 
# stumps_boot_models_step = readRDS(paste(sourcepath, "stumps_step_", family$family, "_cn_", cn,
#                                         '_lambda_', penalty[['step']][['stumps']],
#                                         "_bootstrap_", B[['stumps']], sep = ""))

################################################
##### Create the curve development figures #####
################################################
splines_after_iter_smooth = c(5,15,50,250, splines_model_optimal_test_smooth)
splines_after_iter_step   = c(5,15,50,250, splines_model_optimal_test_step)
stumps_after_iter_smooth  = c(5,15,50,250, stumps_model_optimal_test_smooth)
stumps_after_iter_step    = c(5,15,50,250, stumps_model_optimal_test_step)

### Splines smooth
file_name = paste(paste("splines_smooth_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['smooth']][['splines']], '_boot_', B[['splines']],
                        '_iter_', sep = ""),
                  paste(splines_after_iter_smooth, collapse = '_'), sep = '')

GAMBoost_plot_predictor_contribution(splines_model_smooth,
                                     type = 'smooth', 
                                     boot_models = splines_boot_models_smooth,
                                     cn = training_data_smooth$cn,
                                     after_iter = splines_after_iter_smooth,
                                     save_filename = file_name,
                                     xmin = -1,
                                     xmax = 1,
                                     ymin = -4.65,
                                     ymax = 4.65,
                                     lwd = 1,
                                     color = TRUE) 

file_name = paste(paste("splines_smooth_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['smooth']][['splines']], '_boot_', B[['splines']],
                        '_iter_', sep = ""),
                  paste(splines_after_iter_smooth, collapse = '_'), "_colorless2", sep = '')

GAMBoost_plot_predictor_contribution(splines_model_smooth,
                                     type = 'smooth', 
                                     boot_models = splines_boot_models_smooth,
                                     cn = training_data_smooth$cn,
                                     after_iter = splines_after_iter_smooth,
                                     save_filename = file_name,
                                     xmin = -1,
                                     xmax = 1,
                                     ymin = -4.65,
                                     ymax = 4.65,
                                     lwd = 1,
                                     color = FALSE) 


### Splines step
file_name = paste(paste("splines_step_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['step']][['splines']], '_boot_', B[['splines']],
                        '_iter_', sep = ""),
                  paste(splines_after_iter_step, collapse = '_'), sep = '')

GAMBoost_plot_predictor_contribution(splines_model_step,
                                     type = 'step', 
                                     boot_models = splines_boot_models_step,
                                     cn = training_data_step$cn,
                                     after_iter = splines_after_iter_step,
                                     save_filename = file_name,
                                     xmin = -1,
                                     xmax = 1,
                                     ymin = -4,
                                     ymax = 4,
                                     lwd = 1,
                                     color = TRUE) 

file_name = paste(paste("splines_step_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['step']][['splines']], '_boot_', B[['splines']],
                        '_iter_', sep = ""),
                  paste(splines_after_iter_step, collapse = '_'), "_colorless2", sep = '')

GAMBoost_plot_predictor_contribution(splines_model_step,
                                     type = 'step', 
                                     boot_models = splines_boot_models_step,
                                     cn = training_data_step$cn,
                                     after_iter = splines_after_iter_step,
                                     save_filename = file_name,
                                     xmin = -1,
                                     xmax = 1,
                                     ymin = -4,
                                     ymax = 4,
                                     lwd = 1,
                                     color = FALSE) 


### Stump smooth
file_name = paste(paste("stumps_smooth_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['smooth']][['stumps']], '_boot_', B[['stumps']],
                        '_iter_', sep = ""),
                  paste(stumps_after_iter_smooth, collapse = '_'), sep = '')

GAMBoost_stumps_plot_predictor_contribution(stumps_model_smooth,
                                            type = 'smooth', 
                                            boot_models = stumps_boot_models_smooth,
                                            cn = training_data_smooth$cn,
                                            after_iter = stumps_after_iter_smooth,
                                            save_filename = file_name,
                                            xmin = -1,
                                            xmax = 1,
                                            ymin = -4.65,
                                            ymax = 4.65,
                                            lwd = 1,
                                            color = TRUE) 

file_name = paste(paste("stumps_smooth_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['smooth']][['stumps']], '_boot_', B[['stumps']],
                        '_iter_', sep = ""),
                  paste(stumps_after_iter_smooth, collapse = '_'), '_colorless2', sep = '')

GAMBoost_stumps_plot_predictor_contribution(stumps_model_smooth,
                                            type = 'smooth', 
                                            boot_models = stumps_boot_models_smooth,
                                            cn = training_data_smooth$cn,
                                            after_iter = stumps_after_iter_smooth,
                                            save_filename = file_name,
                                            xmin = -1,
                                            xmax = 1,
                                            ymin = -4.65,
                                            ymax = 4.65,
                                            lwd = 1,
                                            color = FALSE) 


### Stumps step
file_name = paste(paste("stumps_step_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['step']][['stumps']], '_boot_', B[['stumps']],
                        '_iter_', sep = ""),
                  paste(stumps_after_iter_step, collapse = '_'), sep = '')

GAMBoost_stumps_plot_predictor_contribution(stumps_model_step,
                                            type = 'step', 
                                            boot_models = stumps_boot_models_step,
                                            cn = training_data_step$cn,
                                            after_iter = stumps_after_iter_step,
                                            save_filename = file_name,
                                            xmin = -1,
                                            xmax = 1,
                                            ymin = -4,
                                            ymax = 4,
                                            lwd = 1,
                                            color = TRUE) 


file_name = paste(paste("stumps_step_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['step']][['stumps']], '_boot_', B[['stumps']],
                        '_iter_', sep = ""),
                  paste(stumps_after_iter_step, collapse = '_'), '_colorless2', sep = '')

GAMBoost_stumps_plot_predictor_contribution(stumps_model_step,
                                            type = 'step', 
                                            boot_models = stumps_boot_models_step,
                                            cn = training_data_step$cn,
                                            after_iter = stumps_after_iter_step,
                                            save_filename = file_name,
                                            xmin = -1,
                                            xmax = 1,
                                            ymin = -4,
                                            ymax = 4,
                                            lwd = 1,
                                            color = FALSE) 

# 
GAMBoost_stumps_plot_predictor_contribution(stumps_model_step,
                                            type = 'step',
                                            boot_models = NULL,# stumps_boot_models_step,
                                            cn = training_data_step$cn,
                                            after_iter = c(2,4,6,8),
                                            save_filename = file_name,
                                            xmin = -1,
                                            xmax = 1,
                                            ymin = -1.25,
                                            ymax = 1.25,
                                            lwd = 1,
                                            color = FALSE)

########################################
##### CREATE PLOT OF THE CI FOR MU #####
########################################
png(paste("GAMBoost_cn3_lambda_4_300_Gaussian_CI_mu.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
  coverage_splines_smooth = GAMBoost_splines_inside_conf_bands(model = splines_model_smooth,
                                                               boot_models = splines_boot_models_smooth,
                                                               cn = training_data_smooth$cn, 
                                                               training_data = training_data_smooth,
                                                               after_iter = splines_model_optimal_test_smooth,
                                                               type = 'smooth',
                                                               colors = TRUE,
                                                               ymin = -8,
                                                               ymax = 10)
  
  coverage_stumps_smooth = GAMBoost_stumps_inside_conf_bands(model = stumps_model_smooth,
                                                             boot_models = stumps_boot_models_smooth,
                                                             cn = training_data_smooth$cn, 
                                                             training_data = training_data_smooth,
                                                             after_iter = stumps_model_optimal_test_smooth,
                                                             type = 'smooth',
                                                             colors = TRUE,
                                                             ymin = -8,
                                                             ymax = 10)
  
  coverage_splines_step = GAMBoost_splines_inside_conf_bands(model = splines_model_step,
                                                             boot_models = splines_boot_models_step,
                                                             cn = training_data_step$cn, 
                                                             training_data = training_data_step,
                                                             after_iter = splines_model_optimal_test_step,
                                                             type = 'step',
                                                             colors = TRUE,
                                                             ymin = -7,
                                                             ymax = 7)
  
  coverage_stumps_step = GAMBoost_stumps_inside_conf_bands(model = stumps_model_step,
                                                           boot_models = stumps_boot_models_step,
                                                           cn = training_data_step$cn, 
                                                           training_data = training_data_step,
                                                           after_iter = stumps_model_optimal_test_step,
                                                           type = 'step',
                                                           colors = TRUE,
                                                           ymin = -7,
                                                           ymax = 7)
}
dev.off()



png(paste("GAMBoost_cn3_lambda_4_300_Gaussian_CI_mu_colorless.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
  coverage_splines_smooth = GAMBoost_splines_inside_conf_bands(model = splines_model_smooth,
                                                               boot_models = splines_boot_models_smooth,
                                                               cn = training_data_smooth$cn, 
                                                               training_data = training_data_smooth,
                                                               after_iter = splines_model_optimal_test_smooth,
                                                               type = 'smooth',
                                                               colors = FALSE,
                                                               ymin = -8,
                                                               ymax = 10)
  
  coverage_stumps_smooth = GAMBoost_stumps_inside_conf_bands(model = stumps_model_smooth,
                                                             boot_models = stumps_boot_models_smooth,
                                                             cn = training_data_smooth$cn, 
                                                             training_data = training_data_smooth,
                                                             after_iter = stumps_model_optimal_test_smooth,
                                                             type = 'smooth',
                                                             colors = FALSE,
                                                             ymin = -8,
                                                             ymax = 10)
  
  coverage_splines_step = GAMBoost_splines_inside_conf_bands(model = splines_model_step,
                                                             boot_models = splines_boot_models_step,
                                                             cn = training_data_step$cn, 
                                                             training_data = training_data_step,
                                                             after_iter = splines_model_optimal_test_step,
                                                             type = 'step',
                                                             colors = FALSE,
                                                             ymin = -7,
                                                             ymax = 7)
  
  coverage_stumps_step = GAMBoost_stumps_inside_conf_bands(model = stumps_model_step,
                                                           boot_models = stumps_boot_models_step,
                                                           cn = training_data_step$cn, 
                                                           training_data = training_data_step,
                                                           after_iter = stumps_model_optimal_test_step,
                                                           type = 'step',
                                                           colors = FALSE,
                                                           ymin = -7,
                                                           ymax = 7)
}
dev.off()

png(paste("GAMBoost_cn3_lambda_4_300_Gaussian_CI_mu_error_bands_colorless.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(2,2), mar = c(3.5,3.5,2,1))
  coverage_splines_smooth = GAMBoost_splines_inside_conf_bands(model = splines_model_smooth,
                                                               boot_models = splines_boot_models_smooth,
                                                               cn = training_data_smooth$cn, 
                                                               training_data = training_data_smooth,
                                                               after_iter = splines_model_optimal_test_smooth,
                                                               type = 'smooth',
                                                               colors = FALSE,
                                                               ymin = -8,
                                                               ymax = 10,
                                                               figure_type = 2)
  
  coverage_stumps_smooth = GAMBoost_stumps_inside_conf_bands(model = stumps_model_smooth,
                                                             boot_models = stumps_boot_models_smooth,
                                                             cn = training_data_smooth$cn, 
                                                             training_data = training_data_smooth,
                                                             after_iter = stumps_model_optimal_test_smooth,
                                                             type = 'smooth',
                                                             colors = FALSE,
                                                             ymin = -8,
                                                             ymax = 10,
                                                             figure_type = 2)
  
  coverage_splines_step = GAMBoost_splines_inside_conf_bands(model = splines_model_step,
                                                             boot_models = splines_boot_models_step,
                                                             cn = training_data_step$cn, 
                                                             training_data = training_data_step,
                                                             after_iter = splines_model_optimal_test_step,
                                                             type = 'step',
                                                             colors = FALSE,
                                                             ymin = -7,
                                                             ymax = 7,
                                                             figure_type = 2)
  
  coverage_stumps_step = GAMBoost_stumps_inside_conf_bands(model = stumps_model_step,
                                                           boot_models = stumps_boot_models_step,
                                                           cn = training_data_step$cn, 
                                                           training_data = training_data_step,
                                                           after_iter = stumps_model_optimal_test_step,
                                                           type = 'step',
                                                           colors = FALSE,
                                                           ymin = -7,
                                                           ymax = 7,
                                                           figure_type = 2)
}
dev.off()
