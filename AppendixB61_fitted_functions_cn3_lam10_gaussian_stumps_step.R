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
num_iter = list('stumps' = 250)
#penalty  = list('splines' = 300, 'stumps' = 4)
penalty  = list('step' = list('stumps' = 10))

# Number of bootstrap samples for each model
B = list('stumps' = 1000)

seed_train = 2020
seed_test = 2030
### Setup: generate data
{
  training_data_step = create_data_combined(n_train, seed_number = seed_train, cn = cn,
                                            type = 'step', family = family)

  testing_data_step  = create_data_combined(n_test, seed_number = seed_test, cn = cn,
                                            type = 'step', family = family)  
}
data_list = list('training_data_step'   = training_data_step,
                 'testing_data_step'    = testing_data_step)

saveRDS(data_list, 
        file = paste("stumps_step_", family$family, "_cn_", cn, 
                     '_lambda_', penalty[['step']][['stumps']],
                     "_bootstrap_", B[['stumps']], sep = ""))



##### Fit the models to both sets of training data #####
{
  ### Fit the stumps
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
  ### Stumps
  stumps_model_deviance_test_step   = rep(NA, num_iter[['stumps']]+1)
  
  for (iter in 1:(num_iter[['stumps']]+1)) {
    cat(sprintf("Iteration %d of %d\n", iter, num_iter[['stumps']]+1))
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
  stumps_model_optimal_test_step = which.min(stumps_model_deviance_test_step) - 1
  stumps_model_optimal_AIC_step = which.min(stumps_model_step$AIC) - 1
  stumps_model_optimal_BIC_step = which.min(stumps_model_step$BIC) - 1
}

#####################################################
##### Crate a figure of the test deviance error #####
#####################################################
### Crate a figure of the test deviance error alongside the scaled AIC and BIC
{
  # SCaled AIC and BIC such that the initial value is the same as test deviance.
  stumps_scaled_AIC_step   = stumps_model_step$AIC * stumps_model_deviance_test_step[1]/stumps_model_step$AIC[1]
  stumps_scaled_BIC_step   = stumps_model_step$BIC * stumps_model_deviance_test_step[1]/stumps_model_step$BIC[1]
  
  stumps_miny_step   = min(stumps_model_deviance_test_step, stumps_scaled_AIC_step, stumps_scaled_BIC_step)
  stumps_maxy_step   = max(stumps_model_deviance_test_step, stumps_scaled_AIC_step, stumps_scaled_BIC_step)
  
  miny_step  = min(stumps_miny_step)
  maxy_step  = min(stumps_maxy_step)
  range_step = maxy_step - miny_step
}

png(paste("GAMBoost_cn3_lambda_10_300_Gaussian_step_stumps_test_deviance.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
  # Plot 4
  plot(0:num_iter[['stumps']], type='n', ylab = "", xlab = "",
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


png(paste("GAMBoost_cn3_lambda_4_300_Gaussian_step_stumps_test_deviance_colorless.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
  # Plot 4
  plot(0:num_iter[['stumps']], type='n', ylab = "", xlab = "",
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

### Then we do the stumps
{
  ## Step
  cat(sprintf("Currently fitting bootstraped versions of Stumps Step\n"))
  # Set up the clusters
  cl = parallel::makeCluster(cores[1]-1) #not to overload your computer
  registerDoSNOW(cl)
  
  # # Create a progress bar
  # pb = txtProgressBar(max = B[['stumps']], style = 3)
  # progress = function(n) setTxtProgressBar(pb, n)
  # opts = list(progress = progress)
  
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
  # # Stop the progress bar
  # close(pb)
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Save bootstrap models
  saveRDS(stumps_boot_models_step, 
          file = paste("stumps_step_", family$family, "_cn_", cn, 
                       '_lambda_', penalty[['step']][['stumps']],
                       "_bootstrap_", B[['stumps']], sep = ""))
}

stumps_boot_models_step = readRDS(paste("stumps_step_", family$family, "_cn_", cn, 
                                        '_lambda_', penalty[['step']][['stumps']],
                                        "_bootstrap_", B[['stumps']], sep = ""))

# data_list = list('training_data_smooth' = training_data_smooth,
#                  'training_data_step'   = training_data_step,
#                  'testing_data_smooth'  = testing_data_smooth,
#                  'testing_data_step'    = testing_data_step,
#                  'stumps_parallel_timer_smooth' = stumps_parallel_timer_smooth,
#                  'stumps_parallel_timer_steps' = stumps_parallel_timer_steps,
#                  'splines_parallel_timer_smooth' = splines_parallel_timer_smooth,
#                  'splines_parallel_timer_step' = splines_parallel_timer_step)
# 
# saveRDS(data_list,
#         file = paste("splines_stumps_smooth_step_", family$family, "_cn_", cn,
#                      '_lambda_Smooth_', penalty[['smooth']][['splines']], "_", penalty[['smooth']][['stumps']],
#                      '_lambda_step_', penalty[['step']][['splines']], "_", penalty[['step']][['stumps']],
#                      "_bootstrap_", B[['splines']], "_",  B[['stumps']], sep = ""))


# If the we have saveded the bootstrap lists, we can load them
# sourcepath = "C:/Users/lars9/OneDrive/Master/R_codes/GAMBoost_final/"
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
stumps_after_iter_step    = c(5,15,50,250, stumps_model_optimal_test_step)

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

#### FOR LOW ITERATIONS #####
stumps_after_iter_step = c(2,4,6,8,11)
file_name = paste(paste("stumps_step_", family$family, "_cn_", cn,
                        '_lambda_', penalty[['step']][['stumps']], '_include_boot_', B[['stumps']],
                        '_iter_', sep = ""),
                  paste(stumps_after_iter_step, collapse = '_'), '_colorless', sep = '')
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


########################################
##### CREATE PLOT OF THE CI FOR MU #####
########################################
png(paste("GAMBoost_cn3_lambda_10_Gaussian_CI_mu.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(1, 1), mar = c(3.5,3.5,2,1))
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



png(paste("GAMBoost_cn3_lambda_10_Gaussian_CI_mu_colorless.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
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

png(paste("GAMBoost_cn3_lambda_10_Gaussian_CI_mu_error_bands_colorless.png", sep = ""), width = 3500, height = 2250, res = 350)
{
  par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
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


png(paste("GAMBoost_cn3_lambda_4_300_Gaussian_step_stumps_test_deviance_and_CI_mu_colorless3.png", sep = ""),
    width = 3500, height = 1250, res = 350)
{
  par(mfrow=c(1,2), mar = c(3.5,3.5,2,1))
  # Plot 4
  plot(0:num_iter[['stumps']], type='n', ylab = "", xlab = "",
       ylim = c(max(0, miny_step-0.1*range_step), maxy_step+0.1*range_step))
  title(xlab = 'Iteration', ylab = 'Test Deviance / Scaled AIC and BIC', line = 2)
  title(main = bquote(plain('Gaussian Stepwise Stumps')), cex.main = 1.13)
  
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
         lty = c(1,2,3), col = 'black', cex = 0.8)

  
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
