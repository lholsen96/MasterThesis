###
source('GAMBoost_stumps.R', local = TRUE)
source('GAMBoost_common.R', local = TRUE)

# rm(list=ls())
# Include libraries
library(GAMBoost)
library(parallel)
library(foreach)
library(snow)
library(doSNOW)

#
nice_colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)
colors = FALSE


B_outer = 5
B_inner = 20
stepno = 250

n_train = 100
n_test = 1000
p = 5

# types 
types = c('smooth', 'step')

# Signal to noise ratios
stnrs = c(1,3,10)

# Set family
# Define possible distributions and function types
families = c('gaussian', 'binomial', 'poisson')

# set penatly lambda
lambda = 300

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
          
          # Fit model
          splines_model = GAMBoost(x        = training_data$X,
                                   y        = training_data$y,
                                   stepno   = stepno,
                                   penalty  = lambda,
                                   family   = family,
                                   calc.hat = TRUE,
                                   calc.se  = TRUE)
          
          # Find best iter based on lowest evaluation error
          splines_model_deviance_test = rep(NA, stepno+1)
          for (iter in 1:(stepno+1)) {
            splines_model_deviance_test[iter] = sum(family$dev.resids(
              testing_data$y, 
              predict(splines_model, newdata = testing_data$X, at.step = (iter-1), type = "response"),
              rep(1, n_test)))
          }
          # Find the optimal test deviance 
          best_iterations[bo] = which.min(splines_model_deviance_test) - 1 # Subtract 1 since intercept is iter 0.
          best_iterations_dev[bo] = min(splines_model_deviance_test)
          
          # Fit a new GAMBoost model to the best iteration to obtain the hat matrix
          splines_model = GAMBoost(x        = training_data$X,
                                   y        = training_data$y,
                                   stepno   = best_iterations[bo],
                                   penalty  = lambda,
                                   family   = family,
                                   calc.hat = TRUE,
                                   calc.se  = TRUE)
          
          ##### Want to create B_inner Bootstrap samples
          ### Parallel bootstraps
          #Find the number of available cores
          cores=detectCores() 
          
          # Set up the clusters
          cl = parallel::makeCluster(cores[1] - 1) # not to overload your computer
          registerDoSNOW(cl)
          
          # Iterate over the B_inner bootstrap samples 
          parallel_timer = system.time({
            splines_boot_models = foreach(b = 1:B_inner, .export = 'GAMBoost') %dopar% {
              
              # The indices for this bootstrap sample
              boot_idx = sample(1:n_train, n_train, replace = TRUE)
              
              # Get the relevant data for b
              boot_X_train = training_data$X[boot_idx, ]
              boot_y_train = training_data$y[boot_idx]
              
              # Fit the model (100 iterations = 1.5 Mb)
              boot_pspline = GAMBoost(x = boot_X_train, 
                                      y = boot_y_train,
                                      stepno   = best_iterations[bo],
                                      penalty  = lambda,
                                      family   = family,
                                      calc.hat = FALSE, # Do not need these as we are only interested in the 
                                      calc.se  = FALSE) # the predicted curves.
            }})
          
          # Stop the cluster
          parallel::stopCluster(cl)
          
          # record the timer
          times[bo,] = parallel_timer[1:3]
          
          # Just to look at the fitted model
          GAMBoost_plot_predictor_contribution(splines_model,
                                               type = type, 
                                               boot_models = splines_boot_models,
                                               cn = training_data$cn,
                                               after_iter = best_iterations[bo],
                                               save_filename = NULL,
                                               xmin = -1,
                                               xmax = 1,
                                               ymin = NULL,
                                               ymax = NULL,
                                               lwd = 1) 
          
          
          # check if the true line is inside the bands
          if (bo == 1) {
            png(paste("splines_", type, "_", family$family, "_stnr_", STNR,
                      '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                      '_CI_mu.png',  sep = ""),
                width = 3500, height = 1500, res = 350)
          }
          par(mfrow = c(1,1), mar = c(3.5,3.5,2,1))
          coverage_temp = GAMBoost_splines_inside_conf_bands(splines_model,
                                                             splines_boot_models,
                                                             cn = training_data$cn, 
                                                             training_data = training_data,
                                                             after_iter = best_iterations[bo],
                                                             type = type,
                                                             colors = FALSE)

          if (bo == 1) {
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
      
      png(paste("splines_", type, "_", family$family, "_stnr_", STNR,
                '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                '_coverage_final.png',  sep = ""),
          width = 3500, height = 1500, res = 350)
      
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
        
        print(c(mean(coverage_approx_temp), mean(coverage_empir_temp)))
      }
      dev.off()
      
      png(paste("splines_", type, "_", family$family, "_stnr_", STNR,
                '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                '_coverage_mu.png',  sep = ""),
          width = 3500, height = 1500, res = 350)
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
      }
      dev.off()
      

        
      saveRDS(list('times' = times, 'training' = training_data, 'testing' = testing_data,
                   'inside_band_app' = inside_band_approximate, 'inside_band_emp'= inside_band_empirical,
                   'best_iters' = best_iterations, 'family' = family, 'type' = type, 'lambda' = lambda,
                   'stnr' = STNR, 'B_out' = B_outer, 'B_in' = B_inner, 'stepno' = stepno,
                   'best_iters_dev' = best_iterations_dev, 'loop_timer' = loop_timer, 
                   'mu_inside_interval_app' = mu_inside_interval_approximate,
                   'mu_inside_interval_emp' = mu_inside_interval_empirical)
        ,
        paste("splines_", type, "_", family$family, "_stnr_", STNR,
                '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                "_final3", sep='')
      )
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
lambdas = c(30, 125, 500)
B_outer = 500
B_inner = 500

# Create a list to store all the results and extract them fromt the files
# Around xxxMb
splines = list()
for (family in families) {
  splines[[family]] = list()
  for (type in types) {
    splines[[family]][[type]] = list()
    for (STNR in stnrs) {
      splines[[family]][[type]][[paste('STNR',STNR,sep = '')]] = list()
      for (lambda in lambdas) {
        splines[[family]][[type]][[paste('STNR',STNR,sep = '')]][[paste('lam',lambda,sep = '')]] =
          readRDS(paste("splines_", type, "_", family, "_stnr_", STNR,
                        '_lambda_', lambda, '_bootinner_', B_inner,
                        '_bootouter_', B_outer, "_final2", sep=''))
      }
    }
  }
}


#####################################################################
##### CREATE FIGURES DESCRIBING THE GAUSSIAN CASE WITH SNTR = 3 #####
#####################################################################
### Splines
splines_gaussian_stnr3_lam125_smooth = splines$gaussian$smooth$STNR3$lam125
splines_gaussian_stnr3_lam125_step   = splines$gaussian$step$STNR3$lam125

### Stumps
setwd("C:/Users/lars9/OneDrive/Master/R_codes/StumpsCoverage/With Bootstrap/")
stumps_gaussian_stnr3_lam4_smooth = 
  readRDS('stumps_smooth_gaussian_stnr_3_lambda_4_bootinner_200_bootouter_200_final2')

stumps_gaussian_stnr3_lam4_step =
  readRDS('stumps_step_gaussian_stnr_3_lambda_4_bootinner_500_bootouter_500_final3')



##### Plot of the coverage of f_j functions
colors = FALSE
png(paste('splines_stumps_gaussian_smoth_step_stnr_3_lambda_125_4_coverage_f_colorless.png',  sep = ""),
    width = 3500, height = 3250, res = 350)
par(mfrow=c(4,5), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
{
  ##### Smooth Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_approx_temp, col = nice_colors[1])
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4])          
    } else {
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_approx_temp, lty = 1, lwd = 1.5, col = 'black')
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_empir_temp,  lty = 2, lwd = 1.5, col = "#676767")
    }
  }
  
  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_approx_temp, col = nice_colors[1])
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4])          
    } else {
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_approx_temp, lty = 1, lwd = 1.5, col = 'black')
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_empir_temp,  lty = 2, lwd = 1.5, col = "#676767")
    }
  }
  
  
  ##### Stepwise Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_approx_temp, col = nice_colors[1])
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4])          
    } else {
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_approx_temp, lty = 1, lwd = 1.5, col = 'black')
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_empir_temp,  lty = 2, lwd = 1.5, col = "#676767")
    }
  }
  
  
  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_approx_temp, col = nice_colors[1])
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4])          
    } else {
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_approx_temp, lty = 1, lwd = 1.5, col = 'black')
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_empir_temp,  lty = 2, lwd = 1.5, col = "#676767")
    }
  }
}
dev.off()




#######################################################################
##### Plot of the coverage of f_j functions with lines and points #####
#######################################################################
colors = FALSE
png(paste('splines_stumps_gaussian_smoth_step_stnr_3_lambda_125_4_coverage_f_lines_points_colorless.png',  sep = ""),
    width = 3500, height = 3250, res = 350)
par(mfrow=c(4,5), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
{
  ##### Smooth Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4], lwd = 0.75) 
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)  

    } else {
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_approx_temp, col = 'black', lwd = 0.75)
      lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_empir_temp, col = "darkgray", lwd = 0.75)
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
            coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  
  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], lwd = 0.75)   
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_approx_temp, col = 'black', lwd = 0.75)
      lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
            coverage_empir_temp, col = "darkgray", lwd = 0.75)
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  
  
  ##### Stepwise Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], lwd = 0.75)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_approx_temp, col = 'black', lwd = 0.75)
      lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_empir_temp, col = "darkgray", lwd = 0.75)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
            coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  
  
  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], lwd = 0.75)  
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_approx_temp, col = 'black', lwd = 0.75)
      lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_empir_temp, col = "darkgray", lwd = 0.75)
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
            coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
}
dev.off()


#############################################################
##### Plot of the coverage of f_j functions with points #####
#############################################################
colors = FALSE
png(paste('splines_stumps_gaussian_smoth_step_stnr_3_lambda_125_4_coverage_f_points_colorless.png',  sep = ""),
    width = 3500, height = 3250, res = 350)
par(mfrow=c(4,5), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
{
  ##### Smooth Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  
  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  
  
  ##### Stepwise Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  
  
  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
}
dev.off()


########################################################################################
##### Plot of the coverage of f_j functions with points. Divide splines and stumps #####
########################################################################################
colors = FALSE
lines = TRUE
png(paste('splines_gaussian_smooth_step_stnr_3_lambda_125_coverage_f_lines_points_colorless.png',  sep = ""),
    width = 3500, height = 1625, res = 350)
par(mfrow=c(2,5), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
{
  ##### Smooth Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      if (lines) {
        lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
               coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
        lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
               coverage_empir_temp,  col = nice_colors[4], lwd = 0.75)
      }
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      if (lines) {
        lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
               coverage_approx_temp, col = 'black', lwd = 0.75)
        lines(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
               coverage_empir_temp, col = "darkgray", lwd = 0.75)
      }
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_smooth$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  


  ##### Stepwise Data
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Splines", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(splines_gaussian_stnr3_lam125_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      if (lines) {
        lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
               coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
        lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
               coverage_empir_temp,  col = nice_colors[4], lwd = 0.75)   
      }
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      if (lines) {
        lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
               coverage_approx_temp, col = 'black', lwd = 0.75)
        lines(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
               coverage_empir_temp, col = "darkgray", lwd = 0.75)
      }
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(splines_gaussian_stnr3_lam125_step$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
}
dev.off()


##### STUMPS
colors = FALSE
lines = TRUE
png(paste('stumps_gaussian_smooth_step_stnr_3_lambda_4_coverage_f_lines_points_colorless.png',  sep = ""),
    width = 3500, height = 1625, res = 350)
par(mfrow=c(2,5), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
{
  ##### Smooth data
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_smooth$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Smooth Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      if (lines) {
        lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
               coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
        lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
               coverage_empir_temp,  col = nice_colors[4], lwd = 0.75)
      }
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      if (lines) {
        lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
               coverage_approx_temp, col = 'black', lwd = 0.75)
        lines(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
               coverage_empir_temp, col = "darkgray", lwd = 0.75)
      }
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_smooth$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
  
  ##### Stepwise data
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_step$B_out
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = 
      apply(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = "Stepwise Stumps", line = 0.5, cex.main = 1.5)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    
    points(jitter(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
           rep(0, n_train), cex =.65, pch ="|", col = "darkgray")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      if (lines) {
        lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
               coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
        lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
               coverage_empir_temp,  col = nice_colors[4], lwd = 0.75) 
      }
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
    } else {
      if (lines) {
        lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
               coverage_approx_temp, col = 'black', lwd = 0.75)
        lines(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
               coverage_empir_temp, col = "darkgray", lwd = 0.75)
      }
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
      points(sort(stumps_gaussian_stnr3_lam4_step$training$X[,s]),
             coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
    }
  }
}
dev.off()






######################################
##### Plot of the coverage of mu #####
######################################
colors = FALSE
lines = TRUE
lwd = 0.75
png(paste('splines_stumps_gaussian_smoth_step_stnr_3_lambda_125_4_coverage_mu_lines_points_colorless.png',  sep = ""),
    width = 3500, height = 2250, res = 350)
par(mfrow=c(2,2), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
{
  ##### Smooth
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_smooth$B_out
  coverage_approx_mu =
    apply(splines_gaussian_stnr3_lam125_smooth$mu_inside_interval_app, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$mu_inside_interval_app[,1])))
  
  coverage_empircal_mu = 
    apply(splines_gaussian_stnr3_lam125_smooth$mu_inside_interval_emp, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_smooth$mu_inside_interval_emp[,1])))
  
  plot(NA, type = 'n', xlim = c(1, n_train), ylim = c(0,1), xlab = '', ylab = '')
  title(xlab = bquote("Observation"), ylab =  bquote("Coverage"~mu), line = 2.25)
  title(main = "Smooth Splines", line = 0.5)
  abline(h = 0.95, lty = 1, col =  "#D3D3D3", lwd = 2)
  
  if (colors) {
    if (lines) {
      lines(1:n_train, coverage_approx_mu,    col = nice_colors[1])
      lines(1:n_train, coverage_empircal_mu,  col = nice_colors[4])
    }
    points(1:n_train, coverage_approx_mu,   col = nice_colors[1], pch = 19)
    points(1:n_train, coverage_empircal_mu, col = nice_colors[4], pch = 19) 
  } else {
    if (lines) {
      lines(1:n_train,  coverage_approx_mu,   col = 'black',   lty = 1,  lwd = lwd)
      lines(1:n_train,  coverage_empircal_mu, col = "darkgray", lty = 1,  lwd = lwd)
    }
    points(1:n_train, coverage_approx_mu,   col = 'black',   pch = 19, cex = 0.6)
    points(1:n_train, coverage_empircal_mu, col = "darkgray", pch = 19, cex = 0.6)
  }

  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_smooth$B_out
  coverage_approx_mu =
    apply(stumps_gaussian_stnr3_lam4_smooth$mu_inside_interval_app, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$mu_inside_interval_app[,1])))
  
  coverage_empircal_mu = 
    apply(stumps_gaussian_stnr3_lam4_smooth$mu_inside_interval_emp, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_smooth$mu_inside_interval_emp[,1])))
  
  plot(NA, type = 'n', xlim = c(1, n_train), ylim = c(0,1), xlab = '', ylab = '')
  title(xlab = bquote("Observation"), ylab =  bquote("Coverage"~mu), line = 2.25)
  title(main = "Smooth Stumps", line = 0.5)
  abline(h = 0.95, lty = 1, col =  "#D3D3D3", lwd = 2)
  
  if (colors) {
    if (lines) {
      lines(1:n_train, coverage_approx_mu,    col = nice_colors[1])
      lines(1:n_train, coverage_empircal_mu,  col = nice_colors[4])
    }
    points(1:n_train, coverage_approx_mu,   col = nice_colors[1], pch = 19)
    points(1:n_train, coverage_empircal_mu, col = nice_colors[4], pch = 19) 
  } else {
    if (lines) {
      lines(1:n_train,  coverage_approx_mu,   col = 'black',   lty = 1,  lwd = lwd)
      lines(1:n_train,  coverage_empircal_mu, col = "darkgray", lty = 1,  lwd = lwd)
    }
    points(1:n_train, coverage_approx_mu,   col = 'black',   pch = 19, cex = 0.6)
    points(1:n_train, coverage_empircal_mu, col = "darkgray", pch = 19, cex = 0.6)
  }
  
  
  
  
  ##### Step
  ### Splines
  B_outer = splines_gaussian_stnr3_lam125_step$B_out
  coverage_approx_mu =
    apply(splines_gaussian_stnr3_lam125_step$mu_inside_interval_app, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$mu_inside_interval_app[,1])))
  
  coverage_empircal_mu = 
    apply(splines_gaussian_stnr3_lam125_step$mu_inside_interval_emp, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(splines_gaussian_stnr3_lam125_step$mu_inside_interval_emp[,1])))
  
  plot(NA, type = 'n', xlim = c(1, n_train), ylim = c(0,1), xlab = '', ylab = '')
  title(xlab = bquote("Observation"), ylab =  bquote("Coverage"~mu), line = 2.25)
  title(main = "Stepwise Splines", line = 0.5)
  abline(h = 0.95, lty = 1, col =  "#D3D3D3", lwd = 2)
  
  if (colors) {
    if (lines) {
      lines(1:n_train, coverage_approx_mu,    col = nice_colors[1])
      lines(1:n_train, coverage_empircal_mu,  col = nice_colors[4]) 
    }
    points(1:n_train, coverage_approx_mu,   col = nice_colors[1], pch = 19)
    points(1:n_train, coverage_empircal_mu, col = nice_colors[4], pch = 19) 
  } else {
    if (lines) {
      lines(1:n_train,  coverage_approx_mu,   col = 'black',   lty = 1,  lwd = lwd)
      lines(1:n_train,  coverage_empircal_mu, col = "darkgray", lty = 1,  lwd = lwd)
    }
    points(1:n_train, coverage_approx_mu,   col = 'black',   pch = 19, cex = 0.6)
    points(1:n_train, coverage_empircal_mu, col = "darkgray", pch = 19, cex = 0.6)
  }
  
  
  ### Stumps
  B_outer = stumps_gaussian_stnr3_lam4_step$B_out
  coverage_approx_mu =
    apply(stumps_gaussian_stnr3_lam4_step$mu_inside_interval_app, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$mu_inside_interval_app[,1])))
  
  coverage_empircal_mu = 
    apply(stumps_gaussian_stnr3_lam4_step$mu_inside_interval_emp, 2, sum, na.rm = TRUE) / 
    (B_outer - sum(is.na(stumps_gaussian_stnr3_lam4_step$mu_inside_interval_emp[,1])))
  
  plot(NA, type = 'n', xlim = c(1, n_train), ylim = c(0,1), xlab = '', ylab = '')
  title(xlab = bquote("Observation"), ylab =  bquote("Coverage"~mu), line = 2.25)
  title(main = "Stepwise Stumps", line = 0.5)
  abline(h = 0.95, lty = 1, col =  "#D3D3D3", lwd = 2)
  
  if (colors) {
    if (lines) {
      lines(1:n_train, coverage_approx_mu,    col = nice_colors[1])
      lines(1:n_train, coverage_empircal_mu,  col = nice_colors[4]) 
    }
    points(1:n_train, coverage_approx_mu,   col = nice_colors[1], pch = 19)
    points(1:n_train, coverage_empircal_mu, col = nice_colors[4], pch = 19) 
  } else {
    if (lines) {
      lines(1:n_train,  coverage_approx_mu,   col = 'black',   lty = 1,  lwd = lwd)
      lines(1:n_train,  coverage_empircal_mu, col = "darkgray", lty = 1,  lwd = lwd)
    }
    points(1:n_train, coverage_approx_mu,   col = 'black',   pch = 19, cex = 0.6)
    points(1:n_train, coverage_empircal_mu, col = "darkgray", pch = 19, cex = 0.6)
  }
}
dev.off()

























############################
##### EARLIER VERSIONS #####
############################

splines$gaussian$smooth$STNR3$lam125

family = gaussian()
lambda = 125
SNTR = 3
par(mfrow=c(2,5), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
for (type in types) {
  temp = splines[[family$family]][[type]][[paste('STNR',STNR,sep = '')]][[paste('lam',lambda,sep = '')]]
  
  # Create the title
  if (temp$training$type == 'smooth') {
    main_tit = "Smooth Splines"
  } else {
    main_tit = "Stepwise Splines"
  }
  
  for (s in 1:p) {
    # Find the coverage for the n_train values
    coverage_approx_temp = apply(temp$inside_band_app[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(temp$inside_band_app[[s]][,1])))
    
    coverage_empir_temp = apply(temp$inside_band_emp[[s]], 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(temp$inside_band_emp[[s]][,1])))
    
    
    plot(NA, type = 'n', xlim = c(-1,1), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("x"[.(s)]), ylab = 'Coverage', line = 2.25)
    if (s == 3) {
      #title(main = bquote('Coverage f'[.(s)]), line = 1)
      title(main = main_tit, line = 0.5, cex.main = 0.95)
      #title(main = paste("    ", main_tit), line = -1, outer = TRUE)
    }
    points(jitter(temp$training$X[,s]), rep(0, n_train), cex =.65, pch ="|", col = "#D3D3D3")
    abline(h = 0.95, lty = 1, col = "#D3D3D3", lwd = 2)
    if (colors) {
      lines(sort(temp$traininga$X[,s]), coverage_approx_temp, col = nice_colors[1])
      lines(sort(training$X[,s]), coverage_empir_temp,  col = nice_colors[4])          
    } else {
      lines(sort(temp$training$X[,s]), coverage_approx_temp, lty = 1, lwd = 1.5, col = 'black')
      lines(sort(temp$training$X[,s]), coverage_empir_temp,  lty = 2, lwd = 1.5, col = "#676767")
    }
  }
}


png(paste("splines_", type, "_", family$family, "_stnr_", STNR,
          '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
          '_coverage_mu_ver2.png',  sep = ""),
    width = 3500, height = 1500, res = 350)
par(mfrow=c(1,2), mar = c(3.5,3.5,2,1), oma = c(0,0,0,0))
for (type in types) {
  temp = splines[[family$family]][[type]][[paste('STNR',STNR,sep = '')]][[paste('lam',lambda,sep = '')]]
  
  # Create the title
  if (temp$training$type == 'smooth') {
    main_tit = "Smooth Splines"
  } else {
    main_tit = "Stepwise Splines"
  }
  
  {
   # par(mfrow=c(1,1), mar = c(3.5,3.5,2,1))
    coverage_approx_mu = apply(temp$mu_inside_interval_app, 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(temp$mu_inside_interval_app[,1])))
    
    coverage_empircal_mu = apply(temp$mu_inside_interval_emp, 2, sum, na.rm = TRUE) / 
      (B_outer - sum(is.na(temp$mu_inside_interval_emp[,1])))
    
    plot(NA, type = 'n', xlim = c(1, n_train), ylim = c(0,1), xlab = '', ylab = '')
    title(xlab = bquote("Observation"), ylab =  bquote("Coverage"~mu), line = 2.25)
    title(main = main_tit, line = 0.65)
    abline(h = 0.95, lty = 1, col = "darkgray", lwd = 2)
    
    if (colors) {
      #lines(1:n_train, coverage_approx_mu,    col = nice_colors[1])
      #lines(1:n_train, coverage_empircal_mu,  col = nice_colors[4])
      points(1:n_train, coverage_approx_mu,   col = nice_colors[1], pch = 19)
      points(1:n_train, coverage_empircal_mu, col = nice_colors[4], pch = 19) 
    } else {
      #lines(1:n_train,  coverage_approx_mu,   col = 'black',   lty = 1,  lwd = 1.5)
      points(1:n_train, coverage_approx_mu,   col = 'black',   pch = 19, cex = 0.5)
      #lines(1:n_train,  coverage_empircal_mu, col = "gray", lty = 2,  lwd = 1.5)
      points(1:n_train, coverage_empircal_mu, col = "gray", pch = 19, cex = 0.5)
    }
  }
}
dev.off()
