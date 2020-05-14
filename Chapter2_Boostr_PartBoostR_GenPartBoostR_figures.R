# rm(list = ls())
setwd("C:/Users/lars9/OneDrive/Master/R_codes/PartBoostR")

source_path = "C:/Users/lars9/OneDrive/Master/R_codes/PartBoostR/"
source(paste(source_path, "PartBoostR_source_code.R", sep=""))


library(MASS)
library(RColorBrewer)


##################################
# Let the true model be y = b1x1 + b2x2 + b3x3 + b4x4^2
beta_true = c(3, -1, 0, 2) 

# Number of parameters 
p = length(beta_true)
sigma_low = 1
sigma_high = 5
n_train = 100
n_test = 1000

# prefix to all the file_names generated in this code.
save_prefix = "Final"

# Define some colors
nice_colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)

# Set seed for reproducibility
set.seed(421)

{
  # Old method for independent data
  X_train = matrix(runif(n_train*p, -1, 1), ncol = p, nrow = n_train)
  
  X_test = matrix(runif(n_test*p, -1, 1), ncol = p, nrow = n_test)
  
  # Generate the normal response data with variance sigma
  y_train_low = rnorm(n_train, mean = X_train %*% beta_true, sd = sigma_low)
  y_test_low = rnorm(n_test, mean = X_test %*% beta_true, sd = sigma_low)
  
  y_train_high = rnorm(n_train, mean = X_train %*% beta_true, sd = sigma_high)
  y_test_high = rnorm(n_test, mean = X_test %*% beta_true, sd = sigma_high)
  
  # Collect them into dataframes for the linear model to work
  training = data.frame(response=y_train_low, x=X_train)
  testing = data.frame(response=y_test_low, x=X_test)
  
  training_high = data.frame(response=y_train_high, x=X_train)
  testing_high = data.frame(response=y_test_high, x=X_test)
  
  # Find the least squares solutions of the unormalzied covariates
  lm_model = lm(response~.-1, data = training)
  lm_model$coefficients
  
  lm_model_high = lm(response~.-1, data = training_high)
  lm_model_high$coefficients
  
  # Look at the training and test error, mse for unscaled version
  lm_mse_train = mean((y_train_low - predict(lm_model, training))^2)
  lm_mse_test  = mean((y_test_low - predict(lm_model, testing))^2)
  c(lm_mse_train, lm_mse_test)
  
  lm_mse_train_high = mean((y_train_high - predict(lm_model_high, training_high))^2)
  lm_mse_test_high  = mean((y_test_high - predict(lm_model_high, testing_high))^2)
  c(lm_mse_train_high, lm_mse_test_high)

}


### Now We can do the modelfitting
# We standardize
X_train_sc = scale(X_train)
# Use the values from X_train to scale the test set. 
X_test_sc = scale(X_test, attr(X_train_sc, "scaled:center"), attr(X_train_sc, "scaled:scale"))

# Collect the normalized sets
training_sc = data.frame(response=y_train_low, x=X_train_sc)
testing_sc = data.frame(response=y_test_low, x=X_test_sc)

training_sc_high = data.frame(response=y_train_high, x=X_train_sc)
testing_sc_high = data.frame(response=y_test_high, x=X_test_sc)

# Find the least squares solutions to the scaled version
lm_model_sc = lm(response~.-1, data = training_sc)
lm_model_sc$coefficients

lm_model_sc_high = lm(response~.-1, data = training_sc_high)
lm_model_sc_high$coefficients

# Save the least squares solution
ls_solutions = lm_model_sc$coefficients
ls_solutions_high = lm_model_sc_high$coefficients

# OLS is not equivariant wrt scaling. Scaling X with a => scaling \beta_hat with 1/a
lm_model_sc$coefficients / attr(X_train_sc, "scaled:scale")
lm_model$coefficients

lm_model_sc_high$coefficients / attr(X_train_sc, "scaled:scale")
lm_model_high$coefficients

# Look at the training and test error, mse
lm_mse_train = mean((y_train_low - predict(lm_model_sc, training_sc))^2)
lm_mse_test  = mean((y_test_low - predict(lm_model_sc, testing_sc))^2)
c(lm_mse_train, lm_mse_test)

# Look at the training and test error, mse
lm_mse_train_high = mean((y_train_high - predict(lm_model_sc_high, training_sc_high))^2)
lm_mse_test_high  = mean((y_test_high - predict(lm_model_sc_high, testing_sc_high))^2)
c(lm_mse_train_high, lm_mse_test_high)

# Fit a BoostR and a componentwise PartBoostR
num_iter = 2000
lambda = 5000
lambda_values = c(1000, 5000, 10000, 25000) # Some values of lambda

# Fit the models with sigma_low
{
  # BoostR
  mod_BoostR_1 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter, lambda=lambda_values[1],
                          return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  mod_BoostR_2 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter, lambda=lambda_values[2],
                            return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  mod_BoostR_3 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter, lambda=lambda_values[3],
                          return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  mod_BoostR_4 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter, lambda=lambda_values[4],
                          return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  # PartBoostR
  mod_PartBoostR_1 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter, 
                                lambda=lambda_values[1], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
  
  mod_PartBoostR_2 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter, 
                                lambda=lambda_values[2], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
  
  mod_PartBoostR_3 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter,
                                lambda=lambda_values[3], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
  
  mod_PartBoostR_4 = PartBoostR(X_train_sc, y_train_low, num_iter = num_iter,
                                lambda=lambda_values[4], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
}

# Fit the models with sigma_high
{
  # BoostR
  mod_BoostR_high_1 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter, lambda=lambda_values[1],
                            return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  mod_BoostR_high_2 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter, lambda=lambda_values[2],
                            return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  mod_BoostR_high_3 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter, lambda=lambda_values[3],
                            return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  mod_BoostR_high_4 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter, lambda=lambda_values[4],
                            return_hatmat = FALSE, V_compulsory = 1:p, print = FALSE)
  
  # PartBoostR
  mod_PartBoostR_high_1 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter, 
                                lambda=lambda_values[1], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
  
  mod_PartBoostR_high_2 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter, 
                                lambda=lambda_values[2], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
  
  mod_PartBoostR_high_3 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter,
                                lambda=lambda_values[3], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
  
  mod_PartBoostR_high_4 = PartBoostR(X_train_sc, y_train_high, num_iter = num_iter,
                                lambda=lambda_values[4], return_hatmat = FALSE,
                                V_compulsory = c(), print = FALSE)
}


# Save the coefficient build up
file_name = paste(save_prefix, "PartBoostR_BoostR_coeff_buildup_lambda_", lambda, sep = "")
plot_coeff_build_up(list(mod_BoostR_2, mod_PartBoostR_2), ls_solutions, save_plot_name = file_name)
plot_coeff_build_up(list(mod_BoostR_2, mod_PartBoostR_2), ls_solutions)

# See how many updates each parameters have recieved
table(mod_PartBoostR_1$update_trace)
table(mod_PartBoostR_2$update_trace)
table(mod_PartBoostR_3$update_trace)
table(mod_PartBoostR_4$update_trace)
# And in which iteration was the first time the param got updated
match(unique(mod_PartBoostR_1$update_trace), mod_PartBoostR_1$update_trace)
match(unique(mod_PartBoostR_2$update_trace), mod_PartBoostR_2$update_trace)
match(unique(mod_PartBoostR_3$update_trace), mod_PartBoostR_3$update_trace)
match(unique(mod_PartBoostR_4$update_trace), mod_PartBoostR_4$update_trace)


#################################################################################
#### Want to show that num_iter is the main tuning parameter and not lambda #####
#################################################################################
# Save the two previous images together
png(paste(save_prefix, "BoostR_PartBoostR_coeff_buildup.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow=c(1,2))
  plot_coeff_build_up(list(mod_BoostR_1, mod_BoostR_2, mod_BoostR_3, mod_BoostR_4), ls_solutions,
                      max_iter = num_iter, title = "BoostR With Different Regularization")
  
  plot_coeff_build_up(list(mod_PartBoostR_1, mod_PartBoostR_2, mod_PartBoostR_3, mod_PartBoostR_4),
                      ls_solutions, max_iter = num_iter, title = "PartBoostR With Different Regularization")
}
dev.off()

# sigma_high
png(paste(save_prefix, "BoostR_PartBoostR_coeff_buildup_sigma_high.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow=c(1,2))
  plot_coeff_build_up(list(mod_BoostR_high_1, mod_BoostR_high_2, mod_BoostR_high_3, mod_BoostR_high_4),
                      ls_solutions_high, max_iter = num_iter, title = "BoostR With Different Regularization")
  
  plot_coeff_build_up(list(mod_PartBoostR_high_1, mod_PartBoostR_high_2, mod_PartBoostR_high_3, mod_PartBoostR_high_4),
                      ls_solutions_high, max_iter = num_iter, title = "PartBoostR With Different Regularization")
}
dev.off()


# # With only 250 to highlight the values
# png(paste(save_prefix, "BoostR_PartBoostR_coeff_buildup_maxiter.png", sep = ""), width = 3250, height = 1500, res = 350)
# {
#   par(mfrow=c(1,2))
#   plot_coeff_build_up(list(mod_BoostR_1, mod_BoostR_2, mod_BoostR_3, mod_BoostR_4), ls_solutions,
#                       max_iter = 250, title = "BoostR With Different Regularization")
#   
#   plot_coeff_build_up(list(mod_PartBoostR_1, mod_PartBoostR_2, mod_PartBoostR_3, mod_PartBoostR_4),
#                       ls_solutions, max_iter = 250, title = "PartBoostR With Different Regularization")
# }
# dev.off()



####################################################
##### Look at the effective degrees of freedom #####
####################################################
png(paste(save_prefix, "PartBoostR_Effective_Degrees_of_Freedom.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow = c(1,1))
  par(mar = c(3.5,3.5,2,1))
  plot(0:num_iter, type = "n", ylim = c(-0.05, p+0.15),
       col =nice_colors[1], main = "Effective Degrees of Freedom",
       xlab = "", ylab = "")
  title(xlab = "Iteration", ylab = "Effective Degrees of Freedom", line = 2)
  abline(h = p, lty = 2, col = 'gray')
  lines(0:num_iter, mod_BoostR_1$df_hat, lty = 1, col = nice_colors[1], type = "s")
  lines(0:num_iter, mod_BoostR_2$df_hat, lty = 1, col = nice_colors[4], type = "s")
  lines(0:num_iter, mod_BoostR_3$df_hat, lty = 1, col = nice_colors[3], type = "s")
  lines(0:num_iter, mod_BoostR_4$df_hat, lty = 1, col = nice_colors[2], type = "s")
  
  lines(0:num_iter, mod_PartBoostR_1$df_hat, lty = 2, col = nice_colors[1], type = "s")
  lines(0:num_iter, mod_PartBoostR_2$df_hat, lty = 2, col = nice_colors[4], type = "s")
  lines(0:num_iter, mod_PartBoostR_3$df_hat, lty = 2, col = nice_colors[3], type = "s")
  lines(0:num_iter, mod_PartBoostR_4$df_hat, lty = 2, col = nice_colors[2], type = "s")
  ll = legend("bottomright", lty = 1, col = nice_colors[c(1,4,3,2)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "PartBoostR", cex = 0.75)
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1,4,3,2)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  legend(ll$rect$left - kk$rect$w, ll$rect$top - (ll$rect$h - kk$rect$h),
         lty = 2, col = nice_colors[c(1,4,3,2)],
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
         title = "BoostR", cex = 0.75)
}
dev.off()

png(paste(save_prefix, "PartBoostR_Effective_Degrees_of_Freedom_sigma_high.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow = c(1,1))
  par(mar = c(3.5,3.5,2,1))
  plot(0:num_iter, type = "n", ylim = c(-0.05, p+0.15),
       col =nice_colors[1], main = "Effective Degrees of Freedom",
       xlab = "", ylab = "")
  title(xlab = "Iteration", ylab = "Effective Degrees of Freedom", line = 2)
  abline(h = p, lty = 2, col = 'gray')
  lines(0:num_iter, mod_BoostR_high_1$df_hat, lty = 1, col = nice_colors[1], type = "s")
  lines(0:num_iter, mod_BoostR_high_2$df_hat, lty = 1, col = nice_colors[4], type = "s")
  lines(0:num_iter, mod_BoostR_high_3$df_hat, lty = 1, col = nice_colors[3], type = "s")
  lines(0:num_iter, mod_BoostR_high_4$df_hat, lty = 1, col = nice_colors[2], type = "s")
  
  lines(0:num_iter, mod_PartBoostR_high_1$df_hat, lty = 2, col = nice_colors[1], type = "s")
  lines(0:num_iter, mod_PartBoostR_high_2$df_hat, lty = 2, col = nice_colors[4], type = "s")
  lines(0:num_iter, mod_PartBoostR_high_3$df_hat, lty = 2, col = nice_colors[3], type = "s")
  lines(0:num_iter, mod_PartBoostR_high_4$df_hat, lty = 2, col = nice_colors[2], type = "s")
  ll = legend("bottomright", lty = 1, col = nice_colors[c(1,4,3,2)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "PartBoostR", cex = 0.75)
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1,4,3,2)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  legend(ll$rect$left - kk$rect$w, ll$rect$top - (ll$rect$h - kk$rect$h),
         lty = 2, col = nice_colors[c(1,4,3,2)],
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
         title = "BoostR", cex = 0.75)
}
dev.off()


############################################################
##### NOW WE WANT TO LOOK AT THE TEST MSE OF OUR MODEL #####
############################################################
# Calculate the test MSE for our models
{
  # sigma_low
  mod_BoostR_1_pred = PartBoostR_predict(mod_BoostR_1, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)
  mod_BoostR_2_pred = PartBoostR_predict(mod_BoostR_2, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)
  mod_BoostR_3_pred = PartBoostR_predict(mod_BoostR_3, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)
  mod_BoostR_4_pred = PartBoostR_predict(mod_BoostR_4, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)
  
  mod_PartBoostR_1_pred = PartBoostR_predict(mod_PartBoostR_1, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)
  mod_PartBoostR_2_pred = PartBoostR_predict(mod_PartBoostR_2, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)
  mod_PartBoostR_3_pred = PartBoostR_predict(mod_PartBoostR_3, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)
  mod_PartBoostR_4_pred = PartBoostR_predict(mod_PartBoostR_4, X_test = X_test_sc, y_test = y_test_low, after_iter = NULL)

  # sigma_high
  mod_BoostR_high_1_pred = PartBoostR_predict(mod_BoostR_high_1, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
  mod_BoostR_high_2_pred = PartBoostR_predict(mod_BoostR_high_2, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
  mod_BoostR_high_3_pred = PartBoostR_predict(mod_BoostR_high_3, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
  mod_BoostR_high_4_pred = PartBoostR_predict(mod_BoostR_high_4, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
  
  mod_PartBoostR_high_1_pred = PartBoostR_predict(mod_PartBoostR_high_1, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
  mod_PartBoostR_high_2_pred = PartBoostR_predict(mod_PartBoostR_high_2, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
  mod_PartBoostR_high_3_pred = PartBoostR_predict(mod_PartBoostR_high_3, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
  mod_PartBoostR_high_4_pred = PartBoostR_predict(mod_PartBoostR_high_4, X_test = X_test_sc, y_test = y_test_high, after_iter = NULL)
}

{ 
  # Look at the smallest training MSE
  aux1 = range(mod_BoostR_1$mse_training)
  aux2 = range(mod_BoostR_2$mse_training)
  aux3 = range(mod_BoostR_3$mse_training)
  aux4 = range(mod_BoostR_4$mse_training)
  
  aux5 = range(mod_PartBoostR_1$mse_training)
  aux6 = range(mod_PartBoostR_2$mse_training)
  aux7 = range(mod_PartBoostR_3$mse_training)
  aux8 = range(mod_PartBoostR_4$mse_training)
  aux9 = lm_mse_train
  range_train = range(aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9)
  range_train_length = range_train[2] - range_train[1]
  
  # Look at the smallest test MSE
  aux1 = range(mod_BoostR_1_pred$mse_loss)
  aux2 = range(mod_BoostR_2_pred$mse_loss)
  aux3 = range(mod_BoostR_3_pred$mse_loss)
  aux4 = range(mod_BoostR_4_pred$mse_loss)
  
  aux5 = range(mod_PartBoostR_1_pred$mse_loss)
  aux3 = range(mod_PartBoostR_2_pred$mse_loss)
  aux7 = range(mod_PartBoostR_3_pred$mse_loss)
  aux8 = range(mod_PartBoostR_4_pred$mse_loss)
  aux9 = lm_mse_test
  range_test =range(aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9)
  range_test_length = range_test[2] - range_test[1]
  
  min_y = range(range_train, range_test)[1]
  max_y = range(range_train, range_test)[2]
}
ylims = c(range_train[1] - 0.1*range_train_length,
          range_train[2] + 0.1*range_train_length)
ylims2 = ylims
#ylims = c(23,29)
#ylims2 = c(24.8, 28.5)

# Plot the MSE development for each iter
png(paste(save_prefix, "PartBoostR_Train_Test_MSE.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))
  # Training MSE PLOT
  plot(1, type="n", main = "Train MSE for PartBoostR",
       xlim = c(0, num_iter), ylim = ylims,
       xlab = "", ylab = "")
  title(xlab = "Iterations", ylab = "MSE", line = 2.3)
  lines(0:num_iter, mod_BoostR_2$mse_training, col = nice_colors[1], lty = 1)
  lines(0:num_iter, mod_PartBoostR_1$mse_training, col = nice_colors[4], lty = 1)
  lines(0:num_iter, mod_PartBoostR_2$mse_training, col = nice_colors[3], lty = 1)
  lines(0:num_iter, mod_PartBoostR_3$mse_training, col = nice_colors[2], lty = 1)
  abline(h = lm_mse_train, col = 'darkgray', lty = 2)
  ll = legend("topright", lty = 1, col = nice_colors[c(4,3,2)], cex = 0.7,
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " "))),
              title = "PartBoostR")
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.7)
  
  legend(ll$rect$left + (ll$rect$w - kk$rect$w),
         ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(1)],
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
         title = "BoostR", cex = 0.7)
  
  # TEST MSE PLOT
  plot(1, type="n", main = "Test MSE for PartBoostR",
       xlim = c(0, num_iter), ylim = ylims2,
       xlab = "", ylab = "")
  title(xlab = "Iterations", ylab = "MSE", line = 2.3)
  lines(0:num_iter, mod_BoostR_2_pred$mse_loss, col = nice_colors[1], lty = 1)
  lines(0:num_iter, mod_PartBoostR_1_pred$mse_loss, col = nice_colors[4], lty = 1)
  lines(0:num_iter, mod_PartBoostR_2_pred$mse_loss, col = nice_colors[3], lty = 1)
  lines(0:num_iter, mod_PartBoostR_3_pred$mse_loss, col = nice_colors[2], lty = 1)
  abline(h = lm_mse_test, col = 'darkgray', lty = 2)
  ll = legend("topright", lty = 1, col = nice_colors[c(4,3,2)], cex = 0.75,
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " "))),
              title = "PartBoostR")
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  legend(ll$rect$left + (ll$rect$w - kk$rect$w),
         ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(1)],
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
         title = "BoostR", cex = 0.75)
}
dev.off()

##################################################
##### Create Regular Ridge regression models #####
##################################################
Ridge_regression = function(X_train, y_train, lambda) {
  # Function that fits a ridge model
  # Returns the ridge solution, df, and mse training error.
  
  B_ridge = solve(t(X_train) %*% X_train + lambda * diag(ncol(X_train))) %*% t(X_train)
  beta_ridge = B_ridge %*% y_train
  df = sum(diag(X_train %*% B_ridge))
  pred_val = X_train %*% beta_ridge
  mse_loss = mean((y_train - pred_val)^2)
  
  return(list(beta = beta_ridge, train_mse = mse_loss, df = df))
}

# Some parameters and variables to store values
{
  num_ridge_mod = 5000
  lambda_values_iter = seq(0, 5000, length.out = num_ridge_mod)
  
  # sigma_low
  df_obtained = rep(NA, num_ridge_mod)
  train_err = rep(NA, num_ridge_mod)
  test_err = rep(NA, num_ridge_mod)
  beta_ridge_mat = matrix(NA, ncol = p, nrow = num_ridge_mod)
  
  # sigma_high
  df_obtained_high = rep(NA, num_ridge_mod)
  train_err_high = rep(NA, num_ridge_mod)
  test_err_high = rep(NA, num_ridge_mod)
  beta_ridge_mat_high = matrix(NA, ncol = p, nrow = num_ridge_mod)
}

# Iterate over all the different penalities for ridge regression
for (idx in 1:num_ridge_mod) {
  ### sigma_low
  temp_ridge = Ridge_regression(X_train_sc, y_train_low, lambda_values_iter[idx])
  df_obtained[idx] = temp_ridge$df
  train_err[idx] = temp_ridge$train_mse
  beta_ridge_mat[idx, ] = temp_ridge$beta
  
  # Calculate test error
  test_err[idx] = mean((X_test_sc %*% temp_ridge$beta - y_test_low)^2)
  
  ### sigma_high
  temp_ridge_high = Ridge_regression(X_train_sc, y_train_high, lambda_values_iter[idx])
  df_obtained_high[idx] = temp_ridge_high$df
  train_err_high[idx] = temp_ridge_high$train_mse
  beta_ridge_mat_high[idx, ] = temp_ridge_high$beta

  # Calculate test error
  test_err_high[idx] = mean((X_test_sc %*% temp_ridge_high$beta - y_test_high)^2)
}


# Plot coef build up for ridge and parboost5000
{
  min_beta = min(ls_solutions)
  max_beta = max(ls_solutions)
  range_beta = max_beta - min_beta 
  min_beta = min_beta - 0.05 * range_beta
  max_beta = max_beta + 0.35 * range_beta
}

###############################################################################
##### Look at the effective degrees of freedom and Coef build up together #####
###############################################################################
png(paste(save_prefix, "Coef_Build_Up_PartBoostR_EDF.png", sep = ""), width = 3250, height = 1500, res = 350)
{ 
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))
 
  # Plot the effective degrees of freedom
  plot(0:num_iter, type = "n", ylim = c(-0.05, p+0.15),
       col =nice_colors[1], main = "Effective Degrees of Freedom",
       xlab = "", ylab = "")
  title(xlab = "Iteration", ylab = "Effective Degrees of Freedom", line = 2)
  abline(h = p, lty = 2, col = 'gray')
  lines(0:num_iter, mod_BoostR_1$df_hat, lty = 1, col = nice_colors[1], type = "s")
  lines(0:num_iter, mod_BoostR_2$df_hat, lty = 1, col = nice_colors[4], type = "s")
  lines(0:num_iter, mod_BoostR_3$df_hat, lty = 1, col = nice_colors[3], type = "s")
  lines(0:num_iter, mod_BoostR_4$df_hat, lty = 1, col = nice_colors[2], type = "s")
  
  lines(0:num_iter, mod_PartBoostR_1$df_hat, lty = 2, col = nice_colors[1], type = "s")
  lines(0:num_iter, mod_PartBoostR_2$df_hat, lty = 2, col = nice_colors[4], type = "s")
  lines(0:num_iter, mod_PartBoostR_3$df_hat, lty = 2, col = nice_colors[3], type = "s")
  lines(0:num_iter, mod_PartBoostR_4$df_hat, lty = 2, col = nice_colors[2], type = "s")
  ll = legend("bottomright", lty = 1, col = nice_colors[c(1,4,3,2)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "PartBoostR", cex = 0.75)
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1,4,3,2)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  legend(ll$rect$left - kk$rect$w, ll$rect$top - (ll$rect$h - kk$rect$h),
         lty = 2, col = nice_colors[c(1,4,3,2)],
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
         title = "BoostR", cex = 0.75)
  
  
  # Plot of coeff build up With ridge and PartBoostR
  plot(NA, type = 'n',
       main = "Coefficient Build-Up", xlim = c(0,p),
       ylim = c(min_beta, max_beta), col = nice_colors[1], xlab = "", ylab = "")
  title(ylab=bquote("Standardized" ~ beta['j']), xlab = "Effective Degrees of Freedom", line = 2)
  
  lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,1], col = nice_colors[1], type = 's')
  lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,2], col = nice_colors[2], type = 's')
  lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,3], col = nice_colors[3], type = 's')
  lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,4], col = nice_colors[4], type = 's')
  
  lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,1], col = nice_colors[1],lty = 2, type = 's')
  lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,2], col = nice_colors[2],lty = 2, type = 's')
  lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,3], col = nice_colors[3],lty = 2, type = 's')
  lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,4], col = nice_colors[4],lty = 2, type = 's')
  
  lines(df_obtained, beta_ridge_mat[,1], col = nice_colors[1], lty = 3, type = 's')
  lines(df_obtained, beta_ridge_mat[,2], col = nice_colors[2], lty = 3, type = 's')
  lines(df_obtained, beta_ridge_mat[,3], col = nice_colors[3], lty = 3, type = 's')
  lines(df_obtained, beta_ridge_mat[,4], col = nice_colors[4], lty = 3, type = 's')
  points(rep(p, p), ls_solutions)
  legend("topleft", legend = c(as.expression(bquote(beta[1] ~ " ")),
                               as.expression(bquote(beta[2] ~ " ")),
                               as.expression(bquote(beta[3] ~ " ")),
                               as.expression(bquote(beta[4] ~ " "))),
         lty = 1, col = nice_colors, cex = 0.75)
  temp = legend("topright", legend = c("PartBoostR", "BoostR", "Ridge"),
                lty = c(1,2,3), col = 'black', cex = 0.75)
}
dev.off()


if (FALSE) {
  png(paste(save_prefix, "PartBoostR_Ridge.png", sep = ""), width = 3250, height = 1500, res = 350)
  {
    par(mfrow = c(1,2))
    par(mar = c(3.5,3.5,2,1))
    plot(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,1], type = 's',
         main = "Coeff. Build Up For PartBoostR and Ridge", xlim = c(0,p),
         ylim = c(min_beta, max_beta), col = nice_colors[1], xlab = "", ylab = "")
    title(ylab=bquote("Standardized" ~ beta['j']), xlab = "Effective Degrees of Freedom", line = 2)
    
    lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,2], col = nice_colors[2], type = 's')
    lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,3], col = nice_colors[3], type = 's')
    lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2$beta_hat[,4], col = nice_colors[4], type = 's')
    
    lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,1], col = nice_colors[1],lty = 2, type = 's')
    lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,2], col = nice_colors[2],lty = 2, type = 's')
    lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,3], col = nice_colors[3],lty = 2, type = 's')
    lines(mod_BoostR_2$df_hat, mod_BoostR_2$beta_hat[,4], col = nice_colors[4],lty = 2, type = 's')
    
    lines(df_obtained, beta_ridge_mat[,1], col = nice_colors[1], lty = 3, type = 's')
    lines(df_obtained, beta_ridge_mat[,2], col = nice_colors[2], lty = 3, type = 's')
    lines(df_obtained, beta_ridge_mat[,3], col = nice_colors[3], lty = 3, type = 's')
    lines(df_obtained, beta_ridge_mat[,4], col = nice_colors[4], lty = 3, type = 's')
    points(rep(p, p), ls_solutions)
    legend("topleft", legend = c(as.expression(bquote(beta[1] ~ " ")),
                                 as.expression(bquote(beta[2] ~ " ")),
                                 as.expression(bquote(beta[3] ~ " ")),
                                 as.expression(bquote(beta[4] ~ " "))),
           lty = 1, col = nice_colors, cex = 0.7)
    temp = legend("topright", legend = c("PartBoostR", "BoostR", "Ridge"),
                  lty = c(1,2,3), col = 'black', cex = 0.7)
    
    # Plot the MSE test with respect to degrees of freedom
    plot(NA, type = 'n', col = nice_colors[1],
         xlab = "", ylab = "", main = "Test MSE for PartBoostR and Ridge", 
         ylim = ylims2, xlim =c(0, p))
    title(xlab = "Effective Degrees of Freedom", ylab = "MSE", line = 2.3)
    
    #lines(mod_BoostR_1$df_hat, mod_BoostR_1_pred$mse_loss, col = nice_colors[1], type = 's')
    lines(mod_BoostR_2$df_hat, mod_BoostR_2_pred$mse_loss, col = nice_colors[1], type = 's')
    #lines(mod_BoostR_3$df_hat, mod_BoostR_3_pred$mse_loss, col = nice_colors[1], type = 's')
    #lines(mod_BoostR_4$df_hat, mod_BoostR_4_pred$mse_loss, col = nice_colors[1], type = 's')
    
    lines(mod_PartBoostR_1$df_hat, mod_PartBoostR_1_pred$mse_loss, col = nice_colors[4], type = 's')
    lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2_pred$mse_loss, col = nice_colors[5], type = 's')
    lines(mod_PartBoostR_3$df_hat, mod_PartBoostR_3_pred$mse_loss, col = nice_colors[2], type = 's')
    lines(mod_PartBoostR_4$df_hat, mod_PartBoostR_3_pred$mse_loss, col = nice_colors[6], type = 's')
    
    
    lines(df_obtained, test_err, col = nice_colors[3], type  = 's')
    
    abline(h = lm_mse_test, lty = 2, col = 'darkgray')
    ll = legend("topright", lty = 1, col = nice_colors[c(4,5,2,6)], cex = 0.7,
                legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                           as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                           as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                           as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
                title = "PartBoostR")
    kk = legend(-100, -100, lty = 1, col = nice_colors[c(1)],
                legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
                title = "BoostR", cex = 0.7)
    
    ll = legend(ll$rect$left + (ll$rect$w - kk$rect$w),
                ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(1)],
                legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
                title = "BoostR", cex = 0.7)
    
    kk = legend(-100, -100, lty = 1, col = nice_colors[c(3)],
                legend = c('Ridge'), cex = 0.7)
    
    legend(ll$rect$left + (ll$rect$w - kk$rect$w),
           ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(3)],
           legend = c('Ridge'), cex = 0.7)
  }
  dev.off()
}

# PLot with test MSE for sigma_high
png(paste(save_prefix, "Test_MSE_sigma_high.png", sep = ""), width = 3250, height = 1500, res = 350)
{ 
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))
  
  # First plot is test MSE of BoostR, PartBoostR wrt num iter
  plot(1, type="n", main = "Test Mean Squared Error",
       xlim = c(0, num_iter), ylim = c(25,30),
       xlab = "", ylab = "")
  title(xlab = "Iterations", ylab = "MSE", line = 2.3)
  
  lines(0:num_iter, mod_BoostR_high_2_pred$mse_loss, col = nice_colors[1], lty = 1, type = 's')
  
  lines(0:num_iter, mod_PartBoostR_high_1_pred$mse_loss, col = nice_colors[4], lty = 1, type = 's')
  lines(0:num_iter, mod_PartBoostR_high_2_pred$mse_loss, col = nice_colors[6], lty = 1, type = 's')
  lines(0:num_iter, mod_PartBoostR_high_3_pred$mse_loss, col = nice_colors[2], lty = 1, type = 's')
  lines(0:num_iter, mod_PartBoostR_high_4_pred$mse_loss, col = nice_colors[5], lty = 1, type = 's')
  
  abline(h = lm_mse_test_high, col = 'darkgray', lty = 2)
  ll = legend("topright", lty = 1, col = nice_colors[c(4,6,2,5)], cex = 0.75,
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "PartBoostR")
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  legend(ll$rect$left + (ll$rect$w - kk$rect$w),
         ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(1)],
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
         title = "BoostR", cex = 0.75)
  
  
  
  # Second plot is test MSE of BoostR, PartBoostR, ridge wrt degrees of freedom
  # Plot the MSE test with respect to degrees of freedom
  plot(NA, type = 'n', col = nice_colors[1],
       xlab = "", ylab = "", main = "Test Mean Squared Error", 
       ylim = c(25,30), xlim =c(0, p))
  title(xlab = "Effective Degrees of Freedom", ylab = "MSE", line = 2.3)
  
  #lines(mod_BoostR_1$df_hat, mod_BoostR_1_pred$mse_loss, col = nice_colors[1], type = 's')
  lines(mod_BoostR_high_2$df_hat, mod_BoostR_high_2_pred$mse_loss, col = nice_colors[1], type = 's')
  #lines(mod_BoostR_3$df_hat, mod_BoostR_3_pred$mse_loss, col = nice_colors[1], type = 's')
  #lines(mod_BoostR_4$df_hat, mod_BoostR_4_pred$mse_loss, col = nice_colors[1], type = 's')
  
  lines(mod_PartBoostR_high_1$df_hat, mod_PartBoostR_high_1_pred$mse_loss, col = nice_colors[4], type = 's')
  lines(mod_PartBoostR_high_2$df_hat, mod_PartBoostR_high_2_pred$mse_loss, col = nice_colors[6], type = 's')
  lines(mod_PartBoostR_high_3$df_hat, mod_PartBoostR_high_3_pred$mse_loss, col = nice_colors[2], type = 's')
  lines(mod_PartBoostR_high_4$df_hat, mod_PartBoostR_high_3_pred$mse_loss, col = nice_colors[5], type = 's')
  
  
  lines(df_obtained, test_err_high, col = nice_colors[3], type  = 's')
  
  abline(h = lm_mse_test_high, lty = 2, col = 'darkgray')
  ll = legend("topright", lty = 1, col = nice_colors[c(4,6,2,5)], cex = 0.75,
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "PartBoostR")
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  ll = legend(ll$rect$left + (ll$rect$w - kk$rect$w),
              ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(3)],
              legend = c('Ridge'), cex = 0.75)
  
  legend(ll$rect$left + (ll$rect$w - kk$rect$w),
         ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(3)],
         legend = c('Ridge'), cex = 0.75)
}
dev.off()



ylims = c(0,6)
# PLot with test MSE for sigma_low
png(paste(save_prefix, "Test_MSE_sigma_low.png", sep = ""), width = 3250, height = 1500, res = 350)
{ 
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))
  
  # First plot is test MSE of BoostR, PartBoostR wrt num iter
  plot(1, type="n", main = "Test Mean Squared Error",
       xlim = c(0, num_iter), ylim = ylims,
       xlab = "", ylab = "")
  title(xlab = "Iterations", ylab = "MSE", line = 2.3)
  
  lines(0:num_iter, mod_BoostR_2_pred$mse_loss, col = nice_colors[1], lty = 1, type = 's')
  
  lines(0:num_iter, mod_PartBoostR_1_pred$mse_loss, col = nice_colors[4], lty = 1, type = 's')
  lines(0:num_iter, mod_PartBoostR_2_pred$mse_loss, col = nice_colors[6], lty = 1, type = 's')
  lines(0:num_iter, mod_PartBoostR_3_pred$mse_loss, col = nice_colors[2], lty = 1, type = 's')
  lines(0:num_iter, mod_PartBoostR_4_pred$mse_loss, col = nice_colors[5], lty = 1, type = 's')
  
  abline(h = lm_mse_test, col = 'darkgray', lty = 2)
  ll = legend("topright", lty = 1, col = nice_colors[c(4,6,2,5)], cex = 0.75,
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "PartBoostR")
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  legend(ll$rect$left - kk$rect$w,
         ll$rect$top, lty = 1, col = nice_colors[c(1)],
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
         title = "BoostR", cex = 0.75)
  
  
  
  # Second plot is test MSE of BoostR, PartBoostR, ridge wrt degrees of freedom
  # Plot the MSE test with respect to degrees of freedom
  plot(NA, type = 'n', col = nice_colors[1],
       xlab = "", ylab = "", main = "Test Mean Squared Error", 
       ylim = ylims, xlim =c(0, p))
  title(xlab = "Effective Degrees of Freedom", ylab = "MSE", line = 2.3)
  
  #lines(mod_BoostR_1$df_hat, mod_BoostR_1_pred$mse_loss, col = nice_colors[1], type = 's')
  lines(mod_BoostR_2$df_hat, mod_BoostR_2_pred$mse_loss, col = nice_colors[1], type = 's')
  #lines(mod_BoostR_3$df_hat, mod_BoostR_3_pred$mse_loss, col = nice_colors[1], type = 's')
  #lines(mod_BoostR_4$df_hat, mod_BoostR_4_pred$mse_loss, col = nice_colors[1], type = 's')
  
  lines(mod_PartBoostR_1$df_hat, mod_PartBoostR_1_pred$mse_loss, col = nice_colors[4], type = 's')
  lines(mod_PartBoostR_2$df_hat, mod_PartBoostR_2_pred$mse_loss, col = nice_colors[6], type = 's')
  lines(mod_PartBoostR_3$df_hat, mod_PartBoostR_3_pred$mse_loss, col = nice_colors[2], type = 's')
  lines(mod_PartBoostR_4$df_hat, mod_PartBoostR_3_pred$mse_loss, col = nice_colors[5], type = 's')
  
  
  lines(df_obtained, test_err, col = nice_colors[3], type  = 's')
  
  abline(h = lm_mse_test, lty = 2, col = 'darkgray')
  ll = legend("topright", lty = 1, col = nice_colors[c(4,6,2,5)], cex = 0.75,
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[4]) ~ " "))),
              title = "PartBoostR")
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  ll = legend(ll$rect$left - kk$rect$w,
              ll$rect$top, lty = 1, col = nice_colors[c(1)],
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  kk = legend(-100, -100, lty = 1, col = nice_colors[c(3)],
              legend = c('Ridge'), cex = 0.75)
  
  legend(ll$rect$left + (ll$rect$w - kk$rect$w),
         ll$rect$top - ll$rect$h, lty = 1, col = nice_colors[c(3)],
         legend = c('Ridge'), cex = 0.75)
}
dev.off()




#############################################################
##### LetS look at the GenPartBoostR and Fisher Scoring #####
#############################################################
# Find the least squares solutions of the unormalzied covariates
lm_model_inter = lm(response~., data = training)
lm_model_inter$coefficients
ls_inter_sol = lm_model_inter$coefficients

lm_model_inter_high = lm(response~., data = training_high)
lm_model_inter_high$coefficients
ls_inter_sol_high = lm_model_inter_high$coefficients

source_path = "C:/Users/lars9/OneDrive/Master/R_codes/Chapter 2 - Ridge regression boosting/"
source(paste(source_path, "GenPartBoostR_Source_Code.R", sep=""))
source(paste(source_path, "FisherScoring_Normal_Data_Source_Code.R", sep=""))


# Fit the GenPartBoostR models
num_iter = 2000
lambda = 5000
lambda_values = c(1000, 5000, 10000, 25000)
# sigma_low
{
  # GenBoostR
  mod_gbr = GenPartBoostR_normal(X_train, y_train_low, rep(0, p+1),
                                 sigma_squared = 1, lambda = lambda_values[2],
                                 num_iter = num_iter, num_iter_init = 10,
                                 V_compulsory = 0:p, TRUE)
  # GenPartBoostR
  mod_gpbr = GenPartBoostR_normal(X_train, y_train_low, rep(0, p+1),
                                  sigma_squared = 1, lambda = lambda_values[2],
                                  num_iter = num_iter, num_iter_init = 10,
                                  V_compulsory = c(), TRUE)
  # FisherScoring
  mod_fish = FisherScoring_normal(X = X_train,
                                  y = y_train_low,
                                  beta_hat = rep(0, p+1),
                                  sigma_squared = 1,
                                  lambda = 0,
                                  num_iter = num_iter,
                                  print = TRUE)
  
}
# sigma_high
{
  # GenBoostR
  mod_gbr_high = GenPartBoostR_normal(X_train, y_train_high, rep(0, p+1),
                                 sigma_squared = 1, lambda = lambda_values[2],
                                 num_iter = num_iter, num_iter_init = 10,
                                 V_compulsory = 0:p, TRUE)
  mod_gbr_high_1 = GenPartBoostR_normal(X_train, y_train_high, rep(0, p+1),
                                      sigma_squared = 1, lambda = lambda_values[1],
                                      num_iter = num_iter, num_iter_init = 10,
                                      V_compulsory = 0:p, TRUE)
  mod_gbr_high_2 = GenPartBoostR_normal(X_train, y_train_high, rep(0, p+1),
                                      sigma_squared = 1, lambda = lambda_values[3],
                                      num_iter = num_iter, num_iter_init = 10,
                                      V_compulsory = 0:p, TRUE)
  
  # GenPartBoostR
  mod_gpbr_high = GenPartBoostR_normal(X_train, y_train_high, rep(0, p+1),
                                  sigma_squared = 1, lambda = lambda,
                                  num_iter = num_iter, num_iter_init = 10,
                                  V_compulsory = c(), TRUE)
  
  mod_gpbr_high_1 = GenPartBoostR_normal(X_train, y_train_high, rep(0, p+1),
                                       sigma_squared = 1, lambda = lambda_values[1],
                                       num_iter = num_iter, num_iter_init = 10,
                                       V_compulsory = c(), TRUE)
  
  mod_gpbr_high_2 = GenPartBoostR_normal(X_train, y_train_high, rep(0, p+1),
                                       sigma_squared = 1, lambda = lambda_values[3],
                                       num_iter = num_iter, num_iter_init = 10,
                                       V_compulsory = c(), TRUE)
  
  # FisherScoring
  mod_fish_high = FisherScoring_normal(X = X_train,
                                  y = y_train_high,
                                  beta_hat = rep(0, p+1),
                                  sigma_squared = 1,
                                  lambda = 0,
                                  num_iter = num_iter,
                                  print = TRUE)
  
}

# See how many updates each parameters have recieved
table(mod_gpbr$update_trace)
table(mod_PartBoostR_2$update_trace)
table(mod_PartBoostR_3$update_trace)
table(mod_PartBoostR_4$update_trace)
# And in which iteration was the first time the param got updated
match(unique(mod_gpbr$update_trace), mod_gpbr$update_trace)
match(unique(mod_PartBoostR_2$update_trace), mod_PartBoostR_2$update_trace)
match(unique(mod_PartBoostR_3$update_trace), mod_PartBoostR_3$update_trace)
match(unique(mod_PartBoostR_4$update_trace), mod_PartBoostR_4$update_trace)


# Make a plot of their coefficient build ups
png(paste(save_prefix, "GenPartBoostR_coeff_build_up.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))
  
  # Take a closer look at the first 250 iterations 
  max_iter = 250
  plot(NA, type = 'n',
       main = "Coeffcient Build-Up", xlim = c(0,max_iter),
       ylim = c(-1.05, 5), col = nice_colors[5], xlab = "", ylab = "")
  title(ylab=bquote(beta['j']), xlab = "Iterations", line = 2)
  
  # The GenBoostR
  lines(0:max_iter, mod_gbr$beta_hat_matrix[1,1:(max_iter+1)], col = nice_colors[5], type = 's', lty = 1)
  lines(0:max_iter, mod_gbr$beta_hat_matrix[2,1:(max_iter+1)], col = nice_colors[1], type = 's', lty = 1)
  lines(0:max_iter, mod_gbr$beta_hat_matrix[3,1:(max_iter+1)], col = nice_colors[2], type = 's', lty = 1)
  lines(0:max_iter, mod_gbr$beta_hat_matrix[4,1:(max_iter+1)], col = nice_colors[3], type = 's', lty = 1)
  lines(0:max_iter, mod_gbr$beta_hat_matrix[5,1:(max_iter+1)], col = nice_colors[4], type = 's', lty = 1)
  
  # The GenPartBoostR
  lines(0:max_iter, mod_gpbr$beta_hat_matrix[1,1:(max_iter+1)], col = nice_colors[5], type = 's', lty = 2)
  lines(0:max_iter, mod_gpbr$beta_hat_matrix[2,1:(max_iter+1)], col = nice_colors[1], type = 's', lty = 2)
  lines(0:max_iter, mod_gpbr$beta_hat_matrix[3,1:(max_iter+1)], col = nice_colors[2], type = 's', lty = 2)
  lines(0:max_iter, mod_gpbr$beta_hat_matrix[4,1:(max_iter+1)], col = nice_colors[3], type = 's', lty = 2)
  lines(0:max_iter, mod_gpbr$beta_hat_matrix[5,1:(max_iter+1)], col = nice_colors[4], type = 's', lty = 2)
  
  # Fisher Scoring 
  lines(0:max_iter, mod_fish$beta_hat_matrix[1,1:(max_iter+1)], col = nice_colors[5], type = 's', lty = 3)
  lines(0:max_iter, mod_fish$beta_hat_matrix[2,1:(max_iter+1)], col = nice_colors[1], type = 's', lty = 3)
  lines(0:max_iter, mod_fish$beta_hat_matrix[3,1:(max_iter+1)], col = nice_colors[2], type = 's', lty = 3)
  lines(0:max_iter, mod_fish$beta_hat_matrix[4,1:(max_iter+1)], col = nice_colors[3], type = 's', lty = 3)
  lines(0:max_iter, mod_fish$beta_hat_matrix[5,1:(max_iter+1)], col = nice_colors[4], type = 's', lty = 3)
  
  # Least squares solution
  points(rep(max_iter, p+1), ls_inter_sol)
  
  # Add legend
  legend("topleft", legend = c(as.expression(bquote(beta[0] ~ " ")),
                               as.expression(bquote(beta[1] ~ " ")),
                               as.expression(bquote(beta[2] ~ " ")),
                               as.expression(bquote(beta[3] ~ " ")),
                               as.expression(bquote(beta[4] ~ " "))),
         lty = 1, col = nice_colors[c(5,1,2,3,4)], cex = 0.7)
  temp = legend("topright", legend = c("Fisher", "GenBoostR", "PartGenBoostR"),
                lty = c(3,1,2), col = 'black', cex = 0.7)

  
  plot(NA, type = 'n',
       main = "Coefficient Build-Up", xlim = c(0,num_iter),
       ylim = c(-1.05, 5), col = nice_colors[5], xlab = "", ylab = "")
  title(ylab=bquote(beta['j']), xlab = "Iterations", line = 2)
  
  # The GenBoostR
  lines(0:num_iter, mod_gbr$beta_hat_matrix[1,], col = nice_colors[5], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[2,], col = nice_colors[1], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[3,], col = nice_colors[2], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[4,], col = nice_colors[3], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[5,], col = nice_colors[4], type = 's', lty = 1)
  
  # The GenPartBoostR
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[1,], col = nice_colors[5], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[2,], col = nice_colors[1], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[3,], col = nice_colors[2], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[4,], col = nice_colors[3], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[5,], col = nice_colors[4], type = 's', lty = 2)
  
  # Fisher Scoring 
  lines(0:num_iter, mod_fish$beta_hat_matrix[1,], col = nice_colors[5], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[2,], col = nice_colors[1], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[3,], col = nice_colors[2], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[4,], col = nice_colors[3], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[5,], col = nice_colors[4], type = 's', lty = 3)
  
  # Least squares solution
  points(rep(num_iter, p+1), ls_inter_sol)
  
  # Add legend
  legend("topleft", legend = c(as.expression(bquote(beta[0] ~ " ")),
                               as.expression(bquote(beta[1] ~ " ")),
                               as.expression(bquote(beta[2] ~ " ")),
                               as.expression(bquote(beta[3] ~ " ")),
                               as.expression(bquote(beta[4] ~ " "))),
         lty = 1, col = nice_colors[c(5,1,2,3,4)], cex = 0.7)
  temp = legend("topright", legend = c("Fisher", "GenBoostR", "PartGenBoostR"),
                lty = c(3,1,2), col = 'black', cex = 0.7)
  
}
dev.off()

######################################
##### LOOK AT THE TRAIN TEST MSE #####
######################################
{
  Z_train = cbind(rep(1, n_train), X_train)
  Z_test = cbind(rep(1, n_test), X_test)
  
  # sigma_low
  ls_inter_train_err = mean((y_train_low - Z_train %*% ls_inter_sol)^2)
  ls_inter_test_err = mean((y_test_low - Z_test %*% ls_inter_sol)^2)
  
  mod_gbr_train_err = rep(NA, num_iter+1)
  mod_gbr_test_err = rep(NA, num_iter+1)
  
  mod_gpbr_train_err = rep(NA, num_iter+1)
  mod_gpbr_test_err = rep(NA, num_iter+1)
  
  mod_fish_train_err = rep(NA, num_iter+1)
  mod_fish_test_err = rep(NA, num_iter+1)
  
  # sigma_high
  ls_inter_train_err_high = mean((y_train_high - Z_train %*% ls_inter_sol_high)^2)
  ls_inter_test_err_high  = mean((y_test_high - Z_test %*% ls_inter_sol_high)^2)
  
  mod_gbr_train_err_high = rep(NA, num_iter+1)
  mod_gbr_test_err_high  = rep(NA, num_iter+1)
  
  mod_gpbr_train_err_high = rep(NA, num_iter+1)
  mod_gpbr_test_err_high  = rep(NA, num_iter+1)
  
  mod_gbr_1_test_err_high  = rep(NA, num_iter+1)
  mod_gbr_2_test_err_high  = rep(NA, num_iter+1)
  
  mod_gpbr_1_test_err_high = rep(NA, num_iter+1)
  mod_gpbr_2_test_err_high = rep(NA, num_iter+1)
  
  mod_fish_train_err_high = rep(NA, num_iter+1)
  mod_fish_test_err_high  = rep(NA, num_iter+1)
}

for (m in 1:(num_iter+1)) {
  ### sigma_low
  # GenBoostR
  temp_beta = mod_gbr$beta_hat_matrix[,m]
  pred_val_train = Z_train %*% temp_beta
  mod_gbr_train_err[m] = mean((pred_val_train - y_train_low)^2)
  pred_val_test = Z_test %*% temp_beta 
  mod_gbr_test_err[m] = mean((pred_val_test - y_test_low)^2)
  
  # GenPartBoostR
  temp_beta = mod_gpbr$beta_hat_matrix[,m]
  pred_val_train = Z_train %*% temp_beta
  mod_gpbr_train_err[m] = mean((pred_val_train - y_train_low)^2)
  pred_val_test = Z_test %*% temp_beta 
  mod_gpbr_test_err[m] = mean((pred_val_test - y_test_low)^2)
  
  # Fisher
  temp_beta = mod_fish$beta_hat_matrix[,m] 
  pred_val_train = Z_train %*% temp_beta
  mod_fish_train_err[m] = mean((pred_val_train - y_train_low)^2)
  pred_val_test = Z_test %*% temp_beta 
  mod_fish_test_err[m] = mean((pred_val_test - y_test_low)^2)
  
  ### sigma_high
  # GenBoostR
  temp_beta_high = mod_gbr_high$beta_hat_matrix[,m]
  pred_val_train_high = Z_train %*% temp_beta_high
  mod_gbr_train_err_high[m] = mean((pred_val_train_high - y_train_high)^2)
  pred_val_test_high = Z_test %*% temp_beta_high
  mod_gbr_test_err_high[m] = mean((pred_val_test_high - y_test_high)^2)
  
  temp_beta_high = mod_gbr_high_1$beta_hat_matrix[,m]
  pred_val_test_high = Z_test %*% temp_beta_high
  mod_gbr_1_test_err_high[m] = mean((pred_val_test_high - y_test_high)^2)
  
  temp_beta_high = mod_gbr_high_2$beta_hat_matrix[,m]
  pred_val_test_high = Z_test %*% temp_beta_high
  mod_gbr_2_test_err_high[m] = mean((pred_val_test_high - y_test_high)^2)
  
  # GenPartBoostR
  temp_beta_high = mod_gpbr_high$beta_hat_matrix[,m]
  pred_val_train_high = Z_train %*% temp_beta_high
  mod_gpbr_train_err_high[m] = mean((pred_val_train_high - y_train_high)^2)
  pred_val_test_high = Z_test %*% temp_beta_high 
  mod_gpbr_test_err_high[m] = mean((pred_val_test_high - y_test_high)^2)
  
  temp_beta_high = mod_gpbr_high_1$beta_hat_matrix[,m]
  pred_val_test_high = Z_test %*% temp_beta_high 
  mod_gpbr_1_test_err_high[m] = mean((pred_val_test_high - y_test_high)^2)
  
  temp_beta_high = mod_gpbr_high_2$beta_hat_matrix[,m]
  pred_val_test_high = Z_test %*% temp_beta_high 
  mod_gpbr_2_test_err_high[m] = mean((pred_val_test_high - y_test_high)^2)
  
  # Fisher
  temp_beta_high = mod_fish_high$beta_hat_matrix[,m] 
  pred_val_train_high = Z_train %*% temp_beta_high
  mod_fish_train_err_high[m] = mean((pred_val_train_high - y_train_high)^2)
  pred_val_test_high = Z_test %*% temp_beta_high 
  mod_fish_test_err_high[m] = mean((pred_val_test_high - y_test_high)^2)
}

ylims = c(0,6)
png(paste(save_prefix, "GenPartBoostR_Fish_MSE.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow = c(1,1))
  par(mar = c(3.5,3.5,2,1))
  
  plot(NA, type = 'n',
       main = "Mean Squared Error", xlim = c(0, num_iter),
       ylim = ylims, xlab = "", ylab = "")
  title(ylab="MSE", xlab = "Iterations", line = 2)
  ll = legend("topright", legend = c("GenPartBoostR", "PartBoostR", "Fisher",
                                     "Least Squares"),
         lty = 1, col = c(nice_colors[c(1,3,4)], 'darkgray'))
  
  kk = legend(-100, -100, lty = 1:2, col = 1,
              legend = c('Test', 'Train'), cex = 1)
  
  legend(ll$rect$left + (ll$rect$w - kk$rect$w),
         ll$rect$top - ll$rect$h,lty = 1:2, col = 1,
         legend = c('Test', 'Train'), cex = 1)
  
  abline(h = ls_inter_train_err, lty = 2, col = 'darkgray')
  lines(mod_gbr_train_err, lty = 2, col = nice_colors[1])
  lines(mod_gpbr_train_err, lty = 2, col = nice_colors[3])
  lines(mod_fish_train_err, lty = 2, col = nice_colors[4])
  
  abline(h = ls_inter_test_err, lty = 1, col = 'darkgray')
  lines(mod_gbr_test_err, lty = 1, col = nice_colors[1])
  lines(mod_gpbr_test_err, lty = 1, col = nice_colors[3])
  lines(mod_fish_test_err, lty = 1, col = nice_colors[4])
}
dev.off()


# Make a plot of their coefficient build ups (sigma_low = 1) together with the test mse (sigma = 5)
png(paste(save_prefix, "GenPartBoostR_coeff_build_up_siglow_MSE_sighigh.png", sep = ""), width = 3250, height = 1500, res = 350)
{
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))

  plot(NA, type = 'n',
       main = "Coefficient Build-Up", xlim = c(0,num_iter),
       ylim = c(-1.05, 5), col = nice_colors[5], xlab = "", ylab = "")
  title(ylab=bquote(beta['j']), xlab = "Iterations", line = 2)
  
  # The GenBoostR
  lines(0:num_iter, mod_gbr$beta_hat_matrix[1,], col = nice_colors[5], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[2,], col = nice_colors[1], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[3,], col = nice_colors[2], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[4,], col = nice_colors[3], type = 's', lty = 1)
  lines(0:num_iter, mod_gbr$beta_hat_matrix[5,], col = nice_colors[4], type = 's', lty = 1)
  
  # The GenPartBoostR
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[1,], col = nice_colors[5], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[2,], col = nice_colors[1], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[3,], col = nice_colors[2], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[4,], col = nice_colors[3], type = 's', lty = 2)
  lines(0:num_iter, mod_gpbr$beta_hat_matrix[5,], col = nice_colors[4], type = 's', lty = 2)
  
  # Fisher Scoring 
  lines(0:num_iter, mod_fish$beta_hat_matrix[1,], col = nice_colors[5], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[2,], col = nice_colors[1], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[3,], col = nice_colors[2], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[4,], col = nice_colors[3], type = 's', lty = 3)
  lines(0:num_iter, mod_fish$beta_hat_matrix[5,], col = nice_colors[4], type = 's', lty = 3)
  
  # Least squares solution
  points(rep(num_iter, p+1), ls_inter_sol)
  
  # Add legend
  legend("topleft", legend = c(as.expression(bquote(beta[0] ~ " ")),
                               as.expression(bquote(beta[1] ~ " ")),
                               as.expression(bquote(beta[2] ~ " ")),
                               as.expression(bquote(beta[3] ~ " ")),
                               as.expression(bquote(beta[4] ~ " "))),
         lty = 1, col = nice_colors[c(5,1,2,3,4)], cex = 0.75)
  temp = legend("topright", legend = c("Fisher", "GenBoostR", "GenPartBoostR"),
                lty = c(3,1,2), col = 'black', cex = 0.75)
  
  
  # Plot the test mse for sigma5
  plot(NA, type = 'n',
       main = "Test Mean Squared Error", xlim = c(0, num_iter),
       ylim = c(25, 30), xlab = "", ylab = "")
  title(ylab="MSE", xlab = "Iterations", line = 2)
  ll = legend("topright", legend = c("GenPartBoostR", "GenBoostR", "Fisher",
                                     "Least Squares"),
              lty = 1, col = c(nice_colors[c(4,1,3)], 'darkgray'), cex = 0.75)
  
  kk = legend(-100, -100, lty = c(2,1,3), col = 1,
              legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                         as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " "))),
              title = "BoostR", cex = 0.75)
  
  legend(ll$rect$left + (ll$rect$w - kk$rect$w),
         ll$rect$top - ll$rect$h, lty = c(2,1,3), col = 1, cex = 0.75,
         legend = c(as.expression(bquote(lambda ~ "=" ~ .(lambda_values[1]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[2]) ~ " ")),
                    as.expression(bquote(lambda ~ "=" ~ .(lambda_values[3]) ~ " "))))

  
  abline(h = ls_inter_test_err_high, lty = 1, col = 'darkgray')
  lines(mod_gbr_test_err_high, lty = 1, col = nice_colors[1])
  lines(mod_gbr_1_test_err_high, lty = 2, col = nice_colors[1])
  lines(mod_gbr_2_test_err_high, lty = 3, col = nice_colors[1])
  
  lines(mod_gpbr_test_err_high, lty = 1, col = nice_colors[4])
  lines(mod_gpbr_1_test_err_high, lty = 2, col = nice_colors[4])
  lines(mod_gpbr_2_test_err_high, lty = 3, col = nice_colors[4])
  
  lines(mod_fish_test_err_high, lty = 1, col = nice_colors[3])
}
dev.off()
