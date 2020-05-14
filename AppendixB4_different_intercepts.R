### File to test the intercept model
library(GAMBoost)
set.seed(2020)
stepno = 0

# Gaussian
family      = gaussian()
y_gaussian  = rnorm(100, 5, 1)
gaussian_v1 = family$linkfun(mean(y_gaussian))
gaussian_v2 = glm(y_gaussian ~ 1, family = family$family)$coeff
gaussian_v3 = GAMBoost(matrix(1:100), y_gaussian,
                       family = family$family, stepno = stepno)
gaussian_v4 = family$linkfun(gaussian_v3$hatmat %*% y_gaussian)[1]

# Binomial
family      = binomial()
y_binomial  = rbinom(100, 1, 0.3)
binomial_v1 = family$linkfun(mean(y_binomial))
binomial_v2 = glm(y_binomial ~ 1, family = family$family)$coeff
binomial_v3 = GAMBoost(matrix(1:100), y_binomial,
                       family = family$family, stepno = stepno)
binomial_v4 = family$linkfun(binomial_v3$hatmat %*% y_binomial)[1]

# Poisson
family     = poisson()
y_poisson  = rpois(100, 4)
poisson_v1 = family$linkfun(mean(y_poisson))
poisson_v2 = glm(y_poisson ~ 1, family = family$family)$coeff
poisson_v3 = GAMBoost(matrix(1:100), y_poisson,
                      family = family$family, stepno = stepno)
poisson_v4 = family$linkfun(poisson_v3$hatmat %*% y_poisson)[1]


# poisson_v3 = GAMBoost(matrix(runif(200), ncol = 2), y_poisson,
#                       family = family$family, stepno = stepno)

# Collect the results
intercepts = matrix(c(gaussian_v1, gaussian_v2, gaussian_v3$beta[[1]], gaussian_v4,
                      binomial_v1, binomial_v2, binomial_v3$beta[[1]], binomial_v4,
                      poisson_v1,  poisson_v2,  poisson_v3$beta[[1]],  poisson_v4),
                    ncol=4, byrow=TRUE)

colnames(intercepts) = c("link of mean", "glm", "GAMBoost", 'GAMBooost: Hat')
rownames(intercepts) = c("Gaussian", "Binomial", "Poisson")
intercepts = as.table(intercepts)
intercepts

#          link of mean        glm   GAMBoost GAMBooost: Hat
# Gaussian    5.1088918  5.1088918  5.1088918      5.1088918
# Binomial   -0.8472979 -0.8472979 -0.8000000     -0.8472979
# Poisson     1.3297240  1.3297240  2.7800000      1.3297240


###
family = gaussian()
y_response = y_gaussian

family = binomial()
y_response = y_binomial

family = poisson()
y_response = y_poisson


png(paste("Differnce_model_hatmat_", family$family, '.png',  sep = ""),
    width = 3500, height = 1750, res = 350)
par(mfrow=c(2,4), mar = c(3.5,3.5,2,1))
for (iter in c(0,1,2,3,5,25,100,500)) {
  mod = GAMBoost(matrix(1:100), y_response, family = family, stepno = iter)
  modmu = family$linkinv(mod$eta[,iter+1])
  hatmu = array(mod$hatmatrix %*% y_response)
  maxy = max(modmu, hatmu)
  miny = min(modmu, hatmu)
  
  plot(modmu, hatmu, main = paste('Iteration:', iter),
       xlim = c(miny, maxy), ylim = c(miny, maxy),
       col = 'darkgray', xlab = '', ylab = '')  
  mtext(side=1, expression(paste('Model: ', hat(mu))), line=2.25, cex = 0.7)
  mtext(side=2, expression(paste('Hat matrix: ', hat(mu))), line=2, cex = 0.7)
  abline(coef = c(0,1), col = 'black', lwd = 1)
  print(sum(family$dev.resids(modmu, hatmu, rep(1, 100))))
}
dev.off()


##### GAMBoost_stumps
mod = GAMBoost_stumps(matrix(1:100), y_response, num_iter = 500, lambda = 2,
                      family = family, print_msg = TRUE)

png(paste("Differnce_model_hatmat_stumps_", family$family, '.png',  sep = ""),
    width = 3500, height = 1750, res = 350)
par(mfrow=c(2,4), mar = c(3.5,3.5,2,1))
for (iter in c(0,1,2,3,5,25,100,500)) {
  modmu = family$linkinv(mod$eta[iter+1,])
  hatmu = array(mod$hat_mat[[iter+1]] %*% y_response)
  maxy = max(modmu, hatmu)
  miny = min(modmu, hatmu)
  
  plot(modmu, hatmu, main = paste('Iteration:', iter),
       xlim = c(miny, maxy), ylim = c(miny, maxy),
       col = 'darkgray', xlab = '', ylab = '')  
  mtext(side=1, expression(paste('Model: ', hat(mu))), line=2.25, cex = 0.7)
  mtext(side=2, expression(paste('Hat matrix: ', hat(mu))), line=2, cex = 0.7)
  abline(coef = c(0,1), col = 'black', lwd = 1)
  print(sum(family$dev.resids(modmu, hatmu, rep(1, 100))))
}
dev.off()


################
family = gaussian()
family = binomial()
family = poisson()

families = c('gaussian', 'binomial', 'poisson')
stnrs = c(0.25, 1, 3, 10)
for (stnr in stnrs) {
  for (family in families) {
    family = switch(family, gaussian = gaussian(), binomial = binomial(), poisson = poisson())
    
    data = create_data_combined(n = 100, seed_number = 123, stnr = stnr,
                                type = 'smooth', family = family)
    
    png(paste("Differnce_model_hatmat_large_data_", family$family, '_stnr_', stnr, '.png',  sep = ""),
        width = 3500, height = 1750, res = 350)
    par(mfrow=c(2,4), mar = c(3.5,3.5,2,1))
    for (iter in c(0,1,2,3,5,25,100,500)) {
      mod = GAMBoost(data$X, data$y, family = family, stepno = iter)
      modmu = family$linkinv(mod$eta[,iter+1])
      hatmu = array(mod$hatmatrix %*% data$y)
      maxy = max(modmu, hatmu)
      miny = min(modmu, hatmu)
      
      # mod$beta[[1]]
      # family$linkfun(mean(y))
      # glm(y ~ 1, family = family$family)$coeff
      
      plot(modmu, hatmu, main = paste('Iteration:', iter),
           xlim = c(miny, maxy), ylim = c(miny, maxy),
           col = 'darkgray', xlab = '', ylab = '')  
      mtext(side=1, expression(paste('Model: ', hat(mu))), line=2.25, cex = 0.7)
      mtext(side=2, expression(paste('Hat matrix: ', hat(mu))), line=2, cex = 0.7)
      abline(coef = c(0,1), col = 'black', lwd = 1)
      
      if (family$family == 'poisson') {
        hatmu[hatmu < 0] = 1e-10
      }
      if (family$family == 'binomial') {
        hatmu[hatmu < 0] = 1e-10
        hatmu[hatmu > 1] = 1-1e-10
      }
      print(sum(family$dev.resids(modmu, hatmu, rep(1, 100))))
    }
    dev.off()
  }
}

##### Similar for GAMBoost_stumps
families = c('gaussian', 'binomial', 'poisson')
stnrs = c(0.25, 1, 3, 10)
for (stnr in stnrs) {
  for (family in families) {
    family = switch(family, gaussian = gaussian(), binomial = binomial(), poisson = poisson())
    
    data = create_data_combined(n = 100, seed_number = 123, stnr = stnr,
                                type = 'smooth', family = family)
    data$cn
    mod = GAMBoost_stumps(data$X, data$y, num_iter = 250, lambda = 2,
                          family = family, print_msg = TRUE)
    
    png(paste("Differnce_model_hatmat_stumps_large_data_", family$family, '_stnr_', stnr, '.png',  sep = ""),
        width = 3500, height = 1750, res = 350)
    par(mfrow=c(2,4), mar = c(3.5,3.5,2,1))
    for (iter in c(0,1,2,3,5,25,100,250)) {
      modmu = family$linkinv(mod$eta[iter+1, ])
      hatmu = array(mod$hat_mat[[iter+1]] %*% data$y)
      maxy = max(modmu, hatmu)
      miny = min(modmu, hatmu)
      
      plot(modmu, hatmu, main = paste('Iteration:', iter),
           xlim = c(miny, maxy), ylim = c(miny, maxy),
           col = 'darkgray', xlab = '', ylab = '')  
      mtext(side=1, expression(paste('Model: ', hat(mu))), line=2.25, cex = 0.7)
      mtext(side=2, expression(paste('Hat matrix: ', hat(mu))), line=2, cex = 0.7)
      abline(coef = c(0,1), col = 'black', lwd = 1)
      
      if (family$family == 'poisson') {
        hatmu[hatmu < 0] = 1e-10
      }
      if (family$family == 'binomial') {
        hatmu[hatmu < 0] = 1e-10
        hatmu[hatmu > 1] = 1-1e-10
      }
      print(sum(family$dev.resids(modmu, hatmu, rep(1, 100))))
    }
    dev.off()
  }
}


