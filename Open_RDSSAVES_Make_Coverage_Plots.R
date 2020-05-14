### FILE TO OPEN THE RDSSAVES GENERATED FROM
### CHAPTER443_COVERAGE_GAMBOOST_SPLINES AND CHAPTER444_COVERAGE_GAMBOOST_STUMPS
### CONTACT LHOLSEN@MATH.UIO.NO FOR OBTAINING MY SAVES

# Some parameters needed to read the files
# All possible cobinations of types, families, SNTRs and lambdas
base_learner = 'stumps' 
types = c('smooth', 'step', 'advancedstep')
stnrs = c(1,3,10)
families = c('gaussian', 'binomial', 'poisson')
lambdas = c(2, 4, 10, 25)
B_outer = 200
B_inner = 200


for (family in families) {
  for (type in types) {
    for (STNR in stnrs) {
      for (lambda in lambdas) {
        temp = readRDS(paste(base_learner, type, "_", family, "_stnr_", STNR,
                             '_lambda_', lambda, '_bootinner_', B_inner,
                             '_bootouter_', B_outer, "_final2", sep=''))
        
        type = temp$type
        family = temp$family
        STNR = temp$stnr
        lambda = temp$lambda
        B_inner = temp$B_in
        B_outer = temp$B_out
        p = ncol(temp$training$X)
        n_train = nrow(temp$training$X)
        inside_band_approximate = temp$inside_band_app
        inside_band_empirical = temp$inside_band_emp
        mu_inside_interval_approximate = temp$mu_inside_interval_app
        mu_inside_interval_empirical = temp$mu_inside_interval_emp
        training_data = temp$training
        height = 1250
        colors = FALSE
        lines = TRUE
        
        png(paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                  '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                  '_coverage_f_lines_points_colorless.png',  sep = ""),
            width = 3500, height = height, res = 350)
        {
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
              if (lines) {
                lines(sort(training_data$X[,s]), coverage_approx_temp, col = nice_colors[1], lwd = 0.75)
                lines(sort(training_data$X[,s]), coverage_empir_temp,  col = nice_colors[4], lwd = 0.75) 
              }
              
              points(sort(training_data$X[,s]), coverage_approx_temp, col = nice_colors[1], pch = 19, cex = 0.6)
              points(sort(training_data$X[,s]), coverage_empir_temp,  col = nice_colors[4], pch = 19, cex = 0.6)          
            } else {
              if (lines) {
                lines(sort(training_data$X[,s]), coverage_approx_temp, lty = 1, lwd = 0.75, col = 'black')
                lines(sort(training_data$X[,s]), coverage_empir_temp,  lty = 1, lwd = 0.75, col = "#676767")
              }
              points(sort(training_data$X[,s]), coverage_approx_temp, col = 'black', pch = 19, cex = 0.6)
              points(sort(training_data$X[,s]), coverage_empir_temp, col = "darkgray", pch = 19, cex = 0.6)
            }
            
            print(c(median(coverage_approx_temp), median(coverage_empir_temp)))
          }
        }
        dev.off()
        
        png(paste("stumps_", type, "_", family$family, "_stnr_", STNR,
                  '_lambda_', lambda, '_bootinner_', B_inner, '_bootouter_', B_outer, 
                  '_coverage_mu_lines_points_colorless.png',  sep = ""),
            width = 3500, height = height, res = 350)
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
            if (lines) {
              lines(1:n_train, coverage_approx_mu,    col = nice_colors[1], lwd = 0.75)
              lines(1:n_train, coverage_empircal_mu,  col = nice_colors[4], lwd = 0.75)
            }
            points(1:n_train, coverage_approx_mu,   col = nice_colors[1], pch = 19, cex = 0.6)
            points(1:n_train, coverage_empircal_mu, col = nice_colors[4], pch = 19, cex = 0.6) 
          } else {
            if (lines) {
              lines(1:n_train,  coverage_approx_mu,   col = 'black',   lty = 1,  lwd = 0.75)
              lines(1:n_train,  coverage_empircal_mu, col = "darkgray", lty = 1,  lwd = 0.75)
            }
            points(1:n_train, coverage_approx_mu,   col = 'black',   pch = 19, cex = 0.6)
            points(1:n_train, coverage_empircal_mu, col = "darkgray", pch = 19, cex = 0.6)
          }
          
          print(c(median(coverage_approx_mu), median(coverage_empircal_mu)))
        }
        dev.off()
      }
    }
  }
}



