#                                                                     
# Written by:                                                         
# --                                                                  
# Lars Henry Berge Olsen                03.03.2020                        
#                                                                     
# Creates figures for splines and trees in chapter 3 in master thesis

################
##### LARS #####
################
rm(list = ls())
set.seed(0)
nice_colors = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)

library(splines)
require(graphics)

xsmall = 0 
xlarge = 4*pi
x = seq(xsmall, xlarge, length.out=50)
y = cos(x) + 0.3*rnorm(length(x))
knots_num = 2
knots_loc = (xlarge - xsmall) / (knots_num + 1) * 1:knots_num

# degree_want = 5
# knots_loc_aug = (xlarge - xsmall) / (knots_num + 1) * -(degree_want-1):(knots_num+degree_want)
# knots_loc_aug
# 
# bs1 = bs(x, degree=2, knots = knots_loc_aug, 
#          Boundary.knots = c(knots_loc_aug[1], knots_loc_aug[length(knots_loc_aug)]))
# matplot(x, bs1, ylim=c(0,max(bs1)), type='l', lwd=2, col=1:10, lty = 1, 
#         xlab="Cubic B-spline basis", ylab="")
#   
# mod1 = lm(y ~ bs1)
# mod1$coefficients
# attr(mod1$terms, "predvars")
# plot( x, y, type="p" );
# lines( x, fitted(mod1) )
# 
# 
# x <- seq(0, 1, by=0.0001)
# spl <- bs(x, degree = 2, knots = c(-0.25, 0, 0.25, 0.5, 0.75, 1, 1.25),
#           Boundary.knots = c(-0.25,1.25))
# spl <- bs(x, degree = 3, knots = c(0, 0.25, 0.5, 0.75, 1),
#           Boundary.knots = c(0,1))
# spl <- bs(x, degree = 3, knots = c(0, 0.2, 0.4, 0.6, 0.8),
#           Boundary.knots = c(-0,1))
# plot(spl[,1]~x, ylim=c(0,max(spl)), type='l', lwd=2, col=nice_colors[1], 
#      xlab="Cubic B-spline basis", ylab="")
# for (j in 2:ncol(spl)) {
#   lines(spl[,j]~x, lwd=2, col=nice_colors[j])
#   j= j+1
# }
# matplot(x, spl, ylim=c(0,max(spl)), type='l', lwd=2, col=nice_colors, lty = 1, 
#         xlab="Cubic B-spline basis", ylab="")
# spl[1,]
# 
# plot(1:10, 1:10, col = 6)
# 
# 
# 
# spl  <- bs(x, degree = 3, knots = c(0, 0.2, 0.4, 0.6, 0.8),
#            Boundary.knots = c(-0,1))
# spl2 <- bs(x, degree = 3, knots = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4),
#            Boundary.knots = c(-0.6, 1.6))
# apply(spl2, 2, sum)
# spl2[,4:11]

##########################
##### PLOT OF BASIS ######
##########################

x <- seq(0, 1, by=0.0001)
spl  <- bs(x, degree = 3, knots = c(0, 0.2, 0.4, 0.6, 0.8),
           Boundary.knots = c(-0,1))
spl2 <- bs(x, degree = 3, knots = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4),
           Boundary.knots = c(-0.6, 1.6))

spl2 = spl2[,4:11] ### Remove extra basises that are always zero since outside intertior

png("CubicBSplineBoundaryKnots.png", width = 3250, height = 1850, res = 350)
{
  par(mfrow = c(2,1))
  par(mar = c(3.5,3.5,2,1))
  matplot(x, spl, ylim=c(0,1.1), type='l', lwd=2, col=nice_colors, lty = 1, 
          xlab="", ylab="")
  lines(x, apply(spl, 1, sum), lty = 2, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 2, col = 'darkgray')
  title(xlab="Cubic B-Spline Basis", line = 2, cex.main = 1.5)
  title(main = "Cubic B-Spline Basis with Four Equally Spaced Interior Knots", line = 0.6)
  
  
  matplot(x, spl2, ylim=c(0,1.1), type='l', lwd=2, col=nice_colors, lty = 1, 
          xlab="", ylab="")
  lines(x, apply(spl2, 1, sum), lty = 2, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 2, col = 'darkgray')
  title(main = "Cubic B-Spline Basis with Equally Spaced Boundary and Interior Knots", line = 0.6)
  title(xlab="Cubic B-Spline Basis", line = 2, cex.main = 1.5)
}
dev.off()


png("CubicBSplineBoundaryKnotsColorless.png", width = 3250, height = 1850, res = 350)
{
  par(mfrow = c(2,1))
  par(mar = c(3.5,3.5,2,1))
  matplot(x, spl, ylim=c(0,1.1), type='l', lwd=1.5, col='black', lty = rep(1:2, 8), 
          xlab="", ylab="")
  lines(x, apply(spl, 1, sum), lty = 4, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 4, col = 'darkgray')
  title(xlab="Cubic B-Spline Basis", line = 2, cex.main = 1.5)
  title(main = "Cubic B-Spline Basis with Four Equally Spaced Interior Knots", line = 0.6)
  
  
  matplot(x, spl2, ylim=c(0,1.1), type='l', lwd=1.5, col='black', lty = rep(1:2, 8), 
          xlab="", ylab="")
  lines(x, apply(spl2, 1, sum), lty = 4, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 4, col = 'darkgray')
  title(main = "Cubic B-Spline Basis with Equally Spaced Boundary and Interior Knots", line = 0.6)
  title(xlab="Cubic B-Spline Basis", line = 2, cex.main = 1.5)
}
dev.off()


#####################################################
##### PLOT UNIFORM KNOTS FOR DIFFERENT DEGREES ######
#####################################################

x = seq(0, 1, by=0.004)
knots = c(0, 0.2, 0.4, 0.6, 0.8, 1)
{
  # Order 1
  spl0 = spline.des(knots, x = x, ord=1, outer.ok = F)$design
  
  # Order 2
  spl1 = bs(x, degree = 1, knots = c(0, 0.2, 0.4, 0.6, 0.8),
            Boundary.knots = c(-0,1))
  
  # Order 3
  spl2 = bs(x, degree = 2, knots = c(0, 0.2, 0.4, 0.6, 0.8, 1),
            Boundary.knots = c(-0.2,1.2))
  spl2 = spl2[,1:7]
  
  # Order 4
  spl3 = bs(x, degree = 3, knots = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2),
            Boundary.knots = c(-0.4,1.4))
  spl3 = spl3[,2:9]
}


png("BSplineDifferentOrder.png", width = 3250, height = 2750, res = 350)
{
  par(mfrow = c(4,1))
  par(mar = c(3.5,3.5,2,1))
  
  # M = 1
  matplot(x, spl0, ylim=c(0,1.1), type='l', lwd=2, col=nice_colors, lty = 1, 
          xlab="", ylab="")
  lines(x, apply(spl0, 1, sum), lty = 2, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 2, col = 'darkgray')
  title(xlab="Constant B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 1", line = 0.6, cex.main=1.5)
  
  # M = 2
  matplot(x, spl1, ylim=c(0,1.1), type='l', lwd=2, col=nice_colors, lty = 1, 
          xlab="", ylab="")
  lines(x, apply(spl1, 1, sum), lty = 2, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 2, col = 'darkgray')
  title(xlab="Linear B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 2", line = 0.6, cex.main=1.5)
  
  # M = 3
  matplot(x, spl2, ylim=c(0,1.1), type='l', lwd=2, col=nice_colors, lty = 1, 
          xlab="", ylab="")
  lines(x, apply(spl2, 1, sum), lty = 2, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 2, col = 'darkgray')
  title(xlab="Quadratic B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 3", line = 0.6, cex.main=1.5)
  
  # M = 4
  matplot(x, spl3, ylim=c(0,1.1), type='l', lwd=2, col=nice_colors, lty = 1, 
          xlab="", ylab="")
  lines(x, apply(spl3, 1, sum), lty = 2, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 2, col = 'darkgray')
  title(xlab="Cubic B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 4", line = 0.6, cex.main=1.5)
}
dev.off()


png("BSplineDifferentOrderColorless.png", width = 3250, height = 2750, res = 350)
{
  par(mfrow = c(4,1))
  par(mar = c(3.5,3.5,2,1))
  
  # M = 1
  matplot(x, spl0, ylim=c(0,1.1), type='l', lwd=1.5, col='black', lty = rep(1:2, 8), 
          xlab="", ylab="")
  lines(x, apply(spl0, 1, sum), lty = 4, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 4, col = 'darkgray')
  title(xlab="Constant B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 1", line = 0.6, cex.main=1.5)
  
  # M = 2
  matplot(x, spl1, ylim=c(0,1.1), type='l', lwd=1.5, col='black', lty = rep(1:2, 8), 
          xlab="", ylab="")
  lines(x, apply(spl1, 1, sum), lty = 4, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 4, col = 'darkgray')
  title(xlab="Linear B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 2", line = 0.6, cex.main=1.5)
  
  # M = 3
  matplot(x, spl2, ylim=c(0,1.1), type='l', lwd=1.5, col='black', lty = rep(1:2, 8), 
          xlab="", ylab="")
  lines(x, apply(spl2, 1, sum), lty = 4, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 4, col = 'darkgray')
  title(xlab="Quadratic B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 3", line = 0.6, cex.main=1.5)
  
  # M = 4
  matplot(x, spl3, ylim=c(0,1.1), type='l', lwd=1.5, col='black', lty = rep(1:2, 8), 
          xlab="", ylab="")
  lines(x, apply(spl3, 1, sum), lty = 4, col = 'darkgray')
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8, 1), lty = 4, col = 'darkgray')
  title(xlab="Cubic B-Spline Basis", line = 2)
  title(main = "The B-Spline Basis of Order 4", line = 0.6, cex.main=1.5)
}
dev.off()

##############################################################
##### USE BSPLINES TO FIT A FUNCTION USING LEAST SQUARES #####
##############################################################
set.seed(123)

xsmall = 0 
xlarge = 4*pi
x = seq(xsmall, xlarge, length.out=50)
y = cos(x) + 0.3*rnorm(length(x))
knots_num = 3
knots_loc = (xlarge - xsmall) / (knots_num + 1) * 1:knots_num

plot_splines = function(x_low, x_up, x_length, degrees, knots_num, noise_scalar = 0.3, seed_num = 123) {
  
  set.seed(seed_num)
  x = seq(x_low, x_up, length.out=x_length)
  y = cos(x) + noise_scalar*rnorm(length(x))
  
  nice_colors = c("#33A02C", "#1F78B4", "#E31A1C", "#FF7F00", "#6A3D9A",
                  "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
                  "#CAB2D6", "#FFFF99", 'black' , 'grey', 2:30)
  
  knots_loc = (x_up - x_low) / (knots_num + 1) * 1:knots_num
  
  knots_loc_aug = lapply(degrees, function(x) 
    (x_up - x_low) / (knots_num + 1) * (-(x-1)):(knots_num+x))
  
  
  fitted_values = matrix(NA, nrow = length(degrees), ncol = length(x))
  if (0 %in% degrees) {
    fitted_values[1, ] = rep(mean(y), length(x)) 
    next_index = 2
  } else {
    next_index = 1
  }
  
  for (j in 1:length(degrees)) {
    degree = degrees[j]
    if (degree != 0) {
      
      bs_temp = bs(x, degree=degree, knots = knots_loc_aug[[j]], 
                   Boundary.knots = c(knots_loc_aug[[j]][1],
                                      knots_loc_aug[[j]][length(knots_loc_aug[[j]])]))
      
      lm_mod = lm(y ~ bs_temp-1)
      print(mean((lm_mod$fitted.values - y)^2))
      
      fitted_values[next_index, ] = fitted(lm_mod)
      
      next_index = next_index + 1
    }
  }
  plot(x, y, type = 'p',  col = 'darkgray', xlab = "", ylab = "", xlim = c(x_low, x_up+1), ylim = c(-1.65, 1.65))
  lines(x, cos(x), col = "black", lwd = 2)
  abline(v = knots_loc, lty = 4, col = 'darkgray')
  matlines(x, t(fitted_values), col = nice_colors, lty = 1, lwd = 2)
  legend("topright", legend = c("Truth", c(degrees)+1), col = c('black', nice_colors),
         lty = 1, title = "Order", cex = 0.75)
  title(xlab = "x", ylab = "f(x)", line = 2)
  title(main = "Fitted B-Splines of Different Orders", line = 0.6, cex.main=1.5)
}


png("FittedCubicSplines.png", width = 3250, height = 1850, res = 350)
 {
  par(mfrow = c(2,1))
  par(mar = c(3.5,3.5,2,1))
  plot_splines(x_low = 0, x_up = 4*pi, x_length = 50, degrees = c(0,1,2,3),
               knots_num = 2, noise_scalar = 0.3, seed_num = 123) 
  plot_splines(x_low = 0, x_up = 4*pi, x_length = 50, degrees = c(0,1,2,3),
               knots_num = 3, noise_scalar = 0.3, seed_num = 123) 
}
dev.off()         


plot_splines_colorless = function(x_low, x_up, x_length, degrees, knots_num, noise_scalar = 0.3, seed_num = 123) {
  
  set.seed(seed_num)
  x = seq(x_low, x_up, length.out=x_length)
  y = cos(x) + noise_scalar*rnorm(length(x))
  
  knots_loc = (x_up - x_low) / (knots_num + 1) * 1:knots_num
  
  knots_loc_aug = lapply(degrees, function(x) 
    (x_up - x_low) / (knots_num + 1) * (-(x-1)):(knots_num+x))
  
  
  fitted_values = matrix(NA, nrow = length(degrees), ncol = length(x))
  if (0 %in% degrees) {
    fitted_values[1, ] = rep(mean(y), length(x)) 
    next_index = 2
  } else {
    next_index = 1
  }
  
  for (j in 1:length(degrees)) {
    degree = degrees[j]
    if (degree != 0) {
      
      bs_temp = bs(x, degree=degree, knots = knots_loc_aug[[j]], 
                   Boundary.knots = c(knots_loc_aug[[j]][1],
                                      knots_loc_aug[[j]][length(knots_loc_aug[[j]])]))
      
      lm_mod = lm(y ~ bs_temp-1)
      print(mean((lm_mod$fitted.values - y)^2))
      
      fitted_values[next_index, ] = fitted(lm_mod)
      
      next_index = next_index + 1
    }
  }
  plot(x, y, type = 'p',  col = 'darkgray', xlab = "", ylab = "", xlim = c(x_low, x_up+1), ylim = c(-1.65, 1.65))
  lines(x, cos(x), col = "darkgray", lwd = 3)
  abline(v = knots_loc, lty = 4, col = 'darkgray')
  matlines(x, t(fitted_values), col = 'black', lty = 1:10, lwd = 1.5)
  legend("topright", legend = c("Truth", c(degrees)+1), col = c('darkgray', rep('black', 10)),
         lty = c(1, 1:10), title = "Order", cex = 0.75)
  title(xlab = "x", ylab = "f(x)", line = 2)
  title(main = "Fitted B-Splines of Different Orders", line = 0.6, cex.main=1.5)
}

png("FittedCubicSplinesColorless.png", width = 3250, height = 1850, res = 350)
{
  par(mfrow = c(2,1))
  par(mar = c(3.5,3.5,2,1))
  plot_splines_colorless(x_low = 0, x_up = 4*pi, x_length = 50, degrees = c(0,1,2,3),
                         knots_num = 2, noise_scalar = 0.3, seed_num = 123) 
  plot_splines_colorless(x_low = 0, x_up = 4*pi, x_length = 50, degrees = c(0,1,2,3),
                         knots_num = 3, noise_scalar = 0.3, seed_num = 123) 
}
dev.off()    
               

###################################
##### CREATE THE TREE FIGURES #####
###################################
# library(tree)
# 
# x_low = -1
# x_up = 1
# x_length = 50
# noise_sd = 0.1
# set.seed(123)
# x = seq(x_low, x_up, length.out=x_length)
# y = x + rnorm(length(x), 0, noise_sd)
# plot(x,y)
# 
# training = data.frame(response=y, x=x)
# 
# tree_mod = tree(response~x, training, control = tree.control(x_length, mincut = 8, minsize = 20))
# summary(tree_mod)
# plot(tree_mod)
# text(tree_mod)
# pred_tree = predict(tree_mod, training)
# plot(x, pred_tree)

# library(rpart)
# rpart_mod = rpart(response~x, training, control =  rpart.control(maxdepth = 4, minsplit = 5))
# summary(rpart_mod)
# plot(rpart_mod)
# text(rpart_mod)
# rpart.plot(rpart_mod)
# prp(rpart_mod)					# Will plot the tree
# prp(rpart_mod, varlen=3)
# fancyRpartPlot(rpart_mod, main = "", sub ="", type = 2, caption = "hei", palettes=c("Oranges", "Greys"))
# pred_rpart = predict(rpart_mod, training)
# plot(x, pred_rpart, type = "s")




##############################################
##### CREATE THE PLOT OF TREE AND STUMPS #####
##############################################
library(rpart)				      # Popular decision tree algorithm
library(rattle)					    # Fancy tree plot
library(RColorBrewer)				# Color selection for fancy tree plot

x_low = 0
x_up = 2
x_length = 51
x = seq(x_low, x_up, length.out=x_length)

noise_sd = 0.3
set.seed(123)
y = x + rnorm(length(x), 0, noise_sd)
training = data.frame(response=y, x=x)


stump_mod = rpart(response~x, training, control =  rpart.control(maxdepth = 1))
tree_mod  = rpart(response~x, training)
tree_mod2 = rpart(response~x, training, control =  rpart.control(maxdepth = 8, minsplit = 2))


fancyRpartPlot(stump_mod, main = "", sub ="", type = 2, palettes=c("Oranges", "Greys"))
fancyRpartPlot(tree_mod,  main = "", sub ="", type = 2, palettes=c("Oranges", "Greys"))
fancyRpartPlot(tree_mod2,  main = "", sub ="", type = 2, palettes=c("Oranges", "Greys"))


stump_pred = predict(stump_mod, training)
tree_pred  = predict(tree_mod,  training)
tree_pred2 = predict(tree_mod2, training)

png("Trees.png", width = 3250, height = 1850, res = 350)
{
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))
  plot(x,y, col = 'darkgray', main = "", xlab = "", ylab = "")
  lines(x, x, lwd = 3, lty = 1, col = 'darkgray')
  lines(x, stump_pred, type = "s", col = nice_colors[1], lwd = 1.5)
  lines(x, tree_pred, type = "s", col = nice_colors[2], lty = 2, lwd = 1.5)
  lines(x, tree_pred2, type = "s", col = nice_colors[3], lty = 4, lwd = 1.5)
  title(xlab="x", ylab = "f(x)", line = 2)
  title(main = "Fitted Trees", line = 0.6, cex.main = 1.22)
  legend("bottomright", legend = c(2,7,15,"Truth"), cex = 0.8,
         lty = c(1,2,4,1), col = c(nice_colors[1:3], 'darkgray'), title="Tree Size")
  
  fancyRpartPlot(tree_mod, main = "Flowchart Diagram         ",
                 sub ="", type = 2, palettes=c("Greens", "Greys"))
}
dev.off()

png("TreesColorless.png", width = 3250, height = 1850, res = 350)
{
  par(mfrow = c(1,2))
  par(mar = c(3.5,3.5,2,1))
  grays = c("#F7F7F7", "#CCCCCC", "#969696", "#636363", "#252525")
  plot(x,y, col = 'darkgray', main = "", xlab = "", ylab = "")
  lines(x, x, lwd = 3, lty=1, col = 'darkgray')
  lines(x, stump_pred, type = "s", col = 1, lwd = 1.5)
  lines(x, tree_pred, type = "s", col = 1, lty = 2, lwd = 1.5)
  lines(x, tree_pred2, type = "s", col = 1, lty = 3, lwd = 1.5)
  title(xlab="x", ylab = "f(x)", line = 2)
  title(main = "Fitted Trees", line = 0.6, cex.main = 1.22)
  legend("bottomright", legend = c(2,7,15,"Truth"), cex = 0.8,
         lty = c(1:3,1), col = c('black','black','black','darkgray'),
         title="Tree Size")
  fancyRpartPlot(tree_mod,  main = "Flowchart Diagram         ",
                 sub ="", type = 2, palettes=c("Greys"))
}
dev.off()

