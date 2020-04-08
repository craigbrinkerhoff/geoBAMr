# utility functions

#' Convert coefficient of variation to sigma parameter of lognormal diistribution
#'
#' @param cv Coefficient of variation
#' @export

cv2sigma <- function (cv) {
  sqrt(log(cv^2 + 1))
}



# functions for getting Q bounds.

#' Minimum across xs of max across time of width
#'
#' @param x a numeric matrix
minmax <- function(x)
  min(apply(x, 1, max))

#' Maximum across xs of min across time of width
#'
#' @param x a numeric matrix
maxmin <- function(x) {
  max(apply(x, 1, min))
}

#Functions for classifying rivers

#'Classify river for expert framework
#'
#'@param Wobs observed widths matrix
classify_func <- function(Wobs) {
  lwbar <- mean(log(Wobs), na.rm=TRUE)#apply(log(Wobs), 1, mean, na.rm = TRUE)
  lwsd <- sd(log(Wobs), na.rm= TRUE)#apply(log(Wobs), 1, sd, na.rm = TRUE)

  maxWidth = 6.5
  classes <- c(2.641197894, 3.143289838,3.376904689,
               3.548754946,3.823191792,4.150094577,
               4.626687029,5.248706866) #median width of each river type

  ifelse(lwbar > maxWidth, 100,
         ifelse(lwsd >= 0.45, 9, which.min(abs(classes-lwbar)))) #100 for big rivers, 9 for width-variable rivers
}

#'Classify river for unsupervised framework
#'
#' @param Wobs observed widths matrix
classify_func_unsupervised <- function(Wobs) {
  width <- median((Wobs)) #note these are not log widths!
  width_sd <- sd((Wobs))

  #One-vs-Rest Logistic regression model
  #test set accuracy rate of 87%
  p_noise <- 1 / (1.0 + exp(-(-3.12907956 + 7.72682055e-04 * width_sd + 1.68285442e-03 * width)))
  p_1 <- 1 / (1.0 + exp(-(1.9719785 + -7.08487404e-04 * width_sd + -9.57359172e-04 * width)))
  p_2 <- 1 / (1.0 + exp(-(-4.0504766  + -1.30605126e-04 * width_sd + 1.60069188e-03 * width)))
  p_3<- 1 / (1.0 + exp(-(-4.34955267 + -6.70236303e-03 * width_sd + 7.59492911e-05 * width)))
  p_4 <- 1 / (1.0 + exp(-(-3.9977568  + 3.24040207e-04 * width_sd + 1.59933716e-03 * width)))
  p_5 <- 1 / (1.0 + exp(-(-4.26058253 + -2.39296957e-03 * width_sd + -1.41460664e-02 * width)))
  p_6 <- 1 / (1.0 + exp(-(-3.86984885 + -2.63864258e-04 * width_sd + -1.68608848e-03 * width)))
  p_7 <- 1 / (1.0 + exp(-(-1.15978639 + 5.68264551e-04 * width_sd + -1.27528508e-01 * width)))

  probs <- data.frame(p_noise, p_1, p_2, p_3, p_4, p_5, p_6, p_7)

  index <- which.max(probs)
  index <- ifelse(index == 1, 101, index-1) #101 for 'noisey' rivers
  return(index)
}
