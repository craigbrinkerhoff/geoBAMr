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
  min(apply(x, 1, max, na.rm=TRUE), na.rm=TRUE)

#' Maximum across xs of min across time of width
#'
#' @param x a numeric matrix
maxmin <- function(x) {
  max(apply(x, 1, min, na.rm=TRUE), na.rm=TRUE)
}

#Functions for classifying rivers----------------------------------------------------------------------

#'Classify river for expert framework
#'
#'@param Wobs observed widths matrix
classify_func <- function(Wobs) {
  lwbar <- mean(log(Wobs), na.rm=TRUE)
  lwsd <- sd(log(Wobs), na.rm= TRUE)

  maxWidth = 6.5
  classes <- c(2.476118144,
               2.864001065,
               3.103015939,
               3.249308032,
               3.284178964,
               3.371669039,
               3.56827873,
               3.664586762,
               3.683922384,
               4.002696788,
               4.031559142,
               4.357733942,
               4.436574004,
               4.921166637,
               5.287893051) #median width of each river type
  index <- ifelse(lwbar > maxWidth, 17, which.min(abs(classes-lwbar))) #17 for big rivers
  index <- ifelse(lwsd >= 0.45, 16, index)  #16 for width-variable rivers, which overrides 'big' rivers
  return(index)
}

#'Classify river for unsupervised framework
#'
#' @param Wobs observed widths matrix
classify_func_unsupervised <- function(Wobs) {
  width <- median((Wobs), na.rm = TRUE) #note these are not log widths!

  #One-vs-Rest Logistic regression model
  #test set accuracy rate of 87%
  p_noise <- 1 / (1.0 + exp(-(-3.11056777 + 0.00189261 * width)))
  p_1 <- 1 / (1.0 + exp(-(1.95763357 + -0.00114445 * width)))
  p_2 <- 1 / (1.0 + exp(-(-4.05274822 + 0.00156913 * width)))
  p_3<- 1 / (1.0 + exp(-(-4.43158034 + -0.00116211 * width)))
  p_4 <- 1 / (1.0 + exp(-(-3.99272449 + 0.00168791* width)))
  p_5 <- 1 / (1.0 + exp(-(-4.28204232 + -0.01493184 * width)))
  p_6 <- 1 / (1.0 + exp(-(-3.87408245 + -0.00175089 * width)))
  p_7 <- 1 / (1.0 + exp(-(-1.16209678 + -0.12664823 * width)))

  probs <- data.frame(p_noise, p_1, p_2, p_3, p_4, p_5, p_6, p_7)

  index <- which.max(probs)
  index <- ifelse(index == 1, 101, index-1) #101 for 'noisey' rivers
  return(index)
}
