#' Title Nonlinear Cross Lag Analysis: autoregression
#'
#' @param xname If cross lagged analysis is used between x and y, 'xname' is the name of x
#' @param yname If cross lagged analysis is used between x and y, 'yname' is the name of y
#' @param data.x1 A numeric variable.
#' @param data.y1 A numeric variable. 'data.x1' and 'data.y1' comes from the first time point
#' @param data.xt A numeric variable.
#' @param data.yt A numeric variable. 'data.xt' and 'data.yt' comes from the another time point
#'
#' @return A dataframe containing the result of autoregression: 'Xt~X1' 'Yt~Y1'
#' @importFrom mgcv gam
#' @export
#' @examples
#' data(test_data1)
#' data(test_data2)
#' clpm_gam_r(xname="ASI",yname = "PWRI",data.x1 = test_data1$ASI,
#'            data.y1 = test_data1$PWRI,data.xt = test_data2$ASI,data.yt = test_data2$PWRI)
clpm_gam_r <- function(xname, yname, data.x1,data.y1,data.xt,data.yt){
  # Ensure that the function parameters are of the correct data types with the main calculations.
  if(!inherits(data.x1,"numeric",which = FALSE))
    stop(" data.x1 must be of class 'numeric' ")
  if(!inherits(data.y1,"numeric",which = FALSE))
    stop(" data.y1 must be of class 'numeric' ")
  if(!inherits(data.xt,"numeric",which = FALSE))
    stop(" data.xt must be of class 'numeric' ")
  if(!inherits(data.yt,"numeric",which = FALSE))
    stop(" data.yt must be of class 'numeric' ")
  if(!inherits(xname,"character",which = FALSE))
    stop(" xname must be of class 'character' ")
  if(!inherits(yname,"character",which = FALSE))
    stop(" yname must be of class 'character' ")
  # Extracte the numeric data from the input parameters for further calculations in the function.
  X1 <- data.x1
  Y1 <- data.y1
  Xt <- data.xt
  Yt <- data.yt
  # Fit a Generalized Additive Model (GAM) to the data with the response variable `Xt` and predictors `Y1` and `X1`.
  fit <- mgcv::gam(Xt ~ s(Y1, k=4) + X1)
  # Storing the summary of the fitted GAM in the variable `s`.
  # This summary typically includes information such as coefficients, standard errors, t-values, p-values, R-squared values, and other diagnostic statistics related to the model fit.
  s <- summary(fit)
  # Calculate the residuals of the model fit.
  r <- Xt - predict(fit)
  # Calculate the Root Mean Residual, which is a measure of the average residual error of the model fit.
  rmr <- sqrt(mean(r^2))

  # Constructing a summary of the GAM fitting process for the relationship between the variables `Xt` and `Y1`.
  result0 <- data.frame(lhs="Xt", op = "~", rhs="X1", beta=s$p.coeff[2], se=s$se[2],
                        t=s$p.t[2], p=s$p.pv[2], x=xname, y=yname, r.sq=s$r.sq, gcv=s$sp.criterion[1], rmr=rmr)
  # Fit another GAM to the data with the response variable `Yt` and predictors `X1` and `Y1`.
  fit <- mgcv::gam(Yt ~ s(X1, k=4) + Y1)
  # Storing the summary of the fitted GAM in the variable `s`. This summary typically includes information such as coefficients, standard errors, t-values, p-values, R-squared values, and other diagnostic statistics related to the model fit.
  s <- summary(fit)
  # Calculate the residuals of the model fit.
  r <- Xt - predict(fit)
  # Calculate the Root Mean Residual.
  rmr <- sqrt(mean(r^2))
  # Create a data frame named `result1` that contains information about the GAM fitting process for the relationship between the response variable `Yt` and the predictor `Y1`.
  result1 <- data.frame(lhs="Yt", op = "~", rhs="Y1", beta=s$p.coeff[2], se=s$se[2],
                        t=s$p.t[2], p=s$p.pv[2], x=xname, y=yname, r.sq=s$r.sq, gcv=s$sp.criterion[1], rmr=rmr)
  # Combine the data frames `result0` and `result1` row-wise into a single data frame named `output`.
  output <- rbind(result0,result1)
  # Modify the column values in the `rst` data frame into certain name.
  rst <- output
  rst$lhs[which(rst$lhs=="Xt")] <- paste0(xname,"_2")
  rst$lhs[which(rst$lhs=="Yt")] <- paste0(yname,"_2")
  rst$rhs[which(rst$rhs=="Y1")] <- paste0(yname,"_1")
  rst$rhs[which(rst$rhs=="X1")] <- paste0(xname,"_1")
  # Assign the string '   -->   ' to the column `op` in the data frame `rst`.
  # Represent an arrow symbol indicating a relationship between two variables.
  rst$op <- '   -->   '
  # Within each iteration, it is concatenating the values in the columns `lhs`, `op`, and `rhs` of that row using the `paste0` function.
  for (i in 1:nrow(rst)) {
    rst$lhs[i] <- paste0(rst$lhs[i],rst$op[i],rst$rhs[i])
  }
  # Remove the 2nd, 3rd, and 4th columns from the 'rst' data frame
  rst <- rst[,-c(2,3,4)]
  # Rename the column named 'lhs' to 'Relation' in the 'rst' data frame
  names(rst)[names(rst) =="lhs"] <-"Relation"
  # Rename the column named 'p' to 'pvalue' in the 'rst' data frame
  names(rst)[names(rst) =="p"] <-"pvalue"
  # Return the modified 'rst' data frame
  return(rst)
}
