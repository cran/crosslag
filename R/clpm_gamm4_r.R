#' Title Cross-lag analysis based on generalized additive mixture models: autoregression
#'
#' @param xname If cross lagged analysis is used between x and y, 'xname' is the name of x
#' @param yname If cross lagged analysis is used between x and y, 'yname' is the name of y
#' @param data.x1 A numeric variable.
#' @param data.y1 A numeric variable. 'data.x1' and 'data.y1' comes from the first time point
#' @param data.xt A numeric variable.
#' @param data.yt A numeric variable. 'data.xt' and 'data.yt' comes from the another time point
#' 'data.x1', 'data.y1','data.xt' and 'data.yt' can be the data processed by the function 'adjust_target()'
#' @param z In the generalized additive mixture model, a random intercept is specified and the random effects are grouped by z.
#' 'z' and 'data.x1', 'data.y1','data.xt','data.yt' are all numeric variables that should ideally correspond to each other.
#' @return A dataframe containing the result of autoregression based on generalized additive mixture model:'Xt~X1' 'Yt~Y1'.
#' @importFrom gamm4 gamm4
#' @export
#' @examples
#' data(test_data1)
#' data(test_data2)
#' result <- clpm_gamm4_r("ASI","PWRI",test_data1$ASI,test_data1$PWRI,
#'                      data.xt = test_data2$ASI,data.yt = test_data2$PWRI,z=test_data1$time)
clpm_gamm4_r <- function(xname, yname, data.x1,data.y1,data.xt,data.yt,z){
  # Check if 'data.x1' is of class 'numeric', if not, stop the function and print an error message
  if(!inherits(data.x1,"numeric",which = FALSE))
    stop(" data.x1 must be of class 'numeric' ")
  # Check if 'data.y1' is of class 'numeric', if not, stop the function and print an error message
  if(!inherits(data.y1,"numeric",which = FALSE))
    stop(" data.y1 must be of class 'numeric' ")
  # Check if 'data.xt' is of class 'numeric', if not, stop the function and print an error message
  if(!inherits(data.xt,"numeric",which = FALSE))
    stop(" data.xt must be of class 'numeric' ")
  # Check if 'data.yt' is of class 'numeric', if not, stop the function and print an error message
  if(!inherits(data.yt,"numeric",which = FALSE))
    stop(" data.yt must be of class 'numeric' ")
  # Check if 'z' is of class 'numeric', if not, stop the function and print an error message
  if(!inherits(z,"numeric",which = FALSE))
    stop(" z must be of class 'numeric' ")
  # Check if 'xname' is of class 'character', if not, stop the function and print an error message
  if(!inherits(xname,"character",which = FALSE))
    stop(" xname must be of class 'character' ")
  # Check if 'yname' is of class 'character', if not, stop the function and print an error message
  if(!inherits(yname,"character",which = FALSE))
    stop(" yname must be of class 'character' ")
  # Assign 'data.x1' to 'X1'
  X1 <- data.x1
  # Assign 'data.y1' to 'Y1'
  Y1 <- data.y1
  # Assign 'data.xt' to 'Xt'
  Xt <- data.xt
  # Assign 'data.yt' to 'Yt'
  Yt <- data.yt
  # Assign 'z' to 'Time'
  Time <- z
  # Combine 'X1', 'Y1', 'Xt', 'Yt', and 'Time' into a single data frame 'adjusted_data'
  adjusted_data <- data.frame("x1"=X1,"xt"=Xt,"y1"=Y1,"yt"=Yt,"time"=Time)
  # Fit a generalized additive mixed model (GAMM) to 'adjusted_data' with 'xt' as the response variable and 's(y1,k=4) + x1' as the predictor variables
  fit <- gamm4::gamm4(xt ~ s(y1,k=4) + x1, data = adjusted_data, random=~(1|time))
  # Get a summary of the fitted model. 's' is a list containing the parameter estimates, standard errors, t values, p values, R squared, and the smoothing parameter selection criterion of the model.
  s <- summary(fit$gam)
  # Calculate the residuals of the model. 'r' is a vector containing the residuals.
  r <- adjusted_data$xt - predict(fit$gam)
  # Calculate the RMSE of the residuals. 'rmr' is a numeric value representing the RMSE.
  rmr <- sqrt(mean(r^2))
  # Create a data frame 'result0' to store the results of the model
  # 'result0' contains the model's parameter estimates, standard errors, t values, p values, R squared, GCV, and RMSE.
  result0 <- data.frame(lhs="Xt", op = "~", rhs="X1", beta=s$p.coeff[2], se=s$se[2],
                        t=s$p.t[2], p=s$p.pv[2], x=xname, y=yname, r.sq=s$r.sq, gcv=s$sp.criterion[1], rmr=rmr)
  # Fit another GAMM to 'adjusted_data' with 'yt' as the response variable and 's(x1,k=4) + y1' as the predictor variables
  fit <- gamm4::gamm4(yt ~ s(x1,k=4) + y1, data = adjusted_data, random=~(1|time))
  # Get a summary of the fitted model.
  s <- summary(fit$gam)
  # Calculate the residuals of the model.
  r <- adjusted_data$yt - predict(fit$gam)
  # Calculate the RMSE of the residuals. 'rmr' is a numeric value representing the RMSE.
  rmr <- sqrt(mean(r^2))
  # Create a data frame 'result1' to store the results of the model. 'result1' contains the model's parameter estimates, standard errors, t values, p values, R squared, GCV, and RMSE.
  result1 <- data.frame(lhs="Yt", op = "~", rhs="Y1", beta=s$p.coeff[2], se=s$se[2],
                        t=s$p.t[2], p=s$p.pv[2], x=xname, y=yname, r.sq=s$r.sq, gcv=s$sp.criterion[1], rmr=rmr)
  # Combine 'result0' and 'result1' from two GAMMs into a single data frame 'output'
  output <- rbind(result0,result1)
  # # Assign 'output' to 'rst'. 'rst' is a data frame containing the results of both models.
  rst <- output
  # Replace name in the 'lhs' column of 'rst' with a certain name,in order to better organize the results
  # Replace 'Xt' in the 'lhs' column of 'rst' with 'xname_2'
  rst$lhs[which(rst$lhs=="Xt")] <- paste0(xname,"_2")
  rst$lhs[which(rst$lhs=="Yt")] <- paste0(yname,"_2")
  rst$rhs[which(rst$rhs=="Y1")] <- paste0(yname,"_1")
  rst$rhs[which(rst$rhs=="X1")] <- paste0(xname,"_1")
  # Create a new column 'op' in 'rst' and assign '   -->   ' to all its values
  rst$op <- '   -->   '
  # Concatenate the values in the 'lhs', 'op_new', and 'rhs' columns
  for (i in 1:nrow(rst)) {
    rst$lhs[i] <- paste0(rst$lhs[i],rst$op[i],rst$rhs[i])
  }
  # Remove the 'op', 'rhs', and 'edf' columns from 'rst'
  rst <- rst[,-c(2,3,4)]
  # Rename the 'lhs' column in 'rst' to 'Relation'
  names(rst)[names(rst) =="lhs"] <-"Relation"
  # Rename the 'p' column in 'rst' to 'pvalue'
  names(rst)[names(rst) =="p"] <-"pvalue"
  # Return 'rst'. 'rst' is a data frame containing the results of both models
  return(rst)
}
