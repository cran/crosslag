#' Title Cross-lag analysis based on generalized additive mixture models
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
#' @return A dataframe containing the result of generalized additive mixture model (consists of a smoothing term, a linear term and a stochastic intercept).
#' @importFrom gamm4 gamm4
#' @export
#' @examples
#' data(test_data1)
#' data(test_data2)
#' result <- clpm_gamm4_c("ASI","PWRI",test_data1$ASI,test_data1$PWRI,
#'                      data.xt = test_data2$ASI,data.yt = test_data2$PWRI,z=test_data1$time)
clpm_gamm4_c <- function(xname, yname, data.x1,data.y1,data.xt,data.yt,z){
  # Check if data.x1 is of class 'numeric'
  if(!inherits(data.x1,"numeric",which = FALSE))
    stop(" data.x1 must be of class 'numeric' ")
  # Check if data.y1 is of class 'numeric'
  if(!inherits(data.y1,"numeric",which = FALSE))
    stop(" data.y1 must be of class 'numeric' ")
  # Check if data.xt is of class 'numeric'
  if(!inherits(data.xt,"numeric",which = FALSE))
    stop(" data.xt must be of class 'numeric' ")
  # Check if data.yt is of class 'numeric'
  if(!inherits(data.yt,"numeric",which = FALSE))
    stop(" data.yt must be of class 'numeric' ")
  # Check if z is of class 'numeric'
  if(!inherits(z,"numeric",which = FALSE))
    stop(" z must be of class 'numeric' ")
  # Check if xname is of class 'character'
  if(!inherits(xname,"character",which = FALSE))
    stop(" xname must be of class 'character' ")
  # Check if yname is of class 'character'
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
  # Get a summary of the fitted model
  s <- summary(fit$gam)
  # Calculate the residuals of the model
  r <- adjusted_data$xt - predict(fit$gam)
  # Calculate the root mean square error (RMSE) of the residuals
  rmr <- sqrt(mean(r^2))
  # Calculate the chi-square statistic divided by the degrees of freedom
  chidf <- s$chi.sq/s$s.table[1]
  # Create a data frame 'result0' to store the results of the model
  result0 <- data.frame(lhs="Xt", op="~", rhs="Y1", edf=s$s.table[1], f=s$s.table[3],
                        p=s$s.table[4], x=xname, y=yname, r.sq=s$r.sq, gcv=s$sp.criterion[1],
                        rmr=rmr, chisqdf=chidf)
  # Fit another GAMM to 'adjusted_data' with 'yt' as the response variable and 's(x1,k=4) + y1' as the predictor variables
  fit <- gamm4::gamm4(yt ~ s(x1,k=4) + y1, data = adjusted_data, random=~(1|time))
  # Get a summary of the fitted model
  s <- summary(fit$gam)
  # Calculate the residuals of the model
  r <- adjusted_data$yt - predict(fit$gam)
  # Calculate the RMSE of the residuals
  rmr <- sqrt(mean(r^2))
  # Calculate the chi-square statistic divided by the degrees of freedom
  chidf <- s$chi.sq/s$s.table[1]
  # Create a data frame 'result1' to store the results of the model
  result1 <- data.frame(lhs="Yt", op="~", rhs="X1", edf=s$s.table[1], f=s$s.table[3],
                        p=s$s.table[4], x=xname, y=yname, r.sq=s$r.sq, gcv=s$sp.criterion[1],
                        rmr=rmr, chisqdf=chidf)
  fit <- gamm4::gamm4(xt ~ s(y1,k=4) + x1, data = adjusted_data, random=~(1|time))
  s <- summary(fit$gam)
  r <- adjusted_data$xt - predict(fit$gam)
  rmr <- sqrt(mean(r^2))
  # Combine 'result0' and 'result1' into a single data frame 'output'
  output <- rbind(result0,result1)
  # Assign 'output' to 'rst'
  rst <- output
  # Replace 'Xt' in the 'lhs' column of 'rst' with 'xname_2'
  rst$lhs[which(rst$lhs=="Xt")] <- paste0(xname,"_2")
  # Replace 'Yt' in the 'lhs' column of 'rst' with 'yname_2'
  rst$lhs[which(rst$lhs=="Yt")] <- paste0(yname,"_2")
  # Replace 'Y1' in the 'rhs' column of 'rst' with 'yname_1'
  rst$rhs[which(rst$rhs=="Y1")] <- paste0(yname,"_1")
  # Replace 'X1' in the 'rhs' column of 'rst' with 'xname_1'
  rst$rhs[which(rst$rhs=="X1")] <- paste0(xname,"_1")
  # Replace all values in the 'op' column of 'rst' with '   -->   '
  rst$op <- '   -->   '
  # For each row in 'rst', concatenate the values in the 'lhs', 'op', and 'rhs' columns
  for (i in 1:nrow(rst)) {
    rst$lhs[i] <- paste0(rst$lhs[i],rst$op[i],rst$rhs[i])
  }
  # Remove the 'op', 'rhs', and 'edf' columns from 'rst'
  rst <- rst[,-c(2,3,4)]
  # Rename the 'lhs' column in 'rst' to 'Relation'
  names(rst)[names(rst) =="lhs"] <-"Relation"
  # Rename the 'p' column in 'rst' to 'pvalue'
  names(rst)[names(rst) =="p"] <-"pvalue"
  # Return the final output of the function `clpm_gamm4_c`.
  return(rst)
}
