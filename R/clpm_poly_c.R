#' Title Nonlinear Cross Lag Analysis based on Polynomial linear regression
#'
#' @param xname If cross lagged analysis is used between x and y, 'xname' is the name of x
#' @param yname If cross lagged analysis is used between x and y, 'xname' is the name of x
#' @param data1 A dataframe containing data of x and y at the same time point.
#' If covariate adjustments are made (adjust=T), data of covariates should also be included.
#' The input 'xname','yname' and 'adjust_name' is written in the same way among the 'data1' column names
#' @param data2 A dataframe containing data of x and y at the another time point.
#' If covariate adjustments are made (adjust=T), data of covariates should also be included.
#' The input 'xname','yname' and 'adjust_name' is written in the same way among the 'data2' column names
#' @param adjust The default value is F.
#' If you want to perform covariate adjustments on x and y based on adjust_target(), please use T.
#' @param adjust_name  The name of the covariate, corresponding to the data contained in data1 or data2.
#' If there are multiple covariates, they need to be saved in a vector.
#' @return A dataframe containing the result of Nonlinear Cross Lag Analysis after polynomial regression
#' @importFrom stats lm
#' @export
#' @examples
#' data(test_data1)
#' data(test_data2)
#' # Not adjusting for covariates
#' result <- clpm_poly_c(xname = "PWRI",yname = "LDL_C",
#'                       data1 = test_data1,data2 = test_data2,adjust = FALSE)
#' # Adjust for covariates
#' result_ad <- clpm_poly_c(xname = "ASI",yname = "PWRI",
#'                          data1 = test_data1,data2 = test_data2,
#'                          adjust = TRUE,adjust_name=c("HDL_C","LDL_C"))
clpm_poly_c <- function(xname, yname, data1,data2,adjust=FALSE,adjust_name=NULL){
  # Check if 'xname' is of class 'character'. If not, stop the function and return an error message
  if(!inherits(xname,"character",which = FALSE))
    stop(" xname must be of class 'character' ")
  # Check if 'xname' is in the names of 'data1'. If not, stop the function and return an error message
  if(!xname%in%names(data1))
    stop(" xname must be in names(data1)")
  # Check if 'xname' is in the names of 'data2'. If not, stop the function and return an error message
  if(!xname%in%names(data2))
    stop(" xname must be in names(data2)")
  # Check if 'yname' is of class 'character'. If not, stop the function and return an error message
  if(!inherits(yname,"character",which = FALSE))
    stop(" yname must be of class 'character' ")
  # Check if 'yname' is in the names of 'data1'. If not, stop the function and return an error message
  if(!yname%in%names(data1))
    stop(" yname must be in names(data1)")
  # Check if 'yname' is in the names of 'data2'. If not, stop the function and return an error message
  if(!yname%in%names(data2))
    stop(" yname must be in names(data2)")
  # Check if 'data1' is of class 'data.frame'. If not, stop the function and return an error message
  if(!inherits(data1,"data.frame",which = FALSE))
    stop(" data1 must be of class 'data.frame' ")
  # Check if 'data2' is of class 'data.frame'. If not, stop the function and return an error message
  if(!inherits(data2,"data.frame",which = FALSE))
    stop(" data2 must be of class 'data.frame' ")
  # Subset 'data1' to select the column named 'xname' and assign it to 'x1_data'
  x1_data <- subset(data1,select=xname)
  # Rename the column in 'x1_data' to 'x1'
  colnames(x1_data) <- "x1"
  # Create a new column 'x1_2' in 'x1_data' that is the square of 'x1'
  x1_data$x1_2 <- x1_data$x1^2
  # Create a new column 'x1_3' in 'x1_data' that is the cube of 'x1'
  x1_data$x1_3 <- x1_data$x1^3
  # Subset 'data1' to select the column named 'yname' and assign it to 'y1_data'
  y1_data <- subset(data1,select=yname)
  # Rename the column in 'y1_data' to 'y1'
  colnames(y1_data) <- "y1"
  # Create a new column 'y1_2' in 'y1_data' that is the square of 'y1'
  y1_data$y1_2 <- y1_data$y1^2
  # Create a new column 'y1_3' in 'y1_data' that is the cube of 'y1'
  y1_data$y1_3 <- y1_data$y1^3
  # If 'adjust' is TRUE, perform covariate adjustment
  if(adjust==T){
    # If 'adjust_name' is NULL, stop the function and return an error message
    if(is.null(adjust_name)){
      stop(" If adjust is 'T', you need to enter adjust_name")
    }else{
      if("FALSE"%in%(adjust_name%in%names(data1)))
        stop("adjust_name must be contained in colname of data1 or data2")
      # Bind the columns 'x1', 'y1', and 'adjust_name' from 'data1' into a new data frame 'd1'
      d1 <- cbind(x1_data$x1,y1_data$y1,subset(data1,select = adjust_name))
      # Rename the columns in 'd1' to 'xname', 'yname', and 'adjust_name'
    colnames(d1) <- c(xname,yname,adjust_name)
    # Adjust the target variables in 'data1' using the function 'adjust_target' and assign the result to 'my_list'
    my_list <- adjust_target(data1,xname,yname,adjust_name)
    # Replace the 'x1' column in 'x1_data' with the first element of 'my_list'
    x1_data$x1 <- my_list[[1]]
    # Replace the 'y1' column in 'y1_data' with the second element of 'my_list'
    y1_data$y1 <- my_list[[2]]
    # Bind the columns 'x1_2', 'y1_2', and 'adjust_name' from 'data1' into a new data frame 'd1'
    d1 <- cbind(x1_data$x1_2,y1_data$y1_2,subset(data1,select = adjust_name))
    # Rename the columns in 'd1' to 'xname', 'yname', and 'adjust_name'
    colnames(d1) <- c(xname,yname,adjust_name)
    # Adjust the target variables in 'data1' using the function 'adjust_target' and assign the result to 'my_list'
    my_list <- adjust_target(data1,xname,yname,adjust_name)
    # Replace the 'x1_2' column in 'x1_data' with the first element of 'my_list'
    x1_data$x1_2 <- my_list[[1]]
    # Replace the 'y1_2' column in 'y1_data' with the second element of 'my_list'
    y1_data$y1_2 <- my_list[[2]]
    # Bind the columns 'x1_3', 'y1_3', and 'adjust_name' from 'data1' into a new data frame 'd1'
    d1 <- cbind(x1_data$x1_3,y1_data$y1_3,subset(data1,select = adjust_name))
    # Rename the columns in 'd1' to 'xname', 'yname', and 'adjust_name'
    colnames(d1) <- c(xname,yname,adjust_name)
    # Adjust the target variables in 'data1' using the function 'adjust_target' and assign the result to 'my_list'
    my_list <- adjust_target(data1,xname,yname,adjust_name)
    # Replace the 'x1_3' column in 'x1_data' with the first element of 'my_list'
    x1_data$x1_3 <- my_list[[1]]
    # Replace the 'y1_3' column in 'y1_data' with the second element of 'my_list'
    y1_data$y1_3 <- my_list[[2]]
    }
  }
  # Subset 'data2' to select the column named 'xname' and assign it to 'xt_data'
  xt_data <- subset(data2,select=xname)
  # Rename the column in 'xt_data' to 'xt'
  colnames(xt_data) <- "xt"
  # Create a new column 'xt_2' in 'xt_data' that is the square of 'xt'
  xt_data$xt_2 <- xt_data$xt^2
  # Create a new column 'xt_3' in 'xt_data' that is the cube of 'xt'
  xt_data$xt_3 <- xt_data$xt^3
  # Subset 'data2' to select the column named 'yname' and assign it to 'yt_data'
  yt_data <- subset(data2,select=yname)
  # Rename the column in 'yt_data' to 'yt'
  colnames(yt_data) <- "yt"
  # Create a new column 'yt_2' in 'yt_data' that is the square of 'yt'
  yt_data$yt_2 <- yt_data$yt^2
  # Create a new column 'yt_3' in 'yt_data' that is the cube of 'yt'
  yt_data$yt_3 <- yt_data$yt^3
  # If 'adjust' is TRUE, perform covariate adjustment to 'xt' and 'yt'
  # The overall purpose is similar to the previous part
  if(adjust==T){
    if(is.null(adjust_name)){
      stop(" If adjust is 'T', you need to enter adjust_name")
    }else{
      if("FALSE"%in%(adjust_name%in%names(data1)))
        stop("adjust_name must be contained in colname of data1 or data2")
    d1 <- cbind(xt_data$xt,yt_data$yt,subset(data2,select = adjust_name))
    colnames(d1) <- c(xname,yname,adjust_name)
    my_list <- adjust_target(d1,xname,yname,adjust_name)
    xt_data$xt <- my_list[[1]]
    yt_data$yt <- my_list[[2]]
    d1 <- cbind(xt_data$xt_2,yt_data$yt_2,subset(data2,select = adjust_name))
    colnames(d1) <- c(xname,yname,adjust_name)
    my_list <- adjust_target(d1,xname,yname,adjust_name)
    xt_data$xt_2 <- my_list[[1]]
    yt_data$yt_2 <- my_list[[2]]
    d1 <- cbind(xt_data$xt_3,yt_data$yt_3,subset(data2,select = adjust_name))

    colnames(d1) <- c(xname,yname,adjust_name)
    my_list <- adjust_target(d1,xname,yname,adjust_name)
    xt_data$xt_3 <- my_list[[1]]
    yt_data$yt_3 <- my_list[[2]]
    }
  }
  # Combine 'x1_data', 'xt_data', 'y1_data', and 'yt_data' into a new data frame 'adjusted_data'
  adjusted_data <- cbind(x1_data,xt_data,y1_data,yt_data)
  # Fit a linear model with 'yt' as the response variable and 'x1', 'x1_2', 'x1_3', and 'y1' as the predictor variables
  fit <- stats::lm(xt~y1+y1_2+y1_3+x1, data = adjusted_data)
  # Get a summary of the fit and assign it to 's'
  s <- summary(fit)
  # Create a data frame 'result1' that contains the coefficients of the fit, the names of the predictor variables, and the R-squared value of the fit
  result0 <- data.frame(lhs="Xt", op="~", rhs="Y1", constant=s$coefficients[1],
                        primary=s$coefficients[2], quadratic=s$coefficients[3],
                        cubic=s$coefficients[4],p1=s$coefficients[17],
                        p2=s$coefficients[18],p3=s$coefficients[19],
                        x=xname, y=yname, r.sq=s$r.squared)
  fit <- stats::lm(yt~x1+x1_2+x1_3+y1, data = adjusted_data)
  s <- summary(fit)
  result1 <- data.frame(lhs="Yt", op="~", rhs="X1", constant=s$coefficients[1],
                        primary=s$coefficients[2], quadratic=s$coefficients[3],
                        cubic=s$coefficients[4],p1=s$coefficients[17],
                        p2=s$coefficients[18],p3=s$coefficients[19],
                        x=xname, y=yname, r.sq=s$r.squared)
  # Combine 'result0' and 'result1' containing the result of two linear model into a new data frame 'output'
  output <- rbind(result0,result1)
  # Replace the entries in the 'lhs' column of 'rst' with the certain name
  rst <- output
  rst$lhs[which(rst$lhs=="Xt")] <- paste0(xname,"_2")
  rst$lhs[which(rst$lhs=="Yt")] <- paste0(yname,"_2")
  rst$rhs[which(rst$rhs=="Y1")] <- paste0(yname,"_1")
  rst$rhs[which(rst$rhs=="X1")] <- paste0(xname,"_1")
  # Replace the entries in the 'op' column of 'rst' with '   -->   '
  rst$op <- '   -->   '
  # concatenate the values in the 'lhs', 'op', and 'rhs' columns and assign the result to the 'lhs' column
  for (i in 1:nrow(rst)) {
    rst$lhs[i] <- paste0(rst$lhs[i],rst$op[i],rst$rhs[i])
  }
  # Remove the 'op', 'rhs', and 'constant' columns from 'rst'
  rst <- rst[,-c(2,3,4)]
  # Rename the 'lhs' column in 'rst' to 'Relation'
  names(rst)[names(rst) =="lhs"] <-"Relation"
  # Rename the 'p1','p2','p3' column in 'rst' to 'pvalue_1','pvalue_2','pvalue_3'
  names(rst)[names(rst) =="p1"] <-"pvalue_1"
  names(rst)[names(rst) =="p2"] <-"pvalue_2"
  names(rst)[names(rst) =="p3"] <-"pvalue_3"
  # Return 'rst'. 'rst' is a data frame containing the results of both models
  return(rst)
}

