#' Title Nonlinear Cross Lag Analysis based on Polynomial linear regression: autoregression
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
#' @return A dataframe containing the result of Nonlinear Cross Lag Analysis after polynomial regression: Xt~X1,Yt~Y1
#' @importFrom stats lm
#' @export
#' @examples
#' data(test_data1)
#' data(test_data2)
#' # Not adjusting for covariates
#' result <- clpm_poly_r(xname = "ASI",yname = "PWRI",
#'                       data1 = test_data1,data2 = test_data2,adjust = FALSE)
#' # Adjust for covariates
#' result_ad <- clpm_poly_r(xname = "ASI",yname = "PWRI",
#'                          data1 = test_data1,data2 = test_data2,
#'                          adjust = TRUE,adjust_name=c("HDL_C","LDL_C"))
clpm_poly_r <- function(xname, yname, data1,data2,adjust=FALSE,adjust_name=NULL){
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
  # Subset 'data1' to select 'xname' and assign the result to 'x1_data'
  x1_data <- subset(data1,select=xname)
  # Rename the column in 'x1_data' to "x1"
  colnames(x1_data) <- "x1"
  # Create a new column 'x1_2' in 'x1_data' that is the square of 'x1'
  x1_data$x1_2 <- x1_data$x1^2
  # Create a new column 'x1_3' in 'x1_data' that is the cube of 'x1'
  x1_data$x1_3 <- x1_data$x1^3
  # Subset 'data1' to select 'yname' and assign the result to 'y1_data'
  y1_data <- subset(data1,select=yname)
  # Rename the column in 'y1_data' to "y1"
  colnames(y1_data) <- "y1"
  # Create a new column 'y1_2' in 'y1_data' that is the square of 'y1'
  y1_data$y1_2 <- y1_data$y1^2
  # Create a new column 'y1_3' in 'y1_data' that is the cube of 'y1'
  y1_data$y1_3 <- y1_data$y1^3
  # If 'adjust' is TRUE,perform covariate adjustment to 'x1' and 'y1'
  if(adjust==T){
    # If 'adjust_name' is NULL, stop the function and display an error message
    if(is.null(adjust_name)){
      stop(" If adjust is 'T', you need to enter adjust_name")
    }else{
      # If "FALSE" is in the result of checking if 'adjust_name' is in the names of 'data1', stop the function and display an error message
      if("FALSE"%in%(adjust_name%in%names(data1)))
        stop("adjust_name must be contained in colname of data1 or data2")
      # Combine 'x1_data$x1', 'y1_data$y1', and the subset of 'data1' selected by 'adjust_name' into a new data frame 'd1'
      d1 <- cbind(x1_data$x1,y1_data$y1,subset(data1,select = adjust_name))
      # Rename the columns in 'd1' to 'xname', 'yname', and 'adjust_name'
      colnames(d1) <- c(xname,yname,adjust_name)
      # Adjust the target using 'data1', 'xname', 'yname', and 'adjust_name' and assign the result to 'my_list'
      my_list <- adjust_target(data1,xname,yname,adjust_name)
      # Replace 'x1_data$x1' with the first element in 'my_list'
      x1_data$x1 <- my_list[[1]]
      # Replace 'y1_data$y1' with the second element in 'my_list'
      y1_data$y1 <- my_list[[2]]
      # Combine 'x1_data$x1_2', 'y1_data$y1_2', and the subset of 'data1' selected by 'adjust_name' into a new data frame 'd1'
      d1 <- cbind(x1_data$x1_2,y1_data$y1_2,subset(data1,select = adjust_name))
      # Rename the columns in 'd1' to 'xname', 'yname', and 'adjust_name'
      colnames(d1) <- c(xname,yname,adjust_name)
      # Adjust the target using 'data1', 'xname', 'yname', and 'adjust_name' and assign the result to 'my_list'
      my_list <- adjust_target(data1,xname,yname,adjust_name)
      # Replace 'x1_data$x1_2' with the first element in 'my_list'
    x1_data$x1_2 <- my_list[[1]]
    # Replace 'y1_data$y1_2' with the second element in 'my_list'
    y1_data$y1_2 <- my_list[[2]]
    # Combine 'x1_data$x1_3', 'y1_data$y1_3', and the subset of 'data1' selected by 'adjust_name' into a new data frame 'd1'
    d1 <- cbind(x1_data$x1_3,y1_data$y1_3,subset(data1,select = adjust_name))
    # Rename the columns in 'd1' to 'xname', 'yname', and 'adjust_name'
    colnames(d1) <- c(xname,yname,adjust_name)
    # Adjust the target using 'data1', 'xname', 'yname', and 'adjust_name' and assign the result to 'my_list'
    my_list <- adjust_target(data1,xname,yname,adjust_name)
    # Replace 'x1_data$x1_3' with the first element in 'my_list'
    x1_data$x1_3 <- my_list[[1]]
    # Replace 'y1_data$y1_3' with the second element in 'my_list'
    y1_data$y1_3 <- my_list[[2]]
    }
  }
  # Subset 'data2' to select 'xname' and assign the result to 'xt_data'
  xt_data <- subset(data2,select=xname)
  # Rename the column in 'xt_data' to "xt"
  colnames(xt_data) <- "xt"
  # Create a new column 'xt_2' in 'xt_data' that is the square of 'xt'
  xt_data$xt_2 <- xt_data$xt^2
  # Create a new column 'xt_3' in 'xt_data' that is the cube of 'xt'
  xt_data$xt_3 <- xt_data$xt^3
  # Subset 'data2' to select 'yname' and assign the result to 'yt_data'
  yt_data <- subset(data2,select=yname)
  # Rename the column in 'yt_data' to "yt"
  colnames(yt_data) <- "yt"
  # Create a new column 'yt_2' in 'yt_data' that is the square of 'yt'
  yt_data$yt_2 <- yt_data$yt^2
  # Create a new column 'yt_3' in 'yt_data' that is the cube of 'yt'
  yt_data$yt_3 <- yt_data$yt^3
  # Check if adjustment is required
  if(adjust==T){
    # Check if adjustment name is provided
    if(is.null(adjust_name)){
      # Stop execution if adjustment name is not provided
      stop(" If adjust is 'T', you need to enter adjust_name")
    }else{
      # Check if adjustment name is in the column names of data1
      if("FALSE"%in%(adjust_name%in%names(data1)))
        # Stop execution if adjustment name is not in the column names of data1
        stop("adjust_name must be contained in colname of data1 or data2")
      # Combine xt_data, yt_data and a subset of data2 into a new dataframe
      d1 <- cbind(xt_data$xt,yt_data$yt,subset(data2,select = adjust_name))
      # Set the column names of the new dataframe
    colnames(d1) <- c(xname,yname,adjust_name)
    # Adjust the target variables in 'd1' using the function 'adjust_target' and assign the result to 'my_list'
    my_list <- adjust_target(d1,xname,yname,adjust_name)
    # Replace the 'xt' column in 'x1_data' with the first element of 'my_list'
    xt_data$xt <- my_list[[1]]
    # Replace the 'yt' column in 'y1_data' with the second element of 'my_list'
    yt_data$yt <- my_list[[2]]
    # Bind the columns 'xt_2', 'yt_2', and 'adjust_name' from 'data2' into a new data frame 'd1'
    # Repeat the above steps for xt_data$xt_2 and yt_data$yt_2
    d1 <- cbind(xt_data$xt_2,yt_data$yt_2,subset(data2,select = adjust_name))
    colnames(d1) <- c(xname,yname,adjust_name)
    my_list <- adjust_target(d1,xname,yname,adjust_name)
    xt_data$xt_2 <- my_list[[1]]
    yt_data$yt_2 <- my_list[[2]]
    # Repeat the above steps for xt_data$xt_3 and yt_data$yt_3
    d1 <- cbind(xt_data$xt_3,yt_data$yt_3,subset(data2,select = adjust_name))
    colnames(d1) <- c(xname,yname,adjust_name)
    my_list <- adjust_target(d1,xname,yname,adjust_name)
    xt_data$xt_3 <- my_list[[1]]
    yt_data$yt_3 <- my_list[[2]]
    }
  }
  # Combine x1_data, xt_data, y1_data, yt_data into a new dataframe
  adjusted_data <- cbind(x1_data,xt_data,y1_data,yt_data)
  # Fit a linear model with xt as the dependent variable
  fit <- stats::lm(xt~y1+y1_2+y1_3+x1,data = adjusted_data)
  # Get the summary of the fitted model
  s <- summary(fit)
  # Get the ANOVA table of the fitted model
  a <- anova(fit)
  # Create a dataframe to store the results
  result0 <- data.frame(lhs="Xt", op="~", rhs="X1",estimate=s$coefficients[5],
                        p=s$coefficients[20], x=xname, y=yname, r.sq=s$r.squared)
  # Fit another linear model with yt as the dependent variable
  fit <- stats::lm(yt~x1+x1_2+x1_3+y1,data = adjusted_data)
  s <- summary(fit)
  a <- anova(fit)
  # Create another dataframe to store the results
  result1 <- data.frame(lhs="Yt", op="~", rhs="Y1", estimate=s$coefficients[5],
                        p=s$coefficients[20], x=xname, y=yname, r.sq=s$r.squared)
  # Combine the two result dataframes
  output <- rbind(result0,result1)
  rst <- output
  # Replace "Xt" with "xname_2" and "Yt" with "yname_2" in the lhs column
  rst$lhs[which(rst$lhs=="Xt")] <- paste0(xname,"_2")
  rst$lhs[which(rst$lhs=="Yt")] <- paste0(yname,"_2")
  # Replace "Y1" with "yname_1" and "X1" with "xname_1" in the rhs column
  rst$rhs[which(rst$rhs=="Y1")] <- paste0(yname,"_1")
  rst$rhs[which(rst$rhs=="X1")] <- paste0(xname,"_1")
  # Replace the operation symbol with an arrow
  rst$op <- '   -->   '
  # Combine the lhs, operation symbol and rhs into a single string
  for (i in 1:nrow(rst)) {
    rst$lhs[i] <- paste0(rst$lhs[i],rst$op[i],rst$rhs[i])
  }
  # Remove unnecessary columns
  rst <- rst[,-c(2,3,4)]
  # To better display the results,rename the columns
  names(rst)[names(rst) =="lhs"] <-"Relation"
  names(rst)[names(rst) =="p"] <-"pvalue"
  # Return the final result
  return(rst)
}
