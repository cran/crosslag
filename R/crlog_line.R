#' Title Linear Cross lagged analysis
#'
#' @param xname If linear cross lagged analysis is used between x and y, 'xname' is the name of x
#' @param yname If linear cross lagged analysis is used between x and y, 'yname' is the name of y
#' @param data.x1 A numeric variable.
#' @param data.y1 A numeric variable. 'data.x1' and 'data.y1' comes from the first time point
#' @param data.xt A numeric variable.
#' @param data.yt A numeric variable. 'data.xt' and 'data.yt' comes from the another time point
#' @importFrom lavaan fitMeasures
#' @return a datafrme containing 'Relation','beta','se','z','pvalue','x','y','rmr' and 'cfi'.
#' Relation: This is a string that represents the relationship between the variables.
#' beta: This is the estimated coefficient of the model.
#' p: This is the p-value of the model coefficient, used to test whether the coefficient is significantly different from 0.
#' x: This is the input xname, representing the name of the x variable.
#' y: This is the input yname, representing the name of the y variable.
#' rmr: This is one of the model fit measures, representing the root mean square residual of the model.
#' cfi: This is one of the model fit measures, representing the comparative fit index of the model.
#' @export
#'
#' @examples
#' data(test_data1)
#' data(test_data2)
#' result <- crlog_line(xname="ASI",yname = "PWRI",data.x1 = test_data1$ASI,
#'                      data.y1 = test_data1$PWRI,data.xt = test_data2$ASI,
#'                      data.yt = test_data2$PWRI)
crlog_line <- function(xname,yname,data.x1,data.y1,data.xt,data.yt){
  # Check if data.x1 is numeric, if not, stop execution
  if(!inherits(data.x1,"numeric",which = FALSE))
    stop(" data.x1 must be of class 'numeric' ")
  # Check if data.y1 is numeric, if not, stop execution
  if(!inherits(data.y1,"numeric",which = FALSE))
    stop(" data.y1 must be of class 'numeric' ")
  # Check if data.xt is numeric, if not, stop execution
  if(!inherits(data.yt,"numeric",which = FALSE))
    stop(" data.xt must be of class 'numeric' ")
  # Check if data.yt is numeric, if not, stop execution
  if(!inherits(data.yt,"numeric",which = FALSE))
    stop(" data.yt must be of class 'numeric' ")
  # Check if xname is a character, if not, stop execution
  if(!inherits(yname,"character",which = FALSE))
    stop(" xname must be of class 'character' ")
  # Check if yname is a character, if not, stop execution
  if(!inherits(yname,"character",which = FALSE))
    stop(" yname must be of class 'character' ")
  # Define a Cross-Lagged Panel Model (CLPM)
  CLPM <- '
  # Estimate the lagged effects between the observed variables.
  fxt + fyt ~ fx1 + fy1

  # Estimate the covariance between the observed variables at the first wave.
  fx1 ~~ fy1 # Covariance

  # Estimate the covariances between the residuals of the observed variables.
  fxt ~~ fyt

  # Estimate the (residual) variance of the observed variables.
  fx1 ~~ fx1 # Variances
  fy1 ~~ fy1
  fxt ~~ fxt
  fyt ~~ fyt
'
      # Assign data.x1 to X1
      X1 <- data.x1
      Y1 <- data.y1
      Xt <- data.xt
      Yt <- data.yt
      # Create a data frame with columns fx1, fy1, fxt, fyt
      df <- data.frame(fx1=X1,fy1=Y1,fxt=Xt,fyt=Yt)
      # Fit the CLPM to the data
      fit <- lavaan(CLPM, data = df)
      # Extract two model fit measures
      temp <-  as.numeric(lavaan::fitMeasures(fit,c('rmr','cfi')))
      # Get the summary of the fitted model
      df <- summary(fit)
      # Get the p-values from the summary
      df <- df$p
      # Add xname and yname to the data frame
      df$x <- xname
      df$y <- yname
      # Add the model fit measures to the data frame
      df$rmr <- temp[1]
      df$cfi <- temp[2]
      # Add an operation symbol to the data frame
      df$op <- '   -->   '
      # Replace "fxt" with "xname_2" and "fyt" with "yname_2" in the lhs column
      df$lhs[which(df$lhs=="fxt")] <- paste0(xname,"_2")
      df$lhs[which(df$lhs=="fyt")] <- paste0(yname,"_2")
      # Replace "fx1" with "xname_1" and "fy1" with "yname_1" in the lhs column
      df$lhs[which(df$lhs=="fx1")] <- paste0(xname,"_1")
      df$lhs[which(df$lhs=="fy1")] <- paste0(yname,"_1")
      # Replace "fxt" with "xname_2" and "fyt" with "yname_2" in the rhs column
      df$rhs[which(df$rhs=="fxt")] <- paste0(xname,"_2")
      df$rhs[which(df$rhs=="fyt")] <- paste0(yname,"_2")
      # Replace "fx1" with "xname_1" and "fy1" with "yname_1" in the rhs column
      df$rhs[which(df$rhs=="fx1")] <- paste0(xname,"_1")
      df$rhs[which(df$rhs=="fy1")] <- paste0(yname,"_1")
      # Combine the lhs, operation symbol and rhs into a single string
      for (i in 1:nrow(df)) {
        df$lhs[i] <- paste0(df$lhs[i],df$op[i],df$rhs[i])
      }
      # Remove unnecessary columns
      df <- df[,-c(2,3,4)]
      # To better display the results,rename the columns
      names(df)[names(df) =="lhs"] <-"Relation"
      names(df)[names(df) =="est"] <-"beta"
      # Return the final result
      return(df)
}


