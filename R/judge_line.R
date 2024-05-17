#' Title Determine a linear relationship between x and y by Restricted Cubic Splines(RCS) and plot it
#' @param xname A character, the name of the x. If there are multiple covariates, they need to be saved in a vector
#' @param yname A character, the name of the y. If there are multiple covariates, they need to be saved in a vector
#' @param Data A dataframe containing data of x and y.
#' The input 'xname' and 'yname' is written in the same way among the 'Data' column names
#' @param parm the number of knots of RCS
#' @return A list containing the information of RCS plots.
#' You can use 'ggpubr::ggarrange()' to plot the list.
#' The list is a list of ggplot objects.
#' Each ggplot object is a plot that represents the relationship between one independent variable (from xname) and one dependent variable (from yname).
#' These plots are based on a linear regression model that uses a restricted cubic spline to handle the independent variable.
#' Each plot includes a line of predicted values (solid blue line) and a confidence interval for these predictions (transparent blue area).
#' The title of each plot includes the p-value for non-linearity, as well as the names of the independent and dependent variables.
#' @export
#' @import rms
#' @import ggplot2
#' @import ggpubr
#' @import lavaan
#' @import utils
#' @examples
#' data(test_data1)
#'
#' # Get RCS plots
#' result <- judge_line(xname="ASI",yname="HDL_C",Data = test_data1,parm = 4)
#' ## Plot an RCS curve graph
#' ggpubr::ggarrange(result[[1]])
#'
#' # Simultaneously determining multiple linear relationships
#' results <- judge_line("ASI",c("HDL_C","LDL_C"),Data = test_data1,parm = 4)
#' ## Plot many RCS curve graphs
#' ## Adjust 'nrow' and 'ncol' according to actual situation
#' ggpubr::ggarrange(plotlist = results, nrow = 1, ncol = 2)

judge_line <- function(xname,yname,Data,parm){
  # Check if 'xname' is a character, if not, stop the function and print an error message
  if(!inherits(xname,"character",which = FALSE))
    stop(" xname must be of class 'character' ")
  # Check if 'yname' is a character, if not, stop the function and print an error message
  if(!inherits(yname,"character",which = FALSE))
    stop(" yname must be of class 'character' ")
  # Check if 'xname' is in the column names of 'Data', if not, stop the function and print an error message
  if("FALSE"%in%(xname%in%names(Data)))
    stop("xname must be contained in colname of Data")
  # Check if 'yname' is in the column names of 'Data', if not, stop the function and print an error message
  if("FALSE"%in%(yname%in%names(Data)))
    stop("yname must be contained in colname of Data")
  # Check if 'Data' is a data frame, if not, stop the function and print an error message
  if(!inherits(Data,"data.frame",which = FALSE))
    stop(" Data must be of class 'data.frame' ")
  # Check if 'xname' has only one element
  if(length(xname)==1){
    # Check if 'yname' has only one element
    if(length(yname)==1){
      # Extract the columns from 'Data' that match 'xname' and 'yname'
      x <- Data[,which(names(Data)==xname)]
      y <- Data[,which(names(Data)==yname)]
      # Combine 'x' and 'y' into a new data frame
      data <- data.frame(x, y)
      # Compute the summary statistics of 'data'
      # 'datad' is an object of class "datadist" that contains these summary statistics.
      datad =rms::datadist(data)
      # Convert 'x' and 'y' to numeric if they are not
      if(!inherits(data$x,"numerics",which = FALSE)){data$x <- as.numeric(data$x)}
      if(!inherits(data$y,"numerics",which = FALSE)){data$y <- as.numeric(data$y)}
      # Set the 'datadist' option to 'datad'
      options(datadist = datad)
      # Fit a linear regression model with 'y' as the response variable and a restricted cubic spline of 'x' as the predictor
      fit <- rms::ols(data[,2] ~ rcs(x,parm),data = data)
      # Perform an ANOVA test on the fitted model
      an <- lavaan::anova(fit)
      # Compute the p-value of the test and format it
      p_value <- ifelse(an[2,5] < 0.05, "< 0.05", sprintf("%0.3f",an[2,5]))
      p_value <- ifelse(an[2,5] < 0.001, "< 0.001", sprintf("%0.3f",an[2,5]))
      # Extract the 'x' variable from 'data'
      x <- data[,1]
      # Compute the predicted values and confidence intervals of the fitted model
      OLS1 <- Predict(fit, x)
      yhat <- OLS1$yhat
      lower <- OLS1$lower
      upper <- OLS1$upper
      # Create a new variable name for the predicted values
      yname <- paste0('Predict_',yname)
      # Create a ggplot object with the predicted values and confidence intervals
      p <- ggplot2::ggplot() +
        geom_line(data=OLS1, aes(x, yhat), linetype = 'solid', linewidth = 1, alpha = 0.7, colour = "#0070b9") +
        geom_ribbon(data=OLS1, aes(x, ymin = lower, ymax = upper),alpha = 0.1,fill="#0070b9") +
        theme_classic() +
        labs(title = paste0("Non-linear p-value = ", p_value, "\n", xname, "~  ", yname), x = xname, y = yname)
      # Store the ggplot object in a list
      p_list <- list(p)
    }
  }
  # Check if 'xname' has more than one element
  if(length(xname)>1){
    # Check if 'yname' has only one element
    if(length(yname)==1){
      # Initialize an empty list 'p_list'
      p_list <- NULL
      # Subset 'Data' to keep only the columns that match 'xname'
      allx <- subset(Data,select = xname)
      # Loop over each element in 'xname'
      for (i in 1:length(xname)) {
        # Extract the i-th column from 'allx'
        x <- allx[,i]
        # Extract the column from 'Data' that matches 'yname'
        y <- Data[,which(names(Data)==yname)]
        # Combine 'x' and 'y' into a new data frame
        data <- data.frame(x, y)
        # Compute the summary statistics of 'data'
        datad = rms::datadist(data)
        # Convert 'x' and 'y' to numeric if they are not
        if(!inherits(data$x,"numerics",which = FALSE)){data$x <- as.numeric(data$x)}
        if(!inherits(data$y,"numerics",which = FALSE)){data$y <- as.numeric(data$y)}
        # Set the 'datadist' option to 'datad'
        options(datadist = datad)
        # Fit a linear regression model with 'y' as the response variable and a restricted cubic spline of 'x' as the predictor
        fit <- rms::ols(data[,2] ~ rcs(x,parm),data = data)
        # Perform an ANOVA test on the fitted model
        an <- lavaan::anova(fit)
        # Compute the p-value of the test and format it
        p_value <- ifelse(an[2,5] < 0.05, "< 0.05", sprintf("%0.3f",an[2,5]))
        p_value <- ifelse(an[2,5] < 0.001, "< 0.001", sprintf("%0.3f",an[2,5]))
        # Extract the 'x' variable from 'data'
        x <- data[,1]
        # Compute the predicted values and confidence intervals of the fitted model
        OLS1 <- Predict(fit, x)
        # Create a new variable name for the predicted values
        yn <- paste0('Predict_',yname)
        # Extract the i-th element from 'xname'
        xn <- xname[i]
        # Create a ggplot object with the predicted values and confidence intervals
        p <- ggplot2::ggplot() +
          geom_line(data=OLS1, aes(x, yhat), linetype = 'solid', linewidth = 1, alpha = 0.7, colour = "#0070b9") +
          geom_ribbon(data=OLS1, aes(x, ymin = lower, ymax = upper),alpha = 0.1,fill="#0070b9") +
          theme_classic() +
          labs(title = paste0("Non-linear p-value = ", p_value, "\n", xn, "~  ", yn), x = xn, y = yn)
        # Append the ggplot object to the list 'p_list'
        p_list <- c(p_list, list(p))
      }
    }
  }
  # Check if 'xname' has only one element
  if(length(xname)==1){
    # Check if 'yname' has more than one element
    if(length(yname)>1){
      # Initialize an empty list 'p_list'
      p_list <- NULL
      # Subset 'Data' to keep only the columns that match 'yname'
      ally <- subset(Data,select = yname)
      # Loop over each element in 'yname'
      for (i in 1:length(yname)) {
        # Extract the i-th column from 'ally'
        y <- ally[,i]
        # Extract the column from 'Data' that matches 'xname'
        x <- Data[,which(names(Data)==xname)]
        # Combine 'x' and 'y' into a new data frame
        data <- data.frame(x, y)
        # Compute the summary statistics of 'data'
        datad = rms::datadist(data)
        # Convert 'x' and 'y' to numeric if they are not
        if(!inherits(data$x,"numerics",which = FALSE)){data$x <- as.numeric(data$x)}
        if(!inherits(data$y,"numerics",which = FALSE)){data$y <- as.numeric(data$y)}
        # Set the 'datadist' option to 'datad'
        options(datadist = datad)
        # Fit a linear regression model with 'y' as the response variable and a restricted cubic spline of 'x' as the predictor
        fit <- rms::ols(data[,2] ~ rcs(x,parm),data = data)
        # Perform an ANOVA test on the fitted model
        an <- lavaan::anova(fit)
        # Compute the p-value of the test and format it
        p_value <- ifelse(an[2,5] < 0.05, "< 0.05", sprintf("%0.3f",an[2,5]))
        p_value <- ifelse(an[2,5] < 0.001, "< 0.001", sprintf("%0.3f",an[2,5]))
        # Extract the 'x' variable from 'data'
        x <- data[,1]
        # Compute the predicted values and confidence intervals of the fitted model
        OLS1 <- Predict(fit, x)
        # Create a new variable name for the predicted values
        yn <- paste0('Predict_',yname)
        # Extract the i-th element from 'xname'
        xn <- xname[i]
        # Create a ggplot object with the predicted values and confidence intervals
        p <- ggplot2::ggplot() +
          geom_line(data=OLS1, aes(x, yhat), linetype = 'solid', linewidth = 1, alpha = 0.7, colour = "#0070b9") +
          geom_ribbon(data=OLS1, aes(x, ymin = lower, ymax = upper),alpha = 0.1,fill="#0070b9") +
          theme_classic() +
          labs(title = paste0("Non-linear p-value = ", p_value, "\n", xn, "~  ", yn), x = xn, y = yn)
        # Append the ggplot object to the list 'p_list'
        p_list <- c(p_list, list(p))
      }
    }
  }
  # Check if both 'xname' and 'yname' have more than one element
  if(length(xname)>1){
    if(length(yname)>1){
      # Subset 'Data' to keep only the columns that match 'xname' and 'yname'
      allx <- subset(Data,select = xname)
      ally <- subset(Data,select = yname)
      # Initialize an empty list 'p_list'
      p_list <- NULL
      # Loop over each element in 'yname' and 'xname'
      for (i in 1:length(yname)) {
        for (j in 1:length(xname)) {
          # Extract the i-th column from 'ally' and j-th column from 'allx'
          y <- ally[,i]
          x <- allx[,j]
          # Combine 'x' and 'y' into a new data frame
          data <- data.frame(x, y)
          # Compute the summary statistics of 'data'
          datad = rms::datadist(data)
          # Convert 'x' and 'y' to numeric if they are not
          if(!inherits(data$x,"numerics",which = FALSE)){data$x <- as.numeric(data$x)}
          if(!inherits(data$y,"numerics",which = FALSE)){data$y <- as.numeric(data$y)}
          # Set the 'datadist' option to 'datad'
          options(datadist = datad)
          # Fit a linear regression model with 'y' as the response variable and a restricted cubic spline of 'x' as the predictor
          fit <- rms::ols(data[,2] ~ rcs(x,parm),data = data)
          # Perform an ANOVA test on the fitted model
          an <- lavaan::anova(fit)
          # Compute the p-value of the test and format it
          p_value <- ifelse(an[2,5] < 0.05, "< 0.05", sprintf("%0.3f",an[2,5]))
          p_value <- ifelse(an[2,5] < 0.001, "< 0.001", sprintf("%0.3f",an[2,5]))
          # Extract the 'x' variable from 'data'
          x <- data[,1]
          # Compute the predicted values and confidence intervals of the fitted model
          OLS1 <- Predict(fit, x)
          # Extract the j-th element from 'xname' and i-th element from 'yname'
          xn <- xname[j]
          yn <- paste0('Predict_',yname[i])
          # Create a ggplot object with the predicted values and confidence intervals
          p <- ggplot2::ggplot() +
            geom_line(data=OLS1, aes(x, yhat), linetype = 'solid', linewidth = 1, alpha = 0.7, colour = "#0070b9") +
            geom_ribbon(data=OLS1, aes(x, ymin = lower, ymax = upper),alpha = 0.1,fill="#0070b9") +
            theme_classic() +
            labs(title = paste0("Non-linear p-value = ", p_value, "\n", xn, "~  ", yn), x = xn, y = yn)
          # Append the ggplot object to the list 'p_list'
          p_list <- c(p_list, list(p))
        }
      }
    }
  }
  # Return the list 'p_list'
  return(p_list)
}
