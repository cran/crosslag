#' Title Generate covariates' residuals
#'
#' @param X A dataset containing y1, y2, and covariates. The row name is the name of x,y and covariate
#' @param y1 If cross lagged analysis is used between x and y, 'y1' is the name of x.
#' @param y2 y2 can be NULL. If cross lagged analysis is used between x and y, 'y2' is the name of y. 'y1' and 'y2' come from the same time point.
#' @param cor_vars  the name of covariate
#'
#' @return a list containing the covariate-adjusted value of y1 and y2 (If y2 is not NULL)
#' @importFrom stats lm
#' @export
#'
#' @examples
#' data(test_data1)
#' result <- adjust_target(X=test_data1,y1="ASI",y2=NULL,cor_vars=c("HDL_C","LDL_C"))
adjust_target <- function(X,y1,y2=NULL,cor_vars){
  #validate the input parameters and ensure they meet the required criteria for further processing within the function
  if(inherits(X,"matrix",which = FALSE)){X <- as.data.frame(X)}
  if("FALSE"%in%(y1%in%names(X)))
    stop(" y1 must be in colname of X")
  if(!inherits(y1,"character",which = FALSE))
    stop(" y1 must be of class 'character' ")
  if(!inherits(cor_vars,"character",which = FALSE))
    stop(" cor_vars must be of class 'character' ")
  if("FALSE"%in%(cor_vars%in%names(X)))
     stop(" cor_vars must be in colname of X")
  #creating a formula string for linear regression modeling
  my_formula <- paste0(y1,'~',paste0(cor_vars,collapse = '+'))
  # fitting a linear regression model
  fit <- lm(formula = stats::as.formula(my_formula),X)
  # assigning the residuals of the linear regression model
  y1_value <- fit$residuals
  # Perform two operations on the variable `y1_value`: Standardize and convert into numerical variables
  y1_value <- as.numeric(scale(y1_value))
  #This block of code is responsible for handling the case when the input variable `y2` is not NULL
  if(!is.null(y2)){
    #checking the validity of the input variable `y2`
    if(!inherits(y2,"character",which = FALSE))
      stop(" y2 must be of class 'character' ")
    if(!y2%in%colnames(X))
      stop(" y2 must be in colname of X")
  my_formula <- paste0(y2,'~',paste0(cor_vars,collapse = '+'))
  fit <- stats::lm(formula = stats::as.formula(my_formula),X)
  # Assign the residuals of the linear regression model
  # The residuals represent the differences between the observed values of the dependent variable and the values predicted by the regression model.
  # These residuals are a measure of how well the model fits the data points.
  y2_value <- fit$residuals
  # Perform two operations on the variable `y2_value`: Standardize and convert into numerical variables
  y2_value <- as.numeric(scale(y2_value))
  #creating a list named `result` with two elements.
  result <- list(y1=y1_value,y2=y2_value)
  # Names being assigned are `y1` and `y2`, which are the names of the variables for which residuals are calculated.
  # This allows for easy identification and access of the elements in the list using these specific names.
  names(result) <- c(y1,y2)
  # Return the final output of the function `adjust_target`.`result` is a list containing the residuals calculated for the dependent variables `y1` and `y2` based on the linear regression
  return(result)
  }else{
    # Return the final output of the function `adjust_target`.`y1_result` is the residuals calculated for the variables `y1` based on the linear regression
    return(y1_value)
  }
}
