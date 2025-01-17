% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clpm_gam_r.R
\name{clpm_gam_r}
\alias{clpm_gam_r}
\title{Title Nonlinear Cross Lag Analysis: autoregression}
\usage{
clpm_gam_r(xname, yname, data.x1, data.y1, data.xt, data.yt)
}
\arguments{
\item{xname}{If cross lagged analysis is used between x and y, 'xname' is the name of x}

\item{yname}{If cross lagged analysis is used between x and y, 'yname' is the name of y}

\item{data.x1}{A numeric variable.}

\item{data.y1}{A numeric variable. 'data.x1' and 'data.y1' comes from the first time point}

\item{data.xt}{A numeric variable.}

\item{data.yt}{A numeric variable. 'data.xt' and 'data.yt' comes from the another time point}
}
\value{
A dataframe containing the result of autoregression: 'Xt~X1' 'Yt~Y1'
}
\description{
Title Nonlinear Cross Lag Analysis: autoregression
}
\examples{
data(test_data1)
data(test_data2)
clpm_gam_r(xname="ASI",yname = "PWRI",data.x1 = test_data1$ASI,
           data.y1 = test_data1$PWRI,data.xt = test_data2$ASI,data.yt = test_data2$PWRI)
}
