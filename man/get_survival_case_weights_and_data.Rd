% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/static_glm.R
\name{get_survival_case_weights_and_data}
\alias{get_survival_case_weights_and_data}
\title{Get data.frame for Discrete Time Survival Models}
\usage{
get_survival_case_weights_and_data(
  formula,
  data,
  by,
  max_T,
  id,
  init_weights,
  risk_obj,
  use_weights = TRUE,
  is_for_discrete_model = TRUE,
  c_outcome = "Y",
  c_weights = "weights",
  c_end_t = "t"
)
}
\arguments{
\item{formula}{\code{\link[survival]{coxph}} like formula with \code{\link[survival]{Surv}(tstart, tstop, event)} on the left hand site of \code{~}.}

\item{data}{\code{data.frame} or environment containing the outcome and covariates.}

\item{by}{interval length of the bins in which parameters are fixed.}

\item{max_T}{end of the last interval interval.}

\item{id}{vector of ids for each row of the in the design matrix.}

\item{init_weights}{weights for the rows in \code{data}. Useful e.g., with skewed sampling.}

\item{risk_obj}{a pre-computed result from a \code{\link{get_risk_obj}}. Will be used to skip some computations.}

\item{use_weights}{\code{TRUE} if weights should be used. See details.}

\item{is_for_discrete_model}{\code{TRUE} if the model is for a discrete hazard model is used like the logistic model.}

\item{c_outcome, c_weights, c_end_t}{alternative names to use for the added columns described in the return section. Useful if you already have a column named \code{Y}, \code{t} or \code{weights}.}
}
\value{
Returns a \code{data.frame} where the following is added (column names will differ if you specified them): column \code{Y} for the binary outcome, column \code{weights} for weights of each row and additional rows if applicable. A column \code{t} is added for the stop time of the bin if \code{use_weights = FALSE}. An element \code{Y} with the used \code{Surv} object is added if \code{is_for_discrete_model = FALSE}.
}
\description{
Function used to get \code{data.frame} with weights for a static fit for survivals.
}
\details{
This function is used to get the \code{data.frame} for e.g. a \code{glm} fit that is comparable to a \code{\link{ddhazard}} fit in the sense that it is a static version. For example, say that we bin our time periods into \code{(0,1]}, \code{(1,2]} and \code{(2,3]}. Next, consider an individual who dies at time 2.5. He should be a control in the the first two bins and should be a case in the last bin. Thus the rows in the final data frame for this individual is \code{c(Y = 1, ..., weights = 1)} and \code{c(Y = 0, ..., weights = 2)} where \code{Y} is the outcome, \code{...} is the covariates and \code{weights} is the weights for the regression. Consider another individual who does not die and we observe him for all three periods. Thus, he will yield one row with \code{c(Y = 0, ..., weights = 3)}.

This function use similar logic as the \code{ddhazard} for individuals with time varying covariates (see the vignette \code{vignette("ddhazard", "dynamichazard")} for details).

If \code{use_weights = FALSE} then the two previously mentioned individuals will yield three rows each. The first individual will have \code{c(Y = 0, t = 1, ..., weights = 1)}, \code{c(Y = 0, t = 2, ..., weights = 1)}, \code{c(Y = 1, t = 3, ..., weights = 1)} while the latter will have three rows \code{c(Y = 0, t = 1, ..., weights = 1)}, \code{c(Y = 0, t = 2, ..., weights = 1)}, \code{c(Y = 0, t = 3, ..., weights = 1)}. This kind of data frame is useful if you want to make a fit with e.g. \code{\link[mgcv]{gam}} function in the \code{mgcv} package as described en Tutz et. al (2016).
}
\examples{
library(dynamichazard)
# small toy example with time-varying covariates
dat <- data.frame(
 id     = c(   1,    1, 2,     2),
 tstart = c(   0,    4, 0,     2),
 tstop  = c(   4,    6, 2,     6),
 event  = c(   0,    1, 0,     0),
 x1     = c(1.09, 1.29, 0, -1.16))

get_survival_case_weights_and_data(
 Surv(tstart, tstop, event) ~ x1, dat, by = 1, id = dat$id)$X
get_survival_case_weights_and_data(
 Surv(tstart, tstop, event) ~ x1, dat, by = 1, id = dat$id,
 use_weights = FALSE)$X

}
\references{
Tutz, Gerhard, and Matthias Schmid. \emph{Nonparametric Modeling and Smooth Effects}. Modeling Discrete Time-to-Event Data. Springer International Publishing, 2016. 105-127.
}
\seealso{
\code{\link{ddhazard}}, \code{\link{static_glm}}
}
