\name{summary.sltime}
\alias{summary.sltime}
\alias{print.summary.sltime}

\title{
S3 method for Summarizing a Super Learner
}
\description{
Return goodness-of-fit indicators of a Super Learner obtained by the function \code{survivalSL}.
}

\usage{
\method{summary}{sltime}(object, newdata=NULL,  method="sl",
ROC.precision=seq(.01,.99,.01), digits=7, pro.time=NULL, ...)
}

\arguments{
  \item{object}{An object returned by the function \code{survivalSL}.}
  \item{method}{A character string with the name of the algorithm included in the SL for which the summary is performed. The default is "sl" for the Super Learner.}
  \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}, the summary is performed on the same subjects of the training sample.}
   \item{ROC.precision}{An optional argument with the percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve. 0 (min) and 1 (max) are not allowed. By default, the precision is \code{seq(.01,.99,.01)}.}
   \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
    \item{pro.time}{This optional value of prognostic time represents the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument times. Not used for the following metrics: "ll", "ibs", and "ibll". Default value is the time at which half of the subjects are still at risk.}
     \item{...}{Additional arguments affecting the summary which are passed from \code{libsl} by default.}
}


\details{
The following metrics can be used : "bs" for the Brier score at the prognostic time \code{pro.time}, "p_ci" and "uno_ci" for the concordance index at the prognostic time \code{pro.time} (Pencina and Uno versions), "ll" for the log-likelihood, "ibs" for the integrated Brier score up to the last observed time in the training data, "ibll" for the integrated binomial log-likelihood up to the last observed time in the training data, "ribs" for the restricted integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted integrated binomial log-likelihood up to the prognostic time \code{pro.time}, "bll" for the binomial log-likelihood, "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.
}

\value{
\item{library}{The name of the library for which to test the performance}
\item{pro.time}{The prognostic time used for evaluation.}
\item{ROC.precision}{The precision values used for the ROC curve computation.}
\item{metrics}{A data frame containing the computed predictive performance metrics.}

}

\seealso{
  \code{\link{survivalSL}}.
}

\examples{
#data("dataDIVAT2")

#formula<-Surv(times,failures) ~ age + hla + retransplant + ecd

#sl1 <- survivalSL(formula, method=c("LIB_AFTgamma", "LIB_PHgompertz"
#,"LIB_AFTllogis"),  metric="auc",
  #data=dataDIVAT2[1:100,],
  #pro.time = 12,
  #cv=3)

# The prognostic capacities from the same training sample
#summary(sl1)
}

\keyword{Summarizing}
