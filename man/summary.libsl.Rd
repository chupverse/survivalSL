\name{summary.libsl}
\alias{summary.libsl}
\alias{print.summary.libsl}

\title{
S3 method for Summarizing a Learner
}
\description{
Return predictive performances of a model or algorithm obtained by a library  of the class \code{libsl}.
}

\usage{
\method{summary}{libsl}(object, newdata=NULL, ROC.precision=seq(.01,.99,.01),
      digits=7, pro.time=NULL, ...)
}

\arguments{
  \item{object}{An object returned by a library of the class \code{libsl}.}
  \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}, the summary is performed on the same subjects of the training sample.}
  \item{ROC.precision}{An optional argument with the percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve. 0 (min) and 1 (max) are not allowed. By default, the precision is \code{seq(.01,.99,.01)}.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
   \item{pro.time}{This optional value of prognostic time represents the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument times. Not used for the following metrics: "ll", "ibs", and "ibll". Default value is the time at which half of the subjects are still at risk.}
  \item{...}{Additional arguments affecting the summary which are passed from \code{libsl} by default.}
}



\value{
\item{library}{The name of the library used for model training.}
\item{pro.time}{The prognostic time used for evaluation.}
\item{ROC.precision}{The precision values used for the ROC curve computation.}
\item{metrics}{A data frame containing the computed predictive performance metrics.}

}


\details{
The following metrics can be used : "bs" for the Brier score at the prognostic time \code{pro.time}, "p_ci" and "uno_ci" for the concordance index at the prognostic time \code{pro.time} (Pencina and Uno versions), "ll" for the log-likelihood, "ibs" for the integrated Brier score up to the last observed time in the training data, "ibll" for the integrated binomial log-likelihood up to the last observed time in the training data, "ribs" for the restricted integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted integrated binomial log-likelihood up to the prognostic time \code{pro.time}, "bll" for the binomial log-likelihood, "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.
}


\examples{
data("dataDIVAT2")

# The training of the gompertz  model with the first 400 patients

formula<-Surv(times,failures) ~ age + hla + retransplant + ecd

data<-dataDIVAT2[1:400,]
model <- LIB_PHgompertz(formula, data=data)

# The prognostic capacities from the same training sample
summary(model)

# The prognostic capacities from a validation of the next 150 patients
# (up to 4 years for several indicators)

#newdata<-dataDIVAT2[401:550,]
#summary(model, pro.time=4, newdata=newdata)


}

\keyword{Summarizing}

