\name{tunePHspline}
\alias{tunePHspline}

\title{
Tune a Survival Regression using the Royston/Parmar Spline Model
}

\description{
This function finds the optimal number of knots of the spline function.
}

\usage{
tunePHspline(formula,
data, cv=10, metric="auc", k=1:4, pro.time=NULL,
seed=NULL, ROC.precision=seq(.01, .99, by=.01))
}

\arguments{
  \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model with the same covariates as in the formula.}
  \item{cv}{The value of the number of folds. The default value is 10.}
   \item{metric}{The loss function or metric. See details. Default metric is Area Under ROC (\code{"auc"}).}
  \item{k}{The number of knots optimized over.}
   \item{pro.time}{This optional value of prognostic time represents the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument times. Not used for the following metrics: "ll", "ibs", and "ibll". Default value is the time at which half of the subjects are still at risk.}
   \item{seed}{A random seed to ensure reproducibility during the cv process. If \code{NULL}, a seed is randomly assigned.}
   \item{ROC.precision}{The percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve. Only used when \code{metric="auc"}. 0 (min) and 1 (max) are not allowed. By default: \code{seq(.01,.99,.01)}.}
 }

\details{
The function runs the \code{flexsurvspline} function of the \code{flexsurv} package. The following metrics can be used : "bs" for the Brier score at the prognostic time \code{pro.time}, "p_ci" and "uno_ci" for the concordance index at the prognostic time \code{pro.time} (Pencina and Uno versions), "ll" for the log-likelihood, "ibs" for the integrated Brier score up to the last observed time in the training data, "ibll" for the integrated binomial log-likelihood up to the last observed time in the training data, "ribs" for the restricted integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted integrated binomial log-likelihood up to the prognostic time \code{pro.time}, "bll" for the binomial log-likelihood, "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.
}

\value{
\item{optimal}{The value of \code{k} that optimizes the cross-validated value of the metric/loss function.}
\item{results}{The data frame with the cross-validated value of the metric/loss function according to \code{k}.}
}

\references{
Royston, P. and Parmar, M. (2002). Flexible parametric proportional-hazards and proportional odds models for censored survival data, with application to prognostic modelling and estimation of treatment effects. Statistics in Medicine 21(1):2175-2197. doi: 10.1002/sim.1203
}

\examples{
data("dataDIVAT2")

# The estimation of the hyperparameters on the first 100 patients

formula<-Surv(times,failures) ~ age + hla + retransplant
tune.model <- tunePHspline(formula=formula,
data=dataDIVAT2[1:100,], cv=3, metric="auc",
k=1:2, pro.time=NULL,seed=123,
ROC.precision=seq(.05, .95, by=.05))

# the estimated nodesize value

 tune.model$optimal
 tune.model$results
}

\keyword{Tunning}

