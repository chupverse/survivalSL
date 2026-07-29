\name{tunePLANN}
\alias{tunePLANN}

\title{
Tune a Survival Neural Network Based on the PLANN Method
}

\description{
This function finds the optimal inter, size, decay, maxit, and MaxNWts parameters for the survival neural network.
}

\usage{
tunePLANN(formula, data, cv=10, inter=1, size=c(2, 4, 6, 8, 10),
decay=c(0.001, 0.01, 0.02, 0.05), maxit=100, MaxNWts=10000, maxtime=NULL,
seed=NULL,metric="auc", pro.time=NULL,
ROC.precision=seq(.01, .99, by=.01))
}

\arguments{
  \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame for training the model with the same covariates as in the formula.}
  \item{cv}{The value of the number of folds. The default value is 10.}
   \item{metric}{The loss function or metric. See details. Default metric is Area Under ROC (\code{"auc"}).}
    \item{inter}{The length of the intervals.}
  \item{size}{The number of units in the hidden layer.}
  \item{decay}{The parameter for weight decay.}
  \item{maxit}{The maximum number of iterations.}
  \item{MaxNWts}{The maximum allowable number of weights.}
  \item{maxtime}{A numeric value with the maximum prognostic. If \code{NULL}, the maximum prognostic time is the highest time observed in the data + 1.}
   \item{pro.time}{This optional value of prognostic time represents the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument times. Not used for the following metrics: "ll", "ibs", and "ibll". Default value is the time at which half of the subjects are still at risk.}
   \item{seed}{A random seed to ensure reproducibility during the cv process. If \code{NULL}, a seed is randomly assigned.}
   \item{ROC.precision}{The percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve. Only used when \code{metric="auc"}. 0 (min) and 1 (max) are not allowed. By default: \code{seq(.01,.99,.01)}.}
 }

 \details{
The function runs the \code{flexsurvspline} function of the \code{flexsurv} package. The following metrics can be used : "bs" for the Brier score at the prognostic time \code{pro.time}, "p_ci" and "uno_ci" for the concordance index at the prognostic time \code{pro.time} (Pencina and Uno versions), "ll" for the log-likelihood, "ibs" for the integrated Brier score up to the last observed time in the training data, "ibll" for the integrated binomial log-likelihood up to the last observed time in the training data, "ribs" for the restricted integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted integrated binomial log-likelihood up to the prognostic time \code{pro.time}, "bll" for the binomial log-likelihood, "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.
}


\value{
\item{optimal}{The value of \code{inter}, \code{size}, \code{decay}, \code{maxit}, and \code{MaxNWts} that optimizes the cross-validated value of the metric/loss function.}
\item{results}{The data frame with the cross-validated value of the metric/loss function according to \code{inter}, \code{size}, \code{decay}, \code{maxit}, and \code{MaxNWts}.}
}

\references{
Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.
}

\examples{
data("dataDIVAT2")

formula<-Surv(times,failures) ~ age + retransplant + ecd
tune.model <- tunePLANN(formula=formula, data=dataDIVAT2[1:125,],
cv=5, inter=c(2,1), size=c(16,32), decay=0.01, maxit=100, MaxNWts=1000,
maxtime=NULL,seed=123,metric="auc",
pro.time=NULL,ROC.precision=seq(.05, .95, by=.05))

tune.model$optimal # the optimal hyperparameters

tune.model$results
}

\keyword{Tunning}
