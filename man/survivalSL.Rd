\name{survivalSL}
\alias{survivalSL}

\title{
Super Learner for Censored Outcomes
}
\description{
This function allows to compute a Super Learner (SL) to predict survival outcomes.
}
\usage{
survivalSL(formula, data, methods, metric="auc", penalty=NULL,
cv=10, param.tune=NULL, pro.time=NULL,
optim.local.min=FALSE, ROC.precision=seq(.01,.99,.01),
param.weights.fix=NULL, param.weights.init=NULL,
seed=NULL, optim.method="Nelder-Mead", maxit=1000,
show_progress=TRUE)
}

\arguments{
  \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{data}{A data frame whose columns correspond to the variables present in the formula.}
  \item{methods}{A vector of characters with the names of the algorithms included in the SL. At least two algorithms have to be included.}
  \item{metric}{The loss function or metric used to estimate the weights of the algorithms in the SL. See details.}
  \item{penalty}{A numerical vector that allows the integration of covariates into the final model after selection (It concerns \code{"LIB_COXaic"}.) or/and allows the covariates not to be penalized (It concerns : \code{"LIB_COXen"} \code{"LIB_COXlasso"} and \code{"LIB_COXridge"}.). We give the value 0 if we want to force the covariate in the model or/and not to be penalized otherwise 1. If \code{NULL}, all covariates undergo the selection process or/and penalization process. }
 \item{cv}{The number of splits for cross-validation. The default value is 10.}
 \item{param.tune}{A list with a length equals to the number of algorithms included in \code{methods}. If \code{NULL}, the tunning parameters are estimated (see details).}
 \item{pro.time}{This optional value of prognostic time represents the maximum delay for which the capacity of the variable is evaluated. The same unit than the one used in the argument times. Not used for the following metrics: "ll", "ibs", and "ibll". Default value is the time at which half of the subjects are still at risk.}
 \item{optim.local.min}{An optional logical value. If \code{TRUE}, the optimization is performed twice to better ensure the estimation of the weights. If \code{FALSE} (default value), the optimization is performed once.}
 \item{ROC.precision}{The percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve. Only used when \code{metric="auc"}. 0 (min) and 1 (max) are not allowed. By default: \code{seq(.01,.99,.01)}.}
  \item{param.weights.fix}{A vector with the parameters of the multinomial logistic regression which generates the weights of the algorithms declared in \code{methods}. When completed, the related parameters are not estimated. The default value is NULL: the parameters are estimated by a \code{cv}-fold cross-validation. See details.}
   \item{param.weights.init}{A vector with the initial values of the parameters of the multinomial logistic regression which generates the weights of the algorithms declared in \code{methods}. The default value is NULL: the initial values are equaled to 0. See details.}
   \item{seed}{A random seed to ensure reproducibility. If \code{NULL}, a seed is randomly assigned.}
   \item{optim.method}{The optimization method used to estimate the weights. It can be either \code{"SANN"} or \code{"Nelder-Mead"}. By default we use \code{Nelder-Mead}.}
   \item{maxit}{The number of iterations during the weight optimization process. By default, it is set to 1000.}
   \item{show_progress}{Parameter to display the progress bar. By default, it is set to \code{TRUE}.}
 }


\details{
Each object of the list declared in \code{param.tune} must have the same name than the names of the \code{methods} included in the SL. If \code{param.tune} = \code{NULL}, \code{survivalSL} has already predefined default grids of tunning parameters for each algorithm in this case. The final tunning parameters are chosen thanks to \code{cv}-fold cross-validation (except for \code{LIB_RSF}, which uses the Out Of Bag observations to select the best hyperparameters based on the optimal value of the chosen metric). The following metrics can be used : "bs" for the Brier score at the prognostic time \code{pro.time}, "p_ci" and "uno_ci" for the concordance index at the prognostic time \code{pro.time} (Pencina and Uno versions), "ll" for the log-likelihood, "ibs" for the integrated Brier score up to the last observed time in the training data, "ibll" for the integrated binomial log-likelihood up to the last observed time in the training data, "ribs" for the restricted integrated Brier score up to the prognostic time \code{pro.time}, "ribll" for the restricted integrated binomial log-likelihood up to the prognostic time \code{pro.time}, "bll" for the binomial log-likelihood, "auc" for the area under the time-dependent ROC curve up to the prognostic time \code{pro.time}.


The following learners are available:
  \tabular{llll}{
  Names \tab Description \tab Package  \cr
  \code{"LIB_AFTgamma"} \tab Gamma-distributed AFT model \tab flexsurv  \cr
  \code{"LIB_AFTggamma"} \tab  Generalized Gamma-distributed AFT model \tab flexsurv  \cr
  \code{"LIB_AFTweibull"} \tab  Weibull-distributed AFT model \tab flexsurv  \cr
  \code{"LIB_PHexponential"} \tab  Exponential-distributed PH model \tab flexsurv  \cr
  \code{"LIB_PHgompertz"} \tab  Gompertz-distributed PH model \tab flexsurv  \cr
  \code{"LIB_PHspline"} \tab  Spline-based PH model \tab flexsurv  \cr
  \code{"LIB_COXall"} \tab  Usual Cox model \tab survival  \cr
  \code{"LIB_COXaic"} \tab  Cox model with AIC-based forward selection\tab MASS  \cr
  \code{"LIB_COXen"} \tab  Elastic Net Cox model \tab glmnet  \cr
  \code{"LIB_COXlasso"} \tab  Lasso Cox model \tab glmnet  \cr
  \code{"LIB_COXridge"} \tab   Ridge Cox model \tab glmnet  \cr
  \code{"LIB_RSF"} \tab Survival Random Forest \tab
  randomForestSRC  \cr
  \code{"LIB_PLANN"} \tab Survival Neural Network \tab survivalPLANN \cr}

The following loss functions for the estimation of the super learner weigths are available (\code{metric}):
\itemize{
  \item Area under the ROC curve (\code{"auc"})
  \item Pencina concordance index  (\code{"p_ci"})
  \item Uno concordance index (\code{"uno_ci"})
  \item Brier score (\code{"bs"})
  \item Binomial log-likelihood (\code{"bll"})
  \item Integrated Brier score (\code{"ibs"})
  \item Integrated binomial log-likelihood (\code{"ibll"})
  \item Restricted integrated Brier score (\code{"ribs"})
  \item Restricted integrated binomial log-Likelihood (\code{"ribll"})
  \item Log-Likelihood (\code{"ll"})

}
}

\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{It corresponds to a matrix with the survival predictions related to the SL.}
\item{FitALL}{It corresponds to a list of matrix with the survival predictions related to each of the learner used for the SL construction.}
\item{formula}{The formula object used for the SL construction.}
\item{data}{The data frame used for learning.}
\item{ROC.precision}{The percentiles (between 0 and 1) of the prognostic variable used for computing each point of the time dependent ROC curve.}
\item{cv}{The number of splits for cross-validation.}
\item{methods}{A vector of characters with the names of the algorithms included in the SL.}
\item{pro.time}{The maximum delay for which the capacity of the variable is evaluated.}
\item{models}{A list with the estimated models/algorithms included in the SL.}
\item{weights}{A list composed by two vectors: the regressions \code{coefficients} of the logistic multinomial regression and the resulting weights' \code{values}.}
\item{metric}{A list composed by two vectors: the loss function used to estimate the weights of the algorithms in the SL and its cross validation value.}
\item{param.tune}{The estimated tunning parameters.}
\item{seed}{The random seed used.}
\item{optim.method}{The optimization method used.}
}

\references{
Polley E and van der Laanet M. Super Learner In Prediction. http://biostats.bepress.com. 2010.
}

\examples{
data("dataDIVAT2")

# The Super Learner based from the first 125 individuals of the data base

formula<-Surv(times,failures) ~ age + retransplant + ecd

sl1 <- survivalSL(formula=formula, data=dataDIVAT2[1:125,],
                  methods=c("LIB_AFTgamma", "LIB_PHgompertz"))

# Individual prediction
pred <- predict(sl1, newdata=data.frame(age=c(52,52),
retransplant=c(1,1), ecd=c(0,1)))

plot(y=pred$predictions$sl[1,], x=pred$times, xlab="Time (years)",
ylab="Predicted survival", col=1, type="l", lty=1, lwd=2, ylim=c(0,1))

lines(y=pred$predictions$sl[2,], x=pred$times, col=2, type="l", lty=1, lwd=2)

legend("topright", col=c(1,2), lty=1, lwd=2, c("Subject #1", "Subject #2"))
}

\keyword{Super Learner}

