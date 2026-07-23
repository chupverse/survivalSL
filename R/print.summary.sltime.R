print.summary.sltime <- function(x, digits = max(4, getOption("digits") - 3), ...)
{
  if (!inherits(x, "summary.sltime"))
    stop("'x' must be a 'summary.sltime' object")
  
  cat("\nSuper Learner Time Survival Model Summary")
  cat("\n==========================================\n")
  
  cat("\nModel Information:")
  cat("\n------------------")
  cat("\n  Method:", x$method)
  if (!is.null(x$pro.time))
    cat("\n  Prediction time (pro.time):", x$pro.time)
  if (!is.null(x$ROC.precision))
    cat("\n  ROC precision:", min(x$ROC.precision), "-", max(x$ROC.precision))
  
  cat("\n\nPerformance Metrics:")
  cat("\n--------------------\n")
  
  print(x$metrics, digits = digits, ...)
  
  invisible(x)
}

