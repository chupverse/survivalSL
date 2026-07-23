print.summary.libsl <- function(x, digits = max(4, getOption("digits") - 3), ...)
{
  if (!inherits(x, "summary.libsl"))
    stop("'x' must be a 'summary.libsl' object")
  
  cat("\nSuper Learner Survival Model Summary")
  cat("\n====================================\n")
  
  # Informations sur le modèle
  cat("\nModel Information:")
  cat("\n------------------")
  cat("\n  Library:", x$library)
  if (!is.null(x$pro.time))
    cat("\n  Prediction time (pro.time):", x$pro.time)
  if (!is.null(x$ROC.precision))
    cat("\n  ROC precision:", min(x$ROC.precision), "-", max(x$ROC.precision))
  
  # Performance Metrics
  cat("\n\nPerformance Metrics:")
  cat("\n--------------------\n")
  
  # Formatage du tableau des métriques
  metrics_df <- x$metrics
  
  # Afficher le tableau
  print(metrics_df, digits = digits, ...)
  
  cat("\n---\n")
  if (any(is.na(metrics_df)))
    cat("\nNA values indicate metrics that could not be computed")
  
  invisible(x)
}