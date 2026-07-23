tuneCOXlasso<- function(formula, data, penalty = NULL, cv = 10, parallel =
                          FALSE, lambda=NULL, seed = NULL){


  if(is.null(seed)){
    seed<-sample(1:1000,1)
  }

  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(lambda)) stop("The 'lambda' argument is required.")



  variables_formula <- all.vars(formula)

  times <- variables_formula[1]
  failures <- variables_formula[2]


  if("." %in% variables_formula){
    vars<-setdiff(names(data),c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))
    variables_formula <- all.vars(formula)
  }



  variables_existent <- all(variables_formula %in% names(data))
  if (!variables_existent) stop("One or more variables from the formula do not exist in the data.")

  rm(variables_existent)

  if(any(sapply(data[,variables_formula],is.character)))stop("Some columns are of type character. Only numeric or factor variables are allowed.")


  all_terms <- attr(terms(formula), "term.labels")
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  if(length(strata_terms) >= 1) stop("The 'glmnet' package does not support the use of 'strata()' in the formula.")

  rm(all_terms,strata_terms)

  is_binary <- all(data[[failures]] %in% c(0, 1))

  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)

  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }

  .y <- Surv(data[[times]], data[[failures]])
  .x <- model.matrix(formula,data)[,-1]

  if(!(is.null(penalty))){

    if(length(penalty)!=length(variables_formula[-c(1,2)]))stop("Penalty length does not equal the number of variables.")
    if(!all(unique(penalty) %in% c(0,1)))stop("Penalty must be numeric and have only 0 or 1.")}



  # --- Function to create stratified folds ---
  create_stratified_folds <- function(status, K = 5, seed = 123){
    set.seed(seed)

    folds <- rep(NA, length(status))

    event_idx  <- which(status == 1)
    censor_idx <- which(status == 0)

    event_idx  <- sample(event_idx)
    censor_idx <- sample(censor_idx)

    folds[event_idx]  <- rep(1:K, length.out = length(event_idx))
    folds[censor_idx] <- rep(1:K, length.out = length(censor_idx))

    return(folds)
  }

  # --- Create stratified folds ---
  data$id <- 1:nrow(data)
  data$folds <- create_stratified_folds(status = data[[failures]], K = cv, seed = seed)

  # --- Check factor variables to ensure all levels appear in training ---
  if(any(sapply(data, is.factor))){
    factor_vars <- names(data)[sapply(data, is.factor)]

    inside <- function(factor, train, valid){
      all(unique(valid[, factor]) %in% unique(train[, factor]))
    }

    check_CVtune <- function(factors, CV){
      result <- unlist(lapply(factors, inside, train = CV$train, valid = CV$valid))
      all(result)
    }

    i <- 0
    success <- FALSE
    while(!success){
      i <- i + 1
      seed <- sample(1:1000, 1)
      set.seed(seed)

      data$folds <- create_stratified_folds(status = data[[failures]], K = cv, seed = seed)
      data$id <- 1:nrow(data)

      CVtune <- lapply(1:cv, function(k){
        train <- data[data$folds != k, ]
        valid <- data[data$folds == k, ]
        t_max_fold <- max(train[train[[failures]] == 1, times])
        list(train = train, valid = valid, t_max_fold = t_max_fold)
      })

      success <- check_CVtune(factor_vars, CVtune)

      if(!success){
        warning(paste("Seed was changed to", seed,
                      "because some factor levels were missing in training folds."))
      }

      if(i >= 3 & !success){
        stop("Some levels of factor variables in the validation set are missing from the training set.")
      }
    }
  }

  # --- Prepare foldid for cv.glmnet ---
  foldid <- data$folds

  if(!(is.null(penalty))) {
    .cv.lasso <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                           nfolds = cv, parallel = parallel, alpha=1,keep=F, foldid = foldid,
                           penalty.factor=penalty,
                           lambda=lambda
    )

  } else{
    .cv.lasso <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                           nfolds = cv, parallel = parallel, alpha=1,keep=F, foldid = foldid,
                           lambda=lambda
    )
  }


  return(list(optimal=list(lambda=.cv.lasso$lambda.min), results = data.frame(lambda=.cv.lasso$lambda, deviance=.cv.lasso$cvm)))
}


