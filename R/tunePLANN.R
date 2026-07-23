tunePLANN <- function(formula, data, cv=10, inter=1, size=c(2, 4, 6, 8, 10), decay=c(0.001, 0.01, 0.02, 0.05), maxit=100, MaxNWts=10000, maxtime=NULL, seed=NULL, metric="auc", pro.time=NULL, ROC.precision=seq(.01, .99, by=.01)){


  if(is.null(seed)){
    seed<-sample(1:1000,1)
  }

  if (missing(formula)) stop("The 'formula' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(inter)) stop("The 'inter' argument is required.")
  if (missing(decay)) stop("The 'decay' argument is required.")
  if (missing(maxit)) stop("The 'maxit' argument is required.")
  if (missing(MaxNWts)) stop("The 'MaxNWts' argument is required.")


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
  ns_terms <- grep("ns\\(", all_terms, value = TRUE)
  bs_terms <- grep("bs\\(", all_terms, value = TRUE)

  if(length(strata_terms) >= 1) stop("The 'survivalPLANN' package does not support the use of 'strata()' in the formula.")


  if((length(ns_terms) >= 1)|(length(bs_terms) >= 1)|(length(all_terms)>(length(variables_formula)-2))){
    vars<-setdiff(variables_formula,c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))
  }

  rm(all_terms,strata_terms,ns_terms,bs_terms)



  is_binary <- all(data[[failures]] %in% c(0, 1))


  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }



  if(metric=="ll"){
    haz_function <- function(surv, times) {
      x <- 1
      result <- sapply(2:length(surv), function(i) {
        value <- surv[i]
        if(value != surv[i-1]) {
          x <<- x + 1
        }
        return(x)
      })

      df <- data.frame(temps = times, value = -log(surv), result = c(1, result))
      df_unique <- df[!duplicated(df$result), ]

      if(nrow(df_unique) > 1) {
        # Calculation of hazard rates for each interval
        diff_1 <- diff(df_unique$temps)
        diff_2 <- diff(df_unique$value)
        taux_intervalles <- diff_2 / diff_1

        # Creating a result vector of the correct length
        bj <- rep(NA, length(times))

        # For each time EXCEPT THE LAST, assign the corresponding rate
        for(i in 1:(length(times)-1)) {
          idx_interval <- findInterval(times[i], df_unique$temps)

          # Ensure that the index is valid
          if(idx_interval > 0 && idx_interval <= length(taux_intervalles)) {
            bj[i] <- taux_intervalles[idx_interval]
          }
        }

        # The last time remains NA (already done by the initialization)

      } else {
        # No changes: NA everywhere
        bj <-  rep(NA, length(times))
      }

      return(bj)
    }

  }


  if(is.null(pro.time)){
    pro.time=median(data[[times]])
  }


  .data_bis<-data
  .time <-unique(sort(c(0,pro.time,data[[times]])))





  create_stratified_folds <- function(status, K = 5, seed=123){

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

  set.seed(seed)
  data$id = 1:nrow(data)

  data$folds <- create_stratified_folds(data[[failures]], K = cv, seed=seed)

  CVtune <- lapply(1:cv, function(i) {

    # create train and valid
    train <- data[data$folds != i, ]
    valid <- data[data$folds == i, ]

    # calculate t_max_fold
    t_max_fold <- max(train[train[[failures]] == 1, times])

    list(
      train = train,
      valid = valid,
      t_max_fold = t_max_fold
    )
  })



  if(any(sapply(data, is.factor))){

    factor_vars <- names(data)[sapply(data, is.factor)]

    # Function to check that all factor levels in validation exist in training
    inside <- function(factor, train, valid){
      all(unique(valid[,factor]) %in% unique(train[,factor]))
    }

    # Function to check all factor variables for one CV split
    check_CVtune <- function(factors, CV){
      result <- unlist(lapply(factors, inside, train = CV$train, valid = CV$valid))
      all(result)
    }

    i <- 0
    success <- FALSE

    while(!success){
      i <- i + 1

      # Generate a random seed for reproducibility
      seed <- sample(1:1000, 1)

      # Create stratified folds based on event status
      data$folds <- create_stratified_folds(status = data[[failures]], K = cv, seed = seed)
      data$id <- 1:nrow(data)

      # Build CVtune list
      CVtune <- lapply(1:cv, function(k){
        train <- data[data$folds != k, ]
        valid <- data[data$folds == k, ]
        t_max_fold <- max(train[train[[failures]] == 1, times])
        list(
          train = train,
          valid = valid,
          t_max_fold = t_max_fold
        )
      })

      # Check that all factor levels are present in training
      success <- check_CVtune(factor_vars, CVtune)

      if(!success){
        warning(paste("The seed has been changed to", seed,
                      "because some factor levels are missing in training folds."))
      }

      # Stop if after 3 attempts it still fails
      if(i >= 3 & !success){
        stop("Certain levels of some factor variables in the validation sample are not present in the training sample. Please check your dataset.")
      }
    }

  }


  t_max_global <- min(unlist(lapply(CVtune,function(x)(return(x$t_max_fold)))))
  .time <- .time[.time <= t_max_global]


  if (is.null(maxtime) || maxtime < max(.time)) {
    maxtime <- max(.time) + 1
  }







  hp<-expand.grid(inter, decay, size, maxit,MaxNWts)
  colnames(hp)<-c("inter", "decay", "size", "maxit","MaxNWts")
  hp_list<-apply(hp, 1, as.list)

  plann_function<-function(x,y){

    model_matrix <- model.matrix(formula, data = rbind(x$train,x$valid))[1:nrow(x$train),-1]

    data_bis<-data.frame(times=x$train[[times]],failures=x$train[[failures]])
    colnames(data_bis)<-c(times,failures)

    vars<-setdiff(names(cbind(data_bis,model_matrix)),c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    .formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))

    .plann <- sPLANN(formula=.formula, data=cbind(data_bis,model_matrix), inter=y$inter, size = y$size, decay = y$decay,  maxit = y$maxit, MaxNWts = y$MaxNWts, pro.time=maxtime)

    newdata <- model.matrix(.formula, rbind(x$valid, x$train))[
      1:nrow(x$valid),
      ,
      drop = FALSE
    ]
    .survivals <- predict(.plann, newdata = data.frame(newdata), newtimes = .time)$predictions
    return(predictions=list(survivals=.survivals, id=x$valid$id))

  }



  y_function<-function(y){
    result<-lapply(CVtune,plann_function,y=y)
    surv_list<-lapply(result,function(x)(return(x$survivals)))
    id_list<- lapply(result,function(x)(return(x$id)))
    survivals<-do.call(rbind,surv_list)
    id<-unlist(id_list)
    return(predictions=list(survivals=survivals, id=id))

  }



  result<-lapply(hp_list,y_function)

  metric_function<-function(x){
    data<-data[x$id, ]
    survivals.matrix<-x$survivals
    hazards.matrix<-NULL
    if(metric=="ll"){
      hazards.matrix<-t(apply(x$survivals,1,haz_function,times=.time))
    }
    resultat<-metrics(metric=metric,formula=formula,data=data,survivals.matrix=survivals.matrix,hazards.matrix=hazards.matrix,prediction.times=.time,pro.time=pro.time,ROC.precision=ROC.precision)
    return(resultat)

  }


  metric_results<-unlist(lapply(result,metric_function))


  if(metric %in% c("bs","ibs","ribs","bll","ibll","ribll")){
    .idx<-which.min(metric_results)
  }else{
    .idx<-which.max(metric_results)
  }

  data<-.data_bis


  return( list(optimal=list(
    inter=hp[.idx,1],
    decay=hp[.idx,2],
    size=hp[.idx,3],
    maxit=hp[.idx,4],
    MaxNWts=hp[.idx,5]

  ), results=cbind(hp,metric_results) ))

}





