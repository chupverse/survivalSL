LIB_COXridge <- function(formula,
                         data, penalty=NULL, lambda){

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

  all_terms <- attr(terms(formula), "term.labels")
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  if(length(strata_terms) >= 1) stop("The 'glmnet' package does not support the use of 'strata()' in the formula.")

  rm(all_terms,strata_terms)


  if(any(sapply(data[,variables_formula],is.character)))stop("Some columns are of type character. Only numeric or factor variables are allowed.")

  is_binary <- all(data[[failures]] %in% c(0, 1))


  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  .x <- model.matrix(formula,data)[,-1]
  .y <- Surv(data[[times]], data[[failures]])

  if(!(is.null(penalty))){

    if(length(penalty)!=length(variables_formula[-c(1,2)]))stop("Penalty length does not equal the number of variables.")
    if(!all(unique(penalty) %in% c(0,1)))stop("Penalty must be numeric and have only 0 or 1.")}


  if(!(is.null(penalty))) {
    #.penalty.factor <- rep(1,length(colnames(.x)))
    #.penalty.factor[which(colnames(.x) %in% var)] <- 0
    .ridge<- glmnet(x = .x, y = .y, lambda = lambda,
                    type.measure = "deviance", family = "cox",cox.ties = "breslow",
                    alpha = 0, penalty.factor = penalty)

  }

  else{
    .ridge <- glmnet(x = .x, y = .y, lambda = lambda,  type.measure = "deviance",
                     family = "cox",cox.ties = "breslow", alpha = 0)
  }


  .lp.ridge <- predict(.ridge, newx = .x)
  .b <- glmnet_basesurv(data[[times]], data[[failures]], .lp.ridge, centered = FALSE)
  .H0 <- data.frame(value = .b$cumulative_base_hazard, time = .b$times)

  .pred <- exp(matrix(exp(.lp.ridge)) %*% t(as.matrix(-1*.H0$value)))

  .survivals<-cbind(rep(1, dim(.pred)[1]), .pred)

  .obj <- list(model=.ridge,
               library="LIB_COXridge",
               formula=formula,
               data=data,
               times=c(0,.H0$time),predictions=.survivals)

  class(.obj) <- "libsl"

  return(.obj)
}


