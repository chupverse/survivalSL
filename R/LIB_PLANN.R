LIB_PLANN <- function(formula,
                      data, inter, size, decay, maxit, MaxNWts, maxtime=NULL){

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


  is_binary <- all(data[[failures]] %in% c(0, 1))


  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }



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


  if(is.null(maxtime)||maxtime<max(data[[times]])){
    maxtime<-max(data[[times]])+1
  }


  model_matrix <- model.matrix(formula, data = data)[,-1]

  data_bis<-data.frame(times=data[[times]],failures=data[[failures]])
  colnames(data_bis)<-c(times,failures)

  vars<-setdiff(names(cbind(data_bis,model_matrix)),c(times,failures))
  .outcome <- paste("Surv(", times, ",", failures, ")")
  .formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))

  .plann <- sPLANN(formula=.formula, data=cbind(data_bis,model_matrix), inter=inter, size = size, decay = decay,  maxit = maxit, MaxNWts = MaxNWts, pro.time=maxtime)

  .time <- c(0,sort(unique(data[[times]])))

  .survivals <- (predict(.plann, newdata=data.frame(model_matrix), newtimes=.time))$predictions


  .obj <- list(model=.plann,
               library="LIB_PLANN",
               formula=formula,
               data=data,
               times=.time,  predictions=.survivals)

  class(.obj) <- "libsl"

  return(.obj)
}



