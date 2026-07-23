predict.libsl <- function(object, newdata=NULL, newtimes=NULL, ...){


  if(!(is.null(newtimes))){
    newtimes<-sort(newtimes)
  }


  if(!(is.null(newdata))){
    variables_formula <- all.vars(object$formula)
    if(!all(variables_formula[-c(1,2)] %in% names(newdata)))stop("Error : some variables of the formula are missing.")
    if(any(sapply(newdata[,variables_formula[-c(1,2)]],is.character)))stop("Error : some columns are of type character. Only numeric or factor variables are allowed.")

  }



  # If the user does not provide newtimes or newdata, by default, they will get the values from the original model.
  if(is.null(newdata) & is.null(newtimes))  {
    .survivals<-object$predictions[,-1]
    .time<-object$times[-1]
    nrow<-nrow(object$data)
    ncol<-length(.time)
  }



  # Case where we don't have `newdata` but have `newtimes`.
  if(is.null(newdata) & !(is.null(newtimes))){
    nrow<-nrow(object$data)
    ncol<-length(newtimes)

    # Using `predict` for functions in the `flexsurv` package"
    if(object$library=="LIB_AFTgamma" | object$library=="LIB_AFTggamma"  | object$library=="LIB_AFTweibull"  |
       object$library=="LIB_AFTllogis" | object$library=="LIB_PHexponential"| object$library=="LIB_PHgompertz" |
       object$library=="LIB_PHspline" ){

      .flex<-object$model
      .data<-object$data
      .time<-newtimes
      .survivalist<-predict(
        .flex,
        newdata=.data,
        type = "survival",
        times = .time
      )

      if(length(newtimes)>1){
        .survivals <- t(sapply(.survivalist$.pred, function(x) x[[2]]))
      }
      else{
        .survivals<-.survivalist$.pred_survival
      }




    }

    # For the following functions, we have step survival functions, so we will use findInterval to determine where the new times fall.

    if(object$library=="LIB_COXen" | object$library=="LIB_COXlasso" | object$library=="LIB_COXridge"
       | object$library=="LIB_RSF" | object$library=="LIB_COXall" | object$library=="LIB_COXaic"){

      .surv<-object$predictions
      .times<-object$times
      .time<-newtimes
      .idx<-findInterval(.time,.times)

      # Extrapolation no allowed
      .idx[.time > max(.times)] <- NA
      .survivals<-.surv[,.idx]


    }

    if(object$library=="LIB_PLANN"){
      data <- model.matrix(object$formula, object$data)[,-1]
      .time<-newtimes
      if(0 %in% .time){
        .survivals <- as.matrix(predict(object$model, newdata = as.data.frame(data), newtimes = .time)$prediction)
      }
      else{
        .survivals <- as.matrix(predict(object$model, newdata = as.data.frame(data), newtimes = .time)$prediction)[,-1]
      }

    }

  }

  # Case where we have newdata but no new times.
  if(!(is.null(newdata)) & is.null(newtimes)){
    nrow<-nrow(newdata)
    ncol<-length(object$times[-1])


    if(object$library=="LIB_AFTgamma" | object$library=="LIB_AFTggamma"  | object$library=="LIB_AFTweibull"  |
       object$library=="LIB_AFTllogis" | object$library=="LIB_PHexponential"| object$library=="LIB_PHgompertz" |
       object$library=="LIB_PHspline" ){

      .flex<-object$model
      .formula<-object$formula
      .variables_formula <- all.vars(.formula)
      .times <- .variables_formula[1]
      .time<-object$times[-1]

      .survivalist<-predict(
        .flex,
        newdata=newdata,
        type = "survival",
        times = .time
      )

      .survivals <- t(sapply(.survivalist$.pred, function(x) x[[2]]))

    }

    if(object$library=="LIB_COXall" | object$library=="LIB_COXaic"){
      # The Cox model is a semi-parametric model, so we can calculate the new linear predictors,
      # but for the non-parametric part, we will use the estimation given by the original model.
      .coxphsurv<-survfit(object$model, newdata = newdata,se.fit = F)
      .sumcoxphsurv<-summary(.coxphsurv, times=object$times[-1])
      .survivals <- t(.sumcoxphsurv$surv)
      .time<-object$times[-1]
    }

    if(object$library=="LIB_COXen" | object$library=="LIB_COXlasso" | object$library=="LIB_COXridge"
    ){
      .x <- model.matrix(object$formula, object$data)[,-1]
      .lp <-predict(object$model, newx = .x)
      .times <- all.vars(object$formula)[1]
      .failures <- all.vars(object$formula)[2]
      if(!(.times %in% names(newdata))|!(.failures %in% names(newdata))){
        new<-data.frame(times=rep(0,nrow(newdata)),failures=rep(0,nrow(newdata)))
        colnames(new)<-c(.times,.failures)
        newdata<-cbind(new,newdata)
      }
      newdata<-newdata[,variables_formula]
      .x_bis <- model.matrix(object$formula, rbind(newdata, object$data[, names(newdata)]))[1:nrow(newdata), -1, drop=FALSE]
      .lp_bis <-predict(object$model, newx = .x_bis)
      .b <- glmnet_basesurv(object$data[[.times]], object$data[[.failures]], .lp, centered = FALSE)
      .H0 <- data.frame(value = .b$cumulative_base_hazard, time = .b$times)
      .survivals <- exp(matrix(exp(.lp_bis)) %*% t(as.matrix(-1*.H0$value)))
      .time<-object$times[-1]

    }

    if(object$library=="LIB_RSF"){
      .pred.rf <- predict(object$model,newdata=newdata)
      .survival <- cbind(rep(1, dim(.pred.rf$survival)[1]), .pred.rf$survival)
      .time.interest <- c(0, .pred.rf$time.interest)
      .time<-object$times[-1]
      .idx=findInterval(.time,.time.interest)
      .survivals<-.survival[,.idx]

    }

    if(object$library=="LIB_PLANN"){

      formula<-object$formula
      .times <- all.vars(formula)[1]
      .failures <- all.vars(formula)[2]
      if(!(.times %in% names(newdata))|!(.failures %in% names(newdata))){
        new<-data.frame(times=rep(0,nrow(newdata)),failures=rep(0,nrow(newdata)))
        colnames(new)<-c(.times,.failures)
        newdata<-cbind(new,newdata)
      }
      newdata<-newdata[,variables_formula]
      data <- model.matrix(formula,rbind(newdata,object$data[,names(newdata)]))[1:nrow(newdata),
                                                                                  ,
                                                                                  drop = FALSE]
      .time<-object$times[-1]
      .survivals <- as.matrix((predict(object$model, newdata = as.data.frame(data), newtimes=.time))$predictions)[,-1]

    }

  }
  if(!(is.null(newdata)) & !(is.null(newtimes))){

    nrow<-nrow(newdata)
    ncol<-length(newtimes)


    if(object$library=="LIB_AFTgamma" | object$library=="LIB_AFTggamma"  | object$library=="LIB_AFTweibull"  |
       object$library=="LIB_AFTllogis" | object$library=="LIB_PHexponential"| object$library=="LIB_PHgompertz" |
       object$library=="LIB_PHspline" ){

      .flex<-object$model
      .time<-newtimes

      .survivalist<-predict(
        .flex,
        newdata=newdata,
        type = "survival",
        times = .time
      )

      if(length(newtimes)>1){
        .survivals <- t(sapply(.survivalist$.pred, function(x) x[[2]]))
      }
      else{
        .survivals<-.survivalist$.pred_survival
      }



    }

    if(object$library=="LIB_COXall" | object$library=="LIB_COXaic"){

      .coxphsurv<-survfit(object$model, newdata = newdata,se.fit = F)
      .sumcoxphsurv<-summary(.coxphsurv, times=object$times)
      .surv <- t(.sumcoxphsurv$surv)
      .time<-newtimes
      .idx=findInterval(.time,object$times)
      .idx[.time > max(object$times)] <- NA
      if(nrow(.surv)==1){
        .survivals<-.surv[.idx]
      }
      else{
        .survivals<-.surv[,.idx]

      }

    }

    if(object$library=="LIB_COXen" | object$library=="LIB_COXlasso" | object$library=="LIB_COXridge"
    ){
      .x <- model.matrix(object$formula, object$data)[,-1]
      .lp <-predict(object$model, newx = .x)
      .times <- all.vars(object$formula)[1]
      .failures <- all.vars(object$formula)[2]
      if(!(.times %in% names(newdata))|!(.failures %in% names(newdata))){
        new<-data.frame(times=rep(0,nrow(newdata)),failures=rep(0,nrow(newdata)))
        colnames(new)<-c(.times,.failures)
        newdata<-cbind(new,newdata)
      }
      newdata<-newdata[,variables_formula]
      .x_bis <- model.matrix(object$formula, rbind(newdata, object$data[, names(newdata)]))[1:nrow(newdata), -1, drop=FALSE]
      .lp_bis <-predict(object$model, newx = .x_bis)
      .b <- glmnet_basesurv(object$data[[.times]], object$data[[.failures]], .lp, centered = FALSE)
      .H0 <- data.frame(value = c(0,.b$cumulative_base_hazard), time = c(0,.b$times))
      .surv <- exp(matrix(exp(.lp_bis)) %*% t(as.matrix(-1*.H0$value)))
      .time<-newtimes
      .idx<-findInterval(.time,.H0$time)
      .idx[.time > max(.H0$time)] <- NA
      if(nrow(.surv)==1){
        .survivals<-.surv[.idx]
      }else{
        .survivals<-.surv[,.idx]
      }



    }

    if(object$library=="LIB_RSF"){
      .pred.rf <- predict(object$model,newdata=newdata)
      .surv <- cbind(rep(1, dim(.pred.rf$survival)[1]), .pred.rf$survival)
      .time.interest <- c(0, .pred.rf$time.interest)
      .time<-newtimes
      .idx=findInterval(.time,.time.interest)
      .idx[.time > max(.time.interest)] <- NA
      .survivals<-.surv[,.idx]

    }

    if(object$library=="LIB_PLANN"){

      formula<-object$formula
      .times <- all.vars(formula)[1]
      .failures <- all.vars(formula)[2]
      if(!(.times %in% names(newdata))|!(.failures %in% names(newdata))){
        new<-data.frame(times=rep(0,nrow(newdata)),failures=rep(0,nrow(newdata)))
        colnames(new)<-c(.times,.failures)
        newdata<-cbind(new,newdata)
      }
      newdata<-newdata[,variables_formula]
      data <- model.matrix(formula,rbind(newdata,object$data[,names(newdata)]))[1:nrow(newdata),
                                                                                  ,
                                                                                  drop = FALSE]
      .time<-newtimes
      if(0 %in% .time){
        .survivals <- as.matrix(predict(object$model, newdata = as.data.frame(data), newtimes = .time)$prediction)
      }
      else{
        .survivals <- as.matrix(predict(object$model, newdata = as.data.frame(data), newtimes = .time)$prediction)[,-1]
      }
    }


  }


  .predictions=matrix(.survivals,ncol=ncol,nrow = nrow)

  if(object$library=="LIB_PLANN"){
    .predictions=.survivals
  }




  return(list(times=.time, predictions=.predictions))
}
