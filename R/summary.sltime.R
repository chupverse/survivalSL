summary.sltime <- function(object, newdata=NULL, method="sl",
                           ROC.precision=seq(.01,.99,.01), digits=7,pro.time=NULL, ...) {


  formula<-object$formula
  variables_formula <- all.vars(formula)
  times <- variables_formula[1]
  failures <- variables_formula[2]

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



  if(is.null(pro.time)) {
    pro.time <- median(object$times)
  }

  time.pred <- unique(sort(c(0,pro.time,object$times)))



  if(is.null(newdata)){

    survivals.matrix <- predict(object, newdata=object$data, newtimes=time.pred)$predictions[[method]]

    if(!(method %in% c("sl","LIB_PLANN","LIB_SNN","LIB_RSF","LIB_COXall","LIB_COXaic","LIB_COXlasso","LIB_COXen","LIB_COXridge"))){
      .flex<-object$models[[method]]$model
      .data<-object$models[[method]]$data
      .time<-time.pred
      .hazlist<-predict(
        .flex,
        newdata=.data,
        type = "haz",
        times = .time
      )
      hazards.matrix <- t(sapply(.hazlist$.pred, function(x) x[[2]]))
    }else{
      hazards.matrix<-t(apply(survivals.matrix,1,haz_function,times=time.pred))
      }

    
      out<-list(metrics=round(  data.frame(
        p_ci = metrics(metric="p_ci", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        uno_ci = metrics(metric="uno_ci", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        auc = metrics(metric="auc", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bs = metrics(metric="bs", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                     hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibs = metrics(metric="ibs", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribs = metrics(metric="ribs", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bll = metrics(metric="bll", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibll = metrics(metric="ibll", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribll = metrics(metric="ribll", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                        hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ll= metrics(metric="ll", formula=object$formula, data=object$data, survivals.matrix=survivals.matrix,
                    hazards.matrix=hazards.matrix,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision)), digits = digits ),
        method=method,
        pro.time=pro.time,
        ROC.precision=ROC.precision) 




  }else{

    survivals.matrix <- predict(object, newdata=newdata, newtimes=time.pred)$predictions[[method]]

    if(!(method %in% c("sl","LIB_PLANN","LIB_SNN","LIB_RSF","LIB_COXall","LIB_COXaic","LIB_COXlasso","LIB_COXen","LIB_COXridge"))){
      .flex<-object$models[[method]]$model
      .time<-time.pred
      .hazlist<-predict(
        .flex,
        newdata=newdata,
        type = "haz",
        times = .time
      )
      hazards.matrix <- t(sapply(.hazlist$.pred, function(x) x[[2]]))
      }else{
      hazards.matrix<-t(apply(survivals.matrix,1,haz_function,times=time.pred))

    }


      out<-list(metrics=round(  data.frame(
        p_ci = metrics(metric="p_ci",formula=object$formula , data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        uno_ci = metrics(metric="uno_ci", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        auc = metrics(metric="auc", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bs = metrics(metric="bs", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                     hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibs = metrics(metric="ibs", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribs = metrics(metric="ribs", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        bll = metrics(metric="bll", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                      hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ibll = metrics(metric="ibll", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                       hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ribll = metrics(metric="ribll", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                        hazards.matrix=NULL,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision),
        ll= metrics(metric="ll", formula=object$formula, data=newdata, survivals.matrix=survivals.matrix,
                    hazards.matrix=hazards.matrix,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision)), digits = digits ),
        method=method,
        pro.time=pro.time,
        ROC.precision=ROC.precision) 




  }
  
  class(out)<-"summary.sltime"
  out



}





