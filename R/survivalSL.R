survivalSL <- function(formula, data, methods, metric="auc", penalty=NULL,
                       cv=10, param.tune=NULL, pro.time=NULL,  optim.local.min=FALSE,
                       ROC.precision=seq(.01,.99,.01), param.weights.fix=NULL,
                       param.weights.init=NULL,
                       seed=NULL, optim.method="Nelder-Mead",maxit=1000,show_progress=TRUE) {



  if (missing(methods)) stop("The 'method' argument is required.")
  if (missing(data)) stop("The 'data' argument is required.")
  if (missing(formula)) stop("The 'formula' argument is required.")





  if(show_progress){
    M <- length(methods)
    total <- M + (M * cv) + M + 3
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    progress <- 0

    update_progress <- function(){
      progress <<- progress + 1
      setTxtProgressBar(pb, progress)
      if(progress == total) close(pb)
    }
  }



  if(!(optim.method %in% c("SANN","Nelder-Mead")))stop("You can either choose the SANN or Nelder-Mead methods.")


  if(is.null(seed)){
    seed<-sample(1:1000,1)
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


  variables_formula <- all.vars(formula)

  times <- variables_formula[1]
  failures <- variables_formula[2]


  if("." %in% variables_formula){
    vars<-setdiff(names(data),c(times,failures))
    .outcome <- paste("Surv(", times, ",", failures, ")")
    formula <- as.formula(paste(.outcome, "~", paste(vars, collapse = " + ")))
    variables_formula <- all.vars(formula)
  }


  is_binary <- all(data[[failures]] %in% c(0, 1))


  if (! is_binary) stop("The 'failures' variable is not coded as 0/1.")

  rm(is_binary)

  # Check if all the variables exist in the dataset.
  variables_existent <- all(variables_formula %in% names(data))
  if (!variables_existent) stop("One or more variables from the formula do not exist in the data.")

  rm(variables_existent)

  if(any(sapply(data[,variables_formula],is.character)))stop("Some columns are of type character. Only numeric or factor variables are allowed.")


  if (any(is.na(data[,variables_formula]))){
    subset_data<-na.omit(data[,variables_formula])
    data<-cbind(subset_data, data[!colnames(data) %in% colnames(subset_data), drop = FALSE])
    warning("Data need to be without NA. NA is removed")
  }


  if(!(is.null(penalty))){

    if(length(penalty)!=length(variables_formula[-c(1,2)]))stop("Penalty length does not equal the number of variables.")
    if(!all(unique(penalty) %in% c(0,1)))stop("Penalty must be numeric and have only 0 or 1.")}



  #####################
  ### Quality tests ###
  #####################

  if(length(methods)<=1)
  { stop("Number of methods need to be greater or equal to 2 to estimate a SuperLearner")   }

  if(length(metric)>1){
    warning(paste0("SuperLearner is currently developped for one metric. Results for metric ",metric[1]))
    metric=metric[1]
  }

  if(min(metric%in%c("uno_ci","p_ci","bs","ibs","ibll","bll", "ribs","ribll","auc","ll"))==0){
    stop("The argument \"metric\" must be Brier score (bs),
         Pencina concordance index (p_ci),
         Uno concordance index (uno_ci),
         integrated Brier score (ibs), the binomial log-likelihood (bll),
         the integrated binomial log-likelihood (ibll), the restricted ibs (ribs),
         the restricted ibll (ribll),
         the area under the ROC curve (auc), or the log-likelihood (ll)")
  }

  if(!is.data.frame(data)){
    stop("The argument \"data\" need to be a data.frame") }




  if (min(data[[times]])<=0){
    stop("Time variable need to be positive")
  }

  if (cv < 3 | !is.numeric(cv)) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }

  .meth_rm=c()
  if(sum(methods %in% "LIB_AFTgamma")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_AFTgamma")[-1])
    warning("SuperLearner can use only one LIB_AFTgamma method. We remove the others.")
  }
  if(sum(methods %in% "LIB_AFTllogis")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_AFTllogis")[-1])
    warning("SuperLearner can use only one LIB_AFTllogis method. We remove the others.")
  }
  if(sum(methods %in% "LIB_AFTggamma")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_AFTggamma")[-1])
    warning("SuperLearner can use only one LIB_AFTggamma method. We remove the others.")
  }
  if(sum(methods %in% "LIB_AFTweibull")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_AFTweibull")[-1])
    warning("SuperLearner can use only one LIB_AFTweibull method. We remove the others.")
  }
  if(sum(methods %in% "LIB_PHexponential")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_PHexponential")[-1])
    warning("SuperLearner can use only one LIB_PHexponential method. We remove the others.")
  }
  if(sum(methods %in% "LIB_PHgompertz")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_PHgompertz")[-1])
    warning("SuperLearner can use only one LIB_PHgompertz method. We remove the others.")
  }
  if(sum(methods %in% "LIB_PHspline")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_PHspline")[-1])
    warning("SuperLearner can use only one LIB_PHspline method. We remove the others.")
  }

  if(sum(methods %in% "LIB_COXaic")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_COXaic")[-1])
    warning("SuperLearner can use only one LIB_COXaic method. We remove the others.")
  }


  if(sum(methods %in% "LIB_COXall")>=2){
    .meth_rm=c(.meth_rm,which(methods=="LIB_COXall")[-1])
    warning("SuperLearner can use only one LIB_COXall method. We remove the others.")
  }

  if(sum(methods %in% "LIB_COXlasso")==1){
    if(!(is.null(param.tune[[which(methods=="LIB_COXlasso")]]))){
      if(!is.list(param.tune[[which(methods=="LIB_COXlasso")]])){
        stop("Argument param.tune for LIB_COXlasso need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="LIB_COXlasso")]])%in%"lambda"))==0){
        stop("Tune parameters for LIB_COXlasso need to have lambda")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_COXlasso")]]$lambda)|
           is.null(param.tune[[which(methods=="LIB_COXlasso")]]$lambda))){
        stop("lambda for LIB_COXlasso need to be a scalar or a vector or NULL")
      }
    }
  }

  if(sum(methods %in% "LIB_COXlasso")>=2){
    if(length(param.tune[which(methods=="LIB_COXlasso")])!=length(unique(param.tune[which(methods=="LIB_COXlasso")]))){
      stop("Tune parameters for LIB_COXlasso methods need to be unique")
    }
    for (i in 1:sum(methods %in% "LIB_COXlasso")){
      if(!(is.null(param.tune[[which(methods=="LIB_COXlasso")[i]]]))){
        if(!is.list(param.tune[[which(methods=="LIB_COXlasso")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th LIB_COXlasso need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="LIB_COXlasso")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th LIB_COXlasso need to have lambda"))
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_COXlasso")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="LIB_COXlasso")[i]]]$lambda))){
          stop(paste("lambda for the ",i,"th LIB_COXlasso need to be a scalar or a vector or NULL"))
        }
      }
    }
  }



  if(sum(methods %in% "LIB_PHspline")==1){
    if(!(is.null(param.tune[[which(methods=="LIB_PHspline")]]))){
      if(!is.list(param.tune[[which(methods=="LIB_PHspline")]])){
        stop("Argument param.tune for LIB_PHspline need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="LIB_PHspline")]])%in%"k"))==0){
        stop("Tune parameters for LIB_PHspline need to have k")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_PHspline")]]$k)|
           is.null(param.tune[[which(methods=="LIB_PHspline")]]$k))){
        stop("lambda for LIB_PHspline need to be a scalar or a vector or NULL")
      }
    }
  }

  if(sum(methods %in% "LIB_PHspline")>=2){
    if(length(param.tune[which(methods=="LIB_PHspline")])!=length(unique(param.tune[which(methods=="LIB_PHspline")]))){
      stop("Tune parameters for LIB_PHspline methods need to be unique")
    }
    for (i in 1:sum(methods %in% "LIB_PHspline")){
      if(!(is.null(param.tune[[which(methods=="LIB_PHspline")[i]]]))){
        if(!is.list(param.tune[[which(methods=="LIB_PHspline")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th LIB_PHspline need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="LIB_PHspline")[i]]])%in%"k"))==0){
          stop(paste("Tune parameters for the ",i,"th LIB_PHspline need to have k"))
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_PHspline")[i]]]$k)|
             is.null(param.tune[[which(methods=="LIB_PHspline")[i]]]$k))){
          stop(paste("lambda for the ",i,"th LIB_PHspline need to be a scalar or a vector or NULL"))
        }
      }
    }
  }


  if(sum(methods %in% "LIB_COXridge")==1){
    if(!(is.null(param.tune[[which(methods=="LIB_COXridge")]]))){
      if(!is.list(param.tune[[which(methods=="LIB_COXridge")]])){
        stop("Argument param.tune for LIB_COXridge need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="LIB_COXridge")]])%in%"lambda"))==0){
        stop("Tune parameters for LIB_COXridge need to have lambda")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_COXridge")]]$lambda)|
           is.null(param.tune[[which(methods=="LIB_COXridge")]]$lambda))){
        stop("lambda for LIB_COXridge need to be a scalar or a vector or NULL")
      }
    }
  }
  if(sum(methods %in% "LIB_COXridge")>=2){
    if(length(param.tune[which(methods=="LIB_COXridge")])!=length(unique(param.tune[which(methods=="LIB_COXridge")]))){
      stop("Tune parameters for LIB_COXridge methods need to be unique")
    }
    for (i in 1:sum(methods %in% "LIB_COXridge")){
      if(!(is.null(param.tune[[which(methods=="LIB_COXridge")[i]]]))){
        if(!is.list(param.tune[[which(methods=="LIB_COXridge")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th LIB_COXridge need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="LIB_COXridge")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th LIB_COXridge need to have lambda"))
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_COXridge")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="LIB_COXridge")[i]]]$lambda))){
          stop(paste("lambda for the ",i,"th LIB_COXridge need to be a scalar or a vector or NULL"))
        }
      }
    }
  }

  if(sum(methods %in% "LIB_COXen")==1){
    if(!(is.null(param.tune[[which(methods=="LIB_COXen")]]))){
      if(!is.list(param.tune[[which(methods=="LIB_COXen")]])){
        stop("Argument param.tune for LIB_COXen need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="LIB_COXen")]])%in%"lambda"))==0){
        stop("Tune parameters for LIB_COXen need to have lambda")
      }
      if(sum((names(param.tune[[which(methods=="LIB_COXen")]])%in%"alpha"))==0){
        stop("Tune parameters for LIB_COXen need to have alpha")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_COXen")]]$lambda)|
           is.null(param.tune[[which(methods=="LIB_COXen")]]$lambda))){
        stop("lambda for LIB_COXen need to be a scalar or a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_COXen")]]$alpha)|
           is.null(param.tune[[which(methods=="LIB_COXen")]]$alpha))){
        stop("alpha for LIB_COXen need to be a scalar or a vector or NULL")
      }
      if(min(param.tune[[which(methods=="LIB_COXen")]]$alpha)<0 | max(param.tune[[which(methods=="LIB_COXen")]]$alpha)>1){
        stop("tune parameters for LIB_COXen alpha need to be in ]0;1[")
      }
    }
  }
  if(sum(methods %in% "LIB_COXen")>=2){
    if(length(param.tune[which(methods=="LIB_COXen")])!=length(unique(param.tune[which(methods=="LIB_COXen")]))){
      stop("Tune parameters for LIB_COXen methods need to be unique")
    }
    for (i in 1:sum(methods %in% "LIB_COXen")){
      if(!(is.null(param.tune[[which(methods=="LIB_COXen")[i]]]))){
        if(!is.list(param.tune[[which(methods=="LIB_COXen")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th LIB_COXen need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="LIB_COXen")[i]]])%in%"lambda"))==0){
          stop(paste("Tune parameters for the ",i,"th LIB_COXen need to have lambda"))
        }
        if(sum((names(param.tune[[which(methods=="LIB_COXen")[i]]])%in%"alpha"))==0){
          stop(paste("Tune parameters for the ",i,"th LIB_COXen need to have alpha"))
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_COXen")[i]]]$lambda)|
             is.null(param.tune[[which(methods=="LIB_COXen")[i]]]$lambda))){
          stop(paste("lambda for the ",i,"th LIB_COXen need to be a scalar or a vector or NULL"))
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_COXen")[i]]]$alpha)|
             is.null(param.tune[[which(methods=="LIB_COXen")[i]]]$alpha))){
          stop(paste("alpha for the ",i,"th LIB_COXen need to be a scalar or a vector or NULL"))
        }
        if(min(param.tune[[which(methods=="LIB_COXen")[i]]]$alpha)<0 | max(param.tune[[which(methods=="LIB_COXen")[i]]]$alpha)>1){
          stop("tune parameters for LIB_COXen alpha need to be in ]0;1[")
        }
      }
    }
  }



  if(sum(methods %in% "LIB_RSF")==1){
    if(!(is.null(param.tune[[which(methods=="LIB_RSF")]]))){
      if(!is.list(param.tune[[which(methods=="LIB_RSF")]])){
        stop("Argument param.tune for LIB_RSF need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="LIB_RSF")]])%in%"nodesize"))==0){
        stop("Tune parameters for LIB_RSF need to have nodesize")
      }
      if(sum((names(param.tune[[which(methods=="LIB_RSF")]])%in%"mtry"))==0){
        stop("Tune parameters for LIB_RSF need to have mtry")
      }
      if(sum((names(param.tune[[which(methods=="LIB_RSF")]])%in%"ntree"))==0){
        stop("Tune parameters for LIB_RSF need to have ntree")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_RSF")]]$nodesize)|
           is.null(param.tune[[which(methods=="LIB_RSF")]]$nodesize))){
        stop("nodesize for LIB_RSF need to be a scalar or a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_RSF")]]$mtry)|
           is.null(param.tune[[which(methods=="LIB_RSF")]]$mtry))){
        stop("mtry for LIB_RSF need to be a scalar or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_RSF")]]$ntree)|
           is.null(param.tune[[which(methods=="LIB_RSF")]]$ntree))){
        stop("ntree for LIB_RSF need to be a scalar or NULL")
      }
    }
  }
  if(sum(methods %in% "LIB_RSF")>=2){
    if(length(param.tune[which(methods=="LIB_RSF")])!=length(unique(param.tune[which(methods=="LIB_RSF")]))){
      stop("Tune parameters for LIB_RSF methods need to be unique")
    }
    for (i in 1:sum(methods %in% "LIB_RSF")){
      if(!is.list(param.tune[[which(methods=="LIB_RSF")[i]]])){
        stop(paste("Argument param.tune for the ",i,"th LIB_RSF need to be a list"))
      }
      if(!(is.null(param.tune[[which(methods=="LIB_RSF")[i]]]))){
        if(sum((names(param.tune[[which(methods=="LIB_RSF")[i]]])%in%"nodesize"))==0){
          stop("Tune parameters for LIB_RSF need to have nodesize")
        }
        if(sum((names(param.tune[[which(methods=="LIB_RSF")[i]]])%in%"mtry"))==0){
          stop("Tune parameters for LIB_RSF need to have mtry")
        }
        if(sum((names(param.tune[[which(methods=="LIB_RSF")[i]]])%in%"ntree"))==0){
          stop("Tune parameters for LIB_RSF need to have ntree")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_RSF")[i]]]$nodesize)|
             is.null(param.tune[[which(methods=="LIB_RSF")[i]]]$nodesize))){
          stop("nodesize for LIB_RSF need to be a scalar or a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_RSF")[i]]]$mtry)|
             is.null(param.tune[[which(methods=="LIB_RSF")[i]]]$mtry))){
          stop("mtry for LIB_RSF need to be a scalar or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_RSF")[i]]]$ntree)|
             is.null(param.tune[[which(methods=="LIB_RSF")[i]]]$ntree))){
          stop("ntree for LIB_RSF need to be a scalar or NULL")
        }
      }
    }
  }





  if(sum(methods %in% "LIB_PLANN")==1){
    if(!(is.null(param.tune[[which(methods=="LIB_PLANN")]]))){
      if(!is.list(param.tune[[which(methods=="LIB_PLANN")]])){
        stop("Argument param.tune for LIB_PLANN need to be a list")
      }
      if(sum((names(param.tune[[which(methods=="LIB_PLANN")]])%in%"inter"))==0){
        stop("Tune parameters for LIB_PLANN need to have inter")
      }
      if(sum((names(param.tune[[which(methods=="LIB_PLANN")]])%in%"size"))==0){
        stop("Tune parameters for LIB_PLANN need to have size")
      }
      if(sum((names(param.tune[[which(methods=="LIB_PLANN")]])%in%"decay"))==0){
        stop("Tune parameters for LIB_PLANN need to have decay")
      }
      if(sum((names(param.tune[[which(methods=="LIB_PLANN")]])%in%"maxit"))==0){
        stop("Tune parameters for LIB_PLANN need to have maxit")
      }
      if(sum((names(param.tune[[which(methods=="LIB_PLANN")]])%in%"MaxNWts"))==0){
        stop("Tune parameters for LIB_PLANN need to have MaxNWts")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")]]$inter)|
           is.null(param.tune[[which(methods=="LIB_PLANN")]]$inter))){
        stop("inter for LIB_PLANN need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")]]$size)|
           is.null(param.tune[[which(methods=="LIB_PLANN")]]$size))){
        stop("size for LIB_PLANN need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")]]$decay)|
           is.null(param.tune[[which(methods=="LIB_PLANN")]]$decay))){
        stop("decay for LIB_PLANN need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")]]$maxit)|
           is.null(param.tune[[which(methods=="LIB_PLANN")]]$maxit))){
        stop("maxit for LIB_PLANN need to be a scalar, a vector or NULL")
      }
      if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")]]$MaxNWts)|
           is.null(param.tune[[which(methods=="LIB_PLANN")]]$MaxNWts))){
        stop("MaxNWts for LIB_PLANN need to be a scalar, a vector or NULL")
      }
    }
  }
  if(sum(methods %in% "LIB_PLANN")>=2){
    if(length(param.tune[which(methods=="LIB_PLANN")])!=length(unique(param.tune[which(methods=="LIB_PLANN")]))){
      stop("Tune parameters for LIB_PLANN methods need to be unique")
    }
    for (i in 1:sum(methods %in% "LIB_PLANN")){
      if(!(is.null(param.tune[[which(methods=="LIB_PLANN")[i]]]))){
        if(!is.list(param.tune[[which(methods=="LIB_PLANN")[i]]])){
          stop(paste("Argument param.tune for the ",i,"th LIB_PLANN need to be a list"))
        }
        if(sum((names(param.tune[[which(methods=="LIB_PLANN")[i]]])%in%"inter"))==0){
          stop("Tune parameters for LIB_PLANN need to have inter")
        }
        if(sum((names(param.tune[[which(methods=="LIB_PLANN")[i]]])%in%"size"))==0){
          stop("Tune parameters for LIB_PLANN need to have size")
        }
        if(sum((names(param.tune[[which(methods=="LIB_PLANN")[i]]])%in%"decay"))==0){
          stop("Tune parameters for LIB_PLANN need to have decay")
        }
        if(sum((names(param.tune[[which(methods=="LIB_PLANN")[i]]])%in%"maxit"))==0){
          stop("Tune parameters for LIB_PLANN need to have maxit")
        }
        if(sum((names(param.tune[[which(methods=="LIB_PLANN")[i]]])%in%"MaxNWts"))==0){
          stop("Tune parameters for LIB_PLANN need to have MaxNWts")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")[i]]]$inter)|
             is.null(param.tune[[which(methods=="LIB_PLANN")[i]]]$inter))){
          stop("inter for LIB_PLANN need to be a scalar or a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")[i]]]$size)|
             is.null(param.tune[[which(methods=="LIB_PLANN")[i]]]$size))){
          stop("size for LIB_PLANN need to be a scalar, a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIPLA_PLANN")[i]]]$decay)|
             is.null(param.tune[[which(methods=="LIPLA_PLANN")[i]]]$decay))){
          stop("decay for LIB_PLANN need to be a scalar, a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")[i]]]$maxit)|
             is.null(param.tune[[which(methods=="LIB_PLANN")[i]]]$maxit))){
          stop("maxit for LIB_PLANN need to be a scalar, a vector or NULL")
        }
        if(!(is.numeric(param.tune[[which(methods=="LIB_PLANN")[i]]]$MaxNWts)|
             is.null(param.tune[[which(methods=="LIB_PLANN")[i]]]$MaxNWts))){
          stop("MaxNWts for LIB_PLANN need to be a scalar, a vector or NULL")
        }
      }
    }
  }



  if(length(.meth_rm)>=1){
    methods=methods[-.meth_rm]
    param.tune=param.tune[-.meth_rm]
  }

  if((max(ROC.precision)==1) | (min(ROC.precision)==0)){
    stop("values for ROC.precision need to be in ]0;1[")
  }


  if(is.null(param.weights.fix)==FALSE & is.null(param.weights.init)==FALSE){
    warning("Weights can not be fix and initial at the same time. SuperLearner ignored initial values")
    param.weights.init<-NULL
  }

  if(is.null(param.weights.fix)==FALSE | is.null(param.weights.init)==FALSE){
    if(is.null(param.weights.fix)==FALSE){
      if(is.numeric(param.weights.fix)==FALSE){
        stop("param.weights.fix need to be numeric")
      }
    }

    if(is.null(param.weights.init)==FALSE){
      if(is.numeric(param.weights.init)==FALSE){
        stop("param.weights.init need to be numeric")
      }

      if(length(param.weights.init)!=(length(methods)-1)){
        stop("wrong lenth for param.weights.init")
      }

    }
  }
  if(is.null(param.weights.fix)==TRUE & is.null(param.weights.init)==TRUE){
    param.weights.init<-rep(0,length(methods)-1)
  }


  if(sum(!(methods %in% c("LIB_COXlasso", "LIB_COXridge", "LIB_RSF", "LIB_COXen",
                          "LIB_AFTweibull","LIB_AFTweibull","LIB_AFTggamma","LIB_AFTgamma",
                          "LIB_PHgompertz","LIB_PHexponential", "LIB_PLANN",
                          "LIB_AFTllogis","LIB_COXaic","LIB_COXall", "LIB_PHspline")))>=1){
    stop("New method is not yet implemented") }



  M <-length((methods))
  N <- length(data[[times]])

  names.meth=c(rep(NA,M))
  for(i in unique(methods)){
    if(length(which(methods==i))==1){
      names.meth[which(methods==i)]=i
    }
    else{
      names.meth[which(methods==i)]=paste0(i,1:length(which(methods==i)))
    }
  }


  if(is.null(pro.time)==T & sum(metric %in% c("p_ci","bs","bll","ribs","ribll","auc","uno_ci"))){
    pro.time=median(data[[times]])
  }

  time.pred <- unique(sort(c(0,pro.time,data[[times]])))
  time.pred_bis <- unique(sort(c(0,pro.time,data[[times]])))

  if("LIB_PLANN" %in% methods){
    maxtime=max(time.pred)+1
  }





  if(is.null(param.tune)==FALSE){
    if(length(param.tune)!=M){
      stop("Param.tune need to have one element per method. Please modifiy param.tune or set it = NULL")
    }
    for (me in 1:M){
      if(is.null(param.tune[[me]])==T & !(methods[me] %in% c("LIB_AFTgamma",
                                                             "LIB_AFTggamma","LIB_AFTweibull","LIB_AFTllogis","LIB_PHexponential","LIB_PHgompertz","LIB_COXall","LIB_COXaic"))){
        if(methods[me] %in%"LIB_PHspline"){
          param.tune[[me]]=list(k=1:4)
        }
        if(methods[me] %in%"LIB_COXen"){
          param.tune[[me]]=list(alpha=seq(.1,.9,.1), lambda=NULL)
        }
        if(methods[me] %in%"LIB_PLANN"){
          param.tune[[me]]=list(inter=1,
                                size=c(2, 4, 6, 8, 10),
                                decay=c(0.001, 0.01, 0.02, 0.05),
                                maxit=100,
                                MaxNWts=10000)
        }
        if(methods[me] %in% "LIB_COXlasso"){
          param.tune[[me]]=list(lambda=NULL)
        }
        if(methods[me] %in%"LIB_COXridge"){
          param.tune[[me]]=list(lambda=NULL)
        }
        if(methods[me] %in%"LIB_RSF"){
          param.tune[[me]]=list(mtry=(length(variables_formula)-2)/2+2,
                                nodesize=c(2, 4, 6, 10, 20, 30, 50, 100),
                                ntree=500)
        }
      }
    }
  }


  if(is.null(param.tune)){
    param.tune=vector("list",M)
    for (me in 1:M){
      if(methods[me] %in%"LIB_COXen"){
        param.tune[[me]]=list(alpha=seq(.1,.9,.1), lambda=NULL)
      }
      if(methods[me] %in%"LIB_PLANN"){
        param.tune[[me]]=list(inter=1,
                              size=c(2, 4, 6, 8, 10),
                              decay=c(0.001, 0.01, 0.02, 0.05),
                              maxit=100,
                              MaxNWts=10000)
      }
      if(methods[me] %in% "LIB_COXlasso"){
        param.tune[[me]]=list(lambda=NULL)
      }
      if(methods[me] %in% "LIB_COXridge"){
        param.tune[[me]]=list(lambda=NULL)
      }
      if(methods[me] %in% "LIB_PHspline"){
        param.tune[[me]]=list(k=1:4)
      }
      if(methods[me] %in% "LIB_RSF"){
        param.tune[[me]]=list(mtry=seq(1,(length(variables_formula)-2)/2+2),
                              nodesize=c(2, 4, 6, 10, 20, 30, 50, 100), ntree=500)
      }
    }
  }



  .model<-vector("list",M)
  names(.model)<-names.meth

  .data_bis<-data
  .tune.optimal<-vector("list",M)
  names(.tune.optimal)<-names.meth

  .tune.results<-vector("list",M)
  names(.tune.results)<-names.meth


  for (me in 1:M){
    if(show_progress){ update_progress()}
    if(methods[me] == "LIB_AFTweibull" ){
      .tune.optimal[[me]]=NA

      .LIB_AFTweibull <- LIB_AFTweibull(formula=formula,data=data)

      .model[[me]]<-.LIB_AFTweibull


      rm(.LIB_AFTweibull) }

    if(methods[me] == "LIB_AFTggamma"){
      .tune.optimal[[me]]=NA

      .LIB_AFTggamma <- LIB_AFTggamma(formula=formula,data=data)

      .model[[me]]<-.LIB_AFTggamma

      rm(.LIB_AFTggamma) }

    if(methods[me] == "LIB_AFTgamma" ){
      .tune.optimal[[me]]=NA

      .LIB_AFTgamma <- LIB_AFTgamma(formula=formula,data=data)

      .model[[me]]<-.LIB_AFTgamma

      rm(.LIB_AFTgamma)  }


    if(methods[me] == "LIB_AFTllogis" ){
      .tune.optimal[[me]]=NA

      .LIB_AFTllogis <- LIB_AFTllogis(formula=formula,data=data)

      .model[[me]]<-.LIB_AFTllogis

      rm(.LIB_AFTllogis)    }

    if(methods[me] == "LIB_PHgompertz" ){
      .tune.optimal[[me]]=NA

      .LIB_PHgompertz <- LIB_PHgompertz(formula=formula,data=data)

      .model[[me]]<-.LIB_PHgompertz

      rm(.LIB_PHgompertz) }

    if(methods[me] == "LIB_PHexponential"){
      .tune.optimal[[me]]=NA

      .LIB_PHexponential <- LIB_PHexponential(formula=formula,data=data)

      .model[[me]]<-.LIB_PHexponential

      rm(.LIB_PHexponential) }

    if(methods[me] == "LIB_COXall"){
      .tune.optimal[[me]]=NA

      .LIB_COXall <- LIB_COXall(formula=formula,data=data)

      .model[[me]]<-.LIB_COXall

      rm(.LIB_COXall) }

    if(methods[me] == "LIB_COXaic"){
      .tune.optimal[[me]]=NA

      .LIB_COXaic <- LIB_COXaic(formula=formula,data=data,penalty=penalty)

      .model[[me]]<-.LIB_COXaic

      rm(.LIB_COXaic) }

    if(methods[me] == "LIB_PHspline"){

      if(is.null(param.tune[[me]]$k)==T | length(param.tune[[me]]$k)>1){

        .tune<- tunePHspline(formula=formula, metric=metric,pro.time=pro.time, data=data, cv=cv, seed=seed,
                             ROC.precision=ROC.precision,
                             k=param.tune[[me]]$k)
        .tune.optimal[[me]]=.tune$optimal
        .tune.results[[me]]=.tune$results
        rm(.tune)  }

      else{ .tune.optimal[[me]]=list(k=param.tune[[me]]$k) }


      .LIB_PHspline <- LIB_PHspline(formula=formula,data=data, k=.tune.optimal[[me]]$k)

      .model[[me]]<-.LIB_PHspline

      rm(.LIB_PHspline)    }

    if(methods[me] == "LIB_COXlasso"){

      if(is.null(param.tune[[me]]$lambda)==T | length(param.tune[[me]]$lambda)>1){

        .tune<- tuneCOXlasso(formula=formula, penalty=penalty, data=data, cv=cv, seed=seed,
                             parallel=FALSE, lambda=param.tune[[me]]$lambda)

        .tune.optimal[[me]]=.tune$optimal
        .tune.results[[me]]=.tune$results
        rm(.tune)  }

      else{ .tune.optimal[[me]]=list(lambda=param.tune[[me]]$lambda) }


      .LIB_COXlasso <- LIB_COXlasso(formula=formula,penalty=penalty, data=data, lambda=.tune.optimal[[me]]$lambda)

      .model[[me]]<-.LIB_COXlasso

      rm(.LIB_COXlasso)  }

    if(methods[me] == "LIB_COXridge"){

      if(is.null(param.tune[[me]]$lambda)==T | length(param.tune[[me]]$lambda)>1){
        .tune<- tuneCOXridge(formula=formula, penalty=penalty, data=data, cv=cv, seed=seed,
                             parallel=FALSE, lambda=param.tune[[me]]$lambda)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]=list(lambda=param.tune[[me]]$lambda)
      }

      .LIB_COXridge <- LIB_COXridge(formula=formula,penalty=penalty, data=data, lambda=.tune.optimal[[me]]$lambda)
      .model[[me]]<-.LIB_COXridge


      rm(.LIB_COXridge)    }

    if(methods[me] == "LIB_COXen"){

      if(length(param.tune[[me]]$alpha)==1 & length(param.tune[[me]]$lambda)==1){
        .tune.optimal[[me]]=list(alpha=param.tune[[me]]$alpha, lambda=param.tune[[me]]$lambda)
      }
      else{
        .tune <- tuneCOXen(formula=formula,penalty=penalty, data=data, cv=cv, seed=seed,
                           parallel=FALSE,
                           alpha=param.tune[[me]]$alpha,
                           lambda=param.tune[[me]]$lambda)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results

        rm(.tune) }

      .LIB_COXen <- LIB_COXen(formula=formula, penalty=penalty, data=data, alpha=.tune.optimal[[me]]$alpha,
                              lambda=.tune.optimal[[me]]$lambda)

      .model[[me]]<-.LIB_COXen

      rm(.LIB_COXen) }

    if (methods[me] == "LIB_RSF"){

      if(length(param.tune[[me]]$nodesize)!=1 | length(param.tune[[me]]$mtry)!=1 |
         length(param.tune[[me]]$ntree)!=1){
        .tune<-tuneRSF(formula=formula, data=data,
                       nodesize=param.tune[[me]]$nodesize,
                       mtry=param.tune[[me]]$mtry,
                       ntree=param.tune[[me]]$ntree,seed=seed)

        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]<-list(nodesize=param.tune[[me]]$nodesize,
                                  mtry=param.tune[[me]]$mtry,
                                  ntree=param.tune[[me]]$ntree)
      }
      .LIB_RSF <-LIB_RSF(formula=formula, data=data,
                         nodesize=.tune.optimal[[me]]$nodesize,
                         mtry=.tune.optimal[[me]]$mtry, ntree=.tune.optimal[[me]]$ntree, seed=seed)

      .model[[me]]<-.LIB_RSF

      rm(.LIB_RSF) }

    if (methods[me] == "LIB_PLANN"){

      if(length(param.tune[[me]]$inter)!=1 | length(param.tune[[me]]$size)!=1 |
         length(param.tune[[me]]$decay)!=1 | length(param.tune[[me]]$maxit)!=1 |
         length(param.tune[[me]]$MaxNWts)!=1){

        .tune <- tunePLANN(formula=formula, data=data, cv=cv,
                           inter=param.tune[[me]]$inter,
                           size=param.tune[[me]]$size,
                           decay=param.tune[[me]]$decay,
                           maxit=param.tune[[me]]$maxit,
                           MaxNWts=param.tune[[me]]$MaxNWts,
                           maxtime=maxtime,
                           seed=seed,metric=metric,pro.time=pro.time,ROC.precision=ROC.precision)
        .tune.optimal[[me]]<-.tune$optimal
        .tune.results[[me]]<-.tune$results
        rm(.tune)
      }
      else{
        .tune.optimal[[me]]<-list(inter=param.tune[[me]]$inter,
                                  size=param.tune[[me]]$size,
                                  decay=param.tune[[me]]$decay,
                                  maxit=param.tune[[me]]$maxit,
                                  MaxNWts=param.tune[[me]]$MaxNWts)
      }

      .LIB_PLANN <-LIB_PLANN(formula=formula, data=data,
                             inter=as.numeric(.tune.optimal[[me]]$inter),
                             size=as.integer(.tune.optimal[[me]]$size),
                             decay=as.numeric(.tune.optimal[[me]]$decay),
                             maxit=as.integer(.tune.optimal[[me]]$maxit),
                             maxtime=maxtime,
                             MaxNWts=as.integer(.tune.optimal[[me]]$MaxNWts))


      .model[[me]]<-.LIB_PLANN


      rm(.LIB_PLANN)    }


  }

  ########################
  ### Cross-Validation  ##
  ########################


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
      t_max_fold = t_max_fold,
      id=data[data$folds == i, "id"]
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
      set.seed(seed)

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
          t_max_fold = t_max_fold,
          id=data[data$folds == k, "id"]
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



  id_vect<- unlist(lapply(CVtune,function(x)(return(x$id))))
  # Calculate the maximum EVENT times per fold
  max_event_time_across_folds <- min(sapply(CVtune, function(fold) {
    # Only take the lines where failures == 1 (event)
    event_times <- fold$train[[times]][fold$train[[failures]] == 1]
    if(length(event_times) > 0) {
      return(max(event_times))
    } else {
      # If there is no event in this fold (rare case), take the maximum of all time
      warning("Fold without events detected, using max time instead")
      return(max(fold$train[[times]]))
    }
  }))


  # Adapt time.pred with this terminal
  time.pred <- time.pred[time.pred <= max_event_time_across_folds]



  CV_one_method<-function(meth, Tune, CV, penalty, formula, time.pred){
    if(show_progress){ update_progress()}
    if(meth == "LIB_AFTweibull"){
      fit<-LIB_AFTweibull(formula=formula, data=CV$train)
      pred=predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_AFTggamma"){
      fit<-LIB_AFTggamma(formula=formula, data=CV$train)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_AFTgamma"){
      fit<-LIB_AFTgamma(formula=formula, data=CV$train)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_AFTllogis"){
      fit<-LIB_AFTllogis(formula=formula, data=CV$train)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_PHgompertz"){
      fit<-LIB_PHgompertz(formula=formula, data=CV$train)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_PHexponential"){
      fit<-LIB_PHexponential(formula=formula, data=CV$train)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_PHspline"){
      fit<-LIB_PHspline(formula=formula, data=CV$train,
                        k=Tune$LIB_PHspline$k)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_COXlasso"){
      fit<-LIB_COXlasso(formula=formula, penalty=penalty, data=CV$train,
                        lambda=Tune$LIB_COXlasso$lambda)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth == "LIB_COXen"){
      fit<-LIB_COXen(formula=formula,penalty = penalty, data=CV$train,
                     alpha=Tune$LIB_COXen$alpha, lambda=Tune$LIB_COXen$lambda)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }
    if(meth =="LIB_COXridge"){
      fit<-LIB_COXridge(formula=formula, penalty=penalty, data=CV$train, lambda=Tune$LIB_COXridge$lambda)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions

    }
    if(meth =="LIB_COXaic"){
      fit<-LIB_COXaic(formula=formula, penalty=penalty, data=CV$train)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions

    }
    if(meth =="LIB_COXall"){
      fit<-LIB_COXall(formula=formula, data=CV$train)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions

    }
    if(meth =="LIB_RSF"){
      fit<-LIB_RSF(formula=formula,  data=CV$train,
                   nodesize=as.numeric(Tune$LIB_RSF$nodesize), mtry=as.numeric(Tune$LIB_RSF$mtry),
                   ntree=as.numeric(Tune$LIB_RSF$ntree),seed=seed)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions

    }

    if(meth =="LIB_PLANN"){
      fit<-LIB_PLANN(formula=formula, data=CV$train,
                     inter=as.numeric(Tune$LIB_PLANN$inter),
                     size=as.numeric(Tune$LIB_PLANN$size),
                     decay=as.numeric(Tune$LIB_PLANN$decay),
                     maxit=as.integer(Tune$LIB_PLANN$maxit),
                     MaxNWts=as.integer(Tune$LIB_PLANN$MaxNWts),
                     maxtime=maxtime)
      pred<-predict(fit, newtimes=time.pred, newdata=CV$valid)$predictions
    }

    return(survivals=pred)

  }



  cv_function<-function(meth, Tune, penalty, formula, time.pred, CVtune){
    result<-lapply(CVtune,CV_one_method, meth=meth, Tune=Tune,penalty=penalty, formula=formula, time.pred=time.pred)
    survivals<-do.call(rbind,result)
    return(list(survivals=survivals))

  }


  cv_all_meth<-list()

  for(i in 1:M){
    Tune<-.tune.optimal[i]
    meth<-methods[i]
    survivals<-((cv_function(meth,Tune,penalty,formula,time.pred,CVtune))$survivals)[order(id_vect),]
    cv_all_meth[[i]]<-survivals

  }



  ########################
  ### Fonction  ##
  ########################

  metric_weight_function<-function(param,list_surv,metric){
    param<-exp(param)/(1+sum(exp(param)))
    param<-c(param,1-sum(param))
    weighted_matrices <- mapply(function(mat, weight) mat * weight,list_surv,param, SIMPLIFY = FALSE)
    survivals.matrix <- Reduce("+", weighted_matrices)
    hazards.matrix<-NULL
    if(metric=="ll"){
      hazards.matrix<-t(apply(survivals.matrix,1,haz_function,times=time.pred))
    }
    result<-metrics(metric=metric,formula=formula, data=.data_bis, survivals.matrix=survivals.matrix, hazards.matrix=hazards.matrix,prediction.times=time.pred,pro.time=pro.time, ROC.precision=ROC.precision)
    if(metric %in% c("uno_ci","p_ci","auc","ll")){result<-(-result)}
    return(result)

  }




  if(is.null(param.weights.fix)==TRUE){
    set.seed(seed)
    if(show_progress){ update_progress()}
    estim<-optim(par=param.weights.init, fn=metric_weight_function, list_surv=cv_all_meth, metric=metric,
                 hessian=F,
                 method=optim.method,control=list(maxit=maxit))

    if(optim.local.min==T){
      start_par=estim$par
      estim<-optim(par=start_par, fn=metric_weight_function, list_surv = cv_all_meth,
                   metric=metric, hessian=F,
                   method=optim.method, control=list(maxit=maxit))
    }

    if(show_progress){ update_progress() }

  }








  ############################
  # Compute Survival from SL #
  ############################

  if(is.null(param.weights.fix)==FALSE){
    estim=list()
    estim$par=param.weights.fix
  }


  rm(cv_all_meth)

  w.sl <- c(exp(c(estim$par,0)) / ( 1+sum(exp(estim$par))) )

  FitALL<- vector("list",M)

  for(me in 1:M){
    if(show_progress){ update_progress()}
    FitALL[[me]]<-predict(object= .model[[me]],newtimes=time.pred_bis)$predictions
  }

  weighted_matrices <- mapply(function(mat, weight) mat * weight, FitALL, w.sl, SIMPLIFY = FALSE)
  survivals <- Reduce("+", weighted_matrices)

  ######################
  # PREPARATION RETURN #
  ######################




  value<-estim$value

  if(metric %in% c("uno_ci","p_ci","auc","ll")){
    value<-(-estim$value)
  }

  data<-.data_bis

  res<-list(times=time.pred_bis,
            predictions=survivals,
            FitALL=FitALL,
            data=data,
            formula=formula,
            ROC.precision=ROC.precision,
            cv=cv,
            pro.time=pro.time,
            methods=methods,
            models=.model,
            weights=list(coefficients=estim$par, values=w.sl),
            metric=list(metric=metric, value=value),
            seed=seed,
            optim.method=optim.method,
            param.tune=list(optimal=.tune.optimal, results=.tune.results))

  class(res) <- "sltime"
  if(show_progress){ update_progress()}

  return(res)
}




