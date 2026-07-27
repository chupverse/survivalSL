metrics <- function(metric, formula=NULL, data=NULL, survivals.matrix=NULL, hazards.matrix=NULL,
                    prediction.times=NULL, object=NULL, pro.time=NULL, ROC.precision=seq(.01, .99, by=.01)){


  # The name of the metric must be provided.
  if (missing(metric)) stop("The 'metric' argument is required.")
  if (! (is.character(metric))) stop("The 'metric' argument must be character.")



  # Check the consistency of the dimensions
  if (!is.null(hazards.matrix)) {
    if (ncol(hazards.matrix) != length(prediction.times)) {
      stop("Number of columns in hazards.matrix must match length(prediction.times).")
    }
  }

  if (!is.null(survivals.matrix)) {
    if (ncol(survivals.matrix) != length(prediction.times)) {
      stop("Number of columns in survivals.matrix must match length(prediction.times).")
    }
  }


  # If there is no object: then the rest is needed!
  if(is.null(object)){
    if (is.null(formula)) stop("Because the 'object' is missing, the 'formula' argument is required.")
    if (is.null(survivals.matrix)) stop("Because the 'object' is missing, the 'survivals.matrix' argument is required.")
    if (is.null(prediction.times)) stop("Because the 'object' is missing, the 'prediction.times' argument is required.")
    if (is.null(data)) stop("Because the 'object' is missing, the 'data' argument is required.")
     if((metric=="ll") & is.null(hazards.matrix))stop("The 'hazards.matrix' argument is required to compute the log-likelihood.")

  }else{
    if (is.null(data)){data<-object$data}
    if (is.null(formula)){formula<-object$formula}
  }

  # We retrieve the names of the variables that interest us:

  variables_formula <- all.vars(formula)

  times <- variables_formula[1]
  failures <- variables_formula[2]


  # We need a pro.time for certain metrics, so we take the median time
  if(is.null(pro.time) & !(metric %in% c("ll","ibll","ibs"))) {pro.time <- median(data[[times]])}


  # The metric must be one of the following:

  if(min(metric %in% c("uno_ci","p_ci","bs","ibs","ibll","bll", "ribs","ribll","auc","ll"))==0){
    stop("The argument \"metric\" must be Brier score (bs),
         Pencina concordance index (p_ci),
         Uno concordance index (uno_ci),
         integrated Brier score (ibs), the binomial log-likelihood (bll),
         the integrated binomial log-likelihood (ibll), the restricted ibs (ribs),
         the restricted ibll (ribll),
         the area under the ROC curve (auc), or the log-likelihood (ll)")
  }



  if((!(is.null(object))) & (is.null(survivals.matrix)|is.null(hazards.matrix))){

    # If there is no prediction.times, we take the timestamps from our object:
    if(is.null(prediction.times)){prediction.times<-sort(unique(c(pro.time,object$times)))}
    # This gives us our matrix with the survival probabilities for each of our individuals:
    survivals.matrix<-predict(object,newdata=data,newtimes=prediction.times)$predictions

    if((as.character(object$model$call[1]) %in% c("flexsurvreg","flexsurvspline"))& metric=="ll"){
      .haz<-predict(
        object$model,
        newdata=data,
        type = "haz",
        times = prediction.times
      )

      hazards.matrix <- t(sapply(.haz$.pred, function(x) x[[2]]))

    }

    if(object$model$call[1]!="flexsurvreg()"& metric=="ll"){
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


      hazards.matrix<-t(apply(survivals.matrix,1,haz_function,times=prediction.times))


    }

  }



    data.times <- data[[times]]  # our times
    data.failures <- data[[failures]] # our events
    obj_surv <- Surv(data.times, data.failures)

    time <- obj_surv[, 1]
    ot <- order(time)
    cens <- obj_surv[ot, 2]
    time <- time[ot]

    .temp <- survfit(Surv(time, cens==0) ~ 1)
    .csurv <- summary(.temp, times=time, extend=TRUE)$surv  # proba of censorship for all data times

    .csurv[.csurv == 0] <- Inf

    if(!(is.null(pro.time))){
      .csurv.pro.time <- summary(.temp, times=pro.time, extend=TRUE)$surv
      .csurv.pro.time[.csurv.pro.time == 0] <- Inf
    }

    survivals.matrix <- survivals.matrix[ot,] # reorder our survivals matrix

    hazards.matrix <- hazards.matrix[ot,] # same for our hazards matrix

    .data_bis<-data
    data<-data[ot,] # same with our data


    switch(metric,
           p_ci={
             status <- ifelse((data[[failures]]==1) & (data[[times]]>pro.time),0,data[[failures]])
             time <- ifelse(data[[times]]>pro.time,pro.time,data[[times]])
             length<-dim(survivals.matrix[,prediction.times>=pro.time])
             if(is.null(length)){
               predicted<-survivals.matrix[,prediction.times>=pro.time]
             }else{
               predicted<-survivals.matrix[,prediction.times>=pro.time][,1]
             }
             permissible <- 0 # comparable pairs
             concord <- 0
             n <- length(time)
             for (i in 1:n) {
               for (j in (i + 1):n) {
                 if(time[i] < time[j] & status[i] == 1 ){
                   permissible<-permissible + 1
                   if(predicted[i] < predicted[j]){
                     concord <- concord + 1
                   }
                   if(predicted[i]==predicted[j]){
                     concord <- concord + 0.5
                   }
                 }

                 if(status[i] == 1 & status[j] == 0 & time[i]==time[j] ){
                   permissible<-permissible + 1
                   if(predicted[i] < predicted[j]){
                     concord <- concord + 1
                   }
                 }



               }

             }
             RET<-ifelse(permissible==0,0,concord/permissible)

           },


           ibs={
             bsc <- sapply(1:length(prediction.times), FUN = function(j)
             {pro.time<-prediction.times[j]
             .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
             .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
             score_risque<-survivals.matrix[,j]
             help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
             help2 <- ifelse(data[[times]] > pro.time ,1,0)
             bs <- mean((0 - score_risque)^2 * help1 * (1/.csurv) + (1 - score_risque)^2 * help2 * (1/.csurv.pro.time))
             bs <- as.numeric(bs)
             return(bs)

             })

             idx <- 2:length(prediction.times)
             RET <- diff(prediction.times) %*% ((bsc[idx - 1] + bsc[idx])/2)
             RET <- RET/diff(range(prediction.times))
             RET <- as.matrix(RET)
           },
           bs={
             length<-dim(survivals.matrix[,prediction.times>=pro.time])
             if(is.null(length)){
               score_risque<-survivals.matrix[,prediction.times>=pro.time]
             }else{
               score_risque<-survivals.matrix[,prediction.times>=pro.time][,1]
             }
             help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
             help2 <- ifelse(data[[times]] > pro.time,1,0)
             bs <- mean((0 - score_risque)^2 * help1 * (1/.csurv) + (1 - score_risque)^2 * help2 * (1/.csurv.pro.time))
             RET <- as.numeric(bs)

           },
           bll={
             length<-dim(survivals.matrix[,prediction.times>=pro.time])
             if(is.null(length)){
               pro.time.surv<-survivals.matrix[,prediction.times>=pro.time]
             }else{
               pro.time.surv<-survivals.matrix[,prediction.times>=pro.time][,1]
             }
             pro.time.surv[which(pro.time.surv==0)]<-10**-7
             pro.time.surv[which(pro.time.surv==1)]<-1-10**-7
             help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
             help2 <- ifelse(data[[times]] > pro.time,1,0)
             bll <- -mean(log(1-pro.time.surv) * help1 * (1/.csurv) + log(pro.time.surv) * help2 * (1/.csurv.pro.time))
             RET <- as.numeric(bll)
           },
           ibll={
             bllc <- sapply(1:length(prediction.times), FUN = function(j)
             {pro.time<-prediction.times[j]
             .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
             .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
             score_risque<-unlist(survivals.matrix[,j])
             score_risque[which(score_risque < 1e-15)]<-1e-7
             score_risque[which(score_risque > 1-1e-15)]<-1-1e-7
             help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
             help2 <- ifelse(data[[times]] > pro.time ,1,0)
             bll <- -mean(log(1-score_risque) * help1 * (1/.csurv) + log(score_risque) * help2 * (1/.csurv.pro.time))
             return(as.numeric(bll)) })
             idx <- 2:length(prediction.times)
             RET <- diff(prediction.times) %*% ((bllc[idx - 1] + bllc[idx])/2)
             RET <- RET/diff(range(prediction.times))
             RET <- as.matrix(RET)
           },

           ribll={
             prediction.times <- prediction.times[prediction.times <= pro.time]
             survivals.matrix<-survivals.matrix[,prediction.times<=pro.time]
             bllc <- sapply(1:length(prediction.times), FUN = function(j)
             {pro.time<-prediction.times[j]
             .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
             .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
             score_risque<-survivals.matrix[,j]
             score_risque[which(score_risque < 1e-15)]<-1e-7
             score_risque[which(score_risque > 1-1e-15)]<-1-1e-7
             help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
             help2 <- ifelse(data[[times]] > pro.time ,1,0)
             bll <- -mean(log(1-score_risque) * help1 * (1/.csurv) + log(score_risque) * help2 * (1/.csurv.pro.time))
             return(as.numeric(bll)) })
             idx <- 2:length(prediction.times)
             RET <- diff(prediction.times) %*% ((bllc[idx - 1] + bllc[idx])/2)
             RET <- RET/diff(range(prediction.times))
             RET <- as.matrix(RET)

           },

           ribs={
             prediction.times <- prediction.times[prediction.times <= pro.time]
             survivals.matrix<-survivals.matrix[,prediction.times<=pro.time]
             bsc <- sapply(1:length(prediction.times), FUN = function(j)
             {pro.time<-prediction.times[j]
             .csurv.pro.time<-summary(.temp, times=pro.time, extend=TRUE)$surv
             .csurv.pro.time<-ifelse(.csurv.pro.time==0,Inf,.csurv.pro.time)
             score_risque<-survivals.matrix[,j]
             help1 <- ifelse(data[[times]] <= pro.time & data[[failures]] == 1,1,0)
             help2 <- ifelse(data[[times]] > pro.time ,1,0)
             bs <- mean((0 - score_risque)^2 * help1 * (1/.csurv) + (1 - score_risque)^2 * help2 * (1/.csurv.pro.time))
             bs <- as.numeric(bs)
             return(bs)

             })

             idx <- 2:length(prediction.times)
             RET <- diff(prediction.times) %*% ((bsc[idx - 1] + bsc[idx])/2)
             RET <- RET/diff(range(prediction.times))
             RET <- as.matrix(RET)
           },
           auc={
             length<-dim(survivals.matrix[,prediction.times>=pro.time])
             if(is.null(length)){
               pro.time.surv<-survivals.matrix[,prediction.times>=pro.time]
             }else{
               pro.time.surv<-survivals.matrix[,prediction.times>=pro.time][,1]
             }
             .data <- data.frame(times=data[[times]], failures=data[[failures]],
                                 variable=(1-pro.time.surv))

             roc <- function(times, failures, variable, confounders, data, pro.time, precision=seq(.01, .99, by=.01)) # FROM RISCA PACKAGE
             {

               # Getting cut-offs from precision sequence
               cut.off <- quantile(data[[variable]], probs=precision, na.rm=TRUE)

               if((max(precision)==1) | (min(precision)==0)){ stop("The cut-off values have to be different from the minimum or the maximum of the variable") }

               data$temp <- data[[failures]] + data[[variable]] + data[[times]]



               form0 <- update.formula(confounders, temp ~ .)

               if((length(data[[failures]]) - summary(glm(form0, data=data))$df.null - 1) > 0) {stop("Missing values are not allowed")}



               # To calculate censorship times for each individual in Ti
               km.c <- summary(survfit(Surv(data[[times]], data[[failures]]==0)~1))
               km.c <- data.frame(time=c(0, km.c$time), surv=c(1, km.c$surv))
               km.c$surv[km.c$surv==0] <- 0.0001
               s.c <- 1 / (sapply(data[[times]], FUN=function(x) {km.c$surv[km.c$time<=x][sum(km.c$time<=x)]}))

               .cox <- eval(parse(text=paste("coxph(Surv(",times ,", ",failures,")",paste(confounders[1],confounders[2]),", data=data)",sep="")))


               if(confounders=="~1"){
                 survfit.object <- survival::survfit(.cox,se.fit=FALSE,conf.int=FALSE)
                 W <- summary(survfit.object,times=pro.time)$surv
               }else{
                 survfit.object <- survival::survfit(.cox, newdata=data,se.fit=FALSE,conf.int=FALSE)
                 W <- summary(survfit.object, times=pro.time)$surv
               }

               se <- function(x) {
                 sum(1 * ((data[[variable]] > cut.off[x]) & (data[[times]] <= pro.time)) * data[[failures]] * s.c *
                       pmin(1/(1 - W),1/(mean(1-W)^2*0.1)))/sum(1 * (data[[times]] <= pro.time) *
                                                                  data[[failures]] * s.c * pmin(1/(1 - W),1/(mean(1-W)^2*0.1)))
               }

               sp <- function(x) {
                 sum(1 * ((data[[variable]] <= cut.off[x]) & (data[[times]] > pro.time)) *
                       pmin(1/W,1/(mean(W)^2*0.1)))/sum(1 * (data[[times]] > pro.time) * pmin(1/W,1/(mean(W)^2*0.1)))
               }

               temp.se <- sapply(1:length(cut.off), FUN = "se")
               temp.sp <- sapply(1:length(cut.off), FUN = "sp")

               .tab <-data.frame(cut.off = cut.off, se = temp.se, sp1 = 1-temp.sp)

               .tab$se[.tab$se > 1] <- NA
               .tab$sp1[.tab$sp1 > 1] <- NA
               .tab.res.temp <- .tab[!is.na(.tab$sp1 + .tab$se), ]
               .tab.res.temp <- .tab.res.temp[order(.tab.res.temp$sp1, .tab.res.temp$se), ]
               .tab.res.temp <- rbind(c(NA, 0, 0), .tab.res.temp, c(NA, 1, 1))
               colnames(.tab.res.temp) <- colnames(.tab)

               .tab$sp <- 1 - .tab$sp1
               .tab <- rbind(c(min(data[[variable]]),1,1,0) ,.tab, c(max(data[[variable]]),0,0,1))



               auc.fun <- function(sens, spec)
               {
                 .tab.res <- data.frame(se=sens, sp=spec)
                 .tab.res <- .tab.res[!is.na(.tab.res$sp + .tab.res$se),]
                 .tab.res$sp1 <- 1-.tab.res$sp
                 .tab.res <- .tab.res[order(.tab.res$sp1, .tab.res$se),]
                 .tab.res <- rbind(c(0,1,0), .tab.res, c(1,0,1))

                 # It's just the formula for the area of a trapezoid
                 return( sum((.tab.res$sp1[2:length(.tab.res$sp1)] -
                                .tab.res$sp1[1:(length(.tab.res$sp1)-1)]) * 0.5 *
                               (.tab.res$se[2:length(.tab.res$se)]+.tab.res$se[1:length(.tab.res$se)-1])) )
               }


               if(dim(.tab.res.temp)[2]>2){
                 .auc <- auc.fun(.tab.res.temp$se, 1-.tab.res.temp$sp1)
               }else{.auc<-NA}


               .tab$J <- .tab$sp + .tab$se - 1

               .obj <- list(table=.tab[,c("cut.off", "se", "sp", "J")], auc = .auc)

               class(.obj) <- "rocrisca"

               return(.obj)
             }



             RET <- roc(times="times", failures="failures", variable="variable",
                        confounders=~1, data=.data,
                        pro.time=pro.time, precision=ROC.precision)$auc
           },

           uno_ci={
             status <- ifelse((data[[failures]]==1) & (data[[times]]>pro.time),0,data[[failures]])
             time <- ifelse(data[[times]]>pro.time,pro.time,data[[times]])
             length<-dim(survivals.matrix[,prediction.times>=pro.time])
             if(is.null(length)){
               predicted<-survivals.matrix[,prediction.times>=pro.time]
             }else{
               predicted<-survivals.matrix[,prediction.times>=pro.time][,1]
             }
             permissible <- 0 # comparable pairs
             concord <- 0
             n <- length(time)
             for (i in 1:n) {
               for (j in (i + 1):n) {
                 if(time[i] < time[j] & status[i] == 1 ){
                   permissible<-permissible + 1/(.csurv[i])**2
                   if(predicted[i] < predicted[j]){
                     concord <- concord + 1/(.csurv[i])**2
                   }
                   if(predicted[i]==predicted[j]){
                     concord <- concord + 0.5*(1/(.csurv[i])**2)
                   }
                 }

                 if(status[i] == 1 & status[j] == 0 & time[i]==time[j] ){
                   permissible<-permissible + 1/(.csurv[i])**2
                   if(predicted[i] < predicted[j]){
                     concord <- concord + 1/(.csurv[i])**2
                   }
                 }



               }

             }
             RET<-ifelse(permissible==0,0,concord/permissible)

           },

           ll={
             # Copy the input dataset
             .Data <- data

             # Create an id column to keep track of the original row order
             # (because we will temporarily reorder the dataset)
             .Data$id <- 1:nrow(.Data)

             # Sort the dataset by observed times
             # This makes it easier to identify identical consecutive times
             .Data <- .Data[order(.Data[[times]]), ]

             # Map each UNIQUE observed time to the corresponding prediction time index
             # findInterval returns the index of the largest prediction.time <= observed time
             .idx <- findInterval(unique(.Data[[times]]), prediction.times)

             # Initialize a counter that will represent the rank of unique times
             .b <- 1

             # Build an index indicating which unique time each observation corresponds to
             # The goal is to map each row to the correct column in the matrices
             .result <- sapply(2:nrow(.Data), function(i){

               value <- .Data[[times]][i]

               # If the time changes compared to the previous observation,
               # we move to the next unique time rank
               if (value != .Data[[times]][i-1]) {
                 .b <<- .b + 1
               }

               return(.b)
             })

             # Add the first observation rank (since the loop starts at row 2)
             .result <- c(1, .result)

             # rank_time indicates which unique observed time corresponds to each row
             # Example:
             # times      = 1.2 1.2 2.7 4 4
             # rank_time  = 1   1   2   3 3
             #
             # This allows us to match each observation with the correct column
             # in the matrices built from unique times.
             .Data$rank_time <- .result

             # Restore the original row order of the dataset
             .Data <- .Data[order(.Data$id), ]

             # Keep only the columns of the survival matrix corresponding to
             # the prediction times matched to the unique observed times
             .survivals.matrix_bis <- survivals.matrix[, .idx]

             # Same operation for the hazard matrix
             .hazards.matrix_bis <- hazards.matrix[, .idx]


             # Extract the survival probability corresponding to each observation's time
             # rank_time tells us which column to use
             .Data$surv <- sapply(1:nrow(.Data), function(i){
               .b <- .survivals.matrix_bis[i, (.Data$rank_time[i])]
               return(.b)
             })

             # Extract the hazard value corresponding to each observation's time
             .Data$haz <- sapply(1:nrow(.Data), function(i){
               .b <- .hazards.matrix_bis[i, (.Data$rank_time[i])]
               return(.b)
             })

             # If hazard is NA, set the failure indicator to 0
             # (the event cannot contribute to the likelihood if hazard is missing)
             .Data[[failures]] <- ifelse(is.na(.Data$haz), 0, .Data[[failures]])

             # Replace missing hazard values with a very small number
             # to avoid numerical issues in log()
             .Data$haz <- ifelse(is.na(.Data$haz), 10^-7, .Data$haz)

             # Replace infinite hazard values with a small number as well
             .Data$haz <- ifelse(is.infinite(.Data$haz), 10^-7, .Data$haz)


             # Compute the log-likelihood contribution:
             # sum( event * log(hazard) + log(survival) )
             RET <- sum(.Data[[failures]] * log(.Data$haz) + log(.Data$surv))

           })
    data<-.data_bis
    return(as.numeric(RET))


}





