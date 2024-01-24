load.libraries <- function(){
  
  # Baselina libraries
  library(tidyverse)
  library(L1pack)
  library(latex2exp)
  library(xtable)
  
  ## Libraries for Bykhovskaya.regression
  #install.packages("nloptr")
  library(nloptr)
  library(matlib)
  library(stats)
  library(expm)
  
}

## Load and preprocess data
load.vol <- function(){
  folder <- "C:/data/"
  vol <- read.csv(paste0(folder,"rk_1m.txt"), header = FALSE)
  colnames(vol) <- unlist(read.csv(paste0(folder,"ids.txt"), header = FALSE))
  preserve <- colnames(read.csv("C:/conectedness_2net_rk.txt", header = TRUE))[-1]
  vol <- vol[,preserve]
  return(vol)
}
load.net <- function(){
  net <- read.csv("C:/matlab/results/LF/conectedness_0network_rk.txt", header = FALSE)
  net <- as.data.frame(t(net))
  days <- read.csv("C:/matlab/data/daysplit.txt", header = FALSE)
  rownames(net) <- as.Date(as.character(days[,1]), format = "%d-%m-%Y")
  rm(days)
  lbls <- as.character(read.csv("C:/results/LF/conectedness_0lbls.txt", header = FALSE))
  cns <- rep("",length(lbls)*length(lbls))
  for(i in 1:length(lbls))for(j in 1:length(lbls)) cns[(i-1)*10+j] <- paste0(lbls[j],",",lbls[i])
  colnames(net) <- cns
  return(net)
}
load.mkt <- function(net){
  mkt <- read.csv("C:/data/MKT.csv", header = TRUE)
  mkt <- mkt%>%mutate(MKT = Mkt.RF+RF,X=parse_date(as.character(X), format = "%Y%m%d"))%>%select(X,MKT)%>%mutate(nMKT=1*(MKT<0))
  rownames(mkt) <- mkt$X
  mkt <- mkt%>%select(MKT,nMKT)
  mkt$MKT <- mkt$MKT / sd(mkt$MKT, na.rm = TRUE)
  mkt <- mkt[rownames(mkt)%in%rownames(net),]
  return(mkt)
}
load.sentiment <- function(net,market.wide = TRUE){
  sentiment <- read.csv("C:/data/sentiment_ecsector_daily_ew.csv")
  rownames(sentiment) <- as.Date(sentiment[,1], format = "%d.%m.%Y")
  sentiment <- sentiment[,-1]
  sentiment <- sentiment[rownames(sentiment)%in%rownames(net),]
  #sentiment <- sentiment%>% mutate(across(where(is.numeric), scale))
  #pSentiment <- nSentiment <- sentiment / sd(unlist(sentiment))
  #pSentiment[pSentiment<0] <- 0
  #nSentiment[pSentiment>0] <- 0
  #nSentiment <- abs(nSentiment)
  if(market.wide){
    sentiment$TOT <- rowSums(sentiment)
    sentiment$TOT <- sentiment$TOT / sd(sentiment$TOT)
    return(sentiment%>%select(TOT))
  }else{
    sentiment <- sentiment / sd(unlist(sentiment))
    return(sentiment)
  }
}
construct.data <- function(net,mkt,sentiment){
  
  # 2.1 Construct baseline data set: first two network lags
  mdl.data <- data.frame(net[rownames(net)[-c(1,2)],colnames(net)[2]])
  colnames(mdl.data)[1] <- paste0("X_",colnames(net)[2],"_1")
  peers <- list()
  for(lag in 1:2){
    for(edge in colnames(net)){
      cmpnents <- unlist(strsplit(edge, split=","))
      peers <- c(peers,cmpnents[1],cmpnents[2])
      if(length(cmpnents)==2&&cmpnents[2]!=cmpnents[1]){
        if(lag==1){
          mdl.data[,paste0("X_",edge,"_",lag)] <- data.frame(net[rownames(net)[-c(1,nrow(net))],edge])
        }else if(lag==2){
          mdl.data[,paste0("X_",edge,"_",lag)] <- data.frame(net[rownames(net)[-c(nrow(net)-1,nrow(net))],edge])
        }
      }
    }
  }
  peers <- unique(unlist(peers))
  n <- length(peers) # n <- (1+sqrt(1+4*ncol(mdl.data)/2))/2
  
  # 2.2 Add (lag-1) triangular peer effects
  col_list <- list()
  for(colmn in colnames(mdl.data)) if(unlist(strsplit(colmn,"_"))[3]=="1") col_list <- c(col_list,colmn)
  col_list <- unlist(col_list)
  for(edge in col_list){
    edge_info <- unlist(strsplit(edge,"_"))
    cmponents <- unlist(strsplit(edge_info[2],","))
    cn <- paste(edge_info[1],edge_info[2],"tr", sep = "_")
    mdl.data[,cn] <- 0
    for(intermediate_peer in peers) if(intermediate_peer != cmponents[1] && intermediate_peer != cmponents[2]){
      c1 <- paste0("X_",cmponents[1],",",intermediate_peer,"_1")
      c2 <- paste0("X_",intermediate_peer,",",cmponents[2],"_1")
      mdl.data[,cn] <- mdl.data[,cn] + sqrt(mdl.data[,c1]*mdl.data[,c2]) / (n - 2)
      #print(paste0(cn," += ",c1," * ",c2))
    }
  }
  cutoffs <- (1:3) * (n*(n-1))
  
  # 2.3 Add sentiment data 
  #mdl.data[,paste0("S_plus_",colnames(pSentiment))] <- pSentiment[-c(1,nrow(net)),]
  #mdl.data[,paste0("S_minus_",colnames(nSentiment))] <- nSentiment[-c(1,nrow(net)),]
  #cutoffs <- c(cutoffs,ncol(mdl.data))
  
  # 2.4 Add quadratic sentiment data
  #mdl.data[,paste0("S2_plus_",colnames(pSentiment))] <- pSentiment[-c(1,nrow(net)),]^2
  #mdl.data[,paste0("S2_minus_",colnames(nSentiment))] <- nSentiment[-c(1,nrow(net)),]^2
  #cutoffs <- c(cutoffs,ncol(mdl.data))
  
  # 2.5 Add interactions with information (sign of MTK return) variable
  #mdl.data[,"I"] <- mkt$nMKT[-c(1,nrow(net))]
  #mdl.data[,paste0("IS_plus_",colnames(pSentiment))] <- mkt$nMKT[-c(1,nrow(net))]*pSentiment[-c(1,nrow(net)),]
  #mdl.data[,paste0("IS_minus_",colnames(nSentiment))] <- mkt$nMKT[-c(1,nrow(net))]*nSentiment[-c(1,nrow(net)),]
  #mdl.data[,paste0("IS2_plus_",colnames(pSentiment))] <- mkt$nMKT[-c(1,nrow(net))]*pSentiment[-c(1,nrow(net)),]^2
  #mdl.data[,paste0("IS2_minus_",colnames(nSentiment))] <- mkt$nMKT[-c(1,nrow(net))]*nSentiment[-c(1,nrow(net)),]^2
  
  # 2.3-2.5 Add market data, sentiment data, and interactions
  mdl.data[,"MKT"] <- mkt$MKT[-c(1,nrow(net))]
  mdl.data <- mdl.data%>%mutate(aMKT = abs(MKT), MKT2 = MKT^2)
  mdl.data[,"S"] <- sentiment$TOT[-c(1,nrow(net))]
  mdl.data <- mdl.data%>%mutate(aS = abs(S), S2 = S^2)
  mdl.data <- mdl.data%>%mutate(MMT_S = MKT*S, MMT_aS = MKT*aS, MMT_S2 = MKT*S2)
  
  # 2.6 Retain only non-null data
  cutoffs <- c(cutoffs,ncol(mdl.data))
  keep <- rowSums(is.na(mdl.data))==0
  
  # Return results
  return(list(mdl.data=mdl.data,cutoffs=cutoffs,keep=keep,peers=peers,n=n))
  
}

## Implementation of Anna Bykhovskaya (2023) Time Series Approach to the Evolution of Networks: Prediction and Estimation, Journal of Business & Economic Statistics, 41:1, 170-183, DOI: 10.1080/07350015.2021.2006669
Bykhovskaya.regression <- function(dta, maxiter = 2.5*10^4, fit.lower.bound = 0, fit.upper.bound = Inf){
  
  # Define data
  y <- dta[,"y"]
  X <- dta
  X[,"y"] <- 1
  X <- X%>%rename("(Intercept)"="y")
  n <- ncol(X)
  T_obs <- nrow(X)
  
  # Define fit function
  fit.f <- function(X,b) return(as.numeric(rowSums(sweep(X, MARGIN=2,b, `*`))))
  
  # Define objective function
  eval_f <- function(b){
    return(sum(abs(y-pmax(rep(0,length(y)),as.numeric(rowSums(sweep(X, MARGIN=2,b, `*`)))))))
  }
  
  # Define constraint function
  eval_g_ineq <- function(b){
    return(c(fit.lower.bound-min(as.numeric(rowSums(sweep(X, MARGIN=2,b, `*`)))), max(as.numeric(rowSums(sweep(X, MARGIN=2,b, `*`))))-fit.upper.bound))
    #return(fit.lower.bound-min(as.numeric(rowSums(sweep(X, MARGIN=2,b, `*`)))))
  }
  
  # Initialize parameters based on standard LTD regression estimates
  mdl2 <- lad(y~., data = dta)
  coef <- as.numeric(mdl2$coefficients)
  se <- sqrt(diag(vcov(mdl2)))
  tstat <- coef / se
  pval <- 2 * pt(abs(tstat), df = nrow(dta) - nrow(mdl2$R) - 2, lower.tail = FALSE)
  res.b <- data.frame(coef, se, tstat, pval)
  rownames(res.b) <- names(mdl2$coefficients)
  coef.b <- coef
  rm(coef,se,tstat,pval)
  
  # Perform optimization
  # See https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  start <- Sys.time()
  optmodel <- nloptr::nloptr(x0=coef.b,                         
                             eval_f = eval_f,
                             eval_g_ineq = eval_g_ineq,
                             lb = rep(-Inf,length(coef.b)),
                             ub = rep(Inf,length(coef.b)), 
                             opts = list("algorithm"="NLOPT_LN_COBYLA", # "NLOPT_LN_BOBYQA"
                                         "xtol_rel" = 1.0e-10,"ftol_abs" = 1.0e-10,
                                         "maxeval" = maxiter))
  end <- Sys.time()
  
  # Coefficients and fit
  coef <- as.numeric(optmodel$solution)
  fit <- fit.f(X, coef)
  
  # Objective and constraint function values
  f <- eval_f(coef)
  g <- eval_g_ineq(coef)
  
  # Coef Cov matrix
  X.adj <- X
  X.adj[fit<0,] <- 0
  M0_hat <- matrix(nrow = n, ncol = n)
  for(i in 1:n) for(j in 1:n) M0_hat[i,j] <- mean(X.adj[,i]*X.adj[,j])
  M0_sq_inv <- chol2inv(chol(sqrtm(M0_hat)))
  
  if(min(mdl2$fitted.values)<fit.lower.bound||max(mdl2$fitted.values)>fit.upper.bound){
    fu <- density(y-fit, kernel = "rectangular")
    i <- sum(fu$x<=0)
    fu0 <- fu$y[i] -fu$x[i]/(fu$x[i+1] - fu$x[i])*(fu$y[i+1] - fu$y[i])
    scale <- 1/2/fu0
  }else scale <- mdl2$scale
  
  #vcov(mdl2)
  #M0_sq_inv%*%M0_sq_inv * 2 * mdl2$scale^2 / T_obs
  #se <- sqrt(diag(M0_sq_inv%*%M0_sq_inv * 2 * mdl2$scale^2 / T_obs))
  #sqrt(diag(vcov(mdl2)))
  
  # se, tstat, pvals
  R <- M0_sq_inv%*%M0_sq_inv * 2 * scale^2 / T_obs
  se <- sqrt(diag(R))
  tstat <- coef / se
  pval <- 2 * pt(abs(tstat), df = nrow(dta) - ncol(dta) - 2, lower.tail = FALSE)
  coef.summary <- data.frame(coef,se,tstat,pval)
  rownames(coef.summary) <- names(mdl2$coefficients)
  
  # Return results
  return(list(
    opt.time=end-start,
    converged=optmodel$iterations<maxiter&&sum(g>0)==0,
    coeficients=list(opmod=coef.summary,benchmod=res.b),
    Cov=R,
    scale=scale,
    fit=fit,
    resid=y-fit,
    SumAbsResid=f,
    IneqCons=g))
}
run.Bykhovskaya.regression <- function(net,mdl.data,keep,n,X.cols = NA){
  
  if(is.na(X.cols[1])) X.cols <- list(full=(1:ncol(mdl.data)))
  
  filename <- "C:/results/coefficients.bench"
  if(file.exists(filename)){
    load("C:/results/actual.bench")
    load("C:/results/fit.bench")
    load(filename)
  }else{
    row <- 0
    actual <- data.frame(matrix(nrow = sum(keep), ncol = n * (n-1)))
    fit <- list()
    fit[["base"]] <- data.frame(matrix(nrow = sum(keep), ncol = n * (n-1)))
    for(mdl.name in names(X.cols)) fit[[mdl.name]] <- data.frame(matrix(nrow = sum(keep), ncol = n * (n-1)))
    results <- list()
    for(edge in colnames(net)){
      cmpnents <- unlist(strsplit(edge, split=","))
      if(length(cmpnents)==2&&cmpnents[2]!=cmpnents[1]){
        
        row <- row + 1
        print(paste0("[",format(Sys.time(),"%a %b %d %X %Y"),"] Estimating baseline LAD model and saving results for [",row,"/",n*(n-1),"] ",edge,".."))
        actual[,row] <- as.numeric(data.frame(net[rownames(net)[-c(1,2)],edge])[keep,])
        colnames(actual)[row] <- edge
        
        # Baseline model
        filename2 <- paste0("results/models/benchmark_",edge)
        if(file.exists(filename2)) load(filename2) else{
          dta <- cbind(actual[,row],mdl.data[keep,grep(edge, colnames(mdl.data))])
          colnames(dta)[1] <- "y"
          mdl <- Bykhovskaya.regression(dta, 10^5, fit.lower.bound = 0, fit.upper.bound = 100)
          save(mdl, file = filename2)
        }
        fit[["base"]][,row] <- mdl$fit
        colnames(fit[["base"]])[row] <- edge
        results[["base"]][[edge]] <- mdl$coeficients$opmod
        #plot(as.numeric(data.frame(net[rownames(net)[-c(1,2)],edge])[keep,]), type='l')
        #lines(mdl$fit, col="blue")
        
        # Other models
        for(mdl.name in names(X.cols)){
          print(paste0("[",format(Sys.time(),"%a %b %d %X %Y"),"] Estimating ",mdl.name," LAD model and saving results for [",row,"/",n*(n-1),"] ",edge,".."))
          filename2 <- paste0("results/models/bench_",mdl.name,"_",edge)
          if(file.exists(filename2)) load(filename2) else {
            dta <- cbind(actual[,row],mdl.data[keep,c(grep(edge, colnames(mdl.data)),X.cols[[mdl.name]])])
            colnames(dta)[1] <- "y"
            mdl <- Bykhovskaya.regression(dta, 10^5, fit.lower.bound = 0, fit.upper.bound = 100)
            save(mdl, file = filename2)
          }
          fit[[mdl.name]][,row] <- mdl$fit
          colnames(fit[[mdl.name]])[row] <- edge
          results[[mdl.name]][[edge]] <- mdl$coeficients$opmod
          #lines(mdl$fit, col="red")
        }
      }
    }
    save(actual, file = "C:/results/actual.bench")
    save(fit, file = "C:/results/fit.bench")
    save(results, file = filename)
  }
  
  return(list(actual=actual,fit=fit,coef.summary=results))
  
}


results.out <- function(actual,fit,results,n){
  
  TRBC.map <- list("X50"="EN",
                   "X51"="MT",
                   "X52"="ID",
                   "X53"="CS_C",
                   "X54"="CS_N",
                   "X55"="FN",
                   "X56"="HC",
                   "X57"="IT",
                   "X59"="UT",
                   "X60"="RE")
  
  # 3.1 General model R^2 and F-stat for S_plus and S_minus coefficients
  R2 <- 0
  T_obs <- length(unlist(c(actual)))
  cutoffs <- list()
  nms <- names(fit)
  for(ctf in nms){
    nR2 <- cor(unlist(c(actual)),unlist(c(fit[[ctf]])))^2
    cutoffs <- c(cutoffs,nrow(results[[ctf]][[1]]))
    print(paste0("R2 for model ",ctf," is: ", sprintf("%.3f", nR2)))
    print(paste0("Difference for model ",ctf," is: ", sprintf("%.3f", nR2-R2), " (",sprintf("%.4f", 1 - pnorm((atanh(sqrt(nR2))-atanh(sqrt(R2)))/sqrt(2/(T_obs - 3)))),")"))
    print(paste0("Adjusted R2 for model ",ctf," is: ",sprintf("%.3f", 1-(1-nR2)*(T_obs - 1)/(T_obs - nrow(results[[ctf]][[1]]) - 1))))
    R2 <- nR2
  }
  cutoffs <- as.numeric(unlist(cutoffs))
  
  # 3.2 F-statistics for each individual equation
  for(i in 2:length(nms)){
    for(j in 1:(i-1)){
      Fstat <- sd(unlist(c(actual - fit[[nms[i-j]]]))) / sd(unlist(c(actual - fit[[nms[i]]])))
      print(paste0("F-stat for model ",i," vs. model ",i-j," is: ", sprintf("%.4f",Fstat), " (", sprintf("%.4f", 1-pf(Fstat, T_obs-1, T_obs-1)), ")"))
    }
  }
  
  # 3.3 Number of significant coefficients at the 5% level in each model
  for (mdl.name in nms) {
    significant.coefs <- 0
    total.coefs <- 0
    for(edge in names(results[[mdl.name]])){
      total.coefs <- total.coefs + nrow(results[[mdl.name]][[edge]])
      significant.coefs <- significant.coefs + sum(results[[mdl.name]][[edge]]$pval<=0.05)
    }
    print(paste0("Total number of significant coefficients (p <= 0.05) in model ",mdl.name," is ",significant.coefs,", i.e. ",sprintf("%.3f",100*significant.coefs/total.coefs),"%."))
  }
  
  # 3.4 Histograms of t-statistics for each coefficient in each model
  for (mdl.name in nms) {
    coef.names <- rownames(results[[mdl.name]][[names(results[[mdl.name]])[1]]])
    edge.names <- names(results[[mdl.name]])
    tstats <- matrix(nrow = length(edge.names), ncol = length(coef.names))
    rownames(tstats) <- edge.names
    colnames(tstats) <- coef.names
    for(edge in edge.names){
      tstats[edge,] <- results[[mdl.name]][[edge]]$tstat
    }
    chart.nrow <- (ncol(tstats)-1)/3
    png(file=paste0("results/hist_tstats_",mdl.name,".png"), width=6, height=2*chart.nrow, units="in", res=600)
    par(mfrow=c(chart.nrow,3), mar = c(4, 2, 4, 2))
    for(i in 2:ncol(tstats)){
      coef.name <- colnames(tstats)[i]
      my_hist <- hist(tstats[,i], breaks = 20, plot = F)
      cls <- ifelse(my_hist$breaks < -1.96, "red", ifelse(my_hist$breaks > 1.96, "forestgreen", "gray"))
      title <- "" #title <- Tex(paste0("t-stats for ",coef.name,""))
      plot(my_hist, border=F, xlab=TeX(paste0("$\\beta_{",i-1,"}$")), ylab="", main = title, cex = 0.5, col = cls)
    }
    dev.off()
  }
  
  # 3.5 Coefficient matrices for each variables vs. each edge in the network
  for (mdl.name in nms) {
    edge.names <- names(results[[mdl.name]])
    coef.names <- rownames(results[[mdl.name]][[names(results[[mdl.name]])[1]]])
    for(cn in coef.names){
      tbl <- data.frame(matrix(nrow = data$n, ncol = data$n))
      colnames(tbl) <- rownames(tbl) <- data$peers
      for(edge in edge.names){
        c1 <- unlist(strsplit(edge,","))[1]
        c2 <- unlist(strsplit(edge,","))[2]
        #stars <- ifelse(results[[mdl.name]][[edge]][cn,"pval"]<0.01,"^{***}",ifelse(results[[mdl.name]][[edge]][cn,"pval"]<0.05,"^{**}",ifelse(results[[mdl.name]][[edge]][cn,"pval"]<0.1,"^{*}","")))
        stars <- ifelse(results[[mdl.name]][[edge]][cn,"pval"]<0.01,"***",ifelse(results[[mdl.name]][[edge]][cn,"pval"]<0.05,"**",ifelse(results[[mdl.name]][[edge]][cn,"pval"]<0.1,"*","")))
        #tbl[c1,c2] <- capture.output(cat("\\begin{tabular}{@{}c@{}}",sprintf("%.4f", results[[mdl.name]][[edge]][cn,"coef"]),stars,"\\\\(",sprintf("%.3f", results[[mdl.name]][[edge]][cn,"se"]),")\\end{tabular}", sep = ""))
        tbl[c1,c2] <- paste0(sprintf("%.4f", results[[mdl.name]][[edge]][cn,"coef"]),stars," (",sprintf("%.3f", results[[mdl.name]][[edge]][cn,"tstat"]),")")
        #tbl[c1,c2] <- results[[mdl.name]][[edge]][coef.name,"t-stat"]
        tbl[c1,c2] <- paste0(sprintf("%.4f", results[[mdl.name]][[edge]][cn,"coef"]),stars)
        
      }
      #print(xtable(tbl, type = "latex", digits = 3, caption = paste0("t-statistics for ",cn), label = paste0("tbl:",cn,"_",unlist(strsplit(coef.name,"_"))[2])), 
      #      file=paste0("C:/results/tstats.",cn,"_",unlist(strsplit(coef.name,"_"))[2],".txt"), 
      #      include.rownames=TRUE)
      print(xtable(tbl, type = "latex", caption = paste0("Coefficients for ",cn), label = paste0("tbl:",mdl.name,"_",cn)), 
            file=paste0("results/result_tbl_",mdl.name,"_",cn,".txt"), 
            include.rownames=TRUE)
    }
  }
  
}

