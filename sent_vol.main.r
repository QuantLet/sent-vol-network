## 0. Load libraries
source("sent_vol.scripts.r")
load.libraries()
summ.stat <- function(series,label=NULL){
  ss <- c(
    count=round(sum(is.finite(series))),avg=round(mean(series,na.rm=T),3),sd=round(sd(series,na.rm=T),3),
    skew=round(skewness(series,na.rm=T),3),kurt=round(kurtosis(series,na.rm=T),3),
    min=round(min(series,na.rm=T),3),q1=round(quantile(series,0.25,na.rm=T),3),q2=round(quantile(series,0.5,na.rm=T),3),
    q3=round(quantile(series,0.75,na.rm=T),3),max=round(max(series,na.rm=T),3)
  )
  if(is.null(label)) return(ss)
  rn <- names(ss)
  ss <- data.frame(ss)
  colnames(ss) <- label
  rownames(ss) <- rn
  return(ss)
}

## 1. Load & process data
# 1.1 Load Financial network data
net <- load.net()
# 1.2 Load Market return data
mkt <- load.mkt(net)
require(moments)
desc.stat <- summ.stat(mkt$MKT)
# 1.3 Load Market sentiment data
sentiment <- load.sentiment(net)
par(mar=c(2, 4, 2, 4) + 0.1)
plot(as.Date(rownames(sentiment)),sentiment[,1], axes=FALSE, type='l',main="",xlab="",ylab="",col="black",ylim=c(-7,7))
box()
#axis(1,pretty(range(as.Date(rownames(sentiment)))))
axis(1, as.Date(rownames(sentiment))[seq(1,nrow(sentiment),length.out=9)], 
     format(as.Date(rownames(sentiment)[seq(1,nrow(sentiment),length.out=9)]), "%b %Y"), cex.axis = .7)
mtext("Investor Sentiment",side=2, line=2)
axis(2, ylim=c(0,1),col="black",las=1)
abline(h=0,lty="dotted")
par(new=TRUE)
plot(as.Date(rownames(sentiment)),mkt$MKT, axes=FALSE,type='l',main="",xlab="",ylab="",col=rgb(1,0,0,0.2),ylim=c(-11,11))
axis(4, ylim=c(0,1),col="black",las=1)
mtext("Market Return",side=4, line=2)
desc.stat <- summ.stat(sentiment$TOT)
# 1.4 Construct SAMPLE data.frame
data <- construct.data(net,mkt,sentiment)
# 1.5 SUMMARY STATISTICS
vol <- load.vol()
summary(vol)
ss <- summ.stat(vol[,1],colnames(vol)[1])
if(ncol(vol)>1) for(j in 2:ncol(vol)) ss <- cbind(ss,summ.stat(vol[,j],colnames(vol)[j]))
xtable(ss, type = "latex", digits = 3)

## 2. Run Bykhovskaya regression models on the network
estimation.results <- run.Bykhovskaya.regression(
  net=net,
  mdl.data=data$mdl.data,
  keep=data$keep,
  n=data$n,
  X.cols = list(market = (data$cutoffs[length(data$cutoffs)-1]+1):(data$cutoffs[length(data$cutoffs)-1]+3),
                sentiment = (data$cutoffs[length(data$cutoffs)-1]+1):data$cutoffs[length(data$cutoffs)]))

## 3. Analyze and output results 
results.out(
  actual=estimation.results$actual, 
  fit = estimation.results$fit, 
  results = estimation.results$coef.summary)





{
  
  # 2.7 Estimate baseline Bykhovskaya (2023) models and save results
  
  # 2.8
  results.out(actual,fit,results,order.reverse=T)
  
  # 2.9 Estimate extended models and save results
  filename <- "C:/results/coefficients"
  if(file.exists(filename)){
    load("C:/results/actual")
    load("C:/results/fit")
    load(filename)
  }else{
    row <- 0
    actual <- data.frame(matrix(nrow = sum(keep), ncol = n * (n-1)))
    fit <- list()
    for(ctf in cutoffs) fit[[as.character(ctf)]] <- data.frame(matrix(nrow = sum(keep), ncol = n * (n-1)))
    results <- list()
    for(edge in colnames(net)){
      cmpnents <- unlist(strsplit(edge, split=","))
      if(length(cmpnents)==2&&cmpnents[2]!=cmpnents[1]){
        row <- row + 1
        print(paste0("[",format(Sys.time(),"%a %b %d %X %Y"),"] Estimating LAD models and saving results for [",row,"/",n*(n-1),"] ",edge,".."))
        actual[,row] <- as.numeric(data.frame(net[rownames(net)[-c(1,2)],edge])[keep,])
        colnames(actual)[row] <- edge
        for(ctf in cutoffs){
          print(paste0("[",format(Sys.time(),"%a %b %d %X %Y"),"] Estimating and saving model with ",ctf," explanatory variables.."))
          # Constructing data set
          dta <- cbind(actual[,row],mdl.data[keep,1:ctf])
          colnames(dta)[1] <- "y"
          # Estimating and saving model
          mdl <- lad(y~., data = dta)
          filename2 <- paste0("C:/results/models/",edge,"_",ctf)
          save(mdl, file = filename2)
          # Saving model fit
          fit[[as.character(ctf)]][,row] <- mdl$fitted.values
          colnames(fit[[as.character(ctf)]])[row] <- edge
          # Saving parameter estimates and associated statistics
          se <- sqrt(diag(vcov(mdl)))
          tstats <- as.numeric(mdl$coefficients) / se
          pvals <- 2 * pt(abs(tstats), df = nrow(dta) - nrow(vcov(mdl)) - 2, lower.tail = FALSE)
          res <- data.frame(as.numeric(mdl$coefficients), se, tstats, pvals)
          colnames(res) <- c("coef","s.e.","t-stat","p-value")
          rownames(res) <- names(mdl$coefficients)
          results[[as.character(ctf)]][[edge]] <- res
        }
      }
    }
    
    save(actual, file = "C:/results/actual")
    save(fit, file = "C:/results/fit")
    save(results, file = filename)
    rm(cmpnents,edge,filename,keep,pvals,R2,row,se,tstats,mdl,res,dta)
  }
  
  # 2.10 
  results.out(actual,fit,results)
}

## 3. Organize and output results

# 3.4 Tables of S_plus and S_minus coefficients
