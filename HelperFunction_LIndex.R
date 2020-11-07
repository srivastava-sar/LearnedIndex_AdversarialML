## @knitr Polyfit_amz

func_fitPoly <- function(X,df_data){
  library(boot)
  print(length(X))
  print(length(date))
  degree=1:10
  cv.error=rep(0,10)
  for (d in degree){
    glm.fit=glm(X~poly(date,d),data=df_data)
    cv.error[d]=cv.glm(df_data,glm.fit)$delta[1]
  }
  print("cv.error")
  print(cv.error)
  
  plot(degree,cv.error,type='b')
}

## @knitr Polyfit_uid

func_fitPoly_uid <- function(X,df_data){
  library(boot)
  print(length(X))
  print(length(date))
  degree=1:10
  cv.error=rep(0,10)
  for (d in degree){
    glm.fit=glm(X~poly(uid,d),data=df_data)
    cv.error[d]=cv.glm(df_data,glm.fit)$delta[1]
  }
  print("cv.error")
  print(cv.error)
  
  plot(degree,cv.error,type='b')
}



## @knitr Polyfit_Visual

func_PolyVis <- function(X,deg,df_keys, df_data){
  QuadFunc=glm(X~poly(df_keys,deg),data=df_data)
  plot(X~df_keys,df_data)
  points(df_keys,fitted(QuadFunc),col='red',pch=2)
  return(QuadFunc)
}


## @knitr key_diff_visual

func_keyDiff <- function(X,df_keys){
  keydiff <- data.frame(diff(as.matrix(df_keys)))
  keyDim <- dim(keydiff)[1]
  plot(keydiff[,1]~X[1:keyDim],type="l", xlab="index", ylab="KeyDiff", main="Plot of Adj. Keys diff vs index") #Adhacent key differences vs index (x-axis)
  #dev.copy(png,'keydiff.png')
  #dev.off()
  return(keydiff)
}

## @knitr GaussFit

func_FitNorm <- function(keyDiff){
  library(fitdistrplus)
  FIT <- fitdist(unlist(keyDiff), "norm")    ## note: it is "norm" not "normal"
  class(FIT)
  print(FIT)
  plot(FIT)    ## use method `plot.fitdist`
  return(FIT)
}



## @knitr CUSUM_Impl2

func_CUSUM2 <- function(FIT,T_factor,datediff,dateDim,y_lim){
  Cp=rep(0,length(datediff[,1]))
  Cm=rep(0,length(datediff[,1]))
  
  mean=FIT$estimate[1]
  K = mean/2
  #Thresh = 5 * FIT$estimate[2]
  
  Thresh=T_factor*FIT$estimate[[2]]
  flagP=FALSE
  flagN=FALSE
  flagEndP=FALSE
  flagEndN=FALSE
  exitpointP=0
  exitpointN=0
  outcomeExitsP<-0
  outcomeExitsN<-0
  
  for (ii in 1:dateDim){
    if (ii == 1){
      Cp[ii]=0
      Cm[ii]=0
    }else{
      Cp[ii]=max(0,datediff[ii,1]-(mean+K)+Cp[ii-1])
      Cm[ii]=max(0,(mean-K)-datediff[ii,1]+Cm[ii-1])
    }
    if (isTRUE(Cp[ii] >= Thresh)){
      if(flagP==FALSE){
        #print(ii)
        flagP=TRUE
        outcomeExitsP <- c(outcomeExitsP, ii)
      }else{
        exitpointP=ii
        flagEndP=FALSE
      }
    }else{
      flagP=FALSE
      if(flagEndP==FALSE){
        #print(exitpointP)
        outcomeExitsP <- c(outcomeExitsP, exitpointP)
        flagEndP=TRUE
      }
    }
    if (isTRUE(Cm[ii] >= Thresh)){
      #print(Cm[ii])
      if(flagN==FALSE){
        #print(ii)
        flagN=TRUE
        outcomeExitsN <- c(outcomeExitsN, ii)
      }else{
        exitpointN=ii
        flagEndN=FALSE
      }
    }else{
      flagN=FALSE
      if(flagEndN==FALSE){
        #print("Negative")
        #print(exitpointN)
        outcomeExitsN <- c(outcomeExitsN, exitpointN)
        flagEndN=TRUE
      }
    }
  }
  #print(outcomeExitsP)
  #length(outcomeExitsP)
  outcome <- cbind(datediff[,1],Cp,Cm)
  
  #print(outcomeExitsN)
  
  head(outcome)
  actualOutcomeP <-0
  actualOutcomeP <- outcomeExitsP[3:length(outcomeExitsP)]
  actualOutcomeP <- cbind(actualOutcomeP,actualOutcomeP)
  
  actualOutcomeN <-0
  actualOutcomeN <- outcomeExitsN[3:length(outcomeExitsN)]
  actualOutcomeN <- cbind(actualOutcomeN,actualOutcomeN)
  
  if(length(outcomeExitsN) > 2 & length(outcomeExitsP) > 2){
    outcomeExits <- append(outcomeExitsP[3:length(outcomeExitsP)],outcomeExitsN[3:length(outcomeExitsN)])
  } else if (length(outcomeExitsN)<=2){
    outcomeExits <- outcomeExitsP[3:length(outcomeExitsP)]
  }else {
    outcomeExits <- outcomeExitsN[3:length(outcomeExitsN)]
  }  
  #Visualizig on CUSUM plot
  par(mfrow = c(2,1))
  plot(outcome[,1]~X[1:dateDim],type="l",col="black",ylim=c(-20,y_lim),main="CUSUM plot based on difference in adjacent keys",xlab="Index",ylab="Adj. Key diff.",sub="Pink=KeyDiff, Green=Cp, Red=Cm, Brown=Start/End")
  points(X[1:dateDim],outcome[,2], col="green",cex=2,pch=20,type = "l")
  points(X[1:dateDim],outcome[,3], col="red",cex=2,pch=20,type = "l")
  if(length(actualOutcomeP)>4){
    points(actualOutcomeP,outcome[actualOutcomeP,2], col="lightgreen",cex=2,pch=20) #Predicted residual positions near the actual based on CUSUM
  }
  if(length(actualOutcomeN)>4){
    points(actualOutcomeN,outcome[actualOutcomeN,3], col="pink",cex=2,pch=20) #Predicted residual positions near the actual based on CUSUM
  }
  abline(h=Thresh)
  return(outcomeExits)
}




## @knitr PlotRes

func_plotRes <- function(Quad2,maxQ,ultraQ){
  plot(Quad2$residuals^2,col="orange",type="line",main = "Residual square plot based on model fitting",xlab = "Index", ylab = "Residual Sq.",sub = "Orange=ResidualSq., Blue=Top Residuals, Black=MaxResidual") #Actual residuals
  points(maxQ,Quad2$residuals[maxQ-1]^2, col="blue",cex=2,pch=20) #Actual maximum residual based on maxQ
  points(ultraMax,Quad2$residuals[ultraMax]^2, col="black",cex=2,pch=20) #UltraMax
  # dev.copy(png,'myplot.png')
  # dev.off()
}



## @knitr Evaluation1

func_eval1 <- function(lastN,keyDim,Quad2,maxQ,outcomeExits,printFlag=TRUE){
  # Residual ordering
  s_SqRes <- data.frame(sqres=Quad2$residuals^2)
  s_SqRes$rank <-  order(order(s_SqRes$sqres, decreasing=TRUE))
  sort(s_SqRes$rank, index.return=TRUE)$ix[1:30]
  
  countSet <- 0
  resultSet <- 0
  relevanceScore <- 0
  maxRekevance <- 0
  for (j in 1:length(maxQ)){
    for (i in 1:length(outcomeExits)){
      if ((maxQ[[j]] > outcomeExits[i]-lastN) && (maxQ[[j]] <= outcomeExits[i] + lastN)){
        countSet <- countSet + 1
        resultSet <- c(resultSet,maxQ[[j]])
        relevanceScore <- relevanceScore + (1/s_SqRes[maxQ[[j]],]$rank)
        #print(s_SqRes[maxQ[[j]],]$rank)
        break
      }
    }
    maxRekevance <- maxRekevance + (1/s_SqRes[maxQ[[j]],]$rank)
  }
  if(printFlag){
    plot(Quad2$residuals^2,col="orange",type="line",main = "Residual square plot with ponts detected using CUSUM",xlab = "Index", ylab = "Residual Sq.") #Actual residuals
    points(resultSet[2:length(resultSet)],Quad2$residuals[resultSet]^2, col="green",cex=2,pch=20) #Actual maximum residual based on maxQ
    # dev.copy(png,'myplot.png')
    # dev.off()
    print("countSet")
    print(countSet)    
  }
  EvalMeasure = relevanceScore / maxRekevance
  return(EvalMeasure)
}


## @knitr Evaluation2

func_eval2 <- function(outcomeExits,dateDim,lastN){
  #No of evaluations:
  NoOfEval <- 2*lastN*length(outcomeExits)/dateDim
  return(NoOfEval)
}



## @knitr CUSUM_Impl3

func_CUSUM3 <- function(FIT,T_factor,datediff,dateDim,y_lim,plotflag=TRUE){
  Cp=rep(0,length(datediff[,1]))
  Cm=rep(0,length(datediff[,1]))

  mean=FIT$estimate[1]
  K = mean/2
  #Thresh = 5 * FIT$estimate[2]
  
  Thresh=T_factor*FIT$estimate[[2]]
  flagP=FALSE
  flagN=FALSE
  flagEndP=FALSE
  flagEndN=FALSE
  exitpointP=0
  exitpointN=0
  outcomeExitsP<-c()
  outcomeExitsN<-c()
  
  for (ii in 1:dateDim){
    if (ii == 1){
      Cp[ii]=0
      Cm[ii]=0
    }else{
      Cp[ii]=max(0,datediff[ii,1]-(mean+K)+Cp[ii-1])
      Cm[ii]=max(0,(mean-K)-datediff[ii,1]+Cm[ii-1])
    }
    if (isTRUE(Cp[ii] >= Thresh)){
      if(flagP==FALSE){
        #print(ii)
        flagP=TRUE
        outcomeExitsP <- c(outcomeExitsP, ii)
      }else{
        exitpointP=ii
        flagEndP=FALSE
      }
    }else{
      flagP=FALSE
      if(flagEndP==FALSE){
        #print(exitpointP)
        #outcomeExitsP <- c(outcomeExitsP, exitpointP)
        flagEndP=TRUE
      }
    }
    if (isTRUE(Cm[ii] >= Thresh)){
      #print(Cm[ii])
      if(flagN==FALSE){
        #print(ii)
        flagN=TRUE
        outcomeExitsN <- c(outcomeExitsN, ii)
      }else{
        exitpointN=ii
        flagEndN=FALSE
      }
    }else{
      flagN=FALSE
      if(flagEndN==FALSE){
        #print("Negative")
        #print(exitpointN)
        #outcomeExitsN <- c(outcomeExitsN, exitpointN)
        flagEndN=TRUE
      }
    }
  }
  #print("positive")
  #print(outcomeExitsP)
  #print(length(outcomeExitsP))
  outcome <- cbind(datediff[,1],Cp,Cm)
  
  #print(outcomeExitsN)

  if(length(outcomeExitsN) > 0 & length(outcomeExitsP) > 0){
    outcomeExits <- list(P=outcomeExitsP,N=outcomeExitsN)
  } else if (length(outcomeExitsP) > 0){
    outcomeExits <- list(P=outcomeExitsP,N=c())
  }else {
    outcomeExits <- list(P=c(),N=outcomeExitsN)
  }  
  
  if(plotflag){
    #print("outcomeExits")
    #print(outcomeExits)
    # Visualizig on CUSUM plot
    par(mfrow = c(1,1))
    plot(outcome[,1]~X[1:dateDim],type="l",col="black",ylim=c(-10,y_lim),main="Reverse CUSUM plot based on difference in adjacent keys",xlab="Index",ylab="Adj. Key diff.")
    points(X[1:dateDim],outcome[,2], col="green",cex=2,pch=20,type = "l")
    points(X[1:dateDim],outcome[,3], col="red",cex=2,pch=20,type = "l")
    if(length(outcomeExits$P)>0){
      points(outcomeExits$P,outcome[outcomeExits$P,2], col="lightgreen",cex=2,pch=20) #Predicted residual positions near the actual based on CUSUM
    }
    if(length(outcomeExits$N)>0){
      points(outcomeExits$N,outcome[outcomeExits$N,3], col="pink",cex=2,pch=20) #Predicted residual positions near the actual based on CUSUM
    }
    legend(0, 2900000, legend=c(expression("Z+"), expression("Z-"),"KeyDiff","Threshold"),
           col=c("green", "red","black","brown"), lty=c(1,1,1,2), cex=0.8)
    abline(h=Thresh,col="brown",type="l",lty=2) 
    dev.copy(png,'2CUSUMa.png')
    dev.off()
  }
  return(outcomeExits)
}


## @knitr CUSUM_Wrapper

func_CUSUM_wrapper <- function(FIT,T_factor,keydiff,keyDim,y_lim,plotFlag=TRUE){
  
  keydiff_rev <- data.frame(Rev(keyDiff,1))
  outcomeExits = func_CUSUM3(FIT,T_factor,keyDiff,keyDim,y_lim,FALSE)
  outcomeExits_rev = func_CUSUM3(FIT,T_factor, keydiff_rev,keyDim,y_lim,plotFlag)
  
  outcomeExits$PDiff <- round((outcomeExits$P + Rev(keyDim-outcomeExits_rev$P,1))/2,0)
  outcomeExits$NDiff <- round((outcomeExits$N + Rev(keyDim-outcomeExits_rev$N,1))/2,0)
  
  actualExits <- c(outcomeExits$PDiff,outcomeExits$NDiff)

  return(actualExits)
}



