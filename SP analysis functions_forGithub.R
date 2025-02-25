# SP analysis functions master

######### FUNCTIONS ##########
addInMean=function(a2,newCol,startCol,endCol,colName){
  meanValList=c()
  for(i in 1:nrow(a2)){
    a2b=a2[i,startCol:endCol]
    meanVal=mean(na.omit(as.numeric(a2b)))
    meanValList=c(meanValList,meanVal)
  }
  if(newCol<=ncol(a2)){
    a2c=cbind(a2[,1:(newCol-1)],temp=meanValList,a2[,newCol:ncol(a2)])
    names(a2c)[newCol]=colName
  }
  else{
    a2c=cbind(a2,temp=meanValList)
    names(a2c)[ncol(a2c)]=colName
  }
  return(a2c)
}
addInPeak=function(a2,newCol,startCol,endCol,colName,direction='min'){
  minValList=c()
  for(i in 1:nrow(a2)){
    a2b=a2[i,startCol:endCol]
    if(direction=='min'){
      minVal=min(as.numeric(a2b),na.rm=TRUE)
    } else{
      minVal=max(as.numeric(a2b),na.rm=TRUE)
    }
    minValList=c(minValList,minVal)
  }
  a2c=cbind(a2[,1:(newCol-1)],temp=minValList,a2[,newCol:ncol(a2)])
  names(a2c)[newCol]=colName
  return(a2c)
}
addInPeakDensity=function(a2,capCol,peakCol,peakDensityCol,peakDensityName='PeakDensity'){
  peakDensity=a2[,peakCol]/a2[,capCol]
  colSize=ncol(a2)
  if(colSize+1==peakDensityCol){
    a2b=cbind(a2,pd=peakDensity)
  } else{
    a2b=cbind(a2[,1:(peakDensityCol-1)],pd=peakDensity,a2[,peakDensityCol:ncol(a2)])
  }
  names(a2b)[peakDensityCol]=peakDensityName
  return(a2b)
}
addInPeakDensityNorm=function(a4,wtMean,peakDensityCol,colNum,name='PeakDensityNorm'){
  PeakDensityNorm=round(100*as.numeric(a4[,peakDensityCol])/wtMean,2)
  colSize=ncol(a4)
  if(colSize+1==colNum){
    a5=cbind(a4,PeakDensityNorm)
  }
  else{
    a5=cbind(a4[,1:(colNum-1)],PeakDensityNorm,a4[,colNum:(ncol(a4))])
  }
  names(a5)[colNum]=name
  return(a5)
}
addSpaces=function(mutFrameForFig){
  last=mutFrameForFig[1,'MutDescription']
  firsttime=TRUE
  blank=mutFrameForFig[1,]
  blank[1,1]='*'
  blank[1,2]='*'
  blank[1,3]=0
  blank[1,4]=0
  blank[1,5]=0
  for(i in 1:nrow(mutFrameForFig)){
    current=mutFrameForFig[i,'MutDescription']
    if(firsttime){
      mutFrameForFig2=mutFrameForFig[1,]
      firsttime=FALSE
    }
    else{
      if(current==last){
        mutFrameForFig2=rbind(mutFrameForFig2,mutFrameForFig[i,])
      }
      else{
        mutFrameForFig2=rbind(mutFrameForFig2,blank,mutFrameForFig[i,])
      }
    }
    last=current
  }
  return(mutFrameForFig2)
}
bulk90adjust = function(i,inactSummary5){
  VhalfInact_WT=inactSummary5[inactSummary5$Mutation=='WT','VHalfInactMean']
  for(row in 1:nrow(i)){
    peakDensityNorm=i[row,'PeakDensityNorm']
    currMutation=i[row,'Mutation']
    VhalfInact_mut=inactSummary5[inactSummary5$Mutation==currMutation,'VHalfInactMean']
    if(length(VhalfInact_mut)==0){
      VhalfInact_mut=VhalfInact_WT
    }
    deltaBoltzmannVal=deltaBoltzmann(k=15.8,VhalfInact_WT,VhalfInact_mut)
    peakDensityNormBulk90=peakDensityNorm*deltaBoltzmannVal
    i[row,'PeakDensityNormBulk90']=peakDensityNormBulk90
  }
  return(i)
}
calcActivation=function(a3,minPeakForAct,errorCutoff,IV_voltages,reversalPotential,firstCol,peakCol,outfile,replaceWithNorm=FALSE){
  pdf(outfile)
  par(mfrow=c(4,4))
  a3b=cbind(a3[,1:(firstCol-1)],VHalfAct=rep(NA,nrow(a3)),VHalfActError=rep(NA,nrow(a3)),VHalfActInclude=rep(NA,nrow(a3)),a3[firstCol:ncol(a3)])
  sweepList=1:length(IV_voltages)
  for(i in 1:nrow(a3)){
    type=a3[i,'Type']
    well=a3[i,'Well']
    allAdjustedPeakVals=as.numeric(a3[i,firstCol:(firstCol+length(IV_voltages)-1)])
    allAdjustedPeakVals2=as.numeric(a3[i,firstCol:(firstCol+29-1)])
    peak=a3[i,peakCol]
    print(well)
    if(!type=='good' | peak>minPeakForAct){
      a3b[i,'VHalfAct']=NA
      a3b[i,'VHalfActError']=NA
      a3b[i,'VHalfActInclude']=FALSE
    } 
    else{
      allExpected=IV_voltages-reversalPotential
      allExpected2=seq(-80,60,5)-reversalPotential
      normalized1=allAdjustedPeakVals/allExpected
      normalized2=allAdjustedPeakVals2/allExpected2
      
      l=length(normalized1)
      toAdjustBy=mean(normalized1[(l-4):l])
      activationratio=normalized1/toAdjustBy
      activationratio2=normalized2/toAdjustBy
      #print(activationratio2)
      if(replaceWithNorm){
        for(j in 1:29){
          a3b[i,j+firstCol+2]=activationratio2[j]
        }
      }
      plot(IV_voltages,activationratio,ylab='Activation fraction',xlab='Voltage (mV)',main=well)
      abline(h=0,lty=2)
      abline(h=1,lty=2)
      #fits ratio to boltzmann equation
      outputB=NA
      errorSD=NA
      include=FALSE
      try({
        fitmodel <- nls(activationratio~(1/(1+exp(((a-IV_voltages)/b)))), start=list(a=-50,b=1))
        result=coef(fitmodel)
        outputA=list(a=result[1],b=activationratio, c=predict(fitmodel,data.frame(voltages=IV_voltages)))
        outputB=result[1]
        bestFit=predict(fitmodel,data.frame(voltages=IV_voltages))
        errorSD=sd(bestFit-activationratio)
        if(errorSD<errorCutoff){
          include=TRUE
        }
      })
      if(!is.na(outputB)){
        lines(IV_voltages,predict(fitmodel))
        text(0,0.5,round(outputB,1))
        textCol='black'
        if(errorSD>errorCutoff){
          textCol='red'
        }
        text(0,0.2,round(errorSD,2),col=textCol)
      }
      a3b[i,'VHalfAct']=round(outputB,1)
      a3b[i,'VHalfActError']=round(errorSD,3)
      a3b[i,'VHalfActInclude']=include
    }
  }
  dev.off()
  return(a3b)
}
calcAndPlotPeakDensityIV=function(a3,mutList,startCol=20,endCol=48){
  result=data.frame()
  count=0
  for(mut in mutList){
    count=count+1
    result[count,'Mutation']=mut
    result[count,'MutDescription']=a3[which(a3$Mutation==mut)[1],'MutDescription']
    
    # Only normals
    print(mut)
    a4=a3[a3$Mutation==mut,]
    data=a4[,startCol:endCol]
    data2=data
    for(i in 1:nrow(data)){
      cap=a4[i,'Capacitance']
      for(j in 1:ncol(data)){
        data2[i,j]=data[i,j]/cap
      }
    }
    dataMean=sapply(data,mean,na.rm=TRUE)
    data2Mean=sapply(data2,mean,na.rm=TRUE)
    library('plotrix')
    dataSE=sapply(data,std.error,na.rm=TRUE)
    data2SE=sapply(data2,std.error,na.rm=TRUE)
    peakDensity=min(data2Mean)
    numSteps=length(data2Mean)
    
    result[count,'peakDensity']=peakDensity
    result[count,'numCells']=nrow(data2)
    for(jj in 1:numSteps){
      result[count,4+jj]=data2Mean[jj]
    }
    for(jj in 1:numSteps){
      result[count,4+numSteps+jj]=data2SE[jj]
    }
  }
  names(result)=c('Mutation','MutDescription','peakDensity','numCells',names(data2Mean),names(data2SE))
  result2=result[order(result$peakDensity),]
  result2$peakDensityNorm_all=result2$peakDensity/result2[result2$Mutation=='WT','peakDensity']
  result2=cbind(result2[,1:3],peakDensityNorm_all=result2[,ncol(result2)],result2[,5:ncol(result2)-1])
  return(result2)
}
calcCOV=function(rfi,colList){
  #coefficient of variation
  covList=c()
  for(row in 1:nrow(rfi)){
    vals=as.numeric(rfi[row,colList])
    cov=abs(sd(vals)/mean(vals))
    covList=c(covList,cov)
  }
  return(covList)
}
calcFancyDelta=function(i3_good,var,varDelta,ratio=FALSE){
  newCol=which(names(i3_good)==var)[1]+1
  runList=i3_good$Plate
  experimentList=remove_ending_letter_vector(runList)
  i3_good$Experiment=experimentList
  experiments=unique(experimentList)
  firsttime=TRUE
  for(experiment in experiments){
    i3_exp=i3_good[i3_good$Experiment==experiment,]
    exp_wt=i3_exp[i3_exp$Mutation=='WT',]
    wt_exp_mean=mean(exp_wt[,var],na.rm=TRUE)
    if(!ratio){
      i3_exp$delta=i3_exp[,var]-wt_exp_mean
    } else{
      i3_exp$delta=i3_exp[,var]/wt_exp_mean
    }
    if(firsttime){
      i4=i3_exp
      firsttime=FALSE
    } else{
      i4=rbind(i4,i3_exp)
    }
  }
  numCols=ncol(i4)
  names(i4)[numCols]=varDelta
  i4=i4[,c(1:(newCol-1),numCols,newCol:(numCols-1))]
  return(i4)
}
calcInactivation=function(inact4,inactCutoff,errorCutoff,Inact_voltages,firstCol,peakCol,outpdf,replaceWithNorm=FALSE,includeK=FALSE){
  pdf(outpdf)
  par(mfrow=c(4,4))
  inact4b=cbind(inact4[,1:(firstCol-1)],VHalfInact=rep(NA,nrow(inact4)),VHalfInactError=rep(NA,nrow(inact4)),VHalfInactInclude=rep(NA,nrow(inact4)),inact4[,firstCol:ncol(inact4)])
  for(i in 1:nrow(inact4)){
    type=inact4[i,'Type']
    well=inact4[i,'Well']
    if(as.numeric(inact4[i,peakCol])>inactCutoff | !type=='good'){
      inact4b[i,'VHalfInact']=NA
      inact4b[i,'VHalfInactError']=NA
      inact4b[i,'VHalfInactInclude']=FALSE
    } 
    else{
      print(well)
      allAdjustedPeakVals=as.numeric(inact4[i,firstCol:(firstCol+length(Inact_voltages)-1)])
      inactivationratio=allAdjustedPeakVals/min(allAdjustedPeakVals)
      if(replaceWithNorm){
        for(j in 1:length(Inact_voltages)){
          inact4b[i,j+firstCol+2]=inactivationratio[j]
        }
      }
      #fits ratio to boltzmann equation
      outputB=NA
      errorSD=NA
      includeInact=FALSE
      try({
        fitmodel <- nls(inactivationratio~(1/(1+exp(((Inact_voltages-a)/b)))), start=list(a=-100,b=1))
        result=coef(fitmodel)
        outputA=list(a=result[1],b=inactivationratio, c=predict(fitmodel,data.frame(voltages=Inact_voltages)))
        outputB=result[1]
        outputC=result[2]
        bestFit=predict(fitmodel,data.frame(voltages=Inact_voltages))
        errors=bestFit-inactivationratio
        errorSD=sd(errors)
        if(errorSD<errorCutoff){
          includeInact=TRUE
        }
      })
      plot(Inact_voltages,inactivationratio,ylab='Activation fraction',xlab='Voltage (mV)',main=well)
      abline(h=0,lty=2)
      abline(h=1,lty=2)
      if(!is.na(outputB)){
        lines(Inact_voltages,predict(fitmodel))
        text(-20,0.8,round(outputB,1))
        textCol='black'
        if(errorSD>errorCutoff){
          textCol='red'
        }
        text(-20,0.5,round(errorSD,2),col=textCol)
      }
      inact4b[i,'VHalfInact']=round(outputB,1)
      inact4b[i,'VHalfInactError']=round(errorSD,3)
      inact4b[i,'VHalfInactInclude']=includeInact
      if(includeK){
        inact4b[i,'VHalfInactK']=round(outputC,2)
      }
    }
  }
  dev.off()
  return(inact4b)
}
calcInactivationTime=function(inacttime,folder,filePrefix,fileSuffix,errorCutoff,pdfname,makeFullPlot=FALSE,stimulusIncluded=FALSE,xlim=c(9,20),addTwo=FALSE){
  pdf(pdfname)
  par(mfrow=c(4,4))
  inacttime2=inacttime
  inacttime2$inacttimeTau=NA
  inacttime2$inacttimeError=NA
  inacttime2$inacttimeInclude=NA
  count=0
  for(well in inacttime$Well){
    count=count+1
    inacttime2[count,'inacttimeInclude']=FALSE
    if(!inacttime[inacttime$Well==well,'PeakInclude']){
      next()
    }
    #well='E01'
    print(well)
    traceFile=paste0(folder,filePrefix,well,fileSuffix)
    trace=read.csv(traceFile,stringsAsFactors=FALSE,skip=3,header=FALSE)
    #print("C")
    if(stimulusIncluded){
      if(!addTwo){
        trace=trace[,c(1,2,4)]
      }
      else{
        trace=trace[,c(1,2,6)]
      }
    }
    names(trace)=c('Count','Time','Current')
    peakVal=min(trace$Current)
    peakCount=which(trace$Current==peakVal)[1]
    peakTime=trace$Time[peakCount]
    trace2=trace[trace$Count>=peakCount,]
    trace$CurrentNorm=trace$Current/peakVal
    trace2$CurrentNorm=trace2$Current/peakVal
    trace$Time2=trace$Time
    trace2$Time2=trace2$Time-peakTime
    time=trace2$Time2/1000
    time_full=trace$Time2/1000
    current=trace2$CurrentNorm
    current_full=trace$CurrentNorm
    if(max(time)>25){
      time=time[time<=25]
      #time_full=time_full[time_full<=25+peakTime/1000]
      current=current[1:length(time)]
      #current_full=current_full[1:length(current_full)]
    }
    plot(time,current,xlab='Time (ms)',ylab='Current (% of peak)',main=well)
    a=NA
    b=NA
    halflife=NA
    timeconstant=NA
    errorSD=NA
    try({
      fitmodel <- nls(current~((1-b)^time)+a, start=list(a=0.1,b=.1))
      result=coef(fitmodel)
      predicted=predict(fitmodel,data.frame(time=time))
      lines(predicted~time,col='red',lwd=2)
      a=result[1]
      b=result[2]
      halflife=log10(0.5-a)/log10(1-b)
      timeconstant=log10(1/exp(1)-a)/log10(1-b)
      abline(v=timeconstant)
      errorSD=sd(predicted-current)
      if(errorSD<errorCutoff){
        include=TRUE
        text(x=15,y=0.6,round(errorSD,2))
        text(x=15,y=0.8,round(timeconstant,2))
        inacttime2[count,'inacttimeInclude']=TRUE
      } else{
        text(x=15,y=0.8,round(timeconstant,2),col='red')
        text(x=15,y=0.6,round(errorSD,2),col='red')
      }
      if(makeFullPlot){
        plot(trace$Time/1000,trace$Current*10^12,type='n',xlim=xlim)
        lines(trace$Time/1000,trace$Current*10^12)
        #plot(time_full,-1*current_full,xlab='Time (ms)',ylab='Current (% of peak)',main=well)
        #lines(predicted~time,col='red',lwd=2)
      }
    })
    inacttime2[count,'inacttimeTau']=timeconstant
    inacttime2[count,'inacttimeError']=errorSD
  }
  par(mfrow=c(1,1))
  dev.off()
  return(inacttime2)
}
calcRFI=function(rfi,firstACol,Brange,pdfname,errorCutoff,RFI_times){
  library(nlsr)
  pdf(pdfname)
  par(mfrow=c(4,4))
  rfi2=cbind(rfi[,1:(firstACol-1)],rfiTau=NA,rfiError=NA,rfi50=NA,rfiInclude=FALSE,rfi[firstACol:ncol(rfi)])
  print(names(rfi2))
  for(row in 1:nrow(rfi)){
    well=rfi[row,'Well']
    include=rfi2[row,'goodRFI']
    if(include){
      print(well)
      peakValues=as.numeric(rfi[row,Brange])
      # Single exponential
      success=FALSE
      tryCatch({
        fitmodel <- nls(peakValues~(-1*exp(-(RFI_times/a))+1), start=list(a=1)) #single exponential
        success=TRUE
      }, error=function(e){})
      
      # Double exponential, couldn't get to work
      #fitmodel=nlxb(peakValues~(-x*exp(-RFI_times/a)-(1-x)*exp(-RFI_times2/b)+1), start=c(a=10,b=50,x=0.9))
      #success=FALSE
      #tryCatch({
      #  fitmodel <- nls(peakValues~(-x*exp(-RFI_times/a)-(1-x)*exp(-RFI_times2/b)+c), start=list(a=10,b=50,c=1,x=0.9)) #double exponential
      #  success=TRUE
      #}, error=function(e){})
      print(success)
      plot(RFI_times,peakValues,xlab='Recovery time (ms)',ylab='Fraction of peak current',main=well,ylim=c(0,1),log='x')
      abline(h=0,lty=2)
      abline(h=1,lty=2)
      if(success){
        result=coef(fitmodel)
        predictedPeak=predict(fitmodel,data.frame(RFI_times=RFI_times))
        lines(RFI_times,predictedPeak)
        text(90,0.4,paste(round(result[[1]],1)))
        errorSD=sd(predictedPeak-peakValues)
        rfi2[row,'rfiError']=errorSD
        textColor='red'
        rfi2[row,'rfiTau']=result[[1]]
        if(errorSD<errorCutoff){
          rfi2[row,'rfiInclude']=TRUE
          textColor='black'
        }
        text(90,0.15,round(errorSD[[1]],2),col=textColor)
        # Calculate time to 50% recovery
        RFI_times_fake2=seq(1,200,by=0.1)
        predictedPeak2=predict(fitmodel,data.frame(RFI_times=RFI_times_fake2))
        minPos=which(abs(predictedPeak2-0.5)==min(abs(predictedPeak2-0.5)))[1]
        rfi50=RFI_times_fake2[minPos]
        if(errorSD<errorCutoff){
          rfi2[row,'rfi50']=rfi50
        }
      } else{
        print('FAIL')
        result=list(NA,NA,NA,NA)
      }
      #plot(RFI_times,peakValues,xlab='Recovery time (ms)',ylab='Fraction of peak current',main=well,ylim=c(0,1),log='x')
      #lines(RFI_times,predictedPeak)
      #text(30,0.4,paste('Tau',round(result,2)))
    }
  }
  dev.off()
  par(mfrow=c(1,1))
  return(rfi2)
}
checkLateFile=function(filename,seqBy,firstVarCol){
  tempnumSweeps=100
  a=read.csv(filename,stringsAsFactors=FALSE,header=TRUE)
  head(a)
  colsToPull=seq(firstVarCol,firstVarCol+(tempnumSweeps-1)*seqBy,seqBy)
  colsToPull2=colsToPull[colsToPull<=ncol(a)]
  a2=a[3,colsToPull2]
  missing=which(!is.finite(as.numeric(a2)))
  numSweeps=length(colsToPull2)
  print(c(missing,numSweeps))
  return(list(names(a),as.character(a[1,])))
}
deltaBoltzmann=function(k,WTinact,mutInact){
  V=seq(-120,-50,by=0.1)
  P = 1/(1 + exp((V-WTinact) / k))
  val120_WT=1/(1 + exp((-120-WTinact) / k)) #0.9866131
  val90_WT=1/(1 + exp((-90-WTinact) / k)) #0.1544653
  percentAvailable_WT=val90_WT/val120_WT #0.1565611
  V_mut=seq(-120,-50,by=0.1)
  P_mut = 1/(1 + exp((V_mut-mutInact) / k))
  #plot(V,P)
  #points(V_mut,P_mut,col='red')
  val120_mut=1/(1 + exp((-120-mutInact) / k)) #0.9996646
  val90_mut=1/(1 + exp((-90-mutInact) / k)) #0.8807971
  
  percentAvailable_mut=val90_mut/val120_mut #0.8810926
  adjustment_mut=percentAvailable_mut/percentAvailable_WT
  return(adjustment_mut)
}
filterOutliers=function(df,var,vals,sdNum,verbose=T){
  # initial filter out big outliers >4SD from vals variable
  #print('MAX VALS:')
  #print(max(vals))
  initial_sd=4
  m=mean(vals,na.rm=TRUE)
  s=sd(vals,na.rm=TRUE)
  minVal=m-initial_sd*s
  maxVal=m+initial_sd*s
  #print(length(vals))
  initialLength=length(vals)
  vals=vals[vals>=minVal & vals<maxVal]
  numOutliersInitial=initialLength-length(vals)
  #print(length(vals))

  # Main curation
  m=mean(vals,na.rm=TRUE)
  s=sd(vals,na.rm=TRUE)
  minVal=m-sdNum*s
  maxVal=m+sdNum*s
  valsToFilter=as.numeric(df[,var])
  #print('AAAA')
  valsToFilter2=valsToFilter[!is.na(valsToFilter) & valsToFilter>minVal*3 & valsToFilter<minVal*3]
  #print('BBBB')
  keepTF=!is.na(valsToFilter) & valsToFilter>=minVal & valsToFilter<maxVal
  df2=df[keepTF,]
  range=max(vals)-min(vals)
  plotFraction=range/(maxVal-minVal)
  breaks=10000
  #if(plotFraction>2){
  #  breaks=100000
  #}
  if(verbose){
    xlim_max=m+sdNum*s*1.5
    xlim_min=m-sdNum*s*1.5
    #print(xlim_max)
    #print('CCCC')
    hist(valsToFilter,breaks=breaks,main=var,xlim=c(xlim_min,xlim_max))
    #print('DDDD')
    #hist(valsToFilter2,breaks=breaks,main=var,xlim=c(0,100))
    #print('MAX VALS:')
    #print(max(vals))
    abline(v=m,col='blue')
    abline(v=minVal,col='red')
    abline(v=maxVal,col='red')
    #print(paste('Mean:',round(m,3),'SD:',round(s,3)))
  }
  numOutliers=nrow(df)-nrow(df2)
  if(verbose){
    print(paste('Initial filter: filtered',numOutliersInitial,'outliers outside range of',round(initial_sd,3),'SD'))
    print(paste('Filtered',numOutliers,'outliers between',round(minVal,3),'and',round(maxVal,3)))
  }
  #print(nrow(df2))
  return(df2)
}
filterWeirdPlates=function(i,wtfilter=TRUE){
  i= filter(i,!str_detect(Plate,"^SP26") & !str_detect(Plate,"^SP40") & !str_detect(Plate,"^SP44"))
  if(wtfilter){
    i[i$Mutation=='WT2' | i$Mutation=='WT3','Mutation']='WT'
  }
  return(i)
}
fixCategories=function(i,mutList,newCategoryList){
  for(i in 1:length(mutList)){
    mut=mutList[i]
    newCategory=newCategoryList[i]
    i[i$Mutation==mut,'MutDescription']=newCategory
  }
  return(i)
}
  getMutDescriptionAndFilter=function(peakSummary2,reorderCols=TRUE){
  for(j in 1:nrow(peakSummary2)){
    mut=peakSummary2[j,'Mutation']
    desc=as.character(i[i$Mutation==mut,'MutDescription'])[1]
    peakSummary2[j,'MutDescription']=desc
  }
  if(reorderCols){
    peakSummary2=peakSummary2[,c(1,6,2,3,4,5)]
  }
  summary(as.factor(peakSummary2$MutDescription))
  peakSummary3=peakSummary2[!peakSummary2$MutDescription=='brs_drug' & !peakSummary2$MutDescription=='dms' & !peakSummary2$MutDescription=='empty' & !peakSummary2$MutDescription=='lpneg' & !peakSummary2$MutDescription=='WT_drug',]
  peakSummary4=peakSummary3[!peakSummary3$Mutation=='WT2' & !peakSummary3$Mutation=='WT3' & !peakSummary3$Mutation=='AirLeak' & !peakSummary3$Mutation=='toss' & !peakSummary3$Mutation=='R893C_double', ]
  return(peakSummary4)
}
getStats=function(a,varName='Peak',goodSealCol='GoodSeal'){
  varName2=paste0('numCells',varName)
  includeName=paste0(varName,'Include')
  mutFrame=data.frame(Mutation=unique(a$Mutation))
  for(row in 1:nrow(mutFrame)){
    mut=as.character(mutFrame[row,'Mutation'])
    print(mut)
    a2=a[a$Mutation==mut,]
    mutFrame[row,'MutDescription']=a2[1,'MutDescription']
    #mutFrame[row,'GoodSeal']=nrow(a2)-sum(a2$Type=='badSeal')
    mutFrame[row,'GoodSeal']=sum(a2[,goodSealCol])
    mutFrame[row,varName2]=sum(as.logical(a2[,includeName]))-sum(as.logical(a2[,includeName]) & a2$Type=='none')
    mutFrame[row,'Observed']=round(100*mutFrame[row,varName2]/mutFrame[row,'GoodSeal'],1)
  }
  return(mutFrame)
}
makeBarplot=function(mutFrame2,y,ySE,ylim,invertTF,numCells,outpdf,space=NA,topOnly=FALSE,col=NULL){
  invertFactor=1
  if(invertTF){
    invertFactor=-1
  }
  toplot=invertFactor*as.numeric(mutFrame2[,y])
  y1=toplot-invertFactor*(as.numeric(mutFrame2[,ySE]))
  y2=toplot+invertFactor*(as.numeric(mutFrame2[,ySE]))
  if(topOnly){
    y1=toplot
  }
  wtVal=invertFactor*mutFrame2[mutFrame2$Mutation=='WT',y]
  textLocation=ylim[1]+.95*(ylim[2]-ylim[1])
  print(textLocation)
  if(is.na(space)){
    x=as.numeric(barplot(toplot,ylim=ylim, las=2, cex.names=0.5))
  } else{
    x=as.numeric(barplot(toplot,ylim=ylim,space=space, las=2, cex.names=0.5))
  }
  pdf(outpdf)
  names(toplot)=mutFrame2$Mutation
  if(is.na(space)){
    barplot(toplot,ylab=y,main=y,ylim=ylim,xlab='Mutation',las=2,cex.names=0.5, col=col)
  } else{
    barplot(toplot,ylab=y,main=y,ylim=ylim,xlab='Mutation',las=2,cex.names=0.5, space=space,col=col)
  }
  segments(x,y1,x,y2)
  segments(0,wtVal,nrow(mutFrame2)*12.7/11,wtVal,lty=2,cex=2)
  scalingFactor=max(x)/length(x)
  if(!is.na(numCells)){
    for(i in 1:nrow(mutFrame2)){
      text(i*scalingFactor-.5,textLocation,mutFrame2[i,numCells])
    }
  }
  dev.off()
}

makeIVplots=function(a3,pdfname,IV_voltages,startCol,endCol,keepCriteria='GoodSeal',xlab='Voltage (uV)',scaleFactor=10^12,vertLines=NA){
  pdf(pdfname)
  par(mfrow=c(4,4))
  for(row in 1:nrow(a3)){
    well=a3[row,'Well']
    plate=a3[row,'Plate']
    title=paste(plate,well)
    print(well)
    includeA=a3[row,keepCriteria]
    if(includeA){
      if(all(is.nan(as.numeric(a3[row,startCol:endCol])))){
        plot(0,0,main=title,xlab='All NaN') #all NaNs, skip
        if(!is.na(vertLines)){
          for(val in vertLines){
            abline(v=val,col='blue')
          }
        }
      }
      else{
        plot(IV_voltages,scaleFactor*a3[row,startCol:endCol],main=title,ylab='Current (pA)',xlab=xlab)
        abline(h=0,lty=2)
        if(!is.na(vertLines)){
          for(val in vertLines){
            abline(v=val,col='blue')
          }
        }
      }
    }
  }
  dev.off()
  par(mfrow=c(1,1))
}
meanWithoutOutliers=function(vals,sdCutoff){
  rawMean=mean(vals)
  rawSD=sd(vals)
  cutoff1=rawMean-sdCutoff*rawSD
  cutoff2=rawMean+sdCutoff*rawSD
  vals2=vals[vals>=cutoff1 & vals<=cutoff2]
  newMean=mean(vals2)
  return(newMean)
}
orderMutFrame=function(mutFrame,orderBy,putFirst=NA,decreasing=FALSE){
  mutFrame2=mutFrame[order(mutFrame[,orderBy],decreasing=decreasing),]
  if(!is.na(putFirst)){
    WTrow=which(mutFrame2$Mutation==putFirst)[1]
    if(WTrow==1){
      return(mutFrame2)
    }
    pre=c()
    if(WTrow>1){
      pre=1:(WTrow-1)
    }
    post=c()
    if(WTrow<nrow(mutFrame2)){
      post=(WTrow+1):nrow(mutFrame2)
    }
    rows=c(WTrow,pre,post)
    mutFrame2=mutFrame2[rows,] # to have wt first
  }
  return(mutFrame2)
}
peak_concordance <- function(df,mutation)
{
  x = df[df$Mutation == mutation & df$PeakInclude,]
  y = names(summary(as.factor(x$Plate)))
  
  iii <- data.frame(matrix(ncol=11,nrow=1))
  colnames(iii)=c("mutation","R1","P1","R2","P2","R3","P3","R4","P4","R5","P5")
  iii[1,1] = mutation
  
  actualCount1=0
  actualCount2=1
  prev_Plate='xxxxx'
  
  # different plates on which the mutation was run
  for (counter in 1:length(y)){
    plate.col = counter*2
    value.col= counter*2+1
    current_PlateE = y[counter]
    current_Plate = substr(current_PlateE, 1,nchar(current_PlateE)-1)
    if(current_Plate==prev_Plate & counter>1){
      next()
    }
    else{
      actualCount1=actualCount1+2
      actualCount2=actualCount2+2
    }
    if(current_PlateE =="SP1" | current_PlateE =="SP2" | 
       current_PlateE =="SP3" | current_PlateE =="SP4" | 
       current_PlateE =="SP5"){
      file.name = paste('~/Dropbox/Andrew-Yuko/Syncropatch/',current_PlateE,
                        "/Analysis/",current_PlateE,"_mutationSummary.csv",sep="")
    }
    else if(current_PlateE == "SP6" | current_PlateE == "SP7"){
      current_Plate="SP6+7"
      file.name = paste('~/Dropbox/Andrew-Yuko/Syncropatch/',current_Plate,
                        "/Analysis/",current_Plate,"_mutationSummary.csv",sep="")
    }
    else if(current_PlateE == "SP8" | current_PlateE == "SP9"){
      current_Plate="SP8+9"
      file.name = paste('~/Dropbox/Andrew-Yuko/Syncropatch/',current_Plate,
                        "/Analysis/",current_Plate,"_mutationSummary.csv",sep="")
    }
    else{
      folder_to_check = as.numeric(substr(current_PlateE,3,nchar(current_PlateE)-1))
      if(folder_to_check < 20)
      {
        file.name = paste('~/Dropbox/Andrew-Yuko/Syncropatch/',current_Plate,
                          "/Analysis/",current_Plate,"_mutationSummary.csv",sep="")
      }
      else{
        file.name = paste('~/Dropbox/Andrew - Ayesha/Syncropatch/',current_Plate,
                          "/Analysis/",current_Plate,"_mutationSummary.csv",sep="")
      }
    }
    df2 <- read.csv(file.name)
    df3 <- df2[df2$Mutation == mutation,]
    if(actualCount1<12){
      iii[1,actualCount1]=current_Plate
      if(current_PlateE =="SP1" | current_PlateE =="SP2" | 
         current_PlateE =="SP3" | current_PlateE =="SP4" | 
         current_PlateE =="SP5"){
        iii[1,actualCount1]=current_PlateE
      }
      iii[1,actualCount2]=df3$PeakDensityNormMean
    }
    prev_Plate=current_Plate
  }
  return(iii)
}
processInitialCsv=function(a,numSweeps,seqBy,firstVarStart,secondVarStart=NA,capacitanceCol=NA,capMin=5e-12,capMax=30e-12,averageCap=FALSE,sealCol,sealMin=0.5*10^9,sealMax=10*10^9,plate=NA,firstVarLetter='A',secondVarLetter='B'){
  # This is a function  that processes...
  names(a)=1:ncol(a)
  a[,1]=substrRight(as.character(a[,1]),3)
  if(is.na(secondVarStart)){
    a2=a[,c(1,2,seq(firstVarStart,ncol(a),by=seqBy))]
    names(a2)=c('Well','GoodSeal',paste(firstVarLetter,1:numSweeps,sep=''))
  } else{
    a2=a[,c(1,2,seq(firstVarStart,ncol(a),by=seqBy),seq(secondVarStart,ncol(a),by=seqBy))]
    names(a2)=c('Well','GoodSeal',paste(firstVarLetter,1:numSweeps,sep=''),paste(secondVarLetter,1:numSweeps,sep=''))
  }
  if(!is.na(capacitanceCol)){
    if(as.logical(averageCap)){
      capFrame=a[,seq(capacitanceCol,ncol(a),by=seqBy)]
      capVals=c()
      for(j in 1:nrow(capFrame)){
        toMean=as.numeric(capFrame[j,])
        toMean=na.omit(toMean)
        capVal=mean(toMean,na.rm=TRUE)
        capVals=c(capVals,capVal)
      }
    }
    else{
      capVals=a[,capacitanceCol]
    }
    a3=cbind(a2[,1:2],capVals,a2[3:ncol(a2)])
    a2=a3
    names(a2)[3]='Capacitance'
    a2$Capacitance[is.nan(a2$Capacitance)]=NA
    a2$GoodSeal[a2$Capacitance<capMin]=FALSE
    a2$GoodSeal[a2$Capacitance>capMax]=FALSE
  }
  else{
    a2=cbind(a2[,1:2],Capacitance=rep(NA,nrow(a2)),a2[3:ncol(a2)])
  }
  if(!is.na(plate)){
    a2=cbind(plate,a2)
    names(a2)[1]="Plate"
    a2$Plate=as.character(a2$Plate)
  }
  a3=cbind(a2[,1:4],a[,sealCol],a2[(4+1):ncol(a2)])
  names(a3)[(4+1)]='Seal'
  a2=a3
  a2$GoodSeal[a2$Seal<sealMin]=FALSE
  a2$GoodSeal[a2$Seal>sealMax]=FALSE
  return(a2)
}
processMissing=function(i,m,result,name,numCols,workingDir='C:/Users/jessa/Dropbox/Andrew-Jessa/Syncropatch/'){
  i=paste0(workingDir,i)
  m=paste0(workingDir,m)
  
  # Process result list
  old_i2=result[[1]]
  empty=result[[2]]

  # Read in data and missing tables
  i2=read.csv(i,stringsAsFactors=FALSE)
  i2=i2[,1:numCols]
  m2=read.csv(m,stringsAsFactors=FALSE)
  # Merge include decisions and add in empty rows with 0 current for missing cells
  for(row in 1:nrow(i2)){
    mut=i2[row,'Mutation']
    includeTF=m2[m2$Mutation==mut,'Include']
    if(is.na(includeTF)){
      includeTF=FALSE
      m2[m2$Mutation==mut,'Include']=FALSE
    }
    i2[row,'Include']=includeTF
  }
  for(row in 1:nrow(m2)){
    mut=m2[row,'Mutation']
    if(m2[row,'Include']){
      whichRow=which(empty$Mutation==mut)
      if(length(whichRow)>0){
        empty[whichRow,'MissingCells']=empty[whichRow,'MissingCells']+m2[row,'MissingCells']
        oldSP=empty[whichRow,'SP']
        newSP=paste(oldSP,';',name,sep='')
        empty[whichRow,'SP']=newSP
      } else{
        empty=rbind(empty,data.frame(Mutation=mut,MissingCells=m2[row,'MissingCells'],SP=name))
      }
    }
  }
  if(!is.character(i2)){
    i2=rbind(old_i2,i2)
  }
  return(list(i2,empty))
}

processMissing2=function(i,m,result,name,numCols){
  # This function aggregates the data by (mutation and experiment) not by mutation
  # Process result list
  old_i2=result[[1]]
  empty=result[[2]]
  
  # Read in data and missing tables
  i2=read.csv(i,stringsAsFactors=FALSE)
  i2=i2[,1:numCols]
  m2=read.csv(m,stringsAsFactors=FALSE)
  for(row in 1:nrow(m2)){
    mut=m2[row,'Mutation']
    spname=name
    m2[row,'Mutation']=paste0(mut,"_",spname)
  }
  for(row in 1:nrow(i2)){
    mut=i2[row,'Mutation']
    spname=name
    i2[row,'Mutation']=paste0(mut,"_",spname)
  }
  # Merge include decisions and add in empty rows with 0 current for missing cells
  for(row in 1:nrow(i2)){
    mut=i2[row,'Mutation']
    includeTF=m2[m2$Mutation==mut,'Include']
    i2[row,'Include']=includeTF[1]
  }
  for(row in 1:nrow(m2)){
    mut=m2[row,'Mutation']
    if(m2[row,'Include']){
      whichRow=which(empty$Mutation==mut)
      if(length(whichRow)>0){
        empty[whichRow,'MissingCells']=empty[whichRow,'MissingCells']+m2[row,'MissingCells']
        oldSP=empty[whichRow,'SP']
        newSP=paste(oldSP,';',name,sep='')
        empty[whichRow,'SP']=newSP
      } else{
        empty=rbind(empty,data.frame(Mutation=mut,MissingCells=m2[row,'MissingCells'],SP=name))
      }
    }
  }
  if(!is.character(i2)){
    i2=rbind(old_i2,i2)
  }
  return(list(i2,empty))
}
processRFI=function(rfi1,rfi2,rfi3,act,cutoff,RFI_times,initial=5,range1=6:12,range2=6:12,range3=6:12,actRange){
  rfi=cbind(act[,actRange],goodRFI=NA,rfi1[,initial],rfi2[,initial],rfi3[,initial],rfi1[,range1],rfi2[,range2],rfi3[,range3])
  f=max(actRange) #last col of act
  print('F:')
  print(f)
  names(rfi)[(f+2):(f+4)]=c('P0_A','P0_B','P0_C')
  names(rfi)[(f+5):(f+25)]=paste('A',1:21,sep='')
  print(names(rfi))
  rfiCOV=calcCOV(rfi,(f+2):(f+4))
  rfi$goodRFI=act$Type=='good' & as.numeric(rfi[,'P0_C'])<cutoff & rfiCOV<.1 #93 total
  sweepList=1:21
  colList=c('01','02','03','04','05','06','07','08','09',10:24)
  rowList=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P')
  count=0
  for(i in (f+5):(f+25)){
    count=count+1
    if(i<(f+12)){
      P0col=f+2
    }
    if(i>=(f+12) & i<(f+19)){
      P0col=f+3
    }
    if(i>=(f+19)){
      P0col=f+4
    }
    rawData=rfi[,i]
    adjustedData=rawData/as.numeric(rfi[,P0col])
    rfi[,i+21]=adjustedData
    names(rfi)[i+21]=paste('B',count,sep='')
  }
  return(rfi)
}
remove_ending_letter <- function(x) {
  if (grepl("[a-zA-Z]$", x)) {
    return(sub("[a-zA-Z]$", "", x))
  } else {
    return(x)
  }
}
remove_ending_letter_vector = function(char_vector){
  result <- lapply(char_vector, remove_ending_letter)
  result <- unlist(result)
  return(result)
}

reorderMutations=function(df,peakOrder,varName){
  newOrder=c()
  for(currMut in peakOrder){
    currRow=which(as.character(df[,varName])==currMut)
    if(length(currRow)>0){
      newOrder=c(newOrder,currRow)
    }
  }
  df2=df[newOrder,]
  row.names(df2)=1:nrow(df2)
  return(df2)
}
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
summarizeByMutationGeneral=function(a3,variable,includeCol,scaleBy=1,missingFrame=NA,roundUnits=1){
  library(plotrix)
  mutList=unique(a3$Mutation)
  mutFrame=data.frame()
  goodTraceName=paste('numCells',variable,sep='')
  goodTraceName2=paste('numCells',variable,'WithMissing',sep='')
  meanName=paste(variable,'Mean',sep='')
  seName=paste(variable,'SE',sep='')
  count=0
  a3[,variable]=as.numeric(a3[,variable])
  for(mut in mutList){
    count=count+1
    a3b=a3[a3$Mutation==mut & !is.na(a3[,variable]) & as.logical(a3[,includeCol]),]
    mutFrame[count,'Mutation']=mut
    mutFrame[count,goodTraceName]=nrow(a3b)
    vars=a3b[,variable]
    vars[is.na(vars)] <- 0
    if(class(missingFrame)=='data.frame'){
      missingRow=which(missingFrame$Mutation==mut)
      if(length(missingRow)>0){
        missingCells=round(missingFrame[missingRow,'MissingCells'],0)
        missingCells[is.na(missingCells)] <- 0
        vars=c(vars,rep(0,missingCells))
        mutFrame[count,goodTraceName2]=length(vars)
      }
    }
    mutFrame[count,meanName]=round(scaleBy*mean(vars,na.rm=TRUE),roundUnits)
    mutFrame[count,seName]=round(scaleBy*std.error(vars,na.rm=TRUE),roundUnits)
  }
  if(class(missingFrame)=='data.frame'){
    print('Adding missing only')
    for(row in 1:nrow(missingFrame)){
      missingMut=missingFrame[row,'Mutation']
      if(sum(missingMut==mutFrame$Mutation,na.rm=TRUE)<1){
      #if(!any(missingMut==mutFrame$Mutation)){
     #   print(missingMut)
        count=count+1
        mutFrame[count,meanName]=0
        mutFrame[count,seName]=0
        mutFrame[count,goodTraceName]=0
        mutFrame[count,goodTraceName2]=missingFrame[row,'MissingCells']
        mutFrame[count,'Mutation']=missingFrame[row,'Mutation']
      }
    }
  }
  mutFrame[!is.finite(mutFrame[,meanName]),meanName]=NA
  mutFrame[!is.finite(mutFrame[,seName]),seName]=NA
  return(mutFrame)
}

filter_stats <- function(df,df2,col_nam){
# takes the original dataframe (df), and the dataframe with conventional filtering (df2)
# also takes a column name (e.g. inacttimeTau) based on which you want to filter 
# spits out a list of mutations for which conventional filtering removes more than 10% of the cells 
# (for any mutation with more than 10 cells)
  df=df[!is.na(df[,col_nam]),]
  pre=summary(as.factor(df$Mutation), maxsum=1000)
  post=summary(as.factor(df2$Mutation), maxsum=1000)
  pre2=data.frame(mut=names(pre),pre=pre)
  post2=data.frame(mut=names(post),post=post)
  library(tidyverse)
  counts=merge(pre2,post2,all=TRUE) %>%
    mutate(outliers=pre-post,
           outlier_percentage = round(100*outliers/pre,0)) %>%
    arrange(outlier_percentage)
  
  counts2 = counts[counts$outlier_percentage >=15 & counts$pre >=15,]
  counts2=counts2[!is.na(counts2$mut),]
  return(counts2)
}


fancy_outlier_removal <- function(df,counts2,var,sd, df_wt){
# this function takes a df of the form count, plate, well, mutation and all the EP parameters
# in the form of 1 qrow per cell. it also takes a list of mutations (counts2), that has more than
# 20% of the outliers filtered out based on the regular criteria 
# 
  list_of_mutations <- unique(df$Mutation)
  firsttime=T
  v = 0 
  for (x in 1:length(list_of_mutations)){
    curr = list_of_mutations[x]
    y = sum(curr == as.character(counts2$mut))
    df3 = df[df$Mutation==curr,]
    if(!y){
      rez = filterOutliers(df3,var,df_wt[,var],3,verbose=F)
    }
    if(y){
      rez = filterOutliers(df3,var,df3[,var],3,verbose=F)
      print(paste("Custom outlier filtration for mutation", curr))
      v = v+1
    }
    if(firsttime){
      rez_final = rez
      firsttime = F
    }
    else{
      rez_final = rbind(rez_final,rez)
    }
  }
  print(paste("Number of mutations fancy filtered:",v))
  return(rez_final)
}

exportSingleCell = function(i, filterFile, filterColumn, filename, empty=NA, peakVar1 = NA, peakVar2 = NA){
  filterVals <- i[,filterColumn]
  i2 = i[filterVals,]
  filterFrame <- read.csv(filterFile, stringsAsFactors = FALSE, header = FALSE)
  i3 = i2[i2$Mutation %in% filterFrame$V1, ]
  i4 = i3
  if(is.data.frame(empty)){
    empty2 = empty[empty$Mutation %in% filterFrame$V1, ]
    for(n in 1:nrow(empty2)){
      missingCellsNumber <- round(empty2[n, "MissingCells"])
      if(missingCellsNumber == 0){
        next
      }
      missingMutation <- empty2[n, "Mutation"]
      for(m in 1:missingCellsNumber){
        newrownumber = nrow(i4)+1
        i4[newrownumber, "Mutation"] = missingMutation
        i4[newrownumber, "PeakDensityNorm"] = 0
        i4[newrownumber, "PeakDensity"] = 0
        i4[newrownumber, "PeakDensityNormSQRT"] = 0
        if(!is.na(peakVar1)){
          i4[newrownumber, peakVar1] = 0
        }
        if(!is.na(peakVar2)){
          i4[newrownumber, peakVar2] = 0
        }
      }
    }
  }
  print(dim(i4))
  i4 <- i4[order(i4$Mutation, i4$Count), ]
  write.csv(i4, filename, row.names=FALSE)
}



exportSingleCell78 = function(i, filterFile, filterColumn, filename, empty=NA, peakVar1 = NA, peakVar2 = NA){
  filterVals <- i[,filterColumn]
  print(dim(i))
  i2 = i[filterVals,]
  filterFrame <- read.csv(filterFile, stringsAsFactors = FALSE, header = FALSE)
  print(dim(i2))
  i3 = i2[i2$Mutation %in% filterFrame$V1, ]
  print(dim(i3))
  i4 = i3
  if(is.data.frame(empty)){
    empty2 = empty[empty$Mutation %in% filterFrame$V1, ]
    for(n in 1:nrow(empty2)){
      missingCellsNumber <- round(empty2[n, "MissingCells"])
      print(n)
      if(missingCellsNumber == 0){
        next
      }
      missingMutation <- empty2[n, "Mutation"]
      #print("check1")
      for(m in 1:missingCellsNumber){
        newrownumber = nrow(i4)+1
        i4[newrownumber, "Mutation"] = missingMutation
        i4[newrownumber, "Peak78DensityNorm"] = 0
        i4[newrownumber, "Peak78Density"] = 0
        i4[newrownumber, "Peak78DensityNormSQRT"] = 0
        if(!is.na(peakVar1)){
          i4[newrownumber, peakVar1] = 0
        }
        if(!is.na(peakVar2)){
          i4[newrownumber, peakVar2] = 0
        }
        #print("check2")
      }
    }
  }
  #print(dim(i4))
  i4 <- i4[order(i4$Mutation, i4$Count), ]
  write.csv(i4, filename, row.names=FALSE)
}

exportSingleCell90bulk = function(i, filterFile, filterColumn, filename, empty=NA, peakVar1 = NA, peakVar2 = NA){
  filterVals <- i[,filterColumn]
  i2 = i[filterVals,]
  filterFrame <- read.csv(filterFile, stringsAsFactors = FALSE, header = FALSE)
  i3 = i2[i2$Mutation %in% filterFrame$V1, ]
  i4 = i3
  if(!is.na(empty)){
    empty2 = empty[empty$Mutation %in% filterFrame$V1, ]
    for(n in 1:nrow(empty2)){
      missingCellsNumber <- round(empty2[n, "MissingCells"])
      if(missingCellsNumber == 0){
        next
      }
      missingMutation <- empty2[n, "Mutation"]
      #print("check1")
      for(m in 1:missingCellsNumber){
        newrownumber = nrow(i4)+1
        i4[newrownumber, "Mutation"] = missingMutation
        i4[newrownumber, "PeakDensityNormBulk90"] = 0
        i4[newrownumber, "PeakDensityNormBulk90SQRT"] = 0
        if(!is.na(peakVar1)){
          i4[newrownumber, peakVar1] = 0
        }
        if(!is.na(peakVar2)){
          i4[newrownumber, peakVar2] = 0
        }
        #print("check2")
      }
    }
  }
  #print(dim(i4))
  i4 <- i4[order(i4$Mutation, i4$Count), ]
  write.csv(i4, filename, row.names=FALSE)
}


remove_terminal_letter = function(my_string){
  if (grepl("[a-zA-Z]$", my_string)) {
    my_string <- substr(my_string, 1, nchar(my_string) - 1)
  }
  return(my_string)
}

fixSP6789names=function(allAct){
  allAct[allAct$Plate=='SP6','Plate']='SP6+7'
  allAct[allAct$Plate=='SP7','Plate']='SP6+7'
  allAct[allAct$Plate=='SP8','Plate']='SP8+9'
  allAct[allAct$Plate=='SP9','Plate']='SP8+9'
  return(allAct)  
}

# Code to perform per-plate sqrt transformation, followed by normalization

sqrt_transform_120 = function(i){
  for(plate in unique(i$Plate)){
    plate2 <- remove_terminal_letter(plate)
    i$Plate[i$Plate == plate] <- plate2
  }
  firsttime=TRUE
  for(plate in unique(i$Plate)){
    print(plate)
    i2 <- i[i$Plate == plate, ]
    i2$PeakDensitySQRTraw <- -1*sqrt(-1*i2$PeakDensity)
    i2_WT <- i2[i2$Mutation == "WT" & i2$PeakInclude & !is.na(i2$PeakDensitySQRTraw), ]
    sqrt_norm_factor <- mean(i2_WT$PeakDensitySQRTraw)
    i2$PeakDensityNormSQRT <- (i2$PeakDensitySQRTraw/sqrt_norm_factor) * 100
    if(firsttime){
      j <- i2
      firsttime = FALSE
    }
    else{
      j <- rbind(j, i2)
    }
  }
  return(j)
}

sqrt_transform_78 = function(i,Positive2Zero=FALSE){
  if(Positive2Zero){
    print('A')
    i[!is.na(i$Peak78Density) & i$Peak78Density>0,'Peak78Density']=0
    print('B')
  }
  for(plate in unique(i$Plate)){
    plate2 <- remove_terminal_letter(plate)
    i$Plate[i$Plate == plate] <- plate2
  }
  firsttime=TRUE
  for(plate in unique(i$Plate)){
    print(plate)
    i2 <- i[i$Plate == plate, ]
    i2$Peak78DensitySQRTraw <- -1*sqrt(-1*i2$Peak78Density)
    i2_WT <- i2[i2$Mutation == "WT" & i2$Peak78Include & !is.na(i2$Peak78DensitySQRTraw), ]
    sqrt_norm_factor <- mean(i2_WT$Peak78DensitySQRTraw)
    i2$Peak78DensityNormSQRT <- (i2$Peak78DensitySQRTraw/sqrt_norm_factor) * 100
    if(firsttime){
      j <- i2
      firsttime = FALSE
    }
    else{
      j <- rbind(j, i2)
    }
  }
  return(j)
}

sqrt_transform_90bulk = function(i,Negative2Zero=FALSE){
  if(Negative2Zero){
    #print('A')
    i[!is.na(i$PeakDensityNormBulk90) & i$PeakDensityNormBulk90<0,'PeakDensityNormBulk90']=0
    #print('B')
  }
  for(plate in unique(i$Plate)){
    plate2 <- remove_terminal_letter(plate)
    i$Plate[i$Plate == plate] <- plate2
  }
  firsttime=TRUE
  for(plate in unique(i$Plate)){
    print(plate)
    i2 <- i[i$Plate == plate, ]
    i2$PeakDensityNormBulk90SQRT <- sqrt(i2$PeakDensityNormBulk90)
    i2_WT <- i2[i2$Mutation == "WT" & i2$PeakInclude & !is.na(i2$PeakDensityNormBulk90SQRT), ]
    sqrt_norm_factor <- mean(i2_WT$PeakDensityNormBulk90SQRT)
    i2$PeakDensityNormBulk90SQRT <- (i2$PeakDensityNormBulk90SQRT/sqrt_norm_factor) * 100
    if(firsttime){
      j <- i2
      firsttime = FALSE
    }
    else{
      j <- rbind(j, i2)
    }
  }
  return(j)
}

convert_z_acmg = function(z_score){
  classes = c()
  for(i in 1:length(z_score)){
    val = z_score[i]
    if(is.na(val)){
      class = NA
      classes = c(classes, class)
      next
    }
    if(val < -4){
      class = "PS3_strong"
    }
    if(val < -3 & val >= -4){
      class = "PS3_moderate"
    }
    if(val < -2 & val >= -3){
      class = "PS3_supporting"
    }
    if(val < -1 & val >= -2){
      class = "BS3_supporting"
    }
    if(val < 1 & val >= -1){
      class = "BS3_moderate"
    }
    if(val < 2 & val >= 1){
      class = "BS3_supporting"
    }
    if(val > 2){
      class = NA
    }
    classes = c(classes, class)
  }
  return(classes)
}

mergeEmpty=function(oldEmpty,empty){
  library(dplyr)
  
  # Add a row index to each data frame
  empty <- empty %>% mutate(row_index = row_number())
  oldEmpty <- oldEmpty %>% mutate(row_index = row_number())
  
  # Perform a full join on the Mutation column
  merged_df <- full_join(empty, oldEmpty, by = "Mutation", suffix = c(".empty", ".oldEmpty"))
  
  # Replace NA values with 0 for MissingCells and empty strings for SP
  merged_df <- merged_df %>%
    mutate(
      MissingCells.empty = ifelse(is.na(MissingCells.empty), 0, MissingCells.empty),
      MissingCells.oldEmpty = ifelse(is.na(MissingCells.oldEmpty), 0, MissingCells.oldEmpty),
      SP.empty = ifelse(is.na(SP.empty), "", SP.empty),
      SP.oldEmpty = ifelse(is.na(SP.oldEmpty), "", SP.oldEmpty)
    )
  
  # Sum the MissingCells values and concatenate the SP values
  merged_df <- merged_df %>%
    mutate(
      MissingCells = MissingCells.empty + MissingCells.oldEmpty,
      SP = paste(SP.oldEmpty, SP.empty, sep = "")
    ) %>%
    select(Mutation, MissingCells, SP, row_index.empty, row_index.oldEmpty) # Select the required columns
  
  # Order the merged dataframe by the row indices
  merged_df <- merged_df %>%
    arrange(row_index.oldEmpty,row_index.empty) %>%
    select(Mutation, MissingCells, SP) # Remove the row_index columns
}

