# Combined Analysis of Individual SCN5A SP Experiments 2/18/25

# NOTE: This analysis file matches the analysis in the O'Neill et al manuscript.

# STEP 1: Combine each parameter from each experiment
# STEP 2: Analyze these combined datasets, perform outlier removal when appropriate
# STEP 3: Merge each parameter for final dataset used in manuscript and variant interpretation

library(stringr)
library(tidyverse)
library(ggplot2)

##### STEP 0: Set working directory #####
# If on Jessa's computer:
dir="C:/Users/jessa/Dropbox/Andrew-Jessa/Syncropatch/"
# If on Andrew's computer:
dir='~/Dropbox/Andrew-Jessa/Syncropatch/'
setwd(dir)

########## STEP 1: Combine each parameter from each experiment ##########
# STEP 1A: Combine SPs together (peak current and voltage of activation)
# STEP 1B: Combine SPs together (voltage of inactivation)
# STEP 1C: Combine SPs together (late current)
# STEP 1D: Combine SPs together (Recovery From Inactivation)
# STEP 1E: Combine SPs together (Inactivation time)
# STEP 1F: Combine SPs together (Peak current -90mV)

##### STEP 1A: Combine SPs together (peak current and voltage of activation) #####
#SP96+
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
result=processMissing('SP97/Analysis/SP97_NaIV_processed3.csv','SP97/Analysis/SP97_missing.csv',result,'SP97',48,dir)
result=processMissing('SP99/Analysis/SP99_NaIV_processed3.csv','SP99/Analysis/SP99_missing.csv',result,'SP99',48,dir)
result=processMissing('SP100/Analysis/SP100_NaIV_processed3.csv','SP100/Analysis/SP100_missing.csv',result,'SP100',48,dir)
allAct=result[[1]]
empty=result[[2]]
allAct=allAct[!allAct$Type=='start',] #remove "start line"
allAct$PeakInclude=as.logical(allAct$PeakInclude)
allAct$VHalfActInclude=as.logical(allAct$VHalfActInclude)
allAct[is.na(allAct$Capacitance),'PeakInclude']=FALSE

#Combining with SP5-95
oldAct=read.csv(paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_act.csv'))
allAct=rbind(oldAct,allAct)
write.csv(allAct, 'CombinedAnalysis/SP5-95/SP5-95_allCells_act.csv',row.names=FALSE,quote=FALSE)
oldEmpty=read.csv(paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_missingCells.csv'))
empty=mergeEmpty(oldEmpty,empty)
write.csv(empty,paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_missingCells.csv'),row.names=FALSE,quote=FALSE)

##### STEP 1B: Combine SPs together (voltage of inactivation) #####
#SP96+
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
result=processMissing('SP97/Analysis/SP97_Inac_processed2.csv','SP97/Analysis/SP97_missing.csv',result,'SP97',61,dir)
result=processMissing('SP99/Analysis/SP99_Inac_processed2.csv','SP99/Analysis/SP99_missing.csv',result,'SP99',61,dir)
result=processMissing('SP100/Analysis/SP100_Inac_processed2.csv','SP100/Analysis/SP100_missing.csv',result,'SP100',61,dir)
allInact=result[[1]]
empty=result[[2]]
allInact=allInact[!allInact$Type=='start',] #remove "start line"
allInact$PeakInclude=as.logical(allInact$PeakInclude)
allInact[is.na(allInact$Capacitance),'PeakInclude']=FALSE
allInact$VHalfInactInclude=as.logical(allInact$VHalfInactInclude)

#Combining with SP5-95
oldInact=read.csv(paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_inact.csv'))
allInact=rbind(oldInact,allInact)
write.csv(allInact, paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_inact.csv'),row.names=FALSE,quote=FALSE)

##### STEP 1C: Combine SPs together (late current) #####
#SP96+
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
result=processMissing('SP97/Analysis/SP97_Late_processed2.csv','SP97/Analysis/SP97_missing.csv',result,'SP97',45,dir)
result=processMissing('SP99/Analysis/SP99_Late_processed2.csv','SP99/Analysis/SP99_missing.csv',result,'SP99',45,dir)
result=processMissing('SP100/Analysis/SP100_Late_processed2.csv','SP100/Analysis/SP100_missing.csv',result,'SP100',45,dir)
allLate=result[[1]]
empty=result[[2]]
allLate=allLate[!allLate$Type=='start',] #remove "start line"
# We had a weird late current protocol for some SP runs in the 60's--these generate extra parameters
# So create a proxy data frame with NAs to align names correctly
proxy <- data.frame(matrix(ncol = 11, nrow = nrow(allLate)))
colnames(proxy) <- c("Peak_30Pre_Late", "Peak_30Post_Late", "Peak_30_Late", "Late200_30Pre_Late", "Late200_30Post_Late", "Late200_30_Late", "Late200_30Ratio_Late", "Late50_30Pre_Late", "Late50_30Post_Late", "Late50_30_Late", "Late50_30Ratio_Late")
allLate2 <- cbind(allLate[,1:45], proxy, Include = allLate[,46])
allLate2$PeakInclude=as.logical(allLate2$PeakInclude)
allLate2[is.na(allLate2$Capacitance),'PeakInclude']=FALSE

#Combining with SP5-95
oldLate=read.csv(paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_late.csv'))
allLate2=rbind(oldLate,allLate2)
write.csv(allLate2, paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_late.csv'),row.names=FALSE,quote=FALSE)

##### STEP 1D: Combine SPs together (Recovery From Inactivation) ##### 
#SP96+
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
result=processMissing('SP97/Analysis/SP97_RecInact_processed2.csv','SP97/Analysis/SP97_missing.csv',result,'SP97',95,dir)
result=processMissing('SP99/Analysis/SP99_RecInact_processed2.csv','SP99/Analysis/SP99_missing.csv',result,'SP99',95,dir)
result=processMissing('SP100/Analysis/SP100_RecInact_processed2.csv','SP100/Analysis/SP100_missing.csv',result,'SP100',95,dir)
allRFI=result[[1]]
empty=result[[2]]
allRFI=allRFI[!allRFI$Type=='start',] #remove "start line"
allRFI$PeakInclude=as.logical(allRFI$PeakInclude)
allRFI[is.na(allRFI$Capacitance),'PeakInclude']=FALSE

#Combining with SP5-95
oldRFI=read.csv('C:/Users/jessa/Dropbox/Andrew-Jessa/Syncropatch/CombinedAnalysis/SP5-95/SP5-95_allCells_rfi.csv')
allRFI=rbind(oldRFI,allRFI)
write.csv(allRFI, paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_RFI.csv'),row.names=FALSE,quote=FALSE)

##### STEP 1E: Combine SPs together (Inactivation time) #####
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
result=processMissing('SP97/Analysis/SP97_Inacttime_processed.csv','SP97/Analysis/SP97_missing.csv',result,'SP97',53,dir)
result=processMissing('SP99/Analysis/SP99_Inacttime_processed.csv','SP99/Analysis/SP99_missing.csv',result,'SP99',53,dir)
result=processMissing('SP100/Analysis/SP100_Inacttime_processed.csv','SP100/Analysis/SP100_missing.csv',result,'SP100',53,dir)
allInacttime=result[[1]]
empty=result[[2]]
allInacttime=allInacttime[!allInacttime$Type=='start',] #remove "start line"
allInacttime$PeakInclude=as.logical(allInacttime$PeakInclude)
allInacttime[is.na(allInacttime$Capacitance),'PeakInclude']=FALSE

#Combining with SP5-95
oldInacttime=read.csv(paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_inacttime.csv'))
allInacttime=rbind(oldInacttime,allInacttime)
write.csv(allInacttime, paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_inacttime.csv'),row.names=FALSE,quote=FALSE)

##### STEP 1F: Combine SPs together (Peak current -90mV) #####
empty=data.frame(Mutation=c('WT'),MissingCells=c(0),SP='',stringsAsFactors=FALSE) #initialize empty data frame
result=list('start',empty)
result=processMissing('SP97/Analysis/SP97_Peak90_processed2.csv','SP97/Analysis/SP97_missing.csv',result,'SP97',61,dir)
result=processMissing('SP99/Analysis/SP99_Peak90_processed2.csv','SP99/Analysis/SP99_missing.csv',result,'SP99',61,dir)
result=processMissing('SP100/Analysis/SP100_Peak90_processed2.csv','SP100/Analysis/SP100_missing.csv',result,'SP100',61,dir)
allPeak90=result[[1]]
empty=result[[2]]
allPeak90=allPeak90[!allPeak90$Type=='start',] #remove "start line"
allPeak90$PeakInclude=as.logical(allPeak90$PeakInclude)
allPeak90[is.na(allPeak90$Capacitance),'PeakInclude']=FALSE

#Combining with SP5-95
oldPeak90=read.csv(paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_allCells_Peak90.csv'))
allPeak90=rbind(oldPeak90,allPeak90)
file_path <- file.path(dir, 'SP5-95_allCells_Peak90.csv')
write.csv(allPeak90, file_path, row.names = FALSE, quote = FALSE)
oldEmpty=read.csv(paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_missingCells2.csv'))
empty=mergeEmpty(oldEmpty,empty)
write.csv(empty,paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_missingCells2.csv'),row.names=FALSE,quote=FALSE)


########## STEP 2: Analyze these combined datasets, perform outlier removal when appropriate ##########


# STEP 2A) Calculate -120 mV peak current (no sqrt transformation)
# STEP 2B) Calculate -120 mV peak current (with sqrt transformation)
# STEP 2C) Calculate Recovery from Inactivation
# STEP 2D) Calculate Inactivation time
# STEP 2E) Calculate Voltage of Activation
# STEP 2F) Calculate Voltage of Inactivation
# STEP 2G) Calculate Late current
# STEP 2H) Calculate -90 mV peak current using bulk inactivation number (no sqrt transformation)
# STEP 2I) Calculate -90 mV peak current using bulk inactivation number (with sqrt transformation)

##### STEP 2A: Calculate -120 mV peak current (no sqrt transformation) #####
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_act.csv',stringsAsFactors=FALSE)
empty=read.csv('CombinedAnalysis/SP5-95/SP5-95_missingCells.csv',stringsAsFactors = FALSE)
peakSummary=summarizeByMutationGeneral(i,'PeakDensityNorm','PeakInclude',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'PeakDensityNormMean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(peakSummary3,paste0(dir,'CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormData.csv'),row.names=FALSE)

##### STEP 2B: Calculate -120 mV peak current (with sqrt transformation) #####
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_act.csv',stringsAsFactors=FALSE)
empty=read.csv('CombinedAnalysis/SP5-95/SP5-95_missingCells.csv',stringsAsFactors = FALSE)
j = sqrt_transform_120(i) #creates 2 variables, PeakDensitySQRTraw and PeakDensityNormSQRT
peakSummary=summarizeByMutationGeneral(j,'PeakDensityNormSQRT','PeakInclude',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'PeakDensityNormSQRTMean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(peakSummary3,'CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormSQRTData.csv',row.names=FALSE)

##### STEP 2C: Calculate Recovery from Inactivation  #####
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_rfi.csv',stringsAsFactors=FALSE)
i=filterWeirdPlates(i)
i_good=i[i$rfiInclude,]
i_fixed=read.csv('CombinedAnalysis/SPfixed/SPfixed_allCells_rfi.csv',stringsAsFactors=FALSE)
i_fixed=filterWeirdPlates(i_fixed)
i_wtgood=i_fixed[i_fixed$rfiInclude & i_fixed$Mutation=='WT',]
i2_good=filterOutliers(i_good,'rfi50',i_wtgood$rfi50,3)
counts2 = filter_stats(i_good,i2_good, "rfi50")
i3_good <- fancy_outlier_removal(i_good,counts2, 'rfi50',3, i_wtgood)
i4=calcFancyDelta(i3_good,'rfi50','rfi50Delta')
RFISummary=summarizeByMutationGeneral(i4,'rfi50','rfiInclude',roundUnits=3)
RFISummaryDelta=summarizeByMutationGeneral(i4,'rfi50Delta','rfiInclude',roundUnits=3)[,c(1,3,4)]
RFISummary=merge(RFISummary,RFISummaryDelta,all=TRUE)
RFISummary2=orderMutFrame(RFISummary,'rfi50DeltaMean','WT',decreasing=TRUE)
RFISummary3=RFISummary2[RFISummary2$numCellsrfi50>=5,]
RFISummary4=getMutDescriptionAndFilter(RFISummary3,reorderCols=FALSE)
RFISummary4=RFISummary4[,c('Mutation','MutDescription','rfi50Mean','rfi50SE','rfi50DeltaMean','rfi50DeltaSE','numCellsrfi50')]
write.csv(RFISummary4,'CombinedAnalysis/SP5-95/SP5-95_RFI_all.csv',row.names=FALSE)

##### STEP 2D: Calculate Inactivation time #####
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_inacttime.csv',stringsAsFactors=FALSE)
i=filterWeirdPlates(i)
i_good=i[!is.na(i$inacttimeError) & !is.na(i$inacttimeTau) & i$inacttimeError<.1 & i$inacttimeTau>0,]
i_fixed=read.csv('CombinedAnalysis/SPfixed/SPfixed_allCells_inacttime.csv',stringsAsFactors=FALSE)
i_fixed=filterWeirdPlates(i_fixed)
i_fixed_good=i_fixed[!is.na(i_fixed$inacttimeError) & !is.na(i_fixed$inacttimeTau) & i_fixed$inacttimeError<.1 & i_fixed$inacttimeTau>0,]
i_wtgood=i_fixed_good[i_fixed_good$Mutation=='WT',]
sum(i_wtgood$inacttimeTau>100)
i_wtgood=i_wtgood[i_wtgood$inacttimeTau<100,]
hist(i_wtgood$inacttimeTau,breaks=50)
i2_good=filterOutliers(i_good,'inacttimeTau',i_wtgood$inacttimeTau,3)
counts2 = filter_stats(i_good,i2_good, "inacttimeTau")
i3_good <- fancy_outlier_removal(i_good,counts2, 'inacttimeTau', 3, i_wtgood)
i4=calcFancyDelta(i3_good,'inacttimeTau','inacttimeTauDelta')
InacttimeSummary=summarizeByMutationGeneral(i4,'inacttimeTau','inacttimeInclude',roundUnits=3)
InacttimeSummaryDelta=summarizeByMutationGeneral(i4,'inacttimeTauDelta','inacttimeInclude',roundUnits=3)[,c(1,3,4)]
InacttimeSummary=merge(InacttimeSummary,InacttimeSummaryDelta,all=TRUE)
InacttimeSummary2=orderMutFrame(InacttimeSummary,'inacttimeTauMean','WT',decreasing=TRUE)
InacttimeSummary3=InacttimeSummary2[InacttimeSummary2$numCellsinacttimeTau>=5,]
InacttimeSummary4=getMutDescriptionAndFilter(InacttimeSummary3,reorderCols=FALSE)
InacttimeSummary4=InacttimeSummary4[,c('Mutation','MutDescription','inacttimeTauMean','inacttimeTauSE','inacttimeTauDeltaMean','inacttimeTauDeltaSE','numCellsinacttimeTau')]
write.csv(InacttimeSummary4,'CombinedAnalysis/SP5-95/SP5-95_inacttime_all.csv',row.names=FALSE)

##### STEP 2E: Calculate Voltage of Activation #####
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_act.csv',stringsAsFactors=FALSE)
i=filterWeirdPlates(i)
i_good=i[i$VHalfActInclude & i$Peak > (-3e-09),]  # Only including cells with peak current <3pA
i_fixed=read.csv('CombinedAnalysis/SPfixed/SPfixed_allCells_act.csv',stringsAsFactors=FALSE)
i_fixed=filterWeirdPlates(i_fixed)
i_fixed_good=i_fixed[i_fixed$VHalfActInclude & i_fixed$Peak > (-3e-09),]  # Only including cells with peak current <3pA
i_wtgood=i_fixed_good[i_fixed_good$Mutation=='WT',]
i2_good=filterOutliers(i_good,'VHalfAct',i_wtgood$VHalfAct,3)
counts2 = filter_stats(i_good,i2_good, "VHalfAct")
i3_good <- fancy_outlier_removal(i_good,counts2, 'VHalfAct',3, i_wtgood)
i4=calcFancyDelta(i3_good,'VHalfAct','VHalfActDelta')
actSummary=summarizeByMutationGeneral(i3_good,'VHalfAct','VHalfActInclude',roundUnits=3)
actSummaryDelta=summarizeByMutationGeneral(i4,'VHalfActDelta','VHalfActInclude',roundUnits=3)[,c(1,3,4)]
actSummary=merge(actSummary,actSummaryDelta,all=TRUE)
actSummary2=orderMutFrame(actSummary,'VHalfActMean','WT',decreasing=TRUE) #500
actSummary3=actSummary2[actSummary2$numCellsVHalfAct>=5,] #438
actSummary4=getMutDescriptionAndFilter(actSummary3,reorderCols=FALSE)
actSummary4=actSummary4[,c('Mutation','MutDescription','VHalfActMean','VHalfActSE','VHalfActDeltaMean','VHalfActDeltaSE','numCellsVHalfAct')]
write.csv(actSummary4,'CombinedAnalysis/SP5-95/SP5-95_act_all.csv',row.names=FALSE)

##### STEP 2F: Calculate Voltage of Inactivation #####
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_inact.csv',stringsAsFactors=FALSE)
i=filterWeirdPlates(i)
i_good=i[i$VHalfInactInclude,]
i_fixed=read.csv('CombinedAnalysis/SPfixed/SPfixed_allCells_inact.csv',stringsAsFactors=FALSE)
i_fixed=filterWeirdPlates(i_fixed)
i_fixed_good=i_fixed[i_fixed$VHalfInactInclude,]
i_wtgood=i_fixed_good[i_fixed_good$Mutation=='WT',]
i2_good=filterOutliers(i_good,'VHalfInact',i_wtgood$VHalfInact,3)
counts2 = filter_stats(i_good,i2_good, "VHalfInact")
i3_good <- fancy_outlier_removal(i_good,counts2, 'VHalfInact',3, i_wtgood)
i4=calcFancyDelta(i3_good,'VHalfInact','VHalfInactDelta')
inactSummary=summarizeByMutationGeneral(i4,'VHalfInact','VHalfInactInclude')
inactSummaryDelta=summarizeByMutationGeneral(i4,'VHalfInactDelta','VHalfInactInclude',roundUnits=3)[,c(1,3,4)]
inactSummary=merge(inactSummary,inactSummaryDelta,all=TRUE)
inactSummary2=orderMutFrame(inactSummary,'VHalfInactMean','WT',decreasing=TRUE)
inactSummary3=inactSummary2[inactSummary2$numCellsVHalfInact>=5,]
inactSummary4=getMutDescriptionAndFilter(inactSummary3,reorderCols=FALSE)
inactSummary4=inactSummary4[,c('Mutation','MutDescription','VHalfInactMean','VHalfInactSE','VHalfInactDeltaMean','VHalfInactDeltaSE','numCellsVHalfInact')]
write.csv(inactSummary4,'CombinedAnalysis/SP5-95/SP5-95_inact_all.csv',row.names=FALSE)

##### STEP 2G: Calculate late current ##### 
peakCutoff=-1000e-12
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_late.csv',stringsAsFactors=FALSE)
i=filterWeirdPlates(i)
i_good=i[i$GoodSeal_Late & i$LateInclude &!is.na(i$GoodSeal_Late) & i$Peak_Late<peakCutoff & i$PeakPost_Late/i$PeakPre_Late<0.1,]
i_fixed=read.csv('CombinedAnalysis/SPfixed/SPfixed_allCells_late.csv',stringsAsFactors=FALSE)
i_fixed=filterWeirdPlates(i_fixed)
i_fixed_good=i_fixed[i_fixed$GoodSeal_Late & i_fixed$LateInclude & (!is.na(i_fixed$GoodSeal_Late))  & i_fixed$Peak_Late<peakCutoff & i_fixed$PeakPost_Late/i_fixed$PeakPre_Late<0.1,]
i_wtgood=i_fixed_good[i_fixed_good$Mutation=='WT',]

# Plot Late current at 50 ms and 200 ms
forplot=i_wtgood[i_wtgood$Late200Ratio_Late>0,]
forplot$Peak_Late2=(10^12)*forplot$Peak_Late
plot(forplot$Peak_Late2,forplot$Late200Ratio_Late,xlab='Peak current (pA)',ylab='late ratio',ylim=c(-.01,.05),pch='.',xlim=c(-10000,0))
abline(h=0)
abline(v=-1000)
plot(forplot$Peak_Late2,forplot$Late50Ratio_Late,xlab='Peak current (pA)',ylab='late ratio',ylim=c(-.01,.05),pch='.',xlim=c(-10000,0))
abline(h=0)
abline(v=-1000)

# Late50
i_wtgood_50=i_wtgood[i_wtgood$Late50Ratio_Late>0,]
i_good_50=i_good[i_good$Late50Ratio_Late>0,]
i2_good_50=filterOutliers(i_good_50,'Late50Ratio_Late',i_wtgood_50$Late50Ratio_Late,3)
counts2 = filter_stats(i_good_50,i2_good_50, "Late50Ratio_Late")
i3_good_50 <- fancy_outlier_removal(i_good_50,counts2, 'Late50Ratio_Late',3, i_wtgood_50)
i4_50=calcFancyDelta(i3_good_50,'Late50Ratio_Late','Late50Ratio_LateDelta')
lateSummary_50=summarizeByMutationGeneral(i4_50,'Late50Ratio_Late','LateInclude',roundUnits=4)
lateSummaryDelta_50=summarizeByMutationGeneral(i4_50,'Late50Ratio_LateDelta','LateInclude',roundUnits=4)[,c(1,3,4)]
lateSummary_50=merge(lateSummary_50,lateSummaryDelta_50,all=TRUE)
lateSummary_50=lateSummary_50[lateSummary_50$numCellsLate50Ratio_Late>=5,]

# Late200
i_wtgood_200=i_wtgood[i_wtgood$Late200Ratio_Late>0,]
i_good_200=i_good[i_good$Late200Ratio_Late>0,]
i2_good_200=filterOutliers(i_good_200,'Late200Ratio_Late',i_wtgood_200$Late200Ratio_Late,3)
counts2 = filter_stats(i_good_200,i2_good_200, "Late200Ratio_Late")
i3_good_200 <- fancy_outlier_removal(i_good_200,counts2, 'Late200Ratio_Late',3, i_wtgood_200)
i4_200=calcFancyDelta(i3_good_200,'Late200Ratio_Late','Late200Ratio_LateDelta')
lateSummary_200=summarizeByMutationGeneral(i4_200,'Late200Ratio_Late','LateInclude',roundUnits=4)
lateSummaryDelta_200=summarizeByMutationGeneral(i4_200,'Late200Ratio_LateDelta','LateInclude',roundUnits=4)[,c(1,3,4)]
lateSummary_200=merge(lateSummary_200,lateSummaryDelta_200,all=TRUE)
lateSummary_200=lateSummary_200[lateSummary_200$numCellsLate200Ratio_Late>=5,]

# RampMin
i_wtgood_ramp=i_wtgood[i_wtgood$RampMinRatio_Late>0,]
i_good_ramp=i_good[i_good$RampMinRatio_Late>0,]
i2_good_ramp=filterOutliers(i_good_ramp,'RampMinRatio_Late',i_wtgood_ramp$RampMinRatio_Late,3)
counts2 = filter_stats(i_good_ramp,i2_good_ramp, "RampMinRatio_Late")
i3_good_ramp <- fancy_outlier_removal(i_good_ramp,counts2, 'RampMinRatio_Late',3, i_wtgood_ramp)
i4_ramp=calcFancyDelta(i3_good_ramp,'RampMinRatio_Late','RampMinRatio_LateDelta')
lateSummary_ramp=summarizeByMutationGeneral(i4_ramp,'RampMinRatio_Late','LateInclude',roundUnits=4)
lateSummaryDelta_ramp=summarizeByMutationGeneral(i4_ramp,'RampMinRatio_LateDelta','LateInclude',roundUnits=4)[,c(1,3,4)]
lateSummary_ramp=merge(lateSummary_ramp,lateSummaryDelta_ramp,all=TRUE)
lateSummary_ramp=lateSummary_ramp[lateSummary_ramp$numCellsRampMinRatio_Late>=5,]

lateSummary2 <- lateSummary_50 %>% merge(lateSummary_200,all=TRUE) %>% merge(lateSummary_ramp,all=TRUE) %>%
  orderMutFrame(.,'Late200Ratio_LateMean','WT',decreasing=TRUE)
lateSummary3=lateSummary2
lateSummary4=getMutDescriptionAndFilter(lateSummary3,reorderCols=FALSE)
lateSummary4=lateSummary4[,c('Mutation','MutDescription','Late50Ratio_LateMean','Late50Ratio_LateSE','Late50Ratio_LateDeltaMean','Late50Ratio_LateDeltaSE','numCellsLate50Ratio_Late','Late200Ratio_LateMean','Late200Ratio_LateSE','Late200Ratio_LateDeltaMean','Late200Ratio_LateDeltaSE','numCellsLate200Ratio_Late','RampMinRatio_LateMean','RampMinRatio_LateSE','RampMinRatio_LateDeltaMean','RampMinRatio_LateDeltaSE','numCellsRampMinRatio_Late')]
lateSummary5=lateSummary4[with(lateSummary4,order(MutDescription,Late200Ratio_LateMean,decreasing=TRUE)),]
row.names(lateSummary5)=1:nrow(lateSummary5)
write.csv(lateSummary5,'CombinedAnalysis/SP5-95/SP5-95_late_all.csv',row.names=FALSE)

##### STEP 2H) Calculate -90 mV peak current using bulk inactivation number (no sqrt transformation) #####
#inactSummary5=read.csv('CombinedAnalysis/SP5-95/SP5-95_inact_all.csv',stringsAsFactors =FALSE)
inactSummary5=read.csv('~/Dropbox/Andrew-Matthew/Syncropatch/CombinedAnalysis/SP5-95/SP5-95_inact_all.csv')
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_allCells_act.csv',stringsAsFactors=FALSE)
i=filterWeirdPlates(i,FALSE)
empty=read.csv('CombinedAnalysis/SP5-95/SP5-95_missingCells.csv',stringsAsFactors = FALSE)
j=bulk90adjust(i,inactSummary5)
peakSummary=summarizeByMutationGeneral(j,'PeakDensityNormBulk90','PeakInclude',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'PeakDensityNormBulk90Mean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(j,'CombinedAnalysis/SP5-95/SP5-95_inact_all_withBulk90.csv',row.names=FALSE)
write.csv(peakSummary3,'CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormBulk90Data.csv',row.names=FALSE)

##### STEP 2I) Calculate -90 mV peak current using bulk inactivation number (with sqrt transformation) #####
i=read.csv('CombinedAnalysis/SP5-95/SP5-95_inact_all_withBulk90.csv',stringsAsFactors=FALSE)
empty=read.csv('CombinedAnalysis/SP5-95/SP5-95_missingCells.csv',stringsAsFactors = FALSE)
j = sqrt_transform_90bulk(i,Negative2Zero = TRUE)
peakSummary=summarizeByMutationGeneral(j,'PeakDensityNormBulk90SQRT','PeakInclude',missingFrame=empty)
peakSummary2=orderMutFrame(peakSummary,'PeakDensityNormBulk90SQRTMean','WT',decreasing=TRUE)
peakSummary3=getMutDescriptionAndFilter(peakSummary2)
write.csv(peakSummary3,'CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormBulk90SQRTData.csv',row.names=FALSE)

########### STEP 3: Make merged file of summary data ##########

peakSummaryX=read.csv('CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormData.csv',stringsAsFactors = FALSE)
peakSummarySQRTX=read.csv('CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormSQRTData.csv',stringsAsFactors = FALSE)
peakSummarySQRTX=peakSummarySQRTX[,c(1,5,6)]
actSummaryX=read.csv('CombinedAnalysis/SP5-95/SP5-95_act_all.csv',stringsAsFactors = FALSE)
actSummaryX=actSummaryX[,-2]
InacttimeSummaryX=read.csv('CombinedAnalysis/SP5-95/SP5-95_inacttime_all.csv',stringsAsFactors=FALSE)
InacttimeSummaryX=InacttimeSummaryX[,-2]
inactSummaryX=read.csv('CombinedAnalysis/SP5-95/SP5-95_inact_all.csv',stringsAsFactors = FALSE)
inactSummaryX=inactSummaryX[,-2]
RFISummaryX=read.csv('CombinedAnalysis/SP5-95/SP5-95_RFI_all.csv',stringsAsFactors =FALSE)
RFISummaryX=RFISummaryX[,-2]
lateSummaryX=read.csv('CombinedAnalysis/SP5-95/SP5-95_late_all.csv',stringsAsFactors=FALSE)
lateSummaryX=lateSummaryX[,-2]
peakSummary90BulkX=read.csv('CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormBulk90Data.csv',stringsAsFactors = FALSE)
peakSummary90BulkX=peakSummary90BulkX[,-2]
peakSummary90BulkSQRTX=read.csv('CombinedAnalysis/SP5-95/SP5-95_NaIV_peakNormBulk90SQRTData.csv',stringsAsFactors = FALSE)
peakSummary90BulkSQRTX=peakSummary90BulkSQRTX[,c(1,5,6)]

ab1=merge(peakSummaryX, peakSummarySQRTX, all = TRUE)
ab2=merge(ab1,peakSummary90BulkX,all=TRUE)
ab3=merge(ab2,peakSummary90BulkSQRTX,all=TRUE)
ab4=merge(ab3,actSummaryX,all=TRUE)
ab5=merge(ab4,InacttimeSummaryX,all=TRUE)
ab6=merge(ab5,inactSummaryX,all=TRUE)
ab7=merge(ab6,RFISummaryX,all=TRUE)
ab8=merge(ab7,lateSummaryX,all=TRUE)
ab8[,c('Mutation','MutDescription')]
allEP=ab8[!is.na(ab8$Mutation),]

# Add in Z-scores based on benign distribution from VCCRI/VUMC collaboration [NOTE: not final!]
benign_mean <- 94.2
benign_sd <- 10.1 
allEP$z_score_120peak <- (allEP$PeakDensityNormSQRTMean-benign_mean) / benign_sd
allEP$z_se_120peak <- (allEP$PeakDensityNormSQRTSE-benign_sd) / benign_sd
allEP$ACMG_120peak <- convert_z_acmg(allEP$z_score_120peak)

benign_mean <- 95.50
benign_sd <- 17.48
allEP$z_score_90peak <- (allEP$PeakDensityNormBulk90SQRTMean-benign_mean) / benign_sd
allEP$z_se_90peak <- (allEP$PeakDensityNormBulk90SQRTSE-benign_sd) / benign_sd
allEP$ACMG_90peak <- convert_z_acmg(allEP$z_score_90peak)

summary(as.factor(allEP$MutDescription))
write.csv(allEP,'CombinedAnalysis/SP5-95/SP5-95_allEP_v3.csv',row.names=FALSE)

# Make summary of allEP table including ACMG criteria for subsets of variants
allEP <- read.csv('CombinedAnalysis/SP5-95/SP5-95_allEP_v3.csv')

# Walsh cohort+benign+path controls
Walsh_variants <- read.csv('CombinedAnalysis/VariantGroups/Walsh_Variants.csv', header = FALSE)$V1
Benign_variants <- read.csv('CombinedAnalysis/VariantGroups/Benign_variants.csv', header = FALSE)$V1
Pathogenic_variants <- read.csv('CombinedAnalysis/VariantGroups/Pathogenic_variants.csv', header = FALSE)$V1
WalshBP_variants=unique(c(Walsh_variants,Benign_variants,Pathogenic_variants))
allEP_subset <- allEP[allEP$Mutation %in% WalshBP_variants,]
write.csv(allEP_subset, 'CombinedAnalysis/SP5-95/SP5-95_allEP_WalshBP_final.csv',row.names=FALSE)
sum(allEP_subset$z_score_90peak<= -4) #42 strong
