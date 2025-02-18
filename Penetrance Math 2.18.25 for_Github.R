# Code to analyze SCN5A-BrS penetrance using method of McGurk et al 2023 AJHG (PMID 37652022)

##### CORE FUNCTION TO CALCULATE PENETRANCE  #####
# Note: the following function is adopted from Github repository cited in PMID 37652022
penetrance <- function(x_D, n_D, x_AgD, n_AgD, x_A, n_A){
  set.seed(28061971)
  digits <- 6   
  alpha <- 0.05  
  p_D <- x_D / n_D
  p_AgD <- x_AgD / n_AgD
  p_A <- x_A / n_A
  d <- 0.5   
  p_D_wod <- (x_D + d) / (n_D + d)
  p_AgD_wod <- (x_AgD + d) / (n_AgD + d)
  p_A_wod <- (x_A + d) / (n_A + d)
  log_AR <- log(p_D_wod * p_AgD_wod / p_A_wod) + 1/2 * ((1 / p_A_wod) * (1 - p_A_wod) / (n_A + d) - 
                                                          (1 / p_D_wod) * (1 - p_D_wod) / (n_D + d) - 
                                                          (1 / p_AgD_wod) * (1 - p_AgD_wod) / (n_AgD + d))
  Var_log_AR <- (1 / p_D_wod) * (1 - p_D_wod) / (n_D + d) + 
    (1 / p_AgD_wod) * (1 - p_AgD_wod) / (n_AgD + d) + 
    (1 / p_A_wod) * (1 - p_A_wod) / (n_A + d)
  log_LCI <- log_AR - qnorm(1 - alpha / 2) * sqrt(Var_log_AR)
  log_UCI <- log_AR + qnorm(1 - alpha / 2) * sqrt(Var_log_AR)
  penetrance <- pmin(1,pmax(0,exp(log_AR)))
  log_lCI <- pmin(1,pmax(0,exp(log_LCI)))
  log_uCI <- pmin(1,pmax(0,exp(log_UCI)))
  my_list <- list("penetrance" = penetrance, "lci" = log_lCI, "uci" = log_uCI)
  return(my_list)
}

##### OTHER FUNCTIONS #####
calcPenetrance=function(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN){
  missenseCount=length(missenseList)
  Cases_total=c(walsh_pLOF_cases)
  gnomAD.Count=c(gnomAD_AC_pLOF)
  gnomAD_AF=c(gnomAD_AF_pLOF)
  for(binCount in 1:missenseCount){
    currMissenseDF=missenseList[[binCount]]
    bin_cases=sum(currMissenseDF$walsh_total)
    bin_gnomAD_AC=sum(currMissenseDF$gnomAD_v4.1_AC)
    bin_gnomAD_AF=sum(currMissenseDF$gnomAD.v4.1_AF)
    Cases_total=c(Cases_total,bin_cases)
    gnomAD.Count=c(gnomAD.Count,bin_gnomAD_AC)
    gnomAD_AF=c(gnomAD_AF,bin_gnomAD_AF)
  }
  Bin.Number=0:missenseCount
  Case_AF=Cases_total/walsh_AN
  df_brs=data.frame(Bin.Number=Bin.Number,BinDescription=BinDescription,Cases_total=Cases_total,gnomAD.Count=gnomAD.Count,Case_AF=Case_AF,gnomAD_AF=gnomAD_AF)
  df_brs$gnomAD_AN=df_brs$gnomAD.Count/df_brs$gnomAD_AF
  binResults=list()
  for(binCount in 0:missenseCount){
    bin=penetrance(wilde_brugadaCases,wilde_brugadaN,df_brs[binCount+1,'Cases_total'],walsh_AN,df_brs[binCount+1,'gnomAD.Count'],df_brs[binCount+1,'gnomAD_AN']) 
    binResults[[binCount + 1]] <- bin
  }
  # Aggregate results into 1 table
  penetranceTable=data.frame()
  print(length(binResults))
  for(i in 1:(missenseCount+1)){
    currBin=binResults[[i]]
    penetranceTable[i,'bin']=i-1
    penetranceTable[i,'penetrance']=currBin[[1]]
    penetranceTable[i,'lci']=currBin[[2]]
    penetranceTable[i,'uci']=currBin[[3]]
  }
  penetranceTable$description=df_brs$BinDescription
  return(list(penetranceTable,df_brs))
}
plotPenetrance=function(penetranceTable,df_brs,pdfNamePen,includePLof,main){
  numBins=nrow(df_brs)
  # Doing math for right side of axis tick marks (odds ratio)
  ccc=1:5000
  orrr=100*(804760-ccc)/(ccc*(3335-100))
  ccc_val1=ccc[which.min(abs(orrr - 1))]
  ccc_val10=ccc[which.min(abs(orrr - 10))]
  ccc_val50=ccc[which.min(abs(orrr - 50))]
  ccc_val100=ccc[which.min(abs(orrr - 100))]
  ccc_val200=ccc[which.min(abs(orrr - 200))]
  ccc_val300=ccc[which.min(abs(orrr - 300))]
  ccc_val400=ccc[which.min(abs(orrr - 400))]
  ccc_val500=ccc[which.min(abs(orrr - 500))]
  ccc_val600=ccc[which.min(abs(orrr - 600))]
  ccc_val700=ccc[which.min(abs(orrr - 700))]
  pen1=(1/2000)*(804760/3335)*100/ccc_val1
  pen10=(1/2000)*(804760/3335)*100/ccc_val10
  pen50=(1/2000)*(804760/3335)*100/ccc_val50
  pen100=(1/2000)*(804760/3335)*100/ccc_val100
  pen200=(1/2000)*(804760/3335)*100/ccc_val200
  pen300=(1/2000)*(804760/3335)*100/ccc_val300
  pen400=(1/2000)*(804760/3335)*100/ccc_val400
  pen500=(1/2000)*(804760/3335)*100/ccc_val500
  pen600=(1/2000)*(804760/3335)*100/ccc_val600
  pen700=(1/2000)*(804760/3335)*100/ccc_val700
  
  if(!includePLof){
    numBins=numBins-1
    penetranceTable=penetranceTable[-1,]
    df_brs=df_brs[-1,]
  }
  
  # Plot penetrance estimates by Z-score
  pdf(pdfNamePen)
  plot(1:numBins,rep(0,numBins),type='n',xlab='Peak current bin',ylab='BrS Penetrance',ylim=c(0,0.4),main=main)
  print(penetranceTable)
  for(i in 1:numBins){
    value=penetranceTable[i,'penetrance']
    points(i,value,pch=20,cex=2)
    lci=penetranceTable[i,'lci']
    uci=penetranceTable[i,'uci']
    segments(i,lci,i,uci,lwd=1.5)
    segments(i-.1,lci,i+.1,lci,lwd=1.5)
    segments(i-.1,uci,i+.1,uci,lwd=1.5)
  }
  # Tick marks for OR
  abline(h=0)
  segments(numBins+0.13,pen1,numBins+0.2,pen1)
  segments(numBins+0.13,pen10,numBins+0.2,pen10)
  segments(numBins+0.13,pen50,numBins+0.2,pen50)
  segments(numBins+0.13,pen100,numBins+0.2,pen100)
  segments(numBins+0.13,pen200,numBins+0.2,pen200)
  segments(numBins+0.13,pen300,numBins+0.2,pen300)
  segments(numBins+0.13,pen400,numBins+0.2,pen400)
  segments(numBins+0.13,pen500,numBins+0.2,pen500)
  segments(numBins+0.13,pen600,numBins+0.2,pen600)
  segments(numBins+0.13,pen700,numBins+0.2,pen700)

  # OR conversion calculation
  # Average allele number in gnomADv4.1:
  gnomAD_AN=missense$gnomAD_v4.1_AC/missense$gnomAD.v4.1.AF
  gnomAD_AN=gnomAD_AN[is.finite(gnomAD_AN)]
  gnomAD_AN_avg=mean(gnomAD_AN) #1,609,520
  gnomAD_AN_over2=gnomAD_AN_avg/2 #804760
  
  # Plot penetrance estimates by odds ratios
  # We calculated the odds ratio (OR) associated with each variant class according to the formula (a/b)/(c/d), 
  # where a = BrS cases with variant, b = BrS cases without variant, c = gnomAD controls with variant, and 
  # d = gnomAD controls without variant. Because the allele number varied for different variants in gnomAD, 
  # the average allele number was calculated over all relevant variant types (missense, frameshift, nonsense, 
  # and splice site) and divided by 2 to obtain a count of sequenced gnomAD participants to use in OR calculations, 
  # following a previously published approach.
  df_brs$OR_a=df_brs$Cases_total
  df_brs$OR_b=walsh_sampleSize-df_brs$Cases_total
  df_brs$OR_c=df_brs$gnomAD.Count
  df_brs$OR_d=df_brs$gnomAD_AN/2 - df_brs$gnomAD.Count
  df_brs$OR=(df_brs$OR_a/df_brs$OR_b)/(df_brs$OR_c/df_brs$OR_d)
  barplot(df_brs$OR~df_brs$Bin.Number,xlab='Bin number',ylab='Odds Ratio')
  dev.off()
  
  print('Odds ratios:')
  print(df_brs$OR)
}
convert1to3=function(aa){
  one_digit_abb <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y','*','-')
  three_letter_abb <- c('Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile', 'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr', 'Val', 'Trp', 'Tyr','Ter','del')
  convertDF=data.frame(oneLetter=one_digit_abb,threeLetters=three_letter_abb)
  aa2=convertDF[convertDF$oneLetter==aa,'threeLetters']
  return(aa2)
}




##### MAIN PENETRANCE CODE -- BrS penetrance math from summary table #####

##### CONSTANT BITS -- RUN ONCE #####
# BRUGADA
# Mizusawa and Wilde 2012: Metaanalysis of 26 studies--prevalance of Type I BrS pattern is 0.05% (1 in 2000).
# https://www.ahajournals.org/doi/full/10.1161/CIRCEP.111.964577
# n=388,237. Don't have exact numerator.
# 388237*.0005=194.1-->round to 194.

walsh_sampleSize=3335
walsh_AN = 3335*2
wilde_brugadaCases=194
wilde_brugadaN=388237

# Calculate penetrance for predicted loss of function variants
walsh_pLOF_cases = 130

# gnomad v4.1.0 for MANE Select transcript ENST00000423572.7 / NM_000335.5 -- click on pLOF only
gnomAD_AC_pLOF = 232
gnomAD_AF_pLOF = 0.00014423
gnomAD_AN_pLOF = gnomAD_AC_pLOF/gnomAD_AF_pLOF #1608542
bin_LOF=penetrance(wilde_brugadaCases,wilde_brugadaN,walsh_pLOF_cases,walsh_AN,gnomAD_AC_pLOF,gnomAD_AN_pLOF) # not as bad as missense variants! similar to what we saw earlier

##### CALCULATING PENETRANCE BY CURRENT DENSITY #####
# functional_1 = CD_90
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4,] #246/252
missense_bin1=missense[missense$functional_1< (-4),]
missense_bin2=missense[missense$functional_1< (-3) & missense$functional_1>= (-4),]
missense_bin3=missense[missense$functional_1< (-2) & missense$functional_1>= (-3),]
missense_bin4=missense[missense$functional_1< (-1) & missense$functional_1>= (-2),]
missense_bin5=missense[missense$functional_1>= (-1),]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5)
sapply(missenseList,nrow)
#100  25  21  27  73
BinDescription=c('pLOF','Z_under_-4','Z_under_-3','Z_under_-2','Z_under_-1','Z_over_-1')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_functional_v7.pdf',FALSE,'Functional data (Z-score)')
#bin  penetrance         lci         uci description
#2   1 0.149139226 0.116751724 0.190511180  Z_under_-4
#3   2 0.038266957 0.027233871 0.053769807  Z_under_-3
#4   3 0.014475857 0.009698042 0.021607499  Z_under_-2
#5   4 0.016053915 0.011435551 0.022537451  Z_under_-1
#6   5 0.006688034 0.005280231 0.008471181   Z_over_-1
#[1] "Odds ratios:"
#[1] 318.44942  77.77435  29.21927  32.55114  13.80646

##### CALCULATING PENETRANCE BY CURRENT DENSITY (no boundary) #####
# Fx_evidence_noBoundary = CD_90 class (no boundary)
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4,] #246/252
summary(as.factor(missense$Fx_evidence_noBoundary))
#BS3_moderate BS3_supporting   PS3_moderate    PS3_strong1    PS3_strong2    PS3_strong3 PS3_supporting 
#73             27             25             23             35             42             21 
missense_bin1=missense[missense$Fx_evidence_noBoundary=="PS3_strong3",]
missense_bin2=missense[missense$Fx_evidence_noBoundary=="PS3_strong2",]
missense_bin3=missense[missense$Fx_evidence_noBoundary=="PS3_strong1",]
missense_bin4=missense[missense$Fx_evidence_noBoundary=="PS3_moderate",]
missense_bin5=missense[missense$Fx_evidence_noBoundary=="PS3_supporting",]
missense_bin6=missense[missense$Fx_evidence_noBoundary=="BS3_supporting",]
missense_bin7=missense[missense$Fx_evidence_noBoundary=="BS3_moderate",]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5,missense_bin6,missense_bin7)
sapply(missenseList,nrow)
#42 35 23 25 21 27 73
BinDescription=c('pLOF','Z_under_-6','Z_under_-5','Z_under_-4','Z_under_-3','Z_under_-2','Z_under_-1','Z_over_-1')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_functional_noBoundary_v7.pdf',FALSE,'Functional data (ACMG)')
#bin  penetrance         lci         uci description
#2   1 0.245213702 0.159350485 0.377342810  Z_under_-6
#3   2 0.190786004 0.133326205 0.273009339  Z_under_-5
#4   3 0.077865756 0.053413275 0.113512530  Z_under_-4
#5   4 0.038266957 0.027233871 0.053769807  Z_under_-3
#6   5 0.014475857 0.009698042 0.021607499  Z_under_-2
#7   6 0.016053915 0.011435551 0.022537451  Z_under_-1
#8   7 0.006688034 0.005280231 0.008471181   Z_over_-1
#[1] "Odds ratios:"
#[1] 501.11177 392.37431 158.22618  77.77435  29.21927  32.55114  13.80646

##### CALCULATING PENETRANCE BY FINAL CLASSIFICATION--AM CALIBRATED POST #####
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4,] #246/252
missense_bin1=missense[missense$AM_calibrated_post=='P',]
missense_bin2=missense[missense$AM_calibrated_post=='LP',]
missense_bin3=missense[missense$AM_calibrated_post=='VUS',]
missense_bin4=missense[missense$AM_calibrated_post=='LB',]
missense_bin5=missense[missense$AM_calibrated_post=='B',]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5)
sapply(missenseList,nrow)
#15 111 113   6   1
BinDescription=c('pLOF','P','LP','VUS','LB','B')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_AM_Post_Class_v7.pdf',FALSE,'Final classification')
#bin   penetrance          lci        uci description
#2   1 0.1809711903 0.1242729775 0.26353735           P
#3   2 0.0871287680 0.0691319138 0.10981068          LP
#4   3 0.0099150087 0.0080418795 0.01222443         VUS
#5   4 0.0023892400 0.0011945175 0.00477889          LB
#6   5 0.0009287204 0.0001847567 0.00466842           B
#[1] "Odds ratios:"
#[1] 140.573813 370.820903 184.953547  20.841489   4.781855   1.729458



##### CALCULATING PENETRANCE BY ALPHAMISSENSE ONLY #####
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4 & !missense$Computational_AM=='indel',] #238/252
summary(as.factor(missense$Computational_AM))
#BP4_-3   BP4_moderate BP4_supporting  none       PP3_+3   PP3_moderate     PP3_strong  PP3_supporting
#8             17             15       64          29             37             53        23
missense_bin1=missense[missense$Computational_AM=='PP3_strong',]
missense_bin2=missense[missense$Computational_AM=='PP3_+3',]
missense_bin3=missense[missense$Computational_AM=='PP3_moderate',]
missense_bin4=missense[missense$Computational_AM=='PP3_supporting',]
missense_bin5=missense[missense$Computational_AM=='none',]
missense_bin6=missense[missense$Computational_AM=='BP4_supporting',]
missense_bin7=missense[missense$Computational_AM=='BP4_moderate',]
missense_bin8=missense[missense$Computational_AM=='BP4_-3',]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5,missense_bin6,missense_bin7,missense_bin8)
sapply(missenseList,nrow)
# 53 29 37 23 56 15 17  8
BinDescription=c('pLOF','PP3_strong','PP3_+3','PP3_moderate','PP3_supporting','none','BP4_supporting','BP4_moderate','BP4_-3')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_AMonly_v7.pdf',FALSE, 'AlphaMissense')
# bin  penetrance         lci         uci    description
# 2   1 0.082150185 0.061257727 0.110168190     PP3_strong
# 3   2 0.067359002 0.049375156 0.091893081         PP3_+3
# 4   3 0.023332202 0.017509141 0.031091853   PP3_moderate
# 5   4 0.023171266 0.016676886 0.032194715 PP3_supporting
# 6   5 0.013316492 0.010292900 0.017228279           none
# 7   6 0.007285236 0.004566897 0.011621602 BP4_supporting
# 8   7 0.005071264 0.003120395 0.008241813   BP4_moderate
# 9   8 0.003203379 0.001656665 0.006194154         BP4_-3
#[1] "Odds ratios:"
#[1] 169.241946 137.960108  47.699893  47.073135  27.360128  14.656832  10.194548   6.416076

##### CALCULATING PENETRANCE BY REVEL ONLY #####
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4 & !missense$Computational_AM=='indel',] #238/252
summary(as.factor(missense$Computational_REVEL))
#BP4_supporting           none         PP3_+3   PP3_moderate     PP3_strong PP3_supporting 
#2                       44             32             32            108             20 
missense_bin1=missense[missense$Computational_REVEL=='PP3_strong',]
missense_bin2=missense[missense$Computational_REVEL=='PP3_+3',]
missense_bin3=missense[missense$Computational_REVEL=='PP3_moderate',]
missense_bin4=missense[missense$Computational_REVEL=='PP3_supporting',]
missense_bin5=missense[missense$Computational_REVEL=='none',]
missense_bin6=missense[missense$Computational_REVEL=='BP4_supporting',]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5,missense_bin6)
sapply(missenseList,nrow)
# [1] 108  32  32  20  44   2
BinDescription=c('pLOF','PP3_strong','PP3_+3','PP3_moderate','PP3_supporting','none','BP4_supporting')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_REVELonly_v7.pdf',FALSE,'REVEL')
#bin  penetrance         lci         uci    description
#2   1 0.059965276 0.048617519 0.073961699     PP3_strong
#3   2 0.012497655 0.009173894 0.017025635         PP3_+3
#4   3 0.014331875 0.010673062 0.019244958   PP3_moderate
#5   4 0.010189002 0.006825355 0.015210308 PP3_supporting
#6   5 0.006662699 0.004895520 0.009067792           none
#7   6 0.117979536 0.020328202 0.684722182 BP4_supporting
#[1] "Odds ratios:"
#[1] 128.85943  25.41208  29.21325  20.55751  13.53207 236.24350


##### CALCULATING PENETRANCE BY HOTSPOT #####
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4,] #246/252
summary(as.factor(missense$hotspot))
#none   PM1_moderate PM1_supporting 
#69            138             39 
missense_bin1=missense[missense$hotspot=='PM1_moderate',]
missense_bin2=missense[missense$hotspot=='PM1_supporting',]
missense_bin3=missense[missense$hotspot=='none',]
missenseList=list(missense_bin1,missense_bin2,missense_bin3)
sapply(missenseList,nrow)
#[1] 138  39  69
BinDescription=c('pLOF','PM1_moderate','PM1_supporting','none')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_hotspot_v7.pdf',FALSE,'Hotspot data')
# bin  penetrance         lci         uci    description
# 2   1 0.076106754 0.061903230 0.093569237   PM1_moderate
# 3   2 0.056110605 0.040012505 0.078685401 PM1_supporting
# 4   3 0.006242454 0.004967948 0.007843929           none
#[1] "Odds ratios:"
# [1] 165.41997 114.29053  12.92876

# Try for sex-specific?



### SUPP FIGURE ###
##### CALCULATING PENETRANCE BY FINAL CLASSIFICATION--REVEL CALIBRATED POST #####
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4,] #246/252
missense_bin1=missense[missense$REVEL_calibrated_post=='P',]
missense_bin2=missense[missense$REVEL_calibrated_post=='LP',]
missense_bin3=missense[missense$REVEL_calibrated_post=='VUS',]
missense_bin4=missense[missense$REVEL_calibrated_post=='LB',]
missense_bin5=missense[missense$REVEL_calibrated_post=='B',]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5)
sapply(missenseList,nrow)
#  19 118 108   1   0
BinDescription=c('pLOF','P','LP','VUS','LB','B')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_REVEL_Post_Class_v7.pdf',FALSE,'Final classification')
#bin   penetrance          lci         uci description
#2   1 0.1575808184 0.1126272772 0.220476912           P
#3   2 0.0720179782 0.0575182234 0.090172973          LP
#4   3 0.0079984421 0.0064669054 0.009892688         VUS
#5   4 0.0009287204 0.0001847567 0.004668420          LB
#6   5          NaN          NaN         NaN           B
#[1] "Odds ratios:"
#[1] 324.473604 153.011378  16.747876   1.729458        NaN


##### CALCULATING PENETRANCE BY FINAL CLASSIFICATION--REVEL SUPPORTING POST #####
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4,] #246/252
missense_bin1=missense[missense$REVEL_supporting_post=='P',]
missense_bin2=missense[missense$REVEL_supporting_post=='LP',]
missense_bin3=missense[missense$REVEL_supporting_post=='VUS',]
missense_bin4=missense[missense$REVEL_supporting_post=='LB',]
missense_bin5=missense[missense$REVEL_supporting_post=='B',]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5)
sapply(missenseList,nrow)
#  12 102 130   2   0
BinDescription=c('pLOF','P','LP','VUS','LB','B')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_REVEL_Supporting_Post_Class_v7.pdf',FALSE,'Final classification')
#bin  penetrance          lci         uci description
#2   1 0.176128782 0.1199606790 0.258595968           P
#3   2 0.140989329 0.1089252425 0.182492052          LP
#4   3 0.009394373 0.0076868553 0.011481189         VUS
#5   4 0.001931205 0.0006649092 0.005609117          LB
#6   5         NaN          NaN         NaN           B
#[1] "Odds ratios:"
#[1] 360.346520 297.720415  19.911542   3.823909        NaN

##### CALCULATING PENETRANCE BY FINAL CLASSIFICATION--ALPHAMISSENSE SUPPORTING POST #####
gnomad4.1=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/gnomad_v4.1_pLOF_forR.csv',stringsAsFactors=FALSE)
missense=read.csv('~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/functionalData_penetrance_forR_v7.csv')
missense=missense[missense$gnomAD.v4.1_AF<1e-4,] #246/252
missense_bin1=missense[missense$AM_supporting_post=='P',]
missense_bin2=missense[missense$AM_supporting_post=='LP',]
missense_bin3=missense[missense$AM_supporting_post=='VUS',]
missense_bin4=missense[missense$AM_supporting_post=='LB',]
missense_bin5=missense[missense$AM_supporting_post=='B',]
missenseList=list(missense_bin1,missense_bin2,missense_bin3,missense_bin4,missense_bin5)
sapply(missenseList,nrow)
#  11 101 132   1   1
BinDescription=c('pLOF','P','LP','VUS','LB','B')
result=calcPenetrance(missenseList,BinDescription,walsh_pLOF_cases,gnomAD_AC_pLOF,walsh_AN)
penetranceTable=result[[1]]
df_brs=result[[2]]
plotPenetrance(penetranceTable,df_brs,'~/Dropbox/Andrew-Matthew/SCN5A SP ACMG VUS (Bezzina)/WarePenetrance/BrS_Penetrance_Mcgurk_AM_Supporting_Post_Class_v7.pdf',FALSE,'Final classification')
#bin   penetrance          lci        uci description
#2   1 0.2094130544 0.1375439723 0.31883496           P
#3   2 0.1459909411 0.1128847808 0.18880627          LP
#4   3 0.0094415466 0.0077294052 0.01153294         VUS
#5   4 0.0049462766 0.0013787954 0.01774422          LB
#6   5 0.0009287204 0.0001847567 0.00466842           B
#[1] "Odds ratios:"
#[1] 427.545682 308.771115  20.023964   9.677787   1.729458