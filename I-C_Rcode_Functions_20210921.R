#Trade-off between competition ability and invulnerability to predation in marine microbes – contrasting of protist grazing and viral lysis effects
#Jinny Wu Yang, Feng-Hsun Chang, Yi-Chun Yeh, An-Yi Tsai, Fuh-Kwo Shiah, Kuo-Ping Chiang, Gwo-Ching Gong, Chih-hao Hsieh
#jinnyang@umich.edu

##### Rarefaction ####
Rarefaction_Overall_Or_T12=function(table,SampleID,Overall){
  if (Overall==FALSE){  
    SampleID=SampleID[-1,]
    SamplePool="T12"
  }else{
    SampleID=SampleID
    SamplePool="T0AndT12"
  }
  
  #Select samples based on SampleID list
  for (j in 1:ncol(SampleID)){
    SS=SampleID[,j]
    TableEachSample=matrix(0,nrow(table),length(SS))
    for(k in 1:length(SS)){
      t=is.na(as.character(SS)[k])
      if (t==T) {
        TableEachSample[,k]=rep(0,nrow(table))
      }else{
        TableEachSample[,k]=table[,which(colnames(table)==as.character(SS)[k])]
      }
    }
    colnames(TableEachSample)=as.character(SS)
    rownames(TableEachSample)=rownames(table)
    
    ##remove unavailable sample with NA
    if(all(!is.na(colnames(TableEachSample)))==F){
      FinalEachTable=TableEachSample[,-which(is.na(colnames(TableEachSample)))] 
    }else{
      FinalEachTable=TableEachSample
    }
    
    ##Rareified based on the minimu number of reads in the selected sample pool.
    OTU = otu_table(FinalEachTable, taxa_are_rows = TRUE)
    raretable=rarefy_even_depth(OTU, sample.size = min(colSums(FinalEachTable)))
    
    write.table(raretable,file=paste("table_ASV",colnames(SampleID)[j],SamplePool,paste(min(colSums(FinalEachTable)),"txt",sep="."),sep="_"),quote=F)
  }
  
}


##### Match and select sample between flow cytometry enumeration and sequencing data ####
SelectFlowCounts_T12=function(tableD30,T12_30){
  DCruise_NotClean  = strsplit(colnames(tableD30), "Or2", fixed = TRUE)
  DCruise = sapply(strsplit(sapply(strsplit(colnames(tableD30), "Or2", fixed = TRUE), "[", 2), "St", fixed = TRUE), "[", 1)
  
  DStation_NotClean = strsplit(colnames(tableD30), "St0", fixed = TRUE)
  DStation=sapply(strsplit(sapply(DStation_NotClean, "[", 2), "R", fixed = TRUE),"[", 1)
  
  DDilution_NotYet = sapply(strsplit(colnames(tableD30), "D", fixed = TRUE),"[", 3)
  DDilution = as.numeric(DDilution_NotYet)*25
  
  DReplicate_NotClean = sapply(strsplit(colnames(tableD30), "R", fixed = TRUE),"[", 2)
  DReplicate = sapply(strsplit(DReplicate_NotClean, "D", fixed = TRUE),"[", 1)
 
   #Pick the right total bacterial abundance with given sampleID
  FlowCountTableEachSample_D30_T12_EachASV=numeric()
  for ( i in 1:length(DStation)){
    FlowCountTable = T12_30[which(T12_30$Station == DStation[i] 
                                  & T12_30$Cruise== DCruise [i]
                                  & T12_30$Replicate == DReplicate[i]
                                  & T12_30$Dilution == DDilution[i]),]
    
    FlowCountTableEachSample_D30_T12_EachASV[i]=FlowCountTable$TB_D 
  }
  names(FlowCountTableEachSample_D30_T12_EachASV)=colnames(tableD30)
  
  return(FlowCountTableEachSample_D30_T12_EachASV)
  
}
SelectFlowCounts_T0=function(tableT0,T0_30){
  DCruise_NotClean  = strsplit(colnames(tableT0), "Or2", fixed = TRUE)
  DCruise = sapply(strsplit(sapply(strsplit(colnames(tableT0), "Or2", fixed = TRUE), "[", 2), "St", fixed = TRUE), "[", 1)
  
  DStation_NotClean = strsplit(colnames(tableT0), "St0", fixed = TRUE)
  DStation = sapply(strsplit(sapply(DStation_NotClean, "[", 2), "", fixed = TRUE),"[", 1)
  
  FlowCountTableEachSample_D30_T0_EachASV=numeric()
  FlowCountTable = T0_30[which(T0_30$Station == DStation 
                               & T0_30$Cruise== DCruise),]
  
  FlowCountTableEachSample_D30_T0_EachASV=FlowCountTable$TB_D #Picking total bacterial count
  
  ##Naming colnames
  D=c("D1","D2","D3","D4")
  Name=numeric()
  for (i in 1:4){
    Name[i]=paste(colnames(tableT0),D[i],sep="_")
  }
  names(FlowCountTableEachSample_D30_T0_EachASV)=rep(Name,2)
  return(FlowCountTableEachSample_D30_T0_EachASV)
  
}


##### Function_TradeOffPlot ####
TradeOffPLot=function(x,y,Xlim,Ylim,XlimLabels,YlimLabels,XLAB,YLAB,label,StationCruise){
  plot(x, y, col="black", pch = 1, cex = 1, ylim = Ylim, xlim = Xlim, xlab="",xaxt='n',yaxt='n',ylab="")
  axis(side = 1, at = XlimLabels, labels = FALSE, tck = -0.03)
  axis(1, at = XlimLabels, line = -1, labels = XLAB, cex.axis = 0.5, tick = F)
  axis(side= 2,at = YlimLabels, labels = F,tck = -0.03)
  axis(side= 2,at = YlimLabels, line=-0.5,labels = YLAB, cex.axis = 0.6, tick = F)
  
  mtext(label, line=-1.25, cex=0.75, at=Xlim[1]+0.5)
  
  l=lm(y~x)
  if (summary(l)$coef[2,4]<=0.05){
    LTY=1
    lwd=2
  }else{
    LTY=2
    lwd=1
  }
  abline(l,lty=LTY,col="blue",lwd=lwd)
  
  Cor=cor(x,y)
  Corp=cor.test(x,y)$p.value
  if (Corp<=0.001){
    Corp="<0.001"
  }else{
    Corp=Corp
  }
  
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(r) == MYVALUE), list(MYVALUE = format(Cor, digits = 2)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(Corp, digits = 2)))[2]
  #legend('bottomleft', legend = rp, bty = 'n',cex=0.7)
  legend('bottomright', legend = rp, bty = 'n',cex=0.7)
  legend('topright',StationCruise, bty = 'n',cex=0.85)
}

#### Fit rank-normalize RAD to zipfbrot model and calculate RAD decay coefficient ####
zipfbrot.gamma=function(table){
  PV=numeric()
  PVName=numeric()
  for(i in 1:ncol(table)){
    mod = rad.zipfbrot(table[,i])
    PV[i] = mod$coefficients[2]
    PVName[i]=colnames(table)[i]
  }
  names(PV)=colnames(table)
  names(PV)=sapply(strsplit(sapply(strsplit(PVName,"F02Or2"), "[", 2) ,"R"), "[", 1)
  return(PV)
}

#### Diversity under viral effect ####
Virus_gamma=function(rl){
  S=c("Cr2014St01","Cr2014St09","Cr2057St01","Cr2057St09","Cr2109St01","Cr2109St09","Cr2164St01")
  Virus_gamma=list()
  for (j in 1:length(S)){
    Data=rl[which(rl$Sample==S[j]),]
    all=matrix(0,4,4)
    
    for (i in 1:4){
      D=Data[which(Data$Dilution==i),]
      DD=D[which(D$Treatment=="D"),]$gamma
      DD30=D[which(D$Treatment=="D30"),]$gamma
      
      NewD=c(DD30[1]-DD[1],DD30[1]-DD[2],DD30[2]-DD[1],DD30[2]-DD[2])
      all[i,]=NewD
      colnames(all)=c("r1","r2","r3","r4")
      rownames(all)=c("D1","D2","D3","D4")
    }
    
    Virus_gamma[[j]]=all
  }
  names(Virus_gamma)=S
  mVirus_gamma=melt(Virus_gamma)
  colnames(mVirus_gamma)=c("Dilution","Replicates","Gamma","Sample")
  return(mVirus_gamma)
}
Virus_evenness=function(rl){
  S=c("Cr2014St01","Cr2014St09","Cr2057St01","Cr2057St09","Cr2109St01","Cr2109St09","Cr2164St01")
  Virus_evenness=list()
  for (j in 1:length(S)){
    Data=rl[which(rl$Sample==S[j]),]
    all=matrix(0,4,4)
    
    for (i in 1:4){
      D=Data[which(Data$Dilution==i),]
      DD=D[which(D$Treatment=="D"),]$Evenness
      DD30=D[which(D$Treatment=="D30"),]$Evenness
      
      NewD=c(DD30[1]-DD[1],DD30[1]-DD[2],DD30[2]-DD[1],DD30[2]-DD[2])
      all[i,]=NewD
      colnames(all)=c("r1","r2","r3","r4")
      rownames(all)=c("D1","D2","D3","D4")
    }
    
    Virus_evenness[[j]]=all
  }
  names(Virus_evenness)=S
  mVirus_evenness=melt(Virus_evenness)
  colnames(mVirus_evenness)=c("Dilution","Replicates","Evenness","Sample")
  return(mVirus_evenness)
}
Virus_richness=function(rl){
  S=c("Cr2014St01","Cr2014St09","Cr2057St01","Cr2057St09","Cr2109St01","Cr2109St09","Cr2164St01")
  Virus_richness=list()
  for (j in 1:length(S)){
    Data=rl[which(rl$Sample==S[j]),]
    all=matrix(0,4,4)
    
    for (i in 1:4){
      D=Data[which(Data$Dilution==i),]
      DD=D[which(D$Treatment=="D"),]$Richness
      DD30=D[which(D$Treatment=="D30"),]$Richness
      
      NewD=c(DD30[1]-DD[1],DD30[1]-DD[2],DD30[2]-DD[1],DD30[2]-DD[2])
      all[i,]=NewD
      colnames(all)=c("r1","r2","r3","r4")
      rownames(all)=c("D1","D2","D3","D4")
    }
    
    Virus_richness[[j]]=all
  }
  names(Virus_richness)=S
  mVirus_richness=melt(Virus_richness)
  colnames(mVirus_richness)=c("Dilution","Replicates","Richness","Sample")
  return(mVirus_richness)
}

#### LMM ####
LMMEst=function(x,y,r){
  library(nlme) 
  ms=lme(y~x,random= ~1|r, na.action=na.omit)
  LMM=c(summary(ms)[]$tTable[2,],fixef(ms))
  names(LMM)=c("Value" ,"Std. Error","DF","t_value","p_value","Fix_intercept","Fix_slope")  
  LMM=as.data.frame(t(LMM))
  return(LMM)
}

#### Plot Predation-diversity relationship ####
DiversityPlot_Combine_Protists=function(l,x,Tr,Ylim,LMMTable){
  COL=c("red","tomato2","blue","midnightblue","springgreen4","chartreuse3","maroon","tan4","black")
  PCH=c(1,20,5,18,2,17,8,13)
  
  plot(as.numeric(l[[1]][which(l[[1]]$Treatment==Tr),]$Dilution),
       as.numeric(t(l[[1]][which(l[[1]]$Treatment==Tr),][x])), 
       ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
       cex.axis=0.7, ylab="", pch=PCH[1], col=COL[1],cex=0.8)
  
  ll=lm(as.numeric(t(l[[1]][which(l[[1]]$Treatment==Tr),][x]))~
          as.numeric(l[[1]][which(l[[1]]$Treatment==Tr),]$Dilution))
  
  if (summary(ll)$coeff[2,4]<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  
  abline(ll, lty=LTY,col=COL[1])
  
  for (i in 2:7){
    points(as.numeric(l[[i]][which(l[[i]]$Treatment==Tr),]$Dilution),
           as.numeric(t(l[[i]][which(l[[i]]$Treatment==Tr),][x])), 
           ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
           cex.axis=0.7, ylab="", pch=PCH[i], col=COL[i],cex=0.8)
    
    ll=lm(as.numeric(t(l[[i]][which(l[[i]]$Treatment==Tr),][x]))~
            as.numeric(l[[i]][which(l[[i]]$Treatment==Tr),]$Dilution))
    
    if (summary(ll)$coeff[2,4]<=0.05){
      LTY=1
    }else{
      LTY=2
    }
    abline(ll, lty=LTY,col=COL[i])
  }
  
  axis(side = 1,cex.axis=0.8, at = c(1,2,3,4), labels = c("25%","50%","75%", "100%"), tck = -0.05)
  
  if (LMMTable$p_value<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  if (LMMTable$p_value<=0.001){
    p="<0.001"
  }else{
    p=LMMTable$p_value
  }
  rp = vector('expression',1)
  rp[1] = substitute(expression(italic(LMM-p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(p, digits = 2)))[2]
  

  abline(c(LMMTable$Fix_intercept,LMMTable$Fix_slope), lty=LTY,col="black", lwd=3)
  legend('topright', legend = rp, bty = 'n',cex=0.8)
}
DiversityPlot_Combine_Virus=function(V,Ylim,LMMTable){
  COL=c("red","tomato2","blue","midnightblue","springgreen4","chartreuse3","maroon","tan4","black")
  PCH=c(1,20,5,18,2,17,8,13)
  
  plot(as.numeric(V[[1]]$Dilution),
       as.numeric(t(V[[1]][,3])), 
       ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
       cex.axis=0.7, ylab="", pch=PCH[1], col=COL[1],cex=0.8)
  
  ll=lm(as.numeric(t(V[[1]][,3]))~
          as.numeric(V[[1]]$Dilution))
  
  if (summary(ll)$coeff[2,4]<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  
  abline(ll, lty=LTY,col=COL[1])
  
  for (i in 2:7){
    points(as.numeric(V[[i]]$Dilution),
           as.numeric(t(V[[i]][,3])), 
           ylim = Ylim, xlim=c(1,4), xlab="",xaxt='n', 
           cex.axis=0.7, ylab="", pch=PCH[i], col=COL[i],cex=0.8)
    
    ll=lm(as.numeric(t(V[[i]][,3]))~
            as.numeric(V[[i]]$Dilution))
    
    if (summary(ll)$coeff[2,4]<=0.05){
      LTY=1
    }else{
      LTY=2
    }
    abline(ll, lty=LTY,col=COL[i])
  }
  
  axis(side = 1,cex.axis=0.8, at = c(1,2,3,4), labels = c("25%","50%","75%", "100%"), tck = -0.05)
  
  if (LMMTable$p_value<=0.001){
    p="<0.001"
  }else{
    p=LMMTable$p_value
  }
  if (LMMTable$p_value<=0.05){
    LTY=1
  }else{
    LTY=2
  }
  rp = vector('expression',1)
  rp[1] = substitute(expression(italic(LMM-p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(p, digits = 2)))[2]
  
  
  abline(c(LMMTable$Fix_intercept,LMMTable$Fix_slope), lty=LTY,col="black", lwd=3)
  legend('topright', legend = rp, bty = 'n',cex=0.8)
}

#### Permutation and Shapiro-Wilk test ####
Permutation_shapiro_test=function(Slope,IGR){
  library("lme4")
  Observed=summary(lm(Slope~IGR))$coef[2,1] 
  Null_Slope_IGR=matrix(0,1000,3)
  colnames(Null_Slope_IGR)=c("Intercept","Estimates","p-value")
  for (j in 1:1000){
    nSlope=Slope[sample(1:length(Slope))] 
    l=lm(nSlope ~ IGR)
    Null_Slope_IGR[j,1]=summary(l)$coef[1,1] #
    Null_Slope_IGR[j,2]=summary(l)$coef[2,1] #Estimates
    Null_Slope_IGR[j,3]=summary(l)$coef[2,4] # 
  }
  Final=numeric()
  Final[1]=shapiro.test(Null_Slope_IGR[,2])$statistic
  Final[2]=shapiro.test(Null_Slope_IGR[,2])$p.value
  z = (Observed-mean(Null_Slope_IGR[,2]))/sd(Null_Slope_IGR[,2])
  Final[3]=z 
  p = 2*pnorm(-abs(z)) # two-tail p-value
  Final[4]=p
  names(Final)=c("shapiro","shapiro_pvalue","zscore","p-value")
  return(Final)
}
