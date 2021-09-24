source("I-C_Rcode_Functions_20210921.R")

######### Part 1: Prepare ASV tables for analysis ########
#### 1. Prepare non-rarefied table for each experiment ####
library("gdata")
library("phyloseq")
RawTable=read.table("table_ASV.txt",header=T, row.names=1) #Input raw ASV table
SampleID=read.xls("SampleName.xlsx", sheet=1,header=T) #Input sample name list

##Remove singletons from ASV table
table_Single=RawTable
table_Single[table_Single>=1]=1 ##Remove rows that all zero
table=RawTable[which(rowSums(table_Single)!=1),]

write.table(table,"NoSingletons_ASVTable.txt",quote = F)

table=read.table("NoSingletons_ASVTable.txt",header=T, row.names=1) #Input raw ASV table
MinRA_all=numeric()
for (i in 1:ncol(SampleID)){
  SS=SampleID[,i]
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
    NonRarefiedEachTable=TableEachSample[,-which(is.na(colnames(TableEachSample)))] 
  }else{
    NonRarefiedEachTable=TableEachSample
  }
  NonRarefiedEachTableNoRow0=NonRarefiedEachTable[rowSums(NonRarefiedEachTable)!=0,] #Remove ASV with all 0
  
  ##Output non-rarefied ASV table
  write.table(NonRarefiedEachTableNoRow0,file=paste("NonRarfied_table_ASV",paste(colnames(SampleID)[i],"txt",sep="."),sep="_"),quote=F)

  ##Calculate smallest relative abundance for each experiment 
  MinRA=min(rowSums(NonRarefiedEachTableNoRow0/sum(NonRarefiedEachTableNoRow0))) 
  MinRA_all[i]=MinRA
  }
MinRA_all[i]=MinRA
names(MinRA_all)=colnames(SampleID)
write.csv(MinRA_all, "Smallest_NonRarefied_RA.csv") ##Output smallest relative abundance for each experiment

#### 2. Do rarefaction and calculate relative abundance for each experiment####
  SampleID=read.xls("SampleName.xlsx", sheet=1,header=T) #Input sample name list
  table=read.table("NoSingletons_ASVTable.txt",header=T, row.names=1) #Input raw ASV table
  
  Rarefaction_Overall_Or_T12(table,SampleID,TRUE) # Do rarefaction and output rarefied ASV tables with function "Rarefaction_Overall_Or_T12"
  # When Overall=FALSE, it remove T0 and will be rarefied based on the minimum number of reads from T12 sample; 
  # Otherwise, when Overall=TRUE, it will rarefied based on the minimum number of reads from a T0 and T12 pool.
  
#### 3. Replace 0 in rarefied-table with non-rarefied relative abundance or adding 10-8 if remain 0 before rarefaction####
  NonRarefiedTable_list=list.files(pattern="NonRarfied_table_ASV.*.txt") #Non-rarefied ASV tables list
  RarefiedTable_list=list.files(pattern="table_ASV.*._T0AndT12_.*.txt") #Input rarefied ASV tables list
  s = strsplit(RarefiedTable_list, "\\_|\\ | ")
  rarefaction=as.numeric(sapply(strsplit(sapply(s, "[", 5), ".", fixed = TRUE), "[", 1))
  CruiseStation=sapply(s, "[", 3)
  for (i in 1:length(RarefiedTable_list)){
    RarefiedTable=read.table(RarefiedTable_list[i],header=T,row.names = 1) #Input rarefied ASV tables
    RarefiedTable_RA=RarefiedTable/rarefaction[i] #Calculate relative abundance
    NonRarefiedTable=read.table(NonRarefiedTable_list[i],header=T,row.names = 1) #Input non-rarefied ASV tables
    NonRarefiedTable_RA=NonRarefiedTable/sum(NonRarefiedTable) #Calculate relative abundance
    
    FinalRA=RarefiedTable_RA
    for (j in 1:ncol(RarefiedTable_RA)){
      for (k in 1:nrow(RarefiedTable_RA)){
        # If relative abundance is 0 in rarefied ASV table, replaced it with the relative abundance from non-rarefied ASV table
        if (FinalRA[k,j]==0){
            FinalRA[k,j]=NonRarefiedTable_RA[which(rownames(NonRarefiedTable_RA)==rownames(RarefiedTable_RA)[k]),j]
          }
        # If relative abundance is still 0 in non-rarefied table, replace it with 0.00000001
        if (FinalRA[k,j]==0) {
          FinalRA[k,j]=0.00000001
        }  
      }
    }
    rownames(FinalRA)=rownames(RarefiedTable_RA)
    colnames(FinalRA)=colnames(RarefiedTable_RA)
    write.table(FinalRA,file=paste("FinalRA_T0AndT12",paste(colnames(SampleID)[i],"txt",sep="."),sep="_"),quote=F)
    }
  
  
  
  #################################### End of Part 1 ########################################
  
 #### Part 2: Growth-based competitiveness-invulnerability trade-off ####
  #### 1. Calculate the Competitiveness and invulnerability of each bacterial ASV ####
  library(gdata)
  library(lavaan)
  library(lme4)
  library(nlme)
  library(MASS)
  
  #Input Flow cytometery enumeration data
  T0=read.xls("FlowCytometry_Enumeration_Data.xlsx",sheet=1,header=T,row.names=1)
  T0_30=read.xls("FlowCytometry_Enumeration_Data.xlsx",sheet=3,header=T)
  T12=read.xls("FlowCytometry_Enumeration_Data.xlsx",sheet=2,header=T)
  T12_30=read.xls("FlowCytometry_Enumeration_Data.xlsx",sheet=4,header=T)
 
  #Input sample list 
  table_list = list.files(pattern="FinalRA_T0AndT12.*.txt") 
  s = strsplit(table_list, "\\_|\\ | ")
  CruiseStation=sapply(strsplit(sapply(s, "[", 3),".", fixed = TRUE),"[",1)
  
  for (i in 1:length(table_list)){
    DALL=read.table(table_list[i],header=T,row.names = 1) #Input 
  
  #Extract Sample at T0
  tableT0=as.matrix(DALL[,1])
  rownames(tableT0)=rownames(DALL)
  colnames(tableT0)=colnames(DALL)[1]
  
  #Split tables based on predation dilution treatment 
  DALL_NoT0=DALL[,-1] 
  tableD02=DALL_NoT0[,which(sapply(strsplit(colnames(DALL_NoT0), "S"), "[", 1)=="D")] 
  tableD30=DALL_NoT0[,which(sapply(strsplit(colnames(DALL_NoT0), "S"), "[", 1)=="D30")] 
  
  #Pick and organize flow cytometry data to match to the relative abundance data 
  FlowT12_D02 = SelectFlowCounts_T12(tableD02,T12)
  FlowT12_D30 = SelectFlowCounts_T12(tableD30,T12_30)
  
  FlowT0_D02 = SelectFlowCounts_T0(tableT0,T0)
  FlowT0_D30 = SelectFlowCounts_T0(tableT0,T0_30)
  
  HypothesisI_data=data.frame()
  HypothesisI_data=matrix(0,nrow(DALL),7)
  
  for (j in 1:nrow(DALL)){
    ##Calculate Net growth rate of each bacterial ASV 
    T12_ForEachASV_D02=FlowT12_D02*(tableD02[j,])
    T12_ForEachASV_D30=FlowT12_D30*(tableD30[j,])
    
    T0_ForEachASV_D02=t(FlowT0_D02)*(tableT0[j])
    rownames(T0_ForEachASV_D02)=rownames(tableT0)[j]
    T0_ForEachASV_D30=t(FlowT0_D30)*(tableT0[j])
    rownames(T0_ForEachASV_D30)=rownames(tableT0)[j]
    
    Net_D02=log((T12_ForEachASV_D02/T0_ForEachASV_D02)/12) #per hour
    Net_D30=log((T12_ForEachASV_D30/T0_ForEachASV_D30)/12) #per hour
    
    ##Prepare dilution factor based on sample name
    S_D02=strsplit(colnames(Net_D02), "D")
    SS_D02=sapply(S_D02, "[", 3)
    SS_D02=SS_D02[!is.na(SS_D02)]
    
    S_D30=strsplit(colnames(Net_D30), "D")
    SS_D30=sapply(S_D30, "[", 3)
    SS_D30=SS_D30[!is.na(SS_D30)]
    
    ##Linear regression of Net growth rate versus dilution factors
    l_D02=lm(as.numeric(Net_D02)~as.numeric(SS_D02))
    p_D02=summary(l_D02)$coefficients[2,4]
    l_D30=lm(as.numeric(Net_D30)~as.numeric(SS_D30))
    p_D30=summary(l_D30)$coefficients[2,4]
    
    #calculate the relative abundance for each ASV at T0
    T0_RA=as.numeric(tableT0[j,])
    
    # Take the mean of the whole community absolute abundance at T0 when without predation manipulation (100% predation) as the whole community abundance at T0
    Total_AbsA_D02=mean(FlowT0_D02[which(sapply(strsplit(names(FlowT0_D02),"_"), "[", 2)=="D4")])
    Total_AbsA_D30=mean(FlowT0_D30[which(sapply(strsplit(names(FlowT0_D30),"_"), "[", 2)=="D4")])
    T0_AbsA=T0_RA*mean(c(Total_AbsA_D02,Total_AbsA_D30))
    
    HypothesisI_data[j,]=c(T0_RA,
                           T0_AbsA,
                           coef(l_D02)[2],p_D02,
                           coef(l_D30)[2],p_D30,
                           as.numeric(coef(l_D30)[1]))
  }
  colnames(HypothesisI_data)=c("T0_RA",
                               "T0_AbsA",
                               "Slope_D02","Linear_Pvalue_D02",
                               "Slope_D30","Linear_Pvalue_D30",
                               "IGR")
  rownames(HypothesisI_data)=rownames(DALL)
  write.table(HypothesisI_data, quote=F, 
              file = paste("HypothesisI_variables", CruiseStation[i], 
                           paste(rarefaction[i],"txt",sep="."), sep="_"))
  
  }
  
  
  #### 2. Plot Competitiveness-Invulnerability trade-off ####
  HI_list = list.files(pattern="HypothesisI_variables.*.txt") 
  CruiseStation=sapply(strsplit(HI_list, "\\_|\\ | "), "[", 3)
  
  YLAB=c(T,F,F,F,T,F,F,F)
  
  par(mfrow = c(4,4)) # 2-by-2 grid of plots
  par(oma = c(4, 4, 2, 2)) # make room for the overall x and y axis titles
  par(mar = c(0, 0, 0, 0)) # make the plots be closer together
  
  ### A.Growth-based Competitiveness-Invulnerability trade-off ####
  ## Protists-caused mortality (Slope_D02) ~ Predation-free growth rate (IGR) ####
  Xlim=c(-20,20)
  XlimLabels=c(-20,-15,-10,-5,0,5,10,15,20)
  Ylim=c(-6,6)
  YlimLabels=c(-6,-4,-2,0,2,4,6)
  label=c("a","b","c","d","e","f","g")
  XLAB=c(F,F,F,T,F,F,F,F)
  
  for (i in 1:7){
    t=read.table(HI_list[i],header=T,row.names = 1)
    x=t$IGR
    y=t$Slope_D02
    TradeOffPLot(x,y,Xlim,Ylim,XlimLabels,YlimLabels,XLAB[i],YLAB[i],label[i],CruiseStation[i])
  }  
  # make a empty plot
  plot(x,y, pch = NA,xlab="",xaxt='n',yaxt='n',ylab="") 
  
  ## Viral-caused mortality (SlopeD30-Slope_D02) ~ Predation-free growth rate (IGR) ####
  label=c("h","i","j","k","l","m","n")
  XLAB=c(F,F,F,T,T,T,T,T)
  for (i in 1:7){
    t=read.table(HI_list[i],header=T,row.names = 1)
    x=t$IGR
    y=t$Slope_D30-t$Slope_D02
    TradeOffPLot(x,y,Xlim,Ylim,XlimLabels,YlimLabels,XLAB[i],YLAB[i],label[i],CruiseStation[i])
  }  
  
  plot(x,y, pch = NA,xlab="",xaxt='n',yaxt='n',ylab="") # make a empty plot
  
  mtext("Invulnerability to",side=2, outer=T,cex=0.8, line=2.5, at=0.5)
  mtext("Protist grazing",side=2, outer=T,cex=0.7, line=1.25, at=0.875)
  mtext("Protist grazing",side=2, outer=T,cex=0.7, line=1.25, at=0.625)
  mtext("Viral lysis",side=2, outer=T,cex=0.7, line=1.25, at=0.375)
  mtext("Viral lysis",side=2, outer=T,cex=0.7, line=1.25, at=0.125)
  mtext("Predation-free growth rate",side=1, outer=T,cex=0.8, line=1)
  

  ### B.Density-based Competitiveness-Invulnerability trade-off ####
  ## Protists-caused mortality (Slope_D02) ~ Relative abundance at T0 ####
  Xlim=c(-3,4)
  XlimLabels=c(-3,-2,-1,0,1,2,3,4)
  Ylim=c(-6,6)
  YlimLabels=c(-6,-4,-2,0,2,4,6)
  label=c("a","b","c","d","e","f","g")
  XLAB=c(F,F,F,T,F,F,F,F)
  for (i in 1:7){
    t=read.table(HI_list[i],header=T,row.names = 1)
    x=log10(t$T0_AbsA)
    y=t$Slope_D02
    TradeOffPLot(x,y,Xlim,Ylim,XlimLabels,YlimLabels,XLAB[i],YLAB[i],label[i],CruiseStation[i])
  }  
  # make a empty plot
  plot(x,y, pch = NA,xlab="",xaxt='n',yaxt='n',ylab="") 
  
  ## Viral-caused mortality ~ Relative abundance at T0 ####
  label=c("h","i","j","k","l","m","n")
  XLAB=c(F,F,F,T,T,T,T,T)
  for (i in 1:7){
    t=read.table(HI_list[i],header=T,row.names = 1)
    x=log10(t$T0_AbsA)
    y=t$Slope_D30-t$Slope_D02
    TradeOffPLot(x,y,Xlim,Ylim,XlimLabels,YlimLabels,XLAB[i],YLAB[i],label[i],CruiseStation[i])
  }  
  
  plot(x,y, pch = NA,xlab="",xaxt='n',yaxt='n',ylab="") # make a empty plot
  
  mtext("Invulnerability to",side=2, outer=T,cex=0.8, line=2.5, at=0.5)
  mtext("Protist grazing",side=2, outer=T,cex=0.7, line=1.25, at=0.875)
  mtext("Protist grazing",side=2, outer=T,cex=0.7, line=1.25, at=0.625)
  mtext("Viral lysis",side=2, outer=T,cex=0.7, line=1.25, at=0.375)
  mtext("Viral lysis",side=2, outer=T,cex=0.7, line=1.25, at=0.125)
  mtext("log10-formed initial abundance",side=1, outer=T,cex=0.8, line=1)
  #################################### End of Part 2 ########################################
  
  
  
  ########################### PART 3: Calculate bacterial diversity under predation effect ################################################
  #### 1. Rank-based normalization #### 
  # Rank-normalized ASV tables are for obtaining RAD and Evenness
  library("gdata")
  library("RADanalysis")
  
  SampleID=read.xls("SampleName.xlsx", sheet=1,header=T) #Input Sample list
  SampleID=SampleID[-1,] #remove sample at T0
  Sample=colnames(SampleID)
  table=read.table("NoSingletons_ASVTable.txt",header=T,row.names=1) #Input raw ASV table
  
  for (j in 1:length(Sample)){
    SS=SampleID[,which(colnames(SampleID)==Sample[j])]
    SS=SS[!is.na(SS)] # Remove NA
    
    #Select sample table based on the given list
    TableEachSample=matrix(0,nrow(table),length(SS))
    for(k in 1:length(SS)){
      TableEachSample[,k]=table[,which(colnames(table)==as.character(SS)[k])]
    }
    colnames(TableEachSample)=as.character(SS)
    rownames(TableEachSample)=rownames(table)
    
    ### Finding minimum ranking among samples (smallest number of ASV among samples)
    TTable=as.data.frame(TableEachSample)
    TTable[TTable>=1]=1 #Making ASV number of reads > 1 become 1
    
    ### Normalizes an abundance table to the desired number of ranks 
    nrads=RADnormalization_matrix(input = TableEachSample, max_rank = min(colSums(TTable)),average_over = 1000, sample_in_row = F)
    
    nrad_mat=t(nrads$norm_matrix)
    colnames(nrad_mat)=colnames(TableEachSample)
    write.table(nrad_mat,file=paste("RankFixed_table",Sample[j],paste(min(colSums(TTable)),"txt",sep="."),sep="_"),quote=F)
     }
  
  #### 2. Calculate RAD decay coefficient and Evenness from rank-normalized RAD ####
  library(vegan)
  library(reshape)
  library("lme4")
  library(nlme)
  library(viridis)
  ### Calculate bacterial RAD and Evenness under protists and protists-viral effect ####
  temp = list.files(pattern="RankFixed_table.*.txt") ## Input sample list with T0 and T12 sample pool
  for (i in 1:length(temp)){
    table=read.table(temp[i],header=T,row.names=1)
    PVName=colnames(table)
    ## Calculate RAD decay coefficient with rank with rank-normalized table
    gamma=zipfbrot.gamma(table)
    gamma=gamma*(-1) # covert to positive value as decay coefficient
    EachSampleName=sapply(strsplit(temp,"_"), "[", 3)[i]
    colEachSampleName=rep(EachSampleName,ncol(table))
    
    ## Calculate Evenness with rank-normalized table
    H = diversity(t(table))
    J = H/log(specnumber(t(table)))
    
    ##Prepare and organize data
    Treatment=sapply(strsplit(PVName,"S"), "[", 1)
    Replicates=sapply(strsplit(sapply(strsplit(PVName,"D"), "[", 2),"R"), "[", 2)
    Dilution=sapply(strsplit(sapply(strsplit(PVName,"R"), "[", 2),"D"), "[", 2)
    MaxRank=as.numeric(strsplit(sapply(strsplit(temp,"_"), "[", 3),".txt"))
    
    finalTable=data.frame(colEachSampleName, gamma,J,Treatment,Dilution,Replicates)
    colnames(finalTable)=c("Sample","gamma","Evenness","Treatment","Dilution","Replicates")
    
    write.table(finalTable,file=paste("RankFixed_RAD_Evenness",paste(EachSampleName,"txt",sep="."),sep="_"),quote=F)
  }
  
  ### Calculate bacterial RAD and Evenness under viral effect ####
  temp = list.files(path=".", pattern="RankFixed_RAD_Evenness*")
  l=list()
  for (i in 1:length(temp)){
    l[[i]]=read.table(temp[i],header=T,row.names=1) 
  }
  rl=rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]])
  
  write.table(Virus_gamma(rl),"VirususEffect_Gamma.txt", quote=F)
  write.table(Virus_evenness(rl),"VirususEffect_Evenness.txt", quote=F)
  
  #### 2. Calculate Richness from size-normalized ASV table ####
  ## A. rarefied ASV table with smallest number of reads at T12 for each experiment ##
  SampleID=read.xls("SampleName.xlsx", sheet=1,header=T) #Input sample name list
  table=read.table("NoSingletons_ASVTable.txt",header=T, row.names=1) #Input raw ASV table
  
  # When Overall=FALSE, it remove T0 and will be rarefied based on the minimum number of reads from T12 sample; 
  # Otherwise, when Overall=TRUE, it will rarefied based on the minimum number of reads from a T0 and T12 pool.
  Rarefaction_Overall_Or_T12(table,SampleID,FALSE) 
  
  ### Calculate bacterial richness under protists and protists-viral effect ####
  temp = list.files(pattern="table_ASV.*._T12_.*.txt") ## Input rarefied ASV table at T12
  for (i in 1:length(temp)){
    t=read.table(temp[i], header=T, row.names = 1)
    Obs_R=as.data.frame(t)
    Obs_R[Obs_R>=1]=1
    data2=Obs_R[complete.cases(Obs_R),]
    Richness=colSums(data2)
    Dilution=as.numeric(sapply(strsplit(names(Richness), "D", fixed = TRUE), "[", 3))
    Treatment=sapply(strsplit(names(Richness), "S", fixed = TRUE), "[", 1)
    Sample=rep(sapply(strsplit(temp[i], "_", fixed = TRUE), "[", 3),length(Treatment))
    final=data.frame(Richness,Dilution,Treatment,Sample)
    write.table(final,quote=F,paste("SizeFixed_Richness",sapply(strsplit(temp, "_", fixed = TRUE), "[", 3)[i],".txt",sep=""))
  }
  
  ### Calculate bacterial richness under viral effect ####
  temp = list.files(path=".", pattern="SizeFixed_Richness*")
  l=list()
  for (i in 1:length(temp)){
    l[[i]]=read.table(temp[i],header=T,row.names=1) 
  }
  rl=rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]])
  write.table(Virus_richness(rl),"VirususEffect_Richness.txt", quote=F)
  
  ########################### End of Part 3 ################################################
  
  
  ########################### PART 4: Linear-mixed model ################################################
  library(reshape2)
  #### A. Hypothesis I: Competitiveness ~ invulnerability relationship ####
  ##Organize and recombine the data
  temp = list.files(pattern="HypothesisI_variables.*.txt") 
  l=list()
  for (i in 1:length(temp)){
    table=read.table(temp[i], header=T,row.names=1)
    Sample=sapply(strsplit(temp,"_"), "[", 3)[i]
    Station=sapply(strsplit(Sample,"St"), "[", 2)
    l[[i]]=data.frame(table,rep(Sample,nrow(table)),rep(Station,nrow(table)))
    colnames(l[[i]])=c(colnames(table),"Sample","Station")
    names(l)[i]=Sample
  }
  rl=rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]])
  
  #Do LMM with function "LM_Estimates"
  LMM_Estimates=rbind(LMMEst(rl$IGR,rl$Slope_D02,rl$Sample),
                      LMMEst(rl$IGR,rl$Slope_D30-rl$Slope_D02,rl$Sample),
                      LMMEst(log10(rl$T0_AbsA),rl$Slope_D02,rl$Sample),
                      LMMEst(log10(rl$T0_AbsA),rl$Slope_D30-rl$Slope_D02,rl$Sample))
  
  colnames(LMM_Estimates)=c("Value" ,"Std. Error","DF","t-value","p-value","Fix_intercept","Fix_slope")     
  rownames(LMM_Estimates)=c("Growth_Protists","Growth_Viruses",
                            "Density_Protists","Density_Viruses")
  write.csv(LMM_Estimates,"HypothesisI_LMM.csv", quote=F) # Output LMM result for Hypothesis I
  
  
  #### B. Hypothesis II: Diversity ~ Predation dilution factors relationship ####
  #Prepare and recombine data
  temp = list.files(path=".", pattern="RankFixed_RAD_Evenness_*")
  l=list()
  for (i in 1:length(temp)){
    l[[i]]=read.table(temp[i],header=T,row.names=1) 
  }
  rl=rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]])
  
  DD=rl[which(rl$Treatment=="D"),]
  DD30=rl[which(rl$Treatment=="D30"),]
  DD$Dilution=as.numeric(DD$Dilution)
  
  #Prepare data of bacterial richness response to viral lysis effect
  rtemp = list.files(path=".", pattern="SizeFixed_Richness_*")
  ll=list()
  for (i in 1:length(rtemp)){
    ll[[i]]=read.table(rtemp[i],header=T,row.names=1) 
  }
  rrl=rbind(ll[[1]],ll[[2]],ll[[3]],ll[[4]],ll[[5]],ll[[6]],ll[[7]])
  rDD=rrl[which(rrl$Treatment=="D"),]
  rDD30=rrl[which(rrl$Treatment=="D30"),]
  rDD$Dilution=as.numeric(rDD$Dilution)
  
  VirususEffect_Richness=read.table("VirususEffect_Richness.txt",header=T, row.names = 1)
  
  #Do LMM with function "LM_Estimates"
  LMM_Estimates=rbind(LMMEst(DD30$Dilution,DD30$gamma,DD30$Sample),
                      LMMEst(DD$Dilution,DD$gamma,DD$Sample),
                      LMMEst(as.numeric(Virus_gamma(rl)$Dilution),Virus_gamma(rl)$Gamma,as.factor(Virus_gamma(rl)$Sample)),
                      LMMEst(DD30$Dilution,DD30$Evenness,DD30$Sample),
                      LMMEst(DD$Dilution,DD$Evenness,DD$Sample),
                      LMMEst(as.numeric(Virus_evenness(rl)$Dilution),Virus_evenness(rl)$Evenness,as.factor(Virus_evenness(rl)$Sample)),
                      LMMEst(rDD30$Dilution,rDD30$Richness,rDD30$Sample),
                      LMMEst(rDD$Dilution,rDD$Richness,rDD$Sample),
                      LMMEst(as.numeric(VirususEffect_Richness$Dilution),VirususEffect_Richness$Richness,VirususEffect_Richness$Sample)
  )
  rownames(LMM_Estimates)=c("Combine_gamma","Protists_gamma","Virus_gamma",
                            "Combine_evenness","Protists_evenness","Virus_evenness",
                            "Combine_richness","Protists_richness","Virus_richness")
  write.csv(LMM_Estimates,"HypothesisII_LMM.csv", quote=F)#Output LMM result for Hypothesis II

  ########################### End of Part 4 ################################################
  
  
  
  ########################### PART 5: Plot Predation ~ Diversity relationship ################################################
  par(mfrow = c(1, 3)) # 1-by-3 grid of plots
  par(oma = c(3, 3, 2, 2)) # make room axis titles
  par(mar = c(1, 2, 0, 0)) # adjust width between plots 
  
  StationCruise=c("2014AprSt1","2014AprSt9","2014OctSt1","2014OctSt9",
                  "2015JulSt1","2015JulSt9","2016MaySt1")
  temp = list.files(path=".", pattern="RankFixed_RAD_Evenness_*")
  l=list()
  for (i in 1:length(temp)){
    l[[i]]=read.table(temp[i],header=T,row.names=1) 
  }
  LMMTable=read.csv("HypothesisII_LMM.csv",header=T,row.names=1)

  #### A. Plot predation ~ RAD decay coefficient relationship ####
  Ylim=c(0.5,3)
  x=2
  Tr="D30"
  DiversityPlot_Combine_Protists(l,2,Tr,Ylim,LMMTable[1,])
  mtext("Protists-viruses-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("a", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  Tr="D"
  DiversityPlot_Combine_Protists(l,2,Tr,Ylim,LMMTable[2,])
  mtext("Protists-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("b", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  Ylim=c(-2.5,2)
  VirususEffect_Gamma=read.table("VirususEffect_Gamma.txt", header=T, row.names=1)
  V=split(VirususEffect_Gamma,VirususEffect_Gamma$Sample)
  DiversityPlot_Combine_Virus(V,Ylim,LMMTable[3,])
  mtext("Viruses-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("c", outer=F,side=3, line=-1, at=1, cex=0.8)
  mtext("Predation dilution factors", outer=T,side=1, line=1, cex=0.8)
  mtext("RAD decay coefficient", outer=T, side=2, line=0, cex=0.8)
  
  COL=c("red","tomato2","blue","midnightblue","springgreen4","chartreuse3","maroon","tan4","black")
  PCH=c(1,20,5,18,2,17,8,13)
  
  legend(0.9,-1.2, 
         xpd = "NA",
         legend = StationCruise, 
         pch=PCH, col=COL, cex=0.5)
  
  #### B. Plot Predation ~ Evenness relationship ####
  Ylim=c(0.2,1)
  x=3
  Tr="D30"
  DiversityPlot_Combine_Protists(l,x,Tr,Ylim,LMMTable[4,])
  mtext("Protists-viruses-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("a", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  Tr="D"
  DiversityPlot_Combine_Protists(l,x,Tr,Ylim,LMMTable[5,])
  mtext("Protists-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("b", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  Ylim=c(-0.5,0.5)
  VirususEffect_Evenness=read.table("VirususEffect_Evenness.txt", header=T, row.names=1)
  V=split(VirususEffect_Evenness,VirususEffect_Evenness$Sample)
  DiversityPlot_Combine_Virus(V,Ylim,LMMTable[6,])
  mtext("Viruses-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("c", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  mtext("Predation dilution factors", outer=T,side=1, line=1, cex=0.8)
  mtext("Evenness", outer=T, side=2, line=0, cex=0.8)
  
  legend(-4.23,-0.215, 
         xpd = "NA",
         legend = StationCruise, 
         pch=PCH, col=COL, cex=0.5)
  
  ### C. Plot Predation ~ Richness relationship ####
  temp = list.files(path=".", pattern="SizeFixed_Richness_*")
  ll=list()
  for (i in 1:length(temp)){
    ll[[i]]=read.table(temp[i],header=T,row.names=1) 
  }
  
  Ylim=c(10,600)
  x=1
  Tr="D30"
  DiversityPlot_Combine_Protists(ll,x,Tr,Ylim,LMMTable[7,])
  mtext("Protists-viruses-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("a", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  Tr="D"
  DiversityPlot_Combine_Protists(ll,x,Tr,Ylim,LMMTable[8,])
  mtext("Protists-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("b", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  Ylim=c(-200,400)
  VirususEffect_Richness=read.table("VirususEffect_Richness.txt", header=T, row.names=1)
  V=split(VirususEffect_Richness,VirususEffect_Richness$Sample)
  DiversityPlot_Combine_Virus(V,Ylim,LMMTable[9,])
  mtext("Viruses-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("c", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  mtext("Predation dilution factors", outer=T,side=1, line=1, cex=0.8)
  mtext("Richness", outer=T, side=2, line=0, cex=0.8)
  
  legend(-4.23,-30, 
         xpd = "NA",
         legend = StationCruise, 
         pch=PCH, col=COL, cex=0.5)
  
  ########################### End of Part 5 ################################################
  
  
  
  ########################### PART 6: Bray-Curtis ################################################
  library(vegan)
  library("lme4")
  temp = list.files(pattern="table_ASV_.*._T0AndT12_.*.txt") 
  
  l=list()
  for (i in 1:length(temp)){
    t=read.table(temp[i], header=T, row.names = 1)
    x=vegdist(t(t), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) #Calculate Bray-Curtis dissimilarity between each T12 and T0
    tt=as.matrix(x)[,1]
    Bray_Curtis=tt[-1] #Remove T0
    Dilution=as.numeric(sapply(strsplit(names(Bray_Curtis), "D", fixed = TRUE), "[", 3))
    Treatment=sapply(strsplit(names(Bray_Curtis), "S", fixed = TRUE), "[", 1)
    Sample=rep(sapply(strsplit(temp[i], "_", fixed = TRUE), "[", 3),length(Bray_Curtis))
    final=data.frame(Bray_Curtis,Dilution,Treatment,Sample)
    
    l[[i]]=final
  }

  ##Analyze with LMM analysis
  rl=rbind(l[[1]],l[[2]],l[[3]],l[[4]],l[[5]],l[[6]],l[[7]])
  DD=rl[which(rl$Treatment=="D"),]
  DD30=rl[which(rl$Treatment=="D30"),]
  LMM_Protists=LMMEst(DD$Dilution,DD$Bray_Curtis,DD$Sample)
  LMM_Combine=LMMEst(DD30$Dilution,DD30$Bray_Curtis,DD30$Sample)
  
  #### Plotting ####
  COL=c("red","tomato2","blue","midnightblue","springgreen4","chartreuse3","maroon","tan4","black")
  PCH=c(1,20,5,18,2,17,8,13)
  
  par(mfrow = c(1, 2)) 
  par(oma = c(3, 3, 2, 2)) 
  par(mar = c(2, 2, 0, 0)) 
  
  Ylim=c(0.4,1)
  x=1
  Tr="D"
  DiversityPlot_Combine_Protists(l,x,Tr,Ylim,LMM_Protists)
  mtext("Protists-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("a", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  Tr="D30"
  DiversityPlot_Combine_Protists(l,x,Tr,Ylim,LMM_Combine)
  mtext("Protists-viruses-manipulated", outer=F,side=3, line=0, cex=0.7)
  mtext("b", outer=F,side=3, line=-1, at=1, cex=0.8)
  
  mtext("Predation dilution facotrs", outer=T,side=1, line=0, cex=0.8)
  mtext("Bray-Curtis disimilarity", outer=T, side=2, line=0, cex=0.8)

  Sample=sapply(strsplit(temp, "_", fixed = TRUE), "[", 3)
  legend(0.88,0.565, 
         xpd = "NA",
         legend = Sample, 
         pch=PCH, col=COL, cex=0.4)
  ########################### End of Part 6 ################################################
  
  
  ########################### PART 7: Permutation test ################################################
  temp = list.files(pattern="HypothesisI_variables.*.txt") #Input variable table list
  
  #### 1. Permutation test on Predation-free growth rate versus Protists-caused mortality ####
  Final=matrix(0,7,4)
  for (i in 1:7){
    table=read.table(temp[i],header=T,row.names = 1)
    Final[i,]=Permutation_shapiro_test(table$Slope_D02,table$IGR)
  }
  colnames(Final)=c("Shapiro-Wilk","Shapiro-Wilk_p-value","z-score","p-value")
  rownames(Final)=sapply(strsplit(temp, "_", fixed = TRUE), "[", 3)
  write.csv(Final,"Permutation_shapiro_Protist-IGR_RegresionSlope.csv")
  
  #### 2. Permutation test on Predation-free growth rate versus Viral-caused mortality ####
  Final=matrix(0,7,4)
  for (i in 1:7){
    table=read.table(temp[i],header=T,row.names = 1)
    Final[i,]=Permutation_shapiro_test(table$Slope_D30-table$Slope_D02,table$IGR)
  }
  colnames(Final)=c("Shapiro-Wilk","Shapiro-Wilk_p-value","z-score","p-value")
  rownames(Final)=sapply(strsplit(temp, "_", fixed = TRUE), "[", 3)
  write.csv(Final,"Permutation_shapiro_Virus-IGR_RegresionSlope.csv")
  ########################### End of Part 7 ################################################
  
  ######### Part 8: Examine potential bottle effect - Ranking shift between T0 and T12 without prediction manipulation ########
  library("gdata")
  SampleID=read.xls("SampleName.xlsx", sheet=1,header=T) #Input sample name list
  AllTable=read.table("NoSingletons_ASVTable.txt",header=T,row.names = 1)
  
  #Extract relative abundance at T0
  T0=rowMeans(AllTable[,which(sapply(strsplit(colnames(AllTable), "T"), "[", 2)=="0")] )
  
  Table_D4=AllTable[,which(sapply(strsplit(colnames(AllTable), "D"), "[", 3)=="4")] 
  D4=rowMeans(Table_D4) #take row means to represent relative abundance T T12
  
  T0D4=data.frame(T0,D4) #Combine two ranking
  
  #### Plotting: ASV ranking at T0 versus T12 ####
  plot(NULL, xlim=c(-100,2000), ylim=c(-100,2500),
       col="white", type="n",
       xlab="", ylab="",
       xaxt='n', yaxt='n')
  CEX=log10(rowSums(T0D4)+1)/2
  
  points(rank(-T0D4$T0),rank(-T0D4$D4),
         cex=CEX,
         ylim=c(1,max(rank(-T0D4$D4))),
         xlim=c(1,max(rank(-T0D4$T0))),
         col="black")
  rank=c(1,500,1000,1500,2000,2500)
  axis(side = 1, at = rank, labels = FALSE, tck = -0.03)
  axis(side = 2, at = rank, labels = FALSE, tck = -0.03)
  axis(1, at = rank, line = -1, labels = rank, cex.axis = 0.45, tick = F)
  axis(2, at = rank, line = -0.5, labels = rank, cex.axis = 0.45, tick = F)
  mtext("ASV ranking at T0", side=1, line=1, cex = 0.8)
  mtext("ASV ranking at T12", side= 2, line= 2, cex=0.8)
  mtext("(full predation)", side= 2, line= 1.3, cex=0.7)
  legend(-200,2800, pch=1, 
         pt.cex=min(CEX),
         cex=0.4,
         legend="Min:1.4e-07",
         horiz=TRUE, bty='n', xpd=TRUE, inset=c(0,-0.1))
  legend(500,2800, pch=1, 
         pt.cex=mean(CEX),
         cex=0.4,
         legend="Mean: 3.5e-04",
         horiz=TRUE, bty='n', xpd=TRUE, inset=c(0,-0.1))
  legend(1250,2870, pch=1, 
         pt.cex=max(CEX),
         cex=0.4,
         legend="",
         horiz=TRUE, bty='n', xpd=TRUE, inset=c(0,-0.1))
  mtext("Max: 1.2e-01", side=3, at=1620, cex=0.4)
 
  #### Calculate 10% ranking change between T0 and D4 ####
  #Calculate the ranking difference between T0 and T12
  T0D4=T0D4[rowSums(T0D4)!=0,]
  RankDif=abs(rank(-T0D4$T0)/max(rank(-T0D4$T0))
              -rank(-T0D4$D4)/max(rank(-T0D4$D4)))
  
  tt=data.frame(rowSums(T0D4),RankDif)
  tt_Top100=tt[which(rank(-tt$rowSums.T0D4.)<=100),] #Pick the top-100 ASVs (based on the sum of T0 and T12 samples)
  sum(tt_Top100$rowSums.T0D4.)/sum(tt$rowSums.T0D4.) #Calculate the % abundance of these top-100 ASVs
  length(which(tt_Top100$RankDif<=0.10)) #Calculate how many ASV ranking shift <10% in the top-100 

  ########################### End of Part 8 ################################################
  
  
  
