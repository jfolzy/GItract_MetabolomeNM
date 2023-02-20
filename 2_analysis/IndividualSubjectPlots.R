#######################################
# Individual Subject Clustering by Day/ Timepoint
#######################################
rm(list=ls())
source("0_Configuration.R")

#LoadData 
MetDF <- readRDS(logNormMetDF_C_path) # has NAs

MetOntology <- readRDS(MetOntology_path)

SampleMetaData <- readRDS(SampleMetaData_path)
FoodWithinHours <- readRDS(paste0(CleanDataPath,"FoodWithin12Hours.RDS"))

#All Capsules ###################################################
SubsetRows <- which(SampleMetaData$`SampleType`== "Capsule")
#####################

MetDF <- MetDF[SubsetRows,]
SampleMetaData <- SampleMetaData[SubsetRows,]
FoodWithinHours <- FoodWithinHours[SubsetRows,]
#bind food data for all capsules
SampleMetaData <- cbind(SampleMetaData,FoodWithinHours)

#Unique Subject/ Swallow Time
SubjectTimepoint <- paste(SampleMetaData$`Subject Number`, SampleMetaData$`Swallow Date Time`)
SampleMetaData <- cbind(SampleMetaData,SubjectTimepoint)

TimepointDF <- as.data.frame(unique(SubjectTimepoint))
TimepointDF[,2:3] <- NA
for (i in 1:nrow(TimepointDF)){
  TimepointDF[i,2] <- length(grep(TimepointDF[i,1],SubjectTimepoint))
  TProws <- grep(TimepointDF[i,1],SubjectTimepoint)
  TimepointDF[i,3] <- length(unique(SampleMetaData[TProws,6]))
  }
colnames(TimepointDF) <- c("Timepoint","CapsulesTaken","DifferentCapTypes")

#Sample meta data column for count of samples in timepoint
TimepointCount <- as.data.frame(SubjectTimepoint)
for (i in 1:nrow(TimepointCount)){
  TimepointCount[i,2] <- length(grep(TimepointCount[i,1],SubjectTimepoint))
}
TimepointCount[,2] <- as.numeric(TimepointCount[,2])

SampleMetaData$TimepointCount <- TimepointCount[,2]
colnames(SampleMetaData)[25] <- "Timepoint Count"

#Subsetting Sample Meta Data by count of greater than 2 capsules per timepoint

SamplesToRemove <- which(SampleMetaData$`Timepoint Count` < 2)

MetDF_S <- MetDF[-SamplesToRemove,]
SampleMetaData_S <- SampleMetaData[-SamplesToRemove,]

#For loop for making dendrograms #################

#Variables for storing significants for KW tests
SigRaw_CAPS <- NA
SigFDR_CAPS <- NA
SigRaw_TIMES <- NA
SigFDR_TIMES <- NA
SigNums <- as.data.frame(matrix(nrow=15,ncol=4))
colnames(SigNums) <- c("SigCapType_raw","SigCapType_FDR","SigTime_raw","SigTime_FDR")
SigMets <- as.data.frame(matrix(NA,nrow=2000,ncol=30))

for (J in 1:15){

OneSubjectRows <- which(SampleMetaData_S$`Subject Number` == paste0("Subject ",J))

OneSubject_MetDF <- MetDF_S[OneSubjectRows,]
OneSubject_MetaData <- SampleMetaData_S[OneSubjectRows,]

#Remove Metabolites with more than 50% missing values
OneSubject_MetDF[OneSubject_MetDF == 0] <- NA
ProportionMissingValues <- 0
for (i in 1:ncol(OneSubject_MetDF)){
  ProportionMissingValues[i] <- sum(is.na(OneSubject_MetDF[,i]))/length(rownames(OneSubject_MetDF))
}
OneSubject_MetDF <- OneSubject_MetDF[,which(ProportionMissingValues < 0.5)]
OneSubject_MetDF[is.na(OneSubject_MetDF)] <- 0

#Correlation based Clustering
cormat <- cor(t(OneSubject_MetDF), method = "spearman")
s <- cormat
diag(s) <- 0
s[is.na(s)] <- 0
hc <- hclust(as.dist(1-s), method="ward.D2")

hcd <- as.dendrogram(hc)

#Generate MetaData for timepoints

OneSubject_MetaData <- OneSubject_MetaData[order(OneSubject_MetaData$`Swallow Date Time`),]
OneSTimepointsDF <- as.data.frame(unique(OneSubject_MetaData$SubjectTimepoint))

for (i in 1:nrow(OneSTimepointsDF)){
OneSTimepointsDF[i,2] <- paste0("Timepoint ",i)
}

#Set Timepoint numeric column
OneSubject_MetaData$CapsuleSet <- NA
for (i in 1:nrow(OneSubject_MetaData)){
SET <- which(OneSubject_MetaData$SubjectTimepoint[i]==OneSTimepointsDF[,1])
OneSubject_MetaData$CapsuleSet[i] <- OneSTimepointsDF[SET,2]
  }

#Set Color Column for Sets
SetColorScheme <- as.data.frame(cbind(paste("Timepoint",1:12),brewer.pal(12,"Dark2")))
OneSubject_MetaData$SetColor <- NA
for (i in 1:nrow(OneSubject_MetaData)){
  CLR <- which(OneSubject_MetaData$CapsuleSet[i]== SetColorScheme[,1])
  OneSubject_MetaData$SetColor[i] <- SetColorScheme[CLR,2]
}

#Set Color Column for Capsules
CapColorScheme <- as.data.frame(cbind(paste(1:4),CapType))
OneSubject_MetaData$CapColor <- NA
for (i in 1:nrow(OneSubject_MetaData)){
  CLR <- which(OneSubject_MetaData$`Capsule Type`[i]== CapColorScheme[,1])
  OneSubject_MetaData$CapColor[i] <- CapColorScheme[CLR,2]
}

#Make data frame for capsule versus tiempoint comparison
CapVtime <- cbind(OneSubject_MetaData$CapsuleSet,OneSubject_MetaData$`Sample Type`,OneSubject_MetDF)
colnames(CapVtime)[1:2] <- c("Capsule Set","Capsule Type")

#Set to order of dendrogram
DendOrder <- as.data.frame(cbind(labels(hcd),1:length(labels(hcd))))
OneSubject_MetaData$DendOrder <- NA
for (i in 1:nrow(OneSubject_MetaData)){
  Ord <- which(OneSubject_MetaData$`Sample ID`[i]== DendOrder[,1])
  OneSubject_MetaData$DendOrder[i] <- as.numeric(DendOrder[Ord,2])
}

Dplot <- OneSubject_MetaData[order(OneSubject_MetaData$DendOrder),]

#Plot single subject
hcd %>%
  set("labels_col", Dplot$SetColor) %>% # change color
  set("labels", Dplot$CapsuleSet) %>% # Set Text Labels
  set("labels_cex", 0.6) %>% # Change size --- FUll plot 0.3 -- subset 0.6
  #set("branches_k_color", k = 4) %>% #Branch color and k
  set("leaves_pch", 15) %>%  # node point type
  set("leaves_cex", 1.75) %>%  # node point size --- 0.9 All Capsules -- subset 1.0
  set("branches_lwd", 1.7) %>%
  set("leaves_col", Dplot$SetColor)%>%
  plot(main = paste0("Dendogram of Samples from Subject ",J)) 
colored_bars(colors=Dplot$CapColor,sort_by_labels_order = FALSE)


####           CapType_V_Timepoint                   ##########################

KruskalWallisPs_CAPS <- 0
KruskalWallisPs_TIMES <- 0
for (i in 3:ncol(CapVtime)){
  KWtestCAPS <- kruskal.test(CapVtime[,i]~CapVtime$`Capsule Type`)
  KruskalWallisPs_CAPS[i-2] <- KWtestCAPS$p.value
  
  KWtestTIMES <- kruskal.test(CapVtime[,i]~CapVtime$`Capsule Set`)
  KruskalWallisPs_TIMES[i-2] <- KWtestTIMES$p.value
}

adjustedPs_CAPS <- p.adjust(KruskalWallisPs_CAPS, method= "BH")
adjustedPs_TIMES <- p.adjust(KruskalWallisPs_TIMES, method= "BH")

SigNums[J,1] <- length(which(KruskalWallisPs_CAPS < 0.05))
SigNums[J,2] <- length(which(adjustedPs_CAPS < 0.05))
SigNums[J,3] <- length(which(KruskalWallisPs_TIMES < 0.05))
SigNums[J,4] <- length(which(adjustedPs_TIMES < 0.05))

colnames(SigMets)[(J*2)-1] <- paste("Subject",J,"Sig_CAPS")
colnames(SigMets)[(J*2)] <- paste("Subject",J,"Sig_TIMES")

SigMets_CAPS <- colnames(CapVtime)[which(KruskalWallisPs_CAPS < 0.05)+2]
SigMets_TIMES <- colnames(CapVtime)[which(KruskalWallisPs_TIMES < 0.05)+2]

SigMets[1:length(SigMets_CAPS),((J*2)-1)] <- SigMets_CAPS
SigMets[1:length(SigMets_TIMES),(J*2)] <- SigMets_TIMES

}

#Number of Significants plots
RawSigNums <- SigNums[,c(1,3)]
SubjectNumber <- paste("Subject",1:15)
RawSigNums <- cbind(SubjectNumber,RawSigNums)
colnames(RawSigNums) <- c("Subject Number","Significant Metabolites between Capsule Types","Significant Metabolites between Capsule Sets")
df2plot <- melt(RawSigNums)
df3plot <- melt(RawSigNums)

ggplot(df3plot, aes(fill=`variable`, y=value, x=`Subject Number`)) + 
#  geom_bar(position="fill", stat="identity", width = 0.75)+ #For filled stacked bar chart
  geom_bar(stat="identity", width = 0.75)+
  labs(x="",
       y="Proportion of significant metbolites",
       title= "Metabolites significantly changed by capsule type or by capsule set")  +
  scale_fill_manual(values=c("darkcyan","goldenrod2"))+
  scale_x_discrete(name="Subject",
                   limits=c("Subject 1", "Subject 2", "Subject 3", "Subject 4",
                            "Subject 5", "Subject 6", "Subject 7", "Subject 8",
                            "Subject 9", "Subject 10", "Subject 11", "Subject 12",
                            "Subject 13", "Subject 14", "Subject 15"))+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5),
  )

# Metabolite Types Significantly changed by cap or time and by subject

SigMetOntologies <- SigMets
for (i in 1:ncol(SigMets)){
  for (j in which(!is.na(SigMets[,i]))){
  MetOntRow <- which(SigMets[j,i] == MetOntology[,1])
  SigMetOntologies[j,i] <- MetOntology[MetOntRow,5]
}
}

SigMetOnts_CAPS <- SigMetOntologies[,grep("CAPS",colnames(SigMetOntologies))]
SigMetOnts_TIMES <- SigMetOntologies[,grep("TIMES",colnames(SigMetOntologies))]

colnames(SigMetOnts_CAPS) <- c(paste("Subject",1:15))
colnames(SigMetOnts_TIMES) <- c(paste("Subject",1:15))

#Calculate all unique ontologies from sig mets
AllUniqueOnts <- NA
for (i in 1:30){
SubUni <- unique(SigMetOntologies[,i])
AllUniqueOnts <- unique(SubUni,AllUniqueOnts)
}
OntCols <- c("#201923","#fcff5d","#7dfc00","#0ec434","#228c68",
"#8ad8e8","#235b54","#29bdab","#3998f5","#37294f","#277da7","#3750db",
"#f22020","#991919","#ffcba5","#e68f66","#c56133","#96341c","#632819",
"#ffc413","#f47a22","#2f2aa0","#b732cc","#772b9d","#f07cab","#d30b94",
"#edeff3","#c3a5b4","#946aa2","#5d4c86")


# CAPS differences ontology Stacked Barchart ###########################
CAPS_CountDF <- as.data.frame(AllUniqueOnts[-19])
for (i in 1:15){
  for (j in 1:nrow(CAPS_CountDF)){
OntCunt <- length(grep(CAPS_CountDF[j,1],SigMetOnts_CAPS[,i]))
CAPS_CountDF[j,i+1] <- OntCunt
  }
}

colnames(CAPS_CountDF) <- c("Chemical Class",paste("Subject",1:15))

df2plot <- melt(CAPS_CountDF)

ggplot(df2plot, aes(fill=`Chemical Class`, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity", width = 0.75)+
  scale_fill_manual(values = OntCols)+
  labs(x="",
       y="Proportion of chemical class",
       title= "Metabolites Significantly Different Between Capsule Types")  +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
  )

# TIME stacked bar chart ontologies ############

TIMES_CountDF <- as.data.frame(AllUniqueOnts[-19])
for (i in 1:15){
  for (j in 1:nrow(TIMES_CountDF)){
    OntCunt <- length(grep(TIMES_CountDF[j,1],SigMetOnts_TIMES[,i]))
    TIMES_CountDF[j,i+1] <- OntCunt
  }
}

colnames(TIMES_CountDF) <- c("Chemical Class",paste("Subject",1:15))

df2plot <- melt(TIMES_CountDF)

ggplot(df2plot, aes(fill=`Chemical Class`, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity", width = 0.75)+
  scale_fill_manual(values = OntCols)+
  labs(x="",
       y="Proportion of chemical class",
       title= "Metabolites Significantly Different Between Capsule Sets (Sampling Timepoints)")  +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
  )

