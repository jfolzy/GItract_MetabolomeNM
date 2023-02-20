#######################################
# Ontologies , VIPs
#######################################
rm(list=ls())
source("0_Configuration.R")

#Load Data
MetOntology <- readRDS(MetOntology_path)
SampleMetaData <- readRDS(SampleMetaData_path)
FoodWithinHours <- readRDS(paste0(CleanDataPath,"FoodWithin6Hours.RDS"))
MetDFvip <- readRDS(FiltLogNormMetDF_C_path)
MetDF <- readRDS(FiltLogNormMetDF_C_path)

#All Capsules 
SubsetRows <- which(SampleMetaData$`SampleType`== "Capsule")

MetDF <- MetDF[SubsetRows,]
SampleMetaData <- SampleMetaData[SubsetRows,]
FoodWithinHours <- FoodWithinHours[SubsetRows,]
SampleMetaData <- cbind(SampleMetaData,FoodWithinHours)

### Get VIP order for metabolites with highest VIP scores when separated by cluster

#Impute NAs with 1/10 of minimum value
TenthOfMinimumValues <- 0
for (i in 1:ncol(MetDFvip)){
  TenthOfMinimumValues[i] <- 0.1*min(MetDFvip[which(!is.na(MetDFvip[,i])),i])
}
impute.filt.data <- MetDFvip
for (i in 1:ncol(MetDFvip)){
  impute.filt.data[which(is.na(MetDFvip[,i])),i] <- TenthOfMinimumValues[i]
}

#Center and Scale
scaled.log.filt.data <- scale(impute.filt.data)
df2plot <- as.data.frame(scaled.log.filt.data)

#Set variable of Sample Type to Group plot 
######################
for (i in c(3,12)){

SampleType <- SampleMetaData[,i]

#MetDFvip.PLSDA <- opls(df2plot,SampleType)
VIPscores <- getVipVn(opls(df2plot,SampleType))
VIPscoresDF <- as.data.frame(VIPscores)

VIPorderedPlot <- as.data.frame(cbind(rownames(VIPscoresDF)[order(VIPscoresDF$VIPscores,decreasing = TRUE)], VIPscoresDF[order(VIPscoresDF$VIPscores, decreasing = TRUE),]))

if(i==3){
SubjectVIPs <- VIPorderedPlot
}else{
CapsuleVIPs <- VIPorderedPlot
}

}

##Filling in ontologies ##############################
AllOntologies <- "NA"
for (i in 1:ncol(MetDF)){
OntRow <-  which(colnames(MetDF)[i] == MetOntology[,1])
  AllOntologies[i] <- as.character(MetOntology[OntRow,6]) #13 for chemrich classes
}
SubjectVIPs$MetaboliteOntology <- "NA"
  for (i in 1:nrow(SubjectVIPs)){
    OntRow <-  which(SubjectVIPs[i,1] == MetOntology[,1])
    SubjectVIPs$MetaboliteOntology[i] <- as.character(MetOntology[OntRow,6])
  }
CapsuleVIPs$MetaboliteOntology <- "NA"
for (i in 1:nrow(CapsuleVIPs)){
  OntRow <-  which(CapsuleVIPs[i,1] == MetOntology[,1])
  CapsuleVIPs$MetaboliteOntology[i] <- as.character(MetOntology[OntRow,6])
}

OntologyCount <- as.data.frame(unique(AllOntologies))
OntologyCount$Count <- "NA"
OntologyCount$TopSubjectVIPs <- "NA"
OntologyCount$TopCapsuleTypeVIPs <- "NA"
for (i in 1:nrow(OntologyCount)){
  OntologyCount$Count[i] <-  length(which(OntologyCount[i,1] == AllOntologies))
  OntologyCount$TopSubjectVIPs[i] <- length(which(OntologyCount[i,1] == SubjectVIPs$MetaboliteOntology[1:100]))
  OntologyCount$TopCapsuleTypeVIPs[i] <- length(which(OntologyCount[i,1] == CapsuleVIPs$MetaboliteOntology[1:100]))
  }

###### Filtering ontologies and plot
RsToRmv <- NA
for (i in 1:nrow(OntologyCount)){
  if(as.numeric(OntologyCount[i,2]) < 8 | sum(as.numeric(OntologyCount[i,3:4])) == 0){RsToRmv <- c(RsToRmv, i)}
}
RsToRmv <- RsToRmv[-1]

zOtherCount <- c(sum(as.numeric(OntologyCount[RsToRmv,2])),sum(as.numeric(OntologyCount[RsToRmv,3])),sum(as.numeric(OntologyCount[RsToRmv,4])))
OntologyCount <- OntologyCount[-RsToRmv[-1],]
#OntologyCount[RsToRmv,1] <- "zOther"
OntologyCount <- rbind(OntologyCount,c("zOther",zOtherCount))

OntologyCount[,2] <- as.numeric(OntologyCount[,2])
OntologyCount[,3] <- as.numeric(OntologyCount[,3])
OntologyCount[,4] <- as.numeric(OntologyCount[,4])
colnames(OntologyCount) <- c("Chemical Class","All Metabolites","Subject VIPs","Capsule VIPs")

OntologyCount <- OntologyCount[,-2]

df2plot <- melt(OntologyCount)

#Set Colors
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cpal <- col_vector[sample(1:60, nrow(df2plot), replace=F)]

ggplot(df2plot, aes(fill=`Chemical Class`, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity", width = 0.75)+
#  geom_bar(stat="identity", width = 0.75)+
  scale_fill_manual(values = cpal)+
  labs(x="",
       y="Proportion of chemical class",
       title= "Metabolites important to differentiate subject and capsule types")  +
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.5),
    )


