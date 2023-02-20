#######################################
# Clustering - VIP - Heatmap
#######################################
rm(list=ls())
source("0_Configuration.R")

#non filtered data
SubjectColorScheme <- readRDS(SubjectColorScheme_path)
SampleTypeColorScheme <- readRDS(SampleTypeColorScheme_path)

MetDF <- readRDS(FiltLogNormMetDF_C_path)
MetDF[is.na(MetDF)] <- 0

SampleMetaData <- readRDS(SampleMetaData_path)
FoodWithinHours <- readRDS(paste0(CleanDataPath,"FoodWithin6Hours.RDS"))

# Subset For Clustering
GroupOI <- "Capsule all"

#All Capsules ###################################################
SubsetRows <- which(SampleMetaData$`SampleType`== "Capsule")
#####################

MetDF <- MetDF[SubsetRows,]
SampleMetaData <- SampleMetaData[SubsetRows,]
FoodWithinHours <- FoodWithinHours[SubsetRows,]
SampleMetaData <- cbind(SampleMetaData,FoodWithinHours)

#################################################################

cormat <- cor(t(MetDF), method = "spearman")
s <- cormat
diag(s) <- 0
s[is.na(s)] <- 0
hc <- hclust(as.dist(1-s), method="ward.D2") # ward method provide better clusters.

ClusterNum <- cutree(hc, 4)

clusterdf <- data.frame(SampleID = rownames(MetDF), 
                        Cluster = ClusterNum,
                        stringsAsFactors = F)
colnames(clusterdf) <- c("Sample ID","Cluster")

SampleMetaData$ClusterNum <- c(clusterdf[,2])

########################################
### Get VIP order for metabolites with highest VIP scores when separated by cluster
########################################

# Load Data For Analysis
FiltLogNormMetDF_C <- readRDS(FiltLogNormMetDF_C_path)
MetDFvip <- FiltLogNormMetDF_C

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
FVariables <- colnames(SampleMetaData)[c(3,12:23)]

SampleMetaData[,14:23] <- ifelse(SampleMetaData[,14:23]=="Yes", "black", "white")

FVN <- which(colnames(SampleMetaData)==FVariables[1])
SampleType <- SampleMetaData[1:nrow(df2plot),FVN]

######################

opls.PCA1 <- opls(df2plot)
plot(opls.PCA1, typeVc = "x-score", parAsColFcVn = SampleType)


#MetDFvip.PLSDA <- opls(df2plot,SampleType)
VIPscores <- getVipVn(opls(df2plot,SampleType))
VIPscoresDF <- as.data.frame(VIPscores)

VIPorderedPlot <- as.data.frame(cbind(rownames(VIPscoresDF)[order(VIPscoresDF$VIPscores,decreasing = TRUE)], VIPscoresDF[order(VIPscoresDF$VIPscores, decreasing = TRUE),]))

#########
# VIP Ordered Heatmap 
#########
NumTopMets <- 80

for (i in 1:ncol(MetDF)){
  MetDF[,i] <- as.data.frame(rescale(MetDF[,i], to=c(-1,1)))
}

tVIPs <- as.data.frame(matrix(nrow=274,ncol=NumTopMets))
for (i in 1:NumTopMets){
mdfCol <- which(VIPorderedPlot[i,1] == colnames(MetDF))
tVIPs[,i] <- MetDF[,mdfCol]
colnames(tVIPs)[i] <- colnames(MetDF)[mdfCol]
}
rownames(tVIPs)<- rownames(MetDF)

# Data Frame with correct Column Order and Row Order (Both Correlation Hierarchically Clustered)

cormatMs <- cor(tVIPs, method = "spearman")
sMs <- cormatMs
diag(sMs) <- 0
sMs[is.na(sMs)] <- 0
hcMs <- hclust(as.dist(1-sMs), method="ward.D2") # ward method provide better clusters.

hcdMs <- as.dendrogram(hcMs)
MetaboliteOrder <- labels(hcdMs)
MetabOrder <- as.data.frame(MetaboliteOrder,)

#Rearrange Columns of VIP df
MetOrderedVIPs <- select(tVIPs,MetabOrder[,1])
#Rearrange Rows of df
SampleOrder <- SampleMetaData[,c(1,FVN)]
SampleOrder <- arrange(SampleOrder,SampleMetaData[,2])

ColMetOrderedVIPs <- MetOrderedVIPs
rownames(ColMetOrderedVIPs) <- 1:nrow(ColMetOrderedVIPs)
for (i in 1:nrow(MetOrderedVIPs)){
  ColMetOrderedVIPs[i,] <- MetOrderedVIPs[which(rownames(MetOrderedVIPs) == SampleOrder[i,1]),]
rownames(ColMetOrderedVIPs)[i] <- SampleOrder[i,1]
}

#Add column order 
SamplePOrder <- as.factor(1:274)
HM2plot <- cbind(SamplePOrder,ColMetOrderedVIPs)

#Melt data
HeatMap2plot <- melt(HM2plot)
colnames(HeatMap2plot) <- c("SampleOrder","Metabolite","Metabolite Abundance")

# Build HeatMap

ggplot(HeatMap2plot, aes(x=`SampleOrder`, y=`Metabolite`, fill= `Metabolite Abundance`))+
    geom_tile()+
    coord_fixed()+
    scale_fill_gradient2(low= "#112047",
                         mid= "#ccd5dd",
                         high = "#f44323")+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=5),
          axis.text.x = element_text(angle = 90, size = 5),
          legend.text = element_text())

# Build Variable Row
Legg <- as.data.frame(matrix(ncol=2,nrow=nrow(HM2plot)))
colnames(Legg) <- c("Order","Fixed Variable")
Legg[,1] <- 1:nrow(HM2plot)
Legg[,2] <- SampleOrder[,2]

#Create a custom color scale and reorder to match legend order...
myColors <- SubjectColorScheme[,1]
myColors <- myColors[c(1,10,11,12,13,14,15,2,3,4,5,6,7,8,9)]

ggplot(Legg, aes(x=`Order`, y=1, fill= `Fixed Variable`))+
  geom_tile()+
  coord_fixed()+
  scale_fill_manual(name="Subject",values=myColors)+
  labs(title="Subjects of Study")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(angle = 90, size = 5),
        legend.text = element_text())

### Dendogram to test cluster distinguish
hcd <- as.dendrogram(hc)

DendOrder <- as.data.frame(labels(hcd),)
DendOrder[,2] <-  c(1:length(labels(hcd)))
colnames(DendOrder) <- c("Sample ID","LabelOrder")
MergedSampleMetaData <- merge(SampleMetaData[1:nrow(DendOrder),],DendOrder, sort=TRUE)

MergedSampleMetaData$Cluster <- clusterdf$Cluster

DendLabels <- as.data.frame(MergedSampleMetaData[order(MergedSampleMetaData$LabelOrder, decreasing=FALSE),c(3,12,13,7,14:23,25,1)])
colnames(DendLabels) <- c("Subject In Order","Capsule Type In Order", "Lunch or Dinner","Swallow Time",paste(colnames(SampleMetaData)[14:23]),"Cluster","SampleID")

#Set Subject Colors
DendLabels2plot <- DendLabels
for (i in 1:nrow(DendLabels)){
  DendLabels2plot[i,1] <- SubjectColorScheme[which(rownames(SubjectColorScheme) == DendLabels[i,1]),]
  DendLabels2plot[i,2] <- SampleTypeColorScheme[which(rownames(SampleTypeColorScheme) == DendLabels[i,2]),]

  }
DendLabels2plot[,3] <- ifelse(DendLabels$`Lunch or Dinner`=="Lunch", "gold", "black")
colnames(DendLabels2plot)[1:3] <- c("Subject Color","Capsule Color","Lunch or Dinner")

DendLabels2plot[5:14] <- ifelse(DendLabels[,5:14]=="Yes", "black", "white")
DendLabels2plot[15] <- ifelse(DendLabels[,15]==3, "green", "white")
