#######################################
# Correlation To Caffeine
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source("0_Configuration.R")

#Data
MetDF <- readRDS(FiltLogNormMetDF_CS_path)
SampleMetaData <- readRDS(SampleMetaData_path)

SubsetRows <- 1:274 #All samples

MetDF <- MetDF[SubsetRows,]
SampleMetaData <- SampleMetaData[SubsetRows,]

# Correlation Matrix
MetDF_CorObj <- rcorr(as.matrix(MetDF), type="spearman")
CorMatrix <- as.data.frame(MetDF_CorObj$r)
PvalMatrix <- as.data.frame(MetDF_CorObj$P)

#Build table with top correlating metabolites to each metabolite
TopCoeffTable <- as.data.frame(sapply(1:ncol(CorMatrix),function(i){
  head(rownames(CorMatrix)[order(rank(-CorMatrix[,i]))],100)
}))

colnames(TopCoeffTable) <- colnames(CorMatrix)
rownames(TopCoeffTable) <- c(1:100)

#Select Metabolite for table
MetToCorr <- "Caffeine"

OneCorrTable <- as.data.frame(cbind(TopCoeffTable[MetToCorr],round(as.numeric(CorMatrix[head(order(rank(-CorMatrix[,MetToCorr])),100),MetToCorr]),2)))
pvals <- NULL
for (i in 1:100){TP <- PvalMatrix[which(OneCorrTable[i,1] == rownames(PvalMatrix)),which(MetToCorr == colnames(PvalMatrix))]
pvals <- c(pvals, TP)}
pvals <- p.adjust(pvals, method = "BH", n=ncol(MetDF))
OneCorrTable <- cbind(OneCorrTable,pvals)
colnames(OneCorrTable) <- c("Metabolite",paste("Correlation to",MetToCorr),"p-value(Spearman, BH FDR corrected)")

MetCol <- which(colnames(CorMatrix)== MetToCorr)

Met1 <- grep("Caffeine",colnames(MetDF))
Met2 <- grep("Theobromine",colnames(MetDF))

df2plot <- as.data.frame(cbind(MetDF[,Met1],MetDF[,Met2]))
colnames(df2plot) <- c(colnames(MetDF)[Met1], colnames(MetDF)[Met2])

Met12_lm <- lm(df2plot[,2] ~ df2plot[,1])

ggplot(df2plot, aes(x=df2plot[,1], y=df2plot[,2]))+
  #geom_point()+
  geom_beeswarm(groupOnX = TRUE, alpha=0.7, size=2, color="brown")+
  geom_abline(slope= coef(Met12_lm)[["df2plot[, 1]"]],
              intercept= coef(Met12_lm)[["(Intercept)"]],
              size=1.5, color="black")+
  labs(x=paste0("log10(",colnames(MetDF)[Met1]," Peak Hieght)"),
       y=paste0("log10(",colnames(MetDF)[Met2]," Peak Hieght)"))

