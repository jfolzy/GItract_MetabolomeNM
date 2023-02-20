#######################################
# Linear Mixed Modeling
#######################################
rm(list=ls())
source("0_Configuration.R")

#LoadData 
Filt_MetDF <- readRDS(FiltLogNormMetDF_C_path) # has NA values

SampleMetaData <- readRDS(SampleMetaData_path)
FoodWithinHours <- readRDS(paste0(CleanDataPath,"FoodWithin6Hours.RDS"))

FoodAddn <- as.data.frame(matrix(nrow = nrow(SampleMetaData)-nrow(FoodWithinHours), ncol = ncol(FoodWithinHours),data = "No"))
colnames(FoodAddn) <- colnames(FoodWithinHours)
FoodWithinHoursAddn <- rbind(FoodWithinHours, FoodAddn)

MergedSampleMetaData <- cbind(SampleMetaData,FoodWithinHoursAddn)

BeginOrEnd <- ifelse(MergedSampleMetaData[,6] == "1" | MergedSampleMetaData[,6] == "2"  ,"Capsule 1 or 2", "Capsule 3 or 4")
MergedSampleMetaData$BeginOrEnd <- BeginOrEnd

MetDF <- Filt_MetDF
MetDF[is.na(MetDF)] <- 0

#remove metabolites with all zero values for capsules
ToRmv <- NA
for (i in 1:ncol(MetDF)){
  if(sum(MetDF[,i])<= 0){ToRmv <- c(ToRmv, i)}
}
ToRmv <- ToRmv[-1]
rm(ToRmv)

#Scale from -1 to 1
for (i in 1:ncol(MetDF)){
MetDF[,i] <- as.data.frame(rescale(MetDF[,i], to=c(-1,1)))
}

df <- cbind(MergedSampleMetaData[1:274,],MetDF)
df[,6] <- as.numeric(df[,6])

colnames(df)[(ncol(df)-ncol(MetDF)+1):ncol(df)] <- paste0("metabolite_", 1:ncol(MetDF))
colnames(df) <- gsub(" ", "_", colnames(df))


vars <- colnames(df[,(ncol(df)-ncol(MetDF)+1):ncol(df)])

#All Variables
FixedVs <- c("BeginOrEnd","veggie","meat_or_egg_or_fish","rice_or_pasta_or_bread_or_grain","coffee_or_tea","dessert","beer_or_wine","dairy","fruit","Antibiotics")

# If more than 1 random variable then need to modify the results table for loop
RandomVs <- c("Subject_Number")

results <- as.data.frame(matrix(nrow = length(vars), ncol=length(FixedVs)*2+length(RandomVs)*2+1))
colnames(results) <- c("metabolite", paste0("rnd_",RandomVs,"_variance"), paste0("rnd_",RandomVs,"_stdev"), paste0("fixed_",FixedVs,"_est"), paste0("fixed_",FixedVs,"_pval"))

#Calculate Results Table Columns for the loop
CL1 <- c(2:(1+length(RandomVs)))
CL2 <- c((max(CL1)+1):(max(CL1)+length(RandomVs)))
CL3 <- c((max(CL2)+1):(max(CL2)+length(FixedVs)))
CL4 <- c((max(CL3)+1):(max(CL3)+length(FixedVs)))

var="metabolite_1"
form = as.formula(paste0(paste0(var, " ~ "), paste(FixedVs, collapse=' + ' ), paste0("+ (1|",RandomVs, ")")))
print(form)

for (i in 1:length(vars)){
  var <- vars[i]
  form = as.formula(paste0(paste0(var, " ~ "), paste(FixedVs, collapse=' + ' ), paste0("+ (1|",RandomVs, ")")))
  model <- summary(lmer(form, data = df))
  
  results[i, 1] <- colnames(MetDF)[i]
  results[i, CL1] <- model$varcor$Subject_Number[1]
  results[i, CL2] <- "NA" #model$varcor Can't get random variable variance stdev
  results[i, CL3] <- as.numeric(model$coefficients[2:nrow(model$coefficients), 1]) # row 4 of coefficients is the food_of_interest row, column 1 is the regression estimate and column 5 is the p value
  results[i, CL4] <- model$coefficients[2:nrow(model$coefficients), 5]
  }

# Volcano Plots

for (i in 1:10) {
    
FV2P <- FixedVs[i]
df2plot <- results[,c(1,grep(FV2P, colnames(results)))]
colnames(df2plot) <- c("Metabolites","Fold Change","p-value")  

#convert fold changes from log slope to ... 
df2plot$Transformed <- exp(df2plot[,2]) ########################
df2plot$Padjusted <- p.adjust(df2plot[,3], method = "BH")


FCcutoff <- 0.20 
pvalueCutoff <- 0.05
df2plot$Importance = ifelse(df2plot$`p-value` < pvalueCutoff & abs(df2plot$`Fold Change`) >= FCcutoff, 
                      ifelse(df2plot$`Fold Change` > FCcutoff ,'Up','Down'),
                      'Stable')
df2plot$FDRSignificant = ifelse(df2plot$`Padjusted` < pvalueCutoff, 
                            "Significant with FDR correction",'Not significant with FDR correction')

df2plot$ToLabel <- NA
df2plot$ToLabel[which(df2plot$Importance != "Stable")] <- df2plot[which(df2plot$Importance != "Stable"),1]

p <- ggplot(data= df2plot, aes(x= `Fold Change`, y= -log10(`p-value`), color=`Importance`, label= `ToLabel`)) +
  geom_point(alpha=0.6, size=4.5, aes(shape=`FDRSignificant`)) +
  geom_text_repel() +
  scale_shape_manual(values=c(19,18))+
  scale_color_manual(values=c("blue", "grey","red"))+
  #xlim(c(-4.5, 4.5)) +
  geom_vline(xintercept=c(-FCcutoff,FCcutoff),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(`pvalueCutoff`) ,lty=4,col="black",lwd=0.8) +
  labs(x="Effect Size Coefficient",
       y="-log10 (adjusted p-value)",
       title=paste0("Metabolites significantly effected by ",gsub("_"," ",FV2P)))  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank())
print(p)

}

