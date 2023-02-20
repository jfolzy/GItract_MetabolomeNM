#######################################
# MS meta data stats
#######################################

rm(list=ls())
source("0_Configuration.R")

SampleMetaData <- readRDS(SampleMetaData_path)

#PH 

df2plot <- as.data.frame(as.factor(SampleMetaData$`Sample Type`)[which("Capsule" == SampleMetaData$SampleType)])
df2plot[,2] <- SampleMetaData$pH[which("Capsule" == SampleMetaData$SampleType)]
colnames(df2plot) <- c("Capsule Type","pH")
my_comparisons <- list(c('Capsule 1','Capsule 2'),c('Capsule 1','Capsule 3'),c('Capsule 1','Capsule 4'),c('Capsule 2','Capsule 3'),c('Capsule 2','Capsule 4'),c('Capsule 3','Capsule 4'))

ggplot(df2plot, aes(y=df2plot$pH, x=`Capsule Type`, fill=`Capsule Type`))+
  geom_boxplot(outlier.shape = NA, alpha = 0.9)+
  geom_jitter(width=0.06,fill="black",alpha=0.5,size=2)+
  scale_fill_manual(values=CapTypeAndStoolColors) +
  stat_compare_means(method = "wilcox.test", 
                     paired = FALSE, 
                     aes(label = ..p.signif..), 
                     comparisons=my_comparisons) +
  labs(y='Sample pH',x="Capsule Type",
       title= "pH of samples by capsule type") +
  theme_bw()+
  theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank())
