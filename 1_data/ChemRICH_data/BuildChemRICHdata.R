#######################################
# Build ChemRICH data
#######################################
rm(list=ls())
#setwd("~/Desktop/GItract15S-main")
setwd("C:/Users/jakef/Box/GItract15S-main")
# configure directories, load libraries and base functions
source("0_Configuration.R")
LMMresult_Filt6h <- read_excel(paste0(Tables_dir,"LinMixedModAllVsBeginEndNoNuts4_FilteredData6hour_NoFDR.xlsx"))
LMMresult_Filt6h <- read_excel(paste0(Tables_dir,"LinMixedMod_SEXage_results_nonCorrected.xlsx"))

# To Predict MESH build xlsx file with "compound_name","pubchem_id","smiles" and save to ChemRICH data folder
# https://sites.google.com/view/chemrich/predict-mesh-classes?authuser=0

setwd("C:/Users/jakef/Box/GItract15S-main/1_data/ChemRICH_data")

# THESE THREE LINES WILL RECOMPUTE MESH CLASSES, THIS WAS DONE THEN MANUALLY CURATED
#source("https://raw.githubusercontent.com/barupal/ChemRICH/master/predict_mesh_chemical_class.R")
#load.ChemRICH.Packages()
#predict_mesh_classes(inputfile = "chemrich_input_mesh_prediction_all.xlsx")

#Check smiles with probs ####
checkSmiles <- function(input = "inputfile.xlsx") {
  ndf <- readxl::read_xlsx(input)
  fps <- lapply(1:nrow(ndf), function(x) {
    rcdk::parse.smiles(ndf$smiles[x])
  })
  charvec <- sapply(fps, nchar)
  paste0("Lines with an incorrect SMILES codes are : ", paste(as.integer(which(charvec==4)), collapse = ","))
}
checkSmiles("chemrich_input_mesh_prediction_all.xlsx")

# Load ChemRICH plotting function
run_chemrich_chemical_classes <- function(inputfile = "nameofthefile") {
  
  ndf <- data.frame(inputfile)
  ndf$pvalue <- as.numeric(ndf$pvalue)
  ndf$effect_size <- as.numeric(ndf$effect_size)
  ndf$edirection <- "up"
  ndf$efs <- 1
  
  if(length(which(ndf$effect_size < 0)) >0) { # if regression models
    ndf$edirection[which(ndf$effect_size < 0)] <- "down"
    ndf$edirection[which(ndf$pvalue > 0.1)] <- "no change"
    ndf$efs[which(ndf$effect_size < 0)] <- 1/abs(ndf$effect_size[which(ndf$effect_size < 0)])
    ndf$efs[which(ndf$effect_size > 1)] <- abs(ndf$effect_size[which(ndf$effect_size > 1)])
    ndf$efs [which(ndf$pvalue > 0.1)] <- 1
    
  } else { # if student test
    ndf$edirection[which(ndf$effect_size < 1)] <- "down"
    ndf$edirection[which(ndf$pvalue > 0.1)] <- "no change"
    ndf$efs[which(ndf$effect_size < 1)] <- 1/ndf$effect_size[which(ndf$effect_size < 1)]
    ndf$efs[which(ndf$effect_size > 1)] <- ndf$effect_size[which(ndf$effect_size > 1)]
    ndf$efs [which(ndf$pvalue > 0.1)] <- 1
  }
  
  ndf$xlogp <- as.numeric(sapply(ndf$smiles, function(x)  { rcdk::get.xlogp(rcdk::parse.smiles(x)[[1]]) }))
  
  clusterids <- names(which(table(ndf$set)>2))
  clusterids <- clusterids[which(clusterids!="")]
  cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
    cl.member <- which(ndf$set==x)
    if( length(which(ndf$pvalue[cl.member]<.10)) >0 ){
      pval.cl.member <- ndf$pvalue[cl.member]
      p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
      p.test.results$p.value
    } else {
      1
    }
  })
  
  cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
  #clusterdf <- data.frame(name=clusterids[which(cluster.pvalues!=10)],pvalues=cluster.pvalues[which(cluster.pvalues!=10)], stringsAsFactors = F)
  clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)
  
  clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
    dfx <- ndf[which(ndf$set==x),]
    dfx$compound_name[which.min(dfx$pvalue)]
  })
  
  altrat <- sapply(clusterdf$name, function (k) {
    length(which(ndf$set==k & ndf$pvalue<0.1))/length(which(ndf$set==k))
  })
  
  uprat <-sapply(clusterdf$name, function (k) {
    length(which(ndf$set==k & ndf$pvalue<0.1 & ndf$edirection == "up"))/length(which(ndf$set==k & ndf$pvalue<0.10))
  })
  
  clust_s_vec <- sapply(clusterdf$name, function (k) {
    length(which(ndf$set==k))
  })
  
  clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(ndf$set==k & ndf$pvalue<0.10))})
  
  clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(ndf$set==k & ndf$pvalue<0.05 & ndf$edirection == "up"))})
  
  clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(ndf$set==k & ndf$pvalue<0.05 & ndf$edirection == "down"))})
  
  clusterdf$upratio <- uprat
  clusterdf$altratio <- altrat
  clusterdf$csize <- clust_s_vec
  clusterdf <- clusterdf[which(clusterdf$csize>2),]
  clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")
  
  clusterdf$xlogp <- as.numeric(sapply(clusterdf$name, function(x) {  median(ndf$xlogp[which(ndf$set==x)]) })) ##
  
  clusterdf$Compounds <- sapply(clusterdf$name, function(x) {
    dfx <- ndf[which(ndf$set==x),]
    paste(dfx$compound_name,collapse="<br>")
  })
  
  clustdf <- clusterdf[which(clusterdf$pvalues!=1),]
  
  #################################################
  ########## Impact Visualization Graph ###########
  #################################################
  
  clustdf.alt.impact <- clustdf[which(clustdf$pvalues<0.05 & clustdf$csize>1 & clustdf$alteredMetabolites>1) ,]
  clustdf.alt.impact <- clustdf.alt.impact[order(clustdf.alt.impact$xlogp),]
  clustdf.alt.impact$order <- order(clustdf.alt.impact$xlogp)
  clustdf.alt.impact$logPval <- -log(clustdf.alt.impact$pvalues)
  
  p2 <- ggplot(clustdf.alt.impact,aes(x=xlogp,y=-log(pvalues)))
  
  p2 <- p2 + geom_point(aes(size=csize, color=upratio)) +
    #labs(subtitle = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
    scale_color_gradient(low = "blue", high = "red", limits=c(0,1))+
    scale_size(range = c(5, 30)) +
    scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+4  )) +
    scale_x_continuous(" Lipophilicity (xlogp) ") +
    theme_bw() +
    labs(title = paste0("ChemRICH cluster plot for ",FixedVs[j])) +
    geom_label_repel(aes(label = name), color = "gray20",data=subset(clustdf.alt.impact, csize>2),force = 8,max.overlaps =20)+
    theme(
      plot.title = element_text(face="bold", size=30,hjust = 0.5),
      axis.title.x = element_text(face="bold", size=20),
      axis.title.y = element_text(face="bold", size=20, angle=90),
      panel.grid.major = element_blank(), # switch off major gridlines
      panel.grid.minor = element_blank(), # switch off minor gridlines
      legend.position = "none", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
      legend.title = element_blank(), # switch off the legend title
      legend.text = element_text(size=12),
      legend.key.size = unit(1.5, "lines"),
      legend.key = element_blank(), # switch off the rectangle around symbols in the legend
      legend.spacing = unit(.05, "cm"),
      axis.text.x = element_text(size=10,angle = 0, hjust = 1),
      axis.text.y = element_text(size=15,angle = 0, hjust = 1)
    )
   
  print(p2) #jf
  
#  read_pptx() %>%
#    add_slide(layout = "Title and Content", master = "Office Theme") %>%
#    ph_with(dml(ggobj = p2), location = ph_location(type = "body",width=10, height=8,left = 0, top = 0)) %>%
#    print(target = paste0("chemrich_class_impact_plot.pptx")) %>%
#    invisible()
  
#  ggsave(paste0("chemrich_class_impact_plot",FixedVs[2],".png"), p2,height = 8, width = 12, dpi=300)
  
#  cat(paste0("chemrich_class_impact_plot.pptx"," has been created.\n"))
  
  ## Export the result table.
  clustdf.e <- clusterdf[order(clusterdf$pvalues),]
  clustdf.e$pvalues <- signif(clustdf.e$pvalues, digits = 2)
  clustdf.e$adjustedpvalue <- signif(clustdf.e$adjustedpvalue, digits = 2)
  clustdf.e$upratio <- signif(clustdf.e$upratio, digits = 1)
  clustdf.e$altratio <- signif(clustdf.e$altratio, digits = 1)
  clustdf.e <- clustdf.e[,c("name","csize","pvalues","adjustedpvalue","keycpdname","alteredMetabolites","upcount","downcount","upratio","altratio")]
  names(clustdf.e) <- c("Cluster name","Cluster size","p-values","FDR","Key compound","Altered metabolites","Increased","Decreased","Increased ratio","Altered Ratio")
  #df1$TreeLabels <- treeLabels
  ndf$pvalue <- signif(ndf$pvalue, digits = 2)
  ndf$efs <- signif(ndf$efs, digits = 2)
  ndf$FDR <- signif(  p.adjust(ndf$pvalue), digits = 2)
  l <- list("ChemRICH_Results" = clustdf.e, "Compound_ChemRICH" = ndf )
  openxlsx::write.xlsx(l, file = paste0("chemRICH_class_results_SEXage_",FixedVs[j],".xlsx"), asTable = TRUE)
  cat(paste0("chemRICH_class_results.xlsx", " has been saved.\n"))

}

#Load MESH names and chemical groups
Classes <- read_excel("MeSh_Prediction_Results_curated.xlsx")

#bdf <- LMMresult_Filt12h
#bdf <- LMMresult_NonFilt6h
bdf <- LMMresult_Filt6h

FixedVs <- c("BeginOrEnd","veggie","meat_or_egg_or_fish","rice_or_pasta_or_bread_or_grain","coffee_or_tea","dessert","beer_or_wine","dairy","fruit","Antibiotics","Sex","Age")

# Build Table with compound name order pvalue effect size set
ChemRIChInputFile <- data.frame(compound_name = Classes$compound_name, 
                                smiles = Classes$smiles, 
                                pvalue = NA, 
                                effect_size = NA,
                                set = Classes$MeSH_Class)

# ChemRICH loop

#pdf(file="ChemRICHclusts_May2Updated_WithoutDipeptides.pdf", width=12, height=7)

for (j in 1:length(FixedVs)){
#  for (j in 1){
    
Vc <- grep(FixedVs[j],colnames(bdf))
bdf_1V <- as.data.frame(bdf[,c(2,Vc)])
bdf_1V$RichEffect <- lapply(bdf_1V[,2],function(x){if(x>0){2}else{0.5}})

for (i in 1:nrow(bdf)){
Cinn <- which(bdf_1V[i,1] == ChemRIChInputFile$compound_name)
ChemRIChInputFile[Cinn,3] <- bdf_1V[i,3] #pvalue
ChemRIChInputFile[Cinn,4] <- as.numeric(bdf_1V[i,4]) #ChemRICH effect size
}

### Run ChemRICH
run_chemrich_chemical_classes(ChemRIChInputFile)

}

#dev.off()

