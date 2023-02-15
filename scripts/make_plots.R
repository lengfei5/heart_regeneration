##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: make plots
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Jan 23 16:03:26 2023
##########################################################################
##########################################################################
make_plot_umap_scRNA = FALSE
if(make_plot_umap_scRNA){
  annots = readRDS('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds')
  annots = annots[, c(1:7, 16)]
  
  write.table(annots, file = paste0("/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/", 
                                    "geneAnnotation_curated.geneSymbol.toUse.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  ##########################################
  # CM umap  
  ##########################################
  
  aa1 <- subset(aa,  subtypes %in%  c( 
    
    "B_cells_(FOXO1)", "B_cells_(SIGLEC11)",             "B_cells_Prol"     ,               
    "CM_Atria" ,
    
    "CM_Atria_Tagln", 
    "CM_IS"          ,             
    "CM_OFT"         ,        "CM_PM_(HCN4)" ,
    
    "CM_Prol_1"       ,             "CM_Prol_2"      ,              "CM_Prol_3"           ,       
    "CM_Prol_IS" ,
    
    "CM_ven_(Cav3_1)"        ,       "CM_ven_(Robo2)"       ,                                      "EC" ,
    
    "EC_(CEMIP)"      ,              "EC_(LHX6)"    ,                "EC_(NOS3)" ,                   "EC_(WNT4)" ,
    
    "EC_IS_(IARS1)"    ,              "EC_IS_(LOX)"   ,                "EC_IS_Prol"  ,                    "EC_Prol" ,
    
    "FB_(PKD1)"  ,                  "FB_(TNXB)"  ,                  "FB_(VWA2)" ,               "FB_IS_(TFPI2)" ,
    
    "FB_IS_(TNC)"  ,                    "FB_Prol",               "Megakeryocytes"  ,           "Mo/Macs_(FAXDC2)" ,
    
    "Mo/Macs_(SNX22)" ,                "Mo/Macs_Prol" ,            "Mo/Macs_resident"  ,                 "Neu_(DYSF)" ,
    
    "Neu_(IL1R1)" ,                    "Neuronal", "Proliferating_Megakeryocytes",            "Proliferating_RBC",
    
    "RBC"  ,                    "T_cells" ))
  
  
}


##########################################
# test linear discriminant analysis (LDA) 
##########################################
Test_LDA = FALSE
if(Test_LDA){
  library(MASS)
  library(ggplot2)
  
  #attach iris dataset to make it easy to work with
  attach(iris)
  
  #view structure of dataset
  str(iris)
  
  #scale each predictor variable (i.e. first 4 columns)
  iris[1:4] <- scale(iris[1:4])
  
  #find mean of each predictor variable
  apply(iris[1:4], 2, mean)
  
  #find standard deviation of each predictor variable
  apply(iris[1:4], 2, sd) 
  
  #make this example reproducible
  set.seed(1)
  
  #Use 70% of dataset as training set and remaining 30% as testing set
  sample <- sample(c(TRUE, FALSE), nrow(iris), replace=TRUE, prob=c(0.7,0.3))
  train <- iris[sample, ]
  test <- iris[!sample, ] 
  
  #fit LDA model
  model <- lda(Species~., data=train)
  
  #view model output
  model
  #define data to plot
  
  #use LDA model to make predictions on test data
  predicted <- predict(model, train)
  
  head(predicted$class)
  #view posterior probabilities for first six observations in test set
  head(predicted$posterior)
  
  #view linear discriminants for first six observations in test set
  head(predicted$x)
  
  lda_plot <- cbind(train, predicted$x)
  
  #create plot
  ggplot(lda_plot, aes(LD1, LD2)) +
    geom_point(aes(color = Species))
  
}
