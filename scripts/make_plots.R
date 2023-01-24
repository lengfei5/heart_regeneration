##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: make plots
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Jan 23 16:03:26 2023
##########################################################################
##########################################################################

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

