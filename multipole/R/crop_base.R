#!/usr/bin/env Rscript

source('/home/crg17/Desktop/scripts/multipole/R/preprocessing.R')
Load_Install_Packages(c("doParallel"))

#hearts <- c(paste0("0",1:9),10:20)
#hearts <- c("01")

#registerDoParallel(cores=20)
#foreach(i=1:20) %dopar% {
#for(i in c(1:length(hearts))){
#Crop_base_from_tags(hearts[i],"RR","BiV_FEC_w0_h0",F)
#Crop_base_from_tags(hearts[i],"RR","BiV_FEC_w5_h33_retagged",F)
#Crop_base_from_tags(hearts[i],"RR","BiV_FEC_w5_h70_retagged",F)
#Crop_base_from_tags(hearts[i],"RR","BiV_FEC_w5_h100_retagged",F)
#if(hearts[i] != 21)
#Crop_base_from_tags(hearts[i],"HF","BiV_FEC_w5_h70_scar")
#}

#system("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"Base cropped in healthy finished\"")

hearts <- c(paste0("0",1:9),10:20,22:24)
#registerDoParallel(cores=20)
#foreach(i=1:length(hearts)) %dopar% {
for(i in c(5,14)){
#Crop_base_from_tags(hearts[i],"HF","BiV_FEC_w0_h0",F)
Crop_base_from_tags(hearts[i],"HF","BiV_FEC_w5_h33_retagged",F)                                               
#Crop_base_from_tags(hearts[i],"HF","BiV_FEC_w5_h70_retagged",F)                                               
#Crop_base_from_tags(hearts[i],"HF","BiV_FEC_w5_h100_retagged",F)
#if(hearts[i] != 21)
#Crop_base_from_tags(hearts[i],"HF","BiV_FEC_w5_h70_scar",F)
}

system("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"Base cropped in HF finished\"")

#CreateMonopoles(which_cases = "HF",SA_folder="default_nobase")

#system("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"Monopoles finished, running dipoles...\"")

#registerDoParallel(cores=20)
#foreach(i=1:24) %dopar% {
#CreateDipoles(SA_folder='default_nobase',which_cases='HF',heart=hearts[i])
#}

#system("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"Dipoles finished\"")


#registerDoParallel(cores=20)
#foreach(i=1:24) %dopar% {
#system(paste0("/home/crg17/Desktop/scripts/multipole/bin/build_HAC_table.o /data/SA_multipole/default_nobase/HF/",hearts[i]," /media/crg17/\"Seagate Backup Plus Drive\"/CT_cases/HF_case",hearts[i],"/meshing/1000um/BiV"))
#}

#system("/home/crg17/Desktop/scripts/multipole/bin/build_HAC_table_RV.o default_nobase HF")


#system("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"Building tables finished\"") 
