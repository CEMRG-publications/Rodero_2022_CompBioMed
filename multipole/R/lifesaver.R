#!/usr/bin/env Rscript

source('/home/crg17/Desktop/scripts/multipole/R/multipolar_pipeline.R')
Run_pipeline("HF","AT_table")

#hearts <- c(paste0("0",1:9),10:24)
#hearts <- c("01","10","20")

#registerDoParallel(cores=20)
#foreach(i=1:24) %dopar% {
#for(i in c(1:length(hearts))){
#Remove_base(heart=hearts[i],which_cases="HF",SA_folder="default")
#}

#system("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"Bases cropped finished\"")

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
