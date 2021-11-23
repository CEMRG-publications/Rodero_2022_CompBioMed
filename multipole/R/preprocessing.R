Correct_HFLA_leads <- function(foldername,casenumber, flag_debugging = FALSE,
                               root_directory = "/media/crg17/Seagate Expansion Drive/SA_multipole/"){
  
  source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/leads_operations.R")
  
  # print(paste0("Reading ","/data/SA_multipole/",foldername,"/HF/",casenumber,"/multipole_AT1090.dat"))
  multipole_AT1090 <- read.csv(paste0(root_directory,foldername,"/HF/",casenumber,"/multipole_AT1090.dat"), sep="", stringsAsFactors=FALSE)
  multipole_AT1090$lead <- ascii2bin_lead(bin2ascii_lead(multipole_AT1090$lead))
  multipole_AT1090_corrected <- multipole_AT1090
  
  # The AT activating the base should be shorter
  
  if(multipole_AT1090$LA[multipole_AT1090$lead == ascii2bin_lead(1)] < multipole_AT1090$LA[multipole_AT1090$lead == ascii2bin_lead(8)]){
    for(i in c(1:nrow(multipole_AT1090_corrected))){
      multipole_AT1090_corrected[i,"LA"] <- multipole_AT1090[which(multipole_AT1090$lead == intToUtf8(rev(utf8ToInt(multipole_AT1090$lead[i])))),"LA"]
    }
  }
  
  write.table(multipole_AT1090_corrected,file=paste0(root_directory,foldername,"/HF/",casenumber,"/multipole_AT1090_leadcorrected.dat"),quote = FALSE,sep = " ",dec = ".",row.names = FALSE,col.names = TRUE)
}

#' @description Script to create the monopole activation files from individual
#' electrodes simulations. It just move them and merged them with the BiV RV
#' apex.
#' 
#' @param which_cases RR/h or HF for reverse remodelling or heart failure
#' respectively. "Both" to run both cases.
#' @param SA_folder Name of the folder of the sensitivity analysis. It will be
#' the output folder name and the input folder name if ekbatch=FALSE or the 
#' preffix if ekbatch=TRUE.
#' @param ekbatch If TRUE uses the output from using ekbatch to run the
#' simulations. Otherwise it works as it did before discovering it.
#' @param flag_debugging If TRUE it prints all the files reading and writing.
#' 
#' @return It creates folders in the /data/SA_multipole folder with the
#' different monopolar activations.
CreateMonopoles <- function(which_cases,SA_folder,ekbatch=TRUE, midseptum = FALSE,
                            flag_debugging = FALSE, with_scar = FALSE,
                            root_directory = "/media/crg17/Seagate Expansion Drive/SA_multipole/"){
  if(which_cases == "both"){
    CreateMonopoles(which_cases = "RR", SA_folder = SA_folder,
                    ekbatch = ekbatch , midseptum = midseptum,
                    flag_debugging = flag_debugging, root_directory = root_directory,
                    with_scar = with_scar)
    CreateMonopoles(which_cases = "HF", SA_folder = SA_folder,
                    ekbatch = ekbatch , midseptum = midseptum,
                    flag_debugging = flag_debugging, with_scar = with_scar,
                    root_directory = root_directory)
  }
  else{
  source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/common_functions.R")
  Load_Install_Packages(c("icesTAF","jjb"))
  
  if(midseptum == FALSE){
    bivpacing_name <- "RV_endo.apex"
  }
  else{
    bivpacing_name <- "midseptum"
  }
  
  if(which_cases == "RR"){
    which_cases <- "h"
  }
  if(which_cases == "HF"){
    num_cases <- 24
  }
  else if(which_cases == "h"){
    num_cases <- 20
  }
  
  if(!ekbatch){
    vein_names <- c("AN","AL","LA","PL","PO")
  }
  else{
    vein_names <- c("AN","AL","LA","IL","IN")
  }
  
  hearts <- c(paste0("0",c(1:9)),10:num_cases)
  
  if(with_scar){
    hearts <- hearts[-c(13,21)]
  }
  
  for(heart in hearts){
    if(!ekbatch){
      file2read <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                          which_cases,"_case",heart,"/simulations/multipole/","
                          eikonal_",SA_folder,"/BiV.",bivpacing_name,"/",
                          "vm_act_seq.dat")
   }
    else{
      file2read <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                          which_cases,"_case",heart,"/simulations/multipole/",
                          "init_files/",SA_folder,
                          "_BiV.",bivpacing_name,"_elemwise.dat")
    }
    Debug_message(file2read,"r",flag_debugging)
    AT_apex <- read.table(file2read)
    
    for(vein in vein_names){
      for(electrode in c(1:8)){
        if(!ekbatch){
          file2read <- paste0("/media/crg17/Seagate Backup Plus Drive/",
                              "CT_cases/",which_cases,"_case",heart,
                              "/simulations/multipole/eikonal_",SA_folder,"/",
                              vein,"_",electrode,"/vm_act_seq.dat")
          }
        else{
          file2read <- paste0("/media/crg17/Seagate Backup Plus Drive/",
                              "CT_cases/",which_cases,"_case",heart,
                              "/simulations/multipole/init_files/",SA_folder,
                              "_",vein,"_",electrode,"_elemwise.dat")
        }
        directory <- paste0(root_directory,SA_folder,"/",which_cases,"/",
                            heart)
        outname <- paste0(which_cases,heart,"_",vein,"_",
                          ascii2bin_lead(electrode),".dat")
        if(!(file.exists(paste0(directory,"/",outname)))){
          Debug_message(file2read,"r",flag_debugging)
          AT <- read.table(file2read)
          
          AT_res <- pmin(AT_apex,AT)
          
          
          
          jjb::mkdir(directory, r = TRUE)
          
          Debug_message(paste0(directory,"/",outname),"w",flag_debugging)
          write.table(AT_res,file=paste0(directory,"/",outname),quote = FALSE,
                      row.names = FALSE,col.names = FALSE)
        }
      }
    }
  }
  }
}

#' @description Script to create the dipoles from _all_ the files in the
#' monopoles folder.
#' 
#' @param SA_folder Name of the folder of the sensitivity analysis. It will
#' be input and output folder.
#' @param which_cases RR/h or HF for reverse remodelling or heart failure
#' respectively. "Both" to run both cases.
#' @param heart Number of the heart to read the monopoles from (written with
#' two ciphers). Use "all" to run it in all the hearts.
#' @param ekbatch Boolean. If true, the veins are named with the new naming, 
#' otherwise we used the previous names.
#' @param flag_debugging If TRUE it prints all the files reading and writing.
#' 
#' @return All the files with the bipoles in the SA_folder.

CreateDipoles <- function(SA_folder,which_cases,heart,ekbatch=TRUE,
                          flag_debugging = FALSE, with_scar = FALSE,
                          root_directory = "/media/crg17/Seagate Expansion Drive/SA_multipole/"){
  source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/common_functions.R")
  if(which_cases=="both"){
    CreateDipoles(SA_folder,which_cases = "RR",heart,ekbatch,flag_debugging,
                  root_directory = root_directory, with_scar = with_scar)
    CreateDipoles(SA_folder,which_cases = "HF",heart,ekbatch,flag_debugging,
                  root_directory = root_directory, with_scar = with_scar)
  }
  else{
  if(heart=="all"){
    source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/common_functions.R")
    heart_list <- c(paste0("0",1:9),10:20)
    if(which_cases=="HF"){
      heart_list <- c(heart_list,21:24)
      if(with_scar){
        heart_list <- heart_list[-c(13,21)]
      }
    }
    
    # for(h in heart_list){
    registerDoParallel(cores=20)
    
    foreach(i=1:length(heart_list)) %dopar% {
      CreateDipoles(SA_folder,which_cases,heart = heart_list[i],ekbatch,
                    flag_debugging, root_directory = root_directory)
    }
  }
  else{
  
  if(which_cases == "RR"){
    which_cases <- "h"
  }
  
  if(flag_debugging){
    print(paste0("Working on case ",heart))
  }
  files_vec <-  list.files(paste0(root_directory,SA_folder,"/",which_cases,"/",heart))
  if(ekbatch){
    vein_names <- c("AN","AL","LA","IL","IN")
  }
  else{
    vein_names <- c("AN","AL","LA","PL","PO")
  }
  for (vein in vein_names) {
    print(paste0("Vein ",vein))
    files_vein <- files_vec[which(grepl(vein,files_vec))]
    indices <- combn(length(files_vein),2)
    for(i in c(1:ncol(indices))){

      file2read <- paste0(root_directory,SA_folder,"/",which_cases,"/",heart,"/",files_vein[indices[1,i]])
      Debug_message(file2read,"r",flag_debugging)
      AT1 <- read.table(file2read)
      
      file2read <- paste0(root_directory,SA_folder,"/",which_cases,"/",heart,"/",files_vein[indices[2,i]])
      Debug_message(file2read,"r",flag_debugging)
      AT2 <- read.table(file2read)
      
      AT_res <- pmin(AT1,AT2)
      
      directory <- paste0(root_directory,SA_folder,"/",which_cases,"/",heart)
      name1 <- substr(files_vein[indices[1,i]],nchar(files_vein[indices[1,i]])-11,nchar(files_vein[indices[1,i]])-4)
      name2 <- substr(files_vein[indices[2,i]],nchar(files_vein[indices[2,i]])-11,nchar(files_vein[indices[2,i]])-4)
      
      outname <- paste0(which_cases,heart,"_",vein,"_",merge_design(c(name1,name2)),".dat")
      
      Debug_message(paste0(directory,"/",outname),"w",flag_debugging)
      write.table(AT_res,file=paste0(directory,"/",outname),quote = FALSE,row.names = FALSE,col.names = FALSE)
    }
  }
  }
  }
}



#' @description Sets the values of the activation times corresponding to the base as 10e4.
#' @details The base is defined as whatever outside of the AHA map in the LV and whatever over 0.9 in the UVC_Z for the RV.
#' 
#' @param heart Number of the case with two digits.
#' @param which_cases Healthy ("h" or "RR") or heart failure ("HF").
#' @param SA_folder Name of the folder of the simulations (after the eikonal_ preffix).
#' @param all_subfolders Default TRUE. If false, only does the RV apex.
#' @param flag_debugging Boolean if true prints whatever is reading and writing.
#' @param num_cores Number of cores used for the parallellisation step.
#' 
#' @return Writes a folder with the same name as the input but with the suffix "_nobase"

Crop_base_from_AT <- function(heart,which_cases,SA_folder,all_subfolders = TRUE,
                              flag_debugging = FALSE,num_cores = 20){
  Load_Install_Packages("jjb")
  
  if(which_cases == "RR")
    which_cases <- "h"
  
  if(flag_debugging)
    print(paste0("Reading /media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/meshing/1000um/BiV/UVC/COORDS_Z.dat"))
  
  COORDS_Z <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/meshing/1000um/BiV/UVC/COORDS_Z.dat"), header = FALSE)
  COORDS_Z <- COORDS_Z[[1]]
  
  if(flag_debugging)
    print(paste0("Reading /media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/meshing/1000um/BiV/BiV_AHA.dat"))
  
  BiV_AHA <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/meshing/1000um/BiV/BiV_AHA.dat"), header = FALSE)
  BiV_AHA <- BiV_AHA[[1]]
  
  if(flag_debugging)
    print(paste0("Reading /media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/meshing/1000um/BiV/UVC/COORDS_V.dat"))
  
  COORDS_V <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/meshing/1000um/BiV/UVC/COORDS_V.dat"), header = FALSE)
  COORDS_V <- COORDS_V[[1]]  
  
  subfolder_names <- c("BiV.RV_endo.apex")
  
  if(all_subfolders){
    subfolder_names <- c(subfolder_names,
                         "AN_1","AN_2","AN_3","AN_4","AN_5","AN_6","AN_7","AN_8",
                         "AL_1","AL_2","AL_3","AL_4","AL_5","AL_6","AL_7","AL_8",
                         "LA_1","LA_2","LA_3","LA_4","LA_5","LA_6","LA_7","LA_8",
                         "PL_1","PL_2","PL_3","PL_4","PL_5","PL_6","PL_7","PL_8",
                         "PO_1","PO_2","PO_3","PO_4","PO_5","PO_6","PO_7","PO_8")
  }
  
  full_AT <- as.data.frame(matrix(nrow = length(COORDS_Z), ncol = length(subfolder_names)))
  
  registerDoParallel(cores=num_cores)
  foreach(i=1:length(subfolder_names)) %dopar% {
    # for(i in c(1:length(subfolder_names))){
    if(flag_debugging)
      print(paste0("Reading /media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/simulations/multipole/eikonal_",SA_folder,"/",subfolder_names[i],"/vm_act_seq.dat"))
    
    full_AT[,subfolder_names[i]] <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/simulations/multipole/eikonal_",SA_folder,"/",subfolder_names[i],"/vm_act_seq.dat"), header = FALSE)
    
    # Outside of the AHA map
    full_AT[(BiV_AHA == 0)&(COORDS_V == -1),subfolder_names[i]] <- 1e+100
    # Over 0.9 of the UVC_Z
    full_AT[(COORDS_Z > 0.9)&(COORDS_V == 1),subfolder_names[i]] <- 1e+100
    
    jjb::mkdir(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/simulations/multipole/eikonal_",SA_folder,"_nobase/",subfolder_names[i]), r = TRUE)
    if(flag_debugging)
      print(paste0("Writing /media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/simulations/multipole/eikonal_",SA_folder,"_nobase/",subfolder_names[i],"/vm_act_seq.dat"))
    
    write.table(full_AT[,subfolder_names[i]],paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/simulations/multipole/eikonal_",SA_folder,"_nobase/",subfolder_names[i],"/vm_act_seq.dat"),row.names = FALSE,col.names = FALSE,quote = FALSE)
  }
}

#' @description Sets the values of the tags corresponding to the base to -1.
#' @details The base is defined as whatever outside of the AHA map in the LV 
#' and whatever over ~0.9 in the UVC_Z for the RV.
#' 
#' @param heart Number of the case with two digits.
#' @param which_cases Healthy ("h" or "RR") or heart failure ("HF").
#' @param BiV_name Name of the elem file located in the FEC folder.
#' @param UVC_folder Name of the folder for the UVCs. In the case of the RR
#' cases, corresponds to the most outer folder.
#' @param output_suffix Suffix of the elem file that will be written. Name of 
#' the file will be the original BiV_name + this suffix.
#' @param flag_debugging Boolean if true prints whatever is reading and writing.
#' 
#' @return Writes an elem file with the input name and the suffix "_nobase".
Crop_base_from_tags <- function(heart, which_cases, BiV_name,
                                UVC_folder = "UVC", output_suffix = "nobase",
                                flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr"))
  
  options(stringsAsFactors = F)
  
  if(which_cases == "RR")
    which_cases <- "h"
  
  common_directory <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                             which_cases,"_case",heart,"/meshing/1000um/BiV")

  BiV_elem <- paste0(common_directory, "/FEC/", BiV_name, ".elem") %>%
              Read_table(., skip = 1, header = F, sep = " ",
                       flag_debugging = flag_debugging)
  
  colnames(BiV_elem) <- c("Tt", "p1", "p2", "p3", "p4", "tag")
  
  if(which_cases == "HF"){

    COORDS_Z <- paste0(common_directory, "/", UVC_folder, "/COORDS_Z.pts") %>%
                Read_table(., header = F, skip = 1,
                           flag_debugging = flag_debugging)

    COORDS_V <- paste0(common_directory,"/", UVC_folder, "/COORDS_V.dat") %>%
                Read_table(., header = F, flag_debugging = flag_debugging)
  }
  else{
    COORDS_Z <- paste0(common_directory,"/", UVC_folder,
                       "/UVC/COORDS_Z.pts") %>%
                Read_table(., header = F, skip = 1,
                           flag_debugging = flag_debugging)
    
    COORDS_V <- paste0(common_directory,"/", UVC_folder,
                       "/UVC/COORDS_V.dat") %>%
                Read_table(., header = F, flag_debugging = flag_debugging)
  }

  COORDS_Z <- COORDS_Z[[1]]
  COORDS_V <- COORDS_V[[1]]
  
  BiV_AHA <- paste0(common_directory,"/BiV_AHA.dat") %>%
             Read_table(., header = F, flag_debugging = flag_debugging)
  
  BiV_AHA <- BiV_AHA[[1]]
  
  # Everything is pointwise except for the elem file
  
  #Outside of the AHA map. Indices TRUE/FALSE(pointwise)
  idx_to_change_LV <- which((BiV_AHA == 0)&(COORDS_V == -1)) 
  idx_to_change_RV <- which((COORDS_Z > 0.87)&(COORDS_V == 1))
  idx_to_change <- c(idx_to_change_LV,idx_to_change_RV) %>%
                  unique(.)
  # The indices start with 1 but in the points start with 0
  idx_to_change <- idx_to_change - 1
  
  elem_to_change <- which(BiV_elem$p1 %in% idx_to_change)
  elem_to_change <- c(elem_to_change,which(BiV_elem$p2 %in% idx_to_change))
  elem_to_change <- c(elem_to_change,which(BiV_elem$p3 %in% idx_to_change))
  elem_to_change <- c(elem_to_change,which(BiV_elem$p4 %in% idx_to_change)) %>%
                    unique(.)
  
  BiV_elem$tag[elem_to_change] <- -1
  
  
  BiV_towrite <- paste(BiV_elem$Tt,BiV_elem$p1,BiV_elem$p2,BiV_elem$p3,
                       BiV_elem$p4,BiV_elem$tag) %>%
                  as.data.frame(.) %>%
                  rbind(nrow(.),.)
  
  paste0(common_directory, "/meshes/", BiV_name, "_", output_suffix, ".elem") %>%
  Write_table(BiV_towrite, file = ., quote = F, col.names = F, row.names = F,
              sep = " ", flag_debugging = flag_debugging)
}

Create_first_electrode_single_vein <- function(heart,which_cases,vein,flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr","magrittr"))
  
  
  if(as.numeric(heart) < 10){
    heart <- paste0("0",as.numeric(heart))
  }
  if(which_cases == "RR" || which_cases == "h"){
    which_cases <- "h"
    UVC_dir <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                      which_cases,"_case",heart,"/meshing/1000um/BiV/UVC/UVC")
  }
  else if(which_cases == "HF"){
    UVC_dir <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                      which_cases,"_case",heart,"/meshing/1000um/BiV/UVC")
  }
  
  if(flag_debugging)
    Debug_message(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                  which_cases,"_case",heart,"/meshing/1000um/BiV/BiV_AHA.dat"),"r")
  
  BiV_AHA <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                        which_cases,"_case",heart,"/meshing/1000um/BiV/BiV_AHA.dat"),header = F)
  BiV_AHA <- BiV_AHA[[1]]
  
  if(flag_debugging)
    Debug_message(paste0(UVC_dir,"/COORDS_PHI.dat"),"r")
  
  COORDS_PHI <- read.table(paste0(UVC_dir,"/COORDS_PHI.dat"),header = F)
  COORDS_PHI <- COORDS_PHI[[1]]
  
  if(flag_debugging)
    Debug_message(paste0(UVC_dir,"/COORDS_Z.dat"),"r")
  
  COORDS_Z <- read.table(paste0(UVC_dir,"/COORDS_Z.dat"),header = F)
  COORDS_Z <- COORDS_Z[[1]]
  
  if(flag_debugging)
    Debug_message(paste0(UVC_dir,"/COORDS_RHO.dat"),"r")
  
  COORDS_RHO <- read.table(paste0(UVC_dir,"/COORDS_RHO.dat"),header = F)
  COORDS_RHO <- COORDS_RHO[[1]]
  
  # To avoid the discontinuity on -pi/pi
  original_PHI <- COORDS_PHI
  neg_PHI <- which(COORDS_PHI < 0)
  COORDS_PHI[neg_PHI] <- COORDS_PHI[neg_PHI] + 2*max(COORDS_PHI)
  
  df_total <- data.frame(BiV_AHA,COORDS_Z,COORDS_PHI,COORDS_RHO,original_PHI)
  #Indices start in 0
  rownames(df_total) <- as.numeric(rownames(df_total)) - 1

  
  # We want it epicardial
  
  new_df <- df_total[df_total$COORDS_RHO == 1,]
  PHI_eps <- 0.05
  
  
    if(vein == "AN"){
      new_df <- new_df[new_df$BiV_AHA == 6,]
      
      PHI_centre <- 0
      PHI_lower_eps <- 0
      PHI_upper_eps <- PHI_eps
    }
  if(vein == "AL"){
    new_df <- new_df[new_df$BiV_AHA == 6,]
  
    PHI_centre <- 0.5
    PHI_lower_eps <- PHI_eps
    PHI_upper_eps <- PHI_eps
  }
  if(vein == "LA"){
    new_df <- new_df[new_df$BiV_AHA == 5,] # Could be changed 
    
    PHI_centre <- 0
    PHI_lower_eps <- 0
    PHI_upper_eps <- PHI_eps
  }
  if(vein == "IL"){
    new_df <- new_df[new_df$BiV_AHA == 5,]
    
    PHI_centre <- 0.5
    PHI_lower_eps <- PHI_eps
    PHI_upper_eps <- PHI_eps
  }
  if(vein == "IN"){
    new_df <- new_df[new_df$BiV_AHA == 5,]
    
    PHI_centre <- 1
    PHI_lower_eps <- PHI_eps
    PHI_upper_eps <- 0
  }
  # We define the starting bands
  Z_centre <- 0.9
  Z_eps <- 0.05
  
  min_Z <- min(new_df$COORDS_Z)
  max_Z <- max(new_df$COORDS_Z)
  height <- max_Z - min_Z
   
  min_PHI <- min(new_df$COORDS_PHI)
  max_PHI <- max(new_df$COORDS_PHI)
  width <- max_PHI - min_PHI
    
    Z_greater_than <- min_Z + Z_centre*height
    Z_minor_than <- min_Z + (Z_centre+Z_eps)*height
    
    PHI_greater_than <- min_PHI + (PHI_centre - PHI_lower_eps)*width
    PHI_minor_than <- min_PHI + (PHI_centre + PHI_upper_eps)*width

    # We want it to be in the most upper part
    new_df_temp <- new_df[new_df$COORDS_Z >= Z_greater_than,]
    
    # If not, we reduce how high it is
    while(nrow(new_df_temp) == 0){
      Z_centre <- Z_centre - 0.01
      Z_minor_than <- min_Z + (Z_centre+Z_eps)*height
      Z_greater_than <- min_Z + Z_centre*height
      
      new_df_temp <- new_df[new_df$COORDS_Z >= Z_greater_than,]
      
    }
    
    new_df <- new_df_temp
    
    # We want it to be in the Z band
    new_df_temp <- new_df[new_df$COORDS_Z <= Z_minor_than,]
    
    # If not, we increase the tolerance
    
    while(nrow(new_df_temp) == 0){
      Z_centre <- Z_centre + 0.01
      Z_minor_than <- min_Z + (Z_centre+Z_eps)*height
  
      new_df_temp <- new_df[new_df$COORDS_Z <= Z_minor_than,]
      
    }
    
    new_df <- new_df_temp
    
    
    # We want it to be in the PHI band
    new_df_temp <- new_df[new_df$COORDS_PHI >= PHI_greater_than,]
    new_df_temp <- new_df_temp[new_df_temp$COORDS_PHI <= PHI_minor_than,]
    
    # If not, we increase how wide the band it is
    while(nrow(new_df_temp) == 0){
      PHI_eps <- PHI_eps + 0.05
      
      PHI_lower_eps <- (PHI_lower_eps > 0)*(PHI_lower_eps + PHI_eps)
      PHI_upper_eps <- (PHI_upper_eps > 0)*(PHI_upper_eps + PHI_eps)
      
      PHI_greater_than <- min_PHI + (PHI_centre - PHI_lower_eps)*width
      PHI_minor_than <- min_PHI + (PHI_centre + PHI_upper_eps)*width
      
      new_df_temp <- new_df[new_df$COORDS_PHI >= PHI_greater_than,]
      new_df_temp <- new_df_temp[new_df_temp$COORDS_PHI <= PHI_minor_than,]
      
      }
    
    new_df <- new_df_temp
    
    # If there are several optins, we take the one closer to the original
    # PHI coordinate:

  idx <- abs(new_df$COORDS_PHI - (min_PHI + PHI_centre*width)) %>%
         which.min(.) %>%
         extract2(rownames(new_df),.) %>%
         as.numeric(.)
    
  return(idx)

}

Create_n_electrodes <- function(heart,which_cases,vein,distance=7500,
                                num_elec=8,action="w1by1",flag_debugging = FALSE){
  
  if(as.numeric(heart) < 10){
    heart <- paste0("0",as.numeric(heart))
  }
  if(which_cases == "RR" || which_cases == "h"){
    which_cases <- "h"
    UVC_dir <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                     which_cases,"_case",heart,"/meshing/1000um/BiV/UVC/UVC")
  }
  else if(which_cases == "HF"){
    UVC_dir <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                      which_cases,"_case",heart,"/meshing/1000um/BiV/UVC")
  }
  
  if(flag_debugging)
    Debug_message(paste0(UVC_dir,"/COORDS_PHI.dat"),"r")
  
  COORDS_PHI <- read.table(paste0(UVC_dir,"/COORDS_PHI.dat"),header = F)
  COORDS_PHI <- COORDS_PHI[[1]]
  
  if(flag_debugging)
    Debug_message(paste0(UVC_dir,"/COORDS_Z.dat"),"r")
  
  COORDS_Z <- read.table(paste0(UVC_dir,"/COORDS_Z.dat"),header = F)
  COORDS_Z <- COORDS_Z[[1]]
  
  if(flag_debugging)
    Debug_message(paste0(UVC_dir,"/COORDS_RHO.dat"),"r")
  
  COORDS_RHO <- read.table(paste0(UVC_dir,"/COORDS_RHO.dat"),header = F)
  COORDS_RHO <- COORDS_RHO[[1]]
  
  # To avoid the discontinuity on -pi/pi
  original_PHI <- COORDS_PHI
  neg_PHI <- which(COORDS_PHI < 0)
  COORDS_PHI[neg_PHI] <- COORDS_PHI[neg_PHI] + 2*max(COORDS_PHI)
  
  if(flag_debugging)
    Debug_message(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                         which_cases,"_case",heart,"/meshing/1000um/BiV/BiV.pts"),
                  "r")
  
  BiV_pts <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                        which_cases,"_case",heart,"/meshing/1000um/BiV/BiV.pts")
                        ,skip = 1)
  colnames(BiV_pts) <- c("x","y","z")
  rownames(BiV_pts) <- as.numeric(rownames(BiV_pts))-1
  
  df_original <- data.frame(BiV_pts,COORDS_PHI,COORDS_Z,COORDS_RHO)
  # We restrict to be in the epicardium
  
  df_original <- df_original[df_original$COORDS_RHO == 1,]
  
  electrode_vtx <- Create_first_electrode_single_vein(heart=heart,which_cases=which_cases,vein=vein)
  if(action == "w1by1"){
    if(flag_debugging)
      Debug_message(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                    which_cases,"_case",heart,"/simulations/multipole/electrodes/"
                    ,vein,"_1.vtx"),"w")
    write.table(c("1","intra",electrode_vtx),
                paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                which_cases,"_case",heart,"/simulations/multipole/electrodes/",
                vein,"_1.vtx"),quote = F,sep = "\n",row.names = F,col.names = F)
  }
  
  
  electrodes_vector <- c(electrode_vtx)
  
  elec_created <- 1
  
  while(elec_created < num_elec){

    electrode_pts <- BiV_pts[as.character(electrode_vtx),]
    electrode_PHI <- COORDS_PHI[electrode_vtx+1]
    electrode_Z <- COORDS_Z[electrode_vtx+1]
    
    PHI_greater_than <-  electrode_PHI - 0.05
    PHI_minor_than <-  electrode_PHI + 0.05
    
    # We restrict to have be in the same PHI band 
    df_total <- df_original[df_original$COORDS_PHI >= PHI_greater_than,]
    df_total <- df_total[df_total$COORDS_PHI <= PHI_minor_than,]
    
    # We restrict to be in a lower position Z-wise
    df_total <- df_total[df_total$COORDS_Z < electrode_Z,]
    
    df_distances <- apply(df_total,1,function(x) (x-electrode_pts)^2 %>%
                                                  sum(.) %>% sqrt(.))
    
    df_total$distances <- df_distances
    
    tolerance <- 1000 # 1 mm
    
    df_total <- df_total[df_total$distances <= (distance + tolerance),]
    df_total <- df_total[df_total$distances >= (distance - tolerance),]
    
    
    electrode_vtx <- abs(df_total$COORDS_PHI - electrode_PHI) %>%
      which.min(.) %>%
      extract2(rownames(df_total),.) %>%
      as.numeric(.)
    
    electrodes_vector <- c(electrodes_vector,electrode_vtx)
    elec_created <- elec_created + 1
    
    if(action == "w1by1"){
      if(flag_debugging)
        Debug_message(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                             which_cases,"_case",heart,"/simulations/multipole/electrodes/"
                             ,vein,"_",elec_created,".vtx"),"w")
      write.table(c("1","intra",electrode_vtx),
                  paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                  which_cases,"_case",heart,"/simulations/multipole/electrodes/",
                  vein,"_",elec_created,".vtx"),quote = F,sep = "\n",
                  row.names = F,col.names = F)
    }
  }
  
  if(action == "return")
    return(electrodes_vector)
  else if(action == "writeall"){
    if(flag_debugging)
      Debug_message(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                           which_cases,"_case",heart,"/simulations/multipole/electrodes/"
                           ,vein,".vtx"),"w")
    write.table(c(num_elec,"intra",electrode_vtx),
                paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                       which_cases,"_case",heart,"/simulations/multipole/electrodes/",
                       vein,".vtx"),quote = F,sep = "\n",row.names = F,col.names = F)
  }
}

#' @description Function to create a file with the vertex in the border of AHA
#' zones 8 and 9, in the middle w.r.t. the apico-basal distance.
#' 
#' @param heart Number of the heart with two digits.
#' @param which_cases h/RR or HF for healthy/reverse remodelled or heart failure
#' respectively.
#' @param flag_debugging If TRUE, it prints whatever is reading and writing.

Create_RV_midseptum <- function(heart, which_cases, flag_debugging = FALSE){
  
  if(as.numeric(heart) < 10){
    heart <- paste0("0", as.numeric(heart))
  }
  
  if(which_cases == "RR" || which_cases == "h"){
    which_cases <- "h"
    BiV_dir <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                      which_cases,"_case",heart,"/meshing/1000um/BiV")
    UVC_dir <- paste0(BiV_dir,"/UVC/UVC")
  }
  if(which_cases == "HF"){
    BiV_dir <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                      which_cases,"_case",heart,"/meshing/1000um/BiV")
    UVC_dir <-  paste0(BiV_dir,"/UVC")
  }
  
  COORDS_Z <- Read_table(paste0(UVC_dir,"/COORDS_Z.dat"), header = FALSE,
                         flag_debugging = flag_debugging)
  COORDS_Z <- COORDS_Z[[1]]

  COORDS_PHI <- Read_table(paste0(UVC_dir,"/COORDS_PHI.dat"), header = FALSE,
                           flag_debugging = flag_debugging)
  COORDS_PHI <- COORDS_PHI[[1]]

  COORDS_RHO <- Read_table(paste0(UVC_dir,"/COORDS_RHO.dat"), header = FALSE,
                           flag_debugging = flag_debugging)
  COORDS_RHO <- COORDS_RHO[[1]]

  BiV_AHA <- Read_table(paste0(BiV_dir,"/BiV_AHA.dat"), header = FALSE)
  BiV_AHA <- BiV_AHA[[1]]


  AHA_9_idx <- which(BiV_AHA == 9)
  epi_idx <- which(COORDS_RHO == 1)
  epi_9_idx <- intersect(AHA_9_idx, epi_idx)

  target_Z <- COORDS_Z[epi_9_idx]
  target_PHI <- COORDS_PHI[epi_9_idx]

  max_PHI <- max(target_PHI)
  min_PHI <- min(target_PHI)
  range_PHI <- max_PHI - min_PHI

  goal_PHI <- max_PHI
  tol_PHI <- 0.2*range_PHI

  band_PHI_idx <- which(abs(target_PHI - goal_PHI) < tol_PHI)
  band_PHI_values <- target_PHI[band_PHI_idx]

  band_PHI_global_idx <- which(COORDS_PHI %in% band_PHI_values)

  band_PHI_Z_values <- COORDS_Z[band_PHI_global_idx]

  max_Z <- max(band_PHI_Z_values)
  min_Z <- min(band_PHI_Z_values)
  range_Z <- max_Z - min_Z
  goal_Z <- min_Z + range_Z/2
  tol_Z <- 0.2*range_Z

  band_Z_idx <- which(abs(target_Z - goal_Z) < tol_Z)
  band_Z_values <- target_Z[band_Z_idx]

  band_Z_global_idx <- which(COORDS_Z %in% band_Z_values)

  target_area <- intersect(band_PHI_global_idx,band_Z_global_idx)

  target_area_Z_values <- COORDS_Z[target_area]
  target_area_PHI_values <- COORDS_PHI[target_area]

  goal_function <- abs(target_area_Z_values - goal_Z) *
                   abs(target_area_PHI_values - goal_PHI)

  Z_miseptum_value <- target_area_Z_values[which.min(goal_function)]

  midseptum_idx <- which(COORDS_Z == Z_miseptum_value)
  
  Write.table(c("1", "intra", midseptum_idx - 1),
              paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                     which_cases,"_case",heart,"/simulations/multipole",
                     "/electrodes/BiV.midseptum.vtx"), quote = F, sep = "\n",
              row.names = F, col.names = F, flag_debugging = flag_debugging)
  
  
}

#' @description Function to convert an activation time file node-wise to
#' elem-wise. It's based on a template file because the backslashes in R give
#' too much problem. It takes around 6 sec per file. 
#' 
#' @param SA_folder Name of the sensitivity analysis folder (preffix for the 
#' files).
#' @param which_cases h or HF for healthy or heart failure.
#' @param heart Number of the heart case.
#' 
#' @return Returns the same file element-wise with the same name and the suffix
#' _elemwise.

AT_node2elem <- function(SA_folder,which_cases,heart){
  
  template <- readLines(paste0("/home/crg17/Desktop/scripts/multipole/sh/",
                               "AT_node2elem_template.sh"))
  output_file <- paste0("/home/crg17/Desktop/scripts/multipole/sh/",
                        "AT_node2elem_",SA_folder,"_",which_cases,"_",
                        heart,".sh")
  
  template[3] <- paste0("SA_folder=\"",SA_folder,"\"")
  template[4] <- paste0("which_cases=\"",which_cases,"\"")
  template[5] <- paste0("heart=\"",heart,"\"")
  
  writeLines(template,output_file)
  
  system(paste0("chmod +x ",output_file))
  
  system(output_file)
  
  system(paste0("rm ",output_file))
}

#' @description Function to convert the BiV apex activation time file node-wise
#' to elem-wise. It's based on a template file because the backslashes in R give
#' too much problem. It takes around 6 sec per file. 
#' 
#' @param SA_folder Name of the sensitivity analysis folder (preffix for the 
#' files).
#' @param which_cases h or HF for healthy or heart failure.
#' @param heart Number of the heart case.
#' 
#' @return Returns the same file element-wise with the same name and the suffix
#' _elemwise.
AT_node2elem_BiV <- function(SA_folder, which_cases, heart, midseptum = FALSE){
  
  template <- readLines(paste0("/home/crg17/Desktop/scripts/multipole/sh/",
                               "AT_node2elem_BiV_template.sh"))
  output_file <- paste0("/home/crg17/Desktop/scripts/multipole/sh/",
                        "AT_node2elem_BiV",SA_folder,"_",which_cases,"_",
                        heart,".sh")
  
  template[3] <- paste0("SA_folder=\"",SA_folder,"\"")
  template[4] <- paste0("which_cases=\"",which_cases,"\"")
  template[5] <- paste0("heart=\"",heart,"\"")
  template[6] <- paste0("midseptum=\"",tolower(midseptum),"\"")
  
  writeLines(template,output_file)
  
  system(paste0("chmod +x ",output_file))
  
  system(paste0("bash ",output_file))
  
  system(paste0("rm ",output_file))
}

#' @description Function to create a pseudo-UVC transmural coordinate. It sets
#' a value of 4 in the points of the epicardium and 0 otherwise. 
#' 
#' @param which_cases "RR", "HF" or "both" for reverse remodelled, heart
#' failure or both, respectively.
#' @param flag_debugging If TRUE, prints whatever is reading and writing.
Extract_epi_RHO <- function(which_cases, flag_debugging = FALSE){
  
  if(which_cases == "both"){
    Extract_epi_RHO(which_cases = "RR", flag_debugging = flag_debugging)
    Extract_epi_RHO(which_cases = "HF", flag_debugging = flag_debugging)
  }
  else{
    heart_list <- c(paste0("0",1:9),10:20)
  
    if(which_cases == "RR"){
      which_cases <- "h"
    }
    if(which_cases == "HF"){
      heart_list <- c(heart_list,21:24)
    }
    
    for(heart in heart_list){
    path2epi <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                       which_cases, "_case", heart, "/meshing/1000um/BiV")
    epiname <- "BiV.epi.surf.vtx"
    
    path2RHO <- paste0(path2epi,"/UVC")
    
    if(which_cases == "h"){
      path2RHO <- paste0(path2RHO,"/UVC")
    }
    
    RHOname <- "COORDS_RHO.dat"
    
    epi_vtx <- Read_table(file = paste0(path2epi,"/",epiname), quote = "\"", 
                          comment.char = "", skip = 2,
                          flag_debugging = flag_debugging)
    epi_vtx <- epi_vtx$V1
    epi_vtx <- epi_vtx + 1
    
    RHO_dat <- Read_table(file = paste0(path2RHO, "/", RHOname), quote = "\"",
                          comment.char = "", flag_debugging = flag_debugging)
    RHO_dat <- RHO_dat$V1
    
    RHO_epi <- RHO_dat
    RHO_epi[-epi_vtx] <- 0
    RHO_epi[epi_vtx] <- 4
    
    Write_table(x = RHO_epi, file = paste0(path2RHO,"/COORDS_RHO_epi.dat"),
                flag_debugging = flag_debugging, row.names = FALSE,
                col.names = FALSE)
  
    }
  }
  
}

#' @description Function to map the thickness from the epicardium (coming from
#' the matlab script scar_thinning.m) of a HF heart to the whole myocardium.
#' 
#' @param case_number Number of the HF heart.
#' @param flag_debugging If TRUE prints whatever is reading and writing.
#' 
#' @return Writes a .dat file with the pointwise thickness. 

Map_thickness_transmurally <- function(case_number, flag_debugging = FALSE){
  
  Load_Install_Packages(c("svMisc","dplyr"))
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  
  path2biv <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",
                      case_number,"/meshing/1000um/BiV")

  
  # Comes from the matlab script "scar_thinning.m"
  BiV_thickness <- Read_table(paste0(path2biv,"/BiV_thickness.dat"), quote="\"",
                              comment.char="", flag_debugging = flag_debugging)
  BiV_thickness <- as.vector(BiV_thickness$V1)
  
  UVC <- Read_table(paste0(path2biv,"/UVC/COMBINED_COORDS_Z_RHO_PHI_V.dat"), 
                    quote="\"", comment.char="",
                    flag_debugging = flag_debugging)
  colnames(UVC) <- c("Z","RHO","PHI","V")
  
  
  # We find the indices of the surface, so with possitive thickness
  mesh_idx <- c(1:length(BiV_thickness))
  # surface_idx <- which(BiV_thickness != -1)
  surface_idx <- which(BiV_thickness != -1)
  
  wherediditgetthickness <- c(1:length(BiV_thickness))
  
  p <- progress_estimated(length(mesh_idx[!(mesh_idx%in%surface_idx)]))
  for (i in mesh_idx[!(mesh_idx%in%surface_idx)]) {
    if(UVC[i,"V"] == -1){
      p$tick()
      p$print()
      diffPHI <- abs(UVC[surface_idx,"PHI"] - UVC[i,"PHI"])
      diffZ <- abs(UVC[surface_idx,"Z"] - UVC[i,"Z"])
      
      order_PHI <- order(diffPHI, decreasing = FALSE)
      order_Z <- order(diffZ, decreasing = FALSE)
      
      # We find the closest indices in the Z direction
      # closest_Z_idx <- order_Z[1:100]
      # closest_Z_idx <- order_Z[1:5]
      
      # Their diffPHI values, to see how close they are in terms of PHI
      # diffPHI_closest_Z <- diffPHI[closest_Z_idx]
      
      # mean_vec <- apply(cbind(diffPHI,diffZ),1,'mean')
      # criterion_vec <- apply(cbind(order_PHI,order_Z),1,'prod')
      
      # idx_in_surface <- which.min(mean_vec)[1]
      # idx_in_surface <- which.min(criterion_vec)[1]
      # 
      # We want the index where that diffPHI is minimum
      # CORRECT THE [1]
      # idx_in_surface <- which(abs(diffPHI - min(diffPHI_closest_Z)) < 1e-5)[1]
      idx_in_surface <- Find_first_common_idx(order_PHI, order_Z)
      
      idx_in_mesh <- surface_idx[idx_in_surface][1]
      
      BiV_thickness[i] <- BiV_thickness[idx_in_mesh]
      wherediditgetthickness[i] <- idx_in_mesh
      
    }
  }
  
  Write_table(BiV_thickness, file = paste0(path2biv,
                                          "/BiV_thickness_transmural_diag.dat"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              flag_debugging = flag_debugging)
  
  Write_table(wherediditgetthickness, file = paste0(path2biv,
                                           "/wherediditgetthickness_diag.dat"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              flag_debugging = flag_debugging)
}

Map_thickness_transmurally_pts <- function(case_number, flag_debugging = FALSE){
  
  Load_Install_Packages(c("svMisc","dplyr"))
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  
  path2biv <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",
                     case_number,"/meshing/1000um/BiV")
  
  
  # Comes from the matlab script "scar_thinning.m"
  BiV_thickness <- Read_table(paste0(path2biv,"/BiV_thickness.dat"), quote="\"",
                              comment.char="", flag_debugging = flag_debugging)
  BiV_thickness <- as.vector(BiV_thickness$V1)
  
  UVC <- Read_table(paste0(path2biv,"/UVC/COMBINED_COORDS_Z_RHO_PHI_V.dat"), 
                    quote="\"", comment.char="",
                    flag_debugging = flag_debugging)
  colnames(UVC) <- c("Z","RHO","PHI","V")
  
  BiV_pts <- Read_table(paste0(path2biv,"/meshes/BiV_FEC_w5_h33_retagged_noPVTV.pts"),
                        quote="\"",skip =1,
                        comment.char="", flag_debugging = flag_debugging)
  colnames(BiV_pts) <- c("x","y","z")
  # We find the indices of the surface, so with possitive thickness
  mesh_idx <- c(1:length(BiV_thickness))
  surface_idx <- which(UVC$RHO == 1)
  surface_coords <- BiV_pts[surface_idx,]
  
  i_lv <- which(UVC$V == -1)
  
  # i_epi = find(rho(i_lv)==1);
  # i_endo = find(rho(i_lv)==0);
  
  
  # h = -1*ones(size(pts,1),1);
  h <- BiV_thickness
  
  
  p <- progress_estimated(length(i_lv[!(i_lv%in%surface_idx)]))
  # for (i in mesh_idx[!(mesh_idx%in%surface_idx)]) {
  for(i in i_lv[!(i_lv%in%surface_idx)]){
    # if(UVC[i,"V"] == -1){
      p$tick()
      p$print()
      
      i_coord <- BiV_pts[i_lv[i],]
      
      # v1 = pts(i_lv(i),:);
      # v2 = pts(i_lv(i_endo),:);
      # x = v1 - v2;
      # dist = sqrt(sum(x.^2,2));
      dist <- apply(surface_coords, 1, function(x) dist(rbind(i_coord,x))) 
      h[i_lv[i]] = h[surface_idx[which.min(dist)]]
  }
  
  Write_table(BiV_thickness, file = paste0(path2biv,
                                           "/BiV_thickness_transmural_pts.dat"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              flag_debugging = flag_debugging)
}

#' @description Function to select the points which thickness is less than 5mm,
#' discarding the apex (AHA area 17) and the base (AHA area 0).
#' 
#' @param case_number Number of the HF heart.
#' @param flag_debugging If TRUE prints whatever is reading and writing.
#' 
#' @return A file with 1 if the point is scar tissue and 0 otherwise.
Select_scar_thickness <- function(case_number, scar_thickness_mm = 5, flag_debugging = FALSE){

  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  
  path2biv <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",
                      case_number,"/meshing/1000um/BiV")
  
  # From Map_thickness_transmurally
  BiV_thickness <- Read_table(paste0(path2biv,"/BiV_thickness_transmural_diag.dat"),
                              quote="\"", comment.char="",
                              flag_debugging = flag_debugging)
  BiV_thickness <- as.vector(BiV_thickness$V1)
  
  BiV_AHA <- Read_table(paste0(path2biv,"/BiV_AHA.dat"), quote="\"",
                        comment.char="", flag_debugging = flag_debugging)
  BiV_AHA <- as.vector(BiV_AHA$V1)
  
  BiV_septum <- Read_table(paste0(path2biv, "/BiV.rvsept.surf.vtx"), quote="\"",
                           comment.char="", flag_debugging = flag_debugging)
  ll <- dim(BiV_septum)[1]
  septum_vtx <- as.double(BiV_septum$V1[3:ll])
  
  thin_vtx <- which(BiV_thickness < 1000*scar_thickness_mm)
  
  noapex_vtx <- which(BiV_AHA < 17)
  nobase_vtx <- which(BiV_AHA > 0)
  
  chosen_vtx_septum <- Reduce(intersect, list(thin_vtx,noapex_vtx,nobase_vtx))
  chosen_vtx <- chosen_vtx_septum
  # chosen_vtx <- chosen_vtx_septum[!(chosen_vtx_septum %in% septum_vtx)]
  
  is_scar <- rep(0,length(BiV_AHA))
  is_scar[chosen_vtx] <- 1
  
  Write_table(is_scar,file = paste0(path2biv, "/BiV_scar_", toString(scar_thickness_mm), "mm.dat"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              flag_debugging = flag_debugging)
}

#' @description Function to change the tag file to include tag 100 wherever
#' there is scar. Comes from removing the header of the mesh with
#' tail -n +2 file.elem > file_noheader.elem. The BiV_scar_5mm.dat should
#' also be mapped to elemwise.
#' 
#' @param case_number Number of the HF heart.
#' @param meshname_noheader Name of the elem file without the first line to add
#' the new tag. The file has to end in "_noheader" but the input to the function
#' is the string right before.
#' @param flag_debugging If TRUE prints whatever is reading and writing.  
#' 
#' @return A file with the same name as the meshname_noheader but with _5mm_ in
#' the middle.
#' 
#' @note To add the header later:
#' filename=./BiV_FEC_w5_h33_retagged_noPVTV_scar5mm; 
#' cp $filename"_noheader.elem" $filename".elem";
#'  sed -i '1s/^/'$(wc -l < $filename".elem")'\n/' $filename".elem"
Change_scar_tag <- function(case_number, meshname_noheader, scar_thickness_mm = 5,
                            flag_debugging = FALSE){

  path2biv <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",
                     case_number,"/meshing/1000um/BiV")
  
  BiV_noheader <- Read_table(paste0(path2biv, "/meshes/", meshname_noheader,
                                    "_noheader.elem"),
                             quote="\"", comment.char="",
                             flag_debugging = flag_debugging)
  colnames(BiV_noheader) <- c("el","p1","p2","p3","p4","tag")
  
  BiV_scar <- Read_table(paste0(path2biv,"/BiV_scar_", toString(scar_thickness_mm),"mm_elem.dat"),
                         quote="\"", comment.char="",
                         flag_debugging = flag_debugging)
  
  BiV_noheader$tag[which(BiV_scar > 0)] <- 100
  
  Write_table(BiV_noheader, file = paste0(path2biv, "/meshes/", meshname_noheader,
                                          "_scar",toString(scar_thickness_mm),"mm_noheader.elem"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              flag_debugging = flag_debugging)
  
}

