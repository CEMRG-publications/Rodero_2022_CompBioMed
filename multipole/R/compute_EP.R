#' @description Script to compute EP metrics as the TAT duration, the AT1090
#' and the AT090 for a biventricular mesh from the vector data.
#' 
#' @param BiV_volume Vector which components are the volumes of each element
#' @param AT_vec_elem Vector which components are the activation times averaged 
#' per element. You can use meshtool to generate it.
#' @param BiV_tags Vector with the tag of each element.
#' @param ventricle Ventricle where the EP metrics are computed from. Options 
#' are "LV", "RV" or "both". Default is both.
#' 
#' @return List which first element is the total activation time, the second is 
#' the AT1090 and the third is AT090. 
#' 
#' @details  Only functional for BiV activations
#' 

Compute_EP_metrics_BiV <- function(BiV_volume, AT_vec_elem, BiV_tags, 
                                   ventricle = "both"){
  
  # Retag the FEC layer
  BiV_tags[BiV_tags == 25] <- 1
  BiV_tags[BiV_tags == 26] <- 2
  
  # Extract only the ventricle you want
  if(ventricle == "LV"){
    BiV_volume <- BiV_volume[which(BiV_tags == 1)]
    AT_vec_elem <- AT_vec_elem[which(BiV_tags == 1)]
  }
  else if(ventricle == "RV"){
    BiV_volume <- BiV_volume[which(BiV_tags == 2)]
    AT_vec_elem <- AT_vec_elem[which(BiV_tags == 2)]
  }
  else if(ventricle == "both"){
    BiV_volume <- BiV_volume[c(which(BiV_tags == 1), which(BiV_tags == 2))]
    AT_vec_elem <- AT_vec_elem[c(which(BiV_tags == 1), which(BiV_tags == 2))]  
  }
  
  # For both ventricles
  
  TAT <- max(AT_vec_elem) - min(AT_vec_elem)
  

  # We reorder the volumes putting first the ones that activate before

  vol_sorted_by_AT <- BiV_volume[order(AT_vec_elem)]

  # In percentage

  vol_perc_sorted_by_AT <- 100*vol_sorted_by_AT/sum(vol_sorted_by_AT)

  idx_10_in_sorted <- which(cumsum(vol_perc_sorted_by_AT) >= 10)[1]
  idx_90_in_sorted <- (which(cumsum(vol_perc_sorted_by_AT) >= 90)[1])-1

  AT1090 <- sort(AT_vec_elem)[idx_90_in_sorted] - 
            sort(AT_vec_elem)[idx_10_in_sorted]
  
  AT090 <- sort(AT_vec_elem)[idx_90_in_sorted] - min(AT_vec_elem)

  return(list(TAT,AT1090,AT090))
  
}

#' @description Script to compute EP metrics as the TAT duration and the AT1090 
#' for a biventricular mesh reading from file
#' 
#' @param path2volumefile Path to the .dat volume file. Without the last slash.
#' @param volumefilename File with the volumes of each element. Include the 
#' extension in the name.
#' @param path2ATfile Path to the file with the activation times averaged per 
#' element. You can use meshtool to generate it.
#' @param ATfilename Activation file name with the activation times averaged per
#'  element. You can use meshtool to generate it.
#' @param path2tagfile Path to the file with the tags file of the mesh. Can be
#'  extracted with meshtool extract tags.
#' @param tagfilename Name of the tags file name
#' @param ventricle Ventricle where the EP metrics are computed from. Options 
#' are "LV", "RV" or "both". Default is both.
#' @param flag_debugging If TRUE, it prints whatever is reading and writing.
#' 
#' @return Same output as Compute_EP_metrics_BiV. 
#' 
#' @details Only functional for BiV activations
#' 
Compute_EP_metrics_BiV_from_file <- function(path2volumefile,volumefilename,
                                             path2ATfile,ATfilename,
                                             path2tagfile,tagfilename,
                                             ventricle = "both", 
                                             flag_debugging=FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr"))
  
  
  BiV_volume <- paste0(path2volumefile,"/",volumefilename) %>%
                Read_table(., flag_debugging = flag_debugging)
  BiV_volume <- BiV_volume$V1
  
  
  AT_vec_elem <- paste0(path2ATfile,"/",ATfilename) %>%
                 Read_table(., flag_debugging = flag_debugging)
  AT_vec_elem <- AT_vec_elem$V1
  
  BiV_tags <- paste0(path2tagfile,"/",tagfilename) %>%
              Read_table(., flag_debugging = flag_debugging)
  BiV_tags <- BiV_tags$V1
  
  res_list <- Compute_EP_metrics_BiV(BiV_volume = BiV_volume,
                                     AT_vec_elem = AT_vec_elem,
                                     BiV_tags = BiV_tags,
                                     ventricle = ventricle)
  
  return(res_list)

}


#' @description Script to compute the synchronisation metrics of a BiV mesh from
#'  vector data.
#' 
#' @param AT_vec_elem Vector which components are the activation times averaged 
#' per element. You can use meshtool to generate it.
#' @param BiV_tags Vector with the tag of each element.
#' @param UVC_RHO_elem Transmural UVC coordinate interpolated on the elements.
#' 
#' @return List which first element is the total ventricular electrical 
#' uncoupling (VEU), the second is the mean VEU and the third is the mean VEU
#' but taking into account only the ATs in the epicardium. 
#' 
#' @details  Only functional for BiV activations
Compute_VEU <- function(AT_vec_elem, BiV_tags, UVC_RHO_elem){
  
  # AT_vec_elem <- AT_vec_elem[which(AT_vec_elem < 10e3)]
  # BiV_tags <- BiV_tags[which(AT_vec_elem < 10e3)]
  
  BiV_tags[BiV_tags == 25] <- 1
  BiV_tags[BiV_tags == 26] <- 2
  
  AT_LV <- AT_vec_elem[which(BiV_tags == 1)]
  AT_RV <- AT_vec_elem[which(BiV_tags == 2)]
  UVC_RHO_LV <- UVC_RHO_elem[which(BiV_tags == 1)]
  UVC_RHO_RV <- UVC_RHO_elem[which(BiV_tags == 2)]
  
  AT_LV_epi <- AT_LV[which(UVC_RHO_LV > 0)]
  AT_RV_epi <- AT_RV[which(UVC_RHO_RV > 0)]
  
  LV_TAT <- max(AT_LV) - min(AT_LV)
  RV_TAT <- max(AT_RV) - min(AT_RV)
  
  LV_TAT_epi <- max(AT_LV_epi) - min(AT_LV_epi)
  RV_TAT_epi <- max(AT_RV_epi) - min(AT_RV_epi)
  
  LV_mean <- mean(AT_LV)
  RV_mean <- mean(AT_RV)
  
  LV_mean_epi <- mean(AT_LV_epi)
  RV_mean_epi <- mean(AT_RV_epi)
  
  VEU_TAT <- LV_TAT - RV_TAT
  VEU_mean <- LV_mean - RV_mean
  VEU_epi <- LV_mean_epi - RV_mean_epi
  
  return(list(VEUTAT = VEU_TAT, VEUmean = VEU_mean, VEUepi = VEU_epi,
              LVTAT = LV_TAT, LVTATepi = LV_TAT_epi, LVmean = LV_mean,
              LVmeanepi = LV_mean_epi, RVTAT = RV_TAT, RVTATepi = RV_TAT_epi,
              RVmean = RV_mean, RVmeanepi = RV_mean_epi))
}

#' @description Script to compute the synchronisation metrics of a BiV mesh from
#'  a file.
#' 
#' @param path2ATfile Path to the file with the activation times averaged per 
#' element. You can use meshtool to generate it.
#' @param ATfilename Activation file name with the activation times averaged per
#'  element. You can use meshtool to generate it.
#' @param path2tagfile Path to the file with the tags file of the mesh. Can be 
#' extracted with meshtool extract tags.
#' @param tagfilename Name of the tags file name
#' @param flag_debugging If false, prints whatever is reading and writing.
#' 
#' @return List which first element is the total ventricular electrical 
#' uncoupling (VEU) and the second is the mean VEU. 
#' 
#' @details  Only functional for BiV activations
#' 
Compute_VEU_from_file <- function(path2ATfile,ATfilename,path2tagfile,
                                  tagfilename, path2RHO, RHOname,
                                  flag_debugging = FALSE){
  
  AT_vec_elem <- Read_table(paste0(path2ATfile,"/",ATfilename), 
                            flag_debugging = flag_debugging)
  AT_vec_elem <- AT_vec_elem$V1
  
  BiV_tags <- Read_table(paste0(path2tagfile,"/",tagfilename),
                         flag_debugging = flag_debugging)
  BiV_tags <- BiV_tags$V1
  
  UVC_RHO_elem <- Read_table(paste0(path2RHO, "/", RHOname),
                             flag_debugging = flag_debugging)
  UVC_RHO_elem <- UVC_RHO_elem$V1
  
  result <- Compute_VEU(AT_vec_elem = AT_vec_elem, BiV_tags = BiV_tags,
                        UVC_RHO_elem = UVC_RHO_elem)
  
  return(result)
}

#' @description Script to show the average +- SD of the QRS duration of the BiV.RV_endo.apex folders.
#' 
#' @param CV_folder Name of the folder of the simulations 
#' (with the eikonal_preffix). If ekbatch=TRUE, it corresponds to the preffix
#' of the .dat files in the init-files folder. 
#' @param which_cases Healthy ("h" or "RR") or heart failure ("HF").
#' @param flag_debugging Boolean if TRUE prints whatever is reading and writing.
#' @param ekbatch Boolean. If TRUE uses the init-file folder.
#' @return Prints the average QRS duration +- SD. 

Print_QRS <- function(CV_folder,which_cases="HF", flag_debugging = FALSE,
                      ekbatch=TRUE){
  if(which_cases == "RR"){
    which_cases <- "h"
  }
  
  if(which_cases == "h"){
    num_cases <- 20
  }
  else if(which_cases == "HF")
    num_cases <- 24
  
  QRS_vec <- c()
  for(heart in c(paste0("0",1:9),10:num_cases)){
    if(!ekbatch){
      vm_act_seq <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/simulations/multipole/",CV_folder,"/BiV.RV_endo.apex","/vm_act_seq.dat"),
                    header = FALSE)
    }
    else{
      vm_act_seq <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart,"/simulations/multipole/init_files/",CV_folder,"_BiV.RV_endo.apex.dat"),
                               header = FALSE)
    }
    QRS_vec <- c(QRS_vec,max(vm_act_seq[vm_act_seq < 10e3])-min(vm_act_seq))
  }
  print(sprintf("%i \U00B1 %i",round(mean(QRS_vec)),round(sd(QRS_vec))))
}

#' @description Function to create the files with tables specified in the 
#' output option from the monopoles and dipoles files. 
#' 
#' @param SA_folder Name of the input/output parent folder.
#' @param which_cases RR/h, HF or both for reverse remodelled/healthy, heart
#' failure or both, respectively. 
#' @param heart Number of the heart. If "all", it run all of the cases (20 for
#' RR, 24 for HF). When heart=="all" it seems that it does not work. Run the 
#' loop outside this function.
#' @param with_scar If TRUE, it reads the BiV_tags corresponding to the case 
#' with scar.
#' @param output... If TRUE, it writes those output metrics.
#' @param flag_debugging Boolean if TRUE prints whatever is reading and writing.
#' 
#' @return As many files as TRUE flags are set in the output option are created
#' called multipole_metric.dat in the input folder with the corresponding 
#' metrics.
#' 
#' @note Not working if which_cases == all, does well RR cases but not HF. Run
#' them separately.
Write_EP_files <- function(SA_folder,which_cases,heart,
                           output_TAT = TRUE, output_AT1090 = TRUE,
                           output_AT090 = TRUE, output_VEUTAT = TRUE,
                           output_VEUmean = TRUE, output_LVTAT = TRUE, 
                           flag_debugging = FALSE,
                           root_directory = "/media/crg17/Seagate Expansion Drive/SA_multipole/",
                           with_scar = F){
  
  source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/common_functions.R")
  Load_Install_Packages("doParallel")
  
  if(which_cases == "both"){
    Write_EP_files(SA_folder = SA_folder, which_cases = "RR", heart = heart,
                   output_TAT = output_TAT,
                   output_AT1090 = output_AT1090, output_AT090 = output_AT090,
                   output_VEUTAT = output_VEUTAT, 
                   output_VEUmean = output_VEUmean, output_LVTAT = output_LVTAT,
                   flag_debugging = flag_debugging, root_directory = root_directory,
                   with_scar = with_scar)
    Write_EP_files(SA_folder = SA_folder, which_cases = "HF", heart = heart,
                   output_TAT = output_TAT,
                   output_AT1090 = output_AT1090, output_AT090 = output_AT090,
                   output_VEUTAT = output_VEUTAT,
                   output_VEUmean = output_VEUmean, output_LVTAT = output_LVTAT,
                   flag_debugging = flag_debugging, root_directory = root_directory,
                   with_scar = with_scar)
  }
    else{
      if(heart=="all"){
        heart_list <- c(paste0("0",1:9),10:20)
        if(which_cases=="HF"){
          heart_list <- c(heart_list,21:24)
          if(with_scar){
            heart_list <- heart_list[-c(13,21)]
          }
        }
        registerDoParallel(cores=20)

        foreach(i=1:length(heart_list)) %dopar% {

          Write_EP_files(SA_folder = SA_folder, which_cases = "RR",
                         heart = heart_list[i], 
                         output_TAT = output_TAT, output_AT1090 = output_AT1090,
                         output_AT090 = output_AT090,
                         output_VEUTAT = output_VEUTAT, 
                         output_VEUmean = output_VEUmean,
                         output_LVTAT = output_LVTAT,
                         flag_debugging = flag_debugging,
                         root_directory = root_directory,
                         with_scar = with_scar)
        }
      }
      
      
    if(which_cases == "RR"){
      which_cases <- "h"
    }
    
    path2biv <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                       which_cases,"_case",heart,"/meshing/1000um/BiV")
    

    tagfilename <- paste0("BiV_tags_",SA_folder,".dat")
    
    we_can_skip <- T
    
    if(output_TAT && !(file.exists(paste0(root_directory,SA_folder,"/",
                                          which_cases,"/",heart,
                                          "/multipole_TAT.dat")))){
      we_can_skip <- F
    }
    if(output_AT1090 && !(file.exists(paste0(root_directory,SA_folder,"/",
                                          which_cases,"/",heart,
                                          "/multipole_AT1090.dat")))){
      we_can_skip <- F
    }
    if(output_AT090 && !(file.exists(paste0(root_directory,SA_folder,"/",
                                          which_cases,"/",heart,
                                          "/multipole_AT090.dat")))){
      we_can_skip <- F
    }
    if(output_VEUTAT && !(file.exists(paste0(root_directory,SA_folder,"/",
                                            which_cases,"/",heart,
                                            "/multipole_VEUTAT.dat")))){
      we_can_skip <- F
    }
    if(output_VEUmean && !(file.exists(paste0(root_directory,SA_folder,"/",
                                            which_cases,"/",heart,
                                            "/multipole_VEUmean.dat")))){
      we_can_skip <- F
    }
    if(output_LVTAT && !(file.exists(paste0(root_directory,SA_folder,"/",
                                            which_cases,"/",heart,
                                            "/multipole_LVTAT.dat")))){
      we_can_skip <- F
    }
    
    if(!we_can_skip){
      
      if(output_TAT || output_AT1090 || output_LVTAT || output_AT090){
        biv_volume <- paste0(path2biv,"/BiV_mesh_volume.dat") %>%
                      Read_table(file = ., flag_debugging = flag_debugging)
        biv_volume <- biv_volume$V1
      }
  
      biv_tags <- paste0(path2biv,"/",tagfilename) %>%
                  Read_table(file = ., flag_debugging = flag_debugging)
      biv_tags <- biv_tags$V1
  
  
  
      file_list <- list.files(path=paste0(root_directory,SA_folder,"/",
                                          which_cases,"/",heart))
      file_list <- file_list[which(toupper(substr(file_list,0,1)) == "H")]
  
      dataset_TAT <- as.data.frame(matrix(nrow = 0, ncol = 6))
      colnames(dataset_TAT) <- c("lead","AN","AL","LA","IL","IN")
      
        dataset_AT1090 <- dataset_TAT
        dataset_AT090 <- dataset_TAT
        
        if(output_LVTAT){
          dataset_LVTAT <- dataset_TAT
        }
      
      if(output_VEUTAT || output_VEUmean){
        dataset_VEUTAT <- dataset_TAT
        dataset_VEUmean <- dataset_TAT
      }
      for (i in 1:length(file_list)){
  
        vm_act_seq <- paste0(root_directory,SA_folder,"/",which_cases
        ,"/", heart,"/",file_list[i]) %>%
                      Read_table(., flag_debugging = flag_debugging)
        vm_act_seq <- vm_act_seq[[1]]
  
        lead_name <- file_list[i] %>% substr(.,nchar(.)-11,nchar(.)-4)
        vein_name <- file_list[i] %>% substr(.,nchar(.)-14,nchar(.)-13)
  
        # we check if that lead is already in the dataset
        rowtoinsert <- which(dataset_TAT$lead == lead_name)
  
        if(length(rowtoinsert) == 0){
          dataset_TAT[nrow(dataset_TAT)+1,1] <- lead_name
          dataset_AT1090[nrow(dataset_AT1090)+1,1] <- lead_name
          dataset_AT090[nrow(dataset_AT090)+1,1] <- lead_name
          
          if(output_VEUTAT || output_VEUmean){
            dataset_VEUTAT[nrow(dataset_VEUTAT)+1,1] <- lead_name
            dataset_VEUmean[nrow(dataset_VEUmean)+1,1] <- lead_name
          }
          
          if(output_LVTAT){
            dataset_LVTAT[nrow(dataset_LVTAT)+1,1] <- lead_name
          }
          rowtoinsert <- nrow(dataset_TAT)
        }
        
        if(output_TAT || output_AT1090 || output_AT090){
          res_list <- Compute_EP_metrics_BiV(biv_volume,vm_act_seq,
                                           biv_tags,ventricle = "both")
        }
        if(output_TAT){
          dataset_TAT[rowtoinsert,vein_name] <- res_list[[1]]
        }
        if(output_AT1090){
          dataset_AT1090[rowtoinsert,vein_name] <- res_list[[2]]
        }
        if(output_AT090){
          dataset_AT090[rowtoinsert,vein_name] <- res_list[[3]]
        }
        if(output_LVTAT){
          res_list <- Compute_EP_metrics_BiV(biv_volume,vm_act_seq,
                                             biv_tags,ventricle = "LV")
          dataset_LVTAT[rowtoinsert,vein_name] <- res_list[[1]]
        }
        
        if(output_VEUTAT || output_VEUmean){
          res_list <- Compute_VEU(AT_vec_elem = vm_act_seq, BiV_tags = biv_tags)
        }
        if(output_VEUTAT){
          dataset_VEUTAT[rowtoinsert,vein_name] <- res_list[[1]]
        }
        if(output_VEUmean){
          dataset_VEUmean[rowtoinsert,vein_name] <- res_list[[2]]
        }
      }
      
      if(output_TAT){
        paste0(root_directory,SA_folder,"/",which_cases,"/",heart,
               "/multipole_TAT.dat") %>%
        Write_table(dataset_TAT,.,col.names = TRUE,row.names = FALSE,
                    quote = FALSE, flag_debugging = flag_debugging)
      }
      if(output_AT1090){
        paste0(root_directory,SA_folder,"/",which_cases,"/",heart,
               "/multipole_AT1090.dat") %>%
        Write_table(dataset_AT1090,.,col.names = TRUE,row.names = FALSE,
                    quote = FALSE, flag_debugging = flag_debugging)
      }
      if(output_AT090){
        paste0(root_directory,SA_folder,"/",which_cases,"/",heart,
               "/multipole_AT090.dat") %>%
          Write_table(dataset_AT090,.,col.names = TRUE,row.names = FALSE,
                      quote = FALSE, flag_debugging = flag_debugging)
      }
      if(output_LVTAT){
        paste0(root_directory,SA_folder,"/",which_cases,"/",heart,
               "/multipole_LVTAT.dat") %>%
          Write_table(dataset_LVTAT,.,col.names = TRUE,row.names = FALSE,
                      quote = FALSE, flag_debugging = flag_debugging)
      }
      if(output_VEUTAT){
        paste0(root_directory,SA_folder,"/",which_cases,"/",heart,
               "/multipole_VEUTAT.dat") %>%
        Write_table(dataset_VEUTAT,.,col.names = TRUE,row.names = FALSE,
                    quote = FALSE, flag_debugging = flag_debugging)
      }
      if(output_VEUmean){
        paste0(root_directory,SA_folder,"/",which_cases,"/",heart,
               "/multipole_VEUmean.dat") %>%
        Write_table(dataset_VEUmean,.,col.names = TRUE,row.names = FALSE,
                      quote = FALSE, flag_debugging = flag_debugging)
      }
    }
    }
}

#' @description Function to create the files with values of the TAT duration, 
#' AT1090, AT090, VEU TAT, VEU mean and LVTAT for the RV apex activation. 
#' 
#' @param SA_folder Name of the input/output parent folder.
#' @param which_cases RR/h, HF or both for reverse remodelled/healthy, heart
#' failure or both, respectively. 
#' @param with_scar If TRUE, it reads the BiV_tags corresponding to the case with 
#' scar.
#' @param midseptum If TRUE reads the midseptum apex.
#' @param flag_debugging Boolean if TRUE prints whatever is reading and writing.
#' 
#' @return A file called multipole_RVapex.dat in the parent folder of all the
#' cases whose columns are each one of the metrics.
Write_EP_files_RV <- function(SA_folder,which_cases,
                              midseptum = FALSE,
                              with_scar = FALSE,
                              flag_debugging = FALSE,
                              root_directory = "/media/crg17/Seagate Expansion Drive/SA_multipole/"){
  
  if(which_cases == "both"){
    Write_EP_files_RV(SA_folder = SA_folder, which_cases = "RR",
                      with_scar = with_scar,
                      flag_debugging = flag_debugging,
                      root_directory = root_directory)
    Write_EP_files_RV(SA_folder = SA_folder, which_cases = "HF",
                      with_scar = with_scar,
                      flag_debugging = flag_debugging,
                      root_directory = root_directory)
  }
  
  else{
  if(which_cases == "RR")
    which_cases <- "h"
  
  heart_list <- c(paste0("0",1:9),10:20)
  if(which_cases == "HF"){
    heart_list <- c(heart_list,21:24)
    if(with_scar){
      heart_list <- heart_list[-c(13,21)]
    }
  }
  
  multipole_RV <- as.data.frame(matrix(nrow = length(heart_list), ncol = 14))
  rownames(multipole_RV) <- heart_list
  colnames(multipole_RV) <- c("TAT", "AT1090", "AT090", "VEUTAT", "VEUmean",
                              "VEUepi", "LVTAT", "LVTATepi", "LVmean", 
                              "LVmeanepi", "RVTAT", "RVTATepi", "RVmean", 
                              "RVmeanepi")
  
  volumefilename <- "BiV_mesh_volume.dat"
  if(!midseptum){
    ATfilename <- paste0(SA_folder,"_BiV.RV_endo.apex_elemwise.dat")
  }
  else{
    ATfilename <- paste0(SA_folder,"_BiV.midseptum_elemwise.dat")
  }
  RHOname <- "COORDS_RHO_epi_elemwise.dat"
  
  tagfilename <- paste0("BiV_tags_",SA_folder,".dat")

  
  for(heart in heart_list){
    path2volumefile <- paste0("/media/crg17/Seagate Backup Plus Drive/",
                              "CT_cases/",which_cases,"_case",heart,"/meshing/",
                              "1000um/BiV")
    path2ATfile <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",
                          which_cases,"_case",heart,"/simulations/multipole/",
                          "init_files")
    path2RHO <- paste0(path2volumefile,"/UVC")
    
    if(which_cases == "h"){
      path2RHO <- paste0(path2RHO,"/UVC")
    }

    
    res_list <- Compute_EP_metrics_BiV_from_file(path2volumefile = path2volumefile,
                                                 volumefilename = volumefilename,
                                                 path2ATfile = path2ATfile,
                                                 ATfilename = ATfilename,
                                                 path2tagfile = path2volumefile,
                                                 tagfilename = tagfilename,
                                                 ventricle = "both",
                                                 flag_debugging = flag_debugging)
    multipole_RV[heart, "TAT"] <- res_list[[1]]
    multipole_RV[heart, "AT1090"] <- res_list[[2]]
    multipole_RV[heart, "AT090"] <- res_list[[3]]
    
    res_list <- Compute_VEU_from_file(path2ATfile = path2ATfile, 
                                      ATfilename = ATfilename,
                                      path2tagfile = path2volumefile,
                                      tagfilename = tagfilename,
                                      path2RHO = path2RHO,
                                      RHOname = RHOname,
                                      flag_debugging = flag_debugging)
    
    multipole_RV[heart, "VEUTAT"] <- res_list$VEUTAT
    multipole_RV[heart, "VEUmean"] <- res_list$VEUmean
    multipole_RV[heart, "VEUepi"] <- res_list$VEUepi
    multipole_RV[heart, "LVTAT"] <- res_list$LVTAT
    multipole_RV[heart, "LVTATepi"] <- res_list$LVTATepi
    multipole_RV[heart, "LVmean"] <- res_list$LVmean
    multipole_RV[heart, "LVmeanepi"] <- res_list$LVmeanepi
    multipole_RV[heart, "RVTAT"] <- res_list$RVTAT
    multipole_RV[heart, "RVTATepi"] <- res_list$RVTATepi
    multipole_RV[heart, "RVmean"] <- res_list$RVmean
    multipole_RV[heart, "RVmeanepi"] <- res_list$RVmeanepi
   
  }
  
  paste0(root_directory,SA_folder,"/",which_cases,
                       "/multipole_RVapex.dat") %>%
  Write_table(multipole_RV,file = .,quote = FALSE,col.names = TRUE,
              row.names = TRUE)
  }
  }


Visualise_AT1090 <- function(heart){
  
  
  
  warning("The elem file must be without the header and the activation file must be element-wise (use meshtool)")
  
  aux <- read.table("~/Desktop/transfer/vm_act_seq_11HC.dat", quote="\"", comment.char="")
  AT_12 <- as.vector(aux$V1)
  rm(aux)
  aux <- read.table("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case05/meshing/1000um/BiV/BiV_mesh_volume.dat", quote="\"", comment.char="")
  BiV_mesh_volume <- as.vector(aux$V1)
  
  # We sort the volumes according to the ones activated first
  
  idx2sort <- order(AT_12)
  vol_sorted <- BiV_mesh_volume[idx2sort]
  
  # We get the percentage of the total volume that it's each element activated
  
  tot_vol = sum(vol_sorted)
  
  accumul_vol = 100*cumsum(vol_sorted)/tot_vol
  
  # We recover the original sorting
  
  idx2unsort <- match(AT_12,sort(AT_12))
  
  final_vec <- accumul_vol[idx2unsort]
  
  write.table(final_vec,"/home/crg17/Desktop/transfer/HF05_percentages_AT.dat",row.names = FALSE,col.names = FALSE)
}

