#'  @param which_cases RR/h or HF for reverse remodelling or heart failure
#' respectively. "both" to run both cases.
#'  @param step "crop_base", "node2elem", "monodipoles", "monopoles", "dipoles",
#'  "AT_table","AT_table_RV".
#'  @param SA_folder Name of the sensitivity analysis folder.
#'  @param output_... Only relevant if step = AT_table. Prints only the results
#'  corresponding to the metrics which flags are set to TRUE.
#'  @param flag_debugging If TRUE, prints everything read and written.
Run_pipeline <- function(which_cases, step, SA_folder = "default_HF",
                         output_TAT = TRUE, output_AT1090 = FALSE,
                         output_AT090 = TRUE, output_LVTAT = FALSE,
                         output_VEUTAT = FALSE, output_VEUmean = FALSE,
                         midseptum = FALSE,
                         with_scar = FALSE,
                         flag_debugging = FALSE,
                         root_directory = "/media/crg17/Seagate Expansion Drive/SA_multipole/"){
  ptm <- proc.time()
  
  source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/common_functions.R")
  source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/preprocessing.R")
  source("/home/crg17/Desktop/KCL_projects/MPP/multipole/R/compute_EP.R")
  Load_Install_Packages(c("doParallel","dplyr"))

  if(which_cases == "RR"){
    which_cases <- "h"
  }
  
  if(which_cases == "h"){
    hearts <- c(paste0("0",1:9),10:20)
  }
  else if(which_cases == "HF"){
    hearts <- c(paste0("0",1:9),10:24)
    if(with_scar){
      hearts <- hearts[-c(13,21)]
    }
  }
  
  if(step == "crop_base"){
    if(which_cases == "both"){
      Run_pipeline("RR","crop_base",flag_debugging = flag_debugging)
      Run_pipeline("HF","crop_base",flag_debugging = flag_debugging)
    }
      else{

    registerDoParallel(cores=20)
    
    foreach(i=1:length(hearts)) %dopar% {
    # for(i in c(1:length(hearts))){
    Crop_base_from_tags(hearts[i],which_cases,"BiV_FEC_w5_h70_retagged", 
                        UVC_folder = "UVC_PM", output_suffix = "noPVTV",
                        flag_debugging = flag_debugging)
    }
    
    Send_mail(paste0("Bases cropped in ",which_cases,"."))
      }
  }
  if(step == "node2elem"){
    if(which_cases == "both"){
      Run_pipeline(which_cases = "h", step = step, SA_folder = SA_folder,
                   flag_debugging = flag_debugging, midseptum = midseptum)
      Run_pipeline(which_cases = "HF", step = step, SA_folder = SA_folder,
                   flag_debugging = flag_debugging, midseptum = midseptum)
    }
    else{
      registerDoParallel(cores=20)
      
      foreach(i=1:length(hearts)) %dopar% {
        AT_node2elem(SA_folder = SA_folder, which_cases = which_cases,
                     heart = hearts[i])
        AT_node2elem_BiV(SA_folder = SA_folder, which_cases = which_cases,
                         heart = hearts[i], midseptum = midseptum)
      }
      
      paste0("Mapped activation times from pointwise to elemwise in ",
             which_cases) %>% Send_mail(.)
    }
  }
  if(step == "monodipoles" || step == "monopoles"){
    CreateMonopoles(which_cases, SA_folder = SA_folder, midseptum = midseptum,
                    flag_debugging = flag_debugging, with_scar = with_scar)
    total_time <- proc.time() - ptm
    Send_mail(paste0("Monopoles created in ",which_cases," cases. Total time of",
                     total_time))
  }
  if(step == "monodipoles" || step == "dipoles"){
    
    if(which_cases != "both"){
    # registerDoParallel(cores=20)
    
    # foreach(i=1:length(hearts)) %dopar% {
    # for(i in 1:length(hearts)){
    CreateDipoles(SA_folder = SA_folder, which_cases = which_cases,
                  heart = "all", flag_debugging = flag_debugging, with_scar = with_scar)
    # }
    total_time <- proc.time() - ptm
    Send_mail(paste0("Dipoles created in ",which_cases," cases. Total time of",
                     total_time))
              }
    
    if(which_cases == "both"){
      CreateDipoles(SA_folder = SA_folder, which_cases = which_cases,
                    heart = "all",flag_debugging = flag_debugging)
      total_time <- proc.time() - ptm
      Send_mail(paste0("Dipoles created in all the cases. Total time of",
                       total_time))
    }
  }
  if(step == "AT_table" || step == "AT_table_RV"){
    if(which_cases == "both"){
      Run_pipeline(which_cases = "RR", step = step, SA_folder = SA_folder,
                   midseptum = midseptum,
                   output_TAT = output_TAT, output_AT1090 = output_AT1090,
                   output_AT090 = output_AT090, output_LVTAT = output_LVTAT,
                   output_VEUTAT = output_VEUTAT, 
                   output_VEUmean = output_VEUmean,
                   flag_debugging = flag_debugging,
                   with_scar = with_scar)
      Run_pipeline(which_cases = "HF", step = step, SA_folder = SA_folder,
                   midseptum = midseptum,
                   output_TAT = output_TAT, output_AT1090 = output_AT1090,
                   output_AT090 = output_AT090, output_LVTAT = output_LVTAT,
                   output_VEUTAT = output_VEUTAT, 
                   output_VEUmean = output_VEUmean,
                   flag_debugging = flag_debugging, with_scar = with_scar)
    }
    else{
      
      if(step == "AT_table" || step == "AT_table_RV"){
      
      if(step == "AT_table"){  
      registerDoParallel(cores=20)

      #foreach(i=1:length(hearts)) %dopar% {
      for(i in c(1:length(hearts))){
    Write_EP_files(SA_folder = SA_folder, which_cases, heart = hearts[i],
                   output_TAT = output_TAT,
                   output_AT1090 = output_AT1090, output_AT090 = output_AT090,
                   output_LVTAT = output_LVTAT, output_VEUTAT = output_VEUTAT,
                   output_VEUmean = output_VEUmean,
                   flag_debugging = flag_debugging,
                   with_scar = with_scar)
      }
      }
      
      Write_EP_files_RV(SA_folder = SA_folder, which_cases = which_cases,
                        with_scar = with_scar,
                        flag_debugging = flag_debugging,
                        midseptum = midseptum,
                        with_scar = with_scar)
      total_time <- proc.time() - ptm
      Send_mail(paste0("AT tables created in ", which_cases,". Total time of",
                       total_time))
    }
    }
  }
}