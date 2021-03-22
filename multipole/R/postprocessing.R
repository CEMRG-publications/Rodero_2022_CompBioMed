
#' @description Function to find the optimal lead configuration give the AT
#' files for the whole cohort.
#' 
#' @param cohort "RR","HF" or "RRHF" for reverse remodelled, heart failure or
#' both respectively.
#' @param optimal_idx If 1, we get the configuration that its optimal for the
#' maximum number of patients; if 2, it's the second configuration and so on.
#' @param metric_option "TAT", "AT1090", "AT090", "VEUTAT", "VEUmean", "LVTAT",
#' for each of the metrics. 
#' @param SA_folder Name of the folder of the simulations.
#' @param flag_debugging Flag (TRUE or FALSE) to print the name and location of
#' the files read and written.
#' 
#' @return The best quadripole for the whole cohort.

Print_Global_Optimal_AT <- function(cohort = "HF", optimal_idx = 1,
                                    metric_option, SA_folder = "FEC_70_values", 
                                    vein = "ALL", flag_debugging = FALSE){
  
  

  # We find the optimal quadripoles
  optimal_data <- Find_Optimal_Quadripole_Optimising(metric_option = metric_option,
                                                     vein = vein,
                                                     which_cases = cohort,
                                                     method = "max", 
                                                     SA_folder = SA_folder,
                                                     version = 4,
                                                     flag_debugging = flag_debugging)
  # We extract the first one (optimal_idx)
  optimal_ascii <- optimal_data$Quadripolar[optimal_idx]
  
  # We extract all the possible bipoles in binary
  cyphers <- strsplit(optimal_ascii,"")[[1]]
  combs <- combn(4,2)
  vec_bipoles_from_quadpoles <- c()
  for(i in c(1:ncol(combs))){
    vec_bipoles_from_quadpoles[length(vec_bipoles_from_quadpoles) + 1] <-
    paste0(cyphers[combs[1,i]],cyphers[combs[2,i]])
  }
  bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
  
  if(vein == "ALL"){
  # We read the table with the AT for this specific vein
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "AN",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder,
                         normalising = FALSE, flag_debugging = flag_debugging)
  
  # We keep 
  df_final <- aux[[1]][bipoles_bin,]
  if(cohort == "RRHF"){
    df_final <- cbind(df_final,aux[[2]][bipoles_bin,])
  }
  rownames(df_final) <- paste0("x",rownames(df_final))
  
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "AL",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder, 
                         normalising = FALSE, flag_debugging = flag_debugging)
  first_df <- aux[[1]][bipoles_bin,]
  if(cohort == "RRHF"){
    first_df <- cbind(first_df,aux[[2]][bipoles_bin,])
  }
  df_final <- rbind(df_final,first_df)
  
  rownames(df_final) <- paste0("x",rownames(df_final))
  
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "LA",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder,
                         normalising = FALSE, flag_debugging = flag_debugging)
  first_df <- aux[[1]][bipoles_bin,]
  if(cohort == "RRHF")
    first_df <- cbind(first_df,aux[[2]][bipoles_bin,])
  df_final <- rbind(df_final,first_df)
  rownames(df_final) <- paste0("x",rownames(df_final))
  
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "IL",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder, 
                         normalising = FALSE, flag_debugging = flag_debugging)
  first_df <- aux[[1]][bipoles_bin,]
  if(cohort == "RRHF")
    first_df <- cbind(first_df,aux[[2]][bipoles_bin,])
  df_final <- rbind(df_final,first_df)
  rownames(df_final) <- paste0("x",rownames(df_final))
  
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "IN",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder,
                         normalising = FALSE, flag_debugging = flag_debugging)
  first_df <- aux[[1]][bipoles_bin,]
  if(cohort == "RRHF")
    first_df <- cbind(first_df,aux[[2]][bipoles_bin,])
  df_final <- rbind(df_final,first_df)
  rownames(df_final) <- paste0("x",rownames(df_final))
  }
  else{
    # We read the table with the AT for this specific vein
    aux <- Read_dataframes(metric_option = metric_option, lead_option = vein,
                           which_cases = cohort, version = 4,
                           SA_folder = SA_folder, 
                           normalising = FALSE, flag_debugging = flag_debugging)
    
    # We keep 
    df_final <- aux[[1]][bipoles_bin,]
    if(cohort == "RRHF"){
      df_final <- cbind(df_final,aux[[2]][bipoles_bin,])
    }
    rownames(df_final) <- paste0("x",rownames(df_final))
  }
  
  
  
  final_result <- mapply(min,df_final) %>% as.data.frame(.)
  colnames(final_result) <- optimal_ascii
  
  return(final_result)
}

#' @description Finds the optimal bipoles from the AT_table, without clustering
#' them into quadripoles.
#' 
#' @param cohort "RR", "HF" or "RRHF" for each one of the cohorts or altogether.
#' @param metric_option "TAT", "AT1090", "AT090", "VEUTAT", "VEUmean", "LVTAT",
#' for each of the metrics.
#' @param SA_folder Folder name of the sensitivity analysis case.
#' @param vein Name of the vein where we compute the optimal. If "ALL", it 
#' computes the optimal value using all the veins at once.
#' @param flag_debugging If TRUE, prints whatever is reading.
#' 
#' @return List whose first component is the optimal activation times for each
#' case and second component is the name of the bipolar lead.

Print_Individual_Optimal_AT <- function(cohort, metric_option,
                                        SA_folder, vein = "ALL",
                                        flag_debugging = FALSE){
  

  if(vein == "ALL"){
  
  #----- Anterior
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "AN",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder,
                           normalising = FALSE, flag_debugging = flag_debugging)
  if(cohort == "RRHF"){
    bipoles_df <- cbind(aux[[1]],aux[[2]])
  }
  else{
    bipoles_df <- aux[[1]]
  }
  
  min_AT <- mapply(min,bipoles_df)
  
  min_names <- rownames(bipoles_df)[mapply(which.min, bipoles_df)]
  names(min_names) <- names(min_AT)

  #----- Anterolateral
  
  aux <- Read_dataframes(metric_option = metric_option,
                         lead_option = "AL", which_cases = cohort, version = 4,
                         SA_folder = SA_folder,
                         normalising = FALSE, flag_debugging = flag_debugging)
  
  if(cohort == "RRHF"){
    bipoles_df <- cbind(aux[[1]],aux[[2]])
  }
  else{
    bipoles_df <- aux[[1]]
  }
  
  min_aux <- mapply(min,bipoles_df)
  min_names_aux <- rownames(bipoles_df)[mapply(which.min,bipoles_df)]
  
  # We check if we need to update or not
  old_AT_new_AT <- min_AT <= min_aux 
  
  for(i in c(1:length(old_AT_new_AT))){
    if(!old_AT_new_AT[i]){
      min_names[i] <- min_names_aux[i]
    }
  }
  
  min_AT <- apply(cbind(min_AT,min_aux),1,min)
  
  #----- Lateral
  
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "LA",
                         which_cases = cohort, version = 4, 
                         SA_folder = SA_folder, 
                         normalising = FALSE, flag_debugging = flag_debugging)
  
  if(cohort == "RRHF"){
    bipoles_df <- cbind(aux[[1]],aux[[2]])
  }
  else{
    bipoles_df <- aux[[1]]
  }
  
  min_aux <- mapply(min,bipoles_df)
  
  min_names_aux <- rownames(bipoles_df)[mapply(which.min,bipoles_df)]
  
  old_AT_new_AT <- min_AT <= min_aux
  
  for(i in c(1:length(old_AT_new_AT))){
    if(!old_AT_new_AT[i]){
      min_names[i] <- min_names_aux[i]
    }
  }
  
  min_AT <- apply(cbind(min_AT,min_aux),1,min)
  
  #----- inferolateral
  
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "IL",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder, 
                         normalising = FALSE, flag_debugging = flag_debugging)
  if(cohort == "RRHF"){
    bipoles_df <- cbind(aux[[1]],aux[[2]])
  }
  else{
    bipoles_df <- aux[[1]]
  }
  
  min_aux <- mapply(min,bipoles_df)
  
  min_names_aux <- rownames(bipoles_df)[mapply(which.min,bipoles_df)]
  
  old_AT_new_AT <- min_AT <= min_aux
  
  for(i in c(1:length(old_AT_new_AT))){
    if(!old_AT_new_AT[i]){
      min_names[i] <- min_names_aux[i]
    }
  }
  
  min_AT <- apply(cbind(min_AT,min_aux),1,min)
  
  #----- inferior
  
  aux <- Read_dataframes(metric_option = metric_option, lead_option = "IN",
                         which_cases = cohort, version = 4,
                         SA_folder = SA_folder, 
                         normalising = FALSE, flag_debugging = flag_debugging)
  if(cohort == "RRHF"){
    bipoles_df <- cbind(aux[[1]],aux[[2]])
  }
  else{
    bipoles_df <- aux[[1]]
  }
  
  min_aux <- mapply(min,bipoles_df)
  
  min_names_aux <- rownames(bipoles_df)[mapply(which.min,bipoles_df)]
  
  old_AT_new_AT <- min_AT <= min_aux
  
  for(i in c(1:length(old_AT_new_AT))){
    if(!old_AT_new_AT[i]){
      min_names[i] <- min_names_aux[i]
    }
  }
  
  min_AT <- apply(cbind(min_AT,min_aux),1,min)
  } 
  
  else{
    aux <- Read_dataframes(metric_option = metric_option, lead_option = vein,
                           which_cases = cohort, version = 4,
                           SA_folder = SA_folder,
                           normalising = FALSE, flag_debugging = flag_debugging)
    if(cohort == "RRHF"){
      bipoles_df <- cbind(aux[[1]],aux[[2]])
    }
    else{
      bipoles_df <- aux[[1]]
    }
    
    min_AT <- mapply(min,bipoles_df)
    
    min_names <- rownames(bipoles_df)[mapply(which.min, bipoles_df)]
    names(min_names) <- names(min_AT)
  }
  
  return(list(min_AT,bin2ascii_lead(min_names)))
}

#' @description Function to print the TAT, AT1090, VEUTAT and VEUmean and
#' improvement taking into account personalisation or a common global optimal
#' lead accross the different cohorts.
#' 
#' @param SA_folder Name of the folder of the sensitivity analysis case.
#' @param flag_debugging If TRUE, prints whatever is reading and writing.
#' 
#' @return All the designs, global and personalised for all the cohorts.

Print_Statistics_AT <- function(SA_folder, vein = "ALL", 
                                flag_debugging = FALSE){
  
  source('/home/crg17/Desktop/scripts/multipole/R/common_functions.R')
  source('/home/crg17/Desktop/scripts/multipole/R/optimality_functions.R')
  
  # We read the RV apex:
    RR_baseline <- Read_csv(file = paste0("/data/SA_multipole/",SA_folder,
                                          "/h/multipole_RVapex.dat"), sep = "",
                            flag_debugging = flag_debugging)
    HF_baseline <- Read_csv(file = paste0("/data/SA_multipole/",SA_folder,
                                          "/HF/multipole_RVapex.dat"), sep = "",
                            flag_debugging = flag_debugging)
  

  rownames(RR_baseline) <- c(paste0("RR",rownames(RR_baseline)))
  rownames(HF_baseline) <- c(paste0("HF",rownames(HF_baseline)))
  RRHF_baseline <- rbind(RR_baseline,HF_baseline)
  
  # Global optimal, TAT, for all the cohorts
  RR_global_TAT <- Print_Global_Optimal_AT(cohort = "RR", metric_option = "TAT",
                                           SA_folder = SA_folder, vein = vein,
                                           flag_debugging = flag_debugging)

  HF_global_TAT <- Print_Global_Optimal_AT(cohort = "HF", metric_option = "TAT",
                                           SA_folder = SA_folder, vein = vein,
                                           flag_debugging = flag_debugging)

  RRHF_global_TAT <- Print_Global_Optimal_AT(cohort = "RRHF",
                                             metric_option = "TAT",
                                             SA_folder = SA_folder, vein = vein,
                                             flag_debugging = flag_debugging)

  # Global optimal, AT1090, for all the cohorts

  RR_global_AT1090 <- Print_Global_Optimal_AT(cohort = "RR",
                                              metric_option = "AT1090",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)

  HF_global_AT1090 <- Print_Global_Optimal_AT(cohort = "HF",
                                              metric_option = "AT1090",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)

  RRHF_global_AT1090 <- Print_Global_Optimal_AT(cohort = "RRHF",
                                                metric_option = "AT1090",
                                                SA_folder = SA_folder,
                                                vein = vein,
                                                flag_debugging = flag_debugging)
  # Global optimal, AT090, for all the cohorts
  
  RR_global_AT090 <- Print_Global_Optimal_AT(cohort = "RR",
                                              metric_option = "AT090",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)
  
  HF_global_AT090 <- Print_Global_Optimal_AT(cohort = "HF",
                                              metric_option = "AT090",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)
  
  RRHF_global_AT090 <- Print_Global_Optimal_AT(cohort = "RRHF",
                                                metric_option = "AT090",
                                                SA_folder = SA_folder,
                                                vein = vein,
                                                flag_debugging = flag_debugging)
  
  # Global optimal, VEUTAT, for all the cohorts
  RR_global_VEUTAT <- Print_Global_Optimal_AT(cohort = "RR",
                                              metric_option = "VEUTAT",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)

  HF_global_VEUTAT <- Print_Global_Optimal_AT(cohort = "HF",
                                              metric_option = "VEUTAT",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)

  RRHF_global_VEUTAT <- Print_Global_Optimal_AT(cohort = "RRHF",
                                                metric_option = "VEUTAT",
                                                SA_folder = SA_folder,
                                                vein = vein,
                                                flag_debugging = flag_debugging)

  # Global optimal, VEUmean, for all the cohorts

  RR_global_VEUmean <- Print_Global_Optimal_AT(cohort = "RR",
                                              metric_option = "VEUmean",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)

  HF_global_VEUmean <- Print_Global_Optimal_AT(cohort = "HF",
                                              metric_option = "VEUmean",
                                              SA_folder = SA_folder,
                                              vein = vein,
                                              flag_debugging = flag_debugging)

  RRHF_global_VEUmean <- Print_Global_Optimal_AT(cohort = "RRHF",
                                                metric_option = "VEUmean",
                                                SA_folder = SA_folder,
                                                vein = vein,
                                                flag_debugging = flag_debugging)
  
  # Global optimal, LVTAT, for all the cohorts
  RR_global_LVTAT <- Print_Global_Optimal_AT(cohort = "RR", 
                                           metric_option = "LVTAT",
                                           SA_folder = SA_folder, vein = vein,
                                           flag_debugging = flag_debugging)
  
  HF_global_LVTAT <- Print_Global_Optimal_AT(cohort = "HF", 
                                           metric_option = "LVTAT",
                                           SA_folder = SA_folder, vein = vein,
                                           flag_debugging = flag_debugging)
  
  RRHF_global_LVTAT <- Print_Global_Optimal_AT(cohort = "RRHF",
                                             metric_option = "LVTAT",
                                             SA_folder = SA_folder, vein = vein,
                                             flag_debugging = flag_debugging)
  
  
  # Personalised optimal, TAT, for all the cohorts
  
  RR_personalised_TAT <- 
    Print_Individual_Optimal_AT(cohort = "RR",metric_option = "TAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  HF_personalised_TAT <- 
    Print_Individual_Optimal_AT(cohort = "HF",metric_option = "TAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  RRHF_personalised_TAT <- 
    Print_Individual_Optimal_AT(cohort = "RRHF",metric_option = "TAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  # Personalised optimal, AT1090, for all the cohorts
  
  RR_personalised_AT1090 <- 
    Print_Individual_Optimal_AT(cohort = "RR",metric_option = "AT1090",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  HF_personalised_AT1090 <-
    Print_Individual_Optimal_AT(cohort = "HF",metric_option = "AT1090",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  RRHF_personalised_AT1090 <-
    Print_Individual_Optimal_AT(cohort = "RRHF",metric_option = "AT1090",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  # Personalised optimal, AT090, for all the cohorts
  
  RR_personalised_AT090 <- 
    Print_Individual_Optimal_AT(cohort = "RR",metric_option = "AT090",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  HF_personalised_AT090 <-
    Print_Individual_Optimal_AT(cohort = "HF",metric_option = "AT090",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  RRHF_personalised_AT090 <-
    Print_Individual_Optimal_AT(cohort = "RRHF",metric_option = "AT090",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  # Personalised optimal, VEUTAT, for all the cohorts
  
  RR_personalised_VEUTAT <- 
    Print_Individual_Optimal_AT(cohort = "RR",metric_option = "VEUTAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  HF_personalised_VEUTAT <-
    Print_Individual_Optimal_AT(cohort = "HF",metric_option = "VEUTAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  RRHF_personalised_VEUTAT <-
    Print_Individual_Optimal_AT(cohort = "RRHF",metric_option = "VEUTAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  # Personalised optimal, VEUmean, for all the cohorts
  
  RR_personalised_VEUmean <- 
    Print_Individual_Optimal_AT(cohort = "RR",metric_option = "VEUmean",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  HF_personalised_VEUmean <-
    Print_Individual_Optimal_AT(cohort = "HF",metric_option = "VEUmean",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  RRHF_personalised_VEUmean <-
    Print_Individual_Optimal_AT(cohort = "RRHF",metric_option = "VEUmean",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  # Personalised optimal, LVTAT, for all the cohorts
  
  RR_personalised_LVTAT <- 
    Print_Individual_Optimal_AT(cohort = "RR",metric_option = "LVTAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  HF_personalised_LVTAT <-
    Print_Individual_Optimal_AT(cohort = "HF",metric_option = "LVTAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  RRHF_personalised_LVTAT <-
    Print_Individual_Optimal_AT(cohort = "RRHF",metric_option = "LVTAT",
                                SA_folder = SA_folder,
                                flag_debugging = flag_debugging)
  
  #---------------------- Prints
  #
  #
  #
  #--------------- Baseline AT
  #
  Print_baseline <- function(baseline,baseline_name){
    
    cohort_name <- sub("_.*","",baseline_name)
    
    sprintf("Baseline %s TAT: %i \U00B1 %i ms",
            cohort_name,
            round(mean(baseline$TAT)),
            round(sd(baseline$TAT))) %>%
      print(.)
    
    sprintf("Baseline %s AT1090: %i \U00B1 %i ms",
            cohort_name,
            round(mean(baseline$AT1090)),
            round(sd(baseline$AT1090))) %>%
      print(.)
    
    sprintf("Baseline %s AT090: %i \U00B1 %i ms",
            cohort_name,
            round(mean(baseline$AT090)),
            round(sd(baseline$AT090))) %>%
      print(.)
    
    sprintf("Baseline %s VEUTAT: %i \U00B1 %i ms",
            cohort_name,
            round(mean(baseline$VEUTAT)),
            round(sd(baseline$VEUTAT))) %>%
      print(.)
    
    sprintf("Baseline %s VEUmean: %i \U00B1 %i ms",
            cohort_name,
            round(mean(baseline$VEUmean)),
            round(sd(baseline$VEUmean))) %>%
      print(.)
    
    sprintf("Baseline %s LVTAT: %i \U00B1 %i ms",
            cohort_name,
            round(mean(baseline$LVTAT)),
            round(sd(baseline$LVTAT))) %>%
      print(.)
  }
  
  Print_baseline(RR_baseline,"RR_baseline")
  Print_baseline(HF_baseline,"HF_baseline")
  Print_baseline(RRHF_baseline,"RRHF_baseline")
  
  
  #------------- Global optimal AT
  
  Print_global <- function(global,var_name){
    
    cohort_name <- sub("_.*","",var_name)
    metric_name <- sub(".*_","",var_name)
    
    sprintf("Optimal global %s %s: %i \U00B1 %i ms",
            cohort_name,
            metric_name,
            round(mean(global[,1])),
            round(sd(global[,1]))) %>%
      print(.)
    
  }
  
  Print_global(RR_global_TAT,"RR_global_TAT")
  Print_global(RR_global_AT1090,"RR_global_AT1090")
  Print_global(RR_global_AT090,"RR_global_AT090")
  Print_global(RR_global_VEUTAT,"RR_global_VEUTAT")
  Print_global(RR_global_VEUmean,"RR_global_VEUmean")
  Print_global(RR_global_LVTAT,"RR_global_LVTAT")
  Print_global(HF_global_TAT,"HF_global_TAT")
  Print_global(HF_global_AT1090,"HF_global_AT1090")
  Print_global(HF_global_AT090,"HF_global_AT090")
  Print_global(HF_global_VEUTAT,"HF_global_VEUTAT")
  Print_global(HF_global_VEUmean,"HF_global_VEUmean")
  Print_global(HF_global_LVTAT,"HF_global_LVTAT")
  Print_global(RRHF_global_TAT,"RRHF_global_TAT")
  Print_global(RRHF_global_AT1090,"RRHF_global_AT1090")
  Print_global(RRHF_global_AT090,"RRHF_global_AT090")
  Print_global(RRHF_global_VEUTAT,"RRHF_global_VEUTAT")
  Print_global(RRHF_global_VEUmean,"RRHF_global_VEUmean")
  Print_global(RRHF_global_LVTAT,"RRHF_global_LVTAT")
  
  #-------------- Global Delta
  #
  Print_global_delta <- function(baseline,global,var_name){
    cohort_name <- sub("_.*","",var_name)
    metric_name <- sub(".*_","",var_name)
    
    v1 <- (baseline[,metric_name]-global[,1]) %>% mean(.) %>% round(.)
    v2 <- (baseline[,metric_name]-global[,1]) %>% sd(.) %>% round(.)
    v3 <- (100*(1-global[,1]/baseline[,metric_name])) %>% mean(.) %>% round(.) 
    v4 <- (100*(1-global[,1]/baseline[,metric_name])) %>% sd(.) %>% round(.)  
    
    sprintf("Global %s d%s: %i \U00B1 %i ms (%i%% \U00B1 %i%%)",
            cohort_name,metric_name,v1,v2,v3,v4) %>%
      print(.)
  }
  
  Print_global_delta(RR_baseline,RR_global_TAT,"RR_global_TAT")
  Print_global_delta(RR_baseline,RR_global_AT1090,"RR_global_AT1090")
  Print_global_delta(RR_baseline,RR_global_AT090,"RR_global_AT090")
  Print_global_delta(RR_baseline,RR_global_VEUTAT,"RR_global_VEUTAT")
  Print_global_delta(RR_baseline,RR_global_VEUmean,"RR_global_VEUmean")
  Print_global_delta(RR_baseline,RR_global_LVTAT,"RR_global_LVTAT")
  Print_global_delta(HF_baseline,HF_global_TAT,"HF_global_TAT")
  Print_global_delta(HF_baseline,HF_global_AT1090,"HF_global_AT1090")
  Print_global_delta(HF_baseline,HF_global_AT090,"HF_global_AT090")
  Print_global_delta(HF_baseline,HF_global_VEUTAT,"HF_global_VEUTAT")
  Print_global_delta(HF_baseline,HF_global_VEUmean,"HF_global_VEUmean")
  Print_global_delta(HF_baseline,HF_global_LVTAT,"HF_global_LVTAT")
  Print_global_delta(RRHF_baseline,RRHF_global_TAT,"RRHF_global_TAT")
  Print_global_delta(RRHF_baseline,RRHF_global_AT1090,"RRHF_global_AT1090")
  Print_global_delta(RRHF_baseline,RRHF_global_AT090,"RRHF_global_AT090")
  Print_global_delta(RRHF_baseline,RRHF_global_VEUTAT,"RRHF_global_VEUTAT")
  Print_global_delta(RRHF_baseline,RRHF_global_VEUmean,"RRHF_global_VEUmean")
  Print_global_delta(RRHF_baseline,RRHF_global_LVTAT,"RRHF_global_LVTAT")
  
  #--------- Optimal leads
  
  paste0("Optimal lead design for RR, according to TAT: ",
         colnames(RR_global_TAT)) %>% print(.)
  paste0("Optimal lead design for RR, according to AT1090: ",
         colnames(RR_global_AT1090)) %>% print(.)
  paste0("Optimal lead design for RR, according to AT090: ",
         colnames(RR_global_AT090)) %>% print(.)
  paste0("Optimal lead design for RR, according to VEUTAT: ",
         colnames(RR_global_VEUTAT)) %>% print(.)
  paste0("Optimal lead design for RR, according to VEUmean: ",
         colnames(RR_global_VEUmean)) %>% print(.)
  paste0("Optimal lead design for RR, according to LVTAT: ",
         colnames(RR_global_LVTAT)) %>% print(.)
  
  paste0("Optimal lead design for HF, according to TAT: ",
         colnames(HF_global_TAT)) %>% print(.)
  paste0("Optimal lead design for HF, according to AT1090: ",
         colnames(HF_global_AT1090)) %>% print(.)
  paste0("Optimal lead design for HF, according to AT090: ",
         colnames(HF_global_AT090)) %>% print(.)
  paste0("Optimal lead design for HF, according to VEUTAT: ",
         colnames(HF_global_VEUTAT)) %>% print(.)
  paste0("Optimal lead design for HF, according to VEUmean: ",
         colnames(HF_global_VEUmean)) %>% print(.)
  paste0("Optimal lead design for HF, according to LVTAT: ",
         colnames(HF_global_LVTAT)) %>% print(.)
  
  paste0("Optimal lead design for RRHF, according to TAT: ",
         colnames(RRHF_global_TAT)) %>% print(.)
  paste0("Optimal lead design for RRHF, according to AT1090: ",
         colnames(RRHF_global_AT1090)) %>% print(.)
  paste0("Optimal lead design for RRHF, according to AT090: ",
         colnames(RRHF_global_AT090)) %>% print(.)
  paste0("Optimal lead design for RRHF, according to VEUTAT: ",
         colnames(RRHF_global_VEUTAT)) %>% print(.)
  paste0("Optimal lead design for RRHF, according to VEUmean: ",
         colnames(RRHF_global_VEUmean)) %>% print(.)
  paste0("Optimal lead design for RRHF, according to LVTAT: ",
         colnames(RRHF_global_LVTAT)) %>% print(.)  
  #---------------------- Personalised optimal AT
  #
  #---------- RR
  
  sprintf("Personalised RR TAT: %i \U00B1 %i ms",
          round(mean(RR_personalised_TAT[[1]])),
          round(sd(RR_personalised_TAT[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RR AT1090: %i \U00B1 %i ms",
          round(mean(RR_personalised_AT1090[[1]])),
          round(sd(RR_personalised_AT1090[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RR AT090: %i \U00B1 %i ms",
          round(mean(RR_personalised_AT090[[1]])),
          round(sd(RR_personalised_AT090[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RR VEUTAT: %i \U00B1 %i ms",
          round(mean(RR_personalised_VEUTAT[[1]])),
          round(sd(RR_personalised_VEUTAT[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RR VEUmean: %i \U00B1 %i ms",
          round(mean(RR_personalised_VEUmean[[1]])),
          round(sd(RR_personalised_VEUmean[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RR LVTAT: %i \U00B1 %i ms",
          round(mean(RR_personalised_LVTAT[[1]])),
          round(sd(RR_personalised_LVTAT[[1]]))) %>%
    print(.)
  
  #---------- HF
  
  sprintf("Personalised HF TAT: %i \U00B1 %i ms",
          round(mean(HF_personalised_TAT[[1]])),
          round(sd(HF_personalised_TAT[[1]]))) %>%
    print(.)
  
  sprintf("Personalised HF AT1090: %i \U00B1 %i ms",
          round(mean(HF_personalised_AT1090[[1]])),
          round(sd(HF_personalised_AT1090[[1]]))) %>%
    print(.)  
  
  sprintf("Personalised HF AT090: %i \U00B1 %i ms",
          round(mean(HF_personalised_AT090[[1]])),
          round(sd(HF_personalised_AT090[[1]]))) %>%
    print(.)
  
  sprintf("Personalised HF VEUTAT: %i \U00B1 %i ms",
          round(mean(HF_personalised_VEUTAT[[1]])),
          round(sd(HF_personalised_VEUTAT[[1]]))) %>%
    print(.)
  
  sprintf("Personalised HF VEUmean: %i \U00B1 %i ms",
          round(mean(HF_personalised_VEUmean[[1]])),
          round(sd(HF_personalised_VEUmean[[1]]))) %>%
    print(.)
  
  sprintf("Personalised HF LVTAT: %i \U00B1 %i ms",
          round(mean(HF_personalised_LVTAT[[1]])),
          round(sd(HF_personalised_LVTAT[[1]]))) %>%
    print(.)
  
  #---------- RRHF
  
  sprintf("Personalised RRHF TAT: %i \U00B1 %i ms",
          round(mean(RRHF_personalised_TAT[[1]])),
          round(sd(RRHF_personalised_TAT[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RRHF AT1090: %i \U00B1 %i ms",
          round(mean(RRHF_personalised_AT1090[[1]])),
          round(sd(RRHF_personalised_AT1090[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RRHF AT090: %i \U00B1 %i ms",
          round(mean(RRHF_personalised_AT090[[1]])),
          round(sd(RRHF_personalised_AT090[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RRHF VEUTAT: %i \U00B1 %i ms",
          round(mean(RRHF_personalised_VEUTAT[[1]])),
          round(sd(RRHF_personalised_VEUTAT[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RRHF VEUmean: %i \U00B1 %i ms",
          round(mean(RRHF_personalised_VEUmean[[1]])),
          round(sd(RRHF_personalised_VEUmean[[1]]))) %>%
    print(.)
  
  sprintf("Personalised RRHF LVTAT: %i \U00B1 %i ms",
          round(mean(RRHF_personalised_LVTAT[[1]])),
          round(sd(RRHF_personalised_LVTAT[[1]]))) %>%
    print(.)
  
  #-------------- Personalised Delta
  #
  Print_personalised_delta <- function(baseline,personalised,var_name){
    cohort_name <- sub("_.*","",var_name)
    metric_name <- sub(".*_","",var_name)
    
    
    v1 <- (baseline[,metric_name]-personalised[[1]]) %>%
      mean(.) %>% round(.)
    v2 <- (baseline[,metric_name]-personalised[[1]]) %>%
      sd(.) %>% round(.)
    v3 <- (100*(1-personalised[[1]]/baseline[,metric_name])) %>%
      mean(.) %>% round(.) 
    v4 <- (100*(1-personalised[[1]]/baseline[,metric_name])) %>%
      sd(.) %>% round(.)  
    
    sprintf("Personalised %s d%s: %i \U00B1 %i ms (%i%% \U00B1 %i%%)",
            cohort_name,metric_name,v1,v2,v3,v4) %>%
      print(.)
  }
  
  Print_personalised_delta(RR_baseline,RR_personalised_TAT,
                           "RR_personalised_TAT")
  Print_personalised_delta(RR_baseline,RR_personalised_AT1090,
                           "RR_personalised_AT1090")
  Print_personalised_delta(RR_baseline,RR_personalised_AT090,
                           "RR_personalised_AT090")
  Print_personalised_delta(RR_baseline,RR_personalised_VEUTAT,
                           "RR_personalised_VEUTAT")
  Print_personalised_delta(RR_baseline,RR_personalised_VEUmean,
                           "RR_personalised_VEUmean")
  Print_personalised_delta(RR_baseline,RR_personalised_LVTAT,
                           "RR_personalised_LVTAT")
  Print_personalised_delta(HF_baseline,HF_personalised_TAT,
                           "HF_personalised_TAT")
  Print_personalised_delta(HF_baseline,HF_personalised_AT1090,
                           "HF_personalised_AT1090")
  Print_personalised_delta(HF_baseline,HF_personalised_AT090,
                           "HF_personalised_AT090")
  Print_personalised_delta(HF_baseline,HF_personalised_VEUTAT,
                           "HF_personalised_VEUTAT")
  Print_personalised_delta(HF_baseline,HF_personalised_VEUmean,
                           "HF_personalised_VEUmean")
  Print_personalised_delta(HF_baseline,HF_personalised_LVTAT,
                           "HF_personalised_LVTAT")
  Print_personalised_delta(RRHF_baseline,RRHF_personalised_TAT,
                           "RRHF_personalised_TAT")
  Print_personalised_delta(RRHF_baseline,RRHF_personalised_AT1090,
                           "RRHF_personalised_AT1090")
  Print_personalised_delta(RRHF_baseline,RRHF_personalised_AT090,
                           "RRHF_personalised_AT090")
  Print_personalised_delta(RRHF_baseline,RRHF_personalised_VEUTAT,
                           "RRHF_personalised_VEUTAT")
  Print_personalised_delta(RRHF_baseline,RRHF_personalised_VEUmean,
                           "RRHF_personalised_VEUmean")
  Print_personalised_delta(RRHF_baseline,RRHF_personalised_LVTAT,
                           "RRHF_personalised_LVTAT")
  
  #--------- Optimal leads
  #
  Print_personalised_optimal_leads <- function(global,personalised,var_name){
    person_bip_ascii <- personalised[[2]]
    global_quad_ascii <- colnames(global)
    
    not_in_quapole<- sapply(person_bip_ascii,
                            function(x) !isbipoleinquadpole(x,global_quad_ascii)) %>%
      as.vector(.)
    
    cohort_name <- sub("_.*","",var_name)
    metric_name <- sub(".*_","",var_name)
    
    if(sum(not_in_quapole) == 0){
      sprintf(paste0("All the optimal bipoles for the %s cohort according to %s",
                     " are already included in the global optimal quadripole"),
              cohort_name,metric_name) %>%
      print(.)
    }
    else{
      sprintf(paste0("The electrodes not included in the global optimal lead", 
                     " for the %s cohort according to %s is/are %s"),
              cohort_name,metric_name,
              paste(unique(person_bip_ascii[not_in_quapole]), collapse = ", ")) %>%
        print(.)
    }
  }
  
  Print_personalised_optimal_leads(RR_global_TAT,RR_personalised_TAT,
                                   "RR_personalised_TAT")
  Print_personalised_optimal_leads(RR_global_AT1090,RR_personalised_AT1090,
                                   "RR_personalised_AT1090")
  Print_personalised_optimal_leads(RR_global_AT090,RR_personalised_AT090,
                                   "RR_personalised_AT090")
  Print_personalised_optimal_leads(RR_global_VEUTAT,RR_personalised_VEUTAT,
                                   "RR_personalised_VEUTAT")
  Print_personalised_optimal_leads(RR_global_VEUmean,RR_personalised_VEUmean,
                                   "RR_personalised_VEUmean")
  Print_personalised_optimal_leads(RR_global_LVTAT,RR_personalised_LVTAT,
                                   "RR_personalised_LVTAT")
  Print_personalised_optimal_leads(HF_global_TAT,HF_personalised_TAT,
                                   "HF_personalised_TAT")
  Print_personalised_optimal_leads(HF_global_AT1090,HF_personalised_AT1090,
                                   "HF_personalised_AT1090")
  Print_personalised_optimal_leads(HF_global_AT090,HF_personalised_AT090,
                                   "HF_personalised_AT090")
  Print_personalised_optimal_leads(HF_global_VEUTAT,HF_personalised_VEUTAT,
                                   "HF_personalised_VEUTAT")
  Print_personalised_optimal_leads(HF_global_VEUmean,HF_personalised_VEUmean,
                                   "HF_personalised_VEUmean")
  Print_personalised_optimal_leads(HF_global_LVTAT,HF_personalised_LVTAT,
                                   "HF_personalised_LVTAT")
  Print_personalised_optimal_leads(RRHF_global_TAT,RRHF_personalised_TAT,
                                   "RRHF_personalised_TAT")
  Print_personalised_optimal_leads(RRHF_global_AT1090,RRHF_personalised_AT1090,
                                   "RRHF_personalised_AT1090")
  Print_personalised_optimal_leads(RRHF_global_AT090,RRHF_personalised_AT090,
                                   "RRHF_personalised_AT090")
  Print_personalised_optimal_leads(RRHF_global_VEUTAT,RRHF_personalised_VEUTAT,
                                   "RRHF_personalised_VEUTAT")
  Print_personalised_optimal_leads(RRHF_global_VEUmean,RRHF_personalised_VEUmean,
                                   "RRHF_personalised_VEUmean")
  Print_personalised_optimal_leads(RRHF_global_LVTAT,RRHF_personalised_LVTAT,
                                   "RRHF_personalised_LVTAT")
  

  
  return(list(RR_baseline = RR_baseline,HF_baseline = HF_baseline,
              RRHF_baseline = RRHF_baseline, RR_global_TAT = RR_global_TAT,
              RR_global_AT1090 = RR_global_AT1090,
              RR_global_AT090 = RR_global_AT090,
              RR_global_VEUTAT = RR_global_VEUTAT,
              RR_global_VEUmean = RR_global_VEUmean,
              RR_global_LVTAT = RR_global_LVTAT,
              HF_global_TAT = HF_global_TAT,
              HF_global_AT1090 = HF_global_AT1090,
              HF_global_AT090 = HF_global_AT090, 
              HF_global_VEUTAT = HF_global_VEUTAT,
              HF_global_VEUmean = HF_global_VEUmean,
              HF_global_LVTAT = HF_global_LVTAT,
              RRHF_global_TAT = RRHF_global_TAT,
              RRHF_global_AT1090 = RRHF_global_AT1090,
              RRHF_global_AT090 = RRHF_global_AT090,
              RRHF_global_VEUTAT = RRHF_global_VEUTAT,
              RRHF_global_VEUmean = RRHF_global_VEUmean,
              RRHF_global_LVTAT = RRHF_global_LVTAT,
              RR_personalised_TAT = RR_personalised_TAT,
              RR_personalised_AT1090 = RR_personalised_AT1090,
              RR_personalised_AT090 = RR_personalised_AT090,
              RR_personalised_VEUTAT = RR_personalised_VEUTAT,
              RR_personalised_VEUmean = RR_personalised_VEUmean,
              RR_personalised_LVTAT = RR_personalised_LVTAT,
              HF_personalised_TAT = HF_personalised_TAT,
              HF_personalised_AT1090 = HF_personalised_AT1090,
              HF_personalised_AT090 = HF_personalised_AT090,
              HF_personalised_VEUTAT = HF_personalised_VEUTAT,
              HF_personalised_VEUmean = HF_personalised_VEUmean,
              HF_personalised_LVTAT = HF_personalised_LVTAT,
              RRHF_personalised_TAT = RRHF_personalised_TAT,
              RRHF_personalised_AT1090 = RRHF_personalised_AT1090,
              RRHF_personalised_AT090 = RRHF_personalised_AT090,
              RRHF_personalised_VEUTAT = RRHF_personalised_VEUTAT,
              RRHF_personalised_VEUmean = RRHF_personalised_VEUmean,
              RRHF_personalised_LVTAT = RRHF_personalised_LVTAT
              ))
}


#' @description Function to print several results concerning the optimal 
#' configurations of MPP leads.
#' 
#' @param SA_folder Name of the folder of the sensitivity analysis case.
#' @param vein Name of the lead to check the optimal configuration. If "ALL", 
#' it loops over all the possible cases.
#' @param metric_option Name of the metric we optimise against.
#' @param num_decimals Number of decimal numbers to print.
#' @param return_flag If TRUE, it returns baselines, AT durations with
#' the global optimal and AT durations with the personalised optimal 
#' configurations for all the cohorts.
#' @param output... If TRUE, it prints that particular value.  
#' @param flag_debugging If TRUE, prints whatever is reading and writing.
#' 
#' @return Check return_flag parameter.

Print_results_metric <- function(SA_folder, vein = "ALL", metric_option,
                                 num_decimals = 0, return_flag = FALSE,
                                 output_baseline = TRUE, output_global = TRUE,
                                 output_dglobal = TRUE, output_lead = TRUE,
                                 output_pers = TRUE, output_dpers = TRUE,
                                 output_notinlead = TRUE,
                                 flag_debugging = FALSE){
  source('/home/crg17/Desktop/scripts/multipole/R/common_functions.R')
  source('/home/crg17/Desktop/scripts/multipole/R/optimality_functions.R')
  
  if(vein == "ALL"){
    for(vein_name in c("AN","AL","LA","IL","IN")){
    print(paste0("Results for ",vein_name," position:"))
    Print_results_metric(SA_folder = SA_folder, vein = vein_name,
                         metric_option = metric_option,
                         num_decimals = num_decimals, return_flag = return_flag,
                         output_baseline = output_baseline,
                         output_global = output_global,
                         output_dglobal = output_dglobal,
                         output_lead = output_lead, output_pers = output_pers,
                         output_dpers = output_dpers,
                         output_notinlead = output_notinlead,
                         flag_debugging = flag_debugging)
    }
  }
  
    else{
    # We read the RV apex:
      RR_baseline <- Read_csv(file = paste0("/data/SA_multipole/",SA_folder,
                                            "/h/multipole_RVapex.dat"), sep = "",
                              flag_debugging = flag_debugging)
      HF_baseline <- Read_csv(file = paste0("/data/SA_multipole/",SA_folder,
                                            "/HF/multipole_RVapex.dat"), sep = "",
                              flag_debugging = flag_debugging)
    

    rownames(RR_baseline) <- c(paste0("RR",rownames(RR_baseline)))
    rownames(HF_baseline) <- c(paste0("HF",rownames(HF_baseline)))
    
    RRHF_baseline <- rbind(RR_baseline,HF_baseline)
    
    # Global optimal for all the cohorts
    RR_global <- Print_Global_Optimal_AT(cohort = "RR",
                                         metric_option = metric_option,
                                         SA_folder = SA_folder, vein = vein,
                                         flag_debugging = flag_debugging)
    
    HF_global <- Print_Global_Optimal_AT(cohort = "HF",
                                         metric_option = metric_option,
                                         SA_folder = SA_folder, vein = vein,
                                         flag_debugging = flag_debugging)
    
    RRHF_global <- Print_Global_Optimal_AT(cohort = "RRHF",
                                           metric_option = metric_option,
                                           SA_folder = SA_folder, vein = vein,
                                           flag_debugging = flag_debugging)
    
    # Personalised optimal, for all the cohorts
    
    RR_personalised <- 
      Print_Individual_Optimal_AT(cohort = "RR", metric_option = metric_option,
                                  SA_folder = SA_folder, vein = vein,
                                  flag_debugging = flag_debugging)
    
    HF_personalised <- 
      Print_Individual_Optimal_AT(cohort = "HF", metric_option = metric_option,
                                  SA_folder = SA_folder, vein = vein,
                                  flag_debugging = flag_debugging)
    
    RRHF_personalised <- 
      Print_Individual_Optimal_AT(cohort = "RRHF", metric_option = metric_option,
                                  SA_folder = SA_folder, vein = vein,
                                  flag_debugging = flag_debugging)
    
    
    #---------------------- Prints
    #
    #
    #
    #--------------- Baseline AT
    #
    #
    
    if(output_baseline){
    
    Print_baseline <- function(baseline, cohort, metric_option, num_decimals){
      
  
      sprintf("Baseline %s %s: %g \U00B1 %g ms",
              cohort,
              metric_option,
              round(mean(baseline[, metric_option]), num_decimals),
              round(sd(baseline[, metric_option])), num_decimals) %>%
        print(.)
      
    }
    
    Print_baseline(baseline = RR_baseline, cohort = "RR",
                   metric_option = metric_option, num_decimals = num_decimals)
    Print_baseline(baseline = HF_baseline, cohort = "HF",
                   metric_option = metric_option, num_decimals = num_decimals)
    Print_baseline(baseline = RRHF_baseline, cohort = "RRHF",
                   metric_option = metric_option, num_decimals = num_decimals)
    }
    
    #------------- Global optimal AT
    
    if(output_global){
    Print_global <- function(global, cohort, metric_option, num_decimals){
      
      
      sprintf("Optimal global %s %s: %g \U00B1 %g ms",
              cohort,
              metric_option,
              round(mean(global[,1]), num_decimals),
              round(sd(global[,1])), num_decimals) %>%
        print(.)
      
    }
    
    Print_global(global = RR_global, cohort = "RR",
                 metric_option = metric_option, num_decimals = num_decimals)
    Print_global(global = HF_global, cohort = "HF",
                 metric_option = metric_option, num_decimals = num_decimals)
    Print_global(global = RRHF_global, cohort = "RRHF",
                 metric_option = metric_option, num_decimals = num_decimals)
    }
    #-------------- Global Delta
    #
    if(output_dglobal){
    Print_global_delta <- function(baseline, global, cohort, metric_option,
                                   num_decimals){
      
      v1 <- (baseline[,metric_option]-global[,1]) %>% mean(.) %>%
        round(., num_decimals)
      v2 <- (baseline[,metric_option]-global[,1]) %>% sd(.) %>%
        round(., num_decimals)
      v3 <- (100*(1-global[,1]/baseline[,metric_option])) %>% mean(.) %>%
        round(., num_decimals) 
      v4 <- (100*(1-global[,1]/baseline[,metric_option])) %>% sd(.) %>%
        round(., num_decimals)  
      
      sprintf("Global %s d%s: %g \U00B1 %g ms (%g%% \U00B1 %g%%)",
              cohort,metric_option,v1,v2,v3,v4) %>%
        print(.)
    }
    
    Print_global_delta(baseline = RR_baseline, global = RR_global,
                       cohort = "RR", metric_option = metric_option,
                       num_decimals = num_decimals)
    Print_global_delta(baseline = HF_baseline, global = HF_global,
                       cohort = "HF", metric_option = metric_option,
                       num_decimals = num_decimals)
    Print_global_delta(baseline = RRHF_baseline, global = RRHF_global,
                       cohort = "RRHF", metric_option = metric_option,
                       num_decimals = num_decimals)
    }
    #--------- Optimal leads
    if(output_lead){
    paste0("Optimal lead design for RR, according to ",metric_option,": ",
           colnames(RR_global)) %>% print(.)  
    paste0("Optimal lead design for HF, according to ",metric_option,": ",
           colnames(HF_global)) %>% print(.)
    paste0("Optimal lead design for RRHF, according to ",metric_option,": ",
           colnames(RRHF_global)) %>% print(.)
    }
    #---------------------- Personalised optimal AT
    #
  if(output_pers){
    sprintf(paste0("Personalised RR ",metric_option,": %g \U00B1 %g ms"),
            round(mean(RR_personalised[[1]]), num_decimals),
            round(sd(RR_personalised[[1]])), num_decimals) %>%
      print(.)
    
    sprintf(paste0("Personalised HF ",metric_option,": %g \U00B1 %g ms"),
            round(mean(HF_personalised[[1]]), num_decimals),
            round(sd(HF_personalised[[1]])), num_decimals) %>%
      print(.)
    
    sprintf(paste0("Personalised RRHF ",metric_option,": %g \U00B1 %g ms"),
            round(mean(RRHF_personalised[[1]]), num_decimals),
            round(sd(RRHF_personalised[[1]])), num_decimals) %>%
      print(.)
  }
    #-------------- Personalised Delta
    #
    if(output_dpers){
    Print_personalised_delta <- function(baseline, personalised, cohort,
                                         metric_option, num_decimals){
      
      
      v1 <- (baseline[,metric_option]-personalised[[1]]) %>%
        mean(.) %>% round(., num_decimals)
      v2 <- (baseline[,metric_option]-personalised[[1]]) %>%
        sd(.) %>% round(., num_decimals)
      v3 <- (100*(1-personalised[[1]]/baseline[,metric_option])) %>%
        mean(.) %>% round(., num_decimals) 
      v4 <- (100*(1-personalised[[1]]/baseline[,metric_option])) %>%
        sd(.) %>% round(., num_decimals)  
      
      sprintf("Personalised %s d%s: %g \U00B1 %g ms (%g%% \U00B1 %g%%)",
              cohort,metric_option,v1,v2,v3,v4) %>%
        print(.)
    }
    
    Print_personalised_delta(baseline = RR_baseline,
                             personalised = RR_personalised, cohort = "RR",
                             metric_option = metric_option,
                             num_decimals = num_decimals) 
    Print_personalised_delta(baseline = HF_baseline,
                             personalised = HF_personalised, cohort = "HF",
                             metric_option = metric_option,
                             num_decimals = num_decimals)
    Print_personalised_delta(baseline = RRHF_baseline,
                             personalised = RRHF_personalised, cohort = "RRHF",
                             metric_option = metric_option,
                             num_decimals = num_decimals)
    }
    #--------- Optimal leads
    #
    if(output_notinlead){
    Print_personalised_optimal_leads <- function(global, personalised, cohort,
                                                 metric_option){
      person_bip_ascii <- personalised[[2]]
      global_quad_ascii <- colnames(global)
      
      not_in_quapole<- sapply(person_bip_ascii,
                              function(x) !isbipoleinquadpole(x,global_quad_ascii)) %>%
        as.vector(.)
      
  
      if(sum(not_in_quapole) == 0){
        sprintf(paste0("All the optimal bipoles for the %s cohort according to %s",
                       " are already included in the global optimal quadripole"),
                cohort,metric_option) %>%
          print(.)
      }
      else{
        sprintf(paste0("The electrodes not included in the global optimal lead", 
                       " for the %s cohort according to %s is/are %s"),
                cohort,metric_option,
                paste(unique(person_bip_ascii[not_in_quapole]), collapse = ", ")) %>%
          print(.)
      }
    }
    
    Print_personalised_optimal_leads(global = RR_global,
                                     personalised = RR_personalised,
                                     cohort = "RR",
                                     metric_option = metric_option)
    Print_personalised_optimal_leads(global = HF_global,
                                     personalised = HF_personalised,
                                     cohort = "HF",
                                     metric_option = metric_option)
    Print_personalised_optimal_leads(global = RRHF_global,
                                     personalised = RRHF_personalised,
                                     cohort = "RRHF",
                                     metric_option = metric_option)
    
    }
    if(return_flag){
      return(list(RR_baseline = RR_baseline, HF_baseline = HF_baseline,
                  RRHF_baseline = RRHF_baseline, RR_global = RR_global,
                  HF_global = HF_global,
                  RRHF_global = RRHF_global,
                  RR_personalised = RR_personalised,
                  HF_personalised = HF_personalised,
                  RRHF_personalised = RRHF_personalised
                  ))
    }
    }
  
}

Last_point_activated <- function(heart,which_cases){
  if(which_cases == 'RR')
    which_cases <- 'h'
  
  elem_path <- paste0("/media/crg17/Seagate\\ Backup\\ Plus\\ Drive/CT_cases/",which_cases,"_case",heart)
  
  elem_file_old <- paste0(elem_path,"/meshing/1000um/BiV/FEC/BiV_FEC_w5_h70.elem")
  elem_file_new <- paste0(elem_path,"/meshing/1000um/BiV/FEC/BiV_noheader.elem")
  
  cmd <- paste0("tail -n +2 ",elem_file_old," > ",elem_file_new)
  system(cmd)
  
  elem_path <- paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/",which_cases,"_case",heart)
  elem_file_new <- paste0(elem_path,"/meshing/1000um/BiV/FEC/BiV_noheader.elem")
  vm_act_seq_name <- paste0(elem_path,"/simulations/multipole/eikonal_default/BiV.RV_endo.apex/vm_act_seq.dat")
  
  BiV_elem <- read.table(elem_file_new)
  vm_act_seq <- read.table(vm_act_seq_name)
  
  last_vertex <- which(vm_act_seq==max(vm_act_seq))-1
  
  p1 <- which(BiV_elem$V2 == last_vertex)
  p2 <- which(BiV_elem$V3 == last_vertex)
  p3 <- which(BiV_elem$V4 == last_vertex)
  p4 <- which(BiV_elem$V5 == last_vertex)
  
  if(length(p1) > 0){
    last_tag <- BiV_elem$V6[p1]
  }
  else if(length(p2) > 0){
    last_tag <- BiV_elem$V6[p2]
  }
  else if(length(p3) > 0){
    last_tag <- BiV_elem$V6[p3]
  }
  else if(length(p4) > 0){
    last_tag <- BiV_elem$V6[p4]
  }
  
  print(sprintf(paste0("Case %s: Last vertex is %d for which last tag is %i"),heart,last_vertex,last_tag))
  
}

#' @description Function to get the reduction in AT090 w.r.t. the RV apex
#' activation in the HF cohort when pacing a single electrode in the LV wall in
#' any vein and any position.
#' 
#' @param output "%" or "ms" for the type of output (relative or absolute 
#' respectively).
#' 
#' @return A dataframe with fields "Vein", "Patient" (the number of the 
#' subject) and "Reduction" with the reduction in % or ms in AT090.

BaselineSingleElectrodeReduction <- function(cohort = "HF",
                                             output = "%",
                                             SA_folder = "default_HF_noPVTV",
                                             flag_debugging = FALSE){
  
  if(cohort == "RRHF"){
    single_El_RR <- BaselineSingleElectrodeReduction(cohort = "RR",
                                                     output = output,
                                                     SA_folder = SA_folder,
                                                     flag_debugging = flag_debugging)
    single_El_RR$Patient <- paste0("RR",single_El_RR$Patient)
    single_El_HF <- BaselineSingleElectrodeReduction(cohort = "HF",
                                                     output = output,
                                                     SA_folder = SA_folder,
                                                     flag_debugging = flag_debugging)
    single_El_HF$Patient <- paste0("HF",single_El_HF$Patient)
    
   single_El <- rbind(single_El_RR, single_El_HF)
  }
  else{
  if(cohort == "RR"){
    cohort <- "h"
  }
  RVapex <- Read_csv(file = paste0("/data/SA_multipole/",SA_folder,
                                   "/",cohort,"/multipole_RVapex.dat"), sep = "",
                     flag_debugging = flag_debugging) %>%
    pull(., var = "AT090")
  if(cohort == "HF"){
    names(RVapex) <- c(paste0("0",1:9),10:24)
  }
  else if(cohort == "h"){
    names(RVapex) <- c(paste0("0",1:9),10:20)
  }
  
  single_El <- as.data.frame(matrix(ncol = 3, nrow = 0))
  colnames(single_El) <- c("Patient","Vein","Reduction")
  
  
  for(heart in names(RVapex)){
    singleheart <- paste0("/data/SA_multipole/",SA_folder,"/",cohort,"/",heart,
                          "/multipole_AT090.dat") %>%
      Read_table(., header = TRUE, flag_debugging = flag_debugging)
    
    num_elec <- singleheart$lead %>% bin2ascii_lead() %>% ascii2bin_lead() %>%
                sum_string()
    one_electrode <- singleheart[which(num_elec == 1),2:6]
    
    reduction_vector <- RVapex[heart] - c(t(one_electrode) %>% as.matrix())
    
    if(output == "%"){
      reduction_vector <- 100 * reduction_vector / RVapex[heart]
    }
    
    veins <- c("AN","AL","LA","IL","IN")
      single_El <- rbind(single_El,
                         data.frame(Vein = rep(veins,8),
                                       Patient = rep(heart,40),
                                       Reduction = reduction_vector
                         )
                    )

    
  }
  }
  return(single_El)
  
}

#' @description Function to get the reduction (in % or ms) of all the options
#' of 1-electrode and 2-electrode pacings from a MPP design.

MPPDesignReduction <- function(cohort = "HF",
                           SA_folder = "default_HF_noPVTV",
                           MPP_design = "1348",
                           output = "%",
                           flag_debugging = FALSE){
  
  if(cohort == "RRHF"){
    MPP_1_2_RR <- MPPDesignReduction(cohort = "RR",
                                     SA_folder = SA_folder,
                                     MPP_design = MPP_design,
                                     output = output,
                                     flag_debugging = flag_debugging)
    MPP_1_2_RR$Patient <- paste0("RR",MPP_1_2_RR$Patient)
    MPP_1_2_HF <- MPPDesignReduction(cohort = "HF",
                                     SA_folder = SA_folder,
                                     MPP_design = MPP_design,
                                     output = output,
                                     flag_debugging = flag_debugging)
    MPP_1_2_HF$Patient <- paste0("HF",MPP_1_2_HF$Patient)
    
    MPP_1_2 <- rbind(MPP_1_2_RR, MPP_1_2_HF)
  }
  else{
  
  optimal_reductions_df <- c() 
  optimal_bipoles_df <- c()
  global_reductions_df <- c()
  
  
  cyphers <- strsplit(MPP_design,"")[[1]]
  
  #We extract all the bipoles from the optimal quadpole
  combs <- combn(4,2)
  
  vec_bipoles_from_quadpoles <- c()
  
  for(i in c(1:ncol(combs))){
    vec_bipoles_from_quadpoles[length(vec_bipoles_from_quadpoles) + 1] <-
      paste0(cyphers[combs[1,i]],cyphers[combs[2,i]])
  }
    
    MPP_1_2 <- as.data.frame(matrix(ncol = 3, nrow = 0))
    colnames(MPP_1_2) <- c("Patient","Vein","Reduction")
    patient_names <- c(paste0("0",1:9),10:20)
    
    if(cohort == "HF"){
      patient_names <- c(patient_names,21:24)
      if(SA_folder == "scar"){
        patient_names <- patient_names[-c(13,21)]
      }
    }
      
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein,
                                            which_cases = cohort, version = 4,
                                            output = output,
                                            SA_folder = SA_folder,
                                            flag_debugging = flag_debugging)
      if(cohort == "HF"){
        normalised_df <- normalised_cohorts$normalised_HF
      }
      else{
        normalised_df <- normalised_cohorts$normalised_RR
      }
      
      # Get only bipoles
      num_elec <- rownames(normalised_df) %>%
        bin2ascii_lead(.) %>%
        ascii2bin_lead(.)
      
      triquad_idx <- which(sum_string(num_elec) > 2)
      if(length(triquad_idx) > 0)
        normalised_df <- normalised_df[-triquad_idx,]
      normalised_ATs <- normalised_df
      
      normalised_ATs <- normalised_ATs[rownames(normalised_ATs) %in%
                                         ascii2bin_lead(vec_bipoles_from_quadpoles),]

        MPP_1_2 <- rbind(MPP_1_2,
                         data.frame(Vein = rep(vein,length(normalised_ATs)),
                                    Patient = rep(patient_names,nrow(normalised_ATs)),
                                    Reduction = c(t(normalised_ATs))))
    }
  }
    
    return(MPP_1_2)
}



#' @description Function to get the maximum reduction in AT090 w.r.t. the RV
#' apex activation in the HF cohort when pacing a single or two electrodes in
#' the LV wall.
#' 
#' @param output "%" or "ms" for the type of output (relative or absolute 
#' respectively).
#' 
#' @return A dataframe with fields "Vein", "Patient" (the number of the 
#' subject) and "Reduction" with the reduction in % or ms in AT090.

Optimal12ElectrodeReduction <- function(cohort = "HF", output = "%",
                                             SA_folder = "default_HF_noPVTV",
                                             flag_debugging = FALSE){
  
  hearts <- c(paste0("0",1:9),10:20)
  if(cohort == "HF"){
    hearts <- c(hearts,21:24)
    if(SA_folder == "scar"){
      hearts <- hearts[-c(13,21)]
    }
  }
  
  preoptimal <- All12ElectrodeReduction(cohort = cohort,
                                        metric_option = metric_option,
                                        output = output,
                                        SA_folder = SA_folder,
                                        flag_debugging = flag_debugging)
  
  optimal_df <- as.data.frame(matrix(ncol = length(hearts), nrow = 5))
  colnames(optimal_df) <- hearts
  rownames(optimal_df) <- c("AN", "AL", "LA", "IL", "IN")
  
  for(heart in colnames(optimal_df)){
    for(vein in rownames(optimal_df)){
      subset <- filter(preoptimal, Patient == heart & Vein == vein) 
      optimal_df[vein,heart] <- max(subset$Reduction)
    }
  }
  
  return(optimal_df)
  
}

#' @description Function that returns the value of the HAC heatmap as a data 
#' frame, with the option of returning in ms. To convert it back to matrix use
#' the function pivot_wider from tidyr.
All12ElectrodeReduction <- function(cohort = "HF", metric_option = "AT090",
                                    output = "%", 
                                    SA_folder = "default_HF_noPVTV",
                                    flag_debugging = FALSE){
  Load_Install_Packages("dplyr")
  
  if(cohort == "RR"){
    cohort <- "h"
  }
  
  RVapex <- Read_csv(file = paste0("/data/SA_multipole/",SA_folder,
                                      "/",cohort,"/multipole_RVapex.dat"),
                     sep = "", flag_debugging = flag_debugging) %>%
    pull(., var = metric_option)
 
  names(RVapex) <- c(paste0("0",1:9),10:length(RVapex))
  
  if(SA_folder == "scar"){
    names(RVapex) <- c(paste0("0",1:9),10:12,14:20,22:24)
  }
  
  optimal_df <- as.data.frame(matrix(ncol = 4, nrow = 0))
  colnames(optimal_df) <- c("Patient","Vein","Design","Reduction")
  preoptimal <- optimal_df
  
  for(heart in names(RVapex)){
    singleheart <- paste0("/data/SA_multipole/",SA_folder,"/",cohort,"/",heart,
                          "/multipole_",metric_option,".dat") %>%
      Read_table(., header = TRUE, flag_debugging = flag_debugging)
    
    num_elec <- singleheart$lead %>% bin2ascii_lead() %>% ascii2bin_lead() %>%
      sum_string()
    one_two_electrodes <- singleheart[which(num_elec < 3),2:6]
    
    reduction_vector <- RVapex[heart] - c(t(one_two_electrodes) %>%
                                               as.matrix())
    
    if(output == "%"){
      reduction_vector <- 100 * reduction_vector / RVapex[heart]
    }
    
    veins <- c("AN","AL","LA","IL","IN")
    preoptimal <- rbind(preoptimal,
                        data.frame(Vein = rep(veins,36),
                                   Patient = rep(heart,length(reduction_vector)),
                                   Design = bin2ascii_lead(singleheart$lead),
                                   Reduction = reduction_vector
                        )
    )
    
    
  }
  
  
  return(preoptimal)
  
}

#' @description Function to replicate the sensitivity analysis.
#' 
#' @param analysis "max_change" to get the maximum value of the HAC map for 
#' all the veins, and the difference with the default parameters.
#' "num_MPP_designs for a dataframe of the number of optimal MPP designs for
#' all the veins. 
#' "name_MPP_designs" for a dataframe of the name of the optimal MPP designs
#' across all the veins.
SensitivityAnalysis <- function(analysis, alpha = 0.01, flag_debugging = FALSE,
                                cohort, RV_midseptum = FALSE, scar = FALSE,
                                output = "%"){
  SA_folders <- c("default_HF_noPVTV",
                  "CV_007","CV_08",
                  "FEC_100","FEC_70",
                  "kFEC_10","kFEC_7",
                  "kxf_029","kxf_1")
  
  if(RV_midseptum){
    SA_folders <- c("default_HF_noPVTV",
                    "RVelec_midseptum")
  }
  if(scar){
    SA_folders <- c("default_HF_noPVTV",
                    "scar")
  }
  
  if(analysis == "max_change"){
    differences <- as.data.frame(matrix(nrow = (length(SA_folders)-1), ncol = 5))
    colnames(differences) <- c("SA_case", "Reduction_abs", "Difference_abs",
                               "Reduction_rel","Difference_rel")
    
    default_abs <- All12ElectrodeReduction(SA_folder = "default_HF_noPVTV",
                                           cohort = "HF",
                                           metric_option = "AT090",
                                           output = "ms")
    max_value_abs <- default_abs$Reduction %>% max()
    default_rel <- All12ElectrodeReduction(SA_folder = "default_HF_noPVTV",
                                           cohort = "HF",
                                           metric_option = "AT090",
                                           output = "%")
    max_value_rel <- default_rel$Reduction %>% max()
    
    for(i in c(1:(length(SA_folders)-1))){
      SA_case <- SA_folders[i+1]
      
      SA_case_abs <- All12ElectrodeReduction(SA_folder = SA_case,
                                             cohort = "HF",
                                             metric_option = "AT090",
                                             output = "ms")
      SA_case_rel <- All12ElectrodeReduction(SA_folder = SA_case,
                                             cohort = "HF",
                                             metric_option = "AT090",
                                             output = "%")
      
      differences[i,"SA_case"] <- SA_case
      differences[i,"Reduction_abs"] <- SA_case_abs$Reduction %>% max()
      differences[i,"Difference_abs"] <- abs(differences[i,"Reduction_abs"] - 
                                               max_value_abs)
      differences[i,"Reduction_rel"] <- SA_case_rel$Reduction %>% max()
      differences[i,"Difference_rel"] <- abs(differences[i,"Reduction_rel"] - 
                                               max_value_rel)
    }
    
    return(differences)
  }
  
  if(analysis == "num_MPP_designs"){
    num_MPP_designs <- as.data.frame(matrix(nrow = (length(SA_folders)),
                                            ncol = 4))
    colnames(num_MPP_designs) <- c("SA_case","HF","RR","RRHF")
    
    for(i in c(1:length(SA_folders))){
      SA_case <- SA_folders[i]
      num_MPP_designs[i, "SA_case"] <- SA_case
      
      HF_vec <-  c()
      RR_vec <-  c()
      RRHF_vec <-  c()
      for(vein in c("AN","AL","LA","IL","IN")){
        aux_HF <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                  vein = vein,
                                                  which_cases = "HF", 
                                                  method = "max",
                                                  SA_folder = SA_case, 
                                                  version = 4,
                                                  response = 100)
        HF_vec <- c(HF_vec, aux_HF$Quadripolar)
        
        aux_RR <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                  vein = vein,
                                                  which_cases = "RR", 
                                                  method = "max",
                                                  SA_folder = SA_case, 
                                                  version = 4,
                                                  response = 100)
        RR_vec <- c(RR_vec, aux_RR$Quadripolar)
        
        aux_RRHF <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                  vein = vein,
                                                  which_cases = "RRHF", 
                                                  method = "max",
                                                  SA_folder = SA_case, 
                                                  version = 4,
                                                  response = 100)
        RRHF_vec <- c(RRHF_vec, aux_RRHF$Quadripolar)
        }

      
      num_MPP_designs[i, "HF"] <- HF_vec %>% unique() %>% length()
      num_MPP_designs[i, "RR"] <- RR_vec %>% unique() %>% length()
      num_MPP_designs[i, "RRHF"] <- RRHF_vec %>% unique() %>% length()
    }
    
    return(num_MPP_designs)
  }
  
  if(analysis == "name_MPP_designs"){
    
    
    if(!scar){
      cohorts_vec <- c("HF","RR","RRHF")
    }
    else{
      cohorts_vec <- c("HF")
    }
    
    name_MPP_designs <- as.data.frame(matrix(nrow = (length(SA_folders)),
                                             ncol = length(cohorts_vec) + 1))
    colnames(name_MPP_designs) <- c("SA_case",cohorts_vec)
    
    for(i in c(1:length(SA_folders))){
      SA_case <- SA_folders[i]
      
      name_MPP_designs[i,"SA_case"] <- SA_case
      
      for(cohort in cohorts_vec){
      optimal_df <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                       vein = "ALL", 
                                                       which_cases = cohort,
                                                       method = "max",
                                                       SA_folder = SA_case,
                                                       version = 4,
                                                       flag_debugging = flag_debugging)
      
      name_MPP_designs[i,cohort] <- optimal_df$Quadripolar[1]
      }
    }
    
    return(name_MPP_designs)
  }
  
  if(analysis == "personalised_optimal"){

    Load_Install_Packages("tidyr")
    
    SA_df <- as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(SA_df) <- c("SA_case","Reduction")
    for(i in c(1:length(SA_folders))){
      SA_optimal <- Optimal12ElectrodeReduction(cohort = cohort, output = "%",
                                                SA_folder = SA_folders[i],
                                                flag_debugging = flag_debugging)
      
      SA_df <- rbind(SA_df,data.frame(Reduction = c(t(SA_optimal)),
                                      SA_case = rep(SA_folders[i],
                                                    length(c(t(SA_optimal))))))
    }
   test <- pairwise.wilcox.test(x = SA_df$Reduction, g = SA_df$SA_case,
                                p.adjust.method = "bonferroni", paired = FALSE,
                                alternative = "two.sided", alpha = alpha)
   
   means_vec <- pivot_wider(data = SA_df, names_from = SA_case,
                            values_from = Reduction) %>% 
     apply(., 2, function(x) mean(unlist(x)))
   
   return(list(paste0("alpha: ",alpha),
               test$p.value[,1] < alpha,
               means_vec)
          )
  }
  
  if(analysis == "cohort_optimal"){
    
    SA_df <- as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(SA_df) <- c("SA_case","Reduction")
    for(i in c(1:length(SA_folders))){
      optimal_df <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                       vein = "ALL", 
                                                       which_cases = cohort,
                                                       method = "max",
                                                       SA_folder = SA_folders[i],
                                                       version = 4)
      
      name_MPP_design <- optimal_df$Quadripolar[1]
      
      SA_optimal <- MPPDesignReduction(cohort = cohort,
                                       MPP_design = name_MPP_design,
                                       output = output,
                                       SA_folder = SA_folders[i])
      
      SA_df <- rbind(SA_df,data.frame(Reduction = SA_optimal$Reduction,
                                      SA_case = rep(SA_folders[i],
                                                    length(SA_optimal$Reduction)
                                                    )
                                      )
                     )
    }
    test <- pairwise.wilcox.test(x = SA_df$Reduction, g = SA_df$SA_case,
                                 p.adjust.method = "bonferroni", paired = FALSE,
                                 alternative = "two.sided", alpha = alpha)
    
    means_vec <- pivot_wider(data = SA_df, names_from = SA_case,
                             values_from = Reduction) %>% 
      apply(., 2, function(x) mean(unlist(x)))
    
    return(list(paste0("alpha: ",alpha),
                test$p.value[,1],
                means_vec)
    )
  }
}

#' @description Function to ease the replication of the p-values obtained for
#' the paper.
#' 
#' @param comparison "design" to assess if a given design has a significant 
#' reduction in AT090.
#' "design_vs_personalised" to assess the differences between a given MPP lead
#' design and choosing the personalised value.
#' @param MPP_design A four-digits lead design.
#' @param cohort "HF" or "RR" for heart failure or reverse remodelled, 
#' respectively.
#' @param output "ms" or "%" to express the result in absolute or relative
#' terms.
#' @param alpha Confidence level for the Wilcoxon / Mann-Whitney tests.
#' @param alternative "g", "l" or "t" for grater, less, or two-sided in the 
#' comparison.
#' @param SA_folder Name of the simulations folder
#' @param metric_option Metric to assess. AT090 and TAT are the two more used.
#' @param flag_debugging If TRUE, prints whatever is reading and writing.
#' 
#' @return The p-value and the means of whatever is compared.

Getpvalues <- function(comparison, MPP_design, cohort, output, alpha = 0.01, 
                       alternative, SA_folder = "default_HF_noPVTV",
                       metric_option = "AT090", flag_debugging = FALSE){
  

  if(comparison == "allpacing"){
    df <- All12ElectrodeReduction(cohort = cohort,
                                   metric_option = metric_option,
                                   output = output, SA_folder = SA_folder,
                                   flag_debugging = flag_debugging)
    
    mean_value <- mean(df$Reduction) %>% round(.,2)
    sd_value <- sd(df$Reduction) %>% round(.,2)
    
    test <- wilcox.test(x = df$Reduction, alternative = "greater", mu = 0,
                        exact = FALSE, alpha = alpha)
    
    return(list(paste0("Mean: ", mean_value),
                paste0("SD: ", sd_value),
                paste0("p-value: ",test$p.value)
                )
           )
  }
  if(comparison == "design"){
    design_all <- MPPDesignReduction(cohort = cohort, MPP_design = MPP_design,
                                     output = output,
                                     SA_folder = SA_folder)
    design_vec <- design_all$Reduction
    
    test <- wilcox.test(x = design_vec, alternative = "greater", mu = 0,
                        exact = FALSE, alpha = alpha)
    
    return(list(paste0("p-value: ",test$p.value,2),
                paste0("Mean value: ",
                       round(mean(design_vec),2)),
                paste0("SD: ",
                       round(sd(design_vec),2))
                )
    )
  }
  
  if(comparison == "design_vs_personalised"){
  
    personalised_veins <- Optimal12ElectrodeReduction(cohort = cohort,
                                                      output = output,
                                                      SA_folder = SA_folder)
    personalised_vec <- apply(personalised_veins,2,max)
    
    
    design_all <- MPPDesignReduction(cohort = cohort, MPP_design = MPP_design,
                                     output = output,
                                     SA_folder = SA_folder)
    design_vec <- design_all$Reduction

    test <- wilcox.test(x = personalised_vec, y = design_vec, 
                        alternative = "greater", paired = FALSE, exact = FALSE,
                        alpha = alpha)
    
    return(list(paste0("p-value: ",test$p.value,2),
                paste0("Difference of mean values: ",
                       round(mean(personalised_vec) - mean(design_vec),2))
    )
    )
  }
  
  if(comparison == "global_vs_personalised"){
    
    optimal_df <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                     vein = "ALL", 
                                                     which_cases = "HF",
                                                     method = "max",
                                                     SA_folder = SA_folder,
                                                     version = 4)
    
    name_MPP_design <- optimal_df$Quadripolar[1]
    
    global_optimal_HF <- MPPDesignReduction(cohort = "HF",
                                            MPP_design = name_MPP_design,
                                            output = output,
                                            SA_folder = SA_folder)
    personalised_optimal_HF <- Optimal12ElectrodeReduction(cohort = "HF",
                                                           output = output,
                                                           SA_folder = SA_folder)
    
    if(SA_folder != "scar"){
      optimal_df <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                       vein = "ALL", 
                                                       which_cases = "RR",
                                                       method = "max",
                                                       SA_folder = SA_folder,
                                                       version = 4,
                                                       flag_debugging = flag_debugging)
      
      name_MPP_design <- optimal_df$Quadripolar[1]
      
      
      global_optimal_RR <- MPPDesignReduction(cohort = "RR",
                                              MPP_design = name_MPP_design,
                                              output = output,
                                              SA_folder = SA_folder)
      personalised_optimal_RR <- Optimal12ElectrodeReduction(cohort = "RR",
                                                             output = output,
                                                             SA_folder = SA_folder)
    }
    global_HF_vec <- global_optimal_HF$Reduction
    person_HF_vec <- c(t(personalised_optimal_HF))
    
    if(SA_folder != "scar"){
      global_RR_vec <- global_optimal_RR$Reduction
      person_RR_vec <- c(t(personalised_optimal_RR))
    }
    
    final_df <- as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(final_df) <- c("Reduction", "Strategy")
    
    final_df <- rbind(final_df,data.frame(Reduction = global_HF_vec,
                                          Strategy = rep("HF-based MPP design used in HF",
                                                         length(global_HF_vec))))
    final_df <- rbind(final_df,data.frame(Reduction = person_HF_vec,
                                          Strategy = rep("Personalised MPP design used in HF",
                                                         length(person_HF_vec))))
    if(SA_folder != "scar"){
      final_df <- rbind(final_df,data.frame(Reduction = global_RR_vec,
                                            Strategy = rep("RR-based MPP design used in RR",
                                                           length(global_RR_vec))))
      final_df <- rbind(final_df,data.frame(Reduction = person_RR_vec,
                                            Strategy = rep("Personalised MPP design used in RR",
                                                           length(person_RR_vec))))
    }
    
    test <- pairwise.wilcox.test(x = final_df$Reduction, final_df$Strategy,
                                 p.adjust.method = "bonferroni",
                                 alternative = alternative, alpha = alpha)
    if(SA_folder != "scar"){
      return_list <- list(paste0("alpha: ",alpha),
                          test,
                          paste0("Mean of HF-pop: ",mean(global_HF_vec)),
                          paste0("Mean of HF-pers: ",mean(person_HF_vec)),
                          paste0("Mean of RR-pop: ",mean(global_RR_vec)),
                          paste0("Mean of RR-pers: ",mean(person_RR_vec))
      )
    }
    else{
      return_list <- list(paste0("alpha: ",alpha),
                          test,
                          paste0("Mean of HF-pop: ",mean(global_HF_vec)),
                          paste0("Mean of HF-pers: ",mean(person_HF_vec))
      )
    }
    return(return_list)
  }
  
  if(comparison == "optimal_vs_cohort"){
    optimal_df <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                     vein = "ALL", 
                                                     which_cases = "RR",
                                                     method = "max",
                                                     SA_folder = SA_folder,
                                                     version = 4)
    
    name_MPP_design_RR <- optimal_df$Quadripolar[1]
    optimal_df <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                     vein = "ALL", 
                                                     which_cases = "HF",
                                                     method = "max",
                                                     SA_folder = SA_folder,
                                                     version = 4)
    
    name_MPP_design_HF <- optimal_df$Quadripolar[1]
    
    
    optimal_HF_in_HF <- MPPDesignReduction(cohort = "HF",
                                           MPP_design = name_MPP_design_HF,
                                           output = output,
                                           SA_folder = SA_folder)
    optimal_HF_in_RR <- MPPDesignReduction(cohort = "RR",
                                           MPP_design = name_MPP_design_HF,
                                           output = output,
                                           SA_folder = SA_folder)
    optimal_RR_in_HF <- MPPDesignReduction(cohort = "HF",
                                           MPP_design = name_MPP_design_RR,
                                           output = output,
                                           SA_folder = SA_folder)
    optimal_RR_in_RR <- MPPDesignReduction(cohort = "RR",
                                           MPP_design = name_MPP_design_RR,
                                           output = output,
                                           SA_folder = SA_folder)
    
    final_df <- as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(final_df) <- c("Reduction", "Strategy")
    
    final_df <- rbind(final_df,data.frame(Reduction = optimal_HF_in_HF$Reduction,
                                          Strategy = rep("HF-based MPP design used in HF",
                                                         nrow(optimal_HF_in_HF))))
    final_df <- rbind(final_df,data.frame(Reduction = optimal_HF_in_RR$Reduction,
                                          Strategy = rep("HF-based MPP design used in RR",
                                                         nrow(optimal_HF_in_RR))))
    final_df <- rbind(final_df,data.frame(Reduction = optimal_RR_in_HF$Reduction,
                                          Strategy = rep("RR-based MPP design used in HF",
                                                         nrow(optimal_RR_in_HF))))
    final_df <- rbind(final_df,data.frame(Reduction = optimal_RR_in_RR$Reduction,
                                          Strategy = rep("RR-based MPP design used in RR",
                                                         nrow(optimal_RR_in_RR))))
    
    test <- pairwise.wilcox.test(x = final_df$Reduction, final_df$Strategy,
                                 p.adjust.method = "bonferroni",
                                 alternative = alternative, alpha = alpha)
    return(list(paste0("alpha: ",alpha),
                test,
                paste0("Mean of HF in HF: ",mean(optimal_HF_in_HF$Reduction)),
                paste0("Mean of HF in RR: ",mean(optimal_HF_in_RR$Reduction)),
                paste0("Mean of RR in HF: ",mean(optimal_RR_in_HF$Reduction)),
                paste0("Mean of RR in RR: ",mean(optimal_RR_in_RR$Reduction))
    )
    )
  }
}