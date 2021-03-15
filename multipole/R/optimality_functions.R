
find_optimal_quadripole <- function(vein,which_cases,method, flag_debugging = FALSE){
  
  quadripole <- c()
  if(vein != "ALL"){
    # We get the names of the optimal lead designs in monopoles and dipoles
    bipolar_names_score <- Find_optimal_monodipoles(vein = vein, response = 100, which_cases = which_cases, version = 1, output = "both")
    
    bipolar_names <- bipolar_names_score[[1]]
    bipolar_scores <- bipolar_names_score[[2]]
    
    
    
    #We get the improvement for each design
    bipolar_scores <- c(bipolar_scores[1], bipolar_scores[-1] - bipolar_scores[-length(bipolar_scores)])
    
    
  }
  else{
    hierarchy_veins <- c("LA","IL","AL","IN","AN")
    bipolar_names <- c()
    bipolar_scores <- c()
    for(each_vein in hierarchy_veins){
      # We get the names of the optimal lead designs in monopoles and dipoles
      bipolar_names_score_temp <- Find_optimal_monodipoles(vein = each_vein, response = 100, which_cases = which_cases, version = 1, output = "both")
      
      bipolar_names_temp <- bipolar_names_score_temp[[1]]
      bipolar_scores_temp <- bipolar_names_score_temp[[2]]
      
      
      
      #We get the improvement for each design
      bipolar_scores_temp <- c(bipolar_scores_temp[1], bipolar_scores_temp[-1] - bipolar_scores_temp[-length(bipolar_scores_temp)])
      
      bipolar_names_temp <- bipolar_names_temp[order(bipolar_scores_temp,decreasing = TRUE)]
      bipolar_scores_temp <- sort(bipolar_scores_temp,decreasing = TRUE)
      
      bipolar_names <- c(bipolar_names,bipolar_names_temp)
      bipolar_scores <- c(bipolar_scores,bipolar_scores_temp)
    }
    
    # We have to create a unique data frame 
    
    df_wrapped <- data.frame(names = bipolar_names, scores = bipolar_scores)
    final_df <- data.frame(names = unique(bipolar_names))
    
    for(i in seq(1,nrow(final_df))){
      vector_of_patients <- df_wrapped[df_wrapped$names == final_df$names[i],2]
      #final_df[i,2] <- max(vector_of_patients)
      # The idea is to do the same as above but using a generic function given by the argument method
      body <- paste0(method,"(v)")
      args <- "v"
      eval(parse(text = paste('f <- function(', args, ') { return(' , body , ')}', sep=''))) # The generic function is renamed as f(v)
      
      final_df[i,2] <- f(vector_of_patients)
    }
    rm(bipolar_names,bipolar_scores)
    
    bipolar_names <- as.vector(final_df$names)
    bipolar_scores <- final_df[,2]
    
  }
  
  # We sort them according to the number of patients improved for each bipole because the quadripoles will be chosen based on this.
  
  bipolar_names <- bipolar_names[order(bipolar_scores,decreasing = TRUE)]
  bipolar_scores <- sort(bipolar_scores,decreasing = TRUE)
  
  bipolar_included <- 0*bipolar_scores
  i <- 1
  j <- 2
  
  # While there are bipolars not included in the quadripolar
  
  
  while(sum(bipolar_included) < 2*length(bipolar_included)){
    
    quadripolar_name <- merge_design(c(bipolar_names[i],bipolar_names[j]))
    
    bipolar_included[i] <- 2
    bipolar_included[j] <- 2
    
    #Update the ones included after merging
    for(k in c(1:length(bipolar_names))){
      if(bipolar_included[k] < 2){
        bipolar_included[k] <- check_quadripole_includes_dipole(quadripolar_name,bipolar_names[k])
      }
    }
    
    
    
    keepgoing = TRUE
    while(keepgoing){
      anychange = FALSE
      for(k in c(1:length(bipolar_names))){
        if(((sum_string(quadripolar_name) == 2) && bipolar_included[k] < 2 )|| ((sum_string(quadripolar_name) == 3) && bipolar_included[k] == 1)){
          quadripolar_name <- merge_design(c(quadripolar_name,bipolar_names[k]))
          bipolar_included[k] <- 2
          anychange = TRUE
          break
        }
        
      }
      
      if(sum_string(quadripolar_name) == 4 || (!anychange)){
        keepgoing = FALSE
      }
      
    }
    
    
    #Update the ones included after merging
    for(k in c(1:length(bipolar_names))){
      if(bipolar_included[k] < 2){
        bipolar_included[k] <- check_quadripole_includes_dipole(quadripolar_name,bipolar_names[k])
      }
    }
    # Wether if it's a tripole or not, we add it
    quadripole[length(quadripole) + 1] <- quadripolar_name
    
    # We find the next electrode to take into account
    first_found = FALSE
    second_found = FALSE
    
    for( k in c(1:length(bipolar_included))){
      if(bipolar_included[k] < 2){
        if(first_found == FALSE){
          i <- k
          first_found = TRUE
        }
        else if(first_found == TRUE){
          j <- k
          second_found = TRUE
          break
        }
      }
    }
    
    if(second_found == FALSE){
      if(check_quadripole_includes_dipole(quadripolar_name,bipolar_names[length(bipolar_names)]) < 2){
        quadripole[length(quadripole) + 1] <- bipolar_names[length(bipolar_names)]
      }
      bipolar_included[length(bipolar_included)] <- 2
    }
    
    for(k in c(1:length(bipolar_names))){
      if(bipolar_included[k] < 2){
        bipolar_included[k] <- check_vector_includes_dipole(quadripole,bipolar_names[k])
      }
    }
    
  }
  
  # Finally, some quadripoles might still be bipoles or quadripoles, so we repeat the process
  
  for(i in seq(1,length(quadripole))){
    if(sum_string(quadripole[i]) < 4){
      quadripole[i] <- complete_quadripole(almost_quadripole = quadripole[i], bipolar_names = bipolar_names)
    }
  }
  
  
  quadripole_ascii <- bin2ascii_lead(quadripole)
  
  bipolar_names_ascii <- bin2ascii_lead(bipolar_names)
  return(list("Patient improvement"=bipolar_scores,"Bipolar"=bipolar_names_ascii,"Quadripolar"=unique(quadripole_ascii)))
}

find_optimal_quadripole_combinatorial <- function(vein,which_cases,method, flag_debugging = FALSE){
  
  quadripole <- c()
  if(vein != "ALL"){
    # We get the names of the optimal lead designs in monopoles and dipoles
    bipolar_names_score <- Find_optimal_monodipoles(vein = vein, response = 100, which_cases = which_cases, version = 1, output = "both")
    
    bipolar_names <- bipolar_names_score[[1]]
    bipolar_scores <- bipolar_names_score[[2]]
    
    #We get the improvement for each design
    bipolar_scores <- c(bipolar_scores[1], bipolar_scores[-1] - bipolar_scores[-length(bipolar_scores)])
    
    
  }
  else{
    hierarchy_veins <- c("LA","IL","AL","IN","AN")
    bipolar_names <- c()
    bipolar_scores <- c()
    for(each_vein in hierarchy_veins){
      # We get the names of the optimal lead designs in monopoles and dipoles
      bipolar_names_score_temp <- Find_optimal_monodipoles(vein = each_vein, response = 100, which_cases = which_cases, version = 1, output = "both")
      
      bipolar_names_temp <- bipolar_names_score_temp[[1]]
      bipolar_scores_temp <- bipolar_names_score_temp[[2]]
      
      
      
      #We get the improvement for each design
      bipolar_scores_temp <- c(bipolar_scores_temp[1], bipolar_scores_temp[-1] - bipolar_scores_temp[-length(bipolar_scores_temp)])
      
      bipolar_names_temp <- bipolar_names_temp[order(bipolar_scores_temp,decreasing = TRUE)]
      bipolar_scores_temp <- sort(bipolar_scores_temp,decreasing = TRUE)
      
      bipolar_names <- c(bipolar_names,bipolar_names_temp)
      bipolar_scores <- c(bipolar_scores,bipolar_scores_temp)
    }
    
    # We have to create a unique data frame 
    
    df_wrapped <- data.frame(names = bipolar_names, scores = bipolar_scores)
    final_df <- data.frame(names = unique(bipolar_names))
    
    for(i in seq(1,nrow(final_df))){
      vector_of_patients <- df_wrapped[df_wrapped$names == final_df$names[i],2]
      #final_df[i,2] <- max(vector_of_patients)
      # The idea is to do the same as above but using a generic function given by the argument method
      body <- paste0(method,"(v)")
      args <- "v"
      eval(parse(text = paste('f <- function(', args, ') { return(' , body , ')}', sep=''))) # The generic function is renamed as f(v)
      
      final_df[i,2] <- f(vector_of_patients)
    }
    rm(bipolar_names,bipolar_scores)
    
    bipolar_names <- as.vector(final_df$names)
    bipolar_scores <- final_df[,2]
    
  }
  
  
  # We create all the combinations of 4 electrodes from 1 to 8
  matrix_comb <- t(combn(1:8,4))
  
  vec_comb <- apply(matrix_comb,1,function(x) paste0(x,collapse = ''))
  
  quad_scores <- 0*c(1:length(vec_comb))
  
  # We create a vector with the scores of all the quadripoles
  for(i in seq(1:length(vec_comb))){
    for(j in seq(1:length(bipolar_names))){
      if(check_quadripole_includes_dipole(ascii2bin_lead(vec_comb[i]),bipolar_names[j]) == 2){
        quad_scores[i] = quad_scores[i] + bipolar_scores[j]
      }
    }
  }
  
  # We sort the optimal quadripoles
  quad_combin_sorted <- vec_comb[order(quad_scores, decreasing = TRUE)]
  
  # We initialise a vector saying if the bipoles are taking into account...
  isbipoleincluded <- FALSE*c(1:length(bipolar_names))
  # ... and a vector saying if the quadripole is in the result
  isquadripoleinresult <- FALSE*c(1:length(quad_combin_sorted))
  res_quadpoles <- c()
  
  for(i in seq(1:length(bipolar_names))){ # For each bipole
    if(!(isbipoleincluded[i])){ # If that bipole is not included
      for(j in seq(1:length(quad_combin_sorted))){ # For each quadpole
        if(!(isquadripoleinresult[j])){ # If that quadpole is not already considered
          if(check_quadripole_includes_dipole(ascii2bin_lead(quad_combin_sorted[j]),bipolar_names[i]) == 2){ # If the quadpole includes the bipole
            res_quadpoles[length(res_quadpoles) + 1] <- quad_combin_sorted[j]
            isquadripoleinresult[j] <- TRUE
            isbipoleincluded[i] <- TRUE
            # And we discard all the next bipoles that are in that quadpole
            for(k in c((i+1):(length(bipolar_names)))){
              if(check_quadripole_includes_dipole(ascii2bin_lead(quad_combin_sorted[j]),bipolar_names[k]) == 2){
                isbipoleincluded[k] <- TRUE
              }
            }
            break
          }
        }
      }
    }
  }
  
  # We sort them according to the number of patients improved for each bipole because the quadripoles will be chosen based on this.
  
  bipolar_names <- bipolar_names[order(bipolar_scores,decreasing = TRUE)]
  bipolar_scores <- sort(bipolar_scores,decreasing = TRUE)
  
  bipolar_included <- 0*bipolar_scores
  i <- 1
  j <- 2
  
  
  quadripole_ascii <- res_quadpoles
  
  bipolar_names_ascii <- bin2ascii_lead(bipolar_names)
  return(list("Patient improvement"=bipolar_scores,"Bipolar"=bipolar_names_ascii,"Quadripolar"=unique(quadripole_ascii)))
}

Find_Optimal_Quadripole_Evenly <- function(vein,which_cases,method, flag_debugging = FALSE){
  
  quadripole <- c()
  if(vein != "ALL"){
    # We get the names of the optimal lead designs in monopoles and dipoles
    bipolar_names_score <- Find_optimal_monodipoles(vein = vein, response = 100, which_cases = which_cases, version = 1, output = "both")
    
    bipolar_names <- bipolar_names_score[[1]]
    bipolar_scores <- bipolar_names_score[[2]]
    
    #We get the improvement for each design
    bipolar_scores <- c(bipolar_scores[1], bipolar_scores[-1] - bipolar_scores[-length(bipolar_scores)])
    
    
  }
  else{
    hierarchy_veins <- c("LA","IL","AL","IN","AN")
    bipolar_names <- c()
    bipolar_scores <- c()
    for(each_vein in hierarchy_veins){
      # We get the names of the optimal lead designs in monopoles and dipoles
      bipolar_names_score_temp <- Find_optimal_monodipoles(vein = each_vein, response = 100, which_cases = which_cases, version = 1, output = "both")
      
      bipolar_names_temp <- bipolar_names_score_temp[[1]]
      bipolar_scores_temp <- bipolar_names_score_temp[[2]]
      
      #We get the improvement for each design
      bipolar_scores_temp <- c(bipolar_scores_temp[1], bipolar_scores_temp[-1] - bipolar_scores_temp[-length(bipolar_scores_temp)])
      
      bipolar_names_temp <- bipolar_names_temp[order(bipolar_scores_temp,decreasing = TRUE)]
      bipolar_scores_temp <- sort(bipolar_scores_temp,decreasing = TRUE)
      
      bipolar_names <- c(bipolar_names,bipolar_names_temp)
      bipolar_scores <- c(bipolar_scores,bipolar_scores_temp)
    }
    
    # We have to create a unique data frame 
    
    df_wrapped <- data.frame(names = bipolar_names, scores = bipolar_scores)
    final_df <- data.frame(names = unique(bipolar_names))
    
    for(i in seq(1,nrow(final_df))){
      vector_of_patients <- df_wrapped[df_wrapped$names == final_df$names[i],2]
      #final_df[i,2] <- max(vector_of_patients)
      # The idea is to do the same as above but using a generic function given by the argument method
      body <- paste0(method,"(v)")
      args <- "v"
      eval(parse(text = paste('f <- function(', args, ') { return(' , body , ')}', sep=''))) # The generic function is renamed as f(v)
      
      final_df[i,2] <- f(vector_of_patients)
    }
    rm(bipolar_names,bipolar_scores)
    
    bipolar_names <- as.vector(final_df$names)
    bipolar_scores <- final_df[,2]
    
  }
  
  
  # We create all the combinations of 4 electrodes from 1 to 8
  matrix_comb <- t(combn(1:8,4))
  
  vec_comb <- apply(matrix_comb,1,function(x) paste0(x,collapse = ''))
  
  quad_scores <- 0*c(1:length(vec_comb))
  
  # We create a vector with the scores of all the quadripoles
  for(i in seq(1:length(vec_comb))){
    for(j in seq(1:length(bipolar_names))){
      if(check_quadripole_includes_dipole(ascii2bin_lead(vec_comb[i]),bipolar_names[j]) == 2){
        quad_scores[i] = quad_scores[i] + bipolar_scores[j]
      }
    }
  }
  
  # We sort the optimal quadripoles
  quad_combin_sorted <- vec_comb[order(quad_scores, decreasing = TRUE)]
  # We initialise a vector saying if the bipoles are taking into account...
  isbipoleincluded <- FALSE*c(1:length(bipolar_names))
  # ... and a vector saying if the quadripole is in the result
  isquadripoleinresult <- FALSE*c(1:length(quad_combin_sorted))
  isquadpartialincluded <- isquadripoleinresult
  res_quadpoles <- c()
  
  for(bip_num in c(1:length(bipolar_names))){ # For each bipole
    while(!(isbipoleincluded[bip_num])){ # If that bipole is not included
      for(quad_num in c(1:length(quad_combin_sorted))){ # For each quadpole
        if((prod(isquadpartialincluded))){
          isquadpartialincluded = !(isquadpartialincluded)
        }
        if(!(isquadripoleinresult[quad_num]) && !(isquadpartialincluded[quad_num])){ # If that quadpole is not already considered
          if(check_quadripole_includes_dipole(ascii2bin_lead(quad_combin_sorted[quad_num]),bipolar_names[bip_num]) == 2){ # If the quadpole includes the bipole
            res_quadpoles[length(res_quadpoles) + 1] <- quad_combin_sorted[quad_num]
            isquadripoleinresult[quad_num] <- TRUE
            isquadpartialincluded[quad_num] <- TRUE
            isbipoleincluded[bip_num] <- TRUE
            
            
            # And we discard all the next bipoles that are in that quadpole
            for(k in c((bip_num+1):(length(bipolar_names)))){
              if(check_quadripole_includes_dipole(ascii2bin_lead(quad_combin_sorted[quad_num]),bipolar_names[k]) == 2){
                isbipoleincluded[k] <- TRUE
              }
            }
            for(k in c(1:length(quad_combin_sorted))){
              if(!(isquadpartialincluded[k])){
                for(m in c(1:length(isbipoleincluded))){
                  if((isbipoleincluded[m]) && (check_quadripole_includes_dipole(ascii2bin_lead(quad_combin_sorted[k]),bipolar_names[m]) == 2)){
                    isquadpartialincluded[k] <- TRUE
                  }
                }
              }
            }
            break
          }
        }
        if(quad_num == length(quad_combin_sorted)){
          isquadpartialincluded = 1 + 0*isquadpartialincluded
        }
      }
    }
  }
  
  # We sort them according to the number of patients improved for each bipole because the quadripoles will be chosen based on this.
  
  bipolar_names <- bipolar_names[order(bipolar_scores,decreasing = TRUE)]
  bipolar_scores <- sort(bipolar_scores,decreasing = TRUE)
  
  bipolar_included <- 0*bipolar_scores
  i <- 1
  j <- 2
  
  
  quadripole_ascii <- res_quadpoles
  quadripolar_score <- quad_scores[match(quadripole_ascii,vec_comb)]
  
  bipolar_names_ascii <- bin2ascii_lead(bipolar_names)
  return(list("Patient improvement"=bipolar_scores,"Bipolar"=bipolar_names_ascii,"Quadripolar"=unique(quadripole_ascii),"Quadripolar_Scores"=quadripolar_score))
}

#' @description Gets the optimal dipoles (and monopoles) from a file of 
#' activation times according to a specific metric.
#' 
#' @param metric_option "TAT", "AT1090", "AT090", "VEUTAT", "VEUmean", "TATLV",
#' for each of the metrics. In the latest version, we do not use TATLV anymore.
#' @param vein Vein to check the optimal config. Default is "LA". Other options
#' are "AN", "AL", "IL", "IN".
#' @param response Threshold to include the "best" activations. Less than 100 
#' means that that percentil is included (if 80 then the electrodes that get the
#' 80% of the best response are selected). Default is 100.
#' @param which_cases Condition of the patients. "RR" for reverse remodelled, 
#' "HF" for heart failure. Default is "HF".
#' @param version Number of the version to read the file. Latest is 4 on July 
#' 2020.
#' @param output What to give as output. Options are "number", which returns the
#' number of patients improved with the best designs; "names" return which are 
#' the best designs and "both" returns both outputs. Default is "number".
#' @param SA_folder Name of the folder of the simulations.
#' @param flag_debugging Flag (TRUE or FALSE) to print the name and location of
#' the files read and written.
#' 
#' @return See the output parameter.

Find_optimal_monodipoles <- function(metric_option, vein, response = 100,
                                     which_cases, version, output = "number",
                                     SA_folder, flag_debugging = FALSE){
                           
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr"))
  
  
  # We read the dataframe
  normalised_cohorts <- Read_dataframes(metric_option = metric_option, 
                                        lead_option = vein, 
                                        which_cases = which_cases,
                                        version = version,
                                        SA_folder = SA_folder,
                                        flag_debugging = flag_debugging)
  
  normalised_RR <- normalised_cohorts$normalised_RR
  normalised_HF <- normalised_cohorts$normalised_HF
  
  if(which_cases == "RR" || which_cases == "RRHF"){
    num_elec <- rownames(normalised_RR) %>%
      bin2ascii_lead() %>% ascii2bin_lead() %>% sum_string()
    triquad_idx <- which(num_elec != 2)
    
    if(length(triquad_idx) > 0)
      normalised_RR <- normalised_RR[-triquad_idx,]
  }
  if(which_cases == "HF" || which_cases == "RRHF"){
    num_elec <- rownames(normalised_HF) %>%
      bin2ascii_lead() %>% ascii2bin_lead() %>% sum_string()
    triquad_idx <- which(num_elec != 2)    
    
    if(length(triquad_idx) > 0)
      normalised_HF <- normalised_HF[-triquad_idx,]
  }
  
  if(which_cases == "RRHF"){
    normalised_ATs <- cbind(normalised_RR,normalised_HF)
  }
  if(which_cases == "RR"){
    normalised_ATs <- normalised_RR
  }
  if(which_cases == "HF"){
    normalised_ATs <- normalised_HF
  }
  
  num_patients <- c()
  names_best_leads <- c()
  
  
  while(sum(sum(!is.na(normalised_ATs))) > 0){
    
    patient_quantiles <- apply(normalised_ATs,2,
                               function(x) quantile(x[!is.na(x)],
                                                    probs = response/100.))
    
    binary_df <- apply(normalised_ATs,1,function(x) x >= patient_quantiles)
    binary_df <- t(binary_df)
    
    # For each lead, the number of patients over the percentile specified
    # before.
    patients_improved <- apply(binary_df,1,function(x) sum(x[!is.na(x)])) 
    
    lead2remove <- which.max(patients_improved)
    pos2add <- length(names_best_leads)+1
    names_best_leads[pos2add] <- rownames(normalised_ATs)[lead2remove]
    normalised_ATs[lead2remove,] <- NA # Remove the best lead
    
    patients2remove <- names(which(binary_df[lead2remove,]))
    normalised_ATs[,patients2remove] <- NA #Remove the patients
    
    num_patients <- c(num_patients,max(patients_improved))
    
  }
  
  num_patients <- cumsum(num_patients)
  
  if(output == "number")
    return(num_patients)
  else if(output == "names")
    return(names_best_leads)
  else if(output == "both")
    return(list(names_best_leads,num_patients))
  
}

#' @description Function to find the optimal quadripoles covering the maximum
#' number of patients for which it's optimal. It chooses the minimal number
#' of quadripolar leads.
#' 
#' @param metric_option "TAT", "AT1090", "AT090", "VEUTAT", "VEUmean", "LVTAT"
#' for each of the metrics. In the latest version, we do not use TATLV anymore. 
#' @param vein Vein to check the optimal config. Default is "LA". Other options
#' are "AN", "AL", "IL", "IN".
#' @param which_cases Condition of the patients. "RR" for reverse remodelled, 
#' "HF" for heart failure. Default is "HF".
#' @param method Function to choose the number of patients. Usually is "max",
#' so we optimise to get the maximum number of patients optimised.
#' @param response Threshold to include the "best" activations. Less than 100 
#' means that that percentil is included (if 80 then the electrodes that get the
#' 80% of the best response are selected). Default is 100.
#' @param SA_folder Name of the folder of the simulations.
#' @param version Number of the version to read the file.
#' @param flag_debugging Flag (TRUE or FALSE) to print the name and location of
#' the files read and written.
#' 
#' @return A list whose components are the number of patients improved with
#' the best bipolar designs; the name of the best bipolar designs; the name of
#' the optimal quadripolar designs and the number of patients improved with the 
#' best quadripolar design.

Find_Optimal_Quadripole_Optimising <- function(metric_option, vein, which_cases,
                                               method, response = 100, 
                                               SA_folder, version, 
                                               flag_debugging = FALSE){
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages("dplyr")
  
  quadripole <- c()
  if(vein != "ALL"){
    # We get the names of the optimal lead designs in monopoles and dipoles
    bipolar_names_score <- Find_optimal_monodipoles(metric_option = metric_option,
                                                    vein = vein, 
                                                    response = response, 
                                                    which_cases = which_cases,
                                                    version = version,
                                                    output = "both", 
                                                    SA_folder = SA_folder,
                                                    flag_debugging = flag_debugging)
    
    bipolar_names <- bipolar_names_score[[1]]
    bipolar_scores <- bipolar_names_score[[2]]
    
    #We get the improvement for each design
    n_bs <- length(bipolar_scores)
    bipolar_scores <- c(bipolar_scores[1],
                        bipolar_scores[-1] - bipolar_scores[-n_bs])
    
    
  }
  else{
    hierarchy_veins <- c("LA","IL","AL","IN","AN")
    bipolar_names <- c()
    bipolar_scores <- c()
    for(each_vein in hierarchy_veins){

      # We get the names of the optimal lead designs in monopoles and dipoles
      bipolar_names_score_temp <- Find_optimal_monodipoles(metric_option = metric_option,
                                                           vein = each_vein, 
                                                           response = response,
                                                           which_cases = which_cases,
                                                           version = version, 
                                                           output = "both",
                                                           SA_folder = SA_folder,
                                                           flag_debugging = flag_debugging)
      
      bipolar_names_temp <- bipolar_names_score_temp[[1]]
      bipolar_scores_temp <- bipolar_names_score_temp[[2]]
      #We get the improvement for each design
      n_bs <- length(bipolar_scores_temp)
      bipolar_scores_temp <- c(bipolar_scores_temp[1],
                               bipolar_scores_temp[-1] - bipolar_scores_temp[-n_bs])
      
      bipolar_names_temp <- bipolar_names_temp[order(bipolar_scores_temp,
                                                     decreasing = TRUE)]
      bipolar_scores_temp <- sort(bipolar_scores_temp,decreasing = TRUE)
      
      bipolar_names <- c(bipolar_names,bipolar_names_temp)
      bipolar_scores <- c(bipolar_scores,bipolar_scores_temp)
    }
    
    # We have to create a unique data frame 
    
    df_wrapped <- data.frame(names = bipolar_names, scores = bipolar_scores)
    final_df <- data.frame(names = unique(bipolar_names))
    
    for(i in seq(1,nrow(final_df))){
      vector_of_patients <- df_wrapped[df_wrapped$names == final_df$names[i],2]
      # The idea is to use generic function given by the argument method.
      body <- paste0(method,"(v)")
      args <- "v"
      # The generic function is renamed as f(v)
      eval(parse(text = paste('f <- function(', args, ') { return(' ,
                              body , ')}', sep=''))) 
      
      final_df[i,2] <- f(vector_of_patients)
    }
    rm(bipolar_names,bipolar_scores)
    
    bipolar_names <- as.vector(final_df$names)
    bipolar_scores <- final_df[,2]
    
  }
  bipolar_names_ascii <- bin2ascii_lead(bipolar_names)
  # We create all the combinations of 4 electrodes from 1 to 8
  matrix_comb <- t(combn(1:8,4))
  
  vec_comb <- apply(matrix_comb,1,function(x) paste0(x,collapse = ''))
  
  quad_scores <- 0*c(1:length(vec_comb))
  matrix_quad_includes_bip <- matrix(data = 0,
                                     nrow = length(bipolar_names),
                                     ncol = length(vec_comb))
  
  # We create a vector with the scores of all the quadripoles
  for(i in seq(1:length(vec_comb))){
    for(j in seq(1:length(bipolar_names))){
      dip_in_quad <- ascii2bin_lead(vec_comb[i]) %>%
                      check_quadripole_includes_dipole(.,bipolar_names[j])
      if(dip_in_quad == 2){
        quad_scores[i] = quad_scores[i] + bipolar_scores[j]
        matrix_quad_includes_bip[j,i] <- 1
      }
    }
  }
  
  q <- 1
  possible_solutions <- c()
  

  # Indices to check all the quadpoles
  permutations_q_ciphers <- t(combn(nrow(matrix_comb),1)) 
  for(i in c(1:nrow(permutations_q_ciphers))){
    #All the bipoles will be included if merging all the binary vectors 
    #corresponding to that quadpole have at least a one in all their positions.
    flag_to_break <- matrix_quad_includes_bip[,permutations_q_ciphers[i]] 
    if(prod(flag_to_break) != 0){
      # We add a row to the possible solutions
      possible_solutions <- c(possible_solutions,
                              vec_comb[permutations_q_ciphers[i]])
    }
  }
  q <- q+1
  
  if(length(possible_solutions) != 0){
    possible_solutions <- as.data.frame(possible_solutions)
  }
  
  # print("Choosing the minimal amount of quadpoles needed...")
  # We search for the minimal amount of quadripolar designs that include all 
  # the bipolar designs.
  # We'll stop once we have a solution.
  while(length(possible_solutions) == 0){ 
    # Indices to check all the quadpoles.
    permutations_q_ciphers <- t(combn(nrow(matrix_comb),q)) 
    for(i in c(1:nrow(permutations_q_ciphers))){
      # All the bipoles will be included if merging all the binary vectors 
      # corresponding to that quadpole have at least a one in all their 
      # positions.
      flag_to_break <- matrix_quad_includes_bip[,permutations_q_ciphers[i,1]] 
      for (j in c(2:ncol(permutations_q_ciphers))){
        flag_to_break <- flag_to_break + 
                        matrix_quad_includes_bip[,permutations_q_ciphers[i,j]]
      }
      if(prod(flag_to_break) != 0){
        # We add a row to the possible solutions
        possible_solutions <- rbind(possible_solutions,
                                    vec_comb[permutations_q_ciphers[i,]])
      }
    }
    q <- q+1
  }
  
  
  # Among all the possible solutions we choose the one that maximise the
  #  benefits.
  
  score_possible_solutions <- matrix(nrow=nrow(possible_solutions),
                                     ncol = ncol(possible_solutions))
  
  # Could be optimised but meh
  
  for(i in c(1:nrow(score_possible_solutions))){
    for (j in c(1:ncol(score_possible_solutions))){
      score_possible_solutions[i,j] <- quad_scores[match(possible_solutions[i,j],
                                                         vec_comb)]
    }
  }
  
  cummulative_score <- rowSums(score_possible_solutions)

  # We choose the quadripole(s) with the maximum score
  max_idx <- which(cummulative_score == max(cummulative_score))
  if(q > 2){
    # In case of a draw between several multipolar designs we choose the one 
    # with the maximum benefit on the first one. If still in draw we continue.
    if(length(max_idx) > 1){
      possible_solutions <- possible_solutions[max_idx,]
      score_possible_solutions <- score_possible_solutions[max_idx,]
      for(i in c(1:nrow(possible_solutions))){
        order_scores <- order(score_possible_solutions[i,],decreasing = TRUE)
        possible_solutions[i,] <- possible_solutions[i,order_scores]
        score_possible_solutions[i,] <- score_possible_solutions[i,order_scores]
      }
      
      
      score_possible_solutions_unsorted <- score_possible_solutions
      # score_possible_solutions <- t(apply(score_possible_solutions, 1,
                                          # function(x) sort(x,decreasing = TRUE)))
      
      flag_to_break <- FALSE
      j <- 1
      
      while(!flag_to_break){
        max_idx <- which(score_possible_solutions[,j] ==
                           max(score_possible_solutions[,j]))
        
        if(length(max_idx) == 1){
          flag_to_break <- TRUE
        }
        else{
          j <- j+1
          possible_solutions <- possible_solutions[max_idx,]
          score_possible_solutions <- score_possible_solutions[max_idx,]
          if(j == ncol(score_possible_solutions)){
            flag_to_break = TRUE
            max_idx = 1
          }
        }
      }
      
    }
    else{
      score_possible_solutions_unsorted <- score_possible_solutions
    }
  }

  idx_sorted <- score_possible_solutions_unsorted[max_idx[1],] %>%
                order(., decreasing = TRUE)
  final_quad_solutions <- as.vector(possible_solutions[max_idx[1], idx_sorted])
  
  quad_scores <- score_possible_solutions[max_idx[1],]
  final_quad_solutions <- final_quad_solutions[order(quad_scores,
                                                     decreasing = TRUE)]
  quad_scores <- sort(quad_scores, decreasing = TRUE)
  
  return(list("Patient improvement"=bipolar_scores,
              "Bipolar"=bipolar_names_ascii,
              "Quadripolar"=final_quad_solutions,
              "Quadripolar_Scores"=quad_scores
                )
         )
}

