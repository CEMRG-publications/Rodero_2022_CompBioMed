#' @description Draws the table with the hierarchical agglomeratve clustering
#' of the leads and hearts of either cohort.
#' 
#' @param metric_option "TAT", "AT1090", "AT090", "LVTAT", "VEUTAT", "VEUmean"
#' for each of the metrics. In the latest version, we do not use TATLV anymore.
#' @param lead_option Vein to plot the lead configurations. Possibilities are 
#' "AN", "AL", "LA", "IL", "IN".
#' @param col_clusters Number of clusters for the columns (hearts). In the plot,
#' the different groups will be slightly sepparated.
#' @param row_clusters Analogous to col_clusters but with the rows (lead 
#' designs).
#' @param which_cases Condition of the patients. "RR" for reverse remodelled, 
#' "HF" for heart failure. Default is "RRHF".
#' @param max_bar Maximum value for the colorbar in the plot. Default is 30.
#' @param version Number of the version to read the file.
#' @param SA_folder Name of the folder of the simulations.
#' @param flag_return If TRUE it returns the object for the 
#' heatmap + dendrograms.
#' @param flag_debugging Flag (TRUE or FALSE) to print the name and location of
#' the files read and written.

Hier_clust_multipole <- function(metric_option, lead_option = "IN",
                                 col_clusters = 3, row_clusters = 10, 
                                 which_cases = "RRHF", max_bar = 30, version,
                                 SA_folder, flag_return = FALSE,
                                 flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dendsort","dendextend", "dplyr","seriation",
                          "stringr","BiocManager","circlize"))
  # BiocManager::install("ComplexHeatmap")
  library("ComplexHeatmap")
  
  
  if(metric_option == "AT1090"){
    plottitle <- "AT1090"
  } else if(metric_option == "TAT"){
    plottitle <- "TAT"
  } else if(metric_option == "LVTAT"){
    plottitle <- "TAT LV"
  } else if(metric_option == "VEUTAT"){
    plottitle <- "Total Ventr. Electrical Uncoupling"
  } else if(metric_option == "VEUmean"){
    plottitle <- "Mean Ventr. Electrical Uncoupling"
  } else if(metric_option == "AT090"){
    plottitle <- "AT090"
  }
  
  if(lead_option == "AN"){
    plottitle = paste0(plottitle," reduction pacing in the anterior position")
  } else if(lead_option == "AL"){
    plottitle = paste0(plottitle," reduction pacing in the antero-lateral position")
  } else if(lead_option == "LA"){
    plottitle = paste0(plottitle," reduction pacing in the lateral position")
  } else if(lead_option == "IL"){
    plottitle = paste0(plottitle," reduction pacing in the infero-lateral position")
  } else if(lead_option == "IN"){
    plottitle = paste0(plottitle," reduction pacing in the inferior position")
  }
  
  plottitle = paste0(plottitle,"\n")
  
  # We read the dataframe
  normalised_cohorts <- Read_dataframes(metric_option = metric_option, 
                                        lead_option = lead_option, 
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
  
  hf2color <- colnames(normalised_ATs)
  
  for(i in 1:length(hf2color)){
    if(substr(hf2color[i],2,2) == "F")
      hf2color[i] <- "HF"
    else
      hf2color[i] <- "RR"
  }
  
  
  # Complex heatmap
  
  row_dend <- set(as.dendrogram(dendsort(hclust(dist(normalised_ATs)))), 
                  "branches_lwd", 3)
  row_dend <- color_branches(row_dend, k = row_clusters)
  
  col_dend <- set(as.dendrogram(dendsort(hclust(dist(t(normalised_ATs))))),
                  "branches_lwd", 3)
  col_dend <- color_branches(col_dend, k = col_clusters)
  
  
  o1 = seriate(dist(normalised_ATs), method = "TSP")
  dend = dendsort(hclust(dist(normalised_ATs)))
  
  rownames(normalised_ATs) <- normalised_ATs %>% rownames() %>%
    bin2ascii_lead() %>%
    strsplit(.,",") %>%
    str_pad(.,side = "left",width = 4)
  
  hm=Heatmap(normalised_ATs,
             name = "RR+HF", # Internal name
             
             cluster_rows = row_dend,
             # cluster_rows = as.dendrogram(o1[[1]]),
             clustering_distance_rows = "euclidean",
             clustering_method_rows = "ward.D2",
             #row_dend_reorder = TRUE,
             #row_names_gp = gpar(col = ifelse(leadpos2color == "Medial",
             # "black", "orange")),
             row_split = row_clusters,
             # row_order = get_order(o1),
             # split=x,
             # right_annotation = rowAnnotation(
             #   Position = leadpos2color,
             #   # col = list(Position = c(Apical = "#FF6D00",
             #                             Basal = "#E65100",
             #                             ApicoBasal = "#FFFF00",
             #                             Medial = "#FFB74D")),
             #   col = list(Position = c(Apical = "#E1C62F",
             #                           Basal = "#FFFF00",
             #                           ApicoBasal = "#212121",
             #                           Medial = "#BCAAA4")),
             #   annotation_legend_param = list(
             #     Position = list(
             #       title = "Lead position",
             #       labels = c("Apical", "Basal", 
             #                  "Mixed apico-basal", "Others"),
             # 
             #       position = "topleft"
             #     )
             #   )
             # ),
             
             #cluster_columns = cluster_within_group(normalised_ATs,hf2color),
             cluster_columns = col_dend,
             clustering_distance_columns = "euclidean",
             clustering_method_columns = "ward.D2",
             #column_dend_reorder = TRUE,
             column_names_gp = gpar(col = ifelse(hf2color == "HF",
                                                 "blue", "red"),
                                    fontsize = 18, fontfamily = "Helvetica"),
             column_split = col_clusters,
             #column_order = get_order(o,2),
             # top_annotation = HeatmapAnnotation(
             #   LV_mass = anno_barplot(
             #     bar_width = 1,
             #     height = unit(40,"points"),
             #     c(RR_mass2color,HF_mass2color),
             #     gp = gpar(fill = colorRampPalette(brewer.pal(n = 9,
             #      name = "Purples"))(100))
             #     ),
             #     show_annotation_name = TRUE,
             #   LV_vol = anno_barplot(
             #     bar_width = 1,
             #     height = unit(40,"points"),
             #     c(RR_volume2color,HF_volume2color),
             #     gp = gpar(fill = colorRampPalette(brewer.pal(n = 9,
             #      name = "Greys"))(100))
             #   ),
             #   LV_quotient = anno_barplot(
             #     bar_width = 1,
             #     height = unit(40,"points"),
             #     c(RR_quotient,HF_quotient),
             #     gp = gpar(fill = colorRampPalette(brewer.pal(n = 9,
             #      name = "Purples"))(100))
             #   )
             # ),
             # 
             
             show_column_dend = TRUE,
             show_row_dend = TRUE,
             row_names_gp = gpar(fontsize = 15, fontfamily = "Helvetica"),
             
             col = colorRamp2(seq(0, max_bar, length = 2),
                              c("#EDF8E9", "#005A32")),
             na_col = "black",
             
             heatmap_legend_param = list(
               title = "AT reduction (%)\n",
               legend_height = unit(0.8, "npc"),
               grid_width = unit(0.25, "npc"),
               title_position = "topcenter",
               title_gp = gpar(fontsize = 20),
               labels_gp = gpar(fontsize = 15)
               # LV_mass = list(
               #   title = "LV mass (g)",
               #   direction = "horizonal"
               # )
             ),
             show_heatmap_legend = TRUE,
             
             column_title = plottitle,
             column_title_side = "top",
             column_title_gp = gpar(fontsize = 30, fontface = "bold",
                                    fontfamily = "Helvetica"),
             
             row_title = "\nLead design\n",
             row_title_side = "right",
             row_title_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
             
             width = unit(0.67, "npc"),
             height = unit(0.6, "npc"),
             row_dend_width = unit(0.1,"npc"),
             column_dend_height = unit(0.2, "npc")
             )
  
  # col_fun = colorRamp2(seq(0, 20, length = 2),c("#EDF8E9", "#005A32"))
  # # lgd = Legend(col_fun = col_fun, title = "AT reduction (%)")
  # # draw(lgd, x = unit(0.5, "npc"), y = unit(0.8, "npc"))
  # 
  # lgd = Legend(col_fun = col_fun, title = "foo", 
  # at = c(0, 0.25, 0.5, 0.75, 1))
  
  hm_return <- draw(hm)
  if(flag_return)
    return(hm_return)
  
  
}

#' @description Plots the optimal choice of quadripolar lead design based on
#' the optimal bipolar lead designs, covering all of them with the minimal 
#' amount of quadripoles possible. 
#' 
#' @param which_cases "RR", "HF" or "RRHF" for the corresponding cohort.
#' @param with_lines If FALSE, it only plots the histogram of the optimal
#' bipoles per vein. If TRUE, it also adds lines underneath showing which
#' bipoles are covered by each quadripole.
#' @param metric_option "TAT", "AT1090", "AT090", "VEUTAT", "VEUmean", "LVTAT",
#' for each of the metrics. 
#' @param bipole_from_file If TRUE, it reads the optimal bipolar _designs_ from
#' file, otherwise it computes them.
#' @param method Criterion to choose among the different vein. Default is "max",
#' so it will go always first to the one that maximises the number of patients.
#' @param savefile If TRUE it saves the plot in a file.
#' @param response If less than 100 it takes as "optimal" a wider window.
#' @param sort_lines Vector with the order of the lines at the bottom of the 
#' plot. Sometimes for some reason they are unsorted, so it's the way of sorting
#' them. 
#' @param SA_folder Name of the folder of the simulations.
#' @param flag_debugging Flag (TRUE or FALSE) to print the name and location of
#' the files read and written.

Plot_Dipoles_Fancy <- function(which_cases, with_lines = TRUE, metric_option,
                               bipole_from_file = FALSE, method = "max",
                               savefile = FALSE, response = 100, 
                               sort_lines = c(), SA_folder, 
                               show_percentages = FALSE,
                               flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("scales","dplyr","ggplot2"))
  
  if(bipole_from_file){
    # We read the cohort of bipoles
    optimal_bipole <- read.delim("/data/SA_multipole/FEC_70_values/optimal_bipole.txt")
    # We extract the rest corresponding to which_cases
    optimal_choices <- optimal_bipole[optimal_bipole$Cohort == which_cases,1:3]
    # Convert the bipoles from number to character
    optimal_choices$Bipole <- as.character(optimal_choices$Bipole)
  }
  else{
    optimal_choices <- as.data.frame(matrix(ncol=3,nrow = 0))
    colnames(optimal_choices) <- c("Vein","Patients","Bipole")
    # print("Finding optimal bipoles for each vein...")
    for (vein in c("AN","AL","LA","IL","IN")) {
      res_list <- Find_optimal_monodipoles(vein = vein, response = response,
                                           metric_option = metric_option,
                                           which_cases = which_cases, version=4,
                                           output = "both",
                                           SA_folder = SA_folder,
                                           flag_debugging = flag_debugging)
      
      temp_df <- as.data.frame(matrix(ncol=3,nrow = length(res_list[[1]])))
      colnames(temp_df) <- c("Vein","Patients","Bipole")
      
      temp_df$Vein <- vein
      
      # We get the increment of patients, not the total
      bipolar_scores <- as.vector(res_list[[2]])
      bipolar_scores <- c(bipolar_scores[1], 
                          bipolar_scores[-1] - 
                            bipolar_scores[-length(bipolar_scores)])
      
      temp_df$Patients <- bipolar_scores
      temp_df$Bipole <- bin2ascii_lead(as.vector(res_list[[1]]))
      
      optimal_choices <- rbind(optimal_choices,temp_df)
      rm(res_list)
    }
  }
  if(with_lines){
    # res_list <- find_optimal_quadripole("ALL",which_cases,method)
    #res_list <- find_optimal_quadripole_combinatorial("ALL",which_cases,method)
    #res_list <- Find_Optimal_Quadripole_Evenly("ALL",which_cases,method)
    # print("Finding optimal quadripoles for all the veins at once...")
    res_list <- Find_Optimal_Quadripole_Optimising(vein = "ALL",
                                                   which_cases = which_cases,
                                                   metric_option = metric_option,
                                                   method = method, 
                                                   response = response,
                                                   SA_folder = SA_folder,
                                                   version = 4, 
                                                   flag_debugging = flag_debugging)
    # print(res_list)
    
    optimal_quadpoles <- res_list$Quadripolar
  }
  # print("Setting up the plot...")
  
  # We sort the bipoles (alphabetically)
  all_bipoles <-  sort(unique(optimal_choices$Bipole))
  
  # We create fake bars for when there's no bipole, setting a height of 0.1
  for(vein_name in c("AN","AL","LA","IL","IN")){
    # Extract the data frame corresponding to each vein
    vein_vec <- optimal_choices[optimal_choices$Vein == vein_name,]
    # Extract the bipoles that are not included in this subset
    which_to_add <- all_bipoles[!all_bipoles %in% vein_vec$Bipole]
    # Create dummy cases with height 0.1
    df_to_add <- data.frame(rep(vein_name,length(which_to_add)),
                            rep(0.1,length(which_to_add)),which_to_add)
    names(df_to_add) <- names(optimal_choices)
    optimal_choices <- rbind(optimal_choices,df_to_add)
    rm(list = c("vein_vec","which_to_add","df_to_add"))
  }
  
  # For the title of the plot and number of patients
  if(which_cases == "RRHF"){
    title_str <- "Both cohorts"
    num_patients <- 44
  }
  else if(which_cases == "RR"){
    title_str <- paste0(which_cases," cohort")
    num_patients <- 20
  }
  else if(which_cases == "HF"){
    title_str <- paste0(which_cases," cohort")
    num_patients <- 24
  }
  
  # title_str <- paste0(title_str," grouped by ",method, "\n in the ", 
  #                     response, "% window")
  title_str <- paste0("MPP lead choice in the ",title_str)
  
  # For the lines coordinates
  if(with_lines){
    first_x <- 0.55
    gap_x <- 1
    segment_length <- 0.9
    first_y <- -0.5
    # if(length(optimal_quadpoles) <= 2)
    #   gap_y <- -2
    # else
    gap_y <- -2
    
    x_ini <- c()
    x_end <- c()
    which_quadpole <- c()
    
    # We create the coordinates of each line if the quadpole includes the bipole  
    for(i in seq(1,length(all_bipoles))){
      for(j in seq(1,length(optimal_quadpoles))){
        if(isbipoleinquadpole(all_bipoles[i],optimal_quadpoles[j])){
          x_ini[length(x_ini) + 1] <- first_x + gap_x*(i-1)
          x_end[length(x_end) + 1] <- x_ini[length(x_ini)] + segment_length
          which_quadpole[length(which_quadpole) + 1] <- optimal_quadpoles[j]
        }
      }
    }
    
    # if(!quad_from_files){
    sorted_patients <- res_list$`Patient improvement`[order(res_list$Bipolar)]
    all_bipoles <- sort(res_list$Bipolar)
    # }
    
    # We compute the percentage of patients for which each quadpole is optimal
    perc_bipoles <- c()
    for(i in seq(1:length(optimal_quadpoles))){
      temp_perc <- 0
      for(j in seq(1:length(all_bipoles))){
        if(isbipoleinquadpole(all_bipoles[j],optimal_quadpoles[i]) == 1){
          temp_perc <- temp_perc + sorted_patients[j]
        }
      }
      perc_bipoles[i] <- round(100*(sum(temp_perc)/(num_patients)),2)
    }
    
    perc_quadpoles <- res_list$Quadripolar_Scores %>% 
                      as.vector(.)/(0.01*num_patients) %>%
                      round(.,2)
    perc_bipoles <- round(perc_quadpoles/5,2)
    # We create and auxiliar vector for the text of the legend
    if(show_percentages){
      legend_text <- paste0(optimal_quadpoles," (",perc_bipoles,"%)")
    }
    else{
      legend_text <- optimal_quadpoles
    }
    
    # The y coordinates from up to down will be
    y_coordinates <- first_y + gap_y*(seq(0,length(optimal_quadpoles)))
    if(length(sort_lines) > 0){
      y_coordinates <- y_coordinates[sort_lines]
    }
    # We need the indices of the sorting, from higher to lower
    indices <- rev(order(perc_bipoles))
    
    # The idea is to put the higher percentage lines on a higher position
    y <- c()
    legend_vector <- c()
    
    
    for(i in seq(1,length(which_quadpole))){
      # We check which optimal quadpole is it
      for(j in seq(1,length(optimal_quadpoles))){
        if(which_quadpole[i] == optimal_quadpoles[j]){
          # If the highest value is the 5th quadpole, we want the index 1,
          #  the indices vector will be c(5,...)
          if(length(sort_lines) == 0){
            pos_in_y <- match(j, indices)
          }
          else{
            pos_in_y <- j
          }
          y[i] <- y_coordinates[pos_in_y]
          break
        }
      }
    }
    
    # We create the data frame
    only_quadripole <- which_quadpole
    line_df <- data.frame(x_ini,x_end,y,only_quadripole)
    
    # We need the auxiliar variable for the lines I think
    line_df$row <- seq_len(nrow(line_df))
    
    
    # We add a column for the legend
    for(i in seq(1,nrow(line_df))){
      corresponding_quad <- as.character(line_df$only_quadripole[i])
      pos_in_perc <- match(corresponding_quad,optimal_quadpoles)
      corresponding_perc <- perc_bipoles[pos_in_perc]
      if(show_percentages){
        line_df$Quadripole[i] <- paste0(corresponding_quad," (",
                                        corresponding_perc,"%)")
      }
      else{
        line_df$Quadripole[i] <- corresponding_quad
      }
    }
    
    
    line_df$Quadripole <- factor(line_df$Quadripole,
                                 levels = legend_text[indices])
  }
  
  # colnames(line_df)[ncol(line_df)] <- "Design (% patients for whom it is optimal)"
  
  # We sort the categories so they appear in order in the histogram
  optimal_choices$Vein <- factor(optimal_choices$Vein,
                                 levels = c("AN","AL","LA","IL","IN"))
  
  toplot = ggplot(optimal_choices) + 
    scale_fill_manual(values=c("#FC0000", "#ED7c31", "#FFBF00", "#00B04F",
                               "#00B0F0")) +
    geom_bar(aes(Bipole, Patients, fill=Vein), stat="identity",position='dodge') 
  
  if(with_lines){
    toplot = toplot + geom_segment(data = line_df, aes(x = x_ini, xend = x_end,
                                                       y = y, yend = y,
                                                       group = row, 
                                                       color = Quadripole),
                                   size = 2) +
      # scale_color_manual(values = brewer.pal(length(optimal_quadpoles),"Paired"))
      scale_color_manual(values = c("green","orange","red"))
      
    if(show_percentages){
      toplot <- toplot + labs(color = "Design (% of subjects)")
    }
    else{
      toplot <- toplot + labs(color = "Design")
    }
  }
  # print(num_patients)
  toplot = toplot + theme_classic(base_size = 40) + 
    theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank())  + 
    scale_y_continuous(breaks = seq(0,num_patients,5), 
                       limits = c(min(line_df$y),num_patients),
                       oob = rescale_none) +
    ggtitle(title_str) + xlab("Electrodes activated")
  
  
  print(toplot)
  
  if(savefile)
    ggsave(filename=paste0("/home/crg17/Pictures/",which_cases,"_response_",
                           toString(response),"_",SA_folder,".png"),
           width = 20,height = 10)
  
}

Show_Delta_Optimal_Vein_Quadripole <- function(cohort = "HF", bar_limit = 15, flag_debugging = FALSE){
  
  # Find the optimal quadpoles  
  optimal_quadpoles_global <- Find_Optimal_Quadripole_Optimising(vein = "ALL",
                                                                 which_cases = cohort,method = "max",response=100,SA_folder = "", # deprecated, it was the default one
                                                                 version = 4)
  optimal_quadpoles_global <- optimal_quadpoles_global$Quadripolar # ASCII
  #From those, we keep the best one
  optimal_quadpole <- optimal_quadpoles_global[1]
  cyphers <- strsplit(optimal_quadpole,"")[[1]]
  
  #We extract all the bipoles from the optimal quadpole
  combs <- combn(4,2)
  
  vec_bipoles_from_quadpoles <- c()
  
  for(i in c(1:ncol(combs))){
    vec_bipoles_from_quadpoles[length(vec_bipoles_from_quadpoles) + 1] <- paste0(cyphers[combs[1,i]],cyphers[combs[2,i]])
  }
  
  
  
  if((cohort == "RR") || (cohort == "RRHF")){
    
    delta_df_RR <- c()
    
    for(vein in c("AN","AL","LA","PL","PO")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = "AT1090",
                                            lead_option = vein, 
                                            which_cases = cohort,
                                            version = 4,
                                            SA_folder = "", # It was the default one 
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_RR
      
      
      # Get only bipoles
      triquad_idx <- which(sum_string(ascii2bin_lead(bin2ascii_lead(rownames(normalised_HF)))) != 2)
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- bin2ascii_lead(apply(normalised_ATs, 2, function(x) rownames(normalised_ATs)[as.vector(which.is.max(x))]))
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      
      # If it's in the optimal quadripole that's it, otherwise we take the best that we can
      value_if_not_included <- apply(normalised_ATs[ascii2bin_lead(vec_bipoles_from_quadpoles),],2,max)
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- max_values_bipoles[i] - value_if_not_included[i]
        }
      }
      
      delta_df_RR <- rbind(delta_df_RR,delta_vein)
      
    }
    
    delta_df_RR <- as.data.frame(delta_df_RR)
    colnames(delta_df_RR) <- names(max_values_bipoles)
    rownames(delta_df_RR) <- c("AN","AL","LA","PL","PO")
    
  }
  
  
  if((cohort == "HF") || (cohort == "RRHF")){
    
    delta_df_HF <- c()
    
    for(vein in c("AN","AL","LA","PL","PO")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = "AT1090",
                                            lead_option = vein,
                                            which_cases = cohort, version = 4,
                                            SA_folder = "", # Deprecated, it was the default one
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_HF
      
      
      # Get only bipoles
      triquad_idx <- which(sum_string(ascii2bin_lead(bin2ascii_lead(rownames(normalised_HF)))) != 2)
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- bin2ascii_lead(apply(normalised_ATs, 2, function(x) rownames(normalised_ATs)[as.vector(which.is.max(x))]))
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      
      # If it's in the optimal quadripole that's it, otherwise we take the best that we can
      value_if_not_included <- apply(normalised_ATs[ascii2bin_lead(vec_bipoles_from_quadpoles),],2,max)
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- max_values_bipoles[i] - value_if_not_included[i]
        }
      }
      
      delta_df_HF <- rbind(delta_df_HF,delta_vein)
      
    }
    
    delta_df_HF <- as.data.frame(delta_df_HF)
    colnames(delta_df_HF) <- names(max_values_bipoles)
    rownames(delta_df_HF) <- c("AN","AL","LA","PL","PO")
    
  }
  
  if(cohort == "RR")
    delta_df <- delta_df_RR
  else if(cohort == "HF")
    delta_df <- delta_df_HF
  else if(cohort == "RRHF")
    delta_df <- cbind(delta_df_RR,delta_df_HF)
  
  
  width_npc <- 0.88
  height_npc <- 0.17
  
  width_npc <- width_npc*ncol(delta_df)/44
  #height_npc <- height_npc*ncol(delta_df)/44
  
  
  hm=Heatmap(as.matrix(delta_df),
             rect_gp = gpar(col = 'black', lty = 'dashed'),
             row_names_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
             row_order = rownames(delta_df),
             column_order = colnames(delta_df),
             column_names_gp = gpar(fontsize = 20),
             column_names_rot = 60,
             
             col = colorRamp2(seq(0, bar_limit, length = 2),c("green", "red")),
             na_col = "black",
             
             heatmap_legend_param = list(
               title = "Improvement \nreduction (%)\n",
               legend_height = unit(0.8, "npc"),
               grid_width = unit(0.25, "npc"),
               title_position = "topcenter",
               title_gp = gpar(fontsize = 20),
               labels_gp = gpar(fontsize = 20)
             ),
             show_heatmap_legend = TRUE,
             
             column_title = "Difference of AT improvement between the optimal\n quadripolar design for the cohort and optimal design for each case\n",
             column_title_side = "top",
             column_title_gp = gpar(fontsize = 30, fontface = "bold", fontfamily = "Helvetica"),
             
             row_title = "\nVein position\n",
             row_title_side = "right",
             row_title_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
             
             width = unit(width_npc, "npc"),
             height = unit(height_npc, "npc")
  )
  
  draw(hm)
  return(delta_df)
}

#' @description Function to plot a table where each tile is colorcoded according
#' to the improvement in % of choosing the personalised design versus choosing
#' the global design (per vein). For example, if with the global optimal you 
#' achieve an improvement in AT090 reduction of 20% and with the personalised 
#' you achieve a 21% the output will be 1%.
#' 
#' @param cohort "RR", "HF" or "RRHF" for the corresponding cohort.
#' @param metric_option "TAT", "AT1090", "AT090", "VEUTAT", "VEUmean", "LVTAT",
#' for each of the metrics. 
#' @param SA_folder Name of the folder of the simulations.
#' @param bar_limit Upper bound of the colour scale.
#' @param optimal_idx If 1, we get the configuration that its optimal for the
#' maximum number of patients; if 2, it's the second configuration and so on.
#' @param output The dataframe it plots "df" for the improvement of personalised
#' vs global, "global" to plot the improvement of the optimal global design over
#' the baseline and "personalised" for the analogous plot with the 
#' patient-specific design.
#' @param flag_debugging Flag (TRUE or FALSE) to print the name and location of
#' the files read and written.
#' 
#' @return The dataframe with the values of whatever it has plotted.

Show_Delta_Optimal_Cohort_Quadripole <- function(cohort = "HF", metric_option, 
                                                 SA_folder, bar_limit = 15,
                                                 optimal_idx = 1, output = "df",
                                                 flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr","nnet"))
  
  optimal_reductions_df <- c() 
  optimal_bipoles_df <- c()
  global_reductions_df <- c()
  
  # Find the optimal quadpoles  

  optimal_quadpoles_global <-
    Find_Optimal_Quadripole_Optimising(vein = "ALL", which_cases = cohort,
                                       method = "max", response=100, 
                                       SA_folder = SA_folder, version = 4,
                                       metric_option = metric_option,
                                       flag_debugging = flag_debugging)
  
  optimal_quadpoles_global <- optimal_quadpoles_global$Quadripolar # ASCII
  #From those, we keep the best one
  optimal_quadpole <- optimal_quadpoles_global[optimal_idx]
  print(optimal_quadpole)
  cyphers <- strsplit(optimal_quadpole,"")[[1]]
  
  #We extract all the bipoles from the optimal quadpole
  combs <- combn(4,2)
  
  vec_bipoles_from_quadpoles <- c()
  
  for(i in c(1:ncol(combs))){
    vec_bipoles_from_quadpoles[length(vec_bipoles_from_quadpoles) + 1] <-
      paste0(cyphers[combs[1,i]],cyphers[combs[2,i]])
  }
  
  if((cohort == "RR") || (cohort == "RRHF")){
    
    delta_df_RR <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein, 
                                            which_cases = cohort, version = 4, 
                                            SA_folder = SA_folder, 
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_RR
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
                  bin2ascii_lead(.) %>%
                  ascii2bin_lead(.)
      triquad_idx <- which(sum_string(num_elec) != 2)
      
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
            apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
            bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- max_values_bipoles[i] - value_if_not_included[i]
        }
      }
      
      delta_df_RR <- rbind(delta_df_RR,delta_vein)
      
    }
    
    delta_df_RR <- as.data.frame(delta_df_RR)
    colnames(delta_df_RR) <- names(max_values_bipoles)
    rownames(delta_df_RR) <- c("AN","AL","LA","IL","IN")
    
  }
  
  
  if((cohort == "HF") || (cohort == "RRHF")){
    
    delta_df_HF <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein,
                                            which_cases = cohort, version = 4,
                                            SA_folder = SA_folder,
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_HF
      
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
                  bin2ascii_lead(.) %>%
                  ascii2bin_lead(.)
                  
      triquad_idx <- which(sum_string(num_elec) != 2)
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
            apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
            bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best 
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- max_values_bipoles[i] - value_if_not_included[i]
        }
      }
      
      delta_df_HF <- rbind(delta_df_HF,delta_vein)
      
    }
    
    delta_df_HF <- as.data.frame(delta_df_HF)
    colnames(delta_df_HF) <- names(max_values_bipoles)
    rownames(delta_df_HF) <- c("AN","AL","LA","IL","IN")
    
  }
  
  if(cohort == "RR")
    delta_df <- delta_df_RR
  else if(cohort == "HF")
    delta_df <- delta_df_HF
  else if(cohort == "RRHF")
    delta_df <- cbind(delta_df_RR,delta_df_HF)
  
  if(output == "df"){
    
    width_npc <- 1.05
    height_npc <- 0.29
    
    width_npc <- width_npc*ncol(delta_df)/44
    #height_npc <- height_npc*ncol(delta_df)/44
    
    title_str <- paste0("Improvement of choosing personalised\n",
                        "vs global optimal design")
    
    hm=Heatmap(as.matrix(delta_df),
               # cell_fun = function(j, i, x, y, w, h, col) { 
               #   # add text to each grid
               #   grid.text(round(delta_df[i, j],2), x, y)
               # },
               rect_gp = gpar(col = 'black', lty = 'dashed'),
               row_names_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               row_order = rownames(delta_df),
               column_order = colnames(delta_df),
               column_names_gp = gpar(fontsize = 20),
               column_names_rot = 60,
               
               col = colorRamp2(seq(0, bar_limit, length = 2),
                                c("green", "red")),
               na_col = "black",
               
               heatmap_legend_param = list(
                 title = "Improvement (%)\n",
                 legend_height = unit(0.8, "npc"),
                 grid_width = unit(0.25, "npc"),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = 20),
                 labels_gp = gpar(fontsize = 20)
               ),
               show_heatmap_legend = TRUE,
               
               column_title = title_str,
               column_title_side = "top",
               column_title_gp = gpar(fontsize = 30,
                                      fontface = "bold",
                                      fontfamily = "Helvetica"),
               
               row_title = "\nVein position\n",
               row_title_side = "right",
               row_title_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               
               width = unit(width_npc, "npc"),
               height = unit(height_npc, "npc")
    )
    
    draw(hm)
  }
  else if(output == "global"){
    
    width_npc <- 0.88
    height_npc <- 0.17
    
    width_npc <- width_npc*ncol(global_reductions_df)/44
    #height_npc <- height_npc*ncol(delta_df)/44
    
    title_str <- paste0("Improvement of choosing global\n",
                        "over baseline")
    
    hm=Heatmap(as.matrix(global_reductions_df),
               rect_gp = gpar(col = 'black', lty = 'dashed'),
               row_names_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               row_order = rownames(global_reductions_df),
               column_order = colnames(global_reductions_df),
               column_names_gp = gpar(fontsize = 20),
               column_names_rot = 60,
               
               col = colorRamp2(seq(0, bar_limit, length = 2),
                                c("red", "green")),
               na_col = "black",
               
               heatmap_legend_param = list(
                 title = "Improvement \nreduction (%)\n",
                 legend_height = unit(0.8, "npc"),
                 grid_width = unit(0.25, "npc"),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = 20),
                 labels_gp = gpar(fontsize = 20)
               ),
               show_heatmap_legend = TRUE,
               
               column_title = title_str,
               column_title_side = "top",
               column_title_gp = gpar(fontsize = 30,
                                      fontface = "bold",
                                      fontfamily = "Helvetica"),
               
               row_title = "\nVein position\n",
               row_title_side = "right",
               row_title_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               
               width = unit(width_npc, "npc"),
               height = unit(height_npc, "npc")
    )
    
    draw(hm)
    
    delta_df <- global_reductions_df
  }
  else if(output == "personalised"){
    
    width_npc <- 0.88
    height_npc <- 0.17
    
    width_npc <- width_npc*ncol(optimal_reductions_df)/44
    #height_npc <- height_npc*ncol(delta_df)/44
    
    title_str <- paste0("Improvement of choosing personalised\n",
                        "over baseline")
    
    hm=Heatmap(as.matrix(optimal_reductions_df),
               rect_gp = gpar(col = 'black', lty = 'dashed'),
               row_names_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               row_order = rownames(optimal_reductions_df),
               column_order = colnames(optimal_reductions_df),
               column_names_gp = gpar(fontsize = 20),
               column_names_rot = 60,
               
               col = colorRamp2(seq(0, bar_limit, length = 2),
                                c("red", "green")),
               na_col = "black",
               
               heatmap_legend_param = list(
                 title = "Improvement \nreduction (%)\n",
                 legend_height = unit(0.8, "npc"),
                 grid_width = unit(0.25, "npc"),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = 20),
                 labels_gp = gpar(fontsize = 20)
               ),
               show_heatmap_legend = TRUE,
               
               column_title = title_str,
               column_title_side = "top",
               column_title_gp = gpar(fontsize = 30,
                                      fontface = "bold",
                                      fontfamily = "Helvetica"),
               
               row_title = "\nVein position\n",
               row_title_side = "right",
               row_title_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               
               width = unit(width_npc, "npc"),
               height = unit(height_npc, "npc")
    )
    
    draw(hm)
    
    delta_df <- optimal_reductions_df
  }
  
  return(delta_df)
}


#' @description Function to plot with a single instruction multiple plots.
#' 
#' @param plottype Type of plot to print. For now only input is "HAC", for the
#' hierarchical agglomerative clustering. 
#' @param metric_option "TAT", "AT1090", "AT090" or "LVTAT" for each of the
#'  metrics. "ALL" to use all of them sequentially.
#' @param lead_option "AN", "AL", "LA", "IL", "IN" for each of the lead 
#' positions. "ALL" to use all of them sequentially.
#' @param which_cases "RR", "HF" or "RRHF" for each of the cohorts. "ALL" to use
#' all of them sequentially.
#' @param SA_folder Name of the folder of the simulations.
#' @param flag_debugging If TRUE, it prints whatever is reading and printing.
#' 
#' 
Plot_all <- function(plottype = "HAC", metric_option = "TAT",
                     lead_option = "AN", which_cases = "HF",
                     SA_folder = "default_HF_noPVTV", flag_debugging = FALSE){
  
  if(metric_option == "ALL"){
    for(m_option in c("TAT","LVTAT","AT090","AT1090")){
    Plot_all(plottype = plottype, metric_option = m_option,
             lead_option = lead_option, which_cases = which_cases, 
             SA_folder = SA_folder, flag_debugging = flag_debugging)
    }
  }
  
  if(lead_option == "ALL"){
    for(l_option in c("AN","AL","LA","IL","IN")){
      Plot_all(plottype = plottype, metric_option = metric_option,
               lead_option = l_option, which_cases = which_cases, 
               SA_folder = SA_folder, flag_debugging = flag_debugging)
    }
  }
  
  if(which_cases == "ALL"){
    for(wc_option in c("HF","RR","RRHF")){
      Plot_all(plottype = plottype, metric_option = metric_option,
               lead_option = lead_option, which_cases = wc_option, 
               SA_folder = SA_folder, flag_debugging = flag_debugging)
    }
  }
  
  if(metric_option != "ALL" && lead_option != "ALL" && which_cases != "ALL"){
  if(plottype == "HAC"){
    hm <- Hier_clust_multipole(metric_option = metric_option,
                         lead_option = lead_option, col_clusters = 4,
                         row_clusters = 4, which_cases = which_cases,
                         max_bar = 30, version = 4, SA_folder = SA_folder,
                         flag_return = TRUE,
                         flag_debugging = flag_debugging)
    
    paste0("/home/crg17/Pictures/default_HF_noPVTV/",plottype,"_",which_cases,
           "_",metric_option,"_",lead_option,".png") %>%
      png(.,width = 1920, height = 1003,
          units = "px", pointsize = 12, bg = "white", res = 120)
      draw(hm)
  }
  
  
  dev.off()
  }
  
}

Plot_lack_improvement <- function(metric_option = "AT090",
                                  SA_folder = "default_HF_noPVTV"){
  
  Load_Install_Packages(c("dplyr","corrplot"))
  
  cohorts <- c("HF", "RR", "RRHF")
  df_plot <- as.data.frame(matrix(nrow = 3, ncol = 5))
  rownames(df_plot) <- cohorts 
  
  for(i in c(1:length(cohorts))){
  delta <- Show_Delta_Optimal_Cohort_Quadripole(cohort = cohorts[i],
                                                   metric_option = metric_option,
                                                   SA_folder = SA_folder,
                                                   bar_limit = 1)
  
  df_plot[i,] <- delta %>% apply(., 1, function(x) sum(x < 1))
  df_plot[i,] <- 100* df_plot[i,] / ncol(delta)
  }
  
  colnames(df_plot) <- rownames(delta)
  df_plot <- round(df_plot,2)
  
  png(filename = paste0("/home/crg17/Pictures/pies_no_improvement.png"),
      width = 1800, height = 1000, units = "px", pointsize = 12, bg = "white",
      res = 100)
  cex_value = 3
  corrplot(as.matrix(df_plot),
           "pie",
           is.corr = FALSE,
           cl.lim = c(0,100),
           mar = c(0,0,8,12), #bottom, left, top, right
           tl.cex = cex_value,
           cl.cex = cex_value,
           tl.col = c("#FF0000","#ED7C31","#FFBF00","#00B04F","#00B0F0"),
           tl.srt = 60,
           col = colorRampPalette(c("white","#EDF8E9","#005A32"))(20)
           )
  colnames(df_plot) <- rep("",5)
  
  corrplot(as.matrix(df_plot),
           "pie",
           is.corr = FALSE,
           cl.lim = c(0,100),
           mar = c(0,0,8,12),
           tl.cex = cex_value,
           cl.cex = cex_value,
           tl.col = "black",
           tl.srt = 60,
           col = colorRampPalette(c("white","#EDF8E9","#005A32"))(20),
           add = TRUE
  )
  
  mtext(paste0("Percentage of patients who do not achieve an improvement\n",
               " greater than 1% if the lead design is personalised\n"),
        side = 3,
        cex = cex_value + 1,
        line = -7,
        at = 3)
  dev.off()
}


#' @description Function to calculate the maximum loss (in % of AT) of not
#' choosing the optimal for each vein and choosing the optimal from the quad
#' design.
#' 
ShowDeltaDesignVSCohort <- function(cohort = "RR", metric_option, 
                                    output = "%", SA_folder, bar_limit = 15,
                                    quad_design = "1267", draw_option = "df",
                                    flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr","nnet","ComplexHeatmap","circlize"))
  
  optimal_reductions_df <- c() 
  optimal_bipoles_df <- c()
  global_reductions_df <- c()
  
  optimal_quadpole <- quad_design
  cyphers <- strsplit(optimal_quadpole,"")[[1]]
  
  #We extract all the bipoles from the optimal quadpole
  combs <- combn(4,2)
  
  vec_bipoles_from_quadpoles <- c()
  
  for(i in c(1:ncol(combs))){
    vec_bipoles_from_quadpoles[length(vec_bipoles_from_quadpoles) + 1] <-
      paste0(cyphers[combs[1,i]],cyphers[combs[2,i]])
  }
  
  if((cohort == "RR") || (cohort == "RRHF")){
    
    delta_df_RR <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein, 
                                            which_cases = cohort, version = 4,
                                            output = output,
                                            SA_folder = SA_folder, 
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_RR
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
        bin2ascii_lead(.) %>%
        ascii2bin_lead(.)
      triquad_idx <- which(sum_string(num_elec) != 2)
      
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
        apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
        bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- max_values_bipoles[i] - value_if_not_included[i]
        }
      }
      
      delta_df_RR <- rbind(delta_df_RR,delta_vein)
      
    }
    
    delta_df_RR <- as.data.frame(delta_df_RR)
    colnames(delta_df_RR) <- names(max_values_bipoles)
    rownames(delta_df_RR) <- c("AN","AL","LA","IL","IN")
    
  }
  
  
  if((cohort == "HF") || (cohort == "RRHF")){
    
    delta_df_HF <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein, output = output,
                                            which_cases = cohort, version = 4,
                                            SA_folder = SA_folder,
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_HF
      
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
        bin2ascii_lead(.) %>%
        ascii2bin_lead(.)
      
      triquad_idx <- which(sum_string(num_elec) != 2)
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
        apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
        bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best 
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- max_values_bipoles[i] - value_if_not_included[i]
        }
      }
      
      delta_df_HF <- rbind(delta_df_HF,delta_vein)
      
    }
    
    delta_df_HF <- as.data.frame(delta_df_HF)
    colnames(delta_df_HF) <- names(max_values_bipoles)
    rownames(delta_df_HF) <- c("AN","AL","LA","IL","IN")
    
  }
  
  if(cohort == "RR")
    delta_df <- delta_df_RR
  else if(cohort == "HF")
    delta_df <- delta_df_HF
  else if(cohort == "RRHF")
    delta_df <- cbind(delta_df_RR,delta_df_HF)
  
  if(draw_option == "df"){
    
    width_npc <- 1.05
    height_npc <- 0.29
    
    width_npc <- width_npc*ncol(delta_df)/44
    #height_npc <- height_npc*ncol(delta_df)/44
    
    title_str <- paste0("Improvement of choosing personalised\n",
                        "vs global optimal design")
    
    hm=Heatmap(as.matrix(delta_df),
               # cell_fun = function(j, i, x, y, w, h, col) { 
               #   # add text to each grid
               #   grid.text(round(delta_df[i, j],2), x, y)
               # },
               rect_gp = gpar(col = 'black', lty = 'dashed'),
               row_names_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               row_order = rownames(delta_df),
               column_order = colnames(delta_df),
               column_names_gp = gpar(fontsize = 20),
               column_names_rot = 60,
               
               col = colorRamp2(seq(0, bar_limit, length = 2),
                                c("green", "red")),
               na_col = "black",
               
               heatmap_legend_param = list(
                 title = "Improvement (%)\n",
                 legend_height = unit(0.8, "npc"),
                 grid_width = unit(0.25, "npc"),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = 20),
                 labels_gp = gpar(fontsize = 20)
               ),
               show_heatmap_legend = TRUE,
               
               column_title = title_str,
               column_title_side = "top",
               column_title_gp = gpar(fontsize = 30,
                                      fontface = "bold",
                                      fontfamily = "Helvetica"),
               
               row_title = "\nVein position\n",
               row_title_side = "right",
               row_title_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               
               width = unit(width_npc, "npc"),
               height = unit(height_npc, "npc")
    )
    
    draw(hm)
  }
  else if(draw_option == "global"){
    
    width_npc <- 0.88
    height_npc <- 0.17
    
    width_npc <- width_npc*ncol(global_reductions_df)/44
    #height_npc <- height_npc*ncol(delta_df)/44
    
    title_str <- paste0("Improvement of choosing global\n",
                        "over baseline")
    
    hm=Heatmap(as.matrix(global_reductions_df),
               rect_gp = gpar(col = 'black', lty = 'dashed'),
               row_names_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               row_order = rownames(global_reductions_df),
               column_order = colnames(global_reductions_df),
               column_names_gp = gpar(fontsize = 20),
               column_names_rot = 60,
               
               col = colorRamp2(seq(0, bar_limit, length = 2),
                                c("red", "green")),
               na_col = "black",
               
               heatmap_legend_param = list(
                 title = "Improvement \nreduction (%)\n",
                 legend_height = unit(0.8, "npc"),
                 grid_width = unit(0.25, "npc"),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = 20),
                 labels_gp = gpar(fontsize = 20)
               ),
               show_heatmap_legend = TRUE,
               
               column_title = title_str,
               column_title_side = "top",
               column_title_gp = gpar(fontsize = 30,
                                      fontface = "bold",
                                      fontfamily = "Helvetica"),
               
               row_title = "\nVein position\n",
               row_title_side = "right",
               row_title_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               
               width = unit(width_npc, "npc"),
               height = unit(height_npc, "npc")
    )
    
    draw(hm)
    
    delta_df <- global_reductions_df
  }
  else if(draw_option == "personalised"){
    
    width_npc <- 0.88
    height_npc <- 0.17
    
    width_npc <- width_npc*ncol(optimal_reductions_df)/44
    #height_npc <- height_npc*ncol(delta_df)/44
    
    title_str <- paste0("Improvement of choosing personalised\n",
                        "over baseline")
    
    hm=Heatmap(as.matrix(optimal_reductions_df),
               rect_gp = gpar(col = 'black', lty = 'dashed'),
               row_names_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               row_order = rownames(optimal_reductions_df),
               column_order = colnames(optimal_reductions_df),
               column_names_gp = gpar(fontsize = 20),
               column_names_rot = 60,
               
               col = colorRamp2(seq(0, bar_limit, length = 2),
                                c("red", "green")),
               na_col = "black",
               
               heatmap_legend_param = list(
                 title = "Improvement \nreduction (%)\n",
                 legend_height = unit(0.8, "npc"),
                 grid_width = unit(0.25, "npc"),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = 20),
                 labels_gp = gpar(fontsize = 20)
               ),
               show_heatmap_legend = TRUE,
               
               column_title = title_str,
               column_title_side = "top",
               column_title_gp = gpar(fontsize = 30,
                                      fontface = "bold",
                                      fontfamily = "Helvetica"),
               
               row_title = "\nVein position\n",
               row_title_side = "right",
               row_title_gp = gpar(fontsize = 20, fontfamily = "Helvetica"),
               
               width = unit(width_npc, "npc"),
               height = unit(height_npc, "npc")
    )
    
    draw(hm)
    
    delta_df <- optimal_reductions_df
  }
  
  return(delta_df)
}

#' @description Function to calculate the maximum RELATIVE loss (in % of AT) of
#' not choosing the optimal for each vein and choosing the optimal from the quad
#' design. For example, if the optimal reduces 25% and the design 24%, the 
#' absolute delta would be 1%, while the relative would be 4%.
#' 
  RelativeDeltaDesignVSCohort <- function(cohort = "RR", metric_option, 
                                    SA_folder, 
                                    quad_design = "1267", 
                                    flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr","nnet","ComplexHeatmap","circlize"))
  
  optimal_reductions_df <- c() 
  optimal_bipoles_df <- c()
  global_reductions_df <- c()
  
  
  optimal_quadpole <- quad_design
  cyphers <- strsplit(optimal_quadpole,"")[[1]]
  
  #We extract all the bipoles from the optimal quadpole
  combs <- combn(4,2)
  
  vec_bipoles_from_quadpoles <- c()
  
  for(i in c(1:ncol(combs))){
    vec_bipoles_from_quadpoles[length(vec_bipoles_from_quadpoles) + 1] <-
      paste0(cyphers[combs[1,i]],cyphers[combs[2,i]])
  }
  
  if((cohort == "RR") || (cohort == "RRHF")){
    
    delta_df_RR <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein, 
                                            which_cases = cohort, version = 4, 
                                            SA_folder = SA_folder, 
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_RR
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
        bin2ascii_lead(.) %>%
        ascii2bin_lead(.)
      triquad_idx <- which(sum_string(num_elec) != 2)
      
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
        apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
        bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- 100*(max_values_bipoles[i] - value_if_not_included[i])/max_values_bipoles[i]
        }
      }
      
      delta_df_RR <- rbind(delta_df_RR,delta_vein)
      
    }
    
    delta_df_RR <- as.data.frame(delta_df_RR)
    colnames(delta_df_RR) <- names(max_values_bipoles)
    rownames(delta_df_RR) <- c("AN","AL","LA","IL","IN")
    
  }
  
  
  if((cohort == "HF") || (cohort == "RRHF")){
    
    delta_df_HF <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein,
                                            which_cases = cohort, version = 4,
                                            SA_folder = SA_folder,
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_HF
      
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
        bin2ascii_lead(.) %>%
        ascii2bin_lead(.)
      
      triquad_idx <- which(sum_string(num_elec) != 2)
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
        apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
        bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best 
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- 100*(max_values_bipoles[i] - value_if_not_included[i])/max_values_bipoles[i]
        }
      }
      
      delta_df_HF <- rbind(delta_df_HF,delta_vein)
      
    }
    
    delta_df_HF <- as.data.frame(delta_df_HF)
    colnames(delta_df_HF) <- names(max_values_bipoles)
    rownames(delta_df_HF) <- c("AN","AL","LA","IL","IN")
    
  }
  
  if(cohort == "RR")
    delta_df <- delta_df_RR
  else if(cohort == "HF")
    delta_df <- delta_df_HF
  else if(cohort == "RRHF")
    delta_df <- cbind(delta_df_RR,delta_df_HF)
  
  
  return(delta_df)
}

#' @description Function to calculate the maximum loss (in % of AT) of not
#' choosing the optimal for each vein and choosing the optimal for a different
#' vein. 
#' 
WorstCaseChoicePerVein <- function(which_cases = "HF", metric_option = "AT090",
                                   output = "%", flag_debugging = FALSE){
  veins <- c("AN","AL","LA","IL","LA")
  max_value <- 0
  for(i in 1:5){
  vein <- veins[i]
    
  optimal_designs <- Find_Optimal_Quadripole_Optimising(vein = vein,
                                     which_cases = which_cases,
                                     method = "max", response=100, 
                                     SA_folder = "default_HF_noPVTV",
                                     version = 4,
                                     metric_option = metric_option,
                                     flag_debugging = flag_debugging)
  
  delta <- ShowDeltaDesignVSCohort(cohort = which_cases, output = output,
                                   draw_option = "", 
                                   metric_option = metric_option,
                                   SA_folder = "default_HF_noPVTV",
                                   quad_design = optimal_designs$Quadripolar[[1]])
  max_delta <- delta[-i, ] %>% max() 
  max_value <- max(max_value,max_delta)
  }
  
  return(max_value)
}

WorstCaseChoicePerVeinRelative <- function(which_cases = "HF", metric_option = "AT090",
                                   flag_debugging = FALSE){
  veins <- c("AN","AL","LA","IL","LA")
  max_value <- 0
  for(i in 1:5){
    vein <- veins[i]
    
    optimal_designs <- Find_Optimal_Quadripole_Optimising(vein = vein,
                                                          which_cases = which_cases,
                                                          method = "max", response=100, 
                                                          SA_folder = "default_HF_noPVTV",
                                                          version = 4,
                                                          metric_option = metric_option,
                                                          flag_debugging = flag_debugging)
    
    delta <- RelativeDeltaDesignVSCohort(cohort = which_cases,
                                     metric_option = metric_option,
                                     SA_folder = "default_HF_noPVTV",
                                     quad_design = optimal_designs$Quadripolar[[1]])
    max_delta <- delta[-i, ] %>% max() 
    max_value <- max(max_value,max_delta)
  }
  
  return(max_value)
}

#' @description Function to show the improvement in reduction of metric_option
#' (in % or ms) of using a MPP lead design in a cohort of patients.
MPPImprovement <- function(cohort = "RR", metric_option, 
                                    SA_folder,
                                    quad_design = "1267",
                                    output = "%",
                                    flag_debugging = FALSE){
  
  source("/home/crg17/Desktop/scripts/multipole/R/common_functions.R")
  Load_Install_Packages(c("dplyr","nnet","ComplexHeatmap","circlize"))
  
  optimal_reductions_df <- c() 
  optimal_bipoles_df <- c()
  global_reductions_df <- c()
  
  
  optimal_quadpole <- quad_design
  cyphers <- strsplit(optimal_quadpole,"")[[1]]
  
  #We extract all the bipoles from the optimal quadpole
  combs <- combn(4,2)
  
  vec_bipoles_from_quadpoles <- c()
  
  for(i in c(1:ncol(combs))){
    vec_bipoles_from_quadpoles[length(vec_bipoles_from_quadpoles) + 1] <-
      paste0(cyphers[combs[1,i]],cyphers[combs[2,i]])
  }
  
  if((cohort == "RR") || (cohort == "RRHF")){
    
    delta_df_RR <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein, 
                                            which_cases = cohort, version = 4, 
                                            SA_folder = SA_folder,
                                            output = output,
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_RR
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
        bin2ascii_lead(.) %>%
        ascii2bin_lead(.)
      triquad_idx <- which(sum_string(num_elec) != 2)
      
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
        apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
        bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(NA,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- value_if_not_included[i]
        }
        else{
          delta_vein[i] <- max_values_bipoles[i]
        }
      }
      
      delta_df_RR <- rbind(delta_df_RR,delta_vein)
      
    }
    
    delta_df_RR <- as.data.frame(delta_df_RR)
    colnames(delta_df_RR) <- names(max_values_bipoles)
    rownames(delta_df_RR) <- c("AN","AL","LA","IL","IN")
    
  }
  
  
  if((cohort == "HF") || (cohort == "RRHF")){
    
    delta_df_HF <- c()
    
    for(vein in c("AN","AL","LA","IL","IN")){
      
      # Read the dataframes
      normalised_cohorts <- Read_dataframes(metric_option = metric_option,
                                            lead_option = vein,
                                            which_cases = cohort, version = 4,
                                            output = output,
                                            SA_folder = SA_folder,
                                            flag_debugging = flag_debugging)
      
      normalised_HF <- normalised_cohorts$normalised_HF
      
      
      # Get only bipoles
      num_elec <- rownames(normalised_HF) %>%
        bin2ascii_lead(.) %>%
        ascii2bin_lead(.)
      
      triquad_idx <- which(sum_string(num_elec) != 2)
      if(length(triquad_idx) > 0)
        normalised_HF <- normalised_HF[-triquad_idx,]
      normalised_ATs <- normalised_HF
      
      # Names of optimal bipoles
      optimal_bipoles <- normalised_ATs %>%
        apply(., 2, function(x) rownames(.)[as.vector(which.is.max(x))]) %>%
        bin2ascii_lead(.)
      optimal_bipoles_df <- rbind(optimal_bipoles_df,optimal_bipoles)
      
      # Value of optimal bipoles
      max_values_bipoles <- apply(normalised_ATs,2,max)
      optimal_reductions_df <- rbind(optimal_reductions_df, max_values_bipoles)
      rownames(optimal_reductions_df)[nrow(optimal_reductions_df)] <- vein
      
      # If it's in the optimal quadripole that's it, otherwise we take the best 
      # that we can
      bipoles_bin <- ascii2bin_lead(vec_bipoles_from_quadpoles)
      value_if_not_included <- apply(normalised_ATs[bipoles_bin,],2,max)
      global_reductions_df <- rbind(global_reductions_df, value_if_not_included)
      rownames(global_reductions_df)[nrow(global_reductions_df)] <- vein
      
      delta_vein <- rep(0,length(optimal_bipoles))
      
      for(i in c(1:length(optimal_bipoles))){
        if(!isbipoleinquadpole(optimal_bipoles[i], optimal_quadpole)){
          delta_vein[i] <- value_if_not_included[i]
        }
        else{
          delta_vein[i] <- max_values_bipoles[i]
        }
      }
      
      delta_df_HF <- rbind(delta_df_HF,delta_vein)
      
    }
    
    delta_df_HF <- as.data.frame(delta_df_HF)
    colnames(delta_df_HF) <- names(max_values_bipoles)
    rownames(delta_df_HF) <- c("AN","AL","LA","IL","IN")
    
  }
  
  if(cohort == "RR")
    delta_df <- delta_df_RR
  else if(cohort == "HF")
    delta_df <- delta_df_HF
  else if(cohort == "RRHF")
    delta_df <- cbind(delta_df_RR,delta_df_HF)
  
  
  
  return(delta_df)
}

PlotDensities <- function(curves, output = "%", flag_debugging = FALSE){
  Load_Install_Packages(c("ggplot2"))
  
  if(curves == "all"){
    PlotDensities(curves = "baseline_veins", output = output, flag_debugging = flag_debugging)
    PlotDensities(curves = "baseline_patients", output = output, flag_debugging = flag_debugging)
    PlotDensities(curves = "single_HF_veins", output = output, flag_debugging = flag_debugging)
    PlotDensities(curves = "single_HF_patients", output = output, flag_debugging = flag_debugging)
    PlotDensities(curves = "single_RR_veins", output = output, flag_debugging = flag_debugging)
    PlotDensities(curves = "single_RR_patients", output = output, flag_debugging = flag_debugging)
    PlotDensities(curves = "central", output = output, flag_debugging = flag_debugging)
    
  }
  else{
  if(output == "%"){
    x_bounds <- c(0,30)
    y_bounds <- c(0,0.25)
  }
    if(output == "ms"){
      x_bounds <- c(0,35)
      y_bounds <- c(0,0.16)
    }
  if(curves == "baseline_veins"){
    source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    baseline_df <- BaselineSingleElectrodeReduction(output = output,
                                                    SA_folder = "default_HF_noPVTV",
                                                    flag_debugging = flag_debugging)
    colours_vec <- c("#FC0000", "#ED7c31", "#FFBF00", "#00B04F","#00B0F0")
    
    p<- ggplot(baseline_df, aes(x=Reduction, color=Vein, fill = Vein)) +
      geom_density(alpha = 0.5)+
      labs(title="Density curves of single electrode activation",
           x=paste0("AT090 reduction (",output,")"), y = "Density")
    p <- p + scale_fill_manual(values = colours_vec) 
    p <- p + scale_color_manual(values = colours_vec) 
    p <- p + theme_classic(base_size = 40) 
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank())
    
    if(output == "%"){
      file_name <- "/home/crg17/Pictures/baseline_veins_rel.png"
      }
    else{
      file_name <- "/home/crg17/Pictures/baseline_veins_abs.png"
    }
    p <- p + xlim(x_bounds) + ylim(y_bounds)
    print(p)
    
    ggsave(filename = file_name,
           width = 20,height = 10)
  }
  if(curves == "baseline_patients"){
    source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    baseline_df <- BaselineSingleElectrodeReduction(output = output,
                                                    SA_folder = "default_HF_noPVTV",
                                                    flag_debugging = flag_debugging)

    p<- ggplot(baseline_df, aes(x=Reduction)) +
      labs(title="Density curves of single electrode activation",
           x=paste0("AT090 reduction (",output,")"), y = "Density")
    p <- p + geom_density(aes(group = Patient))
  
    p <- p + theme_classic(base_size = 40) 
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank())
    if(output == "%"){
      file_name <- "/home/crg17/Pictures/baseline_patients_rel.png"
    }
    else{
      file_name <- "/home/crg17/Pictures/baseline_patients_abs.png"
    }
    y_bounds[2] <- 0.85
    p <- p + xlim(x_bounds) + ylim(y_bounds)
    print(p)
    ggsave(filename = file_name,
           width = 20,height = 10)
  }
  if(curves == "single_HF_veins"){
    source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    global_optimal_HF <- MPPDesignReduction(cohort = "HF", output = output)
    colours_vec <- c("#FC0000", "#ED7c31", "#FFBF00", "#00B04F","#00B0F0")
    
    p<- ggplot(global_optimal_HF, aes(x=Reduction, color=Vein, fill = Vein)) +
      geom_density(alpha = 0.5)+
      labs(title="Density curves of global optimal MPP design in HF",
           x=paste0("AT090 reduction (",output,")"), y = "Density")
    p <- p + scale_fill_manual(values = colours_vec) 
    p <- p + scale_color_manual(values = colours_vec) 
    p <- p + theme_classic(base_size = 40) 
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank())
    if(output == "%"){
      file_name <- "/home/crg17/Pictures/single_HF_veins_rel.png"
    }
    else{
      file_name <- "/home/crg17/Pictures/single_HF_veins_abs.png"
    }
    p <- p + xlim(x_bounds) + ylim(y_bounds)
    print(p)
    
    ggsave(filename = file_name,
           width = 20,height = 10)
  }
  if(curves == "single_HF_patients"){
    source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    global_optimal_HF <- MPPDesignReduction(cohort = "HF", output = output)
    
    p <- ggplot(global_optimal_HF, aes(x=Reduction)) +
      labs(title=paste0("Proportion of lead configurations with a reduction of AT090\nusing the cohort-based optimal MPP designs in the HF cohort"),
           x=paste0("AT090 reduction (",output,")"), y = "Lead configurations (%)")
    p <- p + geom_density(aes(group = Patient))
    
    p <- p + theme_classic(base_size = 40) 
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank())
    if(output == "%"){
      file_name <- "/home/crg17/Pictures/single_HF_patients_rel.png"
    }
    else{
      file_name <- "/home/crg17/Pictures/single_HF_patients_abs.png"
    }
    y_bounds[2] <- 0.85
    p <- p + xlim(x_bounds) + ylim(y_bounds)
    
    print(p)
    
    ggsave(filename = file_name,
           width = 20,height = 10)
  }
  if(curves == "single_RR_veins"){
    source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    global_optimal_HF <- MPPDesignReduction(cohort = "RR", output = output)
    colours_vec <- c("#FC0000", "#ED7c31", "#FFBF00", "#00B04F","#00B0F0")
    
    p<- ggplot(global_optimal_HF, aes(x=Reduction, color=Vein, fill = Vein)) +
      geom_density(alpha = 0.5)+
      labs(title="Density curves of global optimal MPP design in RR",
           x=paste0("AT090 reduction (",output,")"), y = "Density")
    p <- p + scale_fill_manual(values = colours_vec) 
    p <- p + scale_color_manual(values = colours_vec) 
    p <- p + theme_classic(base_size = 40) 
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank())
    if(output == "%"){
      file_name <- "/home/crg17/Pictures/single_RR_veins_rel.png"
    }
    else{
      file_name <- "/home/crg17/Pictures/single_RR_veins_abs.png"
    }
    p <- p + xlim(x_bounds) + ylim(y_bounds)
    
    print(p)
    ggsave(filename = file_name,
           width = 20,height = 10)
  }
  if(curves == "single_RR_patients"){
    source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    global_optimal_HF <- MPPDesignReduction(cohort = "RR", output = output)
    
    p<- ggplot(global_optimal_HF, aes(x=Reduction)) +
      labs(title="Density curves of global optimal MPP design in RR",
           x=paste0("AT090 reduction (",output,")"), y = "Density")
    p <- p + geom_density(aes(group = Patient))
    
    p <- p + theme_classic(base_size = 40) 
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank())
    if(output == "%"){
      file_name <- "/home/crg17/Pictures/single_RR_patients_rel.png"
    }
    else{
      file_name <- "/home/crg17/Pictures/single_RR_patients_abs.png"
    }
    
    y_bounds[2] <- 0.85
    p <- p + xlim(x_bounds) + ylim(y_bounds)
    
    print(p)
    ggsave(filename = file_name,
           width = 20,height = 10)
  }
  if(curves == "central"){
    source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    
    HF_design <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                    vein = "ALL",
                                                    which_cases = "HF",
                                                    method = "max",
                                                    SA_folder = "default_HF_noPVTV",
                                                    version = 4)
    
    global_optimal_HF <- MPPDesignReduction(cohort = "HF",
                                            MPP_design = HF_design$Quadripolar[1],
                                            output = output)
    personalised_optimal_HF <- Optimal12ElectrodeReduction(cohort = "HF",
                                                           output = output)
    
    RR_design <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                    vein = "ALL",
                                                    which_cases = "RR",
                                                    method = "max",
                                                    SA_folder = "default_HF_noPVTV",
                                                    version = 4)
    global_optimal_RR <- MPPDesignReduction(cohort = "RR",
                                            MPP_design = RR_design$Quadripolar[1],
                                            output = output)
    personalised_optimal_RR <- Optimal12ElectrodeReduction(cohort = "RR",
                                                           output = output)
    
    
    global_HF_vec <- global_optimal_HF$Reduction
    person_HF_vec <- c(t(personalised_optimal_HF))
    global_RR_vec <- global_optimal_RR$Reduction
    person_RR_vec <- c(t(personalised_optimal_RR))
    
    final_df <- as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(final_df) <- c("Reduction", "Strategy")
    
    final_df <- rbind(final_df,data.frame(Reduction = global_HF_vec,
                                          Strategy = rep("HF-based MPP design used in HF",
                                                         length(global_HF_vec))))
    final_df <- rbind(final_df,data.frame(Reduction = person_HF_vec,
                                          Strategy = rep("\nPersonalised MPP design used in HF\n",
                                                         length(person_HF_vec))))
    final_df <- rbind(final_df,data.frame(Reduction = global_RR_vec,
                                          Strategy = rep("RR-based MPP design used in RR",
                                                         length(global_RR_vec))))
    final_df <- rbind(final_df,data.frame(Reduction = person_RR_vec,
                                          Strategy = rep("\nPersonalised MPP design used in RR\n",
                                                         length(person_RR_vec))))
    p <- ggplot(final_df, aes(x=Reduction)) 
      
    p <- p + labs(title="Proportion of lead configurations with a reduction of AT090\nusing patient and cohort-based optimal MPP designs",
           x=paste0("AT090 reduction (",output,")"), y = "Lead configurations (%)")
    p <- p + geom_line(size = 2, aes(colour = Strategy, linetype = Strategy),
                       stat = "density", position = "identity")
    p <- p + scale_y_continuous(breaks = c(0, .05, .1, .15),
                                labels =c ("0", "5", "10", "15") )
    p <- p + scale_color_manual(values=c("HF-based MPP design used in HF" = "#0000FF", # Blue
                                         "RR-based MPP design used in RR" = "#FF0000", # Red
                                         #"\nSingle-electrode pacing\n" = "#000000", # Black
                                         "\nPersonalised MPP design used in HF\n" = "#92D050", # Light green
                                         "\nPersonalised MPP design used in RR\n" = "#92D050")) # Light green
    p <- p + scale_linetype_manual(values=c("HF-based MPP design used in HF" = "solid",
                                            "RR-based MPP design used in RR" = "dotted",
                                            #"\nSingle-electrode pacing\n" = "solid",
                                            "\nPersonalised MPP design used in HF\n" = "solid",
                                            "\nPersonalised MPP design used in RR\n" = "dotted"))
    p <- p + theme_classic(base_size = 40) 
    p <- p + theme(plot.title = element_text(hjust = 0.1),
                   panel.background = element_blank())
    print(p)
    if(output == "%"){
      file_name <- "/home/crg17/Pictures/central_figure_rel.png"
    }
    else{
      file_name <- "/home/crg17/Pictures/central_figure_abs.png"
    }
    ggsave(filename = file_name,
           width = 20,height = 10)
  }
  if(curves == "optimal_vs_cohort"){
      source("/home/crg17/Desktop/scripts/multipole/R/postprocessing.R")
    
    HF_design <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                    vein = "ALL",
                                                    which_cases = "HF",
                                                    method = "max",
                                                    SA_folder = "default_HF_noPVTV",
                                                    version = 4)
    
      optimal_HF_in_HF <- MPPDesignReduction(cohort = "HF",
                                             MPP_design = HF_design$Quadripolar[1],
                                             output = output)
      optimal_HF_in_RR <- MPPDesignReduction(cohort = "RR",
                                             MPP_design = HF_design$Quadripolar[1],
                                             output = output)
      
      RR_design <- Find_Optimal_Quadripole_Optimising(metric_option = "AT090",
                                                      vein = "ALL",
                                                      which_cases = "RR",
                                                      method = "max",
                                                      SA_folder = "default_HF_noPVTV",
                                                      version = 4)
      
      optimal_RR_in_HF <- MPPDesignReduction(cohort = "HF",
                                             MPP_design = RR_design$Quadripolar[1],
                                             output = output)
      optimal_RR_in_RR <- MPPDesignReduction(cohort = "RR",
                                             MPP_design = RR_design$Quadripolar[1],
                                             output = output)
      
      final_df <- as.data.frame(matrix(nrow = 0, ncol = 2))
      colnames(final_df) <- c("Reduction", "Strategy")
      
      final_df <- rbind(final_df,data.frame(Reduction = optimal_HF_in_HF$Reduction,
                                            Strategy = rep("HF-based MPP design used in HF",
                                                           nrow(optimal_HF_in_HF))))
      final_df <- rbind(final_df,data.frame(Reduction = optimal_HF_in_RR$Reduction,
                                            Strategy = rep("\nHF-based MPP design used in RR\n",
                                                           nrow(optimal_HF_in_RR))))
      final_df <- rbind(final_df,data.frame(Reduction = optimal_RR_in_HF$Reduction,
                                            Strategy = rep("RR-based MPP design used in HF",
                                                           nrow(optimal_RR_in_HF))))
      final_df <- rbind(final_df,data.frame(Reduction = optimal_RR_in_RR$Reduction,
                                            Strategy = rep("\nRR-based MPP design used in RR\n",
                                                           nrow(optimal_RR_in_RR))))
      
      p <- ggplot(final_df, aes(x = Reduction)) 
      
      p <- p + labs(title="Proportion of lead configurations with reduction in AT090\nin population-based optimal MPP designs",
                    x=paste0("AT090 reduction (",output,")"), y = "Lead configurations (%)")
      # p <- p + geom_line(stat = "density", size = 1)
      p <- p + geom_line(size = 2, aes(colour = Strategy, linetype = Strategy),
                            stat = "density", position = "identity")
      p <- p + scale_y_continuous(breaks = c(0, .05, .1, .15),
                                  labels =c ("0", "5", "10", "15") )
      p <- p + scale_color_manual(values=c("HF-based MPP design used in HF" = "#0000FF",
                                           "\nHF-based MPP design used in RR\n" = "#0000FF",
                                           "RR-based MPP design used in HF" = "#FF0000",
                                           "\nRR-based MPP design used in RR\n" = "#FF0000"))
      p <- p + scale_linetype_manual(values=c("HF-based MPP design used in HF" = "solid",
                                              "\nHF-based MPP design used in RR\n" = "dotted",
                                              "RR-based MPP design used in HF" = "solid",
                                              "\nRR-based MPP design used in RR\n" = "dotted"))
      p <- p + theme_classic(base_size = 40) 
      p <- p + theme(plot.title = element_text(hjust = 0.2),
                     panel.background = element_blank())
      print(p)
      if(output == "%"){
        file_name <- "/home/crg17/Pictures/optimal_vs_cohort_rel.png"
      }
      else{
        file_name <- "/home/crg17/Pictures/optimal_vs_cohort_abs.png"
      }
      ggsave(filename = file_name,
             width = 20,height = 10)
    }
  }
}
