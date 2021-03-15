library("RColorBrewer") # For the colours 
#library("gplots")
require("pheatmap")
library("d3heatmap")
library("viridis")
require("ComplexHeatmap")
library("circlize")
library("dendextend")
library("dendsort")
library("seriation")
library("ggplot2")
library("nnet")
require("stringr")


read_dataframes <- function(metric_option, lead_option, which_cases, version, SA_folder, RV_midseptum = FALSE){
  # if(exists("RR_cases")) rm(RR_cases)
  # if(exists("HF_cases")) rm(HF_cases)
  # if(exists("normalised_RR")) rm(normalised_RR)
  # if(exists("normalised_HF")) rm(normalised_HF)
  if(version == 1){
    RR_cases <- read.table(paste0("/data/multipoles_files/baseline/multipole_H_",lead_option,"_",metric_option,".dat"), header=TRUE)
    HF_cases <- read.table(paste0("/data/multipoles_files/baseline/multipole_HF_",lead_option,"_",metric_option,".dat"), header=TRUE)
    
    RR_RVapex <- read.table("/data/multipoles_files/baseline/AT_metric_RVapex.dat", header=TRUE)
    HF_RVapex <- read.table("/data/multipoles_files/baseline/AT_metric_RVapex_HF.dat", header=TRUE)
    
    # Normalising
    
    normalised_RR = 100*(1-mapply('/', RR_cases, RR_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))
    normalised_HF = 100*(1-mapply('/', HF_cases, HF_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))
    
    colnames(normalised_RR) <- c("RR01","RR02","RR03","RR04","RR05","RR06","RR07","RR08","RR09","RR10","RR11","RR12","RR13","RR14","RR15","RR16","RR17","RR18","RR19","RR20")
    rownames(normalised_RR) <- rownames(RR_cases)
    
    colnames(normalised_HF) <- c("HF01","HF02","HF03","HF04","HF05","HF06","HF07","HF08","HF09","HF10","HF11","HF12","HF13","HF14","HF15","HF16","HF17","HF18","HF19","HF20","HF21","HF22","HF23","HF24")
    rownames(normalised_HF) <- rownames(HF_cases)
  
  }
  
  else if(version == 2){
    heart_cases <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24")
    
    if(exists("temp_HF")) rm(temp_HF)
    if(exists("HF_cases")) rm(HF_cases)
    
    temp_HF <- NA*read.table(paste0("/data/SA_multipole/FEC_70_values/HF/",heart_cases[1],"/multipole_",metric_option,".dat"), header = TRUE)
    HF_cases <- temp_HF[,c(1,2)] 
    
    for(i in c(1,length(heart_cases))){
      if(exists("temp_HF")) rm(temp_HF)
      temp_HF <- read.table(paste0("/data/SA_multipole/FEC_70_values/HF/",heart_cases[i],"/multipole_",metric_option,".dat"), header = TRUE)
      HF_cases[,heart_cases[i]] <- temp_HF[,lead_option]
    }
    
    RR_cases <- read.table(paste0("/data/multipoles_files/baseline/multipole_H_",lead_option,"_",metric_option,".dat"), header=TRUE)

    RR_RVapex <- read.table("/data/multipoles_files/baseline/AT_metric_RVapex.dat", header=TRUE)
    HF_RVapex <- read.table("/data/multipoles_files/baseline/AT_metric_RVapex_HF.dat", header=TRUE)
    
    # Normalising
    
    normalised_RR = 100*(1-mapply('/', RR_cases, RR_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))
    normalised_HF = 100*(1-mapply('/', HF_cases, HF_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))
    
    colnames(normalised_RR) <- c("RR01","RR02","RR03","RR04","RR05","RR06","RR07","RR08","RR09","RR10","RR11","RR12","RR13","RR14","RR15","RR16","RR17","RR18","RR19","RR20")
    rownames(normalised_RR) <- rownames(RR_cases)
    
    colnames(normalised_HF) <- paste0("HF",heart_cases)
    rownames(normalised_HF) <- rownames(HF_cases)
    
  }
  
  else if(version == 3){
    
    singlerotlead <- c()
    RR_cases <- c()
    HF_cases <- c()
    # if(exists("singleheart")) rm(singleheart)
    # if(exists("normalised_RR")) rm(normalised_RR)
    if(which_cases == "RR" || which_cases == "RRHF"){
    for(heart in c("01","02","03","04","05","06","07","08","09",10:20)){
      singleheart <- read.table(paste0("/data/SA_multipole/FEC_70_values/h/",heart,"/multipole_",metric_option,".dat"), header = TRUE)
      if(lead_option == "AN")
        singlerotlead <- singleheart[,1]
      if(lead_option == "AL")
        singlerotlead <- singleheart[,2]
      if(lead_option == "LA")
        singlerotlead <- singleheart[,3]
      if(lead_option == "PL")
        singlerotlead <- singleheart[,4]
      if(lead_option == "PO")
        singlerotlead <- singleheart[,5]
      
      RR_cases <- cbind(RR_cases,singlerotlead)    
    }
      RR_RVapex <- read.table("/data/multipoles_files/baseline/AT_metric_RVapex.dat", header=TRUE)
      
      # Normalising
      
      normalised_RR = 100*(1-mapply('/', as.data.frame(RR_cases), RR_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))

      colnames(normalised_RR) <- c("RR01","RR02","RR03","RR04","RR05","RR06","RR07","RR08","RR09","RR10","RR11","RR12","RR13","RR14","RR15","RR16","RR17","RR18","RR19","RR20")
      rownames(normalised_RR) <- rownames(singleheart)
    }
    
    if(which_cases == "HF" || which_cases == "RRHF"){
    # if(exists("singleheart")) rm(singleheart)
    HF_cases <- c()
    for(heart in c("01","02","03","04","05","06","07","08","09",10:24)){
      monodipole <- read.table(paste0("/data/SA_multipole/FEC_70_values/HF/",heart,"/monopole_dipole/multipole_",metric_option,".dat"), header = TRUE)
      triquapole <- read.table(paste0("/data/SA_multipole/FEC_70_values/HF/",heart,"/tripole_quadripole/multipole_",metric_option,".dat"), header = TRUE)
      
      singleheart <- rbind(monodipole,triquapole)

      if(lead_option == "AN")
        singlerotlead <- singleheart[,1]
      if(lead_option == "AL")
        singlerotlead <- singleheart[,2]
      if(lead_option == "LA")
        singlerotlead <- singleheart[,3]
      if(lead_option == "PL")
        singlerotlead <- singleheart[,4]
      if(lead_option == "PO")
        singlerotlead <- singleheart[,5]
      
      HF_cases <- cbind(HF_cases,singlerotlead)    
    }
    
    
    HF_RVapex <- read.table("/data/multipoles_files/baseline/AT_metric_RVapex_HF.dat", header=TRUE)
    
    # Normalising
    
    normalised_HF = 100*(1-mapply('/', as.data.frame(HF_cases), HF_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))
    
    colnames(normalised_HF) <- c("HF01","HF02","HF03","HF04","HF05","HF06","HF07","HF08","HF09","HF10","HF11","HF12","HF13","HF14","HF15","HF16","HF17","HF18","HF19","HF20","HF21","HF22","HF23","HF24")
    rownames(normalised_HF) <- rownames(singleheart)
    }
    

    
  }
  
  else if(version == 4){
    
    singlerotlead <- c()
    RR_cases <- c()
    HF_cases <- c()
    # if(exists("singleheart")) rm(singleheart)
    # if(exists("normalised_RR")) rm(normalised_RR)
    if(which_cases == "RR" || which_cases == "RRHF"){
      for(heart in c("01","02","03","04","05","06","07","08","09",10:20)){
        singleheart <- read.table(paste0("/data/SA_multipole/",SA_folder,"/h/",heart,"/multipole_AT1090_leadcorrected.dat"), header = TRUE)
        singlerotlead <- singleheart[,lead_option]
        
        RR_cases <- cbind(RR_cases,singlerotlead)    
      }
      # print(paste0("/data/SA_multipole/",SA_folder,"/h/multipole_RVapex.dat"))
      if(!RV_midseptum)
      RR_RVapex <- read.table(paste0("/data/SA_multipole/",SA_folder,"/h/multipole_RVapex.dat"), header = TRUE)
      else
        RR_RVapex <- read.table(paste0("/data/SA_multipole/",SA_folder,"/h/AT_metric_midseptum_h.dat"), header = TRUE)
      # Normalising
      
      normalised_RR = 100*(1-mapply('/', as.data.frame(RR_cases), RR_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))
      
      colnames(normalised_RR) <- c("RR01","RR02","RR03","RR04","RR05","RR06","RR07","RR08","RR09","RR10","RR11","RR12","RR13","RR14","RR15","RR16","RR17","RR18","RR19","RR20")
      if(version != 4){
        rownames(normalised_RR) <- rownames(singleheart)
      }
      else{
        rownames(normalised_RR) <- ascii2bin_lead(bin2ascii_lead(singleheart$lead))
      }
    }
    
    if(which_cases == "HF" || which_cases == "RRHF"){
      # if(exists("singleheart")) rm(singleheart)
      HF_cases <- c()
      for(heart in c("01","02","03","04","05","06","07","08","09",10:24)){
        # print(paste0("/data/SA_multipole/",SA_folder,"/HF/",heart,"/multipole_AT1090_leadcorrected.dat"))
        singleheart <- read.table(paste0("/data/SA_multipole/",SA_folder,"/HF/",heart,"/multipole_AT1090_leadcorrected.dat"), header = TRUE)

        singlerotlead <- singleheart[,lead_option]
        
        
        HF_cases <- cbind(HF_cases,singlerotlead)    
      }
      
      # print(paste0("/data/SA_multipole/",SA_folder,"/HF/multipole_RVapex.dat"))
      if(!RV_midseptum)
        HF_RVapex <- read.table(paste0("/data/SA_multipole/",SA_folder,"/HF/multipole_RVapex.dat"), header=TRUE)
      else
        HF_RVapex <- read.table(paste0("/data/SA_multipole/",SA_folder,"/HF/AT_metric_midseptum_HF.dat"), header=TRUE)
      # Normalising
      
      normalised_HF = 100*(1-mapply('/', as.data.frame(HF_cases), HF_RVapex[,as.integer(metric_option == "QRS") + 2*as.integer(metric_option == "TATLV") + 3*as.integer(metric_option == "AT1090")]))
      
      colnames(normalised_HF) <- c("HF01","HF02","HF03","HF04","HF05","HF06","HF07","HF08","HF09","HF10","HF11","HF12","HF13","HF14","HF15","HF16","HF17","HF18","HF19","HF20","HF21","HF22","HF23","HF24")
      if(version != 4){
        rownames(normalised_HF) <- rownames(singleheart)
      }
      else{
        rownames(normalised_HF) <- ascii2bin_lead(bin2ascii_lead(singleheart$lead))
      }
    }
    
  }
  
  if(which_cases == "RRHF")
    return(list("normalised_RR" = normalised_RR, "normalised_HF" = normalised_HF))
  else if(which_cases == "RR")
    return(list("normalised_RR" = normalised_RR))
  else if(which_cases == "HF")
    return(list("normalised_HF" = normalised_HF))
}


hier_clust_multipole <- function(metric_option = "AT1090", lead_option = "PO", col_clusters = 3, row_clusters = 10, draw_hm = FALSE, thresh_plot = 100, which_cases = "RRHF", max_bar = 30, version, output = "number", SA_folder = SA_folder, RV_midseptum = RV_midseptum){

  if(metric_option == "AT1090"){
  plottitle = "AT 10%-90% "
} else if(metric_option == "QRS"){
  plottitle = "QRS duration "
} else if(metric_option == "TATLV"){
  plottitle = "TAT of the LV "
}

if(lead_option == "AN"){
  plottitle = paste0(plottitle,"pacing in the anterior position")
} else if(lead_option == "AL"){
  plottitle = paste0(plottitle,"pacing in the antero-lateral position")
} else if(lead_option == "LA"){
  plottitle = paste0(plottitle,"pacing in the lateral position")
} else if(lead_option == "PL"){
  plottitle = paste0(plottitle,"pacing in the postero-lateral position")
} else if(lead_option == "PO"){
  plottitle = paste0(plottitle,"pacing in the posterior position")
}

plottitle = paste0(plottitle,"\n")

# We read the dataframe
if(exists("normalised_ATs")) rm(normalised_ATs)
normalised_cohorts <- read_dataframes(metric_option = metric_option, lead_option = lead_option, which_cases = which_cases, version = version, SA_folder = SA_folder, RV_midseptum = RV_midseptum)
normalised_RR <- normalised_cohorts$normalised_RR
normalised_HF <- normalised_cohorts$normalised_HF
if(which_cases == "RR" || which_cases == "RRHF"){
  triquad_idx <- which(sum_string(ascii2bin_lead(bin2ascii_lead(rownames(normalised_RR)))) != 2)
  if(length(triquad_idx) > 0)
    normalised_RR <- normalised_RR[-triquad_idx,]
}
if(which_cases == "HF" || which_cases == "RRHF"){
  triquad_idx <- which(sum_string(ascii2bin_lead(bin2ascii_lead(rownames(normalised_HF)))) != 2)
  if(length(triquad_idx) > 0)
    normalised_HF <- normalised_HF[-triquad_idx,]
}
# Coloring


RR_mass <- 1.05*read.table("/media/crg17/Seagate Backup Plus Drive/CT_cases/forall/meshes_vol.dat", header = TRUE)
RR_mass2color <- RR_mass$RV
HF_mass <- 1.05*0.001*read.table("/media/crg17/Seagate Backup Plus Drive/CT_cases/forall/meshes_vol_HF.dat", header = TRUE)
HF_mass2color <- HF_mass$RV

RR_volumes <- read.table("/media/crg17/Seagate Backup Plus Drive/CT_cases/forall/chambers_volumes.txt", header = TRUE)
RR_volume2color <- RR_volumes$RV_mL
HF_volumes <- read.table("/media/crg17/Seagate Backup Plus Drive/CT_cases/forall/chambers_volumes_HF.csv", header = TRUE, sep = ",")
HF_volume2color <- HF_volumes$RV

RR_quotient <- RR_mass2color/RR_volume2color
HF_quotient <- HF_mass2color/HF_volume2color

RR_big <- RR_volume2color > 250 
RR_small <- RR_volume2color < 100
HF_big <- HF_volume2color > 250 
HF_small <- HF_volume2color < 100

if(exists("volumes2color")) rm(volumes2color)
if(which_cases == "RRHF"){
  volumes2color <- sprintf("%s",1 - c(RR_big,HF_big) + c(RR_small,HF_small)) # 0 for big, 1 for medium, 2 for small
}
if(which_cases == "RR"){
  volumes2color <- sprintf("%s",1 - RR_big + RR_small) # 0 for big, 1 for medium, 2 for small
}
if(which_cases == "HF"){
  volumes2color <- sprintf("%s",1 - HF_big + HF_small) # 0 for big, 1 for medium, 2 for small
}

for(i in 1:length(volumes2color)){
  if(volumes2color[i] == "0"){
    volumes2color[i] <- "Big"
  }
  else if(volumes2color[i] == "1"){
    volumes2color[i] <- "Medium"
  }
  else if(volumes2color[i] == "2"){
    volumes2color[i] <- "Small"
  }
}

if(exists("normalised_ATs")) rm(normalised_ATs)
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


mat_col <- data.frame(Size = volumes2color, Diagnosis = hf2color)
rownames(mat_col) <- colnames(normalised_ATs)


# Colors for rows
leadpos2color <- rownames(normalised_ATs)

cut_base = 2
cut_apex = 7

for(j in 1:length(leadpos2color)){
  base_act <- lengths(regmatches(substr(leadpos2color[j],1,cut_base), gregexpr("1", substr(leadpos2color[j],1,cut_base)))) 
  mid_act <- lengths(regmatches(substr(leadpos2color[j],cut_base+1,cut_apex-1), gregexpr("1", substr(leadpos2color[j],cut_base+1,cut_apex-1)))) 
  apex_act <- lengths(regmatches(substr(leadpos2color[j],cut_apex,8), gregexpr("1", substr(leadpos2color[j],cut_apex,8)))) 
  
  
  if(apex_act >= 1 && base_act == 0 && mid_act == 0){
    leadpos2color[j] <- "Apical"
  }
  else if(base_act >= 1 && apex_act == 0 && mid_act == 0){
    leadpos2color[j] <- "Basal"
  }
  else if(apex_act >= 1 && base_act >= 1 && mid_act == 0){
    leadpos2color[j] <- "ApicoBasal"
  }
  else{
    leadpos2color[j] <- "Medial"
  }
}

mat_row <- data.frame(Position = leadpos2color)
rownames(mat_row) <- rownames(normalised_ATs)


# To cut the data set
# 
#  pheatmap(normalised_ATs,
#          cutree_rows = row_clusters,
#          cutree_cols = col_clusters,
# 
#          clustering_distance_rows = "euclidean", # distance measure used in clustering rows. Possible values are "correlation" for Pearson correlation and all the distances supported by dist, such as "euclidean", etc. If the value is none of the above it is assumed that a distance matrix is provided.
#          clustering_distance_cols = "euclidean", # clustering method used. Accepts the same values as hclust
#          clustering_method = "ward.D2", # clustering method used. Accepts the same values as hclust
# 
#          border_color = NA,
#          color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
#          annotation_col = mat_col,
#          annotation_row = mat_row,
#          annotation_colors = my_colour,
#          annotation_legend = TRUE,
# 
#          cellwidth = NA,
#          cellheight = NA,
#          # main = "Improvement (%) of AT 10%\u201390%, anterior", # Title
#          main = paste0("Improvement (%) of ",metric_option,", ",lead_option),
#          fontSize = 10,
#          display_numbers = TRUE,
# 
#          # breaks = seq(min(normalised_ATs), max(normalised_ATs), length.out = 101),
#          breaks = seq(0,20, length.out = 101),
#          legend = TRUE
#          )


# Complex heatmap

row_dend <- as.dendrogram(dendsort(hclust(dist(normalised_ATs))))
row_dend <- color_branches(row_dend, k = row_clusters)

col_dend <- as.dendrogram(dendsort(hclust(dist(t(normalised_ATs)))))
col_dend <- color_branches(col_dend, k = col_clusters)



if(draw_hm == TRUE){
  
  
  o1 = seriate(dist(normalised_ATs), method = "ARSA")
  # o1 = seriate(dist(normalised_ATs), method = "BBURCG")
  # o1 = seriate(dist(normalised_ATs), method = "BBWRCG")
  o1 = seriate(dist(normalised_ATs), method = "GW")
  o1 = seriate(dist(normalised_ATs), method = "GW_average")
  o1 = seriate(dist(normalised_ATs), method = "GW_single")
  o1 = seriate(dist(normalised_ATs), method = "GW_ward")
  o1 = seriate(dist(normalised_ATs), method = "MDS_angle")
  o1 = seriate(dist(normalised_ATs), method = "MDS")
  o1 = seriate(dist(normalised_ATs), method = "Identity")
  o1 = seriate(dist(normalised_ATs), method = "HC_ward")
  # o1 = seriate(dist(normalised_ATs), method = "HC_single")
  # o1 = seriate(dist(normalised_ATs), method = "HC_complete")
  # o1 = seriate(dist(normalised_ATs), method = "HC_average")
  # o1 = seriate(dist(normalised_ATs), method = "HC")
  # o1 = seriate(dist(normalised_ATs), method = "MDS_metric")
  # o1 = seriate(dist(normalised_ATs), method = "MEDS_nonmetric")
  # o1 = seriate(dist(normalised_ATs), method = "OLO")
  # o1 = seriate(dist(normalised_ATs), method = "OLO_average")
  # o1 = seriate(dist(normalised_ATs), method = "OLO_complete")
  # o1 = seriate(dist(normalised_ATs), method = "OLO_single")
  # o1 = seriate(dist(normalised_ATs), method = "OLO_ward")
  # o1 = seriate(dist(normalised_ATs), method = "QAP_2SUM")
  # o1 = seriate(dist(normalised_ATs), method = "Spectral_norm")
  # o1 = seriate(dist(normalised_ATs), method = "Spectral")
  # o1 = seriate(dist(normalised_ATs), method = "SA")
  # o1 = seriate(dist(normalised_ATs), method = "Random")
  # o1 = seriate(dist(normalised_ATs), method = "R2E")
  # o1 = seriate(dist(normalised_ATs), method = "QAP_LS")
  # o1 = seriate(dist(normalised_ATs), method = "QAP_Inertia")
  # o1 = seriate(dist(normalised_ATs), method = "QAP_BAR")
  # o1 = seriate(dist(normalised_ATs), method = "SPIN_NH")
  # o1 = seriate(dist(normalised_ATs), method = "SPIN_STS")
  o1 = seriate(dist(normalised_ATs), method = "TSP")
  # o1 = seriate(dist(normalised_ATs), method = "VAT")
  #o2 = seriate(dist(t(normalised_ATs)), method = "TSP")
  dend = dendsort(hclust(dist(normalised_ATs)))
  
  rownames(normalised_ATs) <- str_pad(strsplit(bin2ascii_lead(rownames(normalised_ATs)),","),side = "left",width = 4)
  hm=Heatmap(normalised_ATs,
             name = "RR+HF", # Internal name
             
             cluster_rows = row_dend,
             # cluster_rows = as.dendrogram(o1[[1]]),
             clustering_distance_rows = "euclidean",
             clustering_method_rows = "ward.D2",
             #row_dend_reorder = TRUE,
             #row_names_gp = gpar(col = ifelse(leadpos2color == "Medial", "black", "orange")),
             row_split = row_clusters,
             # row_order = get_order(o1),
             # split=x,
             # right_annotation = rowAnnotation(
             #   Position = leadpos2color,
             #   # col = list(Position = c(Apical = "#FF6D00", Basal = "#E65100", ApicoBasal = "#FFFF00", Medial = "#FFB74D")),
             #   col = list(Position = c(Apical = "#E1C62F", Basal = "#FFFF00", ApicoBasal = "#212121", Medial = "#BCAAA4")),
             #   annotation_legend_param = list(
             #     Position = list(
             #       title = "Lead position",
             #       labels = c("Apical", "Basal", "Mixed apico-basal", "Others"),
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
             column_names_gp = gpar(col = ifelse(hf2color == "HF", "blue", "red"), fontsize = 18, fontfamily = "Helvetica"),
             column_split = col_clusters,
             #column_order = get_order(o,2),
             # top_annotation = HeatmapAnnotation(
             #   LV_mass = anno_barplot(
             #     bar_width = 1,
             #     height = unit(40,"points"),
             #     c(RR_mass2color,HF_mass2color),
             #     gp = gpar(fill = colorRampPalette(brewer.pal(n = 9, name = "Purples"))(100))
             #     ),
             #     show_annotation_name = TRUE,
             #   LV_vol = anno_barplot(
             #     bar_width = 1,
             #     height = unit(40,"points"),
             #     c(RR_volume2color,HF_volume2color),
             #     gp = gpar(fill = colorRampPalette(brewer.pal(n = 9, name = "Greys"))(100))
             #   ),
             #   LV_quotient = anno_barplot(
             #     bar_width = 1,
             #     height = unit(40,"points"),
             #     c(RR_quotient,HF_quotient),
             #     gp = gpar(fill = colorRampPalette(brewer.pal(n = 9, name = "Purples"))(100))
             #   )
             # ),
             # 
             
             show_column_dend = TRUE,
             show_row_dend = TRUE,
             row_names_gp = gpar(fontsize = 15, fontfamily = "Helvetica"),
             
             col = colorRamp2(seq(0, max_bar, length = 2),c("#EDF8E9", "#005A32")),
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
             column_title_gp = gpar(fontsize = 30, fontface = "bold", fontfamily = "Helvetica"),
             
             row_title = "\nLead design\n",
             row_title_side = "right",
             row_title_gp = gpar(fontsize = 25, fontfamily = "Helvetica"),
             
             width = unit(0.67, "npc"),
             height = unit(0.7, "npc"),
             row_dend_width = unit(0.1,"npc"),
             column_dend_height = unit(0.2, "npc")
  )
  
  # col_fun = colorRamp2(seq(0, 20, length = 2),c("#EDF8E9", "#005A32"))
  # # lgd = Legend(col_fun = col_fun, title = "AT reduction (%)")
  # # draw(lgd, x = unit(0.5, "npc"), y = unit(0.8, "npc"))
  # 
  # lgd = Legend(col_fun = col_fun, title = "foo", at = c(0, 0.25, 0.5, 0.75, 1))
  draw(hm)
  # return(normalised_ATs)
}
if(draw_hm == FALSE){
  # 
  # num_patients = c()
  # original_patients = ncol_prev = ncol(normalised_ATs)
  # patients2remove <- c()
  # leads2remove <- c()
  # rm(sub_matrix)
  # patients2remove <- c()
  # while(sum(sum(normalised_ATs == -100)) < (ncol(normalised_ATs) * nrow(normalised_ATs))){
  #   normalised_ATs[is.na(normalised_ATs)] <- -100
  #   
  #   max_vec <- apply(normalised_ATs,2,'max') # Vector with the max value for each patient
  #   max_names <- apply(normalised_ATs,2,'which.is.max') # Vec. with idx of the lead achieving each of that value
  #   #max_names <- which(normalised_ATs >= 1*max(max_vec))
  #   
  #   first2die <- names(which(max_vec == max(max_vec))) # Which patient does have that maximum?
  #   leads2remove <- c(leads2remove,(max_names[first2die])) # Indices
  #   sub_matrix <- normalised_ATs[leads2remove[length(leads2remove)],]
  #   patients2remove_aux <- names(which(sub_matrix >= (thresh_plot/100)*max_vec & (sub_matrix > -100)))
  #   
  #   #patients2remove_aux = names(which(max_names == leads2remove[length(leads2remove)])) # Names of the patients achieving the maximum in the same leads
  #   patients2remove <- c(patients2remove,which(names(max_names) %in% patients2remove_aux))
  #   patients2remove <- unique(patients2remove)
  #   
  #   normalised_ATs[,patients2remove] <- -100
  #   normalised_ATs[leads2remove,] <- -100
  #   
  #   num_patients = unique(c(num_patients,length(patients2remove)))
  #   
  #   normalised_ATs[normalised_ATs < -10] <- NA
  # }
  # num_patients
  num_patients <- c()
  names_best_leads <- c()


  while(sum(sum(!is.na(normalised_ATs))) > 0){

  patient_quantiles <- apply(normalised_ATs,2,function(x) quantile(x[!is.na(x)], probs = thresh_plot/100.))
  
  binary_df <- apply(normalised_ATs,1,function(x) x >= patient_quantiles)
  binary_df <- t(binary_df)
  
  patients_improved <- apply(binary_df,1,function(x) sum(x[!is.na(x)])) # For each lead, the number of patients over the percentile specified before
  
  lead2remove <- which.max(patients_improved)
  names_best_leads[length(names_best_leads)+1] <- rownames(normalised_ATs)[lead2remove]
  normalised_ATs[lead2remove,] <- NA # Remove the best lead
  
  patients2remove <- names(which(binary_df[lead2remove,]))
  normalised_ATs[,patients2remove] <- NA #Remove the patients
  
  num_patients <- c(num_patients,max(patients_improved))

  }
  
  num_patients <- cumsum(num_patients)
  # if(version == 4){
  #   names_best_leads <- ascii2bin_lead(names_best_leads)
  # }
  
  if(output == "number")
    return(num_patients)
  else if(output == "names")
    return(names_best_leads)
  else if(output == "both")
    return(list(names_best_leads,num_patients))
}
}


Find_optimal_monodipoles <- function(vein="LA",response=100,which_cases="HF",version, output = "number", SA_folder = SA_folder, RV_midseptum = RV_midseptum){

  y<-hier_clust_multipole(lead_option = vein, thresh_plot = response, which_cases = which_cases,version=version, draw_hm = FALSE, output = output, SA_folder = SA_folder, RV_midseptum = RV_midseptum)
  
  return(y)
  
}

merge_design <- function(s_vec){
  
  s_res <- s_vec[1]
  
  for(i_char in c(1:nchar(s_res))){
    for(i_string in c(1:length(s_vec))){
      if(substr(s_vec[i_string],i_char,i_char) == '1'){
        substr(s_res,i_char,i_char) <- '1'
      }    
    }
  }
  
  return(s_res)
}

check_quadripole_includes_dipole <- function(quad,dip){
  
  quad_ascii <- bin2ascii_lead(quad)
  dip_ascii <- bin2ascii_lead(dip)
  
  quad_vec <- strsplit(quad_ascii,"")[[1]]
  dip_vec <- strsplit(dip_ascii,"")[[1]]
  
  count <- 0
  
  for(i in c(1:length(dip_vec))){
    for(j in c(1:length(quad_vec))){
      if(dip_vec[i] == quad_vec[j]){
        count <- count+1
      }
    }
  }
  
  if(count == length(dip_vec)){
    return(2)
  }
  else if(count == 0){
    return(0)
  }
  else{
    return(1)
  }
}

check_vector_includes_dipole <- function(vec,dip){
  result = 0
  for(i in c(1:length(vec))){
    if(check_quadripole_includes_dipole(vec[i],dip) == 2){
      result = 2
      break
    }
    else if(check_quadripole_includes_dipole(vec[i],dip) == 1){
      result = 1
    }
  }
  return(result)
}

sum_string <- function(s){
  sol <- c()
  if(length(s) > 1){
    for(i in c(1:length(s))){
      elem_s <- sum_string(s[i])
      sol <- c(sol,elem_s)
    }
  }
  else{
    sol <- sum(as.integer(strsplit(s,"")[[1]]))
  }
  
  return(sol)
}

bin2ascii_lead <- function(s_bin){
  
  s_ascii <- c()
  
  for(element in c(1:length(s_bin))){
    if(nchar(s_bin[element]) < 8){
      s_bin[element] <- paste(c(rep("0",8-nchar(s_bin[element])),s_bin[element]),collapse = "")
    }
    s_ascii_int <- c()
    for(character in c(1:nchar(s_bin[element]))){
      if(substr(s_bin[element],character,character) == "1")
        s_ascii_int[length(s_ascii_int) + 1] <- character
    }
    s_ascii[element] <- paste(s_ascii_int,collapse="")
  }
  
  return(s_ascii)
}

ascii2bin_lead <- function(s_ascii_supervec){
  
  s_bin_supervec <- c()
  
  for(s_ascii in s_ascii_supervec){
    s_ascii <- toString(s_ascii)
    s_vec <- as.integer(strsplit(s_ascii,"")[[1]])
    
    s_bin_vec <- c(rep("0",8))
    s_bin_vec[s_vec] <- "1"
    
    s_bin_supervec[length(s_bin_supervec) + 1] <- paste(s_bin_vec,collapse = "")
  }
  
  return(s_bin_supervec)
}

complete_quadripole <- function(almost_quadripole, bipolar_names){
  
  i <- 1
  while((sum_string(almost_quadripole) != 4) && (i < length(bipolar_names))){
    final_quadripole <- merge_design(c(almost_quadripole,bipolar_names[i]))
    if(sum_string(final_quadripole) < 4){
      almost_quadripole <- final_quadripole
      i <- i+1
    }
    else if(sum_string(final_quadripole) > 4){
      i <- i+1
    }
    else if(sum_string(final_quadripole) == 4)
      almost_quadripole <- final_quadripole
  }
  
  if(i == length(bipolar_names)){ # We just put it as spaced as possible
      quad_ints <- as.integer(strsplit(bin2ascii_lead(almost_quadripole),"")[[1]])
      dist_vec <- c()
    for(num in c(1:8)){
      dist_vec[length(dist_vec) + 1] <- min(dist(c(num,quad_ints))[1:length(quad_ints)])
    }
      almost_quadripole <- merge_design(c(almost_quadripole,ascii2bin_lead(which.max(dist_vec))))
  }
  
  if(sum_string(almost_quadripole) < 4){
    almost_quadripole <- complete_quadripole(almost_quadripole,bipolar_names)
  }
  
  return(almost_quadripole)
}

find_optimal_quadripole <- function(vein,which_cases,method){
  
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
    hierarchy_veins <- c("LA","PL","AL","PO","AN")
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

find_optimal_quadripole_combinatorial <- function(vein,which_cases,method){
  
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
    hierarchy_veins <- c("LA","PL","AL","PO","AN")
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

Find_Optimal_Quadripole_Evenly <- function(vein,which_cases,method){
  
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
    hierarchy_veins <- c("LA","PL","AL","PO","AN")
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

# Function to get the second or Nth maximum of an array
maxN <- function(x, N=2){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

Find_Optimal_Quadripole_Optimising <- function(vein,which_cases,method,response=100,SA_folder,version, RV_midseptum = RV_midseptum){
  quadripole <- c()
  if(vein != "ALL"){
    # We get the names of the optimal lead designs in monopoles and dipoles
    bipolar_names_score <- Find_optimal_monodipoles(vein = vein, response = response, which_cases = which_cases, version = version, output = "both", SA_folder = SA_folder, RV_midseptum = RV_midseptum)
    
    bipolar_names <- bipolar_names_score[[1]]
    bipolar_scores <- bipolar_names_score[[2]]
    
    #We get the improvement for each design
    bipolar_scores <- c(bipolar_scores[1], bipolar_scores[-1] - bipolar_scores[-length(bipolar_scores)])
    
    
  }
  else{
    hierarchy_veins <- c("LA","PL","AL","PO","AN")
    bipolar_names <- c()
    bipolar_scores <- c()
    for(each_vein in hierarchy_veins){
      # print(paste0("Checking vein ",each_vein,"..."))
      # We get the names of the optimal lead designs in monopoles and dipoles
      bipolar_names_score_temp <- Find_optimal_monodipoles(vein = each_vein, response = response, which_cases = which_cases, version = version, output = "both", SA_folder = SA_folder, RV_midseptum = RV_midseptum)
      
      bipolar_names_temp <- bipolar_names_score_temp[[1]]
      bipolar_scores_temp <- bipolar_names_score_temp[[2]]
      # print(bin2ascii_lead(bipolar_names_temp))
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
  # print("Checking the combinatorial...")
  bipolar_names_ascii <- bin2ascii_lead(bipolar_names)
  # We create all the combinations of 4 electrodes from 1 to 8
  matrix_comb <- t(combn(1:8,4))
  
  vec_comb <- apply(matrix_comb,1,function(x) paste0(x,collapse = ''))
  
  quad_scores <- 0*c(1:length(vec_comb))
  matrix_quad_includes_bip <- matrix(data = 0,nrow = length(bipolar_names), ncol = length(vec_comb))
  
  # We create a vector with the scores of all the quadripoles
  for(i in seq(1:length(vec_comb))){
    for(j in seq(1:length(bipolar_names))){
      if(check_quadripole_includes_dipole(ascii2bin_lead(vec_comb[i]),bipolar_names[j]) == 2){
        quad_scores[i] = quad_scores[i] + bipolar_scores[j]
        matrix_quad_includes_bip[j,i] <- 1
      }
    }
  }
  
  q <- 1
  possible_solutions <- c()
  
  # First iteration
  
  permutations_q_ciphers <- t(combn(nrow(matrix_comb),1)) # Indices to check all the quadpoles
  for(i in c(1:nrow(permutations_q_ciphers))){
    #All the bipoles will be included if merging all the binary vectors corresponding to that
    #quadpole have at least a one in all their positions
    flag_to_break <- matrix_quad_includes_bip[,permutations_q_ciphers[i]] 
    if(prod(flag_to_break) != 0){
      # We add a row to the possible solutions
      possible_solutions <- c(possible_solutions,vec_comb[permutations_q_ciphers[i]])
    }
  }
  q <- q+1
  
  if(length(possible_solutions) != 0){
    possible_solutions <- as.data.frame(possible_solutions)
  }

  # print("Choosing the minimal amount of quadpoles needed...")
  # We search for the minimal amount of quadripolar designs that include all the bipolar designs
  while(length(possible_solutions) == 0){ #We'll stop once we have a solution
    # print(paste0("Checking if with ",q))
    permutations_q_ciphers <- t(combn(nrow(matrix_comb),q)) # Indices to check all the quadpoles
    for(i in c(1:nrow(permutations_q_ciphers))){
      #All the bipoles will be included if merging all the binary vectors corresponding to that
      #quadpole have at least a one in all their positions
      flag_to_break <- matrix_quad_includes_bip[,permutations_q_ciphers[i,1]] 
      for (j in c(2:ncol(permutations_q_ciphers))){
        flag_to_break <- flag_to_break + matrix_quad_includes_bip[,permutations_q_ciphers[i,j]]
      }
      if(prod(flag_to_break) != 0){
        # We add a row to the possible solutions
        possible_solutions <- rbind(possible_solutions,vec_comb[permutations_q_ciphers[i,]])
      }
    }
    q <- q+1
  }
  
  
  # Among all the possible solutions we choose the one that maximise the benefits.
  
  score_possible_solutions <- matrix(nrow=nrow(possible_solutions),ncol = ncol(possible_solutions))
  
  # Could be optimised but meh
  
  for(i in c(1:nrow(score_possible_solutions))){
    for (j in c(1:ncol(score_possible_solutions))){
      score_possible_solutions[i,j] <- quad_scores[match(possible_solutions[i,j],vec_comb)]
    }
  }
  
  cummulative_score <- rowSums(score_possible_solutions)
  # print("Choosing the optimal option...")
  # We choose the quadripole(s) with the maximum score
  max_idx <- which(cummulative_score == max(cummulative_score))
  if(q > 2){
  if(length(max_idx) > 1){ # In case of a draw between several multipolar designs we choose the one with the maximum benefit on the first one. If still in draw we continue.
    possible_solutions <- possible_solutions[max_idx,]
    
    
    for(i in c(1:nrow(possible_solutions))){
      possible_solutions[i,] <- possible_solutions[i,order(score_possible_solutions[i,],decreasing = TRUE)]
    }
    
    score_possible_solutions <- score_possible_solutions[max_idx,]
    score_possible_solutions <- t(apply(score_possible_solutions, 1, function(x) sort(x,decreasing = TRUE)))
    
    flag_to_break <- FALSE
    j <- 1
  
    while(!flag_to_break){
      max_idx <- which(score_possible_solutions[,j] == max(score_possible_solutions[,j]))
      if(length(max_idx) == 1){
        flag_to_break <- TRUE
      }
      else{
        j <- j+1
        if(j == ncol(score_possible_solutions)){
          flag_to_break = TRUE
          max_idx = 1
        }
      }
    }
  }
  }
  
  return(list("Patient improvement"=bipolar_scores,"Bipolar"=bipolar_names_ascii,"Quadripolar"=as.vector(possible_solutions[max_idx[1],]),"Quadripolar_Scores"=score_possible_solutions[max_idx[1],]))
}

isbipoleinquadpole <- function(bipole,quadpole){
  
  bipole <- toString(bipole)
  quadpole <- toString(quadpole)
  
  bip_digits<-strsplit(bipole,"")[[1]]
  quad_digits<-strsplit(quadpole,"")[[1]]
  
  return(prod(bip_digits %in% quad_digits))
}

Plot_Dipoles_Fancy <- function(which_cases = "HF", with_lines = TRUE, bipole_from_file = FALSE,method = "max",savefile=FALSE, response = 100,sort_lines = c(), SA_folder, RV_midseptum = RV_midseptum){
  
  require("scales")
  
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
    for (vein in c("AN","AL","LA","PL","PO")) {
      res_list <- Find_optimal_monodipoles(vein = vein, response = response, which_cases = which_cases, version=4, output = "both", SA_folder = SA_folder, RV_midseptum = RV_midseptum)
      
      temp_df <- as.data.frame(matrix(ncol=3,nrow = length(res_list[[1]])))
      colnames(temp_df) <- c("Vein","Patients","Bipole")
      
      temp_df$Vein <- vein
  
      # We get the increment of patients, not the total
      bipolar_scores <- as.vector(res_list[[2]])
      bipolar_scores <- c(bipolar_scores[1], bipolar_scores[-1] - bipolar_scores[-length(bipolar_scores)])
      
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
      res_list <- Find_Optimal_Quadripole_Optimising("ALL",which_cases,method,response=response,SA_folder=SA_folder,version=4, RV_midseptum = RV_midseptum)
      # print(res_list)

      optimal_quadpoles <- res_list$Quadripolar
  }
  # print("Setting up the plot...")
  
  # We sort the bipoles (alphabetically)
  all_bipoles <-  sort(unique(optimal_choices$Bipole))
  
  # We create fake bars for when there's no bipole, setting a height of 0.1
  for(vein_name in c("AN","AL","LA","PL","PO")){
  # Extract the data frame corresponding to each vein
  vein_vec <- optimal_choices[optimal_choices$Vein == vein_name,]
  # Extract the bipoles that are not included in this subset
  which_to_add <- all_bipoles[!all_bipoles %in% vein_vec$Bipole]
  # Create dummy cases with height 0.1
  df_to_add <- data.frame(rep(vein_name,length(which_to_add)),rep(0.1,length(which_to_add)),which_to_add)
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
  
  title_str <- paste0(title_str," grouped by ",method, "\n in the ", response, "% window")
  
  
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
  
  perc_quadpoles <- round(100*as.vector(res_list$Quadripolar_Scores)/num_patients,2)
  perc_bipoles <- perc_quadpoles/5
  # We create and auxiliar vector for the text of the legend
  legend_text <- paste0(optimal_quadpoles," (",perc_bipoles,"%)")
  
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
        # If the highest value is the 5th quadpole, we want the index 1, the indices vector will be c(5,...)
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
    line_df$Quadripole[i] <- paste0(corresponding_quad," (",corresponding_perc,"%)")
  }
  
  
  line_df$Quadripole <- factor(line_df$Quadripole, levels = legend_text[indices])
  }

# We sort the categories so they appear in order in the histogram
optimal_choices$Vein <- factor(optimal_choices$Vein,levels = c("AN","AL","LA","PL","PO"))
  
  toplot = ggplot(optimal_choices) + 
    geom_bar(aes(Bipole, Patients, fill=Vein), stat="identity",position='dodge') 

  if(with_lines){
    toplot = toplot + geom_segment(data = line_df, aes(x = x_ini, xend = x_end, y = y, yend = y, group = row, color = Quadripole), size = 2) +
    scale_color_manual(values = brewer.pal(length(optimal_quadpoles),"Paired"))
      # scale_color_manual(values = c("green","orange","red"))
  }
  # print(num_patients)
  toplot = toplot +theme_classic(base_size = 40) + theme(plot.title = element_text(hjust = 0.5))  + 
    scale_y_continuous(breaks = seq(0,num_patients,5), limits = c(min(line_df$y),num_patients), oob = rescale_none) + ggtitle(title_str) + xlab("Electrodes activated")

  
  print(toplot)
  
  if(savefile)
    ggsave(filename=paste0("/home/crg17/Pictures/",which_cases,"_response_",toString(response),"_",SA_folder,".png"),width = 20,height = 10)
  
  # return(unique(line_df$Quadripole))
}

Correct_HFLA_leads <- function(foldername,casenumber){
  
  multipole_AT1090 <- read.csv(paste0("/data/SA_multipole/",foldername,"/HF/",casenumber,"/multipole_AT1090.dat"), sep="", stringsAsFactors=FALSE)
  
  multipole_AT1090$lead <- ascii2bin_lead(bin2ascii_lead(multipole_AT1090$lead))
  
  multipole_AT1090_corrected <- multipole_AT1090
  
  for(i in c(1:nrow(multipole_AT1090_corrected))){
    multipole_AT1090_corrected[i,"LA"] <- multipole_AT1090[which(multipole_AT1090$lead == intToUtf8(rev(utf8ToInt(multipole_AT1090$lead[i])))),"LA"]
  }
  
  write.table(multipole_AT1090_corrected,file=paste0("/data/SA_multipole/",foldername,"/HF/",casenumber,"/multipole_AT1090_leadcorrected.dat"),quote = FALSE,sep = " ",dec = ".",row.names = FALSE,col.names = TRUE)
  
  
}

CreateMonopoles <- function(){
  require(jjb)
  
  for(heart in c(paste0("0",c(1:9)),10:20)){
      AT_apex <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart,"/simulations/multipole/eikonal_midseptum/BiV.midseptum/vm_act_seq.dat"))
    for(vein in c("AN","AL","LA","PL","PO")){
      for(electrode in c(1:8)){
      AT <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart,"/simulations/multipole/eikonal/",vein,"_",electrode,"/vm_act_seq.dat"))
      
      AT_res <- pmin(AT_apex,AT)
      
      directory <- paste0("/data/SA_multipole/midseptum/h/",heart)
      outname <- paste0("h",heart,"_",vein,"_",ascii2bin_lead(electrode),".dat")
      
      mkdir(directory, r = TRUE)
      
      write.table(AT_res,file=paste0(directory,"/",outname),quote = FALSE,row.names = FALSE,col.names = FALSE)
      }
    }
  }
  
  for(heart in c(paste0("0",c(1:9)),10:24)){
    AT_apex <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",heart,"/simulations/multipole/eikonal_midseptum/BiV.midseptum/vm_act_seq.dat"))
    for(vein in c("AN","AL","LA","PL","PO")){
      for(electrode in c(1:8)){
        AT <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",heart,"/simulations/multipole/eikonal/",vein,"_",electrode,"/vm_act_seq.dat"))
        
        AT_res <- pmin(AT_apex,AT)
        
        directory <- paste0("/data/SA_multipole/midseptum/HF/",heart)
        outname <- paste0("h",heart,"_",vein,"_",ascii2bin_lead(electrode),".dat")
        
        mkdir(directory, r = TRUE)
        
        write.table(AT_res,file=paste0(directory,"/",outname),quote = FALSE,row.names = FALSE,col.names = FALSE)
      }
    }
  }
}

CreateDipoles <- function(){
  require(jjb)
  
  for(heart in c(paste0("0",c(1:9)),10:20)){
    AT_apex <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart,"/simulations/multipole/eikonal_midseptum/BiV.midseptum/vm_act_seq.dat"))
    for(vein in c("AN","AL","LA","PL","PO")){
      for(electrode in c(1:8)){
        AT <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/h_case",heart,"/simulations/multipole/eikonal/",vein,"_",electrode,"/vm_act_seq.dat"))
        
        AT_res <- pmin(AT_apex,AT)
        
        directory <- paste0("/data/SA_multipole/midseptum/h/",heart)
        outname <- paste0("h",heart,"_",vein,"_",ascii2bin_lead(electrode),".dat")
        
        mkdir(directory, r = TRUE)
        
        write.table(AT_res,file=paste0(directory,"/",outname),quote = FALSE,row.names = FALSE,col.names = FALSE)
      }
    }
  }
  
  for(heart in c(paste0("0",c(1:9)),10:24)){
    AT_apex <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",heart,"/simulations/multipole/eikonal_midseptum/BiV.midseptum/vm_act_seq.dat"))
    for(vein in c("AN","AL","LA","PL","PO")){
      for(electrode in c(1:8)){
        AT <- read.table(paste0("/media/crg17/Seagate Backup Plus Drive/CT_cases/HF_case",heart,"/simulations/multipole/eikonal/",vein,"_",electrode,"/vm_act_seq.dat"))
        
        AT_res <- pmin(AT_apex,AT)
        
        directory <- paste0("/data/SA_multipole/midseptum/HF/",heart)
        outname <- paste0("h",heart,"_",vein,"_",ascii2bin_lead(electrode),".dat")
        
        mkdir(directory, r = TRUE)
        
        write.table(AT_res,file=paste0(directory,"/",outname),quote = FALSE,row.names = FALSE,col.names = FALSE)
      }
    }
  }
}