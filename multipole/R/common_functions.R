#' @description Loads the needed packages and install them if they are not
#' downloaded.
#' @param list.of.packages Vector with the strings of the name of the packages.
Load_Install_Packages <- function(list.of.packages){
  new.packages <- list.of.packages[!(list.of.packages %in% 
                                       installed.packages()[,"Package"])]
  if(length(new.packages)){ 
      install.packages(new.packages)
    }
    flags<-lapply(list.of.packages, require, character.only = TRUE)
    
    return(flags)
}


#' @description Function to read the dataframes of the AT appropriately. Before
#' it was read_dataframes. Last version is version 4, but just in case I'm not
#' deleting previous versions.
#' 
#' @param metric_option "QRS", "AT1090", "VEUTAT" or "VEUmean" to read the 
#' corresponding file. For the case of the VEUs, if normalised, they will be 
#' read on absolute value.
#' @param lead_option Subselects that vein. For version 4 posibilies are "AN",
#' "AL", "LA", "IL", "IN".
#' @param which_cases "RR", "HF" or "RRHF" for reverse remodelled, heart failure
#' or both respectively.
#' @param version Version of the function to read the dataframe. Latest is 4.
#' @param SA_folder Sensitivity analysis folder from which to read the files.
#' @param normalising If TRUE, it normalises the AT over the times of the RV
#' apex.
#' @param output For version = 4, this value can be % or ms for just difference
#' of ms w.r.t. baseline (if output == ms) or normalised (if output == %)
#' @param flag_debugging If TRUE, it prints whatever is reading and writing.
#' 
#' @return The dataframe which columns are each one of the hearts and the rows
#' are each one of the lead configurations.

Read_dataframes <- function(metric_option, lead_option, which_cases, version,
                            SA_folder, normalising = TRUE,
                            output = "%", flag_debugging = FALSE,
                            root_directory = "/media/crg17/Seagate Expansion Drive/SA_multipole/"){
  
  source("/home/crg17/Desktop/scripts/multipole/R/leads_operations.R")
  Load_Install_Packages("dplyr")
  
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
    
    if(which_cases == "RR" || which_cases == "RRHF"){
      for(heart in c("01","02","03","04","05","06","07","08","09",10:20)){
        singleheart <- Read_table(paste0(root_directory,SA_folder,"/h/",
                                         heart,"/multipole_",metric_option,
                                         ".dat"), header = TRUE,
                                  flag_debugging = flag_debugging)
        singlerotlead <- singleheart[,lead_option]
        
        RR_cases <- cbind(RR_cases,singlerotlead)    
      }
        RR_RVapex <- Read_table(paste0(root_directory,SA_folder,
                                       "/h/multipole_RVapex.dat"),
                                header = TRUE, flag_debugging = flag_debugging)
      
      # Normalising
      
      AT_RR <- as.data.frame(RR_cases)
      
      if(substr(metric_option,1,3) != "VEU"){
        # prenorm_RR <- mapply('/', AT_RR, RR_RVapex[,metric_option])
        # normalised_RR <-  100*(1-abs(prenorm_RR))
        prenorm_RR <- mapply('-',RR_RVapex[,metric_option],AT_RR)
        if(output == "ms"){
          normalised_RR <- prenorm_RR
        }
        else if(output == "%"){
          normalised_RR <- 100*t(apply(prenorm_RR,1,
                                       function(x) x/(RR_RVapex[,metric_option])
                                       ))
        }
      } else{
        normalised_RR <- mapply('-', AT_RR, RR_RVapex[,metric_option])
      }
      
      colnames(normalised_RR) <- c("RR01","RR02","RR03","RR04","RR05","RR06",
                                   "RR07","RR08","RR09","RR10","RR11","RR12",
                                   "RR13","RR14","RR15","RR16","RR17","RR18",
                                   "RR19","RR20")
      rownames(normalised_RR) <- singleheart$lead %>% 
                                 bin2ascii_lead() %>% ascii2bin_lead()
                        
      
      colnames(AT_RR) <- colnames(normalised_RR)
      rownames(AT_RR) <- rownames(normalised_RR)
    }
    
    if(which_cases == "HF" || which_cases == "RRHF"){

      HF_cases <- c()
      hearts_vec <- c("01","02","03","04","05","06","07","08","09",10:24)
      
      if(SA_folder == "scar" || SA_folder == "scar_6mm"){
        hearts_vec <- hearts_vec[-c(13,21)]
      }
      for(heart in hearts_vec){
          singleheart <- paste0(root_directory,SA_folder,"/HF/",heart,
                                "/multipole_",metric_option,".dat") %>%
            Read_table(., header = TRUE, flag_debugging = flag_debugging)
        
        singlerotlead <- singleheart[,lead_option]
        
        
        HF_cases <- cbind(HF_cases,singlerotlead)    
      }
      
        HF_RVapex <- paste0(root_directory,SA_folder,
                            "/HF/multipole_RVapex.dat") %>%
          Read_table(., header=TRUE, flag_debugging = flag_debugging)
      
      # Normalising
      AT_HF <- as.data.frame(HF_cases)
      if(substr(metric_option,1,3) != "VEU"){
        # prenorm_HF <- mapply('/', AT_HF, HF_RVapex[,metric_option])
        # normalised_HF <- 100*(1-abs(prenorm_HF))
        prenorm_HF <- mapply('-',HF_RVapex[,metric_option],AT_HF)
        if(output == "ms"){
          normalised_HF <- prenorm_HF
        }
        else if(output == "%"){
          normalised_HF <- 100*t(apply(prenorm_HF,1,
                                       function(x) x/(HF_RVapex[,metric_option])
          ))
        }
      } else{
        normalised_HF <- mapply('-', AT_HF, HF_RVapex[,metric_option])
      }
      
      colnames(normalised_HF) <- paste0("HF",hearts_vec)
      rownames(normalised_HF) <- singleheart$lead %>%
                                 bin2ascii_lead() %>% ascii2bin_lead()

      colnames(AT_HF) <- colnames(normalised_HF)
      rownames(AT_HF) <- rownames(normalised_HF)
      
    }
    
  }
  
  if(which_cases == "RRHF"){
    if(normalising)
      return(list("normalised_RR" = normalised_RR, "normalised_HF" = normalised_HF))
    else
      return(list("normalised_RR" = AT_RR, "normalised_HF" = AT_HF))
  }
  else if(which_cases == "RR"){
    if(normalising)
      return(list("normalised_RR" = normalised_RR))
    else
      return(list("normalised_RR" = AT_RR))
  }
  else if(which_cases == "HF"){
    if(normalising)
      return(list("normalised_HF" = normalised_HF))
    else
      return(list("normalised_HF" = AT_HF))
  }
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

#' @description Function to break from a function without giving an error
#' message.
Stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

#' @description It prints a message of the file reading or writing.
#' @param string Directory and name of the file.
#' @param io "r" if reading the file, "w" if writing it.
#' @param flag_debugging If FALSE it does not run anything.
Debug_message <- function(string,io,flag_debugging=TRUE){
  if(flag_debugging){
  Load_Install_Packages(c("dplyr"))
  
  if(io == "r"){
    paste0("You need the file ",string) %>%
      print(.)
  }
  else if(io == "w"){
    paste0("You are getting the file ",string) %>%
      print(.)
  }
  }
}

#' @description Sends me an email with the subject "Process done" and the input
#' message body.
#' 
#' @param string Body of the e-mail.
Send_mail <- function(string){
  system(paste0("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"",
                string,"\""))
}

#' @description Same as read.csv function with the Debug_message function added.
#' 
#' @param ... Same as the ones in read.csv
#' @param flag_debugging If TRUE, prints whatever is reading.
#' 
#' @return The file read.
Read_csv <- function(file, header = TRUE, sep = ",", quote = "\"",dec = ".",
                     fill = TRUE, comment.char = "", flag_debugging = FALSE,
                     stringsAsFactors = default.stringsAsFactors()){
  source("/home/crg17/Desktop/scripts/4chmodel/R/common_functions.R")
  
  Debug_message(file,"r",flag_debugging = flag_debugging)
  
  file_read <- read.csv(file = file, header = header, sep = sep, quote = quote,
                        dec = dec, fill = fill, comment.char = comment.char,
                        stringsAsFactors = stringsAsFactors)
  return(file_read)
}

#' @description Same as read.table function with the Debug_message function
#'  added.
#' 
#' @param ... Same as the ones in read.table
#' @param flag_debugging If TRUE, prints whatever is reading.
#' 
#' @return The file read.

Read_table <- function(file, header = FALSE, sep = "", quote = "\"'", dec = ".",
                       numerals = c("allow.loss", "warn.loss", "no.loss"),
                       row.names, col.names, as.is = !stringsAsFactors,
                       na.strings = "NA", colClasses = NA, nrows = -1, skip = 0,
                       check.names = TRUE, fill = !blank.lines.skip, 
                       strip.white = FALSE, blank.lines.skip = TRUE, 
                       comment.char = "#", allowEscapes = FALSE, flush = FALSE,
                       stringsAsFactors = default.stringsAsFactors(),
                       fileEncoding = "", encoding = "unknown", text,
                       skipNul = FALSE, flag_debugging = FALSE){
  
  Debug_message(file,"r",flag_debugging = flag_debugging)
  
  file_read <- read.table(file = file, header = header, sep = sep, 
                          quote = quote, dec = dec, numerals = numerals,
                          row.names = row.names, col.names = col.names,
                          as.is = as.is, na.strings = na.strings,
                          colClasses = colClasses, nrows = nrows, skip = skip,
                          check.names = check.names, fill = fill, 
                          strip.white = strip.white,
                          blank.lines.skip = blank.lines.skip, 
                          comment.char = comment.char,
                          allowEscapes = allowEscapes, flush = flush, 
                          stringsAsFactors = stringsAsFactors,
                          fileEncoding = fileEncoding, encoding = encoding,
                          text = text, skipNul = skipNul)
  return(file_read)
}

#' @description Same as write.table function with the Debug_message function
#'  added.
#' 
#' @param ... Same as the ones in write.table
#' @param flag_debugging If TRUE, prints whatever is writing.

Write_table <- function(x, file = "", append = FALSE, quote = TRUE, sep = " ",
                        eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                        col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "", flag_debugging = FALSE){
  
  Debug_message(file,"w",flag_debugging = flag_debugging)
  
  write.table(x = x, file = file, append = append, quote = quote, sep = sep,
              eol = eol, na = na, dec = dec, row.names = row.names, 
              col.names = col.names, qmethod = qmethod,
              fileEncoding = fileEncoding)
}
  

#' @description Extract the name of a variable as a string. Addapted from
#' StackOverflow.
#' 
#' @param variable The variable you want to extract the name from. I think any
#' type of variable is accepted.
#' 
#' @return The name of the variable as a string.
Extract_variable_name <- function(variable){
  variable_symbol <- substitute(variable)
  variable_string <- sapply(variable_symbol,deparse)
  
  return(variable_string)
}

Find_first_common_idx <- function(u,v,i = 1){
  
  element <- u[which(u[1:i]%in%v[1:i])]
  
  if(length(element) == 0){
    element <- Find_first_common_idx(u, v, i = i+1)
  }
    return(element)
  
}
