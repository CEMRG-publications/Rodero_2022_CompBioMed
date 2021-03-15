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

isbipoleinquadpole <- function(bipole,quadpole){
  
  bipole <- toString(bipole)
  quadpole <- toString(quadpole)
  
  bip_digits<-strsplit(bipole,"")[[1]]
  quad_digits<-strsplit(quadpole,"")[[1]]
  
  return(prod(bip_digits %in% quad_digits))
}