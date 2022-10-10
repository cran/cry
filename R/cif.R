#
# This file is part of the cry package
#
# Functions connected to reflections data.

#' Reads and output a CIF file
#'
#' @param filename A character string. The path to a valid CIF file.
#' @param message A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the cif file.
#' @return A named list. Each name correspond to a valid field in the cif.
#'
#' @examples
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"AMS_DATA.cif")
#' lCIF <- readCIF(filename)
#' print(names(lCIF))
#' print(lCIF$INTRO$CELL)
#' print(lCIF$INTRO$HALL)
#' print(lCIF$INTRO$HM)
#' print(lCIF$SYMM)
#'
#' @export
readCIF <- function(filename, message=FALSE){
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  l_list <- grep("loop_",lcif)
  l_list1 <- append(l_list,length(lcif))
  mat<-zoo::rollapply(l_list1, 2,stack)
  ch <- apply(mat, 1, function(x) lcif[(x[1]+1):(x[2]-1)])
  crystal_summary <- lapply(ch, ucheck, pattern="_publ_author_name")
  symmetry <- lapply(ch, ucheck, pattern="_space_group_symop_operation")
  reflection <- lapply(ch, ucheck, pattern="_refln_index_h")
  coordinates <- lapply(ch, ucheck, pattern="_atom_site_fract_x")
  anisotropy <- lapply(ch, ucheck, pattern="_atom_site_aniso_label")
  symbol <- lapply(ch, ucheck, pattern="_atom_type_symbol")
  geom_angle <- lapply(ch, ucheck, pattern="geom_angle_atom_site_label_1")
  geom_distance <- lapply(ch, ucheck, pattern="_geom_bond_atom_site_label_1")
  geom_hbond <- lapply(ch, ucheck, pattern="geom_hbond_atom_site_label_D")
  geom_torsion <- lapply(ch, ucheck, pattern="geom_torsion_atom_site_label_1")

  intro <- r_summ(lcif)
  symm <- ansnull(symmetry)
  reflections1 <- r_reflections1(lcif)
  reflections2 <- if (is.na(nanona(reflection)) == FALSE) r_reflections2(nanona(reflection)) else NA
  reflections <- r_reflcs(reflections1,reflections2)
  coordinate <- if (is.na(nanona(coordinates)) == FALSE) clean(r_spositions(nanona(coordinates))) else NULL
  #coordinate <- r_positions(ansnull(coordinates))
  anisotropies <- if (is.na(nanona(anisotropy)) == FALSE) clean(r_saniso(nanona(anisotropy))) else NULL
  #anisotropies <- r_aniso(ansnull(anisotropy))
  symbolics <- ansnull(symbol)
  angle <- if (is.na(nanona(geom_angle)) == FALSE) clean(r_angle(nanona(geom_angle))) else NULL
  #angle <- r_angle(ansnull(geom_angle))
  distance <- if (is.na(nanona(geom_distance)) == FALSE) clean(r_dist(nanona(geom_distance))) else NULL
  #distance <- r_dist(ansnull(geom_distance))
  hbond <- if (is.na(nanona(geom_hbond)) == FALSE) clean(r_hbond(nanona(geom_hbond))) else NULL
  #hbond <- r_hbond(ansnull(geom_hbond))
  torsion <- if (is.na(nanona(geom_torsion)) == FALSE) clean(r_torsion(nanona(geom_torsion))) else NULL
  #torsion <- r_torsion(ansnull(geom_torsion))
  #CIF=list(INTRO=intro,SYMM=symm,COOR=coordinate,ANISO=anisotropies)
  geome = list(ANGLE=angle,DIST=distance,HBOND=hbond,TORSION=torsion)
  CIF = list(HEADER=intro,SYMM=symm,REFL=reflections,COOR=coordinate,ANISO=anisotropies,SYMB=symbolics,GEOM=geome)
  close(f)
  if (message) {
     if (!is.null(reflections)){
     n <- length(reflections$F_squared_meas)
     f <- as.numeric(reflections$F_squared_meas)
	 msg <- c("\n")
	 msg1 <- c(msg,sprintf("File %s read successfully.\n",filename))
     msg2 <- sprintf("There are %d reflections in this file.\n",n)
     msg3 <- c(msg,"Here is a summary of the observations:\n")
     msg4 <- c("\n")
	 out <- c(msg,msg1,msg2,msg3,msg4)
	 cat(out)
	 print(summary(f))
	 } else {
	 anum <- nrow(coordinate$VAL)
	 msg <- c("\n")
	 msg <- c(msg,sprintf("File %s read successfully.\n",filename))
     msg1 <- sprintf("The file does not contain reflection datablock,
	 please refer corresponding reflection file (fcf or hkl).\n")
     msg2 <- sprintf("There are %d atoms in the molecule.\n",anum)
	 out <- c(msg,msg1,msg2)
     cat(out)
	 }
  }
  return(CIF)
}


### accessory functions ####

stack<-function(x){
  j <- c(x[1],x[2])
  return(j)
}

ucheck <- function(x,pattern){
  r <- unlist(x)
  if (length(grep(pattern,r))>0){
    piece <- r
  } else
  { piece <- NA
  return(piece)
  }
}

ansnull <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- x[!is.na(x)]
  return(out)
  }
}

nanona <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NA
  } else
  { out <- x[!is.na(x)]
  return(out)
  }
}

recheck <- function(r1){
  r2 <- gsub("[:):]","",gsub("[:(:]",",",r1))
  r3 <- as.numeric((strsplit(r2, ",")[[1]])[1])
  return(r3)
}

check <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  r2 <- if (length(grep("[:(:]",r1,value = TRUE)>0) == TRUE) recheck(r1) else r1
  return(r2)
}

check1 <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "'")[[1]])[2] else NA
  return(r1)
}

clean1 <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- nc_type(as.data.frame(x))
  return(out)
  }
}


clean <- function(x){
  co1 <- data.frame(gsub ("[()]","",as.matrix(x),perl=T),stringsAsFactors = FALSE)
  ref <- data.frame(gsub("(?<!\\))(?:\\w+|[^()])(?!\\))","",as.matrix(x),perl=T))
  ref1 <- data.frame(gsub("[()]","",as.matrix(ref),perl=T),stringsAsFactors = FALSE)
  ref1[ref1==""]<-NA
  ref2 <- clean1(ref1)
  col1 <- nc_type(co1)
  return(list(VAL=col1,STD=ref2))
}

reap1 <- function(x){
  if (all(is.na(x)) == TRUE){
    out <- NULL
  } else
  { out <- as.numeric(x)
  return(out)
  }
}

reap <- function(pattern,word){
  r <- grep(pattern, word, value = TRUE)
  r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  v <- as.numeric(gsub ("[()]","",as.matrix(r1),perl=T))
  s <- gsub("(?<!\\))(?:\\w+|[^()])(?!\\))","",as.matrix(r1),perl=T)
  s1 <- gsub("[()]","",as.matrix(s),perl=T)
  s2 <- reap1(s1)
  return(list(VAL=v,STD=s2))
}

nc_type <- function(data){
  count <- as.numeric(ncol(data))
  if (isTRUE(count > 2)) {
    data[] <- lapply(data, function(x) numas(x))
    out <- data
  } else if (count == 2){
    l1_data <- list(data$VAL)
    l_1 <- lapply(l1_data[[1]], function(x) numas(x))
    l2_data <- list(data$KEY)
    l2 <- c(gsub("\\[|\\]" ,"",l2_data[[1]]))
    names(l_1) <- c(l2)
    out <- l_1
    return(out)
  }
}

numas <- function(x){
  data <- x
  out <- (suppressWarnings(as.numeric(data)))
  if (all(is.na(out))== FALSE) {
    out1 <- out
  } else {
    out1 <- as.character(data)
  }
  return(out1)
}

r_reflcs <- function(x,y){
  if (all(is.na(x))== FALSE) {
      out <- x
	  } else {
	  if (all(is.na(y)) == FALSE) {
	  out <- y
	  } else {
	   out <- NULL
	  }
	  return(out)
	 }
   }

r_reflections1 <- function (x){
  data <- unlist(x)
  cu <- grep("_shelx_hkl_file", data)
  l <- length(cu)
  if (l > 0) {
     sc <- grep(";", data)
     sc1 <- sc[sc > cu]
     data <- data[sc1[1]:sc1[2]]
     data <- data[data !=";"]
     data1 <- data[!grepl('_', data)]
     lst <- lapply(split(data1, cumsum(grepl("^V", data1))),
                  function(x) read.table(text=x))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
     colnames(res) <- c("index_h","index_k","index_l","F_squared_meas","F_squared_sigma")
     res <- res
     } else {
	 res <- NA
	 return(res)
	 }
	}

r_reflections2 <- function (x){
  data <- unlist(x)
  nskip <- length((grep("_refln",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_refln",data,value=TRUE))
  l_l <- c(gsub(" ","",l_l))
  colnames(res) <- c(gsub("_refln_","",l_l))
  return(res)
}

r_spositions <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_atom",data,value=TRUE))
  data1 <- data[length(l_l)+1:length(data)]
  data1 <- data1[!is.na(data1)]
  cu <- grep("^\\_",data1)
  if ((length(cu)) >0){
      data2 <- data1[1:cu[1]-1]
    } else {
      data2 <- data1
    }
  #nskip <- length((grep("_atom",data)))
  lst <- lapply(split(data2, cumsum(grepl("^V", data2))),
                function(x) read.table(text=x))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  colnames(res) <- c(gsub("_atom_site_","",l_l))
  return(res)
}


r_saniso <- function(x){
  data <- unlist(x)
  l_l <- c(grep("_atom",data,value=TRUE))
  data1 <- data[length(l_l)+1:length(data)]
  data1 <- data1[!is.na(data1)]
  cu <- grep("^\\_",data1)
  if ((length(cu)) >0){
      data2 <- data1[1:cu[1]-1]
    } else {
      data2 <- data1
    }
  #nskip <- length((grep("_atom",data)))
  lst <- lapply(split(data2, cumsum(grepl("^V", data2))),
                function(x) read.table(text=x))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  colnames(res) <- c(gsub("_atom_site_aniso_","",l_l))
  return(res)
}

r_angle <- function(x){
  data <- unlist(x)
  l_l <- c(grep("_geom_angle",data,value=TRUE))
  data1 <- data[length(l_l)+1:length(data)]
  data1 <- data1[!is.na(data1)]
  cu <- grep("^\\_",data1)
  if ((length(cu)) >0){
      data2 <- data1[1:cu[1]-1]
    } else {
      data2 <- data1
    }
  #nskip <- length((grep("_geom_angle",data)))
  lst <- lapply(split(data2, cumsum(grepl("^V", data2))),
                function(x) read.table(text=x))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  colnames(res) <- c(gsub("_geom_angle_","",l_l))
  return(res)
}

r_dist <- function(x){
  data <- unlist(x)
  l_l <- c(grep("_geom_bond",data,value=TRUE))
  data1 <- data[length(l_l)+1:length(data)]
  data1 <- data1[!is.na(data1)]
  cu <- grep("^\\_",data1)
  if ((length(cu)) >0){
      data2 <- data1[1:cu[1]-1]
    } else {
      data2 <- data1
    }
  #nskip <- length((grep("_geom_bond",data)))
  lst <- lapply(split(data2, cumsum(grepl("^V", data2))),
                function(x) read.table(text=x))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  colnames(res) <- c(gsub("_geom_bond_","",l_l))
  return(res)
}

r_hbond <- function(x){
  data <- unlist(x)
  l_l <- c(grep("_geom_hbond",data,value=TRUE))
  data1 <- data[length(l_l)+1:length(data)]
  data1 <- data1[!is.na(data1)]
  cu <- grep("^\\_",data1)
  if ((length(cu)) >0){
      data2 <- data1[1:cu[1]-1]
    } else {
      data2 <- data1
    }
  #nskip <- length((grep("_geom_hbond",data)))
  lst <- lapply(split(data2, cumsum(grepl("^V", data2))),
                function(x) read.table(text=x))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  colnames(res) <- c(gsub("_geom_hbond_","",l_l))
  return(res)
}

r_torsion <- function(x){
  data <- unlist(x)
  l_l <- c(grep("_geom_torsion",data,value=TRUE))
  data1 <- data[length(l_l)+1:length(data)]
  data1 <- data1[!is.na(data1)]
  cu <- grep("^\\_",data1)
  if ((length(cu)) >0){
      data2 <- data1[1:cu[1]-1]
    } else {
      data2 <- data1
    }
  #nskip <- length((grep("_geom_torsion",data)))
  lst <- lapply(split(data2, cumsum(grepl("^V", data2))),
                function(x) read.table(text=x))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  colnames(res) <- c(gsub("_geom_torsion_","",l_l))
  return(res)
}



r_summ <- function(x){
  data <- unlist(x)
  formula <- check1("_chemical_formula_sum",data)
  a <- reap("_cell_length_a",data)
  b <- reap("_cell_length_b",data)
  c <- reap("_cell_length_c",data)
  al <- reap("_cell_angle_alpha",data)
  be <- reap("_cell_angle_beta",data)
  ga <- reap("_cell_angle_gamma",data)
  v <- reap("_cell_volume",data)
  den <- reap("_exptl_crystal_density_diffrn",data)
  c_sy <- check("_symmetry_cell_setting",data)
  sg_n <- as.numeric(check("_space_group_IT_number",data))
  sg_hall<- check1("_symmetry_space_group_name_Hall",data)
  sg_HM <- check1("_symmetry_space_group_name_H-M",data)
  cell <- list(A=a,B=b,C=c,ALPHA=al,BETA=be,GAMMA=ga)
  #sg <- list(c_sy,sg_n,sg_hall,sg_HM)
  #prop <- list(v,den)
  summ_c <- list(FORMULA=formula,CELL=cell,VOL=v,DEN=den,CrysSys=c_sy,SGN=sg_n,HALL=sg_hall,HM=sg_HM)
  return(summ_c)
}

r_msg <- function(x){
  data <- unlist(x)
  a <- reap("_cell_length_a",data)
  b <- reap("_cell_length_b",data)
  c <- reap("_cell_length_c",data)
  al <- reap("_cell_angle_alpha",data)
  be <- reap("_cell_angle_beta",data)
  ga <- reap("_cell_angle_gamma",data)
  cell <- c(a,b,c,al,be,ga)
  return(cell)
}


