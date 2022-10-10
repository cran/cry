#
# This file is part of the cry package
#
# Functions connected to reflections data.

#' Reads and output an CIF file
#'
#' @param filename A character string. The path to a valid small molecule
#'     reflection file. Typically a hkl or an fcf file.
#' @param message A logical variable. If TRUE (default) the function prints
#'    a message highlighting what is included in the cif file.
#' @return A named list. Each name correspond to a valid field in the fcf or hkl file.
#'
#' @examples
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"1dei-sf.cif")
#' lCIF <- readsm_REFL(filename)
#' print(names(lCIF))
#' print(lCIF$INTRO$CELL)
#' print(lCIF$SYMM)
#' print(lCIF$REFL)
#' @export
readsm_REFL <- function(filename, message=FALSE){
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  l_list <- grep("loop_",lcif)
  l_list1 <- l_find(l_list, length(lcif))
  #l_list1 <- append(l_list,(length(lcif)-1))
  mat<-zoo::rollapply(l_list1, 2,stack)
  ch <- apply(mat, 1, function(x) lcif[(x[1]+1):(x[2]-1)])
  #crystal_summary <- lapply(ch, ucheck, pattern="_audit.revision_id ")
  symmetry <- lapply(ch, ucheck, pattern="_xyz")
  reflection <- lapply(ch, ucheck, pattern="_refln_index_h")


  intro <- r_summ(lcif)
  symm <- if (is.na(nanona(symmetry)) == FALSE) r_symm(nanona(symmetry)) else NULL
  reflections <- if (is.na(nanona(reflection)) == FALSE) clean(r_reflections(nanona(reflection))) else NULL
  CIF = list(HEADER=intro,SYMM=symm,REFL=reflections)
  close(f)
  nrefs <- length(reflections$VAL$F_squared_meas)
  fmeas <- as.numeric(reflections$VAL$F_squared_meas)
  msg <- c("\n")
  if (message) {
    msg <- c(msg,sprintf("File %s read successfully.\n",filename))
    msg2 <- sprintf("There are %d reflections in this file.\n",nrefs)
    msg <- c(msg,msg2)
    msg <- c(msg,"Here is a summary of the observations:\n")
    msg <- c(msg,"\n")
    cat(msg)
    print(summary(fmeas))
  }
  return(CIF)
}



### accessory functions ####

l_find <- function(a,n){
  if (length(a) > 1){
    a1 <- append(a,(n-1))
  } else
  { a1 <- c(2,a,(n-1))
  return(a1)
  }
}

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
  if (isTRUE(count > 2)) { #something wrong in lower version of R
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

r_symm <- function (x){
  data <- unlist(x)
  data <- data[!grepl("_", data)]
  data <- data[data !=" "]
  res <- data[data !=""]
  return(res)
}

r_reflections <- function (x){
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

r_summ <- function(x){
  data <- unlist(x)
  id <- check1("_shelx_title",data)
  a <- reap("_cell_length_a",data)
  b <- reap("_cell_length_b",data)
  c <- reap("_cell_length_c",data)
  al <- reap("_cell_angle_alpha",data)
  be <- reap("_cell_angle_beta",data)
  ga <- reap("_cell_angle_gamma",data)
  sg_n <- as.numeric(check("_symmetry_Int_Tables_number",data))
  sg_hall<- check1("_symmetry_space_group_name_Hall",data)
  sg_HM <- check1("_symmetry_space_group_name_H-M",data)
  co <- check("_shelx_refln_list_code",data)
  F000 <- check("_exptl_crystal_F_000",data)
  hr <- check("_reflns_d_resolution_high",data)
  cell <- list(A=a,B=b,C=c,ALPHA=al,BETA=be,GAMMA=ga)
  summ_c <- list(TITLE=id,CELL=cell,SGN=sg_n,HALL=sg_hall,HM=sg_HM,SHELX_CODE=co,F_000=F000,HIGH_RES=hr)
  return(summ_c)
}
