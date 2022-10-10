#
# This file is part of the cry package
#
# Functions connected to reflections data.

#' Reads and output an mmCIF file
#'
#' @param filename A character string. The path to a valid
#'                 CIF file.
#' @param message A logical variable. If TRUE (default) the
#'                function prints a message highlighting
#'                what is included in the cif file.
#' @return A named list. Each name correspond to a valid
#'         field in the cif.
#'
#' @examples
#' datadir <- system.file("extdata",package="cry")
#' filename <- file.path(datadir,"3syu.cif")
#' lCIF <- readmm_CIF(filename)
#' print(names(lCIF))
#' print(lCIF$HEADER$Entry)
#' print(lCIF$HEADER$Symmtery)
#' print(lCIF$HEADER$CELL)
#' print(lCIF$EXP_DETAILS$CRYSTAL_CON$VAL)
#'
#' @export
readmm_CIF <- function(filename, message=FALSE){
  f <- file(filename)
  lcif <- readLines(f,warn=FALSE)
  l_list <- grep("loop_",lcif)
  l_list1 <- append(l_list,length(lcif))
  h_list <- grep("#",lcif)
  all <- append(l_list1,h_list)
  all <- sort(all)
  #all1 <- all[sapply(1:(length(all)-1), function(i) all[i+1] - all[i]) != 1]
  mat<-zoo::rollapply(all, 2,stack)
  ch <- apply(mat, 1, function(x) lcif[(x[1]+1):(x[2]-1)])
  cellp <- grep("\\b_cell.\\b", lcif, value=TRUE,perl=T)
  cellparam <- r_cellparm(cellp)
  id <- check("_entry.id",lcif)
  c_sy <- check("_symmetry_cell_setting",lcif)
  sg_n <- as.numeric(check("_space_group_IT_number",lcif))
  sg_hall<- check1("_symmetry_space_group_name_Hall",lcif)
  sg_HM <- check1("symmetry.space_group_name_H-M",lcif)

  symm <- grep("\\b_symmetry.", lcif, value=TRUE,perl=T)
  symmetry <- clean(r_symm(symm))

  expetls <- grep("\\b_exptl", lcif, value=TRUE,perl=T)
  exptl <- clean(r_exptl(expetls))

  cry_conds <- lapply(ch, ucheck, pattern="_exptl_crystal_grow.pdbx_details")
  cry_cond <- if (is.na(nanona(cry_conds)) == FALSE) clean(r_cry_con(nanona(cry_conds))) else NULL

  diffractn <- grep("\\b_diffrn", lcif, value=TRUE,perl=T)
  diffr <- clean(r_diff(diffractn))

  refl_block <- lapply(ch, ucheck, pattern="_refln.index_h")
  refl_data <- if (is.na(nanona(refl_block)) == FALSE) clean(r_reflections(nanona(refl_block))) else NULL

  reflections1 <- lapply(ch, ucheck, pattern="\\b_reflns.\\b")
  refl_all <- if (is.na(nanona(reflections1)) == FALSE) clean(r_refl1(nanona(reflections1))) else NULL

  reflections2 <- lapply(ch, ucheck, pattern="\\b_reflns_shell.\\b")
  refl_shell <- if (is.na(nanona(reflections2)) == FALSE) clean(r_refl2(nanona(reflections2))) else NULL

  refinement <- lapply(ch, ucheck, pattern="\\b_refine.\\b")
  refine_all <- if (is.na(nanona(refinement)) == FALSE) clean(r_refine(nanona(refinement))) else NULL

  refinement2 <- lapply(ch, ucheck, pattern="\\b_refine_hist.\\b")
  refine_hist <- if (is.na(nanona(refinement2)) == FALSE) clean(r_refineh(nanona(refinement2))) else NULL

  refinement3 <- lapply(ch, ucheck, pattern="\\b_refine_ls_shell.\\b")
  refine_shell <- if (is.na(nanona(refinement3)) == FALSE) clean(r_refinesh(nanona(refinement3))) else NULL

  refinement4 <- lapply(ch, ucheck, pattern="_pdbx_refine_tls.id")
  refine_tls <- if (is.na(nanona(refinement4)) == FALSE) clean(r_tls1(nanona(refinement4))) else NULL

  refinement5 <- lapply(ch, ucheck, pattern="_pdbx_refine_tls_group.id")
  refine_tls_g <- if (is.na(nanona(refinement5)) == FALSE) clean(r_tls_g(nanona(refinement5))) else NULL

  entities <- lapply(ch, ucheck, pattern="_entity.type")
  entity <- if (is.na(nanona(entities)) == FALSE) clean(r_entity(nanona(entities))) else NULL

  entities_poly <- lapply(ch, ucheck, pattern="_entity_poly.entity_id")
  entity_poly <- if (is.na(nanona(entities_poly)) == FALSE) clean(r_entity_p(nanona(entities_poly))) else NULL

  entities_src <- lapply(ch, ucheck, pattern="_entity_src_")
  entity_src <- if (is.na(nanona(entities_src)) == FALSE) clean(r_entity_s(nanona(entities_src))) else NULL

  sequences <- lapply(ch, ucheck, pattern="_entity_poly_seq.entity_id")
  seq_ent <- if (is.na(nanona(sequences)) == FALSE) clean(r_seqent(nanona(sequences))) else NULL

  strreference <- lapply(ch, ucheck, pattern="_struct_ref.id")
  #u <- length(unique(as.numeric(seq_ent$VAL[,1])))
  str_ref <- if (is.na(nanona(strreference)) == FALSE) clean(r_strref(nanona(strreference))) else NULL

  referenceseq <- lapply(ch, ucheck, pattern="_struct_ref_seq.align_id")
  ref_seq <- if (is.na(nanona(referenceseq)) == FALSE) clean(r_refseq(nanona(referenceseq))) else NULL

  mutantseq <- lapply(ch, ucheck, pattern="_struct_ref_seq_dif.align_id")
  mut_seq <- if (is.na(nanona(mutantseq)) == FALSE) clean(r_mutseq(nanona(mutantseq))) else NULL

  complist <- lapply(ch, ucheck, pattern="_chem_comp.id")
  chem_comp <- if (is.na(nanona(complist)) == FALSE) clean(r_compl(nanona(complist))) else NULL

  asyms_info <- lapply(ch, ucheck, pattern="_struct_asym.id")
  asym_info <- if (is.na(nanona(asyms_info)) == FALSE) clean(r_asym_info(nanona(asyms_info))) else NULL

  conformations <- lapply(ch, ucheck, pattern="_struct_conf.id")
  conf <- if (is.na(nanona(conformations)) == FALSE) clean(r_conf(nanona(conformations))) else NULL

  connets <- lapply(ch, ucheck, pattern="_struct_conn.id")
  cont <- if (is.na(nanona(connets)) == FALSE) clean(r_cont(nanona(connets))) else NULL

  monociss <- lapply(ch, ucheck, pattern="_struct_mon_prot_cis.pdbx_id")
  monocis <- if (is.na(nanona(monociss)) == FALSE) clean(r_monocis(nanona(monociss))) else NULL

  sheets1 <- lapply(ch, ucheck, pattern="_struct_sheet.id")
  sheet1 <- if (is.na(nanona(sheets1)) == FALSE) clean(r_sheet1(nanona(sheets1))) else NULL

  sheets2 <- lapply(ch, ucheck, pattern="_struct_sheet_order.sheet_id")
  sheet2 <- if (is.na(nanona(sheets2)) == FALSE) clean(r_sheet2(nanona(sheets2))) else NULL

  sheets3 <- lapply(ch, ucheck, pattern="_struct_sheet_range.sheet_id")
  sheet3 <- if (is.na(nanona(sheets3)) == FALSE) clean(r_sheet3(nanona(sheets3))) else NULL

  sheets4 <- lapply(ch, ucheck, pattern="_pdbx_struct_sheet_hbond.sheet_id")
  sheet4 <- if (is.na(nanona(sheets4)) == FALSE) clean(r_sheet4(nanona(sheets4))) else NULL

  nscdlim <- lapply(ch, ucheck, pattern="_struct_ncs_dom_lim.dom_id")
  nscdlimt <- if (is.na(nanona(nscdlim)) == FALSE) clean(r_nsclim(nanona(nscdlim))) else NULL

  nscdo <- lapply(ch, ucheck, pattern="_struct_ncs_dom.id")
  nscd <- if (is.na(nanona(nscdo)) == FALSE) clean(r_nscd(nanona(nscdo))) else NULL

  na1 <- lapply(ch, ucheck, pattern="_ndb_struct_conf_na.entry_id")
  na_conf <- if (is.na(nanona(na1)) == FALSE) clean(r_na1(nanona(na1))) else NULL

  na2 <- lapply(ch, ucheck, pattern="_ndb_struct_na_base_pair.i_label_asym_id")
  na_b_int <- if (is.na(nanona(na2)) == FALSE) clean(r_na2(nanona(na2))) else NULL

  na3 <- lapply(ch, ucheck, pattern="_ndb_struct_na_base_pair_step.i_label_asym_id_1")
  na_s_int <- if (is.na(nanona(na3)) == FALSE) clean(r_na3(nanona(na3))) else NULL

  transforms <- lapply(ch, ucheck, pattern="_atom_sites.fract_transf_matrix")
  tranform <- if (is.na(nanona(transforms)) == FALSE) clean(r_trans(nanona(transforms))) else NULL

  coordinates <- lapply(ch, ucheck, pattern="_atom_site.Cartn_x")
  coordinate <- if (is.na(nanona(coordinates)) == FALSE) clean(r_positions(nanona(coordinates))) else NULL

  anisots <- lapply(ch, ucheck, pattern="_atom_site_anisotrop.id")
  anisot <- if (is.na(nanona(anisots)) == FALSE) clean(r_aniso(nanona(anisots))) else NULL

  asym_npo <- lapply(ch, ucheck, pattern="_pdbx_nonpoly_scheme.asym_id")
  asymnonpoly <- if (is.na(nanona(asym_npo)) == FALSE) clean(r_asymnp(nanona(asym_npo))) else NULL

  asym_po <- lapply(ch, ucheck, pattern="_pdbx_poly_seq_scheme.asym_id")
  asympoly <- if (is.na(nanona(asym_po)) == FALSE) clean(r_asymp(nanona(asym_po))) else NULL

  softcits <- lapply(ch, ucheck, pattern="_software.citation_id")
  softcite <- if (is.na(nanona(softcits)) == FALSE) clean(r_softc(nanona(softcits))) else NULL

  close_conct <- lapply(ch, ucheck, pattern="_pdbx_validate_close_contact.id")
  close_c <- if (is.na(nanona(close_conct)) == FALSE) clean(r_closec(nanona(close_conct))) else NULL

  valid_angle <- lapply(ch, ucheck, pattern="_pdbx_validate_rmsd_angle.id")
  val_angle <- if (is.na(nanona(valid_angle)) == FALSE) clean(r_valang(nanona(valid_angle))) else NULL

  valid_tor <- lapply(ch, ucheck, pattern="_pdbx_validate_torsion.id")
  val_tor <- if (is.na(nanona(valid_tor)) == FALSE) clean(r_valtor(nanona(valid_tor))) else NULL

  valid_omg <- lapply(ch, ucheck, pattern="_pdbx_validate_peptide_omega.id")
  val_omg <- if (is.na(nanona(valid_omg)) == FALSE) clean(r_valomg(nanona(valid_omg))) else NULL

  zero_atom <- lapply(ch, ucheck, pattern="_pdbx_unobs_or_zero_occ_atoms.id")
  zo_atom <- if (is.na(nanona(zero_atom)) == FALSE) clean(r_zoatom(nanona(zero_atom))) else NULL

  zero_residue <- lapply(ch, ucheck, pattern="_pdbx_unobs_or_zero_occ_residues.id")
  zo_res <- if (is.na(nanona(zero_residue)) == FALSE) clean(r_zores(nanona(zero_residue))) else NULL

  intro = list(Entry=id,Symmtery=symmetry,CELL=cellparam)
  refn = list(Overall=refine_all,HIST=refine_hist,SHELL=refine_shell,TLS=refine_tls,TLS_Group=refine_tls_g)
  refl_d = list(Overall=refl_all,SHELL=refl_shell)
  ent = list(ENTITY_all=entity,ENTITY_Poly=entity_poly,ENTITY_Source=entity_src)
  bp_d = list(NA_CONF=na_conf,NA_BP_INT=na_b_int,NA_S_INT=na_s_int)
  expr = list(EXPERIMENT=exptl,CRYSTAL_CON=cry_cond,DIFFRACTION=diffr,REFLECTION=refl_d,REFINEMENT=refn)
  str_seq = list(ENTITY=ent,SEQ=seq_ent,STR_Ref=str_ref,SEQ_SOURCE=ref_seq,SEQ_ALT=mut_seq,COMPOSITION=chem_comp,ASYM=asym_info)#
  sheetin = list(SHEETID=sheet1,SHEET_ORDER=sheet2,SHEET_RANGE=sheet3,SHEET_HBOND=sheet4)
  conn_nsc = list(CONFIRMATIONS=conf,CONNECTIONS=cont,PEP_CIS=monocis,SHEET_INFO=sheetin,NSC_LIMIT=nscdlimt,NSC=nscd)
  Str_d = list(SEQ_STR=str_seq,CONN_NSC=conn_nsc,TRANS_INFO=tranform,COOR=coordinate,NA_INFO=bp_d,ANISO=anisot,ASYMNP=asymnonpoly,ASYMP=asympoly)
  val = list(CLOSE_CONT=close_c,VAL_ANGLE=val_angle,VAL_TOR=val_tor,VAL_OMG=val_omg,ZO_atom=zo_atom,ZO_residue=zo_res)
  CIF = list(HEADER=intro,REFL=refl_data,EXP_DETAILS=expr,STRU_DETAILS=Str_d,SOFTWARE=softcite,VAL_DETAILS=val)
  close(f)
 if (message) {
     if (!is.null(refl_data)){
     n <- length(refl_data$VAL$F_meas_au)
     f <- as.numeric(refl_data$VAL$F_meas_au)
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
	 please refer corresponding reflection file (sfcif or mtz).\n")
     msg2 <- sprintf("There are %d atoms in the molecule.\n",anum)
	 out <- c(msg,msg1,msg2)
     cat(out)
	 }
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

r_entity <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_entity.",data,value=TRUE))
  m <- length(l_l)
  n <- length(data)
  data_ex <- data[n]
  data_ex <- (scan(text=data_ex, what='character', quiet=TRUE)[1])
  data1 <- data[m+1:n]
  data1 <- nanona(data1)
  data2 <- scan(text=data1, what='character', quiet=TRUE)
  data2 <- gsub("[;]","",data2)
  data2 <- data2[data2 !=""]
  o <- length(data2)
  data3 <- list(data2)
  list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
  res <- do.call(rbind, list_all)
  res <- as.data.frame(res)
  colnames(res) <- c(gsub("_entity","",l_l))
  return(res)
 }

r_entity_p <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_entity_poly",data,value=TRUE))
  m <- length(l_l)
  c <- check("_entity_poly",data[1])
  c[is.na(c)] <- 0
  c <- as.numeric(c)
  if (c == 1) {
        data1 <- gsub(";", "'", data, fixed = T)
	    data2 <- scan(text=data1, what='character', quiet=TRUE)
	    data2 <- gsub("[\r\n]", "", data2)
	    o <- length(data2)
	    data3 <- list(data2)
	    list_all <- split(data3[[1]], rep(1:m, each = (length(data2)/length(l_l))))
	    res <- do.call(rbind, list_all)
	    res <- as.data.frame(res)
	    colnames(res) <- c("KEY","VAL")
		res <- res
      } else if (c == 0) {
        n <- length(data)
        data1 <- data[m+1:n]
        data1 <- nanona(data1)
        d <- length(grep(";",data1))
        if (d < 1) {
            data1 <- unlist(data1)
	        } else {
	        data1 <- unlist(data1)
	        data1 <- gsub(";", '"', data1, fixed = T)
	        }
        data2 <- scan(text=data1, what='character', quiet=TRUE)
        data2 <- gsub("[\r\n]", "", data2)
        o <- length(data2)
        data3 <- list(data2)
        list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
        res <- do.call(rbind, list_all)
        res <- as.data.frame(res)
        colnames(res) <- c(gsub("_entity_poly","",l_l))
        return(res)
        }
	}


r_entity_s <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_entity_src_",data))))-1)
  if (l < 1) {
     data <- (gsub("_entity_src_","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_entity_src_",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    d <- length(grep(";",data1))
    if (d < 1) {
       data1 <- unlist(data1)
	   } else {
	     data1 <- unlist(data1)
	     data1 <- gsub(";", '"', data1, fixed = T)
	   }
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    data2 <- gsub("[\r\n]", "", data2)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_entity_src_","",l_l))
    return(res)
  }
 }

r_strref <- function (x,u){
  data <- unlist(x)
  c <- check("_struct_ref.",data[1])
  c[is.na(c)] <- 0
  c <- as.numeric(c)
  if (c == 1) {
	  data1 <- gsub(";", "'", data, fixed = T)
	  data2 <- scan(text=data1, what='character', quiet=TRUE)
	  data2 <- gsub("[\r\n]", "", data2)
	  l_l <- c(grep("_struct_ref.",data,value=TRUE))
	  m <- length(l_l)
	  o <- length(data2)
	  data3 <- list(data2)
	  list_all <- split(data3[[1]], rep(1:m, each = (length(data2)/length(l_l))))
	  res <- do.call(rbind, list_all)
	  res <- as.data.frame(res)
	  colnames(res) <- c("KEY","VAL")
	  res <- res
      } else if (c == 0) {
	    l_l <- c(grep("_struct_ref.",data,value=TRUE))
        m <- length(l_l)
        n <- length(data)
        data1 <- data[m+1:n]
        data1 <- nanona(data1)
        data1 <- gsub(";", '"',data1,fixed = T)
		data2 <- scan(text=data1, what='character', quiet=TRUE)
        data2 <- data2[data2 !=""]
        o <- length(data2)
        data3 <- list(data2)
        list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
        res <- do.call(rbind, list_all)
        res <- as.data.frame(res)
        colnames(res) <- c(gsub("_struct_ref.","",l_l))
		return(res)
	  }
	}

r_seqent <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_entity",data))))-1)
  if (l < 1) {
     data <- (gsub("_entity_poly_seq.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_entity_poly_seq.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_entity_poly_seq.",data,value=TRUE))
    colnames(res) <- c(gsub("_entity_poly_seq.","",l_l))
	return(res)
  }
 }

r_refseq <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_ref_seq.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_struct_ref_seq.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_struct_ref_seq.",data,value=TRUE))
    colnames(res) <- c(gsub("_struct_ref_seq.","",l_l))
	return(res)
  }
 }


r_mutseq <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_ref_seq_dif.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_struct_ref_seq_dif.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_struct_ref_seq_dif.",data,value=TRUE))
    colnames(res) <- c(gsub("_struct_ref_seq_dif.","",l_l))
	return(res)
  }
 }

r_compl <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_chem_comp",data,value=TRUE))
  m <- length(l_l)
  n <- length(data)
  data1 <- data[m+1:n]
  data1 <- nanona(data1)
  data2 <- scan(text=data1, what='character',quiet=TRUE)
  lsc <- length(grep("[;]",data2))
  if (lsc <= 1) {
     data2 <- data2
	 } else {
	   if (length(grep("\"", data1, fixed = TRUE)) > 1){
	       data2 <- gsub("[;]","",data2)
           data2 <- data2[data2 !=""]
	    } else {
	      data2 <- scan(text=data2, what='character', sep=";",quiet=TRUE)
          data2 <- data2[data2 !=""]
	    }
		}
  o <- length(data2)
  data3 <- list(data2)
  list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
  res <- do.call(rbind, list_all)
  res <- as.data.frame(res)
  colnames(res) <- c(gsub("_chem_comp","",l_l))
  return(res)
 }

 r_asym_info <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct_asym.",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_asym.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_struct_asym.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_struct_asym.","",l_l))
	return(res)
  }
 }


r_conf <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_conf.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_struct_conf.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_struct_conf.",data,value=TRUE))
    colnames(res) <- c(gsub("_struct_conf.","",l_l))
	return(res)
  }
 }


r_cont <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_struct",data,value=TRUE))
  m <- length(l_l)
  n <- length(data)
  data_ex <- data[n]
  data_ex <- (scan(text=data_ex, what='character', quiet=TRUE)[1])
  data1 <- data[m+1:n]
  data1 <- nanona(data1)
  data2 <- scan(text=data1, what='character', quiet=TRUE)
  data2 <- gsub("[;]","",data2)
  data2 <- data2[data2 !=""]
  o <- length(data2)
  data3 <- list(data2)
  list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
  res <- do.call(rbind, list_all)
  res <- as.data.frame(res)
  colnames(res) <- c(gsub("_struct_conn.","",l_l))
  return(res)
 }

r_monocis <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_mon_prot_cis.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_struct_mon_prot_cis.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_struct_mon_prot_cis.",data,value=TRUE))
    colnames(res) <- c(gsub("_struct_mon_prot_cis.","",l_l))
	return(res)
  }
 }

 r_sheet1 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_sheet.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_struct_sheet.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_struct_sheet.","",l_l))
	return(res)
  }
 }

r_sheet2 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_sheet_order.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_struct_sheet_order.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_struct_sheet_order.","",l_l))
	return(res)
  }
 }

r_sheet3 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_sheet_range.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_struct_sheet_range.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_struct_sheet_range.","",l_l))
	return(res)
  }
 }

r_sheet4 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
     data <- (gsub("_pdbx_struct_sheet_hbond.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_pdbx_struct_sheet_hbond.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_pdbx_struct_sheet_hbond.","",l_l))
	return(res)
  }
 }

r_nsclim <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_ncs_dom_lim.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_struct_ncs_dom_lim.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_struct_ncs_dom_lim.",data,value=TRUE))
    colnames(res) <- c(gsub("_struct_ncs_dom_lim.","",l_l))
	return(res)
  }
 }

r_nscd <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_struct_ncs_dom.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_struct_ncs_dom.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_struct_ncs_dom.",data,value=TRUE))
    colnames(res) <- c(gsub("_struct_ncs_dom.","",l_l))
	return(res)
  }
 }


r_na1 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_ndb_struct",data))))-1)
  if (l < 1) {
     data <- (gsub("_ndb_struct_conf_na.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_ndb_struct_conf_na.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_ndb_struct_conf_na.",data,value=TRUE))
    colnames(res) <- c(gsub("_ndb_struct_conf_na.","",l_l))
	return(res)
  }
 }

 r_na2 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_ndb_struct_na_base_pair.",data))))-1)
  if (l < 1) {
     data <- (gsub("_ndb_struct_na_base_pair.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_ndb_struct_na_base_pair.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_ndb_struct_na_base_pair.","",l_l))
	return(res)
  }
 }


 r_na3 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_ndb_struct_na_base_pair_step.",data))))-1)
  if (l < 1) {
     data <- (gsub("_ndb_struct_na_base_pair_step.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_ndb_struct_na_base_pair_step.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_ndb_struct_na_base_pair_step.","",l_l))
	return(res)
  }
 }

r_trans <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_atom_",data))))-1)
  if (l < 1) {
     data <- (gsub("_atom_sites.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_atom_sites.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_atom_sites.",data,value=TRUE))
    colnames(res) <- c(gsub("_atom_sites.","",l_l))
	return(res)
  }
 }

r_positions <- function (x){
  data <- unlist(x)
  nskip <- length((grep("_atom_site",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_atom_site.",data,value=TRUE))
  colnames(res) <- c(gsub("_atom_site.","",l_l))
  return(res)
}

r_aniso <- function (x){
  data <- unlist(x)
  nskip <- length((grep("_atom_site_",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_atom_site_",data,value=TRUE))
  colnames(res) <- c(gsub("_atom_site_","",l_l))
  return(res)
}

r_asymnp <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_pdbx_nonpoly_scheme.",data,value=TRUE))
  m <- length(l_l)
  n <- length(data)
  data1 <- data[m+1:n]
  data1 <- nanona(data1)
  data2 <- scan(text=data1, what='character', quiet=TRUE)
  o <- length(data2)
  data3 <- list(data2)
  list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
  res <- do.call(rbind, list_all)
  res <- as.data.frame(res)
  colnames(res) <- c(gsub("_pdbx_nonpoly_scheme.","",l_l))
  return(res)
}


r_asymp <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
     data <- (gsub("_pdbx_poly_seq_scheme.","",data))
     lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
     names(lst) <- NULL
     res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_pdbx_poly_seq_scheme.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_pdbx_poly_seq_scheme.","",l_l))
	return(res)
  }
 }


r_softc <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_software",data))))-1)
  if (l < 1) {
      data <- (gsub("_software.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    l_l <- c(grep("_software.",data,value=TRUE))
    m <- length(l_l)
    n <- length(data)
    data1 <- data[m+1:n]
    data1 <- nanona(data1)
    data2 <- scan(text=data1, what='character', quiet=TRUE)
    o <- length(data2)
    data3 <- list(data2)
    list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
    res <- do.call(rbind, list_all)
    res <- as.data.frame(res)
    colnames(res) <- c(gsub("_software.","",l_l))
	return(res)
  }
 }


r_closec <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
      data <- (gsub("_pdbx_validate_close_contact.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_pdbx_validate_close_contact.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                 function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_pdbx_validate_close_contact.",data,value=TRUE))
    colnames(res) <- c(gsub("_pdbx_validate_close_contact.","",l_l))
	return(res)
  }
 }


r_valang <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
      data <- (gsub("_pdbx_validate_rmsd_angle.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_pdbx_validate_rmsd_angle.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                 function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_pdbx_validate_rmsd_angle.",data,value=TRUE))
    colnames(res) <- c(gsub("_pdbx_validate_rmsd_angle.","",l_l))
	return(res)
  }
 }

r_valtor <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
      data <- (gsub("_pdbx_validate_torsion.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_pdbx_validate_torsion.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_pdbx_validate_torsion.",data,value=TRUE))
    colnames(res) <- c(gsub("_pdbx_validate_torsion.","",l_l))
	return(res)
  }
 }

r_valomg <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
      data <- (gsub("_pdbx_validate_peptide_omega.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_pdbx_validate_peptide_omega.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_pdbx_validate_peptide_omega.",data,value=TRUE))
    colnames(res) <- c(gsub("_pdbx_validate_peptide_omega.","",l_l))
	return(res)
  }
 }

r_zoatom <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
      data <- (gsub("_pdbx_unobs_or_zero_occ_atoms.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_pdbx_unobs_or_zero_occ_atoms.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_pdbx_unobs_or_zero_occ_atoms.",data,value=TRUE))
    colnames(res) <- c(gsub("_pdbx_unobs_or_zero_occ_atoms.","",l_l))
	return(res)
  }
 }

r_zores <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_pdbx",data))))-1)
  if (l < 1) {
      data <- (gsub("_pdbx_unobs_or_zero_occ_residues.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_pdbx_unobs_or_zero_occ_residues.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_pdbx_unobs_or_zero_occ_residues.",data,value=TRUE))
    colnames(res) <- c(gsub("_pdbx_unobs_or_zero_occ_residues.","",l_l))
	return(res)
  }
 }

r_cellparm <- function (x){
  data <- unlist(x)
  #data <- data[!grepl("_cell.entry.id", data)]
  data <- data[!grepl("_cell.details", data)]
  data <- (gsub("_cell.","",data))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x ,stringsAsFactors = FALSE))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  res <- suppressWarnings(data.frame(res[1], apply(res[2], 1,numas)))
  colnames(res) <- c("KEY","VAL")
  return(res)
}


r_symm <- function (x){
  data <- unlist(x)
  data <- data[!grepl("_symmetry.entry.id", data)]
  data <- (gsub("_symmetry.","",data))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  res <- suppressWarnings(data.frame(res[1], apply(res[2], 2,numas)))
  return(res)
}

r_exptl <- function (x){
  data <- unlist(x)
  data <- data[!grepl("_details", data)]
  data <- (gsub("_exptl.","",data))
  data <- (gsub("_exptl_","",data))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  return(res)
}

r_cry_con <- function(x){
  data <- unlist(x)
  n <- length(grep("_exptl",data))
  m <- length(data)
  r <- grep("_exptl_crystal_grow.pdbx_details", data, value = TRUE)
  #r1 <- if(length(r) > 0) (strsplit(r, "\\s+")[[1]])[2] else NA
  r1 <- if(length(r) > 0) (strsplit(r, "'")[[1]])[2] else NA
  if (is.na(r1) == TRUE){
	  if (n+1 == m){
      r1 <- data[n+1]
      } else {
      r1 <- paste(data[n+1],data[m])
	  r1 <- strsplit(r1, ";")[[1]][2]
      }
	  }
   res <- strsplit(r1, ",")[[1]]
   res <- as.data.frame(list(res),col.names=c("VAL"))
   return(res)
   }

r_diff <- function (x){
  data <- unlist(x)
  data <- (gsub("_diffrn.","",data))
  data <- (gsub("_diffrn_","",data))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  return(res)
}

r_refl1 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_reflns.",data))))-1)
  if (l < 1) {
      data <- (gsub("_reflns.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_reflns.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_reflns.",data,value=TRUE))
    colnames(res) <- c(gsub("_reflns.","",l_l))
	return(res)
  }
 }

r_refl2 <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_reflns_shell.",data))))-1)
  if (l < 1) {
      data <- (gsub("_reflns_shell.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_reflns_shell.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_reflns_shell.",data,value=TRUE))
    colnames(res) <- c(gsub("_reflns_shell.","",l_l))
	return(res)
  }
 }

r_refine <- function (x){
  data <- unlist(x)
  d <- length(grep(";",data))
  if (d < 1) {
      data <- (gsub("_refine.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
	} else {
	  data1 <- gsub(";", "'", data, fixed = T)
	  data2 <- scan(text=data1, what='character', quiet=TRUE)
	  l_l <- c(grep("_refine.",data,value=TRUE))
	  m <- length(l_l)
	  o <- length(data2)
	  data3 <- list(data2)
	  list_all <- split(data3[[1]], rep(1:m, each = (length(data2)/length(l_l))))
	  res <- do.call(rbind, list_all)
	  res <- as.data.frame(res)
	  colnames(res) <- c("KEY","VAL")
      return(res)
	  }
  }


r_refineh <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_refine_hist.",data))))-1)
  if (l < 1) {
      data <- (gsub("_refine_hist.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_refine_hist.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_refine_hist.",data,value=TRUE))
    colnames(res) <- c(gsub("_refine_hist.","",l_l))
	return(res)
  }
 }

r_refinesh <- function (x){
  data <- unlist(x)
  l <- ((length(data)-length((grep("_refine_ls_shell.",data))))-1)
  if (l < 1) {
      data <- (gsub("_refine_ls_shell.","",data))
      lst <- lapply(split(data, cumsum(grepl("^V", data))),
                   function(x) read.table(text=x ,stringsAsFactors = FALSE,col.names=c("KEY","VAL")))
      names(lst) <- NULL
      res <- do.call(`cbind`, lst)
  } else {
    nskip <- length((grep("_refine_ls_shell.",data)))
    lst <- lapply(split(data, cumsum(grepl("^V", data))),
                  function(x) read.table(text=x,skip=nskip))
    names(lst) <- NULL
    res <- do.call(`cbind`, lst)
    l_l <- c(grep("_refine_ls_shell.",data,value=TRUE))
    colnames(res) <- c(gsub("_refine_ls_shell.","",l_l))
	return(res)
  }
 }


r_tls1 <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_pdbx_refine_tls.",data,value=TRUE))
  l <- ((length(data)-length((grep("_pdbx_refine_tls.",data)))))
  m <- length(l_l)
  n <- length(data)
  if (l < 1) {
      data <- (gsub("_pdbx_refine_tls.","",data))
      data2 <- scan(text=data, what='character', quiet=TRUE)
      data3 <- list(data2)
	  list_all <- split(data3[[1]], rep(1:m, each = (length(data2)/length(l_l))))
	  res <- do.call(rbind, list_all)
	  res <- as.data.frame(res)
	  colnames(res) <- c("KEY","VAL")
      return(res)
	} else {
	  data1 <- data[m+1:n]
      data1 <- nanona(data1)
      data2 <- scan(text=data1, what='character', quiet=TRUE)
      data2 <- gsub("[;]","",data2)
      data2 <- data2[data2 !=""]
      o <- length(data2)
      data3 <- list(data2)
      list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
      res <- do.call(rbind, list_all)
      res <- as.data.frame(res)
      colnames(res) <- c(gsub("_pdbx_refine_tls.","",l_l))
      return(res)
    }
  }


r_tls_g <- function (x){
  data <- unlist(x)
  l_l <- c(grep("_pdbx_refine_tls_group",data,value=TRUE))
  l <- ((length(data)-length((grep("_pdbx_refine_tls.",data)))))
  m <- length(l_l)
  n <- length(data)
  if (l < 1) {
      data <- (gsub("_pdbx_refine_tls_group.","",data))
      data2 <- scan(text=data, what='character', quiet=TRUE)
      data3 <- list(data2)
	  list_all <- split(data3[[1]], rep(1:m, each = (length(data2)/length(l_l))))
	  res <- do.call(rbind, list_all)
	  res <- as.data.frame(res)
	  colnames(res) <- c("KEY","VAL")
      return(res)
	} else {
      data1 <- data[m+1:n]
      data1 <- nanona(data1)
      d <- length(grep(";",data1))
      if (d < 1) {
         data1 <- unlist(data1)
	     } else {
	       data1 <- unlist(data1)
	       data1 <- gsub(";", '"', data1, fixed = T)
	      }
      data2 <- scan(text=data1, what='character', quiet=TRUE)
      data2 <- gsub("[\r\n]", "", data2)
      o <- length(data2)
      data3 <- list(data2)
      list_all <- split(data3[[1]], rep(1:(length(data2)/length(l_l)), each = m))
      res <- do.call(rbind, list_all)
      res <- as.data.frame(res)
      colnames(res) <- c(gsub("_pdbx_refine_tls_group","",l_l))
    return(res)
   }
  }

r_reflections <- function (x){
  data <- unlist(x)
  nskip <- length((grep("_refln",data)))
  lst <- lapply(split(data, cumsum(grepl("^V", data))),
                function(x) read.table(text=x,skip=nskip))
  names(lst) <- NULL
  res <- do.call(`cbind`, lst)
  l_l <- c(grep("_refln",data,value=TRUE))
  colnames(res) <- c(gsub("_refln.","",l_l))
  return(res)
}
