#' parse mutation
#'
#' @param x a string mutation: "-94G>A"
#'
#' @return a vector of different element of the mutation: c("-","94","G","A")
#'
#' @importFrom stringr str_extract
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
parse_mutation <- function(x){
  
  if(!is.character(x)){
    stop("x must be a string like '-94G>A'")
  }
  
  sig <- str_extract(x, "^-") ##################################################################################################
  if(is.na(sig)){
    sig <- str_extract(x, "^\\*")
  }
  pos <- str_extract(x, "[0-9]+")
  mut <- substring(gsub(pos, ">", x),3)
  mut <- str_split_fixed(mut, ">", 2)
  
  c(sig,pos,mut)
}
 
#' parse mutation in exonic region
#'
#' @param x a string mutation: "T22G"
#'
#' @return a vector of different element of the mutation: c("+","22","T","G")
#'
#' @importFrom stringr str_extract
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
parse_mutation_exonic <- function(x){
  if(!is.character(x)){
    stop("x must be a string like 'c.471C>A'")
  }
  
  sig <- ""
  
  pos <- str_extract(x, "[0-9]+")
  mut <- gsub(pos, ">", x)
  mut <- str_split_fixed(mut, ">", 2)
  pos <- gsub("-","",pos)
  
  c(sig,pos,mut)
  
}
#' parse refGene
#'
#' @param x a string refGene mutation: "CCDC116:NM_001331066:exon2:c.T22G:p.S8A"
#'
#' @return a vector of different element of the mutation: c("NM_001010939","c.T22G","+","22","T","G")
#'
#' @importFrom stringr str_split
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
parse_AAChange.refGene <- function(x){

  x <- str_split_fixed(x, ":", 5)

  y <- str_split_fixed(x[,4], "\\.", 5)[,2]

  z <- t(sapply(y, parse_mutation_exonic))

  w <- str_split_fixed(x[,2], "\\_", 2)

  a <- cbind(w[,1],x,z)
  dimnames(a) <- NULL

  a <-  subset(a, select = -c(2,3,4,6) )
  a
}

#' parse refGene
#'
#' @param x a string refGene mutation: "NM_001010939:c.-94G>A"
#'
#' @return a vector of different element of the mutation: c("NM_001010939","c.-94G>A","-","94","G","A")
#'
#' @importFrom stringr str_split
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
parse_GeneDetail.refGene <- function(x){

  if(!is.character(x)){
    stop("x must be a string like 'ENST_001010939:c.-94G>A'")
  }

  # Fix wrong character encoding!?!?!
  x <- gsub("\\\\x3b",";",x)

  x <- str_split(x, ";")[[1]]

  x <- str_split_fixed(x, ":", 2)

  w <- str_split_fixed(x[,1], "\\_", 2)

  y <- str_split_fixed(x[,2], "\\.", 2)[,2] ########################################################################################################
  z <- t(sapply(y, parse_mutation))

  a <- cbind(w[,1],x,z)
  dimnames(a) <- NULL

  a <-  subset(a, select = -c(2) )

  a
}

#' calculate kozak score
#'
#' @param sequence a sequence of 23 nucleotides : "GAATTATTTTATGCTATCATGAT"
#'
#' @return a the value of the kozak similarity score 
#'
#' @importFrom stringr str_split
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
kozak_similarity_score <- function(sequence){
  #0=A, 1=T, 2=G, 3=C, 4=N (Missing)
  weights <- list(
    c(0.04210526, 0.        , 0.03157895, 0.05263158, 0.        ),
    c(0.04210526, 0.05263158, 0.10526316, 0.0625    , 0.        ),
    c(0.03157895, 0.04210526, 0.05263158, 0.07368421, 0.        ),
    c(0.03157895, 0.01052632, 0.04210526, 0.05263158, 0.        ),
    c(0.08421053, 0.07368421, 0.18947368, 0.10526316, 0.        ),
    c(0.04210526, 0.05263158, 0.05263158, 0.08421053, 0.        ),
    c(0.12631579, 0.0625    , 0.12631579, 0.21052632, 0.        ),
    c(0.83157895, 0.12631579, 0.65263158, 0.16842105, 0.        ),
    c(0.15789474, 0.06315789, 0.11578947, 0.2       , 0.        ),
    c(0.21052632, 0.09473684, 0.31578947, 0.51578947, 0.        ),
    c(0.        , 0.        , 0.        , 0.        , 0.        ),
    c(0.        , 0.        , 0.        , 0.        , 0.        ),
    c(0.        , 0.        , 0.        , 0.        , 0.        ),
    c(0.24210526, 0.16666667, 0.53684211, 0.13684211, 0.        ),
    c(0.15789474, 0.09473684, 0.09473684, 0.24210526, 0.        ),
    c(0.05263158, 0.08421053, 0.14736842, 0.09473684, 0.        ),
    c(0.07216495, 0.05263158, 0.10526316, 0.06315789, 0.        ),
    c(0.        , 0.        , 0.        , 0.05263158, 0.        ),
    c(0.05263158, 0.05263158, 0.10526316, 0.09473684, 0.        ),
    c(0.04210526, 0.03157895, 0.05263158, 0.04210526, 0.        ),
    c(0.        , 0.        , 0.        , 0.        , 0.        ),
    c(0.04210526, 0.04210526, 0.08421053, 0.07368421, 0.        ),
    c(0.0625    , 0.04210526, 0.09473684, 0.05263158, 0.        )
  )
  
  
  #sequence = "GAATTATTTTATGCTATCATGAT"
  #sequence <- toupper(sequence)
  sequence <- str_split(sequence,"")
  
  numbers <- c(rep(0, nchar(sequence)))
  
  for(i in 1:length(sequence[[1]])){
    if(sequence[[1]][i] == 'A'){
      numbers[i] <- 1
    }else if(sequence[[1]][i] == 'T'){
      numbers[i] <- 2
    }else if(sequence[[1]][i] == 'G'){
      numbers[i] <- 3
    }else if(sequence[[1]][i] == 'C'){
      numbers[i] <- 4
    }
  }
  
  score <- 0
  for(i in 1:length(sequence[[1]])){
    score <- score + weights[[i]][numbers[i]]
  }
  
  
  max_index <- lapply(weights,which.max)
  max_score <- 0
  for(i in 1:length(weights)){
    max_score <- max_score + (weights[[i]])[max_index[[i]][1]]
  }
  
  #max_score = 3.7368421199999995
  score <- score/max_score
  return(score)  
}



#' get data used by MORFEE
#'
#' @param path is a path to the temporary folder to store downloaded databases
#' @param GRCh is the version of Human reference genome to download (i.e. 37 or 38)
#' @param GENCODE is the version of GENCODE database to download (i.e. 19 or NA)
#' @param force a boolean wheter the function should check if the database is already stored before to download them
#'
#' @importFrom utils download.file
#' @importFrom rtracklayer readGFF
#' @importFrom data.table fread
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_split_fixed
#'
#' @export
#'
get.morfee.data <- function(path=tempdir(), GRCh=37, GENCODE=NA, force=FALSE){

  if(!(GRCh %in% c(37,38))){
    stop('The argument "GRCh" must be equal to "37" or "38"')
  }
  if(!(is.na(GENCODE) | GENCODE==19)){
    stop('The argument "GENCODE" must be equal to "NA" or "19"')
  }

  MORFEE.DATA <- list()

  MORFEE.DATA[["DB_PATH"]] <- paste0(path,"/MORFEE_DB/") # Path to the DB
  MORFEE.DATA[["DB_PATH_ANNOVAR"]] <- paste0(MORFEE.DATA[["DB_PATH"]],"ANNOVAR","/") # Path to the ANNOVAR DB
  MORFEE.DATA[["DB_PATH_GENCODE"]] <- paste0(MORFEE.DATA[["DB_PATH"]],"GENCODE","/") # Path to the GENCODE DB

  MORFEE.DATA[["GRCh"]] <- GRCh

  if(GRCh==37){
    if(is.na(GENCODE)){
      MORFEE.DATA[["GENCODE"]] <- 43
      MORFEE.DATA[["GENCODE_URL"]] <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh37_mapping/"
      MORFEE.DATA[["GENCODE_FILE_ANNOT"]] <- "gencode.v43lift37.annotation.gff3.gz" # Position codon start + strand
      MORFEE.DATA[["GENCODE_FILE_METAD"]] <- "gencode.v43lift37.metadata.RefSeq.gz" # NM code
      MORFEE.DATA[["GENCODE_FILE_SEQUE"]] <- "GRCh37.primary_assembly.genome.fa.gz" # Whole Human Sequence
    }else{
      MORFEE.DATA[["GENCODE"]] <- 19
      MORFEE.DATA[["GENCODE_URL"]] <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/"
      MORFEE.DATA[["GENCODE_FILE_ANNOT"]] <- "gencode.v19.annotation.gff3.gz" # Position codon start + strand
      MORFEE.DATA[["GENCODE_FILE_METAD"]] <- "gencode.v19.metadata.RefSeq.gz" # NM code
      MORFEE.DATA[["GENCODE_FILE_SEQUE"]] <- "GRCh37.p13.genome.fa.gz" # Whole Human Sequence
    }
  }else if(GRCh==38){
    MORFEE.DATA[["GENCODE"]] <- 43
    MORFEE.DATA[["GENCODE_URL"]] <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/"
    MORFEE.DATA[["GENCODE_FILE_ANNOT"]] <- "gencode.v43.annotation.gff3.gz" # Position codon start + strand
    MORFEE.DATA[["GENCODE_FILE_METAD"]] <- "gencode.v43.metadata.RefSeq.gz" # NM code
    MORFEE.DATA[["GENCODE_FILE_SEQUE"]] <- "GRCh38.primary_assembly.genome.fa.gz" # Whole Human Sequence
  }

  MORFEE.DATA[["SEQ_INIT"]] <- Biostrings::DNAStringSet(c("ATG", "ACG", "CTG","GTG","TTG","AAG","AGG","ATC","ATA","ATT")) # START codon
  MORFEE.DATA[["SEQ_STOP"]] <- Biostrings::DNAStringSet(c("TAA", "TAG", "TGA")) # STOP codon

  if(!dir.exists(MORFEE.DATA[["DB_PATH_GENCODE"]])){
    dir.create(MORFEE.DATA[["DB_PATH_GENCODE"]] , recursive = TRUE)
  }

  if(!file.exists(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_ANNOT"]])) | force){
    download.file(paste0(MORFEE.DATA[["GENCODE_URL"]], MORFEE.DATA[["GENCODE_FILE_ANNOT"]]),
                   paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_ANNOT"]]))
  }

  MORFEE.DATA[["GENCODE_ANNOT"]] <- readGFF(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]],
                                                   MORFEE.DATA[["GENCODE_FILE_ANNOT"]]) , version=3)

  if(!file.exists(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_METAD"]])) | force){
    download.file(paste0(MORFEE.DATA[["GENCODE_URL"]], MORFEE.DATA[["GENCODE_FILE_METAD"]]),
                   paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_METAD"]]))
  }

  MORFEE.DATA[["GENCODE_METAD"]] <- fread(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]],
                                                 MORFEE.DATA[["GENCODE_FILE_METAD"]]))

  if(!file.exists(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_SEQUE"]])) | force){
    download.file(paste0(MORFEE.DATA[["GENCODE_URL"]], MORFEE.DATA[["GENCODE_FILE_SEQUE"]]),
                   paste0(MORFEE.DATA[["DB_PATH_GENCODE"]], MORFEE.DATA[["GENCODE_FILE_SEQUE"]]))
  }

  MORFEE.DATA[["GENCODE_SEQ"]] <- readDNAStringSet(paste0(MORFEE.DATA[["DB_PATH_GENCODE"]],
                                                          MORFEE.DATA[["GENCODE_FILE_SEQUE"]]))

  MORFEE.DATA[["GENCODE_SEQ_ORDER"]] <- str_split_fixed(MORFEE.DATA[["GENCODE_SEQ"]]@ranges@NAMES, " ", n=2)[,1]

  return(MORFEE.DATA)
}

#' write a VCF object as an .ods OpenDocument file
#'
#' @param myvcf a VCF object
#' @param file name of the file to be created
#'
#' @importFrom readODS write_ods
#'
#' @export
#'
vcf_2_ods <- function(myvcf, file){

  df.myvcf <- combine.vcf.slot(myvcf)

  write_ods(x=df.myvcf, path=file)

}


#' write a VCF object as an .xlsx Office Open XML file
#'
#' @param myvcf_TIS a VCF object
#' @param mycvf_STOP a VCF object
#' @param file name of the file to be created
#'
#' @importFrom writexl write_xlsx
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_trunc
#' @importFrom stringi stri_length
#'
#' @export
#'
vcf_2_xlsx <- function(myvcf_TIS, mycvf_STOP,mycvf_STOP_creation, file, hg){

  myvcf_list <- list(myvcf_TIS, mycvf_STOP,mycvf_STOP_creation)
  final_df <- data.frame()

  for(i in 1:length(myvcf_list)){
    myvcf <- myvcf_list[[i]]

    df.myvcf <- combine.vcf.slot(myvcf)

    if ((length(df.myvcf$MORFEE_uTIS) + length(df.myvcf$MORFEE_dTIS) + 
            length(df.myvcf$MORFEE_intTIS) + length(df.myvcf$MORFEE_uSTOP) + 
            length(df.myvcf$MORFEE_dSTOP) + length(df.myvcf$MORFEE_intSTOP) +
            length(df.myvcf$MORFEE_uSTOP_creation) + length(df.myvcf$MORFEE_dSTOP_creation) +
            length(df.myvcf$MORFEE_intSTOP_creation)) == 0) {
            next
        }

    df.myvcf$orfSNVs_type <- NA
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_uTIS)] <- "uTIS"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_dTIS)] <- "dTIS"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_intTIS)] <- "intTIS"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_uSTOP)] <- "uSTOP"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_dSTOP)] <- "dSTOP"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_intSTOP)] <- "intSTOP"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_uSTOP_creation)] <- "new_uSTOP"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_dSTOP_creation)] <- "new_dSTOP"
    df.myvcf$orfSNVs_type[!is.na(df.myvcf$MORFEE_intSTOP_creation)] <- "new_intSTOP"


    enst_uTIS  <- unlist(get.enst(df.myvcf$MORFEE_uTIS))
    enst_dTIS  <- unlist(get.enst(df.myvcf$MORFEE_dTIS))
    enst_intTIS  <- unlist(get.enst(df.myvcf$MORFEE_intTIS))
    enst_uSTOP <- unlist(get.enst(df.myvcf$MORFEE_uSTOP))
    enst_dSTOP <- unlist(get.enst(df.myvcf$MORFEE_dSTOP))
    enst_intSTOP <- unlist(get.enst(df.myvcf$MORFEE_intSTOP))
    enst_new_uSTOP <- unlist(get.enst(df.myvcf$MORFEE_uSTOP_creation))
    enst_new_dSTOP <- unlist(get.enst(df.myvcf$MORFEE_dSTOP_creation))
    enst_new_intSTOP <- unlist(get.enst(df.myvcf$MORFEE_intSTOP_creation))
    enst <- paste(enst_uTIS, enst_dTIS,enst_intTIS, enst_uSTOP, enst_dSTOP, enst_intSTOP, enst_new_uSTOP,enst_new_dSTOP,enst_new_intSTOP,sep=";")
    enst <- gsub(";NA","",enst)
    enst <- gsub("NA;","",enst)
    df.myvcf$Transcript <- enst

    nm_uTIS  <- unlist(get.nm(df.myvcf$MORFEE_uTIS))
    nm_dTIS  <- unlist(get.nm(df.myvcf$MORFEE_dTIS))
    nm_intTIS  <- unlist(get.nm(df.myvcf$MORFEE_intTIS))
    nm_uSTOP <- unlist(get.nm(df.myvcf$MORFEE_uSTOP))
    nm_dSTOP <- unlist(get.nm(df.myvcf$MORFEE_dSTOP))
    nm_intSTOP <- unlist(get.nm(df.myvcf$MORFEE_intSTOP))
    nm_new_uSTOP <- unlist(get.nm(df.myvcf$MORFEE_uSTOP_creation))
    nm_new_dSTOP <- unlist(get.nm(df.myvcf$MORFEE_dSTOP_creation))
    nm_new_intSTOP <- unlist(get.nm(df.myvcf$MORFEE_intSTOP_creation))
    nm <- paste(nm_uTIS, nm_dTIS,nm_intTIS, nm_uSTOP, nm_dSTOP, nm_intSTOP, nm_new_uSTOP,nm_new_dSTOP,nm_new_intSTOP,sep=";")
    nm <- gsub(";NA","",nm)
    nm <- gsub("NA;","",nm)
    df.myvcf$RefSeq_transcript <- nm

    orfSNVs_frame_uTIS  <- unlist(get.orfSNVs(df.myvcf$MORFEE_uTIS, type="frame"))
    orfSNVs_frame_dTIS  <- unlist(get.orfSNVs(df.myvcf$MORFEE_dTIS, type="frame"))
    orfSNVs_frame_intTIS  <- unlist(get.orfSNVs(df.myvcf$MORFEE_intTIS, type="frame"))
    orfSNVs_frame_uSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_uSTOP, type="frame"))
    orfSNVs_frame_dSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_dSTOP, type="frame"))
    orfSNVs_frame_intSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_intSTOP, type="frame"))
    orfSNVs_frame_new_uSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_uSTOP_creation, type="frame"))
    orfSNVs_frame_new_dSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_dSTOP_creation, type="frame"))
    orfSNVs_frame_new_intSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_intSTOP_creation, type="frame"))
    orfSNVs_frame <- paste(orfSNVs_frame_uTIS, orfSNVs_frame_dTIS,orfSNVs_frame_intTIS, orfSNVs_frame_uSTOP, orfSNVs_frame_dSTOP, orfSNVs_frame_intSTOP, orfSNVs_frame_new_uSTOP,orfSNVs_frame_new_dSTOP,orfSNVs_frame_new_intSTOP,sep=";")
    orfSNVs_frame <- gsub(";NA","",orfSNVs_frame)
    orfSNVs_frame <- gsub("NA;","",orfSNVs_frame)
    df.myvcf$orfSNVs_frame <- orfSNVs_frame


    orfSNVs_position_uTIS  <- unlist(get.orfSNVs(df.myvcf$MORFEE_uTIS, type="position"))
    orfSNVs_position_dTIS  <- unlist(get.orfSNVs(df.myvcf$MORFEE_dTIS, type="position"))
    orfSNVs_position_intTIS  <- unlist(get.orfSNVs(df.myvcf$MORFEE_intTIS, type="position"))
    orfSNVs_position_uSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_uSTOP, type="position"))
    orfSNVs_position_dSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_dSTOP, type="position"))
    orfSNVs_position_intSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_intSTOP, type="position"))
    orfSNVs_position_new_uSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_uSTOP_creation, type="position"))
    orfSNVs_position_new_dSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_dSTOP_creation, type="position"))
    orfSNVs_position_new_intSTOP <- unlist(get.orfSNVs(df.myvcf$MORFEE_intSTOP_creation, type="position"))
    orfSNVs_position <- paste(orfSNVs_position_uTIS, orfSNVs_position_dTIS, orfSNVs_position_intTIS, orfSNVs_position_uSTOP, orfSNVs_position_dSTOP, orfSNVs_position_intSTOP, orfSNVs_position_new_uSTOP, orfSNVs_position_new_dSTOP, orfSNVs_position_new_intSTOP, sep=";")
    orfSNVs_position <- gsub(";NA","",orfSNVs_position)
    orfSNVs_position <- gsub("NA;","",orfSNVs_position)
    df.myvcf$type_of_generated_ORF <- orfSNVs_position


    start_seq_uTIS  <- unlist(get.start.seq(df.myvcf$MORFEE_uTIS))
    start_seq_dTIS  <- unlist(get.start.seq(df.myvcf$MORFEE_dTIS))
    start_seq_intTIS  <- unlist(get.start.seq(df.myvcf$MORFEE_intTIS))
    start_seq_uSTOP  <- unlist(get.start.seq(df.myvcf$MORFEE_uSTOP))
    start_seq_dSTOP  <- unlist(get.start.seq(df.myvcf$MORFEE_dSTOP))
    start_seq_intSTOP  <- unlist(get.start.seq(df.myvcf$MORFEE_intSTOP))
    start_seq_new_uSTOP  <- unlist(get.start.seq(df.myvcf$MORFEE_uSTOP_creation))
    start_seq_new_dSTOP  <- unlist(get.start.seq(df.myvcf$MORFEE_dSTOP_creation))
    start_seq_new_intSTOP  <- unlist(get.start.seq(df.myvcf$MORFEE_intSTOP_creation))
    start_seq <- paste(start_seq_uTIS, start_seq_dTIS, start_seq_intTIS, start_seq_uSTOP, start_seq_dSTOP, start_seq_intSTOP,start_seq_new_uSTOP, start_seq_new_dSTOP, start_seq_new_intSTOP, sep=";")
    start_seq <- gsub(";NA","",start_seq)
    start_seq <- gsub("NA;","",start_seq)
    df.myvcf$TIS_sequence <- start_seq


    TIS_type_uTIS <- unlist(get.TIS.type(df.myvcf$MORFEE_uTIS))
    TIS_type_dTIS <- unlist(get.TIS.type(df.myvcf$MORFEE_dTIS))
    TIS_type_intTIS <- unlist(get.TIS.type(df.myvcf$MORFEE_intTIS))
    TIS_type_uSTOP <- unlist(get.TIS.type(df.myvcf$MORFEE_uSTOP))
    TIS_type_dSTOP <- unlist(get.TIS.type(df.myvcf$MORFEE_dSTOP))
    TIS_type_intSTOP <- unlist(get.TIS.type(df.myvcf$MORFEE_intSTOP))
    TIS_type_new_uSTOP <- unlist(get.TIS.type(df.myvcf$MORFEE_uSTOP_creation))
    TIS_type_new_dSTOP <- unlist(get.TIS.type(df.myvcf$MORFEE_dSTOP_creation))
    TIS_type_new_intSTOP <- unlist(get.TIS.type(df.myvcf$MORFEE_intSTOP_creation))
    TIS_type <- paste(TIS_type_uTIS, TIS_type_dTIS, TIS_type_intTIS, TIS_type_uSTOP,TIS_type_dSTOP, TIS_type_intSTOP, TIS_type_new_uSTOP, TIS_type_new_dSTOP, TIS_type_new_intSTOP, sep=";")
    TIS_type <- gsub(";NA","",TIS_type)
    TIS_type <- gsub("NA;","",TIS_type)
    df.myvcf$TIS_type <- TIS_type


    modification_type_uTIS  <- unlist(get.modification.type(df.myvcf$MORFEE_uTIS))
    modification_type_dTIS  <- unlist(get.modification.type(df.myvcf$MORFEE_dTIS))
    modification_type_intTIS  <- unlist(get.modification.type(df.myvcf$MORFEE_intTIS))
    modification_type_uSTOP  <- unlist(get.modification.type(df.myvcf$MORFEE_uSTOP))
    modification_type_dSTOP  <- unlist(get.modification.type(df.myvcf$MORFEE_dSTOP))
    modification_type_intSTOP  <- unlist(get.modification.type(df.myvcf$MORFEE_intSTOP))
    modification_type_new_uSTOP  <- unlist(get.modification.type(df.myvcf$MORFEE_uSTOP_creation))
    modification_type_new_dSTOP  <- unlist(get.modification.type(df.myvcf$MORFEE_dSTOP_creation))
    modification_type_new_intSTOP  <- unlist(get.modification.type(df.myvcf$MORFEE_intSTOP_creation))
    modification_type <- paste(modification_type_uTIS, modification_type_dTIS, modification_type_intTIS, modification_type_uSTOP, modification_type_dSTOP, modification_type_intSTOP, modification_type_new_uSTOP, modification_type_new_dSTOP, modification_type_new_intSTOP, sep=";")
    modification_type <- gsub(";NA","",modification_type)
    modification_type <- gsub("NA;","",modification_type)
    df.myvcf$modification_type <- modification_type


    NewAALength_uTIS <- unlist(get.NewAALength(df.myvcf$MORFEE_uTIS))
    NewAALength_dTIS <- unlist(get.NewAALength(df.myvcf$MORFEE_dTIS))
    NewAALength_intTIS <- unlist(get.NewAALength(df.myvcf$MORFEE_intTIS))
    NewAALength_uSTOP <- unlist(get.NewAALength(df.myvcf$MORFEE_uSTOP))
    NewAALength_dSTOP <- unlist(get.NewAALength(df.myvcf$MORFEE_dSTOP))
    NewAALength_intSTOP <- unlist(get.NewAALength(df.myvcf$MORFEE_intSTOP))
    NewAALength_new_uSTOP <- unlist(get.NewAALength(df.myvcf$MORFEE_uSTOP_creation))
    NewAALength_new_dSTOP <- unlist(get.NewAALength(df.myvcf$MORFEE_dSTOP_creation))
    NewAALength_new_intSTOP <- unlist(get.NewAALength(df.myvcf$MORFEE_intSTOP_creation))
    NewAALength <- paste(NewAALength_uTIS, NewAALength_dTIS, NewAALength_intTIS, NewAALength_uSTOP, NewAALength_dSTOP, NewAALength_intSTOP, NewAALength_new_uSTOP, NewAALength_new_dSTOP, NewAALength_new_intSTOP, sep=";")
    NewAALength <- gsub(";NA","",NewAALength)
    NewAALength <- gsub("NA;","",NewAALength)
    df.myvcf$NewAALength <- NewAALength


    StopRelPos_uTIS <- unlist(get.StopRelPos(df.myvcf$MORFEE_uTIS))
    StopRelPos_dTIS <- unlist(get.StopRelPos(df.myvcf$MORFEE_dTIS))
    StopRelPos_intTIS <- unlist(get.StopRelPos(df.myvcf$MORFEE_intTIS))
    StopRelPos_uSTOP <- unlist(get.StopRelPos(df.myvcf$MORFEE_uSTOP))
    StopRelPos_dSTOP <- unlist(get.StopRelPos(df.myvcf$MORFEE_dSTOP))
    StopRelPos_intSTOP <- unlist(get.StopRelPos(df.myvcf$MORFEE_intSTOP))
    StopRelPos_new_uSTOP <- unlist(get.StopRelPos(df.myvcf$MORFEE_uSTOP_creation))
    StopRelPos_new_dSTOP <- unlist(get.StopRelPos(df.myvcf$MORFEE_dSTOP_creation))
    StopRelPos_new_intSTOP <- unlist(get.StopRelPos(df.myvcf$MORFEE_intSTOP_creation))
    StopRelPos <- paste(StopRelPos_uTIS, StopRelPos_dTIS, StopRelPos_intTIS, StopRelPos_uSTOP, StopRelPos_dSTOP, StopRelPos_intSTOP, StopRelPos_new_uSTOP, StopRelPos_new_dSTOP, StopRelPos_new_intSTOP, sep=";")
    StopRelPos <- gsub(";NA","",StopRelPos)
    StopRelPos <- gsub("NA;","",StopRelPos)
    df.myvcf$ORF_size <- StopRelPos


    fixed_annot_uTIS <- unlist(get.FixedAnnot(df.myvcf$MORFEE_uTIS))
    fixed_annot_dTIS <- unlist(get.FixedAnnot(df.myvcf$MORFEE_dTIS))
    fixed_annot_intTIS <- unlist(get.FixedAnnot(df.myvcf$MORFEE_intTIS))
    fixed_annot_uSTOP <- unlist(get.FixedAnnot(df.myvcf$MORFEE_uSTOP))
    fixed_annot_dSTOP <- unlist(get.FixedAnnot(df.myvcf$MORFEE_dSTOP))
    fixed_annot_intSTOP <- unlist(get.FixedAnnot(df.myvcf$MORFEE_intSTOP))
    fixed_annot_new_uSTOP <- unlist(get.FixedAnnot(df.myvcf$MORFEE_uSTOP_creation))
    fixed_annot_new_dSTOP <- unlist(get.FixedAnnot(df.myvcf$MORFEE_dSTOP_creation))
    fixed_annot_new_intSTOP <- unlist(get.FixedAnnot(df.myvcf$MORFEE_intSTOP_creation))
    fixed_annot <- paste(fixed_annot_uTIS, fixed_annot_dTIS, fixed_annot_intTIS, fixed_annot_uSTOP, fixed_annot_dSTOP, fixed_annot_intSTOP, fixed_annot_new_uSTOP, fixed_annot_new_dSTOP, fixed_annot_new_intSTOP, sep=";")
    fixed_annot <- gsub(";NA","",fixed_annot)
    fixed_annot <- gsub("NA;","",fixed_annot)
    df.myvcf$Variation <- fixed_annot

    TIS_geno_pos_uTIS <- unlist(get.TIS_Start(df.myvcf$MORFEE_uTIS))
    TIS_geno_pos_dTIS <- unlist(get.TIS_Start(df.myvcf$MORFEE_dTIS))
    TIS_geno_pos_intTIS <- unlist(get.TIS_Start(df.myvcf$MORFEE_intTIS))
    TIS_geno_pos_uSTOP <- unlist(get.TIS_Start(df.myvcf$MORFEE_uSTOP))
    TIS_geno_pos_dSTOP <- unlist(get.TIS_Start(df.myvcf$MORFEE_dSTOP))
    TIS_geno_pos_intSTOP <- unlist(get.TIS_Start(df.myvcf$MORFEE_intSTOP))
    TIS_geno_pos_new_uSTOP <- unlist(get.TIS_Start(df.myvcf$MORFEE_uSTOP_creation))
    TIS_geno_pos_new_dSTOP <- unlist(get.TIS_Start(df.myvcf$MORFEE_dSTOP_creation))
    TIS_geno_pos_new_intSTOP <- unlist(get.TIS_Start(df.myvcf$MORFEE_intSTOP_creation))
    TIS_geno_pos <- paste(TIS_geno_pos_uTIS, TIS_geno_pos_dTIS, TIS_geno_pos_intTIS, TIS_geno_pos_uSTOP, TIS_geno_pos_dSTOP, TIS_geno_pos_intSTOP, TIS_geno_pos_new_uSTOP, TIS_geno_pos_new_dSTOP, TIS_geno_pos_new_intSTOP, sep=";")
    TIS_geno_pos <- gsub(";NA","",TIS_geno_pos)
    TIS_geno_pos <- gsub("NA;","",TIS_geno_pos)
    df.myvcf$TIS_position <- TIS_geno_pos

    STOP_geno_pos_uTIS <- unlist(get.StopGenoPos(df.myvcf$MORFEE_uTIS))
    STOP_geno_pos_dTIS <- unlist(get.StopGenoPos(df.myvcf$MORFEE_dTIS))
    STOP_geno_pos_intTIS <- unlist(get.StopGenoPos(df.myvcf$MORFEE_intTIS))
    STOP_geno_pos_uSTOP <- unlist(get.StopGenoPos(df.myvcf$MORFEE_uSTOP))
    STOP_geno_pos_dSTOP <- unlist(get.StopGenoPos(df.myvcf$MORFEE_dSTOP))
    STOP_geno_pos_intSTOP <- unlist(get.StopGenoPos(df.myvcf$MORFEE_intSTOP))
    STOP_geno_pos_new_uSTOP <- unlist(get.StopGenoPos(df.myvcf$MORFEE_uSTOP_creation))
    STOP_geno_pos_new_dSTOP <- unlist(get.StopGenoPos(df.myvcf$MORFEE_dSTOP_creation))
    STOP_geno_pos_new_intSTOP <- unlist(get.StopGenoPos(df.myvcf$MORFEE_intSTOP_creation))
    STOP_geno_pos <- paste(STOP_geno_pos_uTIS, STOP_geno_pos_dTIS, STOP_geno_pos_intTIS, STOP_geno_pos_uSTOP, STOP_geno_pos_dSTOP, STOP_geno_pos_intSTOP, STOP_geno_pos_new_uSTOP, STOP_geno_pos_new_dSTOP, STOP_geno_pos_new_intSTOP, sep=";")
    STOP_geno_pos <- gsub(";NA","",STOP_geno_pos)
    STOP_geno_pos <- gsub("NA;","",STOP_geno_pos)
    df.myvcf$STOP_position <- STOP_geno_pos

    CDS_start_pos_uTIS <- unlist(get.CDS_Start(df.myvcf$MORFEE_uTIS))
    CDS_start_pos_dTIS <- unlist(get.CDS_Start(df.myvcf$MORFEE_dTIS))
    CDS_start_pos_intTIS <- unlist(get.CDS_Start(df.myvcf$MORFEE_intTIS))
    CDS_start_pos_uSTOP <- unlist(get.CDS_Start(df.myvcf$MORFEE_uSTOP))
    CDS_start_pos_dSTOP <- unlist(get.CDS_Start(df.myvcf$MORFEE_dSTOP))
    CDS_start_pos_intSTOP <- unlist(get.CDS_Start(df.myvcf$MORFEE_intSTOP))
    CDS_start_pos_new_uSTOP <- unlist(get.CDS_Start(df.myvcf$MORFEE_uSTOP_creation))
    CDS_start_pos_new_dSTOP <- unlist(get.CDS_Start(df.myvcf$MORFEE_dSTOP_creation))
    CDS_start_pos_new_intSTOP <- unlist(get.CDS_Start(df.myvcf$MORFEE_intSTOP_creation))
    CDS_start_pos <- paste(CDS_start_pos_uTIS, CDS_start_pos_dTIS, CDS_start_pos_intTIS, CDS_start_pos_uSTOP, CDS_start_pos_dSTOP, CDS_start_pos_intSTOP, CDS_start_pos_new_uSTOP, CDS_start_pos_new_dSTOP, CDS_start_pos_new_intSTOP, sep=";")
    CDS_start_pos <- gsub(";NA","",CDS_start_pos)
    CDS_start_pos <- gsub("NA;","",CDS_start_pos)
    df.myvcf$CDS_start_position <- CDS_start_pos

    CDS_end_pos_uTIS <- unlist(get.CDS_End(df.myvcf$MORFEE_uTIS))
    CDS_end_pos_dTIS <- unlist(get.CDS_End(df.myvcf$MORFEE_dTIS))
    CDS_end_pos_intTIS <- unlist(get.CDS_End(df.myvcf$MORFEE_intTIS))
    CDS_end_pos_uSTOP <- unlist(get.CDS_End(df.myvcf$MORFEE_uSTOP))
    CDS_end_pos_dSTOP <- unlist(get.CDS_End(df.myvcf$MORFEE_dSTOP))
    CDS_end_pos_intSTOP <- unlist(get.CDS_End(df.myvcf$MORFEE_intSTOP))
    CDS_end_pos_new_uSTOP <- unlist(get.CDS_End(df.myvcf$MORFEE_uSTOP_creation))
    CDS_end_pos_new_dSTOP <- unlist(get.CDS_End(df.myvcf$MORFEE_dSTOP_creation))
    CDS_end_pos_new_intSTOP <- unlist(get.CDS_End(df.myvcf$MORFEE_intSTOP_creation))
    CDS_end_pos <- paste(CDS_end_pos_uTIS, CDS_end_pos_dTIS, CDS_end_pos_intTIS, CDS_end_pos_uSTOP, CDS_end_pos_dSTOP, CDS_end_pos_intSTOP, CDS_end_pos_new_uSTOP, CDS_end_pos_new_dSTOP, CDS_end_pos_new_intSTOP, sep=";")
    CDS_end_pos <- gsub(";NA","",CDS_end_pos)
    CDS_end_pos <- gsub("NA;","",CDS_end_pos)
    df.myvcf$CDS_end_position <- CDS_end_pos

    Kozak_score_uTIS <- unlist(get.kozak.score(df.myvcf$MORFEE_uTIS))
    Kozak_score_dTIS <- unlist(get.kozak.score(df.myvcf$MORFEE_dTIS))
    Kozak_score_intTIS <- unlist(get.kozak.score(df.myvcf$MORFEE_intTIS))
    Kozak_score_uSTOP <- unlist(get.kozak.score(df.myvcf$MORFEE_uSTOP))
    Kozak_score_dSTOP <- unlist(get.kozak.score(df.myvcf$MORFEE_dSTOP))
    Kozak_score_intSTOP <- unlist(get.kozak.score(df.myvcf$MORFEE_intSTOP))
    Kozak_score_new_uSTOP <- unlist(get.kozak.score(df.myvcf$MORFEE_uSTOP_creation))
    Kozak_score_new_dSTOP <- unlist(get.kozak.score(df.myvcf$MORFEE_dSTOP_creation))
    Kozak_score_new_intSTOP <- unlist(get.kozak.score(df.myvcf$MORFEE_intSTOP_creation))
    Kozak_score <- paste(Kozak_score_uTIS, Kozak_score_dTIS, Kozak_score_intTIS, Kozak_score_uSTOP, Kozak_score_dSTOP, Kozak_score_intSTOP, Kozak_score_new_uSTOP, Kozak_score_new_dSTOP, Kozak_score_new_intSTOP, sep=";")
    Kozak_score <- gsub(";NA","",Kozak_score)
    Kozak_score <- gsub("NA;","",Kozak_score)
    df.myvcf$Kozak_score <- Kozak_score

    Kozak_sequence_uTIS <- unlist(get.kozak.sequence(df.myvcf$MORFEE_uTIS))
    Kozak_sequence_dTIS <- unlist(get.kozak.sequence(df.myvcf$MORFEE_dTIS))
    Kozak_sequence_intTIS <- unlist(get.kozak.sequence(df.myvcf$MORFEE_intTIS))
    Kozak_sequence_uSTOP <- unlist(get.kozak.sequence(df.myvcf$MORFEE_uSTOP))
    Kozak_sequence_dSTOP <- unlist(get.kozak.sequence(df.myvcf$MORFEE_dSTOP))
    Kozak_sequence_intSTOP <- unlist(get.kozak.sequence(df.myvcf$MORFEE_intSTOP))
    Kozak_sequence_new_uSTOP <- unlist(get.kozak.sequence(df.myvcf$MORFEE_uSTOP_creation))
    Kozak_sequence_new_dSTOP <- unlist(get.kozak.sequence(df.myvcf$MORFEE_dSTOP_creation))
    Kozak_sequence_new_intSTOP <- unlist(get.kozak.sequence(df.myvcf$MORFEE_intSTOP_creation))
    Kozak_sequence <- paste(Kozak_sequence_uTIS, Kozak_sequence_dTIS, Kozak_sequence_intTIS, Kozak_sequence_uSTOP, Kozak_sequence_dSTOP, Kozak_sequence_intSTOP, Kozak_sequence_new_uSTOP, Kozak_sequence_new_dSTOP, Kozak_sequence_new_intSTOP, sep=";")
    Kozak_sequence <- gsub(";NA","",Kozak_sequence)
    Kozak_sequence <- gsub("NA;","",Kozak_sequence)
    df.myvcf$Kozak_sequence <- Kozak_sequence

    Kozak_strength_uTIS <- unlist(get.kozak.strength(df.myvcf$MORFEE_uTIS))
    Kozak_strength_dTIS <- unlist(get.kozak.strength(df.myvcf$MORFEE_dTIS))
    Kozak_strength_intTIS <- unlist(get.kozak.strength(df.myvcf$MORFEE_intTIS))
    Kozak_strength_uSTOP <- unlist(get.kozak.strength(df.myvcf$MORFEE_uSTOP))
    Kozak_strength_dSTOP <- unlist(get.kozak.strength(df.myvcf$MORFEE_dSTOP))
    Kozak_strength_intSTOP <- unlist(get.kozak.strength(df.myvcf$MORFEE_intSTOP))
    Kozak_strength_new_uSTOP <- unlist(get.kozak.strength(df.myvcf$MORFEE_uSTOP_creation))
    Kozak_strength_new_dSTOP <- unlist(get.kozak.strength(df.myvcf$MORFEE_dSTOP_creation))
    Kozak_strength_new_intSTOP <- unlist(get.kozak.strength(df.myvcf$MORFEE_intSTOP_creation))
    Kozak_strength <- paste(Kozak_strength_uTIS, Kozak_strength_dTIS, Kozak_strength_intTIS, Kozak_strength_uSTOP, Kozak_strength_dSTOP, Kozak_strength_intSTOP, Kozak_strength_new_uSTOP, Kozak_strength_new_dSTOP, Kozak_strength_new_intSTOP, sep=";")
    Kozak_strength <- gsub(";NA","",Kozak_strength)
    Kozak_strength <- gsub("NA;","",Kozak_strength)
    df.myvcf$Kozak_strength <- Kozak_strength

    TIS_c_pos_uTIS <- unlist(get.tis.c.pos(df.myvcf$MORFEE_uTIS))
    TIS_c_pos_dTIS <- unlist(get.tis.c.pos(df.myvcf$MORFEE_dTIS))
    TIS_c_pos_intTIS <- unlist(get.tis.c.pos(df.myvcf$MORFEE_intTIS))
    TIS_c_pos_uSTOP <- unlist(get.tis.c.pos(df.myvcf$MORFEE_uSTOP))
    TIS_c_pos_dSTOP <- unlist(get.tis.c.pos(df.myvcf$MORFEE_dSTOP))
    TIS_c_pos_intSTOP <- unlist(get.tis.c.pos(df.myvcf$MORFEE_intSTOP))
    TIS_c_pos_new_uSTOP <- unlist(get.tis.c.pos(df.myvcf$MORFEE_uSTOP_creation))
    TIS_c_pos_new_dSTOP <- unlist(get.tis.c.pos(df.myvcf$MORFEE_dSTOP_creation))
    TIS_c_pos_new_intSTOP <- unlist(get.tis.c.pos(df.myvcf$MORFEE_intSTOP_creation))
    TIS_c_pos <- paste(TIS_c_pos_uTIS, TIS_c_pos_dTIS, TIS_c_pos_intTIS, TIS_c_pos_uSTOP, TIS_c_pos_dSTOP, TIS_c_pos_intSTOP, TIS_c_pos_new_uSTOP, TIS_c_pos_new_dSTOP, TIS_c_pos_new_intSTOP, sep=";")
    TIS_c_pos <- gsub(";NA","",TIS_c_pos)
    TIS_c_pos <- gsub("NA;","",TIS_c_pos)
    df.myvcf$TIS_c_pos <- TIS_c_pos

    STOP_c_pos_uTIS <- unlist(get.stop.c.pos(df.myvcf$MORFEE_uTIS))
    STOP_c_pos_dTIS <- unlist(get.stop.c.pos(df.myvcf$MORFEE_dTIS))
    STOP_c_pos_intTIS <- unlist(get.stop.c.pos(df.myvcf$MORFEE_intTIS))
    STOP_c_pos_uSTOP <- unlist(get.stop.c.pos(df.myvcf$MORFEE_uSTOP))
    STOP_c_pos_dSTOP <- unlist(get.stop.c.pos(df.myvcf$MORFEE_dSTOP))
    STOP_c_pos_intSTOP <- unlist(get.stop.c.pos(df.myvcf$MORFEE_intSTOP))
    STOP_c_pos_new_uSTOP <- unlist(get.stop.c.pos(df.myvcf$MORFEE_uSTOP_creation))
    STOP_c_pos_new_dSTOP <- unlist(get.stop.c.pos(df.myvcf$MORFEE_dSTOP_creation))
    STOP_c_pos_new_intSTOP <- unlist(get.stop.c.pos(df.myvcf$MORFEE_intSTOP_creation))
    STOP_c_pos <- paste(STOP_c_pos_uTIS, STOP_c_pos_dTIS, STOP_c_pos_intTIS, STOP_c_pos_uSTOP, STOP_c_pos_dSTOP, STOP_c_pos_intSTOP, STOP_c_pos_new_uSTOP, STOP_c_pos_new_dSTOP, STOP_c_pos_new_intSTOP, sep=";")
    STOP_c_pos <- gsub(";NA","",STOP_c_pos)
    STOP_c_pos <- gsub("NA;","",STOP_c_pos)
    df.myvcf$STOP_c_pos <- STOP_c_pos

    STOP_seq_uTIS <- unlist(get.stop.seq(df.myvcf$MORFEE_uTIS))
    STOP_seq_dTIS <- unlist(get.stop.seq(df.myvcf$MORFEE_dTIS))
    STOP_seq_intTIS <- unlist(get.stop.seq(df.myvcf$MORFEE_intTIS))
    STOP_seq_uSTOP <- unlist(get.stop.seq(df.myvcf$MORFEE_uSTOP))
    STOP_seq_dSTOP <- unlist(get.stop.seq(df.myvcf$MORFEE_dSTOP))
    STOP_seq_intSTOP <- unlist(get.stop.seq(df.myvcf$MORFEE_intSTOP))
    STOP_seq_new_uSTOP <- unlist(get.stop.seq(df.myvcf$MORFEE_uSTOP_creation))
    STOP_seq_new_dSTOP <- unlist(get.stop.seq(df.myvcf$MORFEE_dSTOP_creation))
    STOP_seq_new_intSTOP <- unlist(get.stop.seq(df.myvcf$MORFEE_intSTOP_creation))
    STOP_seq <- paste(STOP_seq_uTIS, STOP_seq_dTIS, STOP_seq_intTIS, STOP_seq_uSTOP, STOP_seq_dSTOP, STOP_seq_intSTOP, STOP_seq_new_uSTOP, STOP_seq_new_dSTOP, STOP_seq_new_intSTOP, sep=";")
    STOP_seq <- gsub(";NA","",STOP_seq)
    STOP_seq <- gsub("NA;","",STOP_seq)
    df.myvcf$STOP_sequence <- STOP_seq

    df.myvcf$Ratio_length_pred_obs <- get.Ratio_length_pred_obs(df.myvcf$NewAALength)


    if(hg == 19){
      col2keep <- c("seqnames","start","REF","ALT","Transcript","RefSeq_transcript","Gene.refGene","avsnp150","Variation","Func.ensGene",
                  "orfSNVs_type","orfSNVs_frame","type_of_generated_ORF",
                  "Ratio_length_pred_obs","NewAALength","ORF_size","MORFEE_uTIS","MORFEE_dTIS","MORFEE_intTIS","MORFEE_uSTOP","MORFEE_dSTOP","MORFEE_intSTOP","MORFEE_uSTOP_creation","MORFEE_dSTOP_creation","MORFEE_intSTOP_creation","TIS_sequence", "TIS_type","modification_type","STOP_sequence","CDS_start_position","CDS_end_position","TIS_position","STOP_position","Kozak_sequence","Kozak_score","Kozak_strength","TIS_c_pos","STOP_c_pos",
                  "pLI.refGene","exp_lof.refGene","oe_lof.refGene","oe_lof_lower.refGene","oe_lof_upper.refGene",
                  "gwasCatalog","CLNDN","CLNDISDB","CLNSIG","AF","AF_popmax",
                  "SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred",
                  "MutationAssessor_pred","FATHMM_pred","PROVEAN_pred","MetaSVM_pred","MetaLR_pred",
                  "M.CAP_pred","REVEL_rankscore","MutPred_rankscore","CADD_phred","fathmm.MKL_coding_pred",
                  "integrated_fitCons_score","GERP.._RS","phyloP100way_vertebrate_rankscore","phyloP30way_mammalian_rankscore",
                  "phastCons100way_vertebrate_rankscore","phastCons30way_mammalian_rankscore",
                  "SiPhy_29way_logOdds")# M-CAP_pred, fathmm-MKL_coding_pred, GERP++_RS

    }else if(hg == 38){
        col2keep <- c("seqnames","start","REF","ALT","Transcript","RefSeq_transcript","Gene.ensGene","avsnp150","Variation","Func.ensGene",
                "orfSNVs_type","orfSNVs_frame","type_of_generated_ORF",
                "Ratio_length_pred_obs","NewAALength","ORF_size","MORFEE_uTIS","MORFEE_dTIS","MORFEE_intTIS","MORFEE_uSTOP","MORFEE_dSTOP","MORFEE_intSTOP","MORFEE_uSTOP_creation","MORFEE_dSTOP_creation","MORFEE_intSTOP_creation","TIS_sequence", "TIS_type","modification_type","STOP_sequence","CDS_start_position","CDS_end_position","TIS_position","STOP_position","Kozak_sequence","Kozak_score","Kozak_strength","TIS_c_pos","STOP_c_pos",
                "pLI.refGene","exp_lof.refGene","oe_lof.refGene","oe_lof_lower.refGene","oe_lof_upper.refGene",
                "gwasCatalog","CLNDN","CLNDISDB","CLNSIG","gnomad312_AF","gnomad312_AF_popmax",
                "SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","LRT_pred","MutationTaster_pred",
                "MutationAssessor_pred","FATHMM_pred","PROVEAN_pred","MetaSVM_pred","MetaLR_pred",
                "M.CAP_pred","REVEL_rankscore","MutPred_rankscore","CADD_phred","fathmm.MKL_coding_pred",
                "integrated_fitCons_score","GERP.._RS","phyloP100way_vertebrate_rankscore","phyloP30way_mammalian_rankscore",
                "phastCons100way_vertebrate_rankscore","phastCons30way_mammalian_rankscore",
                "SiPhy_29way_logOdds")# M-CAP_pred, fathmm-MKL_coding_pred, GERP++_RS
    }

    df.ok <- as.data.frame(df.myvcf[,col2keep])
    final_df <- rbind(final_df,df.ok)


}


  for(i in 1:ncol(final_df)){
    temp_aschar <- lapply(final_df[,i], as.character)
    temp_paste <- lapply(temp_aschar, paste, sep="", collapse=",")
    temp_unlist <- unlist(temp_paste)
    final_df[,i] <- gsub("\\\\x3b"," ; ",temp_unlist) # Fix encoding problem
    final_df[,i] <- gsub("Name\\\\x3d","",temp_unlist) # Clean gwasCatalog
  }

  # Workaround the "string exceeds Excel's limit of 32,767 characters"
  df.str.length <- apply(final_df,2,stri_length)
  xlsx_lim <- 32000
  final_df[df.str.length > xlsx_lim] <- str_trunc(final_df[df.str.length > xlsx_lim], xlsx_lim, "right")
  
  write_xlsx(x=final_df, path=file)

}

get.Ratio_length_pred_obs <- function(x){

  z <- sapply(x, get.Ratio_length_pred_obs.meth)

  return(z)
}

get.Ratio_length_pred_obs.meth <- function(x){

  if(is.na(x)){return(NA)}

  x <- str_split_fixed(x,";", n=Inf)
  y <- sapply(x, function(x) eval(parse(text=x)))
  y <- as.numeric(y)
  y.min.id <- which.min(abs(y-1))[1]
  z <- y[y.min.id]

  return(z)
}

get.NewAALength <- function(x){

  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))

  orf <- lapply(morfee_temp, get.NewAALength.meth)

  return(orf)
}

get.NewAALength.meth <- function(y){

  if(ncol(y)==22){  
    length_temp <- paste(y[,10], collapse="; ", sep="")
    length_temp <- gsub("\\(aa\\)","",length_temp)
    length_temp <- gsub("\\[","",length_temp)
    length_temp <- gsub("\\]","",length_temp)
    return(length_temp)
  }else if(ncol(y)==16){
    length_temp <- paste(y[,7], collapse="; ", sep="")
    length_temp <- gsub("\\(aa\\)","",length_temp)
    length_temp <- gsub("\\[","",length_temp)
    length_temp <- gsub("\\]","",length_temp)
    return(length_temp)

  }else{
    return(NA)
  }
}

get.orfSNVs <- function(x, type=c("frame","position")){

  if(!(type %in% c("frame","position"))){
    stop("'type' is not 'frame' nor 'position'")
  }

  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))

  orf <- lapply(morfee_temp, get.orfSNVs.meth, type)

  return(orf)
}

get.orfSNVs.meth <- function(y, type){

  if(ncol(y)==22){   
    if(type=="frame"){
      return(paste(y[,4], collapse="; ", sep=""))
    }else if(type=="position"){
      return(paste(y[,9], collapse="; ", sep=""))
    }
  }else if(ncol(y)==16){
      if(type=="frame"){
      return(paste(y[,4], collapse="; ", sep=""))
    }else if(type=="position"){
      return(paste(y[,6], collapse="; ", sep=""))
    }
  }else{
    return(NA)
  }
}

get.start.seq <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.start.seq.meth)
}

get.start.seq.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,5], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}


get.modification.type <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.modification.type.meth)
}

get.modification.type.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,8], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.TIS.type <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.TIS.type.meth)
}

get.TIS.type.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,7], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.StopRelPos <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.StopRelPos.meth)
}

get.StopRelPos.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,11], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,10], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.FixedAnnot <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.FixedAnnot.meth)
}

get.FixedAnnot.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,2], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,2], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.StopGenoPos <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.StopGenoPos.meth)
}

get.StopGenoPos.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,14], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,9], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.CDS_Start <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.CDS_Start.meth)
}

get.CDS_Start.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,12], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,11], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.CDS_End <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.CDS_End.meth)
}

get.CDS_End.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,15], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,12], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.TIS_Start <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.TIS_Start.meth)
}

get.TIS_Start.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,13], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,8], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}



get.enst <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.enst.meth)
}

get.enst.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,1], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,1], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.kozak.score <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.kozak.score.meth)
}

get.kozak.score.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,17], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.kozak.sequence <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.kozak.sequence.meth)
}

get.kozak.sequence.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,16], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.kozak.strength <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.kozak.strength.meth)
}

get.kozak.strength.meth <- function(y){
  if(ncol(y)==22){  
    return(paste(y[,20], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.tis.c.pos <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.tis.c.pos.meth)
}

get.tis.c.pos.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,18], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,13], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.stop.c.pos <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.stop.c.pos.meth)
}

get.stop.c.pos.meth <- function(y){
  if(ncol(y)==22){   
    return(paste(y[,19], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,14], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.stop.seq <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.stop.seq.meth)
}

get.stop.seq.meth <- function(y){
  if(ncol(y)==22){  
    return(paste(y[,21], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,15], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}

get.nm <- function(x){
  morfee_temp <- lapply(x, function(x) str_split_fixed(x,"\\|", n=Inf)[1,])
  morfee_temp <- lapply(morfee_temp, function(x) str_split_fixed(x,",", n=Inf))
  
  orf <- lapply(morfee_temp, get.nm.meth)
}

get.nm.meth <- function(y){
  if(ncol(y)==22){  
    return(paste(y[,22], collapse="; ", sep=""))
  }else if(ncol(y)==16){
    return(paste(y[,16], collapse="; ", sep=""))
  }else{
    return(NA)
  }
}




#' combine different slots of VCF object
#'
#' @param myvcf a VCF object
#'
#' @return a data.frame
#'
#' @importFrom DelayedArray rowRanges
#'
#' @export
#'
combine.vcf.slot <- function(myvcf){

  if(class(myvcf)!="CollapsedVCF"){
    stop("Class of myvcf must be 'CollapsedVCF'!!!")
  }

  df <- as.data.frame(rowRanges(myvcf))
  df.info <- as.data.frame(info(myvcf))
  df <- cbind(df,df.info)

  df.geno <- as.data.frame(geno(myvcf))
  if(nrow(df.geno)!=0){
    df <- cbind(df,df.geno[df.geno$group_name=="GT",])
  }

 return(df)
}