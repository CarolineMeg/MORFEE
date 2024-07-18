#' vcf annotation with MORFEE
#'
#' @param myvcf_annot an ANNOVAR annotated VCF object to annotate with MORFEE
#' @param MORFEE_DATA data object obtained from get.morfee.data()
#'
#' @importFrom stats na.omit
#' @import foreach
#' @import VariantAnnotation
#' @import Biostrings
#' @import GenomicRanges
#'
#' @export
#'
morfee.annotation <- function(myvcf_annot, MORFEE_DATA){
  myvcf_annot_info <- info(myvcf_annot)
  myvcf_annot_info$MORFEE_uTIS <- NA
  myvcf_annot_info$MORFEE_dTIS <- NA
  myvcf_annot_info$MORFEE_intTIS <- NA
  myvcf_annot_info$MORFEE_uSTOP <- NA
  myvcf_annot_info$MORFEE_dSTOP <- NA
  myvcf_annot_info$MORFEE_intSTOP <- NA
  myvcf_annot_info$MORFEE_uSTOP_creation <- NA
  myvcf_annot_info$MORFEE_dSTOP_creation <- NA
  myvcf_annot_info$MORFEE_intSTOP_creation <- NA
  myvcf_annot_header <- rbind(info(header(myvcf_annot)),
                              data.frame(Number = ".", Type = "String", Description = "New upstream TIS annotation provided by MORFEE", stringsAsFactors = FALSE),
                              data.frame(Number = ".", Type = "String", Description = "New downstream TIS annotation provided by MORFEE", stringsAsFactors = FALSE),
                              data.frame(Number = ".", Type = "String", Description = "New exonic TIS annotation provided by MORFEE", stringsAsFactors = FALSE),
                              data.frame(Number = ".", Type = "String", Description = "Deletion upstream STOP annotation provided by MORFEE", stringsAsFactors = FALSE),
                              data.frame(Number = ".", Type = "String", Description = "Deletion downstream STOP annotation provided by MORFEE", stringsAsFactors = FALSE), 
                              data.frame(Number = ".", Type = "String", Description = "Deletion exonic STOP annotation provided by MORFEE", stringsAsFactors = FALSE), 
                              data.frame(Number = ".", Type = "String", Description = "New upstream STOP annotation provided by MORFEE", stringsAsFactors = FALSE),
                              data.frame(Number = ".", Type = "String", Description = "New downstream STOP annotation provided by MORFEE", stringsAsFactors = FALSE), 
                              data.frame(Number = ".", Type = "String", Description = "New exonic STOP annotation provided by MORFEE", stringsAsFactors = FALSE) )
                              
  rownames(myvcf_annot_header) <- c(rownames(info(header(myvcf_annot))),"MORFEE_uTIS", "MORFEE_dTIS", "MORFEE_intTIS","MORFEE_uSTOP", "MORFEE_dSTOP", "MORFEE_intSTOP","MORFEE_uSTOP_creation","MORFEE_dSTOP_creation","MORFEE_intSTOP_creation")
  info(header(myvcf_annot)) <- myvcf_annot_header
  info(myvcf_annot) <- myvcf_annot_info

  myvcf_annot_info_STOP <- info(myvcf_annot)
  myvcf_annot_info_STOP$MORFEE_uTIS <- NA
  myvcf_annot_info_STOP$MORFEE_dTIS <- NA
  myvcf_annot_info_STOP$MORFEE_intTIS <- NA
  myvcf_annot_info_STOP$MORFEE_uSTOP <- NA
  myvcf_annot_info_STOP$MORFEE_dSTOP <- NA
  myvcf_annot_info_STOP$MORFEE_intSTOP <- NA
  myvcf_annot_info$MORFEE_uSTOP_creation <- NA
  myvcf_annot_info$MORFEE_dSTOP_creation <- NA
  myvcf_annot_info$MORFEE_intSTOP_creation <- NA
  info(myvcf_annot) <- myvcf_annot_info_STOP


  myvcf_annot_info_STOP_creation <- info(myvcf_annot)
  myvcf_annot_info_STOP_creation$MORFEE_uTIS <- NA
  myvcf_annot_info_STOP_creation$MORFEE_dTIS <- NA
  myvcf_annot_info_STOP_creation$MORFEE_intTIS <- NA
  myvcf_annot_info_STOP_creation$MORFEE_uSTOP <- NA
  myvcf_annot_info_STOP_creation$MORFEE_dSTOP <- NA
  myvcf_annot_info_STOP_creation$MORFEE_intSTOP <- NA
  myvcf_annot_info$MORFEE_uSTOP_creation <- NA
  myvcf_annot_info$MORFEE_dSTOP_creation <- NA
  myvcf_annot_info$MORFEE_intSTOP_creation <- NA
  info(myvcf_annot) <- myvcf_annot_info_STOP_creation

  if(MORFEE_DATA[["GRCh"]] == 38){
    MORFEE_DATA[["GENCODE_METAD"]] = rbind(MORFEE_DATA[["GENCODE_METAD"]][1,],MORFEE_DATA[["GENCODE_METAD"]])
    colnames(MORFEE_DATA[["GENCODE_METAD"]]) = c("V1","V2","V3")
  }

  i <- NULL
  for(i in 1:nrow(myvcf_annot_info)){
    my_func <- as.character(myvcf_annot_info[i,"Func.ensGene"])
    if(!(my_func %in% c("UTR5","UTR3","exonic"))){
      message(paste0("Skip variant ",i,": not in UTR5, UTR3 or exonic region"))
      next
    }

    if(my_func == "exonic"){
      my_refgene <- as.character(myvcf_annot_info[i,"AAChange.refGene"][[1]])
      my_ensgene <- as.character(myvcf_annot_info[i,"AAChange.ensGene"][[1]])
    }else{
      my_refgene <- as.character(myvcf_annot_info[i,"GeneDetail.refGene"])
      my_ensgene <- as.character(myvcf_annot_info[i,"GeneDetail.ensGene"])
    }

    if(my_ensgene[[1]]=="."){
      message(paste0("Skip variant ",i,": has no GeneDetail.ensGene annotation"))
      next
    }

    my_gene <- as.character(myvcf_annot_info[i,"Gene.ensGene"])
    if(my_func == "exonic"){
      #my_nm_list <- parse_AAChange.refGene(my_refgene)
      my_enst_list <- parse_AAChange.refGene(my_ensgene)
    }else{
      #my_nm_list <- parse_GeneDetail.refGene(my_refgene)
      my_enst_list <- parse_GeneDetail.refGene(my_ensgene)
    }

    my_snp_pos_geno <- start(ranges(rowRanges(myvcf_annot)))[i]
    my_chr <- as.character(seqnames(rowRanges(myvcf_annot))[i])
    if(length(grep("chr",my_chr))<1){
          my_chr <- paste0("chr",my_chr)
    }
    my_gencode_seque <- MORFEE_DATA[["GENCODE_SEQ"]][which(MORFEE_DATA[["GENCODE_SEQ_ORDER"]]==my_chr)]

    # Loop for each transcript (row)
    for(enst in 1:nrow(my_enst_list)){

      if(!(as.character(my_enst_list[enst,5]) %in% c("A","T","C","G"))){
        message(paste0("Skip variant ",i,": reference allele not A, T, C or G! ",my_enst_list[enst,5]))
#       message(" it could be an indel, but not yet supported!")
        next
      }

      if(!(as.character(my_enst_list[enst,6]) %in% c("A","T","C","G"))){
        message(paste0("Skip variant ",i,": alternative allele not A, T, C or G! ",my_enst_list[enst,6]))
#       message(" it could be an indel, but not yet supported!")
        next
      }

      my_enst <- my_enst_list[enst,1]
      my_snp <- my_enst_list[enst,2]
      my_upordown <- my_enst_list[enst,3]
      my_snp_pos_rel <- as.numeric(my_enst_list[enst,4]) # Position from ATG, ex: 94
      my_nm_id <- grep(strsplit(my_enst,"\\.")[[1]][1], MORFEE_DATA[["GENCODE_METAD"]]$V1)
      my_ref <- as.character(my_enst_list[enst,5])
      my_alt <- as.character(my_enst_list[enst,6])

      my_nm <- MORFEE_DATA[["GENCODE_METAD"]][my_nm_id,2][1]

      if(length(my_nm_id)==0){
        #message(paste0("Skip variant ",i,": annotation for ",my_nm_id," was not found in GENCODE"))
        #next
        my_nm <- "na"
      }

      if(MORFEE_DATA[["GRCh"]]==37){
        my_transcript_id <- grep(paste0(my_enst,"."), MORFEE_DATA[["GENCODE_ANNOT"]]$transcript_id)
      }else if(MORFEE_DATA[["GRCh"]]==38){
        my_transcript_id <- grep(my_enst, MORFEE_DATA[["GENCODE_ANNOT"]]$transcript_id)
      }else{
        stop("Reference Genome unknown")
      }
      gencode_annot_sub <- MORFEE_DATA[["GENCODE_ANNOT"]][my_transcript_id,]
      gencode_annot_sub <- gencode_annot_sub[gencode_annot_sub$seqid==my_chr,]
      gencode_annot_cds <- gencode_annot_sub[gencode_annot_sub$type=="CDS",]
      gencode_annot_exon <- gencode_annot_sub[gencode_annot_sub$type=="exon",]
      gencode_annot_transcript_type <- unique(gencode_annot_exon$transcript_type)

      if(length(gencode_annot_transcript_type)<1){
        message(paste0("Skip variant ",i,": transcript type is not protein coding"))
        next
      }else if(!(gencode_annot_transcript_type %in% "protein_coding")){
        message(paste0("Skip variant ",i,": transcript type is not protein coding"))
        next
      }

      my_init_codon_r <- gencode_annot_sub[gencode_annot_sub$type=="start_codon",]

      if(nrow(my_init_codon_r)==0){
        message(paste0("Skip variant ",i,": associated transcript without start_codon"))
        next
      }

      my_init_codon_start <- as.numeric(my_init_codon_r[,"start"])[1]  
      my_init_codon_end   <- as.numeric(my_init_codon_r[,"end"])[1]  

      my_stop_codon_r <- gencode_annot_sub[gencode_annot_sub$type=="stop_codon",]

      if(nrow(my_stop_codon_r)==0){
        message(paste0("Skip variant ",i,": associated transcript without stop_codon"))
        next
      }

      my_stop_codon_start <- as.numeric(my_stop_codon_r[,"start"])[1]  
      my_stop_codon_end   <- as.numeric(my_stop_codon_r[,"end"])[1]  

      if(my_init_codon_end < my_stop_codon_end){
        my_init_codon_5 <- my_init_codon_start
        my_init_codon_3 <- my_init_codon_end
        my_stop_codon_5 <- my_stop_codon_start
        my_stop_codon_3 <- my_stop_codon_end
      }else{
        my_init_codon_5 <- my_init_codon_end
        my_init_codon_3 <- my_init_codon_start
        my_stop_codon_5 <- my_stop_codon_end
        my_stop_codon_3 <- my_stop_codon_start
      }

      ########################################
      # Test gene orientation
      if(my_init_codon_end < my_stop_codon_end){

        gencode_annot_exon <- gencode_annot_exon[order(as.numeric(gencode_annot_exon$exon_number), decreasing = FALSE),]
        exons_length <- (gencode_annot_exon[,"end"]+1)-gencode_annot_exon[,"start"]
        my_snp_exon <- which(gencode_annot_exon[,"start"] <= my_snp_pos_geno & gencode_annot_exon[,"end"] >= my_snp_pos_geno)
        my_start_exon <- which(gencode_annot_exon[,"start"] <= my_init_codon_5 & gencode_annot_exon[,"end"] >= my_init_codon_5)
        my_stop_exon <- which(gencode_annot_exon[,"start"] <= my_stop_codon_5 & gencode_annot_exon[,"end"] >= my_stop_codon_5)
        CDS_start <- gencode_annot_cds$start[1]
        CDS_end <- gencode_annot_cds$end[nrow(gencode_annot_cds)]
        my_cdna_length_A <- my_cdna_length_B <- my_cdna_length_C <- intron_length <- exon_start <- cdna_length_w_intron <- 0
        exon_start_vec <- exon_end <- c()

        if(length(my_snp_exon)==0){
          message(paste0("Skip variant ",i,": not in an exonic part"))
          next
        }

        # get relative position of exon start and end to determine in which exon any relative position is 
        for(annot_exon in 1:nrow(gencode_annot_exon)){
          exon_start_vec <- append(exon_start_vec, exon_start)
          exon_end <- append(exon_end,exon_start + exons_length[annot_exon])
          cdna_length_w_intron <- cdna_length_w_intron + exons_length[annot_exon]
          exon_start <- cdna_length_w_intron
          }
        exon_df <- data.frame(exon_start_vec, exon_end)
        colnames(exon_df) <- c("start","end")
        exon_df$start <- exon_df$start+1



        #fix annovar nomenclature of variant position 
          if(my_snp_pos_geno < my_init_codon_5){ #if snp in 5utr
            if(my_snp_exon != my_start_exon){
              for(exon in my_snp_exon:(my_start_exon-1)){
                intron_length <- intron_length + (gencode_annot_exon[,"start"][exon + 1] - gencode_annot_exon[,"end"][exon])
              }
            my_snp_pos_rel <- as.integer(my_enst_list[enst,4]) - abs(intron_length) + length(my_snp_exon:(my_start_exon-1))
            }
          }

        # Reference sequence stats
        for(exon_i in 1:nrow(gencode_annot_exon)){
          if(exon_i==1){
            my_cdna   <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
          }else{
            my_cdna_i <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
            my_cdna <- c(my_cdna, my_cdna_i)
          }

          # Calcul position of variant in cDNA
          if(exon_i < my_snp_exon){
            my_cdna_length_A <- my_cdna_length_A + exons_length[exon_i]
          }else if(exon_i == my_snp_exon){
            my_snp_pos_cdna <- my_cdna_length_A + ((my_snp_pos_geno+1) - gencode_annot_exon[exon_i,"start"])
          }

          # Calcul position of reference A(TG) in cDNA
          if(exon_i < my_start_exon){
            my_cdna_length_B <- my_cdna_length_B + exons_length[exon_i]
          }else if(exon_i == my_start_exon){
            my_init_codon_5_cdna <- my_cdna_length_B + ((my_init_codon_5 + 1) - gencode_annot_exon[exon_i,"start"])
          }
            # Calcul position of reference STOP in cDNA
            if(exon_i < my_stop_exon){
                my_cdna_length_C <- my_cdna_length_C + exons_length[exon_i]
            }else if(exon_i == my_stop_exon){
                my_stop_codon_5_cdna <- my_cdna_length_C + ((my_stop_codon_5 + 1) - gencode_annot_exon[exon_i,"start"])
          }
        }

        # Find codon start in reference sequence
          for(j in 1:length(MORFEE_DATA[["SEQ_INIT"]])){

            stats_orig_j <- matchPattern(MORFEE_DATA[["SEQ_INIT"]][[j]], my_cdna)

            if(j==1){
              stats_orig <- stats_orig_j
            }else{
              stats_orig <- c(stats_orig, stats_orig_j)
            }
          }


          # Found all STOP in reference sequence
          for(j in 1:length(MORFEE_DATA[["SEQ_STOP"]])){

            stats_stop_orig_j <- matchPattern(MORFEE_DATA[["SEQ_STOP"]][[j]], my_cdna)

            if(j==1){
              stats_stop_orig <- stats_stop_orig_j
            }else{
              stats_stop_orig <- c(stats_stop_orig, stats_stop_orig_j)
            }
          }


        my_ref_allele <- as.character(subseq(my_cdna, start=my_snp_pos_cdna, end=my_snp_pos_cdna))
        if(my_ref_allele!=my_enst_list[enst,5]){
          message(paste0("Skip variant ",i,": mismatch between alleles"))
          next
        }

        my_ref_ATG <- as.character(subseq(my_cdna, start=my_init_codon_5_cdna, end=(my_init_codon_5_cdna+2)))
        if(my_ref_ATG!="ATG"){
          message(paste0("Skip variant ",i,": reference ATG no detected"))
          next
        }

        my_ref_stop <- as.character(subseq(my_cdna, start=my_stop_codon_5_cdna, end=(my_stop_codon_5_cdna+2)))
        if(my_ref_stop!="TGA" & my_ref_stop!="TAA" & my_ref_stop!="TAG"){
          message(paste0("Skip variant ",i,": reference STOP not detected"))
          next
        }


        # Replace my sequence with SNP
        my_cdna_updated <- replaceLetterAt(my_cdna, my_snp_pos_cdna, my_enst_list[enst,6])

        ########################################
        }else{
        # message("Opposite orientation!")
          #if(MORFEE_DATA[["GRCh"]] == 37){
          #gencode_annot_exon <- gencode_annot_exon[order(as.numeric(gencode_annot_exon$exon_number), decreasing = TRUE),]
          #CDS_start <- gencode_annot_cds$end[nrow(gencode_annot_cds)]
          #CDS_end <- gencode_annot_cds$start[1]
          #}else{
          gencode_annot_exon <- gencode_annot_exon[order(as.numeric(gencode_annot_exon$exon_number), decreasing = FALSE),]
          CDS_start <- gencode_annot_cds$end[1]
          CDS_end <- gencode_annot_cds$start[nrow(gencode_annot_cds)]
          #}
          exons_length <- (gencode_annot_exon[,"end"]+1)-gencode_annot_exon[,"start"]
          my_snp_exon <- which(gencode_annot_exon[,"start"] <= my_snp_pos_geno & gencode_annot_exon[,"end"] >= my_snp_pos_geno)
          my_start_exon <- which(gencode_annot_exon[,"start"] <= my_init_codon_5 & gencode_annot_exon[,"end"] >= my_init_codon_5)
          my_stop_exon <- which(gencode_annot_exon[,"start"] <= my_stop_codon_5 & gencode_annot_exon[,"end"] >= my_stop_codon_5)
          my_cdna_length_A <- my_cdna_length_B <- my_cdna_length_C <- intron_length <- exon_start <- cdna_length_w_intron <- 0
          exon_start_vec <- exon_end <- c()

          if(length(my_snp_exon)==0){
            message(paste0("Skip variant ",i,": not in an exonic part"))
            next
          }

          # get relative position of exon start and end to determine in which exon any relative position is 
          for(annot_exon in 1:nrow(gencode_annot_exon)){
            exon_start_vec <- append(exon_start_vec, exon_start)
            exon_end <- append(exon_end,exon_start + exons_length[annot_exon])
            cdna_length_w_intron <- cdna_length_w_intron + exons_length[annot_exon]
            exon_start <- cdna_length_w_intron
            }
          exon_df <- data.frame(exon_start_vec, exon_end)
          colnames(exon_df) <- c("start","end")
          exon_df$start <- exon_df$start+1



          if(my_snp_pos_geno > my_init_codon_5){ #if snp in 5utr
            if(my_snp_exon != my_start_exon){
              for(exon in my_snp_exon:(my_start_exon-1)){
                intron_length <- intron_length + (gencode_annot_exon[,"start"][exon] - gencode_annot_exon[,"end"][exon + 1])
              }
            my_snp_pos_rel <- as.integer(my_enst_list[enst,4]) - abs(intron_length) + length(my_snp_exon:(my_start_exon-1))
            }
          }

          # Reference sequence stats
          for(exon_i in 1:nrow(gencode_annot_exon)){
            if(exon_i==1){
              my_cdna   <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
              my_cdna   <- reverse(complement(my_cdna))
            }else{
              my_cdna_i <- subseq(my_gencode_seque, start=gencode_annot_exon[exon_i,"start"], end=gencode_annot_exon[exon_i,"end"])[[1]]
              my_cdna_i <- reverse(complement(my_cdna_i))

              my_cdna <- c(my_cdna, my_cdna_i)
            }

            # Calcul position of variant in cDNA
            if(exon_i < my_snp_exon){
              my_cdna_length_A <- my_cdna_length_A + exons_length[exon_i]
            }else if(exon_i == my_snp_exon){
              my_snp_pos_cdna <- my_cdna_length_A + (gencode_annot_exon[exon_i,"end"] - (my_snp_pos_geno-1))
            }

            # Calcul position of reference A(TG) in cDNA
            if(exon_i < my_start_exon){
              my_cdna_length_B <- my_cdna_length_B + exons_length[exon_i]
            }else if(exon_i == my_start_exon){
              my_init_codon_5_cdna <- my_cdna_length_B + (gencode_annot_exon[exon_i,"end"] - (my_init_codon_3+1))
            }

            # Calcul position of reference STOP in cDNA
            if(exon_i < my_stop_exon){
                my_cdna_length_C <- my_cdna_length_C + exons_length[exon_i]
            }else if(exon_i == my_stop_exon){
                my_stop_codon_5_cdna <- my_cdna_length_C + gencode_annot_exon[exon_i,"end"] - ((my_stop_codon_3+1))
            }
          }

          # Find codon start in reference sequence
          for(j in 1:length(MORFEE_DATA[["SEQ_INIT"]])){
            stats_orig_j <- matchPattern(MORFEE_DATA[["SEQ_INIT"]][[j]], my_cdna)
            if(j==1){
              stats_orig <- stats_orig_j
            }else{
              stats_orig <- c(stats_orig, stats_orig_j)
            }
          }

          # Found all STOP in reference sequence
          for(j in 1:length(MORFEE_DATA[["SEQ_STOP"]])){

            stats_stop_orig_j <- matchPattern(MORFEE_DATA[["SEQ_STOP"]][[j]], my_cdna)

            if(j==1){
              stats_stop_orig <- stats_stop_orig_j
            }else{
              stats_stop_orig <- c(stats_stop_orig, stats_stop_orig_j)
            }
          }

          my_ref_allele <- as.character(subseq(my_cdna, start=my_snp_pos_cdna, end=my_snp_pos_cdna))
          if(my_ref_allele!=my_enst_list[enst,5]){
            message(paste0("Skip variant ",i,": mismatch between alleles"))
            next
          }

          my_ref_ATG <- as.character(subseq(my_cdna, start=my_init_codon_5_cdna, end=(my_init_codon_5_cdna+2)))
          if(my_ref_ATG!="ATG"){
            message(paste0("Skip variant ",i,": reference ATG no detected"))
            next
          }

          my_ref_stop <- as.character(subseq(my_cdna, start=my_stop_codon_5_cdna, end=(my_stop_codon_5_cdna+2)))
          if(my_ref_stop!="TGA" & my_ref_stop!="TAA" & my_ref_stop!="TAG"){
            message(paste0("Skip variant ",i,": reference STOP not detected"))
            next
          }

          # Replace my sequence with SNP
          my_cdna_updated <- replaceLetterAt(my_cdna, my_snp_pos_cdna, my_enst_list[enst,6])

        }
        ########################################

        # Find codon start in mutated sequence
        for(j in 1:length(MORFEE_DATA[["SEQ_INIT"]])){
          stats_mut_j <- matchPattern(MORFEE_DATA[["SEQ_INIT"]][[j]], my_cdna_updated) 
          if(j==1){
            stats_mut <- stats_mut_j
          }else{
            stats_mut <- c(stats_mut, stats_mut_j)
          }
        }


        # Found all STOP in mutated sequence
        for(j in 1:length(MORFEE_DATA[["SEQ_STOP"]])){

          stats_stop_mut_j <- matchPattern(MORFEE_DATA[["SEQ_STOP"]][[j]], my_cdna_updated) 

          if(j==1){
            stats_stop_mut <- stats_stop_mut_j
          }else{
            stats_stop_mut <- c(stats_stop_mut, stats_stop_mut_j)
          }
        }

        # Compare start codon in reference and mutated sequences
        new.atg <- ranges(stats_mut)[!c(ranges(stats_mut) %in% ranges(stats_orig)) ,]
        fixed_nomenclature <- paste0("c.",my_enst_list[enst,3],my_snp_pos_rel,my_enst_list[enst,5],">",my_enst_list[enst,6])

        START_sequences <- NA
        #if(length(stats_orig) == length(stats_mut)){ ## to find start codon to start codon substitution ##

        my_cdna_updated_start_seq <- as.character(mapply(subseq,as.character(my_cdna_updated),start(stats_mut),end(stats_mut))) #list of START codon sequence in mutated sequence
        my_cdna_start_seq <- as.character(mapply(subseq,as.character(my_cdna),start(stats_mut),end(stats_mut))) #list of START codon sequence in original sequence
        index <- which(my_cdna_updated_start_seq != my_cdna_start_seq) # index of changed codon
        if(length(index)>0){
          new.atg <- ranges(stats_mut[index])   ## new start codon ##
        } 
        #}


        # Compare stop codons in reference and mutated sequences
        del.stop <- ranges(stats_stop_orig)[!c(ranges(stats_stop_orig) %in% ranges(stats_stop_mut)) ,]
        new.stop <- ranges(stats_stop_mut)[!c(ranges(stats_stop_mut) %in% ranges(stats_stop_orig)) ,]

        if(length(new.atg)>0){
          for(atg in 1:length(new.atg)){
            message("New TIS detected!") ### detection of modification type ###
            start_codon_list <- list(ATG = "canonical_TIS", ACG = "non_canonical_TIS",CTG = "non_canonical_TIS",GTG = "non_canonical_TIS",TTG = "non_canonical_TIS",AAG = "non_canonical_TIS",AGG = "non_canonical_TIS",ATC = "non_canonical_TIS",ATA = "non_canonical_TIS",ATT = "non_canonical_TIS")
            START_sequences <- c(subseq(as.character(my_cdna),start(new.atg)[atg],end(new.atg)[atg]), subseq(as.character(my_cdna_updated),start(new.atg)[atg],end(new.atg)[atg]))
            mut_codon <- start_codon_list[START_sequences[2]]
            if(START_sequences[1] %in% names(start_codon_list)){
              orig_codon <- start_codon_list[subseq(as.character(my_cdna),start(new.atg)[atg],end(new.atg)[atg])]
            }else{
              orig_codon <- "non_TIS"
            }
            modification_info <- paste0(orig_codon,"(",START_sequences[1],")","_to_",mut_codon,"(",START_sequences[2],")")

            new_atg_seq <- START_sequences[2]
            #get kozak sequence and calculate strength
            kozak_seq <- ifelse(start(new.atg)[atg] < 11 || end(new.atg)[atg] > (length(my_cdna_updated)-10), "na", subseq(as.character(my_cdna_updated),start(new.atg)[atg]-10,end(new.atg)[atg]+10))
            kss <- ifelse(start(new.atg)[atg] < 11 || end(new.atg)[atg] > (length(my_cdna_updated)-10), "na", kozak_similarity_score(kozak_seq))

            if(kozak_seq != "na"){
              kozak_3 <- substr(kozak_seq,8,8)
              kozak_4 <- substr(kozak_seq,14,14)
              if(paste0(kozak_3,kozak_4) == "AG" || paste0(kozak_3,kozak_4) == "GG"){
                kozak_strength <- "strong"
              }else if(kozak_3 == "A" ||  kozak_3 == "G" || kozak_4 == "G"){
                kozak_strength <- "moderate"
              }else{
                kozak_strength <- "weak"
              }
            }else{
              kozak_strength <- "na"
            }


            if(my_init_codon_end < my_stop_codon_end){
              my.strand <- "forward"
            }else{
              my.strand <- "reverse"
            }



            new.atg.distance <- my_init_codon_5_cdna - (start(new.atg)[atg])

            test.frame.uatg <- (new.atg.distance%%3)

            if(test.frame.uatg==0){
              in.frame <- "in_frame"
            }else{
              in.frame <- paste0("out_of_frame_(",test.frame.uatg,")")
            }

            # Found all STOP
            for(j in 1:length(MORFEE_DATA[["SEQ_STOP"]])){

              my_stop_j <- matchPattern(MORFEE_DATA[["SEQ_STOP"]][[j]], my_cdna_updated[ start(range(new.atg[atg])) : length(my_cdna_updated) ])

              if(j==1){
                my_stops <- my_stop_j
              }else{
                my_stops <- c(my_stops, my_stop_j)
              }
            }

            # First stop in phase to the new codon start
            my_stops_sort <- sort(start(ranges(my_stops)))
            my_first_stop <- my_stops_sort[my_stops_sort%%3==1][1]
            my_first_stop_cdna <- my_first_stop + start(new.atg[atg]) - 1 #position of first stop in cdna
            ref_tis_to_stop_distance <- my_first_stop - new.atg.distance - 1 #distance from reference ATG to first uSTOP
            orf_size <- my_first_stop + 2
            intron_length <- my_stop_pos_g <- 0

            generated.prot.length <- (my_first_stop-1)/3
            ref.prot.length <- (sum(gencode_annot_cds[,"end"]+1 - gencode_annot_cds[,"start"] ) -3)/3

            if(is.na(my_first_stop)){
              message(paste0("Skip variant ",i,": new TIS detected but no STOP in phase"))
              next
            }

            start_c_pos <- NA
            stop_c_pos <- NA

            if(start(new.atg[atg]) < my_init_codon_5_cdna){ #utr5
              start_c_pos <- paste0("c.",start(new.atg[atg]) - my_init_codon_5_cdna) 
            }else if(start(new.atg[atg]) <= (my_stop_codon_5_cdna+2)){ # exonic
              start_c_pos <- paste0("c.",start(new.atg[atg]) - my_init_codon_5_cdna + 1)
            }else{ #3utr
              start_c_pos <- paste0("c.*",start(new.atg[atg]) - my_stop_codon_5_cdna - 2)
            }

            if(my_first_stop_cdna < my_init_codon_5_cdna){
              stop_c_pos <- paste0("c.",ref_tis_to_stop_distance) 
            }else if(my_first_stop_cdna <= (my_stop_codon_5_cdna+2)){
              stop_c_pos <- paste0("c.",ref_tis_to_stop_distance+1) 
            }else{
              stop_c_pos <- paste0("c.*", my_first_stop_cdna - my_stop_codon_5_cdna - 2)
            }

            #find new stop exonic position - forward
            my_first_stop_exon <- which(exon_df[,"start"] <= my_first_stop_cdna & exon_df[,"end"] >= my_first_stop_cdna)
            
            if(my.strand == "forward"){
              #calculate genomic pos of first nt of uTIS
              uTIS_g_pos <- my_snp_pos_geno + (start(new.atg[atg]) - my_snp_pos_cdna)
              ####calculate distance with introns - forward
              if(my_first_stop_exon != my_snp_exon){
                for(exon in (my_first_stop_exon-1):(my_snp_exon)){
                  intron_length <- intron_length + (gencode_annot_exon[,"start"][exon + 1] - gencode_annot_exon[,"end"][exon]) - 1                  
                  }
                my_stop_pos_g <- uTIS_g_pos + abs(intron_length) + my_first_stop - 1 
                }else{
                my_stop_pos_g <- uTIS_g_pos + my_first_stop - 1
                }
              }else{
                uTIS_g_pos <- my_snp_pos_geno - (start(new.atg[atg]) - my_snp_pos_cdna)
                if(my_first_stop_exon != my_snp_exon){
                  for(exon in (my_first_stop_exon - 1):(my_snp_exon)){
                    intron_length <- intron_length + (((gencode_annot_exon[,"end"][exon + 1] - gencode_annot_exon[,"start"][exon])) + 1)
                    }
                  my_stop_pos_g <- uTIS_g_pos - abs(intron_length) - my_first_stop + 1 
                }else{
                  my_stop_pos_g <- uTIS_g_pos - my_first_stop + 1
                }

              }

            #get orf stop sequence
            STOP_sequence <- subseq(as.character(my_cdna),my_first_stop_cdna,my_first_stop_cdna+2)

            # Determine whether the ORF is overlapping, not overlapping or elongated CD
            if(my_first_stop_cdna <= my_init_codon_5_cdna || start(new.atg)[atg] >= my_stop_codon_5_cdna){
              overlapping.prot <- "not_overlapping"
            }else if( in.frame=="in_frame" & (my_first_stop_cdna > my_init_codon_5_cdna) & my_init_codon_5_cdna > (start(new.atg)[atg]) ){ ## add if atg before reference atg
              overlapping.prot <- "elongated_CDS"
            }else if(my_init_codon_5_cdna <= (start(new.atg)[atg]) & my_stop_codon_5_cdna+2 >= my_first_stop_cdna ){
              overlapping.perc <- (my_first_stop/(ref.prot.length*3))*100
              overlapping.perc.round <- round(overlapping.perc, digits = 2)
              overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
            }else if(my_init_codon_5_cdna <= (start(new.atg)[atg]) & my_stop_codon_5_cdna < my_first_stop_cdna ){
              overlapping.perc <- ((my_stop_codon_5_cdna - (start(new.atg)[atg]))/(ref.prot.length*3))*100
              overlapping.perc.round <- round(overlapping.perc, digits = 2)
              overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
            }else{
              overlapping.perc <- ref_tis_to_stop_distance/(ref.prot.length*3)*100
              overlapping.perc.round <- round(overlapping.perc, digits = 2)
              overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
            }


            if((start(new.atg)[atg]) > my_stop_codon_5_cdna+2){ 
                message( paste("For",my_gene,"-",my_enst,"and",my_snp))
                message(paste0(" - New dTIS detected at: ",-new.atg.distance," from the main ATG!"))
                message( paste(" - new dTIS is",in.frame,"to the main ATG!"))
                message( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
  #             message(paste0(" - DEBUG: i=",i," ; nm=",nm))
  
                message("\n\n")

                # Update myvcf_annot_info
                new_field <- paste( na.omit(c( myvcf_annot_info[i,"MORFEE_dTIS"], 
                                  paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",new.atg.distance,",",in.frame,",",new_atg_seq,",",my_func,",",mut_codon,",",modification_info,",",overlapping.prot,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",orf_size,",",CDS_start,",",uTIS_g_pos,",",my_stop_pos_g,",",CDS_end,",",kozak_seq,",",kss,",",start_c_pos,",",stop_c_pos,",",kozak_strength,",",STOP_sequence,",",my_nm)) )
                                , collapse="|")

                myvcf_annot_info[i,"MORFEE_dTIS"] <- new_field
              }else if((start(new.atg)[atg]) >= my_init_codon_5_cdna & (start(new.atg)[atg]) <= my_stop_codon_5_cdna+2){
                message( paste("For",my_gene,"-",my_enst,"and",my_snp))
                message(paste0(" - New intTIS detected at: ",-new.atg.distance," from the main ATG!"))
                message( paste(" - new intTIS is",in.frame,"to the main ATG!"))
                message( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
  #             message(paste0(" - DEBUG: i=",i," ; nm=",nm))
  
                message("\n\n")

                # Update myvcf_annot_info
                new_field <- paste( na.omit(c( myvcf_annot_info[i,"MORFEE_intTIS"], 
                                  paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",new.atg.distance,",",in.frame,",",new_atg_seq,",",my_func,",",mut_codon,",",modification_info,",",overlapping.prot,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",orf_size,",",CDS_start,",",uTIS_g_pos,",",my_stop_pos_g,",",CDS_end,",",kozak_seq,",",kss,",",start_c_pos,",",stop_c_pos,",",kozak_strength,",",STOP_sequence,",",my_nm)) )
                                , collapse="|")

                myvcf_annot_info[i,"MORFEE_intTIS"] <- new_field
                }else{
                message( paste("For",my_gene,"-",my_enst,"and",my_snp))
                message(paste0(" - New uTIS detected at: ",new.atg.distance," from the main ATG!"))
                message( paste(" - new uTIS is",in.frame,"to the main ATG!"))
                message( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
  #             message(paste0(" - DEBUG: i=",i," ; nm=",nm))
                message("\n\n")

                # Update myvcf_annot_info
                new_field <- paste( na.omit(c( myvcf_annot_info[i,"MORFEE_uTIS"], 
                                  paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",new.atg.distance,",",in.frame,",",new_atg_seq,",",my_func,",",mut_codon,",",modification_info,",",overlapping.prot,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",orf_size,",",CDS_start,",",uTIS_g_pos,",",my_stop_pos_g,",",CDS_end,",",kozak_seq,",",kss,",",start_c_pos,",",stop_c_pos,",",kozak_strength,",",STOP_sequence,",",my_nm)) )
                                , collapse="|")

                myvcf_annot_info[i,"MORFEE_uTIS"] <- new_field
                }
          }
        }# END new ATG
        

        if(length(new.stop)>0){
          for(stop in 1:length(new.stop)){
            message("New STOP detected!")
            new_stop_seq <- c(subseq(as.character(my_cdna_updated),start(new.stop)[stop],end(new.stop)[stop]))

            if(my_init_codon_end < my_stop_codon_end){
              my.strand <- "forward"
            }else{
              my.strand <- "reverse"
            }

            new.stop.distance <- my_init_codon_5_cdna - (start(new.stop)[stop])

            test.frame.stop <- (new.stop.distance%%3)

            if(test.frame.stop==0){
              in.frame <- "in_frame"
            }else{
              in.frame <- paste0("out_of_frame_(",test.frame.stop,")")
            }

            # Found all START
            for(j in 1:length(MORFEE_DATA[["SEQ_INIT"]])){

              my_start_j <- matchPattern(MORFEE_DATA[["SEQ_INIT"]][[j]], my_cdna_updated[ 1 : start(range(new.stop[stop])) ])

              if(j==1){
                my_starts <- my_start_j
              }else{
                my_starts <- c(my_starts, my_start_j)
              }
            }

            # First start in phase to the new codon stop
            my_starts_sort <- sort(start(ranges(my_starts)))
            my_starts_sort_from_stop <- start(new.stop) - my_starts_sort #gives the distance from the stop to the starts instead of the distance from pos 1 to START which allow for the %%3 calculation to work
            #my_first_start <- my_starts_sort_from_stop[my_starts_sort_from_stop%%3==0][length(my_starts_sort_from_stop[my_starts_sort_from_stop%%3==0])]
            my_first_starts <- my_starts_sort_from_stop[my_starts_sort_from_stop%%3==0] #all in-frame starts
            my_first_starts_cdna <- start(new.stop) - my_first_starts #get starts in frame in cdna pos again

            ref_tis_to_stop_distance <- start(new.stop[stop]) - my_init_codon_5_cdna - 1 #distance from reference ATG to first uSTOP

            if(length(my_first_starts) == 0){
              message(paste0("Skip variant ",i,": new STOP detected but no TIS in phase"))
              next
            }

            for(start_in_frame in my_first_starts_cdna){
              
              #starts_between_TIS_and_STOP <- start(stats_stop_mut)[start(stats_stop_mut)>starts_in_frame & start(stats_stop_mut)<start(new.stop)]
              stop_in_frame <- start(stats_stop_mut)[ ((start_in_frame - start(stats_stop_mut)) %%3)==0 ]
              stop_in_frame_upstream <- stop_in_frame[stop_in_frame>start_in_frame & stop_in_frame<start(new.stop)] #only gives stops in between the TIS and the new STOP
              if (length(stop_in_frame_upstream > 0)){
                next
              }

            start_c_pos <- NA
            stop_c_pos <- NA

            if(start_in_frame < my_init_codon_5_cdna){ #utr5
              start_c_pos <- paste0("c.",start_in_frame - my_init_codon_5_cdna) 
            }else if(start_in_frame <= (my_stop_codon_5_cdna+2)){ # exonic
              start_c_pos <- paste0("c.",start_in_frame - my_init_codon_5_cdna + 1)
            }else{ #3utr
              start_c_pos <- paste0("c.*",start_in_frame - my_stop_codon_5_cdna - 2)
            }

            if(start(new.stop) < my_init_codon_5_cdna){
              stop_c_pos <- paste0("c.",start(new.stop) - my_init_codon_5_cdna) 
            }else if(start(new.stop) <= (my_stop_codon_5_cdna+2)){
              stop_c_pos <- paste0("c.",start(new.stop) - my_init_codon_5_cdna +1) 
            }else{
              stop_c_pos <- paste0("c.*", start(new.stop) - my_stop_codon_5_cdna - 2)
            }


              #my_first_start_cdna_pos <- start(new.stop) - start_in_frame
              orf_size <- end(new.stop) - start_in_frame + 1
              intron_length <- my_stop_pos_g <- 0

              generated.prot.length <- (orf_size-3)/3
              ref.prot.length <- (sum(gencode_annot_cds[,"end"]+1 - gencode_annot_cds[,"start"] ) -3)/3

              my_first_start_exon <- which(exon_df[,"start"] <= start_in_frame & exon_df[,"end"] >= start_in_frame)

            if(my.strand == "forward"){
              #calculate genomic pos of first nt of uSTOP
              uSTOP_g_pos <- my_snp_pos_geno + (start(new.stop[stop]) - my_snp_pos_cdna)
              ####calculate distance with introns - forward
              if(my_first_start_exon != my_snp_exon){
                for(exon in (my_snp_exon-1):(my_first_start_exon)){
                  intron_length <- intron_length + (gencode_annot_exon[,"start"][exon + 1] - gencode_annot_exon[,"end"][exon]) - 1
                  }
                my_start_pos_g <- uSTOP_g_pos - abs(intron_length) - orf_size + 3
                }else{
                my_start_pos_g <- uSTOP_g_pos - orf_size + 3
                }
              }else{
                uSTOP_g_pos <- my_snp_pos_geno - (start(new.stop[stop]) - my_snp_pos_cdna)
                if(my_first_start_exon != my_snp_exon){
                  for(exon in (my_snp_exon-1):(my_first_start_exon)){
                    intron_length <- intron_length + (((gencode_annot_exon[,"end"][exon + 1] - gencode_annot_exon[,"start"][exon])) + 1)
                    }
                  my_start_pos_g <- uSTOP_g_pos + abs(intron_length) + orf_size - 3
                }else{
                  my_start_pos_g <- uSTOP_g_pos + orf_size - 3
                }

              }

            if(start(new.stop[stop]) <= my_init_codon_5_cdna || start_in_frame >= my_stop_codon_5_cdna){
              overlapping.prot <- "not_overlapping"
            }else if(my_init_codon_5_cdna <= start_in_frame & my_stop_codon_5_cdna+2 >= start(new.stop[stop]) ){
              overlapping.perc <- (orf_size/(ref.prot.length*3))*100
              overlapping.perc.round <- round(overlapping.perc, digits = 2)
              overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
            }else if(my_init_codon_5_cdna <= start_in_frame & my_stop_codon_5_cdna+2 < start(new.stop[stop]) & start_in_frame <= my_stop_codon_5_cdna){
              overlapping.perc <- ((my_stop_codon_5_cdna - start_in_frame)/(ref.prot.length*3))*100
              overlapping.perc.round <- round(overlapping.perc, digits = 2)
             overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
            }else{
              overlapping.perc <- ref_tis_to_stop_distance/(ref.prot.length*3)*100
              overlapping.perc.round <- round(overlapping.perc, digits = 2)
              overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
            }

            
            if((start(new.stop)[stop]) > my_stop_codon_5_cdna+2){ 
                print(overlapping.prot)
                message( paste("For",my_gene,"-",my_enst,"and",my_snp))
                message(paste0(" - New dSTOP creation detected at: ",-new.stop.distance," from the main ATG!"))
                message( paste(" - new dSTOP is",in.frame,"to the main ATG!"))
                message( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
  #             message(paste0(" - DEBUG: i=",i," ; nm=",nm))

                message("\n\n")

                # Update myvcf_annot_info
                new_field <- paste( na.omit(c( myvcf_annot_info_STOP_creation[i,"MORFEE_dSTOP_creation"], 
                                  paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",new.stop.distance,",",in.frame,",",my_func,",",overlapping.prot,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",my_start_pos_g,",",uSTOP_g_pos,",",orf_size,",",CDS_start,",",CDS_end,",",start_c_pos,",",stop_c_pos,",",new_stop_seq,",",my_nm)) )
                                , collapse="|")

                myvcf_annot_info_STOP_creation[i,"MORFEE_dSTOP_creation"] <- new_field
              }else if((start(new.stop)[stop]) >= my_init_codon_5_cdna & (start(new.stop)[stop]) <= my_stop_codon_5_cdna+2){
                message( paste("For",my_gene,"-",my_enst,"and",my_snp))
                message(paste0(" - New intSTOP creation detected at: ",-new.stop.distance," from the main ATG!"))
                message( paste(" - new intSTOP creation is",in.frame,"to the main ATG!"))
                message( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
  #             message(paste0(" - DEBUG: i=",i," ; nm=",nm))

                message("\n\n")

                # Update myvcf_annot_info
                new_field <- paste( na.omit(c( myvcf_annot_info_STOP_creation[i,"MORFEE_intSTOP_creation"],  
                                  paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",new.stop.distance,",",in.frame,",",my_func,",",overlapping.prot,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",my_start_pos_g,",",uSTOP_g_pos,",",orf_size,",",CDS_start,",",CDS_end,",",start_c_pos,",",stop_c_pos,",",new_stop_seq,",",my_nm)) )
                                , collapse="|")

                myvcf_annot_info_STOP_creation[i,"MORFEE_intSTOP_creation"] <- new_field
                }else{
                message( paste("For",my_gene,"-",my_enst,"and",my_snp))
                message(paste0(" - New uSTOP creation detected at: ",new.stop.distance," from the main ATG!"))
                message( paste(" - new uSTOP creation is",in.frame,"to the main ATG!"))
                message( paste(" - new generated protein has a length of",generated.prot.length,"(aa) vs",ref.prot.length,"(aa)"))
  #             message(paste0(" - DEBUG: i=",i," ; nm=",nm))
                message("\n\n")

                # Update myvcf_annot_info
                new_field <- paste( na.omit(c( myvcf_annot_info_STOP_creation[i,"MORFEE_uSTOP_creation"], 
                                  paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",new.stop.distance,",",in.frame,",",my_func,",",overlapping.prot,",",generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",my_start_pos_g,",",uSTOP_g_pos,",",orf_size,",",CDS_start,",",CDS_end,",",start_c_pos,",",stop_c_pos,",",new_stop_seq,",",my_nm)) )
                                , collapse="|")

                myvcf_annot_info_STOP_creation[i,"MORFEE_uSTOP_creation"] <- new_field
                } 

            }
          }
        }# END new STOP

        if(length(del.stop)>0){

          # Use stats_orig, but could use stats_mut
          uatg <- start(stats_orig)[ c(start(del.stop) - start(stats_orig)) > 0]

          uatg_in_frame <- uatg[((start(del.stop) - uatg) %% 3)==0]

          del.stop.distance <- my_init_codon_5_cdna - (start(del.stop)[1])

          if(length(uatg_in_frame)>0){

            if(my_init_codon_end < my_stop_codon_end){
              my.strand <- "forward"
            }else{
              my.strand <- "reverse"
            }

              print(       "STOP deletion detected!")
              print(       " -  uSTOP deletion in ORF detected!")
              print( paste("For",my_gene,"-",my_enst,"and",my_snp))
              print(paste0(" - Deletion of a uSTOP codon detected at: ",-del.stop.distance," from the main ATG!"))
              print( paste(" --- "   , as.character(        my_cdna[start(del.stop)[1]:end(del.stop)[1]] ),
                           " becomes ",as.character(my_cdna_updated[start(del.stop)[1]:end(del.stop)[1]] ) ))
              print( paste(" --- Gene direction:",my.strand))
#             print(paste0(" - DEBUG: i=",i," ; nm=",nm))

            # several uATG could be present, so the protein length will be different
            for(uatg_i in uatg_in_frame){
            # uatg_i = uatg_in_frame[1]

              # Find next stop in frame with uatg_i
              uatg_i_in_frame <- start(stats_stop_mut)[ ((uatg_i - start(stats_stop_mut)) %%3)==0 ]
              id_ustop <- uatg_i_in_frame > uatg_i

              if(sum(id_ustop)>0){
                first_new_stop <- min( uatg_i_in_frame[id_ustop] )
                if(my_func == "UTR3" || my_func == "exonic" || my_func == "UTR5"){
                    if(first_new_stop < end(del.stop)){
                        message(paste0("Skip variant ",i,": STOP codon in frame found upstream of the deletion"))
                        next
                    }
                }

              }else{
                message(paste0("Skip variant ",i,": uTIS detected but no STOP in phase"))
                next
              }
            
            my_first_uatg_i_exon <- which(exon_df[,"start"] <= uatg_i & exon_df[,"end"] >= uatg_i)
            my_first_new_stop_exon <- which(exon_df[,"start"] <= first_new_stop & exon_df[,"end"] >= first_new_stop)
            del_stop_atg_distance <-  start(del.stop) - uatg_i
            del_stop_new_stop_distance <- first_new_stop - start(del.stop) 
            new_orf_size <- del_stop_atg_distance + del_stop_new_stop_distance + 3
            intron_length_st <- intron_length_sp <- my_new_start_pos_g <- 0


            start_c_pos <- NA
            stop_c_pos <- NA

            if(uatg_i < my_init_codon_5_cdna){ #utr5
              start_c_pos <- paste0("c.",uatg_i - my_init_codon_5_cdna) 
            }else if(uatg_i <= (my_stop_codon_5_cdna+2)){ # exonic
              start_c_pos <- paste0("c.",uatg_i - my_init_codon_5_cdna + 1)
            }else{ #3utr
              start_c_pos <- paste0("c.*",uatg_i - my_stop_codon_5_cdna - 2)
            }

            if(first_new_stop < my_init_codon_5_cdna){
              stop_c_pos <- paste0("c.",first_new_stop - my_init_codon_5_cdna) 
            }else if(first_new_stop <= (my_stop_codon_5_cdna+2)){
              stop_c_pos <- paste0("c.",first_new_stop - my_init_codon_5_cdna +1) 
            }else{
              stop_c_pos <- paste0("c.*", first_new_stop - my_stop_codon_5_cdna - 2)
            }

            ### find genomic pos of uatg
            if(my.strand == "forward"){
              #calculate genomic pos of first nt of deleted STOP
              uSTOP_g_pos <- my_snp_pos_geno + (start(del.stop) - my_snp_pos_cdna)
              ####calculate distance with introns - forward
              if(my_first_uatg_i_exon != my_snp_exon){
                for(exon in (my_snp_exon-1):(my_first_uatg_i_exon)){
                  intron_length_st <- intron_length_st + (gencode_annot_exon[,"start"][exon + 1] - gencode_annot_exon[,"end"][exon]) - 1
                  
                  }
                my_new_start_pos_g <- uSTOP_g_pos - abs(intron_length_st) - del_stop_atg_distance 
                }else{
                my_new_start_pos_g <- uSTOP_g_pos - del_stop_atg_distance 
                }
              }else{
                uSTOP_g_pos <- my_snp_pos_geno - (start(del.stop) - my_snp_pos_cdna)
                if(my_first_uatg_i_exon != my_snp_exon){
                  for(exon in (my_snp_exon-1):(my_first_uatg_i_exon)){
                    intron_length_st <- intron_length_st + (((gencode_annot_exon[,"end"][exon + 1] - gencode_annot_exon[,"start"][exon])) + 1)
                    }
                  my_new_start_pos_g <- uSTOP_g_pos + abs(intron_length_st) + del_stop_atg_distance 
                }else{
                  my_new_start_pos_g <- uSTOP_g_pos + del_stop_atg_distance 
                }

              }

            ### find genomic pos of new_stop
            if(my.strand == "forward"){
              ####calculate distance with introns - forward
              if(my_first_new_stop_exon != my_snp_exon){
                for(exon in (my_first_new_stop_exon-1):(my_snp_exon)){
                  intron_length_sp <- intron_length_sp + (gencode_annot_exon[,"start"][exon + 1] - gencode_annot_exon[,"end"][exon]) - 1
                  
                  }
                my_new_stop_pos_g <- uSTOP_g_pos + abs(intron_length_sp) + del_stop_new_stop_distance 
                }else{
                my_new_stop_pos_g <- uSTOP_g_pos + del_stop_new_stop_distance 
                }
              }else{
                if(my_first_new_stop_exon != my_snp_exon){
                  for(exon in (my_first_new_stop_exon-1):(my_snp_exon)){
                    intron_length_sp <- intron_length_sp + (((gencode_annot_exon[,"end"][exon + 1] - gencode_annot_exon[,"start"][exon])) + 1)
                    }
                  my_new_stop_pos_g <- uSTOP_g_pos - abs(intron_length_sp) - del_stop_new_stop_distance 
                }else{
                  my_new_stop_pos_g <- uSTOP_g_pos - del_stop_new_stop_distance 
                }
              }

              # Compute distance and length
              stop.generated.prot.length <- (first_new_stop-uatg_i)/3
              ref.prot.length <- (sum(gencode_annot_cds[,"end"]+1 - gencode_annot_cds[,"start"] ) -3)/3

              uatg_used <- -(my_init_codon_5_cdna - uatg_i) #distance to main ATG
              stop_used <- (my_init_codon_5_cdna - first_new_stop) #distance from main ATG to uSTOP, used to check for overlap and get overlap size
              #if(my_func == "UTR3" || my_func == "exonic"){
              if(uatg_i >= my_init_codon_5_cdna){ 
                stop_used <- uatg_i - my_stop_codon_5_cdna #used to check for overlap if orf is exonic or downstream
              }


              if(stop_used<0){

                # 1. uATG in frame to ref ATG
                uatg_inframe_refatg <- ((uatg_i-my_init_codon_5_cdna)%%3)==0
                # 2. STOP position used > ATG
                stop_used_downstream <- (first_new_stop > my_init_codon_5_cdna)

                if(uatg_inframe_refatg & stop_used_downstream & uatg_i < my_init_codon_5_cdna){  
                  overlapping.prot <- "elongated_CDS"
                }else if(my_init_codon_5_cdna <= uatg_i & my_stop_codon_5_cdna+2 >= first_new_stop){ #for exonic ORF
                  overlapping.perc <- (new_orf_size/(ref.prot.length*3))*100
                  overlapping.perc.round <- round(overlapping.perc, digits = 2)
                  overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
                }else if(my_init_codon_5_cdna <= uatg_i & my_stop_codon_5_cdna+2 < first_new_stop){ #for 3UTR overlapping ORF
                  overlapping.perc <- (my_stop_codon_5_cdna - uatg_i)/(ref.prot.length*3)*100
                  overlapping.perc.round <- round(overlapping.perc, digits = 2)
                  overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
                }else{ #for 5UTR ORF
                  overlapping.perc <- (-stop_used/(ref.prot.length*3))*100
                  overlapping.perc.round <- round(overlapping.perc, digits = 2)
                  overlapping.prot <- paste0("overlapping_",overlapping.perc.round,"%")
                }

              }else{
                overlapping.prot <- "not_overlapping"
              }

              stop.codon <- as.character(stats_stop_mut[start(stats_stop_mut)==first_new_stop])

              test.frame.ustop <- ((my_init_codon_5_cdna-first_new_stop)%%3)
              if(test.frame.ustop==0){
                in.frame <- "in_frame"
              }else{
                in.frame <- paste0("out_of_frame_(",test.frame.ustop,")")
              }

              print(       " --")
              print( paste(" --- using uTIS at",uatg_used,"to the main ATG!"))
              print(paste0(" --- using STOP (",stop.codon,") at ",-stop_used," to the main ATG!"))
              print( paste(" --- new predicted ORF has a length of",stop.generated.prot.length,"(aa) vs",ref.prot.length,"(aa) for the main protein"))
              print( paste(" --- new predicted ORF is",overlapping.prot,"with the main protein"))

              new_STOP_sequence <- c(subseq(as.character(my_cdna_updated),first_new_stop,first_new_stop+2))


              # Update myvcf_annot_info
              #if(my_func == "UTR5"){
              if(start(del.stop) < my_init_codon_5_cdna){
                  new_field <- paste( na.omit(c( myvcf_annot_info_STOP[i,"MORFEE_uSTOP"],
                                     paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",-del.stop.distance,",",in.frame,",",my_func,",",overlapping.prot,",",stop.generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",my_new_start_pos_g,",",my_new_stop_pos_g,",",new_orf_size,",",CDS_start,",",CDS_end,",",start_c_pos,",",stop_c_pos,",",new_STOP_sequence,",",my_nm)) )
                                 , collapse="|")

              myvcf_annot_info_STOP[i,"MORFEE_uSTOP"] <- new_field
              }else if(start(del.stop)  > my_stop_codon_5_cdna+2){
                  new_field <- paste( na.omit(c( myvcf_annot_info_STOP[i,"MORFEE_dSTOP"],
                                     paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",-del.stop.distance,",",in.frame,",",my_func,",",overlapping.prot,",",stop.generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",my_new_start_pos_g,",",my_new_stop_pos_g,",",new_orf_size,",",CDS_start,",",CDS_end,",",start_c_pos,",",stop_c_pos,",",new_STOP_sequence,",",my_nm)) )
                                 , collapse="|")

              myvcf_annot_info_STOP[i,"MORFEE_dSTOP"] <- new_field
              }else{
                  new_field <- paste( na.omit(c( myvcf_annot_info_STOP[i,"MORFEE_intSTOP"],
                                     paste0(my_enst,":",my.strand,",",fixed_nomenclature,",",-del.stop.distance,",",in.frame,",",my_func,",",overlapping.prot,",",stop.generated.prot.length,"[/",ref.prot.length,"]","(aa)",",",my_new_start_pos_g,",",my_new_stop_pos_g,",",new_orf_size,",",CDS_start,",",CDS_end,",",start_c_pos,",",stop_c_pos,",",new_STOP_sequence,",",my_nm)) )
                                 , collapse="|")

              myvcf_annot_info_STOP[i,"MORFEE_intSTOP"] <- new_field
              }
            }
            cat("\n\n")

          }else{

              message(       "STOP deletion detected!")
              message(       " -  uSTOP deletion detected BUT without an upstream TIS (not in an ORF region)!")
              message( paste("For",my_gene,"-",my_enst,"and",my_snp))
              message(paste0(" - Deletion of a uSTOP codon detected at: ",-del.stop.distance," from the main ATG!"))
              message( paste(" --- ",    as.character(        my_cdna[start(del.stop)[1]:end(del.stop)[1]] ),
                             " becomes ",as.character(my_cdna_updated[start(del.stop)[1]:end(del.stop)[1]] ) ))
              message("\n\n")

          }
        }# END del STOP
      }
  }
  myvcf_annot_2 <- myvcf_annot
  myvcf_annot_3 <- myvcf_annot
  info(myvcf_annot) <- myvcf_annot_info
  info(myvcf_annot_2) <- myvcf_annot_info_STOP
  info(myvcf_annot_3) <- myvcf_annot_info_STOP_creation
  return(list(myvcf_annot,myvcf_annot_2,myvcf_annot_3))
}