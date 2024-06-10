#### Functions for filtering

tumbu <- function(
  x,
  covered_region
) {
  #' Tumor Mutational Burden
  #'
  #' @description Calculates the tumor mutational burden
  #'
  #' @param x dataframe. Table of Mutations
  #' @param sureselect string. Covered region by sequencer in Mb
  #'
  #' @return returns tmb, double. tumor mutational burden
  #'
  #' @note Please make sure that your sequencer is implemented.
  #'
  #' @details The tumor mutational burden is calculated by taking into account
  #' @details every mutation including intronic, up-/downstream and so on.
  #' @details The total number of mutations is divided by the area in MegaBases
  #' @details that the sequencer covers
  covered_region <- as.numeric(covered_region)
  tmb <- nrow(x) / covered_region
  return(list(tmb = tmb, exon_region = covered_region))
}

tmb_ex <- function(x, covered_exons, mode = "T", cov_t) {
  require(GenomicRanges)
  if (mode != "T") {
    tmb <- NULL
  } else {
    bed <- read.delim(covered_exons, header = FALSE, comment.char = "#")
    mani_gr <- GRanges(seqnames = bed$V1, strand = "*",
                       ranges = IRanges(start = bed$V2, end = bed$V3))
    mani_gr <- reduce(mani_gr)
    mut_gr <- GRanges(seqnames = x$Chr, strand = "*",
                      ranges = IRanges(start = x$Start, end = x$Start))
    #
    exon_region <- sum(width(mani_gr))/1000000
    used_exon_region <- exon_region * cov_t
    number_used_mutations_tmb <- length(findOverlaps(mut_gr, mani_gr))
    tmb <- number_used_mutations_tmb/used_exon_region
  }
  return(list(tmb = tmb, exon_region = exon_region, used_exon_region = used_exon_region, number_used_mutations_tmb = number_used_mutations_tmb))
}

covered_region <- function(sureselect, mode = "T") {
  if (mode == "T") {
    bed <- read.delim(sureselect, header = FALSE, comment.char = "#")
    gr <- GRanges(
      seqnames = bed$V1, strand = "*",
      ranges = IRanges(start = bed$V2, end = bed$V3)
    )
    gr <- reduce(gr)
    region <- sum(width(gr)) / 1000000
  } else {
    region <- NULL
  }
  return(region)
}


filt <- function(x, func) {
  #' Filter for function
  #'
  #' @description Filters for functionality
  #'
  #' @param x dataframe. Table of Mutations
  #' @param func string. Describes the functionality
  #'
  #' @return returns x, dataframe. Table of Mutations
  #'
  #' @details The mutations can be filtered by the functionality of the part of
  #' @details gene where it is located. Functionalities like intronic or
  #' @details intergenic can be filtered out.
  idx <- grep(func, x$Variant_Classification)
  if (length(idx) > 0) {
    x <- x[-idx, ]
  }
  return(x)
}

vrz <- function(x, mode, protocol = "Tumor_Normal", manifest){
  #' Extract VAF, Readcounts and Zygosity
  #'
  #' @description Extracts VAF, Readcounts (and Zygosity)
  #'
  #' @param x dataframe. Table of Mutations
  #' @param mode string. "N" for SNPs and Indels, "LOH" for LoH
  #'
  #' @return returns x, dataframe. Table of Mutations with extra columns
  #'
  #' @details Mode: "N"
  #' @details We extract VAF, Readcounts and Zygosity out of the column
  #' @details "Otherinfo".
  #' @details Mode: "LOH"
  #' @details We extract VAF and Readcounts both for Tumor and Normal
  #' @details out of the column "Otherinfo". The dataframe x should contain
  #' @details only mutations calls with a lack of heterozygosity.
  if (nrow(x) == 0) {
    return(x)
  }
  if (mode == "T") {
    x$Variant_Reads <- paste(x$t_alt_count, x$t_depth, sep = "|")
    x$Zygosity <- "hom"
    x$Zygosity[x$t_GT == "0/1"] <- "het"
    x$Type <- "Somatic"
  }
  if (mode == "N") {
    x$Variant_Reads <- paste(x$n_alt_count, x$n_depth, sep = "|")
    x$Zygosity <- "hom"
    x$Zygosity[x$n_GT == "0/1"] <- "het"
    x$Type <- "Germline"
  }
  if (mode == "LOH") {
    x$Count_Normal <- paste(x$n_alt_count, x$n_depth, sep = "|")
    x$Count_Tumor <- paste(x$t_alt_count, x$t_depth, sep = "|")
    x$Type <- "LoH"
  }
  
  x$n_AF <- format(round(100 * as.numeric(x$n_alt_count) / as.numeric(x$n_depth), 2), nsmall = 2)
  x$t_AF <- format(round(100 * as.numeric(x$t_alt_count) / as.numeric(x$t_depth), 2), nsmall = 2)
  
  return(x)
}


mrc <- function(x, min_var_count){
  vrc <- strsplit(x = as.character(x$Variant_Reads), split = "|", fixed = TRUE)
  vrc <- unlist(lapply(vrc, function(x){return(x[[1]])}))
  id_np <- which(as.numeric(vrc) < min_var_count)
  if (length(id_np) > 0){
    x <- x[-id_np, ]
  }
  return(x)
}

actionable <- function(x, column, actionable_genes) {
  if(is.na(actionable_genes)) {
    return(x)
  }
  genes <- read.table(actionable_genes, header = F)
  actionable_matches <- which(get(column, x) %in% genes$V1)
  return(x[actionable_matches,])
}

updateGeneNames <- function(x) {
  if(nrow(x) == 0) {
    return(x)
  }

  newSymbols <- as.character(mget(as.character(mget(c(x$Hugo_Symbol), envir = org.Hs.egALIAS2EG, ifnotfound=NA)), envir = org.Hs.egSYMBOL, ifnotfound=NA))
  x$Hugo_Symbol[newSymbols != "NA"] <- newSymbols[newSymbols != "NA"]
  return(x)
}

### Database Queries
isflag <- function(x, dbfile){
  #' Flags
  #'
  #' @description Find flag genes
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #'
  #' @return returns x, dataframe. Table of Mutations with an extra column
  #'
  #' @details Check Mutations in Table, whether or not they are FLAG genes
  #' @details that means, these genes are FrequentLy mutAted GeneS in public
  #' @details exomes. For further information see
  #' @details Shyr C, Tarailo-Graovac M, Gottlieb M, Lee JJ, van Karnebeek C,
  #' @details Wasserman WW. FLAGS, frequently mutated genes in public exomes.
  #' @details BMC Medical Genomics. 2014;7:64. doi:10.1186/s12920-014-0064-y.
  db <- read.delim(dbfile, header = T, sep = "\t", colClasses = "character")
  myhit <- db[nchar(db[, "FLAGS"]) > 0, "FLAGS"]
  x["is_flag"] <- 0
  idx <- which(x$Hugo_Symbol %in% myhit)
  x$is_flag[idx] <- 1
  return(x)
}

isogtsg <- function(x, dbfile){
  #' Cancer Genes
  #'
  #' @description Find cancer genes
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #'
  #' @return returns x, dataframe. Table of Mutations with two extra columns
  #'
  #' @details Check Mutations in Table, whether or not they are cancergenes
  #' @details that means, if these genes are tumorsuppressor or oncogenes.
  #' @details Two columns are added that contain values 0 or 1 describing
  #' @details whether a gene is a tumorsuppressorgene (1) or not (0).
  #' @details The same procedure is used for oncogenes.
  db <- read.delim(dbfile, header = T, sep = "\t", colClasses = "character")
  x$is_tumorsuppressor <- 0
  ts <- which(db$Is.Tumor.Suppressor.Gene == "Yes")
  idx <- which (x$Hugo_Symbol %in% db$Hugo.Symbol[ts])
  x$is_tumorsuppressor[idx] <- 1

  x$is_oncogene <- 0
  og <- which(db$Is.Oncogene == "Yes")
  idx <- which (x$Hugo_Symbol %in% db$Hugo.Symbol[og])
  x$is_oncogene[idx] <- 1
  return(x)
}

ishs <- function(x, dbfile){
  #' Hotspot Detection for SNVs
  #'
  #' @description Finds Hotspots in the Mutationlist
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #' 
  #' @return returns x, dataframe. Table of Mutations with one extra column
  #'
  #' @details Check Mutations in Table, whether or not they are hotspot
  #' @details mutations that means, if these variants are mutated more
  #' @details frequently. Therefore the exact genomic position is checked.
  #' @details In case as hotspot is detected, the aminoacid change is added to
  #' @details the extra column.
  #' 
  #' @note Please make sure that the columns of the input have the right names:
  #' @note required in db: Hugo_Symbol, Genomic_Position,
  #' @note -Reference_Amino_Acid, Variant_Amino_Acid, Amino_Acid_Position
  #' @note required in x: Hugo_Symbol, Start, HGVSp_Short
  # lisths <- read.delim(dbfile, header = T,
  #                      sep = "\t", colClasses = "character")
  lisths <- readxl::read_xls(path = dbfile, sheet = 1)
  x$is_hotspot <- 0
  idh <- which (x$Hugo_Symbol %in% lisths$Hugo_Symbol)
  phs <- which (lisths$Hugo_Symbol %in% x$Hugo_Symbol)

  if (length(idh) > 0){
    hotspot <- rep(FALSE, times = length(idh))
    # Check for matching Gene Locations and then for matching HGVSp_Short
    for (i in 1:length(idh)){
      for (j in 1:length(phs)){
        gp <- as.character(x$Start[idh][i])
        if (grepl(gp, as.character(lisths$Genomic_Position[phs][j]))){
          if (!is.na(x$HGVSp_Short[idh][i])){
            aac <- as.character(x$HGVSp_Short[idh][i])
            l <- nchar(aac)
            aac <- substr(aac, l, l)
            if (grepl(aac, lisths$Variant_Amino_Acid[phs][j])){
              ea <- substr(lisths$Reference_Amino_Acid[phs][j], 1, 1)
              ap <- lisths$Amino_Acid_Position[phs][j]
              aa <- paste0(ea, ap, aac)
              hotspot[i] <- aa
              if (hotspot[i] != FALSE){
                break
              }
            }
          }
        }
      }
    }
    idh_ret <- which(hotspot != FALSE)
    if (length(idh_ret) > 0){
      idh <- idh[idh_ret]
      aachange <- hotspot[idh_ret]
      if (length(idh != 0)){
        x$is_hotspot[idh] <- aachange
      }
    }
  }
  return(x)
}
isihs <- function(x, dbfile){
  #' Hotspot Detection for InDels
  #'
  #' @description Finds Hotspots in the Mutationlist
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #' 
  #' @return returns x, dataframe. Table of Mutations with one extra column
  #'
  #' @details Check Mutations in Table, whether or not they are hotspot
  #' @details mutations that means, if these variants are mutated more
  #' @details frequently. Therefore the exact genomic position is checked.
  #' @details In case as hotspot is detected, the aminoacid change is added to
  #' @details the extra column.
  #' 
  #' @note Please make sure that the columns of the input have the right names:
  #' @note required in db: Hugo_Symbol, Genomic_Position,
  #' @note -Reference_Amino_Acid, Variant_Amino_Acid, Amino_Acid_Position
  #' @note required in x: Hugo_Symbol, Start, HGVSp_Short
  # lisths <- read.delim(dbfile, header = T,
  #                      sep = "\t", colClasses = "character")
  lisths <- readxl::read_xls(path = dbfile, sheet = 2)

  # list should already habe a hotspot column for snps
  fs <- which(x$Variant_Classification != "Frame_Shift_Del")
  fs2 <- which(x$Variant_Classification[fs] != "fFrame_Shift_Ins")
  fs <- fs[fs2]

  # Check Genes
  idh <- which(x$Hugo_Symbol %in% lisths$Hugo_Symbol)
  phs <- which(lisths$Hugo_Symbol %in% x$Hugo_Symbol)
  idh <- intersect(idh, fs)

  if (length(idh) > 0){
    hotspot <- rep(FALSE, times = length(idh))
    # Check for matching Gene Locations and then for matching HGVSp_Short
    for (i in 1: length(idh)){
      for (j in 1:length(phs)){

        l <- nchar(as.character(lisths$Amino_Acid_Position[phs][j]))
        vaa <- as.character(lisths$Variant_Amino_Acid[phs][j])
        vaa <- unlist(strsplit(vaa, ":"))[1]

        if (is.na(vaa)){
          break
          } else{
          if (nchar(vaa) > l + 4){
            l <- round(l / 2) - 1
            start <- substr(vaa, 2, l + 1)
            end <- substr(vaa, l + 4, 2 * l + 3)
            type <- substr(vaa, 2 * (l + 2), nchar(vaa))
            pos <- paste(start, end, sep = "_")

            aac <- as.character(x$HGVSp_Short[idh[i]])
            a <- unlist(strsplit(aac, ","))
            id <- 5 * c(1:length(a))
            b <- unlist(strsplit(a, ":"))
            ast <- paste(b[id], collapse = ",")

            if (grepl(pos, as.character(ast))
                && (x$Hugo_Symbol[idh][i] == lisths$Hugo_Symbol[phs][j])
                && (grepl(type, as.character(x$HGVSp_Short[idh][i])))
                ){
              hotspot[i] <- paste0("p.", vaa)
              if (hotspot[i] != FALSE){
                break
                }
            }
          } else{
            pos <- substr(vaa, 2, l + 1)
            type <- substr(vaa, l + 2, nchar(vaa))
            aac <- as.character(x$HGVSp_Short[idh][i])

            a <- unlist(strsplit(aac, ","))
            id <- 5 * c(1:length(a))
            b <- unlist(strsplit(a, ":"))
            ast <- paste(b[id], sep = ",", collapse = "")

            if (grepl(pos, as.character(ast))
               && (x$Hugo_Symbol[idh][i] == lisths$Hugo_Symbol[phs][j])
               && grepl(type, as.character(ast))){
              hotspot[i] <- paste0("p.", vaa)
              if (hotspot[i] != FALSE){
                break
              }
            }
          }
        }
      }
    }
    idh_ret <- which(hotspot != FALSE)
    if (length(idh_ret) > 0){
      idh <- idh[idh_ret]
      aachange <- hotspot[idh_ret]
      if (length(idh != 0)){
        x$is_hotspot[idh] <- aachange
      }
    }
  }
  return(x)
}

rvis <- function(x, dbfile){
  #' RVIS Score
  #'
  #' @description Extractn RVIS Score
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #' 
  #' @return returns x, dataframe. Table of Mutations with one extra column
  #'
  #' @details Extract RVIS score out of table for each mutation.
  db <- read.delim(dbfile, header = T, sep = "\t", colClasses = "character",
                   row.names = 1)
  x$rvis <- 0.0
  x$rvis <- as.numeric(db[x$Hugo_Symbol, 1])
  return(x)
}

trgt <- function(x, dbfile){
  #' Target
  #'
  #' @description Check for targeted therapy
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #' 
  #' @return returns x, dataframe. Table of Mutations with one extra column
  #'
  #' @details Check if there is a targeted therapy for the mutation's gene.
  #' @details If there is one possible drugs are written in the extra column.
  db <- read.delim(dbfile, header = T, sep = "\t", colClasses = "character",
                   row.names = NULL)
  x$target <- "."
  for (i in 1:nrow(x)){
    idx <- which(db$Gene == x[i, "Hugo_Symbol"])
    if (length(idx)) { 
    x$target[i] <- db$Examples.of.Therapeutic.Agents[idx]
    }
  }
  return(x)
}

dgidb <- function(x, dbfile){
  #' DGIDB
  #'
  #' @description Check for drug gene interactions
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #' 
  #' @return returns x, dataframe. Table of Mutations with one extra column
  #'
  #' @details Check if there is a interaction between drugs and the mutation's
  #' @detials gene. If there is one, information is written in the extra
  #' @details column.
  db <- read.delim(dbfile, header = T, sep = "\t", colClasses = "character",
                   row.names = NULL)
  x$DGIdb <- NA
  for (i in 1:nrow(x)) {
    idx <- which(db$gene_name  == x[i, "Hugo_Symbol"])
    if (length(idx)){
      x$DGIdb[i] <- paste(db$drug_claim_primary_name[idx], collapse = "; ")
    }
  }
  return(x)
}

oncokb <- function(x, dbfile){
  #' OncoKB - drug notation
  #'
  #' @description Looks for actionable Variants in OncoKB's Database
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #'
  #' @return returns x, dataframe. Table of Mutations with one extra column
  #'
  #' @details Check if there is a interaction between drugs and the mutation's
  #' @detials gene. If there is one, information is written in the extra
  #' @details column.
  db <- read.delim(dbfile, header = T, colClasses = "character")
  x$oncokb <- NA
  idh <- which(x$Hugo_Symbol %in% db$Gene)
  pta <- which(db$Gene %in% x$Hugo_Symbol[idh])
  onkb <- rep(FALSE, times = length(idh))
  if (length(pta)){
    for (j in 1:length(pta)){
      id <- which(x$Hugo_Symbol[idh] == db$Gene[pta[j]])
      if (db$Alteration[pta[j]] == "Oncogenic Mutations"){
        for (i in 1:length(id)){
          out <- ""
          out <- paste0(db$Tumor.Type[pta[j]], ",EL:", db$Level[pta[j]],
                        ":", db$Drugs[pta[j]])
          if (onkb[id[i]] == FALSE){
            onkb[id[i]] <- out
          } else {
            onkb[id[i]] <- paste(onkb[id[i]], ";", out)
          }
        }
      } else if (grepl("Exon", as.character(db$Alteration[j]))
                && grepl("mutations", as.character(db$Alteration[j]))){
        ZAHL <- substr(as.character(db$Alteration[j]), 6, 7)
        ex <- paste0("exon", ZAHL)

        for (i in 1:length(id)){
          if (grepl(ex, x$HGVSp_Short[idh[i]])
              && x$Hugo_Symbol[idh] == db$Gene[pta[j]]){
            out <- ""
            out <- paste0(db$Tumor.Type[pta[j]], ",EL:", db$Level[pta[j]],
                          ":", db$Drugs[pta[j]])
            if (onkb[id[i]] == FALSE){
              onkb[id[i]] <- out
            } else {
              onkb[id[i]] <- paste(onkb[id[i]], ";", out)
            }
          }
        }
      } else if (grepl("Exon", as.character(db$Alteration[pta[j]]))
                && grepl("del", as.character(db$Alteration[pta[j]]))){
        ZAHL <- substr(as.character(db$Alteration[pta[j]]), 6, 7)
        ex <- paste0("exon", ZAHL)
        for (i in 1:length(id)){
          if (grepl(ex, x$HGVSp_Short[i])
              && grepl("del", x$HGVSp_Short[i])
              && x$Hugo_Symbol[idh] == db$Gene[pta[j]]){
            out <- ""
            out <- paste0(db$Tumor.Type[pta[j]], ",EL:", db$Level[pta[j]],
                          ":", db$Drugs[pta[j]])
            if (onkb[id[i]] == FALSE){
              onkb[id[i]] <- out
            } else {
              onkb[id[i]] <- paste(onkb[id[i]], ";", out)
            }
          }
        }
      } else if (grepl("Exon", as.character(db$Alteration[pta[j]]))
                && grepl("ins", as.character(db$Alteration[pta[j]]))){
        ZAHL <- substr(as.character(db$Alteration[pta[j]]), 6, 7)
        ex <- paste0("exon", ZAHL)
        for (i in 1:length(id)){
          if (grepl(ex, x$HGVSp_Short[i])
              && grepl("ins", x$HGVSp_Short[i])
              && x$Hugo_Symbol[idh] == db$Gene[pta[j]]){
            out <- ""
            out <- paste0(db$Tumor.Type[pta[j]], ",EL:", db$Level[pta[j]],
                          ":", db$Drugs[pta[j]])
            if (onkb[id[i]] == FALSE){
              onkb[id[i]] <- out
            } else {
              onkb[id[i]] <- paste(onkb[id[i]], ";", out)
            }
          }
        }
      } else {
        for (i in 1:length(id)){
          if (db$Gene[pta[j]] == x$Hugo_Symbol[i]
              && grepl(db$Alteration[pta[j]], x$HGVSp_Short[i])){
            out <- ""
            out <- paste0(db$Tumor.Type[pta[j]], ",EL:", db$Level[pta[j]],
                          ":", db$Drugs[pta[j]])
            if (onkb[id[i]] == FALSE){
              onkb[id[i]] <- out
            } else {
              onkb[id[i]] <- paste(onkb[id[i]], ";", out)
            }
          }
        }
      }
    }
    idh_ret <- which(onkb != FALSE)
    idh <- idh[idh_ret]
    on_re <- onkb[idh_ret]
    if (length(idh)){
      x$oncokb[idh] <- on_re
    }
  }
  return(x)
}

rare <- function(x, maf = 0.001){
  #' rare
  #'
  #' @description Filters for rare mutations
  #'
  #' @param x dataframe. Table of Mutations
  #'
  #' @return returns x, dataframe. Table of Mutations
  #'
  #' @details Check the MAF-Score of the mutations. gnomAD score
  #' @details has to be lower than 0.001.
  #' @details Only the rare mutations are kept. 
  keep <- c()
  gnomad <- as.numeric(x$MAX_AF)
  if (nrow(x) > 0) {
    for (n in 1:length(gnomad)){
      if (is.na(gnomad[n]) | gnomad[n] <= maf){
        keep <- c(keep, n)
        next
      }
    }
    x <- x[keep, ] 
  }
  return(x)
}

## Condel Functions
condelQuery <- function(chr, start, ref, alt, dbfile){
  #' Condel Query
  #'
  #' @description Annotates Mutation with Condel
  #'
  #' @param chr string. Chromosom
  #' @param start numeric. Startposition
  #' @param ref string. Reference base
  #' @param alt string. Alternative base
  #' @param dbfile string. Filename for Database (textfile)
  #'
  #' @return returns logical. Condel score
  #'
  #' @details Condel predicts the effect of a mutational to the
  #' @details protein structur. The score is generated here.
  if (!(chr %in% c(1:22, "X", "Y"))) return(NA)

  param <- GRanges(c(chr), IRanges(c(start), c(start)))
  if (countTabix(dbfile, param = param) == 0) return(NA)

  res <- scanTabix(dbfile, param = param)
  dff <- Map(function(elt) {
    read.csv(textConnection(elt), sep = "\t", header = FALSE,
             colClasses = c("character"))
    }
    , res)
  dff <- dff[[1]]
  if (sum(dff$V4 == ref & dff$V5 == alt) == 0) return(NA)

  score <- as.numeric(dff$V15[dff$V4 == ref & dff$V5 == alt])
  if (sum(!is.na(score)) == 0) return(NA)

  return(max(score, na.rm = TRUE))
}
addCondel <- function(x, dbfile){
  #' Condel
  #'
  #' @description Add Condel prediction
  #'
  #' @param x dataframe. Table of Mutations
  #' @param dbfile string. Filename for Database (textfile)
  #'
  #' @return returns x dataframe. Table of Mutations with a new column.
  #'
  #' @details To get the Condel prediction of a mutation's effect on the
  #' @details protein structur, the position of the mutation is needed
  #' @details and then compared to its score in the database.
  chr_name <- gsub("chr", "", x$Chromosome)
  start_loc <- x$Start_Position
  if (class(start_loc) == "factor"){
    start_loc <- as.numeric(levels(start_loc))[start_loc]
  } else if (class(start_loc) == "character"){
    start_loc <- as.numeric(start_loc)
  }
  condelscore <- lapply(1:nrow(x), function(i)
    return(condelQuery(chr_name[i],
                       start_loc[i],
                       as.character(x$Reference_Allele[i]),
                       as.character(x$Allele[i]), dbfile)))

  condelscore <- unlist(condelscore)
  condellabel <- rep(NA, length(condelscore))
  isneutral <- condelscore < 0.52
  isneutral[is.na(isneutral)] <- FALSE
  isdel <- condelscore >= 0.52
  isdel[is.na(isdel)] <- FALSE
  condellabel[isneutral] <- "N"
  condellabel[isdel] <- "D"

  x$condel.score <- condelscore
  x$condel.label <- condellabel
  return(x)
}

target_check <- function(input, sureselect) {
  man <- read.delim(file = sureselect, header = FALSE, comment.char = "#")
  gr_man <- makeGRangesFromDataFrame(man[, c(1:3)], keep.extra.columns = FALSE,
                                    ignore.strand = TRUE, seqinfo = NULL, seqnames.field = "V1",
                                    start.field = "V2", end.field = "V3", starts.in.df.are.0based = FALSE)

  gr_x <- makeGRangesFromDataFrame(input, keep.extra.columns = FALSE, ignore.strand = TRUE,
                                seqinfo = NULL, seqnames.field = "Chromosome", start.field = "Start_Position",
                                end.field = "End_Position", starts.in.df.are.0based = FALSE)

  index <- findOverlaps(query = gr_man, subject = gr_x)
  index <- as.data.frame(index)
  output <- input[index[, 2], ]
  return(output)
}

txt2maf <- function(x){
  return(x[, c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Mutation_Status", "HGVSp_Short", "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")])
}

exclude <- function(x, vaf = 5){
  id <- which(as.numeric(x$t_AF) >= vaf)
  if(length(id) > 0) {
    x <- x[id, ]
  }
  return(x)
}

extract_lohs <- function(x) {
  return(x[x$n_GT == "0/1" & x$t_GT != "0/1", ])
}

loh_correction <- function(filt_loh, filt_gd = NULL, protocol = "somaticGermline", vaf = 5, actionable_genes = NA){
  table_id_loss <- which(as.numeric(filt_loh$table$t_AF) < 25 | is.na(as.numeric(filt_loh$table$t_AF)))
  filt_loh_loss <- filt_loh$table[table_id_loss, ]
  
  maf_filt_loh_loss <- filt_loh$maf[table_id_loss, ]
  if (length(table_id_loss) > 0) {
    filt_loh$table <- filt_loh$table[-table_id_loss, ]
    filt_loh$maf <- filt_loh$maf[-table_id_loss, ]
  }
  
  if (protocol == "somaticGermline"){
    filt_loh_loss$Variant_Reads <- filt_loh_loss$Count_Normal
    filt_loh_loss$Zygosity <- rep(x = "het", times = dim(filt_loh_loss)[1])
    id_exclude <- which(filt_loh_loss$n_AF < vaf)
    if (length(id_exclude) > 0){
      filt_loh_loss <- filt_loh_loss[-id_exclude, ]
      maf_filt_loh_loss <- maf_filt_loh_loss[-id_exclude, ]
    }
    id_ex2 <- which(colnames(filt_loh_loss) %in% c("Count_Normal", "Count_Tumor"))
    filt_loh_loss <- filt_loh_loss[, -id_ex2]
    filt_gd$table$Type <- rep("Germline", times = dim(filt_gd$table)[1])
    filt_loh_loss$Type <- rep("LoH", times = dim(filt_loh_loss)[1])
    if(dim(filt_loh_loss)[1] > 0) {
      filt_gd$table <- rbind(filt_gd$table, filt_loh_loss)
    }
    
    filt_gd$maf$Mutation_Status <- rep("Germline", times = dim(filt_gd$maf)[1])
    maf_filt_loh_loss$Mutation_Status <- rep("LoH", times = dim(maf_filt_loh_loss)[1])
    if(dim(maf_filt_loh_loss)[1] > 0) {
      filt_gd$maf <- rbind(filt_gd$maf, maf_filt_loh_loss)
    }

    if (!is.na(actionable_genes) & nrow(filt_gd$maf) > 0 & nrow(filt_gd$table) > 0) {
      filt_gd$maf <- actionable(filt_gd$maf, "Hugo_Symbol", actionable_genes)
      filt_gd$table <- actionable(filt_gd$table, "Hugo_Symbol", actionable_genes)
    }
    filt_gd$table <- filt_gd$table[!duplicated(filt_gd$table), ]
    filt_gd$maf <- filt_gd$maf[!duplicated(filt_gd$maf), ]
  }
  return(list(loh = filt_loh, gd = filt_gd))
}
