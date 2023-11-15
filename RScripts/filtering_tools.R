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
  idx <- grep(func, x$Func.refGene)
  if (length(idx) > 0) {
    x <- x[-idx, ]
  }
  return(x)
}

vrz <- function(x, mode, protocol){
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
  if (mode == "N" | mode == "T"){
    variant_frequency <- c()
    variant_count <- c()
    zygosity <- c()
    for (j in 1:dim(x)[1]) {
      other <- as.character(x[j, "Otherinfo"])
      if (protocol == "somaticGermline" | protocol == "somatic"){
          zygosity <- c(zygosity, substr(other, 1, 3))
          split <- strsplit(other, split = ":", fixed = TRUE)
          variant_frequency <- c(variant_frequency, split[[1]][12])
          variant_count <- c(variant_count, paste(as.numeric(split[[1]][11]),
                                                  as.numeric(split[[1]][9]),
                                                  sep = "|"))
        } else {
          zyg <- NA
          zyg = ifelse(test = grepl("HET=1", other), yes = "het", no = zyg)
          zyg = ifelse(test = grepl("HOM=1", other), yes = "hom", no = zyg)
          zygosity <- c(zygosity, zyg)
          split <- strsplit(other, split = ":", fixed = TRUE)
          variant_frequency <- c(variant_frequency, split[[1]][20])
          variant_count <- c(variant_count, paste(as.numeric(split[[1]][19]),
                                                  as.numeric(split[[1]][17]),
                                                  sep = "|"))
        }
    }
    x <- cbind(x, Variant_Allele_Frequency = variant_frequency,
               Zygosity = zygosity, Variant_Reads = variant_count)
    rownames(x) <- NULL
  } else if (mode == "LOH"){
    vaf.normal <- c()
    vaf.tumor <- c()
    count.normal <- c()
    count.tumor <- c()
    for (j in 1:dim(x)[1]) {
      other <- as.character(x[j, "Otherinfo"])
      split <- strsplit(other, split = ":", fixed = TRUE)
      vaf.normal <- c(vaf.normal, split[[1]][12])
      vaf.tumor <- c(vaf.tumor, split[[1]][18])
      count.normal <- c(count.normal, paste(split[[1]][11], split[[1]][9],
                                            sep = "|"))
      count.tumor <- c(count.tumor, paste(split[[1]][17], split[[1]][15],
                                          sep = "|"))
    } 
    x <- cbind(x, VAF_Normal = vaf.normal, VAF_Tumor = vaf.tumor,
               Count_Normal = count.normal, Count_Tumor = count.tumor)
    vaf.tumor <- as.character(vaf.tumor)
    vaf.tumor2 <- as.numeric(gsub("%", "", vaf.tumor))
  }
  return(x)
}

vrz_gatk <- function(x, mode, protocol = "Tumor_Normal", manifest){
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
  if(dim(x)[1] == 0){
    return(x)
  }
  if (mode == "N" | mode == "T"){
    split_other <- strsplit(x$Otherinfo13, split = ":", fixed = TRUE)
    x$Variant_Allele_Frequency <- unlist(lapply(split_other, function(f) f[3]))
    count_alt <- unlist(lapply(split_other, function(f) f[2]))
    x$Variant_Reads <- paste(unlist(lapply(strsplit(count_alt, split = ",", fixed = TRUE), function(f) f[2])),
                             unlist(lapply(split_other, function(f) f[4])), sep = "|")
    x$Zygosity <- "het"
    x$Zygosity[as.numeric(x$Variant_Allele_Frequency) > 0.75] <- "hom"
    return(x)
  }
  if (mode == "LOH"){
    vaf.normal <- c()
    vaf.tumor <- c()
    count.normal <- c()
    count.tumor <- c()
    for (j in 1:dim(x)[1]) {
      other <- as.character(x[j, "Otherinfo"])
      split <- strsplit(other, split = ":", fixed = TRUE)
      vaf.normal <- c(vaf.normal, split[[1]][12])
      vaf.tumor <- c(vaf.tumor, split[[1]][18])
      count.normal <- c(count.normal, paste(split[[1]][11], split[[1]][9],
                                            sep = "|"))
      count.tumor <- c(count.tumor, paste(split[[1]][17],
                                          (as.numeric(split[[1]][16]) 
                                           + as.numeric(split[[1]][17])),
                                          sep = "|"))
    }
    
    x <- cbind(x, VAF_Normal = vaf.normal, VAF_Tumor = vaf.tumor,
               Count_Normal = count.normal, Count_Tumor = count.tumor)
    vaf.tumor <- as.character(vaf.tumor)
    vaf.tumor2 <- as.numeric(gsub("%", "", vaf.tumor))
  }
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
  idx <- which(x$Gene.refGene %in% myhit)
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
  idx <- which (x$Gene.refGene %in% db$Hugo.Symbol[ts])
  x$is_tumorsuppressor[idx] <- 1

  x$is_oncogene <- 0
  og <- which(db$Is.Oncogene == "Yes")
  idx <- which (x$Gene.refGene %in% db$Hugo.Symbol[og])
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
  #' @note required in x: Gene.refGene, Start, AAChange.refGene
  # lisths <- read.delim(dbfile, header = T,
  #                      sep = "\t", colClasses = "character")
  lisths <- read_xls(path = dbfile, sheet = 1)
  x$is_hotspot <- 0
  idh <- which (x$Gene.refGene %in% lisths$Hugo_Symbol)
  phs <- which (lisths$Hugo_Symbol %in% x$Gene.refGene)

  if (length(idh) > 0){
    hotspot <- rep(FALSE, times = length(idh))
    # Check for matching Gene Locations and then for matching AAChange
    for (i in 1:length(idh)){
      for (j in 1:length(phs)){
        gp <- as.character(x$Start[idh][i])
        if (grepl(gp, as.character(lisths$Genomic_Position[phs][j]))){
          if (!is.na(x$AAChange.refGene[idh][i])){
            aac <- as.character(x$AAChange.refGene[idh][i])
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
  #' @note required in x: Gene.refGene, Start, AAChange.refGene
  # lisths <- read.delim(dbfile, header = T,
  #                      sep = "\t", colClasses = "character")
  lisths <- read_xls(path = dbfile, sheet = 2)

  # list should already habe a hotspot column for snps
  fs <- which(x$ExonicFunc.refGene != "frameshift deletion")
  fs2 <- which(x$ExonicFunc.refGene[fs] != "frameshift insertion")
  fs <- fs[fs2]

  # Check Genes
  idh <- which(x$Gene.refGene %in% lisths$Hugo_Symbol)
  phs <- which(lisths$Hugo_Symbol %in% x$Gene.refGene)
  idh <- intersect(idh, fs)

  if (length(idh) > 0){
    hotspot <- rep(FALSE, times = length(idh))
    # Check for matching Gene Locations and then for matching AAChange
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

            aac <- as.character(x$AAChange.refGene[idh[i]])
            a <- unlist(strsplit(aac, ","))
            id <- 5 * c(1:length(a))
            b <- unlist(strsplit(a, ":"))
            ast <- paste(b[id], collapse = ",")

            if (grepl(pos, as.character(ast))
                && (x$Gene.refGene[idh][i] == lisths$Hugo_Symbol[phs][j])
                && (grepl(type, as.character(x$AAChange.refGene[idh][i])))
                ){
              hotspot[i] <- paste0("p.", vaa)
              if (hotspot[i] != FALSE){
                break
                }
            }
          } else{
            pos <- substr(vaa, 2, l + 1)
            type <- substr(vaa, l + 2, nchar(vaa))
            aac <- as.character(x$AAChange.refGene[idh][i])

            a <- unlist(strsplit(aac, ","))
            id <- 5 * c(1:length(a))
            b <- unlist(strsplit(a, ":"))
            ast <- paste(b[id], sep = ",", collapse = "")

            if (grepl(pos, as.character(ast))
               && (x$Gene.refGene[idh][i] == lisths$Hugo_Symbol[phs][j])
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
  x$rvis <- as.numeric(db[x$Gene.refGene, 1])
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
    idx <- which(db$Gene == x[i, "Gene.refGene"])
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
    idx <- which(db$gene_name  == x[i, "Gene.refGene"])
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
  idh <- which(x$Gene.refGene %in% db$Gene)
  pta <- which(db$Gene %in% x$Gene.refGene[idh])
  onkb <- rep(FALSE, times = length(idh))
  if (length(pta)){
    for (j in 1:length(pta)){
      id <- which(x$Gene.refGene[idh] == db$Gene[pta[j]])
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
          if (grepl(ex, x$AAChange.refGene[idh[i]])
              && x$Gene.refGene[idh] == db$Gene[pta[j]]){
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
          if (grepl(ex, x$AAChange.refGene[i])
              && grepl("del", x$AAChange.refGene[i])
              && x$Gene.refGene[idh] == db$Gene[pta[j]]){
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
          if (grepl(ex, x$AAChange.refGene[i])
              && grepl("ins", x$AAChange.refGene[i])
              && x$Gene.refGene[idh] == db$Gene[pta[j]]){
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
          if (db$Gene[pta[j]] == x$Gene.refGene[i]
              && grepl(db$Alteration[pta[j]], x$AAChange.refGene[i])){
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

gene_name <- function(x){
  #' GeneName
  #'
  #' @description Finds GeneName for genes
  #'
  #' @param x dataframe. Table of Mutations
  #'
  #' @return returns x, dataframe. Table of Mutations with one extra column
  #'
  #' @details Find the whole genename for all the genes with a mutation and
  #' @details add it to an extra column.
  require(org.Hs.eg.db)
  x$GeneName <- rep(NA, times = dim(x)[1])

  sym <- org.Hs.egSYMBOL
  mapped_genes <- mappedkeys(sym)
  y <- unlist(as.list(sym[mapped_genes]))
  doubled1 <- grep(",", x$Gene.refGene)
  doubled2 <- grep(";", x$Gene.refGene)
  doubled <- c(doubled1, doubled2)

  for (i in 1:nrow(x)) {
    if (i %in% doubled) {
      if (i %in% doubled1) {
        gene.symbol <- strsplit(x$Gene.refGene[i], ",")[[1]]
      }
      if (i %in% doubled2) {
        gene.symbol <- strsplit(x$Gene.refGene[i], ";")[[1]]
      }
      entrez.id <- c()
      for (k in gene.symbol) {
        entrez.id <- c(entrez.id, names(y[which(y == k)]))
      }
    } else {
      gene.symbol <- x$Gene.refGene[i]
      entrez.id <- names(y[which(y == gene.symbol)])
    }
    x$GeneName[i] <- paste(unlist(mget(entrez.id, org.Hs.egGENENAME)),
                            collapse = ",")
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
  gnomad <- as.numeric(x$AF_popmax)
  for (n in 1:length(gnomad)){
    if (is.na(gnomad[n]) | gnomad[n] <= maf){
      keep <- c(keep, n)
      next
    }
  }
  x <- x[keep, ]
  return(x)
}

maneselect <- function(x, maneselectfile) {
  #' maneselect
  #'
  #' @description Selects the MANE Select preferred transcript if annovar did output multiple
  #'
  #' @param x dataframe. Table of Mutations
  #' @param maneselect file. File containing maneselect database
  #'
  #' @return returns x, dataframe. Table of Mutations
  mane <- read.table(maneselectfile, comment.char = "", sep = '\t', quote = "", header = T)
  mane <- mane[mane$MANE_status == "MANE Select", ]
  ref <- x$AAChange.refGene
  nms<-lapply(str_split(ref, ",", simplify = F), function(y) str_remove(str_remove(y, ":exon.*"), ".*:"))
  nm <- lapply(nms, function(y) unlist(lapply(y, function(z) mane[z == str_remove(mane$RefSeq_nuc, "\\.[0-9]+"), ]$RefSeq_nuc)))
  nm2<-data.frame(apply(str_split(ref, ",", simplify = T), 2, function(y) str_remove(str_remove(y, ":exon.*"), ".*:")))$X1
  if(!is.null(nm2)) {
    nm[sapply(nm, function(y) length(y) != 1)] <- nm2[sapply(nm, function(y) length(y) != 1)]
  }
  x$CChange <- "."
  prot <- str_remove_all(nm, "\\.[0-9]+")
  if(nrow(x) > 0) {
    for(i in 1:nrow(x)) {
      if(is.na(ref[i])) {
        next
      }
      transcript <- apply(str_split(ref[i], ",", simplify = T), 2, function(y) str_extract(y, paste0(prot[i], ".*")))
      transcript <- transcript[!is.na(transcript)]
      transplit <- str_split(transcript, ":", 5)[[1]]
      x$AAChange[i] <- transplit[length(transplit)]
      x$CChange[i] <- transplit[length(transplit) - 1]
      x$Exon[i] <- str_remove(transplit[2], "exon")
      x$Transcript[i] <- transplit[1]
    }
  }
  
  codon_nr <- str_extract(x$CChange, "[0-9]+")
  codon_ref_alt <- lapply(str_split(str_remove_all(x$CChange, "(c\\.|_|ins[ACGT]+|del|dup)"), "[0-9]+"), function(y) y[y != ""])
  is_snp <- unlist(lapply(codon_ref_alt, function(y) length(y) == 2))
  x$CChange[is_snp] <- paste0("c.", codon_nr[is_snp], unlist(lapply(codon_ref_alt[is_snp], function(y) y[1])), ">", unlist(lapply(codon_ref_alt[is_snp], function(y) y[2])))
  if(any(str_detect(x$CChange, "(del|dup)"))) {
    x[str_detect(x$CChange, "(del|dup)"),]$CChange <- str_remove_all(x[str_detect(x$CChange, "(del|dup)"),]$CChange, "[ACGT]+")
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
  chr_name <- gsub("chr", "", x$Chr)
  start_loc <- x$Start
  if (class(start_loc) == "factor"){
    start_loc <- as.numeric(levels(start_loc))[start_loc]
  } else if (class(start_loc) == "character"){
    start_loc <- as.numeric(start_loc)
  }
  condelscore <- lapply(1:nrow(x), function(i)
    return(condelQuery(chr_name[i],
                       start_loc[i],
                       as.character(x$Ref[i]),
                       as.character(x$Alt[i]), dbfile)))

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
                                seqinfo = NULL, seqnames.field = "Chr", start.field = "Start",
                                end.field = "End", starts.in.df.are.0based = FALSE)

  index <- findOverlaps(query = gr_man, subject = gr_x)
  index <- as.data.frame(index)
  output <- input[index[, 2], ]
  return(output)
}

txt2maf <- function(input, Center = center, refBuild = 'GRCh37', idCol = NULL, id = NULL, sep = "\t", Mutation_Status = c("T", "N","LOH")[1], protocol, snv_vcf, indel_vcf){
  
  #' Text to MAF Converter (cBioPortal Import)
  #' 
  #' @param input data.frame. Input data,frame from the pipeline. Output from the filtering scripts.
  #' @param Center string. Name of the processing Center.
  #' @param refBuild string. Reference genome version, e.g. GRCh37, GRCh38.
  #' @param idCol string. Column name in the input file containing the Patient_ID.
  #' @param id string. Patient_ID.
  #' @param sep string. Separator for the outputfile. Should be "\t" for cBioPortal import.
  #' @param Mutation_Status string. Kind of mutations about to process. Could either be T(umor), N(ormal) or LOH.
  #' @param protocol string. Used protocol either panel/tumor only (panelTumor) or WES / tumor and normal (somaticGerlmine or somatic)
  
  if (protocol == "somaticGermline" | protocol == "somatic"){
    if (Mutation_Status == "T") {
      Mutation_Status <- "Somatic"
    } else if (Mutation_Status == "N") {
      Mutation_Status <- "Germline"
    } else if (Mutation_Status == "LOH") {
      Mutation_Status <- "LoH"
    }
  }
  if (protocol == "panelTumor" | protocol == "tumorOnly"){
    Mutation_Status = "Tumor"
  }
  
  ann <- input
  
  essential.col <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",
                     "AAChange", "CChange", "avsnp150")
  
  for(i in 1:length(essential.col)){
    colId = suppressWarnings(grep(pattern = paste0("^",essential.col[i], "$"),
                                  x <- colnames(ann), ignore.case = TRUE))
    if (length(colId) == 1) {
      colnames(ann)[colId] <- essential.col[i]
    }
  }
  
  if(length(essential.col[!essential.col %in% colnames(ann)]) > 0){
    message("Available fields:")
    print(colnames(ann))
    message(paste0("Missing required field in input file: "))
    print(essential.col[!essential.col %in% colnames(ann)])
    stop()
  }
  
  if(is.null(idCol) & is.null(id)) {
    error('Provide either the column containing the Tumor Sample Barcode or the Tumor Sample Barcode as string!')
    stop()
  }
  
  if(is.null(idCol) & !is.null(id)) {
    ann$Tumor_Sample_Barcode <- paste(as.character(id),"TD",sep = "_")
    ann$Matched_Norm_Sample_Barcode <- paste(as.character(id),"GD",sep = "_")
  }
  
  if(!is.null(idCol)) {
    colnames(ann)[which(colnames(ann) == idCol)] <- "Tumor_Sample_Barcode"
  }
  
  if(is.null(Center)) {
    Center <- NA
  }
  
  ann.mand <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",
                "AAChange", "CChange", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "avsnp150")

  
  ann <- ann[, ann.mand]
  ann$ExonicFunc.refGene <- gsub(pattern = " SNV", replacement = "", x = ann$ExonicFunc.refGene)
  funcSpl <- strsplit(x = as.character(ann$ExonicFunc.refGene), split = ";", fixed = TRUE)
  funcSpl <- sapply(funcSpl, function(l) {l[length(l)]})
  ann$ExonicFunc.refGene <- funcSpl
  
  funcRef <- strsplit(x = as.character(ann$Func.refGene), split = ";", fixed = TRUE)
  funcRef <- sapply(funcRef, function(l) {l[length(l)]})
  ann$Func.refGene <- funcRef
  
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "intronic", yes = "Intron", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "intergenic", yes = "IGR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "downstream", yes = "3'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "upstream", yes = "5'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "splicing", yes = "Splice_Site", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "UTR3", yes = "3'UTR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "UTR5", yes = "5'UTR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "splice_region_variant&synonymous_variant", yes = "Splice_Site", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene %in% c("ncRNA_exonic", "ncRNA_intronic", "ncRNA_UTR3", "ncRNA_UTR5", "ncRNA"), yes = "RNA", no = ann$ExonicFunc.refGene)
  ann.lvls <- c("synonymous", "nonsynonymous", "stopgain", "stoploss", "frameshift insertion", "frameshift deletion", "nonframeshift insertion", "nonframeshift deletion",
                "Intron", "IGR", "Splice_Site", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "unknown", "UNKNOWN", "RNA", "nonframeshift substitution")
  ann.lbls <- c("Silent", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
                "Intron", "IGR", "Splice_Site", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "Unknown", "Unknown", "RNA", "Missense_Mutation")
  names(ann.lbls) <- ann.lvls
  ann$ExonicFunc.refGene <- as.character(ann.lbls[as.character(ann$ExonicFunc.refGene)])
  
  ann.del <- ann[ann$Alt %in% "-",]
  ann <- ann[!ann$Alt %in% "-",]
  
  if(nrow(ann.del) > 0){
    ann.del$var.type <- "DEL"
  } else {
    ann.del$var.type <- rep(NA, times = dim(ann.del)[1])
  }
  
  ann.ins <- ann[ann$Ref %in% "-",]
  ann <- ann[!ann$Ref %in% "-",]
  if(nrow(ann.ins) > 0){
    ann.ins$var.type <- "INS"
  } else {
    ann.ins$var.type <- rep(NA, times = dim(ann.ins)[1])
  }
  
  if(nrow(ann) > 0){
    ann$var.type <- "SNP"
  } else {
    ann$var.type <-  rep(NA, times = dim(ann)[1])
  }
  
  ann <- rbind(ann, ann.del, ann.ins)
   
  # Hugo Gene Symbol
  symbol <- unlist(lapply(strsplit(ann$Gene.refGene, split = ";"), function(x) {x[1]}))
  idx <- which(!is.na(symbol))
  
  # ENTREZ Gene ID
  entrez <- rep(NA, times = length(symbol))
  entrez[idx] <- unlist(lapply(mget(as.character(symbol[idx]), org.Hs.egSYMBOL2EG, ifnotfound = NA), function(x){x[1]}))
  
  # ENSEMBL Gene ID
  ensembl <- rep(NA, times = length(symbol))
  idx <- which(!is.na(entrez))
  ensembl[idx] <- unlist(lapply(mget(as.character(entrez[idx]), org.Hs.egENSEMBL, ifnotfound = NA), function(x){x[1]}))
  
  proteinChange <- as.character(ann$AAChange)
  
  # TxChange
  TxChange <- unlist(lapply(strsplit(x = as.character(ann$CChange), split = ";", fixed = T), function(x) x[1]))
  TxChange <- gsub(TxChange, pattern = " ", replacement = "")
  TxChange[TxChange == ""] <- NA
  
  # VAF
  vaf <- str_split(ifelse("Variant_Reads" %in% colnames(input), input$Variant_Reads, input$Count_Tumor), "\\|")
  t_ref_count <- as.numeric(unlist(lapply(vaf, function(x) x[2]))) - as.numeric(unlist(lapply(vaf, function(x) x[1])))
  t_alt_count <- as.numeric(unlist(lapply(vaf, function(x) x[1])))
  if("Count_Normal" %in% colnames(input)) {
    vaf <- str_split(input$Count_Normal, "\\|")
    n_ref_count <- as.numeric(unlist(lapply(vaf, function(x) x[2]))) - as.numeric(unlist(lapply(vaf, function(x) x[1])))
    n_alt_count <- as.numeric(unlist(lapply(vaf, function(x) x[1])))
  } else {
    n_ref_count <- vector(mode = "character", length = nrow(input))
    n_alt_count <- vector(mode = "character", length = nrow(input))
  }

  ann.maf <- data.table::data.table(Hugo_Symbol = as.character(symbol),
                                    Entrez_Gene_Id = as.character(entrez),
                                    Center = Center,
                                    NCBI_Build = refBuild,
                                    Chromosome = ann$Chr,
                                    Start_Position = ann$Start,
                                    End_Position = ann$End,
                                    Strand = "+",
                                    Variant_Classification = ann$ExonicFunc.refGene,
                                    Variant_Type = ann$var.type,
                                    Reference_Allele = ann$Ref,
                                    Tumor_Seq_Allele1 = ann$Ref,
                                    Tumor_Seq_Allele2 = ann$Alt,
                                    dbSNP_RS = ann$avsnp150,
                                    dbSNP_Val_Status = "",
                                    Tumor_Sample_Barcode = ann$Tumor_Sample_Barcode,
                                    Matched_Norm_Sample_Barcode = ann$Matched_Norm_Sample_Barcode,
                                    Match_Norm_Seq_Allele1 = ann$Ref,
                                    Match_Norm_Seq_Allele2 = ann$Alt,
                                    Tumor_Validation_Allele1 = "",
                                    Tumor_Validation_Allele2 = "",
                                    Match_Norm_Validation_Allele1 = "",
                                    Match_Norm_Validation_Allele2 = "",
                                    Verification_Status = "",
                                    Validation_Status = NA,
                                    Mutation_Status = Mutation_Status,
                                    Sequencing_Phase = "",
                                    Sequencing_Source = "",
                                    Validation_Method = "",
                                    Score = "",
                                    BAM_File = "",
                                    Sequencer = "",
                                    HGVSp_Short = proteinChange,
                                    Amino_Acid_Change = proteinChange,
                                    TxChange = TxChange,
                                    ENSEMBL_Gene_Id = ensembl,
                                    t_ref_count,
                                    t_alt_count,
                                    n_ref_count,
                                    n_alt_count)
  
  return(ann.maf)
}

txt2maf_mutect2 <- function(input, snv_vcf, protocol, Center = center, refBuild = 'GRCh37', idCol = NULL, id = NULL, sep = "\t", Mutation_Status = c("T", "N","LOH")[1]){
  
  #' Text to MAF Converter (cBioPortal Import)
  #' 
  #' @param input data.frame. Input data,frame from the pipeline. Output from the filtering scripts.
  #' @param Center string. Name of the processing Center.
  #' @param refBuild string. Reference genome version, e.g. GRCh37, GRCh38.
  #' @param idCol string. Column name in the input file containing the Patient_ID.
  #' @param id string. Patient_ID.
  #' @param sep string. Separator for the outputfile. Should be "\t" for cBioPortal import.
  #' @param Mutation_Status string. Kind of mutations about to process. Could either be T(umor), N(ormal) or LOH.
  #' @param protocol string. Used protocol either panel/tumor only (panelTumor) or WES / tumor and normal (somaticGerlmine or somatic)
  
  if (protocol == "somaticGermline" | protocol == "somatic"){
    if (Mutation_Status == "T") {
      Mutation_Status <- "Somatic"
    } else if (Mutation_Status == "N") {
      Mutation_Status <- "Germline"
    } else if (Mutation_Status == "LOH") {
      Mutation_Status <- "LoH"
    }
  }
  if (protocol == "panelTumor" | protocol == "tumorOnly" | protocol == "Tumor_only"){
    Mutation_Status = "Tumor"
  }
  
  ann <- input
  
  essential.col <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",
                     "AAChange", "CChange", "avsnp150")
  
  for(i in 1:length(essential.col)){
    colId = suppressWarnings(grep(pattern = paste0("^",essential.col[i], "$"),
                                  x <- colnames(ann), ignore.case = TRUE))
    if (length(colId) == 1) {
      colnames(ann)[colId] <- essential.col[i]
    }
  }
  
  if(length(essential.col[!essential.col %in% colnames(ann)]) > 0){
    message("Available fields:")
    print(colnames(ann))
    message(paste0("Missing required field in input file: "))
    print(essential.col[!essential.col %in% colnames(ann)])
    stop()
  }
  
  if(is.null(idCol) & is.null(id)) {
    error('Provide either the column containing the Tumor Sample Barcode or the Tumor Sample Barcode as string!')
    stop()
  }
  
  if(is.null(idCol) & !is.null(id)) {
    ann$Tumor_Sample_Barcode <- paste(as.character(id),"TD",sep = "_")
    ann$Matched_Norm_Sample_Barcode <- paste(as.character(id),"GD",sep = "_")
  }
  
  if(!is.null(idCol)) {
    colnames(ann)[which(colnames(ann) == idCol)] <- "Tumor_Sample_Barcode"
  }
  
  if(is.null(Center)) {
    Center <- NA
  }
  
  ann.mand <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene",
                "AAChange", "CChange", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "avsnp150")
  
  ann <- ann[, ann.mand]
  ann$ExonicFunc.refGene <- gsub(pattern = " SNV", replacement = "", x = ann$ExonicFunc.refGene)
  funcSpl <- strsplit(x = as.character(ann$ExonicFunc.refGene), split = ";", fixed = TRUE)
  funcSpl <- sapply(funcSpl, function(l) {l[length(l)]})
  ann$ExonicFunc.refGene <- funcSpl
  
  funcRef <- strsplit(x = as.character(ann$Func.refGene), split = ";", fixed = TRUE)
  funcRef <- sapply(funcRef, function(l) {l[length(l)]})
  ann$Func.refGene <- funcRef
  
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "intronic", yes = "Intron", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "intergenic", yes = "IGR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "downstream", yes = "3'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "upstream", yes = "5'Flank", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "splicing", yes = "Splice_Site", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "UTR3", yes = "3'UTR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "UTR5", yes = "5'UTR", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene == "splice_region_variant&synonymous_variant", yes = "Splice_Site", no = ann$ExonicFunc.refGene)
  ann$ExonicFunc.refGene <- ifelse(test = ann$Func.refGene %in% c("ncRNA_exonic", "ncRNA_intronic", "ncRNA_UTR3", "ncRNA_UTR5", "ncRNA"), yes = "RNA", no = ann$ExonicFunc.refGene)
  ann.lvls <- c("synonymous", "nonsynonymous", "stopgain", "stoploss", "frameshift insertion", "frameshift deletion", "nonframeshift insertion", "nonframeshift deletion",
                "Intron", "IGR", "Splice_Site", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "unknown", "UNKNOWN", "RNA", "nonframeshift substitution")
  ann.lbls <- c("Silent", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
                "Intron", "IGR", "Splice_Site", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "Unknown", "Unknown", "RNA", "Missense_Mutation")
  names(ann.lbls) <- ann.lvls
  ann$ExonicFunc.refGene <- as.character(ann.lbls[as.character(ann$ExonicFunc.refGene)])
  
  ann.del <- ann[ann$Alt %in% "-",]
  ann <- ann[!ann$Alt %in% "-",]
  
  if(nrow(ann.del) > 0){
    ann.del$var.type <- "DEL"
  } else {
    ann.del$var.type <- rep(NA, times = dim(ann.del)[1])
  }
  
  ann.ins <- ann[ann$Ref %in% "-",]
  ann <- ann[!ann$Ref %in% "-",]
  if(nrow(ann.ins) > 0){
    ann.ins$var.type <- "INS"
  } else {
    ann.ins$var.type <- rep(NA, times = dim(ann.ins)[1])
  }
  
  if(nrow(ann) > 0){
    ann$var.type <- "SNP"
  } else {
    ann$var.type <-  rep(NA, times = dim(ann)[1])
  }
  
  ann <- rbind(ann, ann.del, ann.ins)
  
  # Hugo Gene Symbol
  symbol <- unlist(lapply(strsplit(ann$Gene.refGene, split = ";"), function(x) {x[1]}))
  idx <- which(!is.na(symbol))
  
  # ENTREZ Gene ID
  entrez <- rep(NA, times = length(symbol))
  entrez[idx] <- unlist(lapply(mget(as.character(symbol[idx]), org.Hs.egSYMBOL2EG, ifnotfound = NA), function(x){x[1]}))
  
  # ENSEMBL Gene ID
  ensembl <- rep(NA, times = length(symbol))
  idx <- which(!is.na(entrez))
  ensembl[idx] <- unlist(lapply(mget(as.character(entrez[idx]), org.Hs.egENSEMBL, ifnotfound = NA), function(x){x[1]}))
  
  # Protein Change
  proteinChange <- as.character(ann$AAChange)
  
  # TxChange
  TxChange <- unlist(lapply(strsplit(x = as.character(ann$CChange), split = ";", fixed = T), function(x) x[1]))
  TxChange <- gsub(TxChange, pattern = " ", replacement = "")
  TxChange[TxChange == ""] <- NA  
  
  # VAF
  vaf <- str_split(ifelse("Variant_Reads" %in% colnames(input), input$Variant_Reads, input$Count_Tumor), "\\|")
  t_ref_count <- as.numeric(unlist(lapply(vaf, function(x) x[2]))) - as.numeric(unlist(lapply(vaf, function(x) x[1])))
  t_alt_count <- as.numeric(unlist(lapply(vaf, function(x) x[1])))
  if("Count_Normal" %in% colnames(input)) {
    vaf <- str_split(input$Count_Normal, "\\|")
    n_ref_count <- as.numeric(unlist(lapply(vaf, function(x) x[2]))) - as.numeric(unlist(lapply(vaf, function(x) x[1])))
    n_alt_count <- as.numeric(unlist(lapply(vaf, function(x) x[1])))
  } else {
    n_ref_count <- vector(mode = "character", length = nrow(input))
    n_alt_count <- vector(mode = "character", length = nrow(input))
  }

  ann.maf <- data.table::data.table(Hugo_Symbol = as.character(symbol),
                                    Entrez_Gene_Id = as.character(entrez),
                                    Center = Center,
                                    NCBI_Build = refBuild,
                                    Chromosome = ann$Chr,
                                    Start_Position = ann$Start,
                                    End_Position = ann$End,
                                    Strand = "+",
                                    Variant_Classification = ann$ExonicFunc.refGene,
                                    Variant_Type = ann$var.type,
                                    Reference_Allele = ann$Ref,
                                    Tumor_Seq_Allele1 = ann$Ref,
                                    Tumor_Seq_Allele2 = ann$Alt,
                                    dbSNP_RS = ann$avsnp150,
                                    dbSNP_Val_Status = "",
                                    Tumor_Sample_Barcode = ann$Tumor_Sample_Barcode,
                                    Matched_Norm_Sample_Barcode = ann$Matched_Norm_Sample_Barcode,
                                    Match_Norm_Seq_Allele1 = ann$Ref,
                                    Match_Norm_Seq_Allele2 = ann$Alt,
                                    Tumor_Validation_Allele1 = "",
                                    Tumor_Validation_Allele2 = "",
                                    Match_Norm_Validation_Allele1 = "",
                                    Match_Norm_Validation_Allele2 = "",
                                    Verification_Status = "",
                                    Validation_Status = NA,
                                    Mutation_Status = Mutation_Status,
                                    Sequencing_Phase = "",
                                    Sequencing_Source = "",
                                    Validation_Method = "",
                                    Score = "",
                                    BAM_File = "",
                                    Sequencer = "",
                                    HGVSp_Short = proteinChange,
                                    Amino_Acid_Change = proteinChange,
                                    TxChange = TxChange,
                                    ENSEMBL_Gene_Id = ensembl,
                                    t_ref_count,
                                    t_alt_count,
                                    n_ref_count,
                                    n_alt_count)
  
  return(ann.maf)
}

exclude <- function(x, vaf = 5){
  variant_freq <- substr(as.character(x$Variant_Allele_Frequency), start = 1,
                         stop = nchar(as.character(x$Variant_Allele_Frequency))-1)
  id <- which(as.numeric(variant_freq) >= vaf)
  if(length(id) > 0) {
    x <- x[id, ]
  }
  return(x)
}

normalize_vaf_gatk <- function(x){
  x$Variant_Allele_Frequency <- as.numeric(x$Variant_Allele_Frequency) * 100
  x$Variant_Allele_Frequency <- paste0(x$Variant_Allele_Frequency, "%")
  return(x)
}

loh_correction <- function(filt_loh, filt_gd = NULL, protocol = "somaticGermline", vaf = 10, actionable_genes = NA){
  filt_loh$table$VAF_Tumor <- as.character(filt_loh$table$VAF_Tumor)
  filt_loh$table$VAF_Normal <- as.character(filt_loh$table$VAF_Normal)
  table_id_loss <- which(as.numeric(substr(filt_loh$table$VAF_Tumor, start = 1,
                                           stop = nchar(filt_loh$table$VAF_Tumor) - 1)) < 20)
  filt_loh_loss <- filt_loh$table[table_id_loss, ]
  filt_loh$table <- filt_loh$table[-table_id_loss, ]
  
  maf_vaf_tumor <-100*as.numeric(filt_loh$maf$t_alt_count) / (as.numeric(filt_loh$maf$t_alt_count) + as.numeric(filt_loh$maf$t_ref_count))
  maf_vaf_normal <- 100*as.numeric(filt_loh$maf$n_alt_count) / (as.numeric(filt_loh$maf$n_alt_count) + as.numeric(filt_loh$maf$n_ref_count))
  maf_id_loss <- which(maf_vaf_tumor < 20)
  maf_filt_loh_loss <- filt_loh$maf[maf_id_loss, ]
  filt_loh$maf <- filt_loh$maf[-maf_id_loss, ]
  
  if (protocol == "somaticGermline"){
    filt_loh_loss$Variant_Reads <- filt_loh_loss$Count_Normal
    filt_loh_loss$Variant_Allele_Frequency <- filt_loh_loss$VAF_Normal
    filt_loh_loss$Zygosity <- rep(x = "het", times = dim(filt_loh_loss)[1])
    id_hom <- which(as.numeric(substr(filt_loh_loss$VAF_Normal, start = 1,
                                      stop = nchar(filt_loh_loss$VAF_Normal) - 1)) > 75)
    id_exclude <- which(as.numeric(substr(filt_loh_loss$VAF_Normal, start = 1,
                                          stop = nchar(filt_loh_loss$VAF_Normal) - 1)) < 10)
    if (length(id_hom) > 0){
      filt_loh_loss$Zygosity[id_hom] <- "hom"
    }
    if (length(id_exclude) > 0){
      filt_loh_loss <- filt_loh_loss[-id_exclude, ]
    }
    id_ex2 <- which(colnames(filt_loh_loss) %in% c("VAF_Normal", "VAF_Tumor", "Count_Normal", "Count_Tumor"))
    filt_loh_loss <- filt_loh_loss[, -id_ex2]
    filt_gd$table$Type <- rep("Germline", times = dim(filt_gd$table)[1])
    filt_loh_loss$Type <- rep("LoH", times = dim(filt_loh_loss)[1])
    if(dim(filt_loh_loss)[1] > 0) {
      filt_gd$table <- rbind(filt_gd$table, filt_loh_loss)
    }
    
    normal_vaf <- 100*as.numeric(maf_filt_loh_loss$n_alt_count) / (as.numeric(maf_filt_loh_loss$n_alt_count) + as.numeric(maf_filt_loh_loss$n_ref_count))
    maf_id_exclude <- which(normal_vaf < 10)
    if (length(maf_id_exclude) > 0){
      maf_filt_loh_loss <- maf_filt_loh_loss[-maf_id_exclude, ]
    }
    filt_gd$maf$Mutation_Status <- rep("Germline", times = dim(filt_gd$maf)[1])
    maf_filt_loh_loss$Mutation_Status <- rep("LoH", times = dim(maf_filt_loh_loss)[1])
    if(dim(maf_filt_loh_loss)[1] > 0) {
      filt_gd$maf <- rbind(filt_gd$maf, maf_filt_loh_loss)
    }

    if (!is.na(actionable_genes) & nrow(filt_gd$maf) > 0 & nrow(filt_gd$table) > 0) {
      filt_gd$maf <- actionable(filt_gd$maf, "Hugo_Symbol", actionable_genes)
      filt_gd$table <- actionable(filt_gd$table, "Gene.refGene", actionable_genes)
    }
  } else {
    filt_gd = filt_gd
  }
  return(list(loh = filt_loh, gd = filt_gd))
}
