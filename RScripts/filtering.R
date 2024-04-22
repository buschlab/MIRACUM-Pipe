#### Filter Variants Function

filtering <- function(
  snpfile,
  path_data,
  path_script,
  mode = "T",
  center = "Freiburg",
  id = id,
  protocol,
  sureselect,
  vaf = 5,
  min_var_count = 20,
  maf = 0.01,
  actionable_genes = NA,
  covered_exons = covered_exons,
  cov_t = 1,
  sureselect_type,
  maneselectfile
  ) {
  #' Filter Variants
  #'
  #' @description Filters the somatic SNPs and InDel for analysis
  #'
  #' @param snpfile dataframe. Table of SNVs
  #' @param path_data string. Directory of the databases
  #' @param path_script string. Directory of required scripts
  #' @param sureselect string. Kind of sequencer
  #' @param mode string. Mode for filtering: "T", "N", "LOH"
  #'
  #' @return returns list of
  #' @return Table dataframe. Filtered table of mutations
  #' @return tmb numerical. Tumor mutational burden
  #'
  #' @note Please make sure that the columns of your input and your databases
  #' @note have the right names.
  #'
  #' @details The mutations are filtered to find the pathogenic mutations.
  #' @details First only the mutations passing all the quality filters are
  #' @details considered. Then the tumor mutational burden is calculated in 
  #' @details Tumor mode ("T"). Afterwards the filter for the functionality of
  #' @details the genetic region is applied and the mutations with the desired
  #' @details exonic functions are chosen. Only rare mutations are kept.
  #' @details In the following some additional information is extraced of the
  #' @details table: in Normal and Tumor mode ("N", "T") Variant Allele
  #' @details Frequency (VAF), Readcounts and Zygosity, in LoH mode ("LOH") VAF
  #' @details and Readcounts both for Tumor and Normal.
  #' @details The full gene names are added to the dataframe.
  #' @details In Addition to that some database queries are done and finally
  #' @details the canonical aminoacid and base changes are added.
  #' @details As last step Condel prediction scores are written into another
  #' @details column.
  require(stringr)
  require(openxlsx)
  require(Rsamtools)
  require(GenomicRanges)

  # Read Data
  x <- read.delim(file = snpfile, header = T, stringsAsFactors = F, comment.char = "#")
  x$HGVSp_Short <- str_replace_all(x$HGVSp_Short, "%3D", "=")
  
  # Filter for LoHs based on VarScan algorithm
  if (mode == "LOH") {
    x <- extract_lohs(x) 
  }
  
  # Filter for actionable genes in Germline
  if (protocol == "somaticGermline" && mode == "N") {
    x <- actionable(x, "Hugo_Symbol", actionable_genes)
  }
  
  # Filter for targeted region in tNGS
  x <- target_check(x, sureselect)

  # Calculate covered region from bed file
  cov_region <- covered_region(sureselect = sureselect, mode = mode)
  

  # Quality Filter
  if (mode == "T") {
    id.pass <- grep("PASS", x$FILTER)
    if (length(id.pass) > 0) {
      x <- x[id.pass, ]
    } else {
      stop("No variant passed quality filter!")
    }
  }
  
  # Extract VAF, Readcounts (and Zygosity)
  x <- vrz(x = x, mode = mode, protocol = protocol)
      
  # VAF Filter
  x <- exclude(x, vaf = vaf)
    
  # Remove Variants with Variant Read Count below 4/20
  x <- mrc(x = x, min_var_count = min_var_count)
    
  # Filter for exonic function
  x <- filt(x, "Intron")
  x <- filt(x, "3'UTR")
  x <- filt(x, "5'UTR")
  x <- filt(x, "RNA")
  x <- filt(x, "IGR")
  x <- filt(x, "3'Flank")
  x <- filt(x, "5'Flank")
  
  # TMB calculation for WES
  if (protocol == "somaticGermline" | protocol == "somatic"){
    if (mode == "T") {
      id_ex <- which(x$BIOTYPE == "protein_coding")
      x_coding <- x[id_ex, ]
      tmb <- tmb_ex(x_coding, covered_exons, mode = "T", cov_t)
    } else {
      tmb <- list(tmb = NULL, exon_region = NULL, used_exon_region = NULL,
                  number_used_mutations_tmb = NULL)
    }
  }
  
  # Filter for exonic function
  test <- as.character(x$Variant_Classification)
  syn.snv <- which(test %in% c("Silent", "Splice_Region"))
  if (length(syn.snv) > 0) {
    x <- x[- syn.snv,]
  }
  
  # Filter for rare mutations
  x <- rare(x, maf = maf)
  
  # TMB calculation for panel (TSO500) and tumorOnly WES; only rare mutations; assumption: rare mutations are "somatic"
  if (protocol == "panelTumor" | protocol == "tumorOnly") {
    if (mode == "T") {
      id_ex <- which(x$BIOTYPE == "protein_coding")
      x_coding <- x[id_ex, ]
      tmb <- tmb_ex(x_coding, covered_exons, mode = "T", cov_t)
    } else {
      tmb <- list(tmb = NULL, exon_region = NULL, used_exon_region = NULL,
                  number_used_mutations_tmb = NULL)
    }
  }

  if (dim(x)[1] != 0) {

    # Include GeneName
    x$Start_Position <- as.numeric(as.character(x$Start_Position))
    x$End_Position <- as.numeric(as.character(x$End_Position))
    x$Hugo_Symbol <- as.character(x$Hugo_Symbol)
    
    # Further Annotation
    # Database Queries
    x <- isflag(x, dbfile = paste(path_data, "flag_genes.txt", sep = "/"))
    x <- isogtsg(x, dbfile = paste(path_data, "cancerGeneList.tsv", sep = "/"))
    x <- ishs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- isihs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- rvis(x, paste(path_data, "RVIS_score.txt", sep = "/"))
    x <- trgt(x, paste(path_data, "TARGET_db.txt", sep = "/"))
    x <- dgidb(x, paste(path_data, "DGIdb_interactions.tsv", sep = "/"))
    x <- oncokb(x, paste(path_data, "oncokb_biomarker_drug_associations.tsv", sep = "/"))
    x <- addCondel(x, paste(path_data, "fannsdb.tsv.gz", sep = "/"))   
    
    # MAF
    out.maf <- txt2maf(x)
    
    return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )

  } else if (mode == "N" | mode == "T") {
    print("No SNVs passed filter!")
    out.maf <- data.frame()
   return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )

  } else if (mode == "LOH") {
    print("No LOH passed filter!")
    out.maf <- data.frame()
    return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )
  }
}
