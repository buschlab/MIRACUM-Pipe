## key_Results
keys <- function(
  mut_sig,
  mutation_analysis_result,
  mutation_analysis_result_gd,
  cnv_analysis_results,
  filt_result_td,
  filt_result_gd,
  med_tmb,
  protocol,
  fusions
) {

  if (sureselect_type %in% c("V5UTR", "V6UTR", "V6")){
    sureselect_type <- paste("Agilent SureSelect", sureselect_type, sep = " ")
  } else if (sureselect_type %in% c("TSO500", "TruSight_Tumor")) {
    sureselect_type <- paste("Illumina", sureselect_type, sep = " ")
  }

  if (protocol == "somaticGermline" | protocol == "somatic") {
    if (mutation_analysis_result$msi < 10) {
      msi_helper <- "MSS"
    } else {
      msi_helper <- "MSI"
    }
    brca_helper <- which(mut_sig$output$Summary$Signature ==  "AC3")
    if (length(brca_helper) == 1 & mut_sig$output$Summary["AC3", 3] > 1.0) {
      brca_helper <- paste0(
        round(
          mut_sig$output$Summary["AC3", 3],
          digits = 1
        ),
        " (", round(mut_sig$output$Summary["AC3", 4], digits = 1),
        ";", round(mut_sig$output$Summary["AC3", 5], digits = 1) , ")"
      )
    } else {
      brca_helper <- "<1%"
    }

    bammatcher_common <- mutation_analysis_result$bam_matcher$FracCommon * 100
    bammatcher_helper <-  paste0("Tumor/Normal passen nicht zusammen (", bammatcher_common, "%).")
    bammatcher_count = mutation_analysis_result$bam_matcher$Same + mutation_analysis_result$bam_matcher$Different
    bammatcher_count_plus = bammatcher_count + max(mutation_analysis_result$bam_matcher$X1het.2sub, mutation_analysis_result$bam_matcher$X1sub.2het)
    if (bammatcher_count <= 20) {
      bammatcher_helper <- "Keine Bestimmung möglich."
    } else if (bammatcher_count <= 100) {
      if (bammatcher_count_plus >= 0.8) {
        bammatcher_helper <- paste0("Tumor/Normal passen zusammen (", bammatcher_common, "%).")
      }
    } else {
      if (mutation_analysis_result$bam_matcher$FracCommon >= 0.8) {
        bammatcher_helper <- paste0("Tumor/Normal passen zusammen (", bammatcher_common, "%).")
      }
    }

  }
  if (protocol == "panelTumor" | protocol == "tumorOnly") {
    if (mutation_analysis_result$msi < 20) {
      msi_helper <- "MSS"
    } else {
      msi_helper <- "MSI"
    }
    brca_helper <- which(mut_sig ==  "AC3")
    if (length(brca_helper) == 1 & mut_sig["AC3", 3]*100 > 1) {
      brca_helper <- paste0(round(mut_sig["AC3", 3]*100, digits = 1))
    } else {
      brca_helper <- "<1%"
    }
  }
  if (protocol == "panelTumor") {
    if (!is.null(fusions)) {
      fus_tmp <- fusions$Fusion_OV$Fusionen
      fus_tmp <- paste(fus_tmp, collapse = ", ")
    } else {
      fus_tmp <- "Keine"
    }
  }

  if (!is.null(med_tmb$med) & !is.null(med_tmb$sd)) {
    tmb_helper <- paste0(
      med_tmb$med, " (",
      as.numeric(med_tmb$sd[1], digits = 2),
      "-", as.numeric(med_tmb$sd[2], digits = 2),
      ")"
    )
  } else if (!is.null(med_tmb$med)) {
    tmb_helper <- med_tmb$med
  } else {
    tmb_helper <- "-"
  }
  if (protocol == "somaticGermline") {
    mut_tab1 <- data.frame(
      Eigenschaften = c(
        "Capture Kit",
        "Abgedeckte Region (total)",
        "Abgedeckte Region (exonisch)",
        paste0("Mutationslast (exonisch, VAF > ", vaf, "%)"),
        paste0("Mittlere TMB der Entit\"at", " (", entity, ")"),
        paste0("Anzahl somatischer Mutationen inkl. LoH (VAF > ", vaf, "%)"),
        paste0("BRCAness (%) inkl. KI", " (VAF > ", vaf, "%)"),
        "Mikrosatelliten Status (Score)",
        "HRD-Score (HRD-LoH|TAI|LST)",
        "bioinformatischer Tumorzellgehalt (%)",
        "bioinformatischer Abgleich",
        "Ploidie",
        "Anzahl CN- Regionen",
        paste0("Anzahl seltener Keimbahnmutationen (VAF > ", germlineVaf, "%)")
      ), Wert = c(
        as.character(sureselect_type),
        paste(
          round(
            x = as.numeric(filt_result_td$covered_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste(
          round(
            x = as.numeric(filt_result_td$exon_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste0(
          round(
            x = filt_result_td$tmb, digits = 2
          ),
          "/Mb", " (", filt_result_td$number_used_mutations_tmb, "/",
          round(x = as.numeric(filt_result_td$used_exon_region),
            digits = 2
          ), ")"
        ),
        tmb_helper,
        as.character(round(
          x = sum(as.numeric(mutation_analysis_result$mut_tab[, 2])),
          digits = 0)
        ),
        brca_helper,
        paste(msi_helper," (", mutation_analysis_result$msi, "%)", sep = ""),
        cnv_analysis_results$hrd$score,
        cnv_analysis_results$purity$purity*100,
        bammatcher_helper,
        cnv_analysis_results$purity$ploidy,
        paste0(
          round(
            x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1],
            digits = 0
          ), " Regionen"
        ),
        as.character(round(
          x = sum(as.numeric(mutation_analysis_result_gd$mut_tab[, 3])),
          digits = 0)
        )
      )
    )
  }
    if (protocol == "somatic") {
    mut_tab1 <- data.frame(
      Eigenschaften = c(
        "Capture Kit",
        "Abgedeckte Region (total)",
        "Abgedeckte Region (exonisch)",
        paste0("Mutationslast (exonisch, VAF > ", vaf, "%)"),
        paste0("Mittlere TMB der Entit\"at", " (", entity, ")"),
        paste0("Anzahl somatischer Mutationen inkl. LoH (VAF > ", vaf, "%)"),
        paste0("BRCAness (%) inkl. KI", " (VAF > ", vaf, "%)"),
        "Mikrosatelliten Status (Score)",
        "HRD-Score (HRD-LoH|TAI|LST)",
        "bioinformatischer Tumorzellgehalt (%)",
        "bioinformatischer Abgleich",
        "Ploidie",
        "Anzahl CN- Regionen"
      ), Wert = c(
        as.character(sureselect_type),
        paste(
          round(
            x = as.numeric(filt_result_td$covered_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste(
          round(
            x = as.numeric(filt_result_td$exon_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste0(
          round(
            x = filt_result_td$tmb, digits = 2
          ),
          "/Mb", " (", filt_result_td$number_used_mutations_tmb, "/",
          round(x = as.numeric(filt_result_td$used_exon_region),
            digits = 2
          ), ")"
        ),
        tmb_helper,
        as.character(round(
          x = sum(as.numeric(mutation_analysis_result$mut_tab[, 2])),
          digits = 0)
        ),
        brca_helper,
        paste(msi_helper," (", mutation_analysis_result$msi, "%)", sep = ""),
        cnv_analysis_results$hrd$score,
        cnv_analysis_results$purity$purity*100,
        bammatcher_helper,
        cnv_analysis_results$purity$ploidy,
        paste0(
          round(
            x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1],
            digits = 0
          ), " Regionen"
        )
      )
    )
  }
  if (protocol == "panelTumor") {
    mut_tab1 <- data.frame(
      Eigenschaften = c(
        "Capture Kit",
        "Abgedeckte Region (total)",
        "Abgedeckte Region (exonisch)",
        paste0("Mutationslast (exonisch, VAF > ", vaf, "%)"),
        paste0("Mittlere TMB der Entit\"at", " (", entity, ")"),
        paste0("Anzahl der Mutationen (VAF > ", vaf, "%)"),
        paste0("BRCAness (%)", " (VAF > ", vaf, "%)"),
        "Mikrosatelliten Status (Score)",
        "HRD-Score (HRD-LoH|TAI|LST)",
        "bioinformatischer Tumozellgehalt (%)",
        "Ploidie",
        "Anzahl CN- Regionen",
        "Fusionen"
      ), Wert = c(
        as.character(sureselect_type),
        paste(
          round(
            x = as.numeric(filt_result_td$covered_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste(
          round(
            x = as.numeric(filt_result_td$exon_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste0(
          round(
            x = filt_result_td$tmb, digits = 2
          ),
          "/Mb", " (", filt_result_td$number_used_mutations_tmb, "/",
          round(
            x = as.numeric(filt_result_td$used_exon_region),
            digits = 2
          ), ")"
        ),
        tmb_helper,
        as.character(round(
          x = sum(as.numeric(mutation_analysis_result$mut_tab[, 3])),
          digits = 0)
        ),
        brca_helper,
        paste(msi_helper," (", mutation_analysis_result$msi, "%)", sep = ""),
        cnv_analysis_results$hrd$score,
        cnv_analysis_results$purity$purity*100,
        cnv_analysis_results$purity$ploidy,
        paste0(
          round(
            x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1],
            digits = 0
          ), " Regionen"
        ),
        fus_tmp
      )
    )
  }
  if (protocol == "tumorOnly") {
    mut_tab1 <- data.frame(
      Eigenschaften = c(
        "Capture Kit",
        "Abgedeckte Region (total)",
        "Abgedeckte Region (exonisch)",
        paste0("Mutationslast (exonisch, VAF > ", vaf, "%)"),
        paste0("Mittlere TMB der Entit\"at", " (", entity, ")"),
        paste0("Anzahl der Mutationen (VAF > ", vaf, "%)"),
        paste0("BRCAness (%)", " (VAF > ", vaf, "%)"),
        "Mikrosatelliten Status (Score)",
        "HRD-Score (HRD-LoH|TAI|LST)",
        "bioinformatischer Tumorzellgehalt (%)",
        "Ploidie",
        "Anzahl CN- Regionen"
      ), Wert = c(
        as.character(sureselect_type),
        paste(
          round(
            x = as.numeric(filt_result_td$covered_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste(
          round(
            x = as.numeric(filt_result_td$exon_region),
            digits = 2
          ), "Mb", sep = ""
        ),
        paste0(
          round(
            x = filt_result_td$tmb, digits = 2
          ),
          "/Mb", " (", filt_result_td$number_used_mutations_tmb, "/",
          round(
            x = as.numeric(filt_result_td$used_exon_region),
            digits = 2
          ), ")"
        ),
        tmb_helper,
        as.character(round(
          x = sum(as.numeric(mutation_analysis_result$mut_tab[, 3])),
          digits = 0)
        ),
        brca_helper,
        paste(msi_helper," (", mutation_analysis_result$msi, "%)", sep = ""),
        cnv_analysis_results$hrd$score,
        cnv_analysis_results$purity$purity*100,
        cnv_analysis_results$purity$ploidy,
        paste0(
          round(
            x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1],
            digits = 0
          ), " Regionen"
        )
      )
    )
  }
  return(mut_tab1)
}

l_gen_nex <- function(df, type = "SNV") {
  if ("Hugo_Symbol" %in% colnames(df)) {
    df$Symbol <- df$Hugo_Symbol
  }
  if ( type == "SNV") {
    vec_gen_nex <- paste0(
      "{\\href{https://www.genomenexus.org/variant/",
      df$Chr, ":g.", df$Start, df$Ref, "\\%3E",
      df$Alt, "}{", df$Symbol, "}}"
    )
  } else if (type == "DEL") {
    vec_gen_nex <- paste0(
      "{\\href{https://www.genomenexus.org/variant/",
      df$Chr, ":g.", df$Start,
      df$Ref, "\\%3Edel}{",
      df$Symbol, "}}"
    )
  } else if (type == "INS") {
    vec_gen_nex <- paste0(
      "{\\href{https://www.genomenexus.org/variant/",
      df$Chr, ":g.", df$Start,
      "_", (df$Start + nchar(as.character(df$Alt))),
     "ins\\%3El", df$Alt, "}{", df$Symbol, "}}"
    )
  }
  return(vec_gen_nex)
}

meta <- function(df) {
  df$HGVSp_Short <- substr(x = df$HGVSp_Short, start = 3, stop = nchar(df$HGVSp_Short))
  hyper_refs <- paste0(
    "https://search.cancervariants.org/\\#", df$Symbol, "\\%20",
    df$HGVSp_Short
  )
  vec_meta <- paste0(
    "{\\href{", hyper_refs, "}{", df$HGVSp_Short, "}}"
  )
  vec_meta <- gsub(
    pattern = "_", replacement = "\\_",
    x = vec_meta, fixed = TRUE
  )
  return(vec_meta)
}

varsome <- function(df, mode) {
  baseURL <- "https://varsome.com/variant/hg19/"
  if (mode == "1") {
    completeVarsomeLink <- paste0(
      "{\\href{", baseURL, df$Chromosome, "\\%3A",df$Start_Position,
      "\\%3A", df$Reference_Allele, "\\%3A", df$Allele, "}{Link}}"
    )
  } else if (mode == "2"){
    completeVarsomeLink <- paste0(
      "{\\href{", baseURL, df$Chromosome, "\\%3A",df$Start_Position,
      "\\%3A", df$Reference_Allele, "\\%3A", df$Allele, "}{", df$Variant_Classification, "}}"
    )
  }
    completeVarsomeLink <- gsub(
      pattern = "-", replacement = "", fixed = TRUE,
      x = completeVarsomeLink
    )
  return(completeVarsomeLink)
}

acmg <- function(df) {
  
  # Split by ,
  df$CLIN_SIG <- str_split(df$CLIN_SIG, ",")
  
  # Remove multi annotations
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "[a-z_]+/[a-z_]+", ""))
  
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^pathogenic$", "5"))
  
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^likely_pathogenic$", "4"))
  
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^conflicting_interpretations_of_pathogenicity$", "4*"))

  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^uncertain_significance$", "3"))
  
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^likely_benign$", "2"))
  
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^benign$", "1"))

  # Remove annotations without acmg relevance
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^other$"))
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^drug_response$"))
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^risk_factor$"))
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^not_provided$"))
  
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace(x, "^$", "0"))
  prio <- lapply(df$CLIN_SIG, function(x) str_replace(x, "4\\*", "4"))
  argmax <- unlist(lapply(prio, function(x) which.max(as.numeric(x))))
  df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace(x, "0", "."))
  
  return(mapply(function(x,y) x[y], df$CLIN_SIG, argmax))
}

revel <- function(df) {
  if("REVEL" %in% colnames(df)) {
    vec_rev <- round(as.numeric(df$REVEL), digits = 1)
  } else {
    vec_rev <- round(as.numeric(df$REVEL_score), digits = 1)
  }
  vec_rev[is.na(vec_rev)] <- "."
  vec_cat <- ifelse(test = vec_rev > 0.5, yes = "D", no = "N")
  vec_cat <- paste0(vec_cat, " (", vec_rev, ")")
  vec_cat[vec_rev == "."] <- "."
  
  return(list(revel = vec_rev, cat = vec_cat))
}

cosmic <- function(df) {
  vec_cos <- str_match_all(df$Existing_variation, "COSV[0-9]+")
  vec_cos <- unlist(lapply(vec_cos, function(x) paste(unlist(x), collapse=', ')))
  vec_cos <- str_remove_all(vec_cos, "COSV")
  vec_cos[nchar(vec_cos) == 0] <- "."
  return(vec_cos)
}

ex_func <- function(df) {
  vec_exf <- df$Variant_Classification

  vec_exf <- gsub(
    pattern = "Nonsense_Mutation",
    replacement = "stopG", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "Nonstop_Mutation",
    replacement = "stopL", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "Translation_Start_Site",
    replacement = "startL", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "Frame_Shift_Del",
    replacement = "fsDel", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "Frame_Shift_Ins",
    replacement = "fsIns", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "In_Frame_Del",
    replacement = "nfsDel", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "In_Frame_Ins",
    replacement = "nfsIns", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "Splice_Site",
    replacement = "splice", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "Missense_Mutation",
    replacement = "nsSNV", x = vec_exf
  )
  
  return(vec_exf)
}

highlight_table <- function(muts_tab, protocol) {
  # Select mutations in tumorsuppressors and oncogenes as well as potentially deleterious
  highlight <- muts_tab[which(
    muts_tab$is_oncogene == 1 | muts_tab$is_tumorsuppressor == 1 |
    muts_tab$CLIN_SIG %in% c(
      "pathogenic",
      "likely_pathogenic",
      "conflicting_interpretations_of_pathogenicity"
    )
  ),
  c(
    "Hugo_Symbol",
    "HGVSp_Short",
    "t_AF",
    "is_tumorsuppressor",
    "is_oncogene",
    "is_hotspot",
    "CLIN_SIG",
    "REVEL",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Allele"
  )]
  colnames(highlight) <- c(
    "Symbol",
    "HGVSp_Short",
    "VAF",
    "TSG",
    "OG",
    "HS",
    "CLIN_SIG",
    "REVEL",
    "Chr",
    "Start",
    "Ref",
    "Alt"
  )
   if (dim(highlight)[1] != 0) {
    ### Interactive links ###
    # Genome Nexus
    highlight$Hugo_Symbol_new <- highlight$Symbol
    if (length(which(highlight$Allele == "-" | highlight$Reference_Allele == "-")) > 0) {
      highlight$Hugo_Symbol_new[- which(
        highlight$Allele == "-" | highlight$Reference_Allele == "-"
      )] <- l_gen_nex(df = highlight[-which(
        highlight$Allele == "-" | highlight$Reference_Allele == "-"
      ), ], type = "SNV") 
      if (length(which(highlight$Alt == "-")) > 0){
        highlight$Hugo_Symbol_new[which(
          highlight$Allele == "-"
        )] <- l_gen_nex(df = highlight[which(
          highlight$Allele == "-"
        ), ], type = "DEL")
      }
      if (length(which(highlight$Ref == "-")) > 0) {
      highlight$Hugo_Symbol_new[which(
        highlight$Reference_Allele == "-"
      )] <- l_gen_nex(df = highlight[which(
        highlight$Reference_Allele == "-"
      ), ], type = "INS")
      }
    } else {
      highlight$Hugo_Symbol_new <- l_gen_nex(df = highlight, type = "SNV")

    }
    # Cancer Consortium Meta-Knowledgebase
    highlight$HGVSp_Short <- meta(df = highlight)

    # VarSome links
    highlight$Varsome <- varsome(df = highlight, mode = "1")

    # VAF
    highlight$VAF <- gsub(pattern = "%", replacement = "", x = highlight$VAF)

    # ClinVar (ACMG)
    highlight$Classification <- acmg(df = highlight)

    # REVEL
    res_revel <- revel(df = highlight)
    highlight$REVEL <- res_revel$revel
    highlight$REVEL_cat <- res_revel$cat

    # Order table
    highlight <- highlight[order(
      as.numeric(highlight$VAF), decreasing = TRUE
    ), , drop = FALSE]
    highlight <- highlight[order(
      as.numeric(highlight$Classification), decreasing = TRUE
    ), , drop = FALSE]
    id_hs <- which(highlight$HS != 0)

    # combined ACMG classification
    highlight$Cat <- highlight$Classification
    highlight$Cat[highlight$Cat == "."] <- "0"
    highlight$Cat[duplicated(highlight$Classification)] <- "."
    highlight$Cat <- gsub(
      pattern = 5, replacement = "Pathogen", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 4, replacement = "Wahrsch. Pathogen", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 3, replacement = "VUS", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 2, replacement = "Wahrsch. Gutartig", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 1, replacement = "Gutartig", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 0, replacement = "Nicht Klassifiziert", x = highlight$Cat
    )

    # cancer genes
    highlight$Cancergene <- "."
    highlight$Cancergene[which(highlight$TSG == 1)] <- "TSG"
    highlight$Cancergene[which(highlight$OG == 1)] <- paste0(
      "OG",
      highlight$Cancergene[which(highlight$OG == 1)]
    )
    highlight$Cancergene <- gsub(
      pattern = "OGTSG",
      replacement = "both",
      x = highlight$Cancergene
    )

    # output
    highlight <- highlight[, c(
      "Cat",
      "Classification",
      "REVEL_cat",
      "Hugo_Symbol_new",
      "HGVSp_Short",
      "VAF",
      "Cancergene",
      "Varsome"
    )]
    colnames(highlight) <- c(
      "Kategorie",
      "ClinVar",
      "REVEL",
      "Gen",
      "AA-Austausch",
      "VAF [\\%]",
      "Cancergene",
      "VarSome"
    )

  } else {
    highlight <- data.frame()
    id_hs <- NULL
  }
    return(list(highlight = highlight, id_hs = id_hs))
}

highlight_detail <- function(muts_tab, Mode = "Tumor", protocol) {
 highlight <- muts_tab
  if (dim(muts_tab)[1] == 0) {
    muts_tab <- NULL
    id_hs <- NULL
  } else {
    # Genome Nexus
    highlight$Hugo_Symbol <- unlist(
      lapply(
        strsplit(
          x = as.character(highlight$Hugo_Symbol), split = ";", fixed = TRUE
        ), function(x) {
          return(x[[1]])
        }
      )
    )
    highlight$Hugo_Symbol_new <- highlight$Hugo_Symbol
    if (length(which(highlight$Allele == "-" | highlight$Reference_Allele == "-")) > 0) {
      highlight$Hugo_Symbol_new[-which(
        highlight$Allele == "-" | highlight$Reference_Allele == "-"
      )] <- l_gen_nex(df = highlight[-which(
        highlight$Allele == "-" | highlight$Reference_Allele == "-"
      ), ], type = "SNV") 
      if (length(which(highlight$Allele == "-")) > 0) {
        highlight$Hugo_Symbol_new[which(
          highlight$Allele == "-"
        )] <- l_gen_nex(df = highlight[which(
          highlight$Allele == "-"
        ), ], type = "DEL") 
      }
      if (length(which(highlight$Reference_Allele == "-")) > 0) {
        highlight$Hugo_Symbol_new[which(
          highlight$Reference_Allele == "-"
        )] <- l_gen_nex(df = highlight[which(
          highlight$Reference_Allele == "-"
        ), ], type = "INS")
      }
    } else {
      highlight$Hugo_Symbol_new <- l_gen_nex(df = highlight, type = "SNV")
    }
    # Cancer Consortium Meta-Knowledgebase
    highlight$HGVSp_Short <- meta(df = highlight)

    # Function
    highlight$Variant_Classification <- ex_func(df = highlight)

    # VarSome links
    highlight$Varsome <- varsome(df = highlight, mode = "2")

    # VAF
    if (Mode %in% c("Tumor", "Germline")) {
      highlight$VAF <- gsub(
        pattern = "%", replacement = "",
        x = highlight$t_AF
      )
      highlight$VAF <- paste0(
        highlight$VAF, " (", highlight$Variant_Reads, ")"
      )
    } else if (Mode == "LoH") {
      highlight$n_AF <- gsub(
        pattern = "%", replacement = "",
        x = highlight$n_AF, fixed = TRUE
      )
      highlight$t_AF <- gsub(
        pattern = "%",replacement = "",
        x = highlight$t_AF, fixed = TRUE
      )
      highlight$n_AF <- paste0(
        highlight$n_AF, " (", highlight$Count_Normal, ")"
      )
      highlight$t_AF <- paste0(
        highlight$t_AF, " (", highlight$Count_Tumor, ")"
      )
    }
    # ClinVar (ACMG)
    highlight$Classification <- acmg(df = highlight)

    # REVEL
    res_revel <- revel(df = highlight)
    highlight$REVEL <- res_revel$revel
    highlight$REVEL_cat <- res_revel$cat

    # COSMIC
    highlight$cosmic <- cosmic(df = highlight)

    # Order table
    if (Mode %in% c("Tumor", "Germline")) {
      highlight <- highlight[order(
        as.numeric(gsub(
          pattern = "%", replacement = "",
          x = highlight$t_AF
        )),
        decreasing = TRUE
      ), , drop = FALSE]
    } else if (Mode == "LoH") {
      highlight <- highlight[order(
        as.numeric(highlight$t_AF), decreasing = TRUE
      ), , drop = FALSE]
    }
    highlight <- highlight[order(
      as.numeric(highlight$Classification), decreasing = TRUE
    ), , drop = FALSE]
    id_hs <- which(highlight$is_hotspot != 0)

    # cancer genes
    highlight$Cancergene <- "."
    highlight$Cancergene[which(highlight$is_tumorsuppressor == 1)] <- "TSG"
    highlight$Cancergene[which(highlight$is_oncogene == 1)] <- "OG"
    highlight$Cancergene[which(
      highlight$is_oncogene == 1 & highlight$is_tumorsuppressor == 1
    )] <- "both"

    # Population frequency
    highlight$MAX_AF <- as.character(format(
      as.numeric(highlight$MAX_AF), scientific = TRUE, digits = 2
    ))
    highlight$MAX_AF[is.na(highlight$MAX_AF)] <- "."

    # output
    if (Mode == "Tumor") {
      muts_tab <- highlight[, c(
        "Hugo_Symbol_new",
        "HGVSp_Short",
        "Varsome",
        "VAF",
        "MAX_AF",
        "Classification",
        "REVEL_cat",
        "cosmic",
        "Cancergene"
      )]
      colnames(muts_tab) <- c(
        "Gen",
        "AA-Austausch",
        "Funktion / VarSome",
        "VAF [\\%] (Coverage)",
        "MAF",
        "ClinVar",
        "REVEL",
        "Cosmic",
        "Cancergene"
      )
    } else if (Mode == "LoH") {
      muts_tab <- highlight[, c(
        "Hugo_Symbol_new",
        "HGVSp_Short",
        "Varsome",
        "t_AF",
        "n_AF",
        "MAX_AF",
        "Classification",
        "REVEL_cat",
        "cosmic",
        "Cancergene"
      )]
      colnames(muts_tab) <- c(
        "Gen",
        "AA-Austausch",
        "Funktion / VarSome",
        "VAF [\\%] (Coverage) Tumor",
        "VAF [\\%] (Coverage) Keimbahn",
        "MAF",
        "ClinVar",
        "REVEL",
        "Cosmic",
        "Cancergene"
      )
    } else if (Mode == "Germline") {
      muts_tab <- highlight[, c(
        "Hugo_Symbol_new",
        "HGVSp_Short",
        "Varsome",
        "VAF",
        "MAX_AF",
        "Classification",
        "REVEL_cat",
        "cosmic",
        "Cancergene"
      )]
      colnames(muts_tab) <- c(
        "Gen",
        "AA-Austausch",
        "Funktion / VarSome",
        "VAF [\\%] (Coverage)",
        "MAF",
        "ClinVar",
        "REVEL",
        "Cosmic",
        "Cancergene"
      )
    }
  }
  return(list(muts_tab = muts_tab, id_hs = id_hs))
 }
 
summary_quality <- function(stats, protocol) {
  if (protocol == "somaticGermline" | protocol == "somatic") {
    q_t <- c()
    if (round(x = sum(
      stats$cover$cov[[2]][, 2] * stats$cover$cov[[2]][, 5]
    ), digits = 2) < 80) {
      q_t[1] <- "Akzeptabel" } else {
        q_t[1] <- "Sehr gut"
      }
    if (round(stats$cover$perc[[2]][1], digits = 2) * 100 > 80) {
      q_t[2] <- "Sehr gut"
    } else { q_t[2] <- "Akzeptabel" }
    if (round(stats$cover$perc[[2]][2], digits = 2) * 100 > 80) {
      q_t[3] <- "Sehr gut"
    } else {q_t[3] <- "Akzeptabel" }
    if (as.numeric(stats$avreads$tin) > 100) {
      q_t[4] <- "Sehr gut"
    } else {q_t[4] <- "Akzeptabel" }
    if (as.numeric(
      stats$qc_check$gc_content[[2]]
    ) >= 40 && as.numeric(stats$qc_check$gc_content[[2]]) <= 60) {
      q_t[5] <- "Sehr gut"
    } else {q_t[5] <- "Akzeptabel"}
    if (round(stats$qc_check$mean_QC[[2]], digits = 2) > 28) {
      q_t[6] <- "Sehr gut"
    } else {q_t[6] <- "Akzeptabel" }
    qp_t <- c(
      round(
        x = sum(
          stats$cover$cov[[2]][, 2] * stats$cover$cov[[2]][, 5]
        ), digits = 2
      ), paste0(round(stats$cover$perc[[2]][1], digits = 2) * 100, "%"),
      paste0(round(stats$cover$perc[[2]][2], digits = 2) * 100, "%"),
      as.numeric(stats$avreads$tin),
      as.numeric(stats$qc_check$gc_content[[2]]),
      round(stats$qc_check$mean_QC[[2]], digits = 2)
    )
    names(qp_t) <- c(
      "Mittlere Coverage",
      "Coverage > 8",
      "Coverage > 40",
      "Insertl\"ange",
      "GC-Anteil",
      "Mittlere Qualit\"at"
    )
    if (length(which(q_t != "Sehr gut") > 0)) {
      warn_t <- paste0(
        names(qp_t)[which(q_t != "Sehr gut")],
        ": ", qp_t[which(q_t != "Sehr gut")]
      )
      q_t1 <- paste0(warn_t, collapse = ", ")
    } else {
      q_t1 <- "Keine."
    }
    q_n <- c()
    if (round(
      x = sum(stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]),
      digits = 2
    ) < 80) {
      q_n[1] <- "Akzeptabel"} else {
        q_n[1] <- "Sehr gut" 
      }
    if (round(stats$cover$perc[[1]][1], digits = 2) * 100 > 80) {
      q_n[2] <- "Sehr gut"
    } else { q_n[2] <- "Akzeptabel"
    }
    if (round(stats$cover$perc[[1]][2], digits = 2) * 100 > 80) {
      q_n[3] <- "Sehr gut"
    } else { q_n[3] <- "Akzeptabel"
    }
    if (as.numeric(stats$avreads$gin) > 100) {
      q_n[4] <- "Sehr gut"
    } else { q_n[4] <- "Akzeptabel"
    }
    if (as.numeric(
      stats$qc_check$gc_content[[1]]
    ) >= 40 && as.numeric(
      stats$qc_check$gc_content[[1]]
    ) <= 60) {
      q_n[5] <- "Sehr gut"
    } else {
      q_n[5] <- "Akzeptabel"
    }
    if (round(stats$qc_check$mean_QC[[1]], digits = 2) > 28) {
      q_n[6] <- "Sehr gut"
    } else {q_n[6] <- "Akzeptabel"
    }
    qp_n <- c(
      round(
        x = sum(stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]),
        digits = 2
      ),
      paste0(round(stats$cover$perc[[1]][1], digits = 2) * 100 , "%"),
      paste0(round(stats$cover$perc[[1]][2], digits = 2) * 100, "%"),
      as.numeric(stats$avreads$gin),
      as.numeric(stats$qc_check$gc_content[[1]]),
      round(stats$qc_check$mean_QC[[1]], digits = 2)
    )
    names(qp_n) <- c(
      "Mittlere Coverage",
      "Coverage > 8",
      "Coverage > 40",
      "Insertl\"ange",
      "GC-Anteil",
      "Mittlere Qualit\"at"
    )
    if (length(which(q_n != "Sehr gut") > 0)) {
      warn_n <- paste0(
        names(qp_n)[which(q_n != "Sehr gut")],
        ": ",
        qp_n[which(q_n != "Sehr gut")]
      )
      q_n1 <- paste0(warn_n, collapse = ", ")
    } else {
      q_n1 <- "Keine."
    }
    tab <- rbind(c("Tumor", q_t1), c("Keimbahn", q_n1))
    colnames(tab) <- c("Probe" , "Auff\"alligkeiten")
    return(tab)
  }
  if (protocol == "panelTumor") {
    q_t <- c()
    if (round(x = sum(
      stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]
    ), digits = 2) < 150) {
      q_t[1] <- "Akzeptabel" } else {
        q_t[1] <- "Sehr gut"
      }
    if (round(stats$cover$perc[[1]][1], digits = 2) * 100 > 90) {
      q_t[2] <- "Sehr gut"
    } else { q_t[2] <- "Akzeptabel" }
    if (round(stats$cover$perc[[1]][2], digits = 2) * 100 > 90) {
      q_t[3] <- "Sehr gut"
    } else {q_t[3] <- "Akzeptabel" }
    if (as.numeric(stats$avreads$tin) > 100) {
      q_t[4] <- "Sehr gut"
    } else {q_t[4] <- "Akzeptabel" }
    if (as.numeric(
      stats$qc_check$gc_content[[1]]
    ) >= 40 && as.numeric(stats$qc_check$gc_content[[1]]) <= 60) {
      q_t[5] <- "Sehr gut"
    } else {q_t[5] <- "Akzeptabel"}
    if (round(stats$qc_check$mean_QC[[1]], digits = 2) > 28) {
      q_t[6] <- "Sehr gut"
    } else {q_t[6] <- "Akzeptabel" }
    qp_t <- c(
      round(
        x = sum(
          stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]
        ), digits = 2
      ), paste0(round(stats$cover$perc[[1]][1], digits = 2) * 100, "%"),
      paste0(round(stats$cover$perc[[1]][2], digits = 2) * 100, "%"),
      as.numeric(stats$avreads$tin),
      as.numeric(stats$qc_check$gc_content[[1]]),
      round(stats$qc_check$mean_QC[[1]], digits = 2)
    )
    names(qp_t) <- c(
      "Mittlere Coverage",
      "Coverage > 50",
      "Coverage > 150",
      "Insertl\"ange",
      "GC-Anteil",
      "Mittlere Qualit\"at"
    )
    if (length(which(q_t != "Sehr gut") > 0)) {
      warn_t <- paste0(
        names(qp_t)[which(q_t != "Sehr gut")],
        ": ", qp_t[which(q_t != "Sehr gut")]
      )
      q_t1 <- paste0(warn_t, collapse = ", ")
    } else {
      q_t1 <- "Keine."
    }
    tab <- rbind(c("Tumor", q_t1))
    colnames(tab) <- c("Probe" , "Auff\"alligkeiten")
    return(tab)
  }
  if (protocol == "tumorOnly") {
    q_t <- c()
    if (round(x = sum(
      stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]
    ), digits = 2) < 80) {
      q_t[1] <- "Akzeptabel" } else {
        q_t[1] <- "Sehr gut"
      }
    if (round(stats$cover$perc[[1]][1], digits = 2) * 100 > 80) {
      q_t[2] <- "Sehr gut"
    } else { q_t[2] <- "Akzeptabel" }
    if (round(stats$cover$perc[[1]][2], digits = 2) * 100 > 80) {
      q_t[3] <- "Sehr gut"
    } else {q_t[3] <- "Akzeptabel" }
    if (as.numeric(stats$avreads$tin) > 100) {
      q_t[4] <- "Sehr gut"
    } else {q_t[4] <- "Akzeptabel" }
    if (as.numeric(
      stats$qc_check$gc_content[[1]]
    ) >= 40 && as.numeric(stats$qc_check$gc_content[[1]]) <= 60) {
      q_t[5] <- "Sehr gut"
    } else {q_t[5] <- "Akzeptabel"}
    if (round(stats$qc_check$mean_QC[[1]], digits = 2) > 28) {
      q_t[6] <- "Sehr gut"
    } else {q_t[6] <- "Akzeptabel" }
    qp_t <- c(
      round(
        x = sum(
          stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]
        ), digits = 2
      ), paste0(round(stats$cover$perc[[1]][1], digits = 2) * 100, "%"),
      paste0(round(stats$cover$perc[[1]][2], digits = 2) * 100, "%"),
      as.numeric(stats$avreads$tin),
      as.numeric(stats$qc_check$gc_content[[1]]),
      round(stats$qc_check$mean_QC[[1]], digits = 2)
    )
    names(qp_t) <- c(
      "Mittlere Coverage",
      "Coverage > 8",
      "Coverage > 40",
      "Insertl\"ange",
      "GC-Anteil",
      "Mittlere Qualit\"at"
    )
    if (length(which(q_t != "Sehr gut") > 0)) {
      warn_t <- paste0(
        names(qp_t)[which(q_t != "Sehr gut")],
        ": ", qp_t[which(q_t != "Sehr gut")]
      )
      q_t1 <- paste0(warn_t, collapse = ", ")
    } else {
      q_t1 <- "Keine."
    }
    tab <- rbind(c("Tumor", q_t1))
    colnames(tab) <- c("Probe" , "Auff\"alligkeiten")
    return(tab)
  }
}

sum_muts <- function(tmp) {
  colnames(tmp) <- c(
    "Mutationstyp",
    "Anzahl",
    "Zygosit\"at",
    "Tumorsuppressoren",
    "Onkogene",
    "Hotspots"
  )
  tmp[c(1, 4), 3] <- "homozygot"
  tmp[c(2, 5), 3] <- "heterozygot"

  return(tmp)
}

cnv_cg <- function(gene_loci, type = "OG") {
  idl <- which(as.numeric(gene_loci$cn) < 1)
  idg <- which(as.numeric(gene_loci$cn) > 4)
  id_d <- which(is.na(as.numeric(gene_loci$cn)))
  if (type == "TSG") {
    colnames_df <- c("TSG", "CN", "CN-Typ", "Status")
  } else {
    colnames_df <- c("OG", "CN", "CN-Typ", "Status")
  }

  if(length(idl) > 0) {
    tmpl <- gene_loci[idl, c(4:6)]
    tmpl$Status <- "Loss"
    tmpl <- tmpl[order(as.numeric(tmpl$cn), decreasing = FALSE), ]
  } else {
    tmpl <- data.frame("", "", "", "")
  }
  colnames(tmpl) <- colnames_df

  if(length(idg) > 0) {
    tmpg <- gene_loci[idg, c(4:6)]
    tmpg$Status <- "Gain"
    tmpg <- tmpg[order(as.numeric(tmpg$cn), decreasing = TRUE), ]
  } else {
    tmpg <- data.frame("", "", "", "")
  }
  colnames(tmpg) <- colnames_df

  if (length(id_d) > 0){
    tmpd <- gene_loci[id_d, c(4:6)]
  } else {
    tmpd <- data.frame("", "", "")
  }
  colnames(tmpd) <- colnames_df[1:3]

  if((tmpg[1,1] != "") & (tmpl[1,1] != "")){
    id_l <- which(tmpl$`CN-Typ` < 1)
    id_g <- which(tmpg$`CN-Typ` > 7)
  }

  return(list(Gains = tmpg, Losses = tmpl, Different = tmpd))
}

cnv_panel <- function(cnv_results) {
  tmp <- cnv_analysis_results$out[, c(
    "Gene",
    "CopyNumber",
    "Status",
    "Type",
    "Cancergene"
  )]
  colnames(tmp) <- c(
    "Gen",
    "Kopien",
    "Status",
    "CN-Typ",
    "Cancergene"
  )
  return(tmp)
}

pathws_cnv <- function(df) {
  if (sum(unlist(lapply(df, function(x){return(dim(x)[1])}))) == 0) {
    return(NULL)
  } else {
    df <- lapply(
      1:length(df),
      function(id) {
        if(nrow(df[[id]])) {
          cbind(df[[id]], names(df)[id])
        }
      }
    )
    table_new <- lapply(df, function(x){
      if (!is.null(x)) {
        new_mat <- as.data.frame(
          matrix(
            nrow = length(unique(x[, 2])),
            ncol = 4,
            data = 0
          )
        )
        colnames(new_mat) <- c("Signalweg", "Status", "Kopien", "Gene")
        new_mat[, 1] <- x[1, 4]
        new_mat[, 3] <- unique(x[, 2])[
          order(unique(x[, 2]), decreasing = FALSE)
        ]
        for (i in 1:dim(new_mat)[1]) {
          new_mat[i, 4] <- paste(
            x[which(as.numeric(x[, 2]) == new_mat[i, 3]), 1],
            collapse = ", "
          )
          new_mat[i, 2] <- ifelse(
            test = as.numeric(new_mat[i, 3]) < 2,
            yes = "Loss",
            no = "Gain"
          )
        }
        return(new_mat)
      }
    })
    df <- data.frame(
      rbind(
        table_new[[1]],
        table_new[[2]],
        table_new[[3]],
        table_new[[4]],
        table_new[[5]]
      )
    )
    df[which(duplicated(df$Signalweg)), 1] <- "."
    df[which(df[, ] == "ddr"), 1] <- "DNA DamResp"
    df[which(df[, ] == "pam"), 1] <- "PI3K-AKT-mTOR"
    df[which(df[, ] == "rme"), 1] <- "RAF-MEK-ERK"
    df[which(df[, ] == "tyk"), 1] <- "Tyr Kin"
    df[which(df[, ] == "cec"), 1] <- "Cell Cycle"
    return(df)
  }
}

pthws_mut <- function(df, protocol) {
  id_to <- which(df$Pathway == "Topart")
  if (length(id_to) != 0) {
    df <- df[-c(id_to:dim(df)[1]), ]
  }
  if (is.null(df) || nrow(df) == 0 || class(df) == "list" || is.null(df)) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    } else if (sum(lapply(df, function(x){return(dim(x)[1])})) == 0){
      return(NULL)
    }
  } else if (class(df) == "data.frame") {
    df <- df[, c(
      "Pathway",
      "Symbol",
      "HGVSp_Short",
      "Variant_Classification",
      "VAF",
      "Reads",
      "MAF",
      "CLIN_SIG",
      "REVEL_score",
      "Existing_variation"
    )]
    colnames(df) <- c(
      "Pathway",
      "Symbol",
      "HGVSp_Short",
      "Variant_Classification",
      "VAF",
      "Reads",
      "MAF",
      "CLIN_SIG",
      "REVEL",
      "Existing_variation"
    )
    # Pathway
    df$Pathway <- gsub(
      pattern = "DNA Damage Response",
      replacement = "DNA DamResp",
      x = df$Pathway
    )
    df$Pathway <- gsub(
      pattern = "Tyrosine Kinases",
      replacement = "Tyr Kin",
      x = df$Pathway
    )
    # HGVSp_Short
    df$HGVSp_Short <- substr(
      x = df$HGVSp_Short, start = 3, stop = nchar(df$HGVSp_Short)
    )
    # VAF
    df$VAF <- gsub(
      pattern = "%", replacement = "", x = df$VAF, fixed = TRUE
    )
    df$VAF <- paste0(df$VAF, " (", df$Reads, ")")
    # HGVSp_Short
    df$HGVSp_Short <- gsub(
      pattern = "_", replacement = "\\_", x = df$HGVSp_Short, fixed = TRUE
    )
    # COSMIC
    df$cosmic <- cosmic(df)
    # Function
    df$Variant_Classification <- ex_func(df)
    # REVEL
    res_revel <- revel(df)
    df$REVEL <- res_revel$revel
    df$REVEL_cat <- res_revel$cat
    # ClinVar (ACMG)
    df$Classification <- acmg(df = df)
    # Population frequency
    df$MAF <- as.character(
      format(as.numeric(df$MAF), scientific = TRUE, digits = 2)
    )
    # output
    df <- df[, c("Pathway", "Symbol", "HGVSp_Short", "Variant_Classification", "VAF", "MAF", "Classification", "REVEL_cat", "cosmic")]
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    colnames(df) <- c(
      "Signalweg",
      "Gen",
      "AA-Austausch",
      "Funktion",
      "VAF [\\%] (Coverage)",
      "MAF",
      "ClinVar",
      "REVEL",
      "Cosmic"
    )
    return(df)
  }
}

topart_mut <- function(df, protocol) {
  id_to <- which(df$Pathway == "Topart")
  if (length(id_to) == 0) {
    return(NULL)
  } else {
    df <- df[id_to:dim(df)[1], c(
      "Pathway",
      "Symbol",
      "HGVSp_Short",
      "Variant_Classification",
      "VAF",
      "Reads",
      "MAF",
      "CLIN_SIG",
      "REVEL_score",
      "Existing_variation"
    )]
    colnames(df) <- c(
      "Pathway",
      "Symbol",
      "HGVSp_Short",
      "Variant_Classification",
      "VAF",
      "Reads",
      "MAF",
      "CLIN_SIG",
      "REVEL",
      "Existing_variation"
    )
    # HGVSp_Short
    df$HGVSp_Short <- substr(
      x = df$HGVSp_Short, start = 3, stop = nchar(df$HGVSp_Short)
    )
    # VAF
    df$VAF <- gsub(pattern = "%", replacement = "", x = df$VAF, fixed = TRUE)

    df$VAF <- paste0(df$VAF, " (", df$Reads, ")")
    # HGVSp_Short
    df$HGVSp_Short <- gsub(
      pattern = "_", replacement = "\\_", x = df$HGVSp_Short, fixed = TRUE
    )
    # COSMIC
    df$cosmic <- cosmic(df)
    # Function
    df$Variant_Classification <- ex_func(df)
    # REVEL
    res_revel <- revel(df)
    df$REVEL <- res_revel$revel
    df$REVEL_cat <- res_revel$cat
    # ClinVar (ACMG)
    df$Classification <- acmg(df = df)
    # Population frequency
    df$MAF <- as.character(format(
      as.numeric(df$MAF), scientific = TRUE, digits = 2
    ))

    # output
    print(colnames(df))
    df <- df[, c("Pathway", "Symbol", "HGVSp_Short", "Variant_Classification", "VAF", "MAF", "Classification", "REVEL_cat", "cosmic")]
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    colnames(df) <- c(
      "Signalweg",
      "Gen",
      "AA-Austausch",
      "Funktion",
      "VAF [\\%] (Coverage)",
      "MAF",
      "ClinVar",
      "REVEL",
      "Cosmic"
    )
    return(df)
  }
}

med_tmb <- function(entity) {
  if (entity == "ACC") {
      med <- 0.82
      sd <- c(0.55, 1.86)
  } else if (entity == "BLCA") {
      med <- 5.26
      sd <- c(2.79, 9.36)
  } else if (entity == "BRCA") {
    med <- 2.63
    sd <- c(0.20, 290.80)
  } else if (entity == "CESC") {
      med <- 3.32
      sd <- c(2.00, 6.07)
  } else if (entity == "CHOL") {
      med <- 1.71
      sd <- c(1.25, 2.80)
  } else if (entity == "COAD") {
      med <- 3.68
      sd <- c(2.67, 5.71)
  } else if (entity == "ESCA") {
      med <- 4.18
      sd <- c(3.16, 6.05)
  } else if (entity == "HNSC") {
      med <- 3.34
      sd <- c(2.05, 5.51)
  } else if (entity == "KIRC") {
      med <- 1.63
      sd <- c(1.21, 2.13)
  } else if (entity == "KIRP") {
      med <- 2.13
      sd <- c(1.34, 2.89)
  } else if (entity == "LIHC") {
      med <- 2.89
      sd <- c(2.03, 3.95)
  } else if (entity == "LUAD") {
      med <- 6.17
      sd <- c(2.47, 12.66)
  } else if (entity == "LUSC") {
      med <- 7.13
      sd <- c(5.01, 10.47)
  } else if (entity == "MESO") {
      med <- 0.95
      sd <- c(0.66, 1.16)
  } else if (entity == "OV") {
      med <- 2.13
      sd <- c(1.50, 3.16)
  } else if (entity == "PAAD") {
      med <- 1.09
      sd <- c(0.76, 1.47)
  } else if (entity == "SKCM") {
      med <- 13.09
      sd <- c(5.91, 28.09)
  } else if (entity == "STAD") {
      med <- 3.5
      sd <- c(2.16, 7.45)
  } else if (entity == "UCEC") {
      med <- 2.63
      sd <- c(1.58, 19.47)
  } else if (entity == "UCS") {
      med <- 1.37
      sd <- c(1.21, 1.79)
  } else if (entity == "UVM") {
      med <- 0.34
      sd <- c(0.26, 0.45)
  } else {
      sd <- NULL
      if (entity == "STA") {
        med <- 3.3
      } else if (entity == "US") {
        med <- 2.6
      } else if (entity == "STR") {
        med <- 2.5
      } else if (entity %in% c(
        "STL", "STU", "BOS", "STS", "STRE",
        "STMPNST", "UPG", "ULS"
      )) {
        med <- 2.5
      } else if (entity == "STMS") {
        med <- 2.2
      } else if (entity == "SIG") {
        med <- 1.8
      } else if (entity == "SG") {
        med <- 1.8
      } else if (entity %in% c(
        "STC", "SSTSFT", "STFS", "UESS",
        "BCS", "STES", "STLS", "STDSRCT",
        "STSS", "STRSA"
      )) {
        med <- 1.7
      } else if (entity == "BC") {
        med <- 1.3
      } else if (entity == "STF") {
        med <- 0.9
      } else if (entity == "MRT") {
        med <- 0.2
      } else {
        med <- NULL
      }
  }
  return(list(med = med, sd = sd))
}
