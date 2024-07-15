#!/usr/bin/env bash

# common settings
readonly CFG_AUTHOR=$(get_config_value common.author "${PARAM_DIR_PATIENT}")
readonly CFG_CENTER=$(get_config_value common.center "${PARAM_DIR_PATIENT}")
readonly CFG_PROTOCOL=$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")
#readonly CFG_SEX=$(get_config_value sex "${PARAM_DIR_PATIENT}")
readonly CFG_ENTITY=$(get_config_value common.entity "${PARAM_DIR_PATIENT}")

# temporary folder
readonly DIR_TMP="$(get_config_value common.dirTmp "${PARAM_DIR_PATIENT}")/${PARAM_DIR_PATIENT}"

readonly CFG_FILE_TUMOR_R1=$(get_config_value common.files.tumor_R1 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_TUMOR_R2=$(get_config_value common.files.tumor_R2 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE_R1=$(get_config_value common.files.germline_R1 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE_R2=$(get_config_value common.files.germline_R2 "${PARAM_DIR_PATIENT}")

readonly CFG_PANEL_FILE_TUMOR=$(get_config_value common.files.panel.tumor "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_FILE_NUMBER=$(get_config_value common.files.panel.numberOfFiles "${PARAM_DIR_PATIENT}")

readonly CFG_FOLDER_RNA=$(get_config_value common.RNA.folder "${PARAM_DIR_PATIENT}")

# folder containing patient output
readonly DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
readonly DIR_WES="${DIR_TARGET}/WES"
readonly DIR_ANALYSES="${DIR_TARGET}/Analyses"
readonly DIR_RNA="${DIR_TARGET}/RNA"
readonly DIR_FUSIONS="${DIR_RNA}/fusioncatcher"

# end paths

# ucsc mysql server
readonly CFG_UCSC_SERVER=$(get_config_value common.ucscServer "${PARAM_DIR_PATIENT}")

# CNV annotation
readonly CFG_CNV_ANNOTATION=$(get_config_value common.cnvAnnotation "${PARAM_DIR_PATIENT}")

## Genome
readonly FILE_GENOME="${DIR_REF}/genome/$(get_config_value reference.genome "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_LENGTH="${DIR_CHROMOSOMES}/$(get_config_value reference.length "${PARAM_DIR_PATIENT}")"
readonly HRD_REF_WIG="${DIR_DATABASE}/$(get_config_value reference.hrdRef "${PARAM_DIR_PATIENT}")"

# depending on measurement machine
## SureSelect (Capture Kit)
readonly CFG_REFERENCE_CAPTUREREGIONS="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.captureRegions "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTUREGENES="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.captureGenes "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_COVEREDREGION="$(get_config_value reference.sequencing.coveredRegion "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTUREREGIONNAME="$(get_config_value reference.sequencing.captureRegionName "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTURECORFACTORS="${DIR_DATABASE}/$(get_config_value reference.sequencing.captureCorFactors "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_COVERED_EXONS="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.coveredExons "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_ACTIONABLEGENES="${DIR_DATABASE}/$(get_config_value reference.sequencing.actionableGenes "${PARAM_DIR_PATIENT}")"

# database for known variants
## dbSNP vcf File
readonly CFG_REFERENCE_DBSNP="${DIR_DBSNP}/$(get_config_value reference.dbSNP "${PARAM_DIR_PATIENT}")"
# END variables

# if no parallel computation, set processes to 1
if [[ -z ${PARALLEL_PROCESSES} ]]; then
  readonly TMP_PROCESSES=1
else
  readonly TMP_PROCESSES="${PARALLEL_PROCESSES}"
fi

# take cpucores/2
readonly CFG_COMMON_CPUCORES=$(($(get_config_value common.cpucores "${PARAM_DIR_PATIENT}")/${TMP_PROCESSES}))

# take memory/2
readonly tmp_memory=$(get_config_value common.memory "${PARAM_DIR_PATIENT}")
readonly CFG_COMMON_MEMORY="$(("${tmp_memory//[^0-9.]/}"/${TMP_PROCESSES}))${tmp_memory//[^a-zA-Z]/}"

# general parameters
readonly CFG_GENERAL_MINBASEQUAL=$(get_config_value general.minBaseQual "${PARAM_DIR_PATIENT}")
readonly CFG_GENERAL_MAFCUTOFF=$(get_config_value general.maf_cutoff "${PARAM_DIR_PATIENT}")
readonly CFG_GENERAL_MINVAF=$(get_config_value general.minVAF "${PARAM_DIR_PATIENT}")
readonly CFG_GENERAL_MINGERMLINEVAF=$(get_config_value general.minGermlineVAF "${PARAM_DIR_PATIENT}")

# mutect2
readonly CFG_MUTECT_CALLABLEDEPTH=$(get_config_value wes.mutect.callableDepth "${PARAM_DIR_PATIENT}")
readonly CFG_MUTECT_PANELOFNORMALS="${DIR_REF}/genome/$(get_config_value wes.mutect.panelOfNormals "${PARAM_DIR_PATIENT}")"
readonly CFG_MUTECT_GERMLINERESOURCE="${DIR_REF}/genome/$(get_config_value wes.mutect.germlineResource "${PARAM_DIR_PATIENT}")"
readonly CFG_MUTECT_GETPILEUPSUMMARIES_KNOWNVARIANTSITES="${DIR_REF}/genome/$(get_config_value wes.mutect.getpileupsummaries.knownVariantSites "${PARAM_DIR_PATIENT}")"

# Panel Parameter
# mutect2
readonly CFG_PANEL_MUTECT_CALLABLEDEPTH=$(get_config_value panel.mutect.callableDepth "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_MUTECT_PANELOFNORMALS="${DIR_REF}/genome/$(get_config_value panel.mutect.panelOfNormals "${PARAM_DIR_PATIENT}")"
readonly CFG_PANEL_MUTECT_GETPILEUPSUMMARIES_KNOWNVARIANTSITES="${DIR_REF}/genome/$(get_config_value wes.mutect.getpileupsummaries.knownVariantSites "${PARAM_DIR_PATIENT}")"

# Fusions
readonly CFG_FUSION_GENES="${DIR_SEQUENCING}/$(get_config_value panel.fusions.fusionGenes "${PARAM_DIR_PATIENT}")"

# Amplifications
readonly CFG_AMPLIFICATION_GENES="${DIR_SEQUENCING}/$(get_config_value panel.amplification.amplificationGenes "${PARAM_DIR_PATIENT}")"

## tumorOnly Parameter
# mutect2
readonly CFG_TUMORONLY_MUTECT_CALLABLEDEPTH=$(get_config_value tumorOnly.mutect.callableDepth "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_MUTECT_PANELOFNORMALS="${DIR_REF}/genome/$(get_config_value tumorOnly.mutect.panelOfNormals "${PARAM_DIR_PATIENT}")"
readonly CFG_TUMORONLY_MUTECT_GERMLINERESOURCE="${DIR_REF}/genome/$(get_config_value tumorOnly.mutect.germlineResource "${PARAM_DIR_PATIENT}")"
readonly CFG_TUMORONLY_MUTECT_GETPILEUPSUMMARIES_KNOWNVARIANTSITES="${DIR_REF}/genome/$(get_config_value wes.mutect.getpileupsummaries.knownVariantSites "${PARAM_DIR_PATIENT}")"

## Tools and paths
# Paths
readonly BIN_JAVA="java -Djava.io.tmpdir=${DIR_TMP} " # path to java

# Pre-Processing
readonly BIN_FASTQC="${DIR_TOOLS}/FastQC/bin/fastqc -t ${CFG_COMMON_CPUCORES} --extract "

readonly BIN_TRIM="${BIN_JAVA} -jar ${DIR_TOOLS}/Trimmomatic/trimmomatic.jar PE -threads ${CFG_COMMON_CPUCORES} -phred33 "
readonly DIR_TRIMMOMATIC_ADAPTER="${DIR_TOOLS}/Trimmomatic/adapters"
readonly BIN_CUT="cut -f1,2,3"

# Alignment
readonly BIN_BWAMEM="${DIR_TOOLS}/bwa-mem2/bwa-mem2 mem -M "
readonly BIN_BWAMEMINDEX="${DIR_TOOLS}/bwa-mem2/bwa-mem2 index"

# SAMTOOLS
readonly BIN_SAMTOOLS="${DIR_TOOLS}/samtools/samtools" # path to samtools
readonly BIN_SAMVIEW="${BIN_SAMTOOLS} view -@ ${CFG_COMMON_CPUCORES} "
readonly BIN_SAMSORT="${BIN_SAMTOOLS} sort -@ ${CFG_COMMON_CPUCORES} "
readonly BIN_SAMRMDUP="${BIN_SAMTOOLS} rmdup "
readonly BIN_SAMINDEX="${BIN_SAMTOOLS} index "
readonly BIN_STATS="${BIN_SAMTOOLS} stats "

# GATK4
readonly BIN_GATK4="${DIR_TOOLS}/gatk4/gatk" # --java-options '-Xmx${CFG_COMMON_MEMORY}'"
readonly BIN_INDEL_REALIGNER="${BIN_GATK4} LeftAlignIndels -R ${FILE_GENOME} "
readonly BIN_BASE_RECALIBRATOR="${BIN_GATK4} BaseRecalibrator -R ${FILE_GENOME} --known-sites ${CFG_REFERENCE_DBSNP} "
readonly BIN_PRINT_READS="${BIN_GATK4} ApplyBQSR -R ${FILE_GENOME} "

# PICARD
readonly BIN_FIX_MATE="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -Dpicard.useLegacyParser=false -jar ${DIR_TOOLS}/picard/picard.jar FixMateInformation "

# VEP
readonly BIN_VEP="${DIR_TOOLS}/vep/vep"
readonly BIN_VCF2MAF="perl ${DIR_TOOLS}/vcf2maf/vcf2maf.pl "

# COVERAGE
readonly BIN_COVERAGE="${DIR_TOOLS}/bedtools2/bin/bedtools coverage -hist -g ${FILE_GENOME}.fai -sorted "

# FREEC
readonly BIN_FREEC="${DIR_TOOLS}/FREEC/bin/freec "
readonly FILE_REFERENCE_MAPPABILITY="${DIR_REF}/mappability/$(get_config_value reference.mappability "${PARAM_DIR_PATIENT}")"

# R
readonly BIN_RSCRIPT=$(command -v Rscript)

# fusioncatcher
readonly FUSIONCATCHER_DB="${DIR_DATABASE}/fusioncatcher/data/current"
readonly BIN_FUSIONCATCHER="/usr/bin/fusioncatcher -p ${CFG_COMMON_CPUCORES} -d ${FUSIONCATCHER_DB} "

# msisnesor2
readonly MSISENSOR2="${DIR_TOOLS}/msisensor2/msisensor2 msi -b ${CFG_COMMON_CPUCORES} -M ${DIR_TOOLS}/msisensor2/models_hg19_GRCh37"

# msisensor-pro
readonly MSISENSOR_PRO="${DIR_TOOLS}/msisensor-pro/binary/msisensor-pro msi -b ${CFG_COMMON_CPUCORES}"
readonly MSISENSOR_PRO_SCAN="${DIR_TOOLS}/msisensor-pro/binary/msisensor-pro scan"
readonly MICROSATELLITE_SITES="${DIR_DATABASE}/$(get_config_value reference.microsatelliteSites "${PARAM_DIR_PATIENT}")"

# HRD
readonly SEQUENZA_UTILS=$(command -v sequenza-utils)
readonly SEQUENZA_WINDOW=$(get_config_value sequenza.window "${PARAM_DIR_PATIENT}")
readonly SEQUENZA_NON_MATCHING_NORMAL="${DIR_REF}/sequenza/$(get_config_value sequenza.nonMatchingNormal "${PARAM_DIR_PATIENT}")"
readonly SEQUENZA_CHROMOSOMES=$(get_config_value sequenza.chromosomes "${PARAM_DIR_PATIENT}")

# BAM-Matcher
readonly BAM_MATCHER="python3 /opt/MIRACUM-Pipe/tools/bam-matcher/bam-matcher.py --reference ${FILE_GENOME} --short-output --do-not-cache"

# export parameters
export CFG_AUTHOR
export CFG_CENTER
export CFG_PROTOCOL
export CFG_ENTITY

export CFG_FILE_TUMOR_R1
export CFG_FILE_TUMOR_R2
export CFG_FILE_GERMLINE_R1
export CFG_FILE_GERMLINE_R2

export CFG_PANEL_FILE_TUMOR
export CFG_PANEL_FILE_NUMBER

export CFG_FOLDER_RNA

export CFG_UCSC_SERVER

export CFG_CNV_ANNOTATION

export DIR_TARGET
export DIR_WES
export DIR_ANALYSES
export DIR_RNA
export DIR_FUSIONS

export FILE_GENOME
export CFG_REFERENCE_LENGTH
export HRD_REF_WIG

export CFG_REFERENCE_CAPTUREREGIONS
export CFG_REFERENCE_CAPTUREGENES
export CFG_REFERENCE_COVEREDREGION
export CFG_REFERENCE_CAPTUREREGIONNAME
export CFG_REFERENCE_CAPTURECORFACTORS
export CFG_REFERENCE_COVERED_EXONS
export CFG_REFERENCE_ACTIONABLEGENES

export CFG_REFERENCE_DBSNP

export CFG_COMMON_CPUCORES
export CFG_COMMON_MEMORY

export CFG_GENERAL_MINBASEQUAL
export CFG_GENERAL_MAFCUTOFF
export CFG_GENERAL_MINVAF
export CFG_GENERAL_MINGERMLINEVAF

export CFG_PANEL_MUTECT_CALLABLEDEPTH
export CFG_PANEL_MUTECT_PANELOFNORMALS
export CFG_PANEL_MUTECT_GERMLINERESOURCE
export CFG_PANEL_MUTECT_GETPILEUPSUMMARIES_KNOWNVARIANTSITES

export CFG_TUMORONLY_MUTECT_CALLABLEDEPTH
export CFG_TUMORONLY_MUTECT_PANELOFNORMALS
export CFG_TUMORONLY_MUTECT_GERMLINERESOURCE
export CFG_TUMORONLY_MUTECT_GETPILEUPSUMMARIES_KNOWNVARIANTSITES

export CFG_MUTECT_CALLABLEDEPTH
export CFG_MUTECT_PANELOFNORMALS
export CFG_MUTECT_GERMLINERESOURCE
export CFG_MUTECT_GETPILEUPSUMMARIES_KNOWNVARIANTSITES

export CFG_FUSION_GENES
export CFG_AMPLIFICATION_GENES

export BIN_JAVA

export BIN_FASTQC
export BIN_TRIM
export DIR_TRIMMOMATIC_ADAPTER
export BIN_CUT

export BIN_BWAMEM
export BIN_BWAMEMINDEX

export BIN_SAMTOOLS
export BIN_SAMVIEW
export BIN_SAMSORT
export BIN_SAMRMDUP
export BIN_SAMINDEX
export BIN_MPILEUP
export BIN_STATS


export BIN_GATK4
export BIN_INDEL_REALIGNER
export BIN_BASE_RECALIBRATOR
export BIN_PRINT_READS

export BIN_FIX_MATE

export BIN_VAR_SCAN
export BIN_SOMATIC
export BIN_PROCESSSOMATIC

export BIN_VEP
export BIN_VCF2MAF

export BIN_COVERAGE

export BIN_FREEC
export FILE_REFERENCE_MAPPABILITY

export BIN_RSCRIPT

export FUSIONCATCHER_DB
export BIN_FUSIONCATCHER

export MSISENSOR2
export MSISENSOR_PRO
export MSISENSOR_PRO_SCAN
export MICROSATELLITE_SITES

export SEQUENZA_UTILS
export SEQUENZA_WINDOW
export SEQUENZA_NON_MATCHING_NORMAL
export SEQUENZA_CHROMOSOMES

export BAM_MATCHER