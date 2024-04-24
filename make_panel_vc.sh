#!/usr/bin/env bash

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)

## load settings
# shellcheck source=common.cfg.sh
. "${DIR_SCRIPT}"/common.cfg.sh


function usage() {
  echo "usage: miracum_pipe.sh -d dir [-h]"
  echo "  -d  dir             specify relative folder of patient"
  echo "  -p                  computing as parallel process"
  echo "  -h                  show this help screen"
  exit 1
}

while getopts d:t:ph option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
  h) usage ;;
  p) readonly PARALLEL_PROCESSES=2 ;;
  \?)
    echo "Unknown option: -$OPTARG" >&2
    exit 1
    ;;
  :)
    echo "Missing option argument for -$OPTARG" >&2
    exit 1
    ;;
  *)
    echo "Unimplemented option: -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# if no patient is defined
if [[ -z "${PARAM_DIR_PATIENT}" ]]; then
  echo "no patient defined."
  echo "--"
  usage
fi


# load patient yaml
readonly CFG_SEX=$(get_config_value sex "${PARAM_DIR_PATIENT}")
if [[ "$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")" = "panel" ]]; then
  readonly CFG_CASE=panelTumor
fi

# check inputs
readonly VALID_SEXES=("XX XY")

if [[ ! " ${VALID_SEXES[@]} " =~ " ${CFG_SEX} " ]]; then
  echo "unknown sex: ${CFG_SEX}"
  echo "use one of the following values: $(join_by ' ' ${VALID_SEXES})"
  exit 1
fi

##################################################################################################################

## load programs
# shellcheck source=programs.cfg.sh
. "${DIR_SCRIPT}/programs.cfg.sh"

##################################################################################################################

[[ -d "${DIR_ANALYSES}" ]] || mkdir -p "${DIR_ANALYSES}"

# names
readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_vc
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td

# temp
readonly mpileup=${DIR_TMP}/${NameD}_mpileup # DIR_WES

# keep
readonly recalbam=${DIR_WES}/${NameTD}_output.sort.rmdup.realigned.fixed.recal.bam
readonly OUTPUT_GZ=${DIR_WES}/${NameTD}_gatk4_mutect2.vcf.gz
readonly OUTPUT_FILTERED_GZ=${DIR_WES}/${NameTD}_gatk4_mutect2_filtered.vcf.gz
readonly OUTPUT=${DIR_WES}/${NameTD}_gatk4_mutect2_filtered
readonly MSI_OUTPUT=${DIR_WES}/${NameTD}_MSI

# Mutect2
${BIN_GATK4} Mutect2 -R ${FILE_GENOME} -I ${recalbam} -O ${OUTPUT_GZ} \
 --callable-depth "${CFG_PANEL_MUTECT_CALLABLEDEPTH}" --intervals "${CFG_REFERENCE_CAPTUREREGIONS}" --min-base-quality-score "${CFG_GENERAL_MINBASEQUAL}" --base-quality-score-threshold "${CFG_GENERAL_MINBASEQUAL}" --panel-of-normals "${CFG_PANEL_MUTECT_PANELOFNORMALS}" --genotype-pon-sites --germline-resource ${CFG_PANEL_MUTECT_GERMLINERESOURCE} --genotype-germline-sites 

# Filter
${BIN_GATK4} FilterMutectCalls -V ${OUTPUT_GZ} -R ${FILE_GENOME} -O ${OUTPUT_FILTERED_GZ} --intervals "${CFG_REFERENCE_CAPTUREREGIONS}" --min-median-base-quality "${CFG_GENERAL_MINBASEQUAL}" --min-allele-fraction "${CFG_GENERAL_MINVAF}"

${BIN_VEP} \
  --offline --cache \
  --assembly GRCh37 \
  --fork ${CFG_COMMON_CPUCORES} \
  --input_file  ${OUTPUT_FILTERED_GZ} \
  --output_file ${OUTPUT}_vep.vcf \
  --vcf \
  --species homo_sapiens \
  --force_overwrite \
  --symbol \
  --numbers \
  --regulatory \
  --canonical \
  --gene_phenotype \
  --af \
  --max_af \
  --af_gnomad \
  --pubmed \
  --allele_number \
  --hgvs \
  --hgvsg \
  --refseq \
  --plugin CADD,${DIR_DATABASE}/vep/CADD_GRCh37/whole_genome_SNVs.tsv.gz,${DIR_DATABASE}/vep/CADD_GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz \
  --plugin REVEL,${DIR_DATABASE}/vep/REVEL/new_tabbed_revel.tsv.gz \
  --no_stats --quiet

${BIN_VCF2MAF} \
  --inhibit-vep \
  --input-vcf ${OUTPUT}_vep.vcf \
  --output-maf ${OUTPUT}_vep.maf \
  --verbose \
  --ref-fasta ${FILE_GENOME} \
  --maf-center "${CFG_CENTER}" \
  --ncbi-build GRCh37 \
  --tumor-id ${NameTD} \
  --retain-fmt GT,AF \
  --retain-ann HGVSg,CADD_PHRED,MAX_AF,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_OTH,gnomADe_SAS,REVEL \
  --custom-enst ${DIR_DATABASE}/refseq_selected.txt

# MSI
${MSISENSOR2} -t "${recalbam}" -o "${MSI_OUTPUT}"

#eo VC
