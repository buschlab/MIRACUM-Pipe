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
if [[ "$(get_config_value common.germline "${PARAM_DIR_PATIENT}")" = "True" ]]; then
  readonly CFG_CASE=somaticGermline
else
  readonly CFG_CASE=somatic
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
readonly NameGD=${CFG_CASE}_${PARAM_DIR_PATIENT}_gd
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td

# input
readonly recalbamGD=${DIR_WES}/${NameGD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
readonly recalbamTD=${DIR_WES}/${NameTD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam

# keep
readonly TD_OUTPUT_GZ=${DIR_WES}/${NameTD}_gatk4_mutect2.vcf.gz
readonly GD_OUTPUT_GZ=${DIR_WES}/${NameTD}_gatk4_haplotype.vcf.gz
readonly TD_OUTPUT_FILTERED_GZ=${DIR_WES}/${NameTD}_gatk4_mutect2_filtered.vcf.gz
readonly TD_OUTPUT=${DIR_WES}/${NameTD}_gatk4_mutect2_filtered
readonly GD_OUTPUT=${DIR_WES}/${NameGD}_gatk4_haplotype
readonly MSI_OUTPUT=${DIR_WES}/${NameD}_MSI
readonly BAMMATCHER_OUTPUT=${DIR_WES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_bam-matcher.txt

# GATK4 Mutect2
# ANNOVAR settings
CODINGARG="--includesnp --onlyAltering --mrnaseq --tolerate"
CONVERTARG="--includeinfo"

# Mutect2
${BIN_GATK4} Mutect2 -R ${FILE_GENOME} -I ${recalbamTD} -I ${recalbamGD} --normal ${NameGD} -O ${TD_OUTPUT_GZ} \
 --callable-depth "${CFG_MUTECT_CALLABLEDEPTH}" --intervals "${CFG_REFERENCE_CAPTUREREGIONS}" --min-base-quality-score "${CFG_GENERAL_MINBASEQUAL}" --base-quality-score-threshold "${CFG_GENERAL_MINBASEQUAL}" --panel-of-normals "${CFG_MUTECT_PANELOFNORMALS}" --genotype-pon-sites --germline-resource ${CFG_MUTECT_GERMLINERESOURCE} --genotype-germline-sites 

# Filter
${BIN_GATK4} FilterMutectCalls -V ${TD_OUTPUT_GZ} -R ${FILE_GENOME} -O ${TD_OUTPUT_FILTERED_GZ} --intervals "${CFG_REFERENCE_CAPTUREREGIONS}" --min-median-base-quality "${CFG_GENERAL_MINBASEQUAL}" --min-allele-fraction "${CFG_GENERAL_MINVAF}"
gunzip "${TD_OUTPUT_FILTERED_GZ}"

# Annovar
${TABLEANNOVAR} "${TD_OUTPUT}.vcf" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" \
 --buildver hg19 --outfile "${TD_OUTPUT}" --operation "${CFG_ANNOVAR_ARGOP}" \
 --nastring . --vcfinput --thread "${CFG_COMMON_CPUCORES}" --maxgenethread "${CFG_COMMON_CPUCORES}" \
 --otherinfo --remove --verbose --polish \
 --convertarg "${CONVERTARG}"
rm "${TD_OUTPUT}.avinput"

# Haplotype Caller
if [[ $CFG_CASE = "somaticGermline" ]]; then
  ${BIN_GATK4} HaplotypeCaller -R ${FILE_GENOME} -I ${recalbamTD} -I ${recalbamGD} -O ${GD_OUTPUT_GZ} \
  --intervals "${CFG_REFERENCE_CAPTUREREGIONS}" --min-base-quality-score "${CFG_GENERAL_MINBASEQUAL}" --base-quality-score-threshold "${CFG_GENERAL_MINBASEQUAL}"
  gunzip "${GD_OUTPUT_GZ}"

  # Annovar
  ${TABLEANNOVAR} "${GD_OUTPUT}.vcf" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" \
  --buildver hg19 --outfile "${GD_OUTPUT}" --operation "${CFG_ANNOVAR_ARGOP}" \
  --nastring . --vcfinput --thread "${CFG_COMMON_CPUCORES}" --maxgenethread "${CFG_COMMON_CPUCORES}" \
  --otherinfo --remove --verbose --polish \
  --convertarg "${CONVERTARG}"
  rm "${GD_OUTPUT}.avinput"
fi

# MSI
if [ ! -f "${MICROSATELLITE_SITES}" ]; then
    echo "${MICROSATELLITE_SITES} does not exist. Generating ..."
    ${MSISENSOR_PRO_SCAN} -d "${FILE_GENOME}" -o "${MICROSATELLITE_SITES}"
fi

${MSISENSOR_PRO} -d ${MICROSATELLITE_SITES} -n "${recalbamGD}" -t "${recalbamTD}" -o "${MSI_OUTPUT}"

${BAM_MATCHER} --bam1 ${recalbamGD} --bam2 ${recalbamTD} --output "${BAMMATCHER_OUTPUT}"

#eo VC