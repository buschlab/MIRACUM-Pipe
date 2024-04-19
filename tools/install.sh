#!/usr/bin/env bash

# variables
readonly VERSION_TRIMMOMATIC="0.39"
readonly VERSION_PICARD="2.27.1"
readonly VERSION_BEDTOOLS="2.30.0"
readonly VERSION_VEP_MAJOR="111"
readonly VERSION_VEP="${VERSION_VEP_MAJOR}.0"

########
readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)


# remove old object files
find ${DIR_SCRIPT} -type f -name '*.o' -delete


##########
# FastQC #
##########
cd ${DIR_SCRIPT}/FastQC
ant build
chmod +x bin/fastqc


####### install from web #######


##########
# picard #
##########
mkdir -p ${DIR_SCRIPT}/picard
cd ${DIR_SCRIPT}/picard

wget https://github.com/broadinstitute/picard/releases/download/${VERSION_PICARD}/picard.jar \
    -O picard.jar

#############
# bedtools2 #
#############
cd ${DIR_SCRIPT}

wget https://github.com/arq5x/bedtools2/releases/download/v${VERSION_BEDTOOLS}/bedtools-${VERSION_BEDTOOLS}.tar.gz \
    -O bedtools2.tar.gz

tar -xzf bedtools2.tar.gz
rm -f bedtools2.tar.gz

cd bedtools2 && make

###############
# Trimmomatic #
###############
cd ${DIR_SCRIPT}

# download new version
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${VERSION_TRIMMOMATIC}.zip \
    -O trimmomatic.zip

# unzip
unzip -o trimmomatic.zip
rm -f trimmomatic.zip

# rename folder and file (neglect version information)
mv Trimmomatic* Trimmomatic
mv Trimmomatic/trimmomatic-${VERSION_TRIMMOMATIC}.jar Trimmomatic/trimmomatic.jar

#######
# VEP #
#######
cd ${DIR_SCRIPT}

wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/${VERSION_VEP}.tar.gz \
    -O vep.tar.gz
tar -xzf vep.tar.gz
rm -f vep.tar.gz
mv ensembl-vep-release-${VERSION_VEP} vep

# #  CADD

# mkdir -p ${DIR_SCRIPT}/../databases/vep/CADD_GRCh37
# wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz -O ${DIR_SCRIPT}/../databases/vep/CADD_GRCh37/whole_genome_SNVs.tsv.gz
# wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz.tbi -O ${DIR_SCRIPT}/../databases/vep/CADD_GRCh37/whole_genome_SNVs.tsv.gz.tbi
# wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz -O ${DIR_SCRIPT}/../databases/vep/CADD_GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz
# wget https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz.tbi -O ${DIR_SCRIPT}/../databases/vep/CADD_GRCh37/gnomad.genomes-exomes.r4.0.indel.tsv.gz.tbi

# # REVEL

# mkdir -p ${DIR_SCRIPT}/../databases/vep/REVEL
# wget https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip -O ${DIR_SCRIPT}/../databases/vep/REVEL/revel-v1.3_all_chromosomes.zip
# cd ${DIR_SCRIPT}/../databases/vep/REVEL
# unzip revel-v1.3_all_chromosomes.zip
# rm revel-v1.3_all_chromosomes.zip
# cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
# sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
# bgzip new_tabbed_revel.tsv
# tabix -f -s 1 -b 2 -e 2 new_tabbed_revel.tsv.gz

# VEP INSTALL

cd ${DIR_SCRIPT}/vep
perl INSTALL.pl -n -a alcfp -s homo_sapiens_refseq -y GRCh37 -g CADD,REVEL

# Index reference fasta

cd ~/.vep/homo_sapiens_refseq/${VERSION_VEP_MAJOR}_GRCh37
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
bgzip --index Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

# VCF2MAF
cd ${DIR_SCRIPT}
wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.21.tar.gz -O vcf2maf.tar.gz
rm vcf2maf.tar.gz
mv vcf2maf-* vcf2maf

###### COMPILE SUBMODULES #######

#########
# FREEC #
#########
# compile FREEC
cd ${DIR_SCRIPT}/FREEC/src && make

# remove object files
rm -f *.o depend.mk

# copy binary
cd ../
mkdir -p bin
chmod +x src/freec
mv src/freec bin

#################
# bam-readcount #
#################
cd ${DIR_SCRIPT}/bam-readcount
mkdir -p build && cd build
cmake ../ && make

# move file
mv ${DIR_SCRIPT}/bam-readcount/build/bin ${DIR_SCRIPT}/bam-readcount/bin
rm -rf ${DIR_SCRIPT}/bam-readcount/build

#######
# bwa-mem2 #
#######
cd ${DIR_SCRIPT}
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 -O bwa-mem2.tar.bz2
tar -xf bwa-mem2.tar.bz2
rm bwa-mem2*.tar.bz2
mv bwa-mem2* bwa-mem2

##########
# htslib #
##########
cd ${DIR_SCRIPT}/htslib
autoheader     # If using configure, generate the header template...
autoconf       # ...and configure script (or use autoreconf to do both)

./configure    # Optional but recommended, for choosing extra functionality


############
# samtools #
############
cd ${DIR_SCRIPT}/samtools

# create config
autoheader
autoconf -Wno-syntax

# configure and make
./configure

# build samtools and htslib
make

rm -f ${DIR_SCRIPT}/samtools/*.o

# htslib is built by samtools
rm -f ${DIR_SCRIPT}/htslib/*.o


# add lib folder system wide
echo "$DIR_SCRIPT/htslib" > /etc/ld.so.conf.d/htslib.conf


#################
# fusioncatcher #
#################
#cd ${DIR_SCRIPT}
#wget http://sf.net/projects/fusioncatcher/files/bootstrap.py -O bootstrap.py
#python bootstrap.py --prefix=${DIR_SCRIPT} -t -y
cd ${DIR_SCRIPT}/fusioncatcher/tools/
./install_tools.sh

##################
# sequenza-utils #
##################
#pip3 install sequenza-utils

##############
# msisensor2 #
##############
cd ${DIR_SCRIPT}/msisensor2
chmod +x msisensor2

#################
# msisensor-pro #
#################
cd ${DIR_SCRIPT}/msisensor-pro
#wget https://github.com/xjtu-omics/msisensor-pro/raw/master/binary/msisensor-pro
#chmod +x msisensor-pro
./INSTALL

############
# agfusion #
############
cd /opt/MIRACUM-Pipe/databases
agfusion download -g hg38

###############
# bam-matcher #
###############
cd ${DIR_SCRIPT}/bam-matcher
pip3 install --upgrade numpy
pip3 install -r requirements.txt

cat >"${DIR_SCRIPT}"/bam-matcher/bam-matcher.conf <<EOI
[VariantCallers]
caller:    gatk4
gatk3:     /opt/MIRACUM-Pipe/tools/gatk/GenomeAnalysisTK.jar
gatk4:     /opt/MIRACUM-Pipe/tools/gatk4/gatk
freebayes: freebayes
samtools:  samtools
varscan:   /opt/MIRACUM-Pipe/tools/bam-matcher/VarScan.jar
java:      java

[ScriptOptions]
DP_threshold: 15
number_of_SNPs:
fast_freebayes: True
VCF_file: /opt/MIRACUM-Pipe/tools/bam-matcher/1kg.exome.highAF.1511.vcf


[VariantCallerParameters]
GATK_MEM:    4
GATK_nt:     1
VARSCAN_MEM: 4

[GenomeReference]
REFERENCE: /opt/MIRACUM-Pipe/assets/references/genome/genome.fa
REF_ALTERNATE:
CHROM_MAP:

[BatchOperations]
CACHE_DIR: ${DIR_TMP}
[Miscellaneous]
EOI

awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${DIR_SCRIPT}/bam-matcher/1kg.exome.highAF.1511.vcf > ${DIR_SCRIPT}/bam-matcher/tmp.vcf
mv ${DIR_SCRIPT}/bam-matcher/tmp.vcf ${DIR_SCRIPT}/bam-matcher/1kg.exome.highAF.1511.vcf

awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${DIR_SCRIPT}/bam-matcher/1kg.exome.highAF.3680.vcf > ${DIR_SCRIPT}/bam-matcher/tmp.vcf
mv ${DIR_SCRIPT}/bam-matcher/tmp.vcf ${DIR_SCRIPT}/bam-matcher/1kg.exome.highAF.3680.vcf

awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${DIR_SCRIPT}/bam-matcher/1kg.exome.highAF.7550.vcf > ${DIR_SCRIPT}/bam-matcher/tmp.vcf
mv ${DIR_SCRIPT}/bam-matcher/tmp.vcf ${DIR_SCRIPT}/bam-matcher/1kg.exome.highAF.7550.vcf
