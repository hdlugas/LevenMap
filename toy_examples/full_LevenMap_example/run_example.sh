#!/bin/bash

# uncomment to load julia if required (e.g. if you're on a HPC cluster); if doing so, you may want to specify version other than 1.12.4
# module load julia/1.12.4

# get path to current working directory
CURRENT_DIRECTORY=$(pwd)

# create output directory
mkdir -p ${CURRENT_DIRECTORY}/data/alignment_output


########## build reference database ########## 
julia ${CURRENT_DIRECTORY}/../../src/BuildReferenceDatabase.jl \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs.fasta \
  --output ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --min_nucleotides 16 \
  --max_nucleotides 120 \
  --max_dist 1




########## extract UMI sequence, filter reads based on length, and align the four samples ##########
echo -e "\nSample 1:"
umi_tools extract \
  --stdin=${CURRENT_DIRECTORY}/data/sample1.fastq \
  --log=${CURRENT_DIRECTORY}/data/sample1_umi_extraction_log_file.txt \
  --stdout=${CURRENT_DIRECTORY}/data/sample1_post_umi_extraction.fastq \
  --extract-method=regex \
  --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'

cutadapt \
  --minimum-length=16 \
  --maximum-length=120 \
  -o ${CURRENT_DIRECTORY}/data/sample1_post_cutadapt.fastq \
  ${CURRENT_DIRECTORY}/data/sample1_post_umi_extraction.fastq > ${CURRENT_DIRECTORY}/data/sample1_cutadapt_log_file.txt

julia --threads 1 ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/sample1_post_cutadapt.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/sample1_alignment.txt \
  --max_dist 1 \
  --threads 1
echo ""



echo -e "\nSample 2:"
umi_tools extract \
  --stdin=${CURRENT_DIRECTORY}/data/sample2.fastq \
  --log=${CURRENT_DIRECTORY}/data/sample2_umi_extraction_log_file.txt \
  --stdout=${CURRENT_DIRECTORY}/data/sample2_post_umi_extraction.fastq \
  --extract-method=regex \
  --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'

cutadapt \
  --minimum-length=16 \
  --maximum-length=120 \
  -o ${CURRENT_DIRECTORY}/data/sample2_post_cutadapt.fastq \
  ${CURRENT_DIRECTORY}/data/sample2_post_umi_extraction.fastq > ${CURRENT_DIRECTORY}/data/sample2_cutadapt_log_file.txt

julia --threads 1 ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/sample2_post_cutadapt.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/sample2_alignment.txt \
  --max_dist 1 \
  --threads 1
echo ""



echo -e "\nSample 3:"
umi_tools extract \
  --stdin=${CURRENT_DIRECTORY}/data/sample3.fastq \
  --log=${CURRENT_DIRECTORY}/data/sample3_umi_extraction_log_file.txt \
  --stdout=${CURRENT_DIRECTORY}/data/sample3_post_umi_extraction.fastq \
  --extract-method=regex \
  --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'

cutadapt \
  --minimum-length=16 \
  --maximum-length=120 \
  -o ${CURRENT_DIRECTORY}/data/sample3_post_cutadapt.fastq \
  ${CURRENT_DIRECTORY}/data/sample3_post_umi_extraction.fastq > ${CURRENT_DIRECTORY}/data/sample3_cutadapt_log_file.txt

julia --threads 1 ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/sample3_post_cutadapt.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/sample3_alignment.txt \
  --max_dist 1 \
  --threads 1
echo ""



echo -e "\nSample 4:"
umi_tools extract \
  --stdin=${CURRENT_DIRECTORY}/data/sample4.fastq \
  --log=${CURRENT_DIRECTORY}/data/sample4_umi_extraction_log_file.txt \
  --stdout=${CURRENT_DIRECTORY}/data/sample4_post_umi_extraction.fastq \
  --extract-method=regex \
  --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'

cutadapt \
  --minimum-length=16 \
  --maximum-length=120 \
  -o ${CURRENT_DIRECTORY}/data/sample4_post_cutadapt.fastq \
  ${CURRENT_DIRECTORY}/data/sample4_post_umi_extraction.fastq > ${CURRENT_DIRECTORY}/data/sample4_cutadapt_log_file.txt

julia --threads 1 ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/sample4_post_cutadapt.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/sample4_alignment.txt \
  --max_dist 1 \
  --threads 1
echo ""



########## consolidate alignment information to create count matrix with rows corresponding to ncRNAs and columns corresponding to samples #########
julia ${CURRENT_DIRECTORY}/../../src/CreateCountMatrix.jl \
  --input_dir ${CURRENT_DIRECTORY}/data/alignment_output
echo ""

rm ${CURRENT_DIRECTORY}/data/*cutadapt*
rm ${CURRENT_DIRECTORY}/data/*umi*


