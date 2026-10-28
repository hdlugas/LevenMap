#!/bin/bash

CURRENT_DIRECTORY=$(pwd)

mkdir -p ${CURRENT_DIRECTORY}/data/alignment_output


julia ${CURRENT_DIRECTORY}/../../src/BuildReferenceDatabase.jl \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs.fasta \
  --output ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --min_nucleotides 16 \
  --max_nucleotides 30 \
  --max_dist 1


echo -e "\nTest case 1: alignment reports edit distance, strand, and multiple optimal alignments (expected: alignment succeeds)"
echo -e "\n FASTQ file of input reads:"
cat ${CURRENT_DIRECTORY}/data/reads_for_ties.fastq

echo -e "\n FASTA file of reference sequences:"
cat ${CURRENT_DIRECTORY}/data/reference_ncRNAs_for_ties.fasta

julia ${CURRENT_DIRECTORY}/../../src/BuildReferenceDatabase.jl \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_for_ties.fasta \
  --output ${CURRENT_DIRECTORY}/data/reference_ncRNAs_for_ties_database \
  --min_nucleotides 16 \
  --max_nucleotides 30 \
  --max_dist 1
julia ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/reads_for_ties.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_for_ties_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_1_alignment.txt \
  --max_dist 1
echo -e "\nAlignment output (with --max_dist 1):"
cat ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_1_alignment.txt


echo -e "\n\n\nTest case 2: parameters match the reference database (expected: alignment succeeds)"
julia ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/reads_in_range.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_2_alignment.txt \
  --max_dist 1
echo -e "\nAlignment successful."


echo -e "\n\n\nTest case 3: --max_dist 2 with a database built with --max_dist 1 (error message expected)"
julia ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/reads_in_range.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_3_alignment.txt \
  --max_dist 2


echo -e "\n\n\nTest case 4: reads of length 14 and 35 with a database covering lengths 16-30 (error message expected)"
julia ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/reads_out_of_range.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_4_alignment.txt \
  --max_dist 1


echo -e "\n\n\nTest case 5: reference database directory does not exist (error message expected)"
julia ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/reads_in_range.fastq \
  --reference ${CURRENT_DIRECTORY}/data/nonexistent_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_5_alignment.txt \
  --max_dist 1


echo -e "\n\n\nTest case 6: reference database directory exists but contains no reference files (error message expected)"
mkdir -p ${CURRENT_DIRECTORY}/data/dummy_database
julia ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/reads_in_range.fastq \
  --reference ${CURRENT_DIRECTORY}/data/dummy_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_6_alignment.txt \
  --max_dist 1


echo -e "\n\n\nTest case 7: negative --max_dist (error message expected)"
julia ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq ${CURRENT_DIRECTORY}/data/reads_in_range.fastq \
  --reference ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database \
  --output ${CURRENT_DIRECTORY}/data/alignment_output/edge_case_7_alignment.txt \
  --max_dist -1




echo -e "\nEdge case testing finished.\n"

rm -r ${CURRENT_DIRECTORY}/data/reference_ncRNAs_database
rm -r ${CURRENT_DIRECTORY}/data/reference_ncRNAs_for_ties_database
rm -r ${CURRENT_DIRECTORY}/data/dummy_database
rm -r ${CURRENT_DIRECTORY}/data/alignment_output

