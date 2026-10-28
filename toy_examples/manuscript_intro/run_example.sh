#!/bin/bash

# uncomment to load dependenties if required (e.g. if you're working on an HPC cluster)
# module load julia
# module load bwa
# module load bowtie
# module load star
# module load blast


CURRENT_DIRECTORY=$(pwd)

mkdir -p ${CURRENT_DIRECTORY}/data/alignment_output

INPUT_FASTQ=${CURRENT_DIRECTORY}/data/toy_sample.fastq
INPUT_FASTA=${CURRENT_DIRECTORY}/data/toy_sample.fasta

REF_FASTA=${CURRENT_DIRECTORY}/data/reference_ncRNAs.fasta
REF_LEVENMAP_DB=${CURRENT_DIRECTORY}/data/reference_ncRNAs_LevenMap_database
REF_STAR_DB=${CURRENT_DIRECTORY}/data/reference_ncRNAs_STAR_database
REF_BLAST_DB=${CURRENT_DIRECTORY}/data/reference_ncRNAs_BLAST_database
REF_BOWTIE_INDEX=${CURRENT_DIRECTORY}/data/reference_ncRNAs_bowtie_index

OUTPUT_LEVENMAP=${CURRENT_DIRECTORY}/data/alignment_output/LevenMap.txt
OUTPUT_BWA_POTLA=${CURRENT_DIRECTORY}/data/alignment_output/BWA_potla.sam
OUTPUT_BOWTIE_M1_L28=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m1_l28.sam
OUTPUT_BOWTIE_M1_L10=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m1_l10.sam
OUTPUT_BOWTIE_M0_L28=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m0_l28.sam
OUTPUT_BOWTIE_M0_L10=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m0_l10.sam
OUTPUT_BOWTIE_POTLA=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_potla.sam
OUTPUT_BOWTIE_MIRDEEP2=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_mirdeep2.sam
OUTPUT_BOWTIE_MIRGE_SET1=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_mirge_set1.sam
OUTPUT_BOWTIE_MIRGE_SET2=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_mirge_set2.sam
OUTPUT_BOWTIE_MIRGE_SET3=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_mirge_set3.sam
OUTPUT_BOWTIE_MIRGE_SET4=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_mirge_set4.sam
OUTPUT_BOWTIE_MIRGE_SET5=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_mirge_set5.sam
OUTPUT_BOWTIE_SRNABENCH0=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_srnabench0.sam
OUTPUT_BOWTIE_SRNABENCH1=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_srnabench1.sam
OUTPUT_BOWTIE_SRNABENCH2=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_srnabench2.sam
OUTPUT_BOWTIE_SHORTSTACK=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_shortstack.sam
OUTPUT_STAR_PREFIX=${CURRENT_DIRECTORY}/data/alignment_output/STAR
OUTPUT_BLAST=${CURRENT_DIRECTORY}/data/alignment_output/BLAST.sam



########## build reference database for each aligner ########## 
julia ${CURRENT_DIRECTORY}/../../src/BuildReferenceDatabase.jl \
  --reference $REF_FASTA \
  --output $REF_LEVENMAP_DB \
  --min_nucleotides 16 \
  --max_nucleotides 120 \
  --max_dist 1

bwa index $REF_FASTA

bowtie-build $REF_FASTA $REF_BOWTIE_INDEX

STAR \
  --runThreadN 1 \
  --runMode genomeGenerate \
  --genomeDir $REF_STAR_DB \
  --genomeFastaFiles $REF_FASTA \
  --limitGenomeGenerateRAM 205254557962

makeblastdb -in $REF_FASTA -dbtype nucl -out $REF_BLAST_DB




########## align the four samples with each aligner (i.e. LevenMap, BWA, Bowtie, BLAST, and STAR ##########
julia --threads auto ${CURRENT_DIRECTORY}/../../src/LevenMap.jl \
  --input_fastq $INPUT_FASTQ \
  --reference $REF_LEVENMAP_DB \
  --output $OUTPUT_LEVENMAP \
  --max_dist 1

bwa aln -l 30 -t 1 $REF_FASTA $INPUT_FASTQ > ${INPUT_FASTQ}_bwa_tmp
bwa samse $REF_FASTA ${INPUT_FASTQ}_bwa_tmp $INPUT_FASTQ > $OUTPUT_BWA_POTLA
rm ${INPUT_FASTQ}_bwa_tmp


bowtie -n 0 -l 30 --norc --best --strata -m 1 $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_POTLA
bowtie -n 0 -e 80 -l 18 -a -m 5 --best --strata $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_MIRDEEP2
bowtie -n 0 --norc $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_MIRGE_SET1
bowtie -n 1 --norc $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_MIRGE_SET2
bowtie -v 1 -a --best --strata --norc $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_MIRGE_SET3
bowtie -v 0 -a --best --strata --norc $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_MIRGE_SET4
bowtie -5 1 -3 2 -v 2 --best --norc $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_MIRGE_SET5
bowtie -v 0 $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_SRNABENCH0
bowtie -v 1 $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_SRNABENCH1
bowtie -v 2 $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_SRNABENCH2
bowtie -v 1 --best --strata -k 50 $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_SHORTSTACK

blastn -query $INPUT_FASTA -db $REF_BLAST_DB -outfmt 17 -out $OUTPUT_BLAST

STAR \
  --runThreadN 1 \
  --genomeDir $REF_STAR_DB \
  --readFilesIn $INPUT_FASTQ \
  --outFileNamePrefix $OUTPUT_STAR_PREFIX \
  --outSAMtype SAM \
  --outSAMunmapped Within \
  --outSAMmultNmax 1 \
  --outSAMattributes Standard

rm ${CURRENT_DIRECTORY}/data/alignment_output/STARLog.final.out
rm ${CURRENT_DIRECTORY}/data/alignment_output/STARLog.progress.out
rm ${CURRENT_DIRECTORY}/data/alignment_output/STARLog.out
rm ${CURRENT_DIRECTORY}/data/alignment_output/STARSJ.out.tab
mv ${CURRENT_DIRECTORY}/data/alignment_output/STARAligned.out.sam ${CURRENT_DIRECTORY}/data/alignment_output/STAR.sam

