#!/bin/bash

# uncomment to load dependenties if required (e.g. if you're working on an HPC cluster)
# module load julia/1.11.3
# module load bwa
# module load bowtie
# module load star/2.7.11b
# module load blast/2.10.1+


CURRENT_DIRECTORY=$(pwd)

mkdir -p ${CURRENT_DIRECTORY}/data/alignment_output

INPUT_FASTQ=${CURRENT_DIRECTORY}/data/toy_sample2.fastq
INPUT_FASTA=${CURRENT_DIRECTORY}/data/toy_sample2.fasta

REF_FASTA=${CURRENT_DIRECTORY}/data/reference_ncRNAs.fasta
REF_LEVENMAP_DB=${CURRENT_DIRECTORY}/data/reference_ncRNAs_levenmap_database
REF_STAR_DB=${CURRENT_DIRECTORY}/data/reference_ncRNAs_STAR_database
REF_BLAST_DB=${CURRENT_DIRECTORY}/data/reference_ncRNAs_BLAST_database
REF_BOWTIE_INDEX=${CURRENT_DIRECTORY}/data/reference_ncRNAs_bowtie_index

OUTPUT_LEVENMAP=${CURRENT_DIRECTORY}/data/alignment_output/levenmap.txt
OUTPUT_BWA_L32=${CURRENT_DIRECTORY}/data/alignment_output/BWA_l32.sam
OUTPUT_BWA_L10=${CURRENT_DIRECTORY}/data/alignment_output/BWA_l10.sam
OUTPUT_BOWTIE_M1_L28=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m1_l28.sam
OUTPUT_BOWTIE_M1_L10=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m1_l10.sam
OUTPUT_BOWTIE_M0_L28=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m0_l28.sam
OUTPUT_BOWTIE_M0_L10=${CURRENT_DIRECTORY}/data/alignment_output/bowtie_m0_l10.sam
OUTPUT_STAR_PREFIX=${CURRENT_DIRECTORY}/data/alignment_output/STAR
OUTPUT_BLAST=${CURRENT_DIRECTORY}/data/alignment_output/BLAST.sam



########## build reference database for each aligner ########## 
julia ${CURRENT_DIRECTORY}/../src/build_reference_database.jl \
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




########## align the four samples with each aligner (i.e. levenmap, BWA, Bowtie, BLAST, and STAR ##########
julia ${CURRENT_DIRECTORY}/../src/levenmap.jl \
  --input_fastq $INPUT_FASTQ \
  --reference $REF_LEVENMAP_DB \
  --output $OUTPUT_LEVENMAP \
  --max_dist 1

bwa aln -l 32 -t 1 $REF_FASTA $INPUT_FASTQ > ${INPUT_FASTQ}_bwa_tmp_l32
bwa aln -l 10 -t 1 $REF_FASTA $INPUT_FASTQ > ${INPUT_FASTQ}_bwa_tmp_l10
bwa samse $REF_FASTA ${INPUT_FASTQ}_bwa_tmp_l32 $INPUT_FASTQ > $OUTPUT_BWA_L32
bwa samse $REF_FASTA ${INPUT_FASTQ}_bwa_tmp_l10 $INPUT_FASTQ > $OUTPUT_BWA_L10
rm ${INPUT_FASTQ}_bwa_tmp_l32
rm ${INPUT_FASTQ}_bwa_tmp_l10

bowtie -n 0 -l 28 --norc --best --strata -m 1 -x $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_M1_L28
bowtie -n 0 -l 10 --norc --best --strata -m 1 -x $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_M1_L10
bowtie -n 0 -l 28 --norc --best -x $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_M0_L28
bowtie -n 0 -l 10 --norc --best -x $REF_BOWTIE_INDEX $INPUT_FASTQ -S $OUTPUT_BOWTIE_M0_L10

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

