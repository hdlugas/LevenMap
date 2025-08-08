# levenmap

Julia-based alignment algorithm designed specifically for aligning small RNA-sequencing data to database of non-coding RNAs (ncRNAs).

## Table of Contents
- [1. Motivation for levenmap](#motivation)
- [2. Installation](#installation)
- [3. Usage](#usage)
   - [3.1. Suggested preprocessing](#preprocessing)
   - [3.2 Build reference database](#build-ref-db)
   - [3.3 Run levenmap](#run-aligner)
   - [3.4 Create count matrix](#create-count-matrix)
- [4. Downstream analysis](#downstream-analysis)

<a name="motivation"></a>
## 1. Motivation for levenmap
Common bioinformatics workflows which process small RNA-sequencing data typically involve using an aligner such as Burrows-Wheeler Aligner (BWA) or Bowtie to align reads to a database of non-coding RNAs (ncRNAs) such as miRbase. In a database of ncRNAs such as miRbase for example, it is often the case that a hairpin microRNA (miRNA) contains an subsequence equivalent to a mature miRNA (e.g. hsa-mir-12117=GTCTCAGT**GAAGTGGAGCACATCAGTGA**AAGGGTGAACTTAACCTTTCACTGGTGTGCTCCATCTCACTCAGAC contains hsa-miR-12117=GAAGTGGAGCACATCAGTGA). If a small RNA-sequencing read is exactly equivalent to a mature miRNA such as hsa-miR-12117 after the unique molecular identifier (UMI) is removed and either BWA or Bowtie is used to align the given read to a database containing both the mature and hairpin miRNA, then the read will often be aligned to the hairpin miRNA even though the read exactly matches the mature miRNA. This is the case because BWA, Bowtie, and other common aligners are designed to align reads to the human genome where (i) reference sequences are not exact subsequences of each other (e.g. chromosome 13 is not an exact subsequence of chromosome 17) and (ii) reference sequences are orders of magnitude longer than the reads produced by the sequencing procedure (e.g. high-throughput sequencing does not produce reads corresponding to an entire chromosome). The alignment tool levenmap is a more biologically-relevant aligner specifically for aligning small RNA-sequencing data to a database of ncRNAs which accounts for the nuances of aligning to short (i.e. reads are on same scale as reference sequences), nested reference sequences. This aligner works by returning the alignment(s) for a given read which have the smallest Levenshtein distance up to a user-specified maximum Levenshtein distance (i.e. a user-specified maximum number of mismatches to allow).

<a name="installation"></a>
## 2. Installation
The Julia packages required to run levenmap are ArgParse, DataStructures, FASTX, Glob, and StringDistances. These dependencies can be installed with the Julia command:
```
using Pkg; Pkg.add(["ArgParse", "FASTX", "DataStructures", "StringDistances", "Glob"])
```

However, if you would like to (i) extract the unique molecular identifier (UMI) sequences from reads and (ii) filter reads based on length (i.e. number of nucleotides) as is recommended by both us and Kapoor et al in "A bioinformatics approach to microRNA-sequencing analysis" ([https://www.sciencedirect.com/science/article/pii/S266591312030131X](https://www.sciencedirect.com/science/article/pii/S266591312030131X)), then umi_tools ([https://github.com/CGATOxford/UMI-tools](https://github.com/CGATOxford/UMI-tools)) and cutadapt ([https://cutadapt.readthedocs.io/en/stable/index.html](https://cutadapt.readthedocs.io/en/stable/index.html)) must also be installed. We recommend creating a conda environment with the necessary dependencies, which can be done by executing the following command once conda is installed on your system, you have downloaded this repository, and have navigated to the directory that contains the environment.yml file:
```
conda env create -f environment.yml
```
To activate the levenmap_env conda environment, run:
```
conda activate levenmap_env
```
To deactivate the levenmap_env conda environement, run:
```
conda deactivate
```

Alternatively, one can create a conda environment, activate it, and install each dependency with:
```
conda create -n levenmap_env
conda activate levenmap_
conda install -c bioconda -c conda-forge umi_tools
conda install bioconda::cutadapt
```

For instructions on installing conda on your system, see [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

<a name="usage"></a>
## 3. Usage
For an example of the complete recommended workflow, clone the repository, navigate to the toy_example directory, activate the levenmap_env conda environment, and run the script example.sh using:
```
./example.sh
```

<a name="preprocessing"></a>
### 3.1. Suggested preprocessing
We recommend performing the first two steps of data preprocessing proposed by Kapoor et al in "A bioinformatics approach to microRNA-sequencing analysis" ([https://www.sciencedirect.com/science/article/pii/S266591312030131X](https://www.sciencedirect.com/science/article/pii/S266591312030131X)) which involve (i) extracting unique molecular identifiers (UMIs) from reads in the original FASTQ files and (ii) filtering reads to keep only reads with a length (i.e. number of nucleotides) in a specified range. If human hairpin and mature miRNAs from miRbase are used as the reference database, then we recommend retaining reads with length 16-120, inclusive, because this is the range of nucleotide sequence length in miRbase's human hairpin and mature miRNAs. The following commands process the original FASTQ file $RAW_FASTQ into the FASTQ file $POST_CUTADAPT_FASTQ which can be used as input for agrepNestAligner:
```
umi_tools extract \
  --stdin=$RAW_FASTQ \
  --log=$UMI_EXTRACTION_LOG_FILE \
  --stdout=$UMI_EXTRACTED_FASTQ \
  --extract-method=regex \
  --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'

cutadapt \
  --minimum-length=16 \
  --maximum-length=120 \
  -o $POST_CUTADAPT_FASTQ \
  $UMI_EXTRACTED_FASTQ > $CUTADAPT_LOG_FILE
```

<a name="build-ref-db"></a>
### 3.2. Build reference database
Building a reference database from the reference FASTA file greatly increases speed and is a necessary step in using levenmap. If the length of the shortest/longest read in the reference FASTA file is S/L and the maximum number of mismatches allowed is N, then a FASTA file with all reads of length M-N to M+N for each M=S,S+1,S+2,...,L-2,L-1,L is created, resulting in (L-S)+1 FASTA files. For a given read in some sample's FASTQ file with length A, the specific FASTA file with all reads of length A-N to A+N is used as the reference to avoid regular expression matching with agrep on reference sequences that are either too short or too long to correspond to the given read. 

To view the usage instructions and parameter descriptions of the script build_reference_database.sh, run:
```
./build_reference_database.sh -h
```

The resulting message is:
```
Usage: ./build_reference_database.sh --min-nucleotides <INT> --max-nucleotides <INT> -x <INT> -r <reference_FASTA_file> -o <output_directory>

Builds reference database corresponding to a reference FASTA file. This step is necessary prior to running agrepNestAligner.

Options:
  --min-nucleotides <INT>                          Integer value of the minimum number of nucleotides each reference sequence must contain to be considered, inclusive. Default = 18.
  --max-nucleotides <INT>                          Integer value of the maximum number of nucleotides each reference sequence must contain to be considered, inclusive. Default = 30.
  -x <maximum number of mismatches allowed>        Integer specifying the maximum number of mismatches allowed during alignment. Default = 0.
  -r <reference_FASTA_file>                        Specifies the FASTA file of reference sequences. Mandatory argument.
  -o <output_directory>                            Specifies the directory the reference database will be written to. Mandatory argument.
  -h                                               Shows this help message.
```

For example, the following writes a reference database to $REFERENCE_ncRNA_DATABASE_DIRECTORY with reference ncRNA FASTA file $REFERENCE_ncRNAs_FASTA:
```
./build_reference_database.sh \
  --min-nucleotides 16 \
  --max-nucleotides 120 \
  -x 1 \
  -r $REFERENCE_ncRNAS_FASTA \
  -o $REFERENCE_ncRNA_DATABASE_DIRECTORY
```

<a name="run-aligner"></a>
### 3.3. Run levenmap
Once the reference database is constructed, agrepNestAligner.sh can be used to align the reads of a given FASTQ file to the reference database. To view the usage instructions and parameter descriptions of the script agrepNestAligner.sh, run:
```
./agrepNestAligner.sh -h
```

The resulting message is:
```
Usage: ./agrepNestAligner.sh -x <INT> -@ <INT> -i <input_FASTQ_file> -r <reference_database_directory> -o <output_directory>

Aligns reads in a given FASTQ file to a reference database of ncRNAs using agrep.

Options:
  -x <maximum number of mismatches allowed>        Integer specifying the maximum number of mismatches allowed during alignment. Default = 0.
  -@ <number of threads>                           Integer specifying the number of threads to utilize. Default = 1.
  -i <input_FASTQ_file>                            Specifies the input FASTQ file. Mandatory argument.
  -r <reference_database_directory>                Path to the directory containing the reference database created from build_reference_database.sh. Mandatory argument.
  -o <output_directory>                            Specifies the directory the output TXT file should be written to. Mandatory argument.
  -h                                               Shows this help message.
```

For example, the following aligns the FASTA file $POST_CUTADAPT_FASTQ to the reference database $REFERENCE_ncRNA_DATABASE_DIRECTORY and writes the output to the directory $ALIGNMENT_OUTPUT_DIRECTORY:
```
./agrepNestAligner.sh \
  -x 1 \
  -@ 1 \
  -i $POST_CUTADAPT_FASTQ \
  -r $REFERENCE_ncRNA_DATABASE_DIRECTORY \
  -o $ALIGNMENT_OUTPUT_DIRECTORY
```

<a name="create-count-matrix"></a>
### 3.4. Create count matrix
Once all samples are aligned with output written to the directory $ALIGNMENT_OUTPUT_DIRECTORY, the alignment information is consolidated into a count matrix with rows corresponding to ncRNAs and columns corresponding to samples. To view the usage instructions and parameter description of the script create_count_matrix.sh, run:
```
./create_count_matrix.sh -h
```

The resulting message is:
```
Usage: ./create_count_matrix.sh -i <input_directory>

Consolidates agrepNestAligner alignment information into count matrix with rows corresponding to ncRNAs and columns corresponding to samples. The tab-delimited text file count_matrix.txt will be written to the input_directory argument.

Options:
  -i <input_directory>                             Specifies the directory the agrepNestAligner output was written to. Mandatory argument.
  -h                                               Shows this help message.
```

For example, the following writes the count matrix to the directory $ALIGNMENT_OUTPUT_DIRECTORY once all agrepNestAligner alignments have been written to $ALIGNMENT_OUTPUT_DIRECTORY:
```
./create_count_matrix.sh -i $ALIGNMENT_OUTPUT_DIRECTORY
```

<a name="downstream-analysis"></a>
## 4. Downstream analysis
Once the count matrix with rows corresponding to ncRNAs and columns corresponding to samples is computed, there are a variety of tools available which perform differential transcription analysis to identify ncRNAs transcribed in significantly different levels among sample groups of interest, most of which use negative binomial regression. Some options include DESeq2, edgeR, or limma.
