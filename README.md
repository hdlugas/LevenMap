# levenmap: alignment of small RNA-sequencing data to database of non-coding RNAs

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
conda activate levenmap_env
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
pip install cutadapt
conda install umi_tools
```

Note that this assumes pip is installed on your system. For instructions on installing conda on your system, see [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

<a name="usage"></a>
## 3. Usage
For an example of the complete recommended workflow, clone the repository, navigate to the toy_example directory, activate the levenmap_env conda environment, and run the script example.sh using:
```
./example.sh
```
The levenmap aligner has been tested and validated on Ubuntu 20.04, Ubuntu 22.04, and CentOS 7.

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
Building a reference database from the reference FASTA file greatly increases speed and is a necessary step in using levenmap. If the length of the shortest/longest read in the reference FASTA file is S/L and the maximum number of mismatches allowed is N, then a FASTA file with all reads of length M-N to M+N for each M=S,S+1,S+2,...,L-2,L-1,L is created, resulting in (L-S)+1 FASTA files. For a given read in some sample's FASTQ file with length A, the specific FASTA file with all reads of length A-N to A+N is used as the reference to avoid computing the Levenshtein distance with reference sequences that are either too short or too long to correspond to the given read.

To view the usage instructions and parameter descriptions of the script build_reference_database.sh, run:
```
julia build_reference_database.jl --help
```

The resulting message is:
```
Help Message for build_reference_database.jl

Usage:
  julia build_reference_database.jl --reference <path> --output <path> [--min_nucleotides <int>] [--max_nucleotides <int>] [--max_dist <float>]

Arguments:
  --reference          Path to FASTA file of reference ncRNAs (required)
  --output             Path to output directory containing processed reference database (required)
  --min_nucleotides    Minimum length (in units of nucleotides) of reference sequences to include in reference database (optional, Int, default = 12)
  --max_nucleotides    Maximum length (in units of nucleotides) of reference sequences to include in reference database (optional, Int, default = 30)
  --max_dist           Maximum Levenshtein distance allowed in alignment (optional, Int, default = 5)
  --help               Show this help message
```

For example, the following writes a reference database to $REFERENCE_ncRNA_DATABASE_DIRECTORY with reference ncRNA FASTA file $REFERENCE_ncRNAs_FASTA:
```
julia build_reference_database.jl \
  --reference $REFERENCE_ncRNAs_FASTA \
  --output $REFERENCE_ncRNA_DATABASE_DIRECTORY \
  --min_nucleotides 16 \
  --max_nucleotides 120 \
  --max_dist 3
```

<a name="run-aligner"></a>
### 3.3. Run levenmap
Once the reference database is constructed, levenmap.jl can be used to align the reads of a given FASTQ file to the reference database. To view the usage instructions and parameter descriptions of the script levenmap.jl, run:
```
julia levenmap.jl --help
```

The resulting message is:
```
Usage:
  julia levenmap.jl --input_fastq <path> --reference <path> --output <path> [--max_dist <float>] [--threads <int>]

Arguments:
  --input_fastq      Path to input FASTQ file (required)
  --reference        Path to reference database output from build_reference_database.jl (required)
  --output           Path to output TXT file (required)
  --max_dist         Maximum Levenshtein distance allowed in alignment (optional, Int, default = 1)
  --threads          Number of threads (optional, Int, default = 1)
  --help             Show this help message
```

For example, the following aligns the FASTA file $POST_CUTADAPT_FASTQ to the reference database $REFERENCE_ncRNA_DATABASE_DIRECTORY allowing for at most 1 mismatch and writes the output alignment text file to $ALIGNMENT_OUTPUT:
```
julia levenmap.jl \
  --input_fastq $POST_CUTADAPT_FASTQ \
  --reference $REFERENCE_ncRNA_DATABASE_DIRECTORY \
  --output ALIGNMENT_OUTPUT \
  --max_dist 1
```

<a name="create-count-matrix"></a>
### 3.4. Create count matrix
Once all samples are aligned with output written to the directory $ALIGNMENT_OUTPUT_DIRECTORY, the alignment information is consolidated into a count matrix with rows corresponding to ncRNAs and columns corresponding to samples. To view the usage instructions and parameter description of the script create_count_matrix.sh, run:
```
julia create_count_matrix.jl --help
```

The resulting message is:
```
Help Message for create_count_matrix.jl

Usage:
  julia create_count_matrix.jl --input <path>

Arguments:
  --input_dir   Path to directory containing alignment output from all samples (required)
  --help    Show this help message
```

For example, the following writes the count matrix to the directory $ALIGNMENT_OUTPUT_DIRECTORY once all agrepNestAligner alignments have been written to $ALIGNMENT_OUTPUT_DIRECTORY:
```
julia create_count_matrix.jl --input_dir $ALIGNMENT_OUTPUT_DIRECTORY
```

<a name="downstream-analysis"></a>
## 4. Downstream analysis
Once the count matrix with rows corresponding to ncRNAs and columns corresponding to samples is computed, there are a variety of tools available which perform differential transcription analysis to identify ncRNAs transcribed in significantly different levels among sample groups of interest. Some options include DESeq2, edgeR, or limma. Alternatively, one can also directly implement their desired hypothesis tests (after normalizing the count matrix) without the use of a package designed specifically for differential expression analysis.
