# Copyright (C) 2014 by Per Unneberg

"""Global variables for snakemakelib.

The variables listed here define the keys for the config dictionary.

"""
##################################################
# General
##################################################

##################################################
# Sections and variables
# 
# Section names follow the convention that if a rules file is named
# dir/subdir/rules.rules, the corresponding variable assignment will
# be DIR_SUBDIR_RULES = "dir.subdir.rules"
#
##################################################

##################################################
# comp
##################################################

##################################################
# comp.settings
##################################################
COMP_SETTINGS = "comp.settings"
AWK = "awk"
CAT = "cat"
GREP = "grep"
GZIP = "gzip"
HEAD = "head"
ZCAT = "zcat"

##################################################
# bio
##################################################

##################################################
# bio.genetics
##################################################

##################################################
# bio.genetics.association.plink
##################################################

BIO_GENETICS_ASSOCIATION_PLINK = "bio.genetics.association.plink"
PLINK = "plink"


##################################################
# bio.ngs
##################################################

##############################
# bio.ngs.settings
##############################
BIO_NGS_SETTINGS = "bio.ngs.settings"
ANNOTATION = "annotation"
ANNOT_LABEL = "annot_label"
BAIT_REGIONS = "bait_regions"
BAM_LIST = "bam_list"
CHR = "chr"
CMD = "cmd"
COMMON_OPTIONS = "common_options"
DB = "db"
DBSNP = "dbsnp"
FASTQ_SUFFIX = "fastq_suffix"
FLOWCELLS = "flowcells"
FLOWCELLRUNS = "flowcellruns"
FLOWCELL_SUFFIX = "flowcell_suffix"
HOME = "home"
INPUTDIR = "inputdir"
JAR = "jar"
JAVA = "java"
JAVA_MEM = "java_mem"
JAVA_TMPDIR = "java_tmpdir"
KNOWN_SITES = "known_sites"
LANES = "lanes"
MAIN = "main"
OPTIONS = "options"
READ1_LABEL = "read1_label"
READ2_LABEL = "read2_label"
READ1_SUFFIX = "read1_suffix"
READ2_SUFFIX = "read2_suffix"
REF = "ref"
RUN_ID_RE = "run_id_re"
SAMPLES = "samples"
SAMPLE_PREFIX = "sample_prefix"
SECTION = "section"
SEQUENCE_CAPTURE = "sequence_capture"
TARGET_REGIONS = "target_regions"
THREADS = "threads"
TRANSCRIPT_ANNOT_GTF = "transcript_annot_gtf"
VCFSUFFIX = "vcfsuffix"


##################################################
# bio.ngs.align
##################################################

##################################################
# bio.ngs.align.bwa
##################################################
BIO_NGS_ALIGN_BWA = "bio.ngs.align.bwa"
MEM = "mem"

##################################################
# bio.ngs.align.bowtie
##################################################
BIO_NGS_ALIGN_BOWTIE = "bio.ngs.align.bowtie"
BOWTIE = "bowtie"
BOWTIE2 = "bowtie2"

##################################################
# bio.ngs.methylseq.bismark
##################################################
BIO_NGS_METHYLSEQ_BISMARK = "bio.ngs.methylseq.bismark"

########################################
# bio.ngs.platform
########################################

##############################
# bio.ngs.platform.illumina
##############################
BIO_NGS_PLATFORM_ILLUMINA = "bio.ngs.platform.illumina"

##################################################
# bio.ngs.popgen
##################################################

##################################################
# bio.ngs.popgen.angsd
##################################################
BIO_NGS_POPGEN_ANGSD = "bio.ngs.popgen.angsd"
ANGSD = "angsd"
CHR_LABEL = "chr_label"
CHR_EXTRA = "chr_extra"
GEM = "gem"
MAPPABILITY_OPTIONS = "mappability_options"
NCHR_OPTION = "nchr_option"
POPULATIONS = "populations"

##################################################
# bio.ngs.qc.picard
##################################################
BIO_NGS_QC_PICARD = "bio.ngs.qc.picard"

PLOTMETRICS = "plotmetrics"

DUPMETRICS_TARGETS = "dupmetrics_targets"
ALIGNMETRICS_TARGETS = "alignmetrics_targets"
MERGE_SAM_TARGETS = "merge_sam_targets"
HSMETRICS_TARGETS = "hsmetrics_targets"
INSERTMETRICS_TARGETS = "insertmetrics_targets"
SORT_SAM_OPTIONS = "sort_sam_options"
MERGE_SAM_OPTIONS = "merge_sam_options"
MERGE_SAM_PREFIX = "merge_sam_prefix"
ADD_OR_REPLACE_READ_GROUPS_OPTIONS = "add_or_replace_read_groups_options"

# Jar program names
BUILD_BAM_INDEX = "BuildBamIndex.jar"
SORT_SAM = "SortSam.jar"
MERGE_SAM_FILES = "MergeSamFiles.jar"
REORDER_SAM = "ReorderSam.jar"
MARK_DUPLICATES = "MarkDuplicates.jar"
CREATE_SEQUENCE_DICTIONARY = "CreateSequenceDictionary.jar"
COLLECT_INSERT_SIZE_METRICS = "CollectInsertSizeMetrics.jar"
COLLECT_ALIGNMENT_SUMMARY_METRICS = "CollectAlignmentSummaryMetrics.jar"
CALCULATE_HS_METRICS = "CalculateHsMetrics.jar"
ADD_OR_REPLACE_READ_GROUPS = "AddOrReplaceReadGroups.jar"

##################################################
# bio.ngs.qc.sequenceprocessing
##################################################
BIO_NGS_QC_SEQUENCEPROCESSING = "bio.ngs.qc.sequenceprocessing"
TRIMMOMATIC_JAR_PROGRAM = "trimmomatic.jar"
THREEPRIME = "threeprime"
FIVEPRIME = "fiveprime"
CUTADAPT = "cutadapt"
TRIMMOMATIC = "trimmomatic"
RESYNCMATES = "resyncmates"
PROCESSING_OPTIONS = "processing_options"
FASTQC = "fastqc"
TRIM_GALORE = "trim_galore"

##################################################
# bio.ngs.tools.gatk
##################################################
BIO_NGS_TOOLS_GATK="bio.ngs.tools.gatk"

# GATK modules
GATK_JAR_PROGRAM = "GenomeAnalysisTK.jar"
UNIFIED_GENOTYPER = "UnifiedGenotyper"
PRINT_READS = "PrintReads"
BASE_RECALIBRATOR = "BaseRecalibrator"
INDEL_REALIGNER = "IndelRealigner"
REALIGNER_TARGET_CREATOR = "RealignerTargetCreator"
VARIANT_EVAL = "VariantEval"
VARIANT_FILTRATION = "VariantFiltration"
READ_BACKED_PHASING = "ReadBackedPhasing"
CLIP_READS = "ClipReads"
SELECT_SNP_VARIANTS = "SelectSnpVariants"

##################################################
# bio.ngs.tools.samtools
##################################################
BIO_NGS_TOOLS_SAMTOOLS = "bio.ngs.tools.samtools"
