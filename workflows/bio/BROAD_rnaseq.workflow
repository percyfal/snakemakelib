rnaseq-sbatch:
	@echo Running rnaseq-sbatch recipe: submit data to queue manager
	@echo All sample targets: $(TARGETS)
	@echo All samples: $(SAMPLES)
	@$(eval NTARGETS=$(shell echo $(SAMPLES) | awk '{print NF}'))
	@i=0; sbatchcounter=0; SAMPLELIST=; JOBNAMES=""; for SMP in $(SAMPLES); do \
	let i=i+1; \
	SAMPLELIST="$$SAMPLELIST$$SMP "; \
	if [ $$(($$i % $(RNASEQ_BATCH_SIZE))) == 0 ] || [ $$i == $(NTARGETS) ]; then  \
	let sbatchcounter=sbatchcounter+1; \
	$(MAKE) rnaseq-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) $(TARGET_NAME) SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=rnaseq-$$sbatchcounter.out SLURM_ERR=rnaseq-$$sbatchcounter.err ;\
	$(MAKE) BROAD_rnaseq-settings SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) $(TARGET_NAME) SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=rnaseq-$$sbatchcounter.out SLURM_ERR=rnaseq-$$sbatchcounter.err > BROAD_rnaseq-$$sbatchcounter.settings; \
	$(MAKE) -pn rnaseq-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) $(TARGET_NAME) SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=rnaseq-$$sbatchcounter.out SLURM_ERR=rnaseq-$$sbatchcounter.err >> BROAD_rnaseq-$$sbatchcounter.settings; \
	JOBNAMES="$$JOBNAMES,rnaseq-$$sbatchcounter"; \
	SAMPLELIST=; \
	fi \
	done ;
clean:
	rm -f $(CLEANTARGETS)
.PHONY: %.log rnaseq all BROAD_rnaseq-settings BROAD_rnaseq-header
BROAD_rnaseq-header:
	@echo -e "\nMakefile.BROAD_rnaseq options"
	@echo "====================="
BROAD_rnaseq-targets:
	@echo -e "\nrnaseq targets"
	@echo "---------------------"
BROAD_rnaseq-settings: BROAD_rnaseq-header print-RNASEQ_BATCH_SIZE
.PRECIOUS: %$(READ1_LABEL)$(FASTQ_SUFFIX).gz %$(READ2_LABEL)$(FASTQ_SUFFIX).gz %.rsem %.tophat2 %.cufflinks %.cufflinks_quant %/$(TOPHAT2_OPTION_OUTPUT_DIR)/accepted_hits.bam
rsem: $(RSEM_TARGETS)
tophat2: $(TOPHAT2_TARGETS)
cufflinks: $(CUFFLINKS_TARGETS)
cufflinks_quant: $(CUFFLINKS_QUANT_TARGETS)
%/$(TOPHAT2_OPTION_OUTPUT_DIR)/accepted_hits.bam: $(TOPHAT2_TARGETS)
	echo "Checkpoint that all tophat2 targets have been completed before running rnaseqc"
rnaseqc_targets: $(RNASEQC_TARGETS)
all: rsem tophat2 cufflinks cufflinks_quant rnaseqc_targets

# -*- snakemake -*-
import os
from snakemakelib.utils import update_config, sml_rules_path

# Start by including the general snakefile
include: os.path.join(sml_rules_path(), 'base_settings.rules')

config_default = { 
	"BROAD_rnaseq" : {
		"RNASEQ_BATCH_SIZE" : "8",
		},
	},
}

config = update_config(config, config_default)

FASTQ_SUFFIX=.P.qtrim.fq
TOPHAT2_OPTION_OUTPUT_DIR=TOPHAT
CUFFLINKS_OPTION_OUTPUT_DIR=CUFFLINKS
RSEM_TARGETS=$(foreach s,$(SAMPLES),$(s)/$(s).rsem)
TOPHAT2_TARGETS=$(foreach s,$(SAMPLES),$(s)/$(s).tophat2)
CUFFLINKS_TARGETS=$(foreach s,$(SAMPLES),$(s)/$(s).cufflinks)
CUFFLINKS_QUANT_TARGETS=$(foreach s,$(SAMPLES),$(s)/$(s).cufflinks_quant)
RNASEQC_TARGETS=$(foreach s,$(SAMPLES),$(s)/$(TOPHAT2_OPTION_OUTPUT_DIR)/accepted_hits.rg.resorted.dup.rnaseqc)
TARGETS=$(RSEM_TARGETS) $(TOPHAT2_TARGETS) $(CUFFLINKS_TARGETS) $(CUFFLINKS_QUANT_TARGETS) $(RNASEQC_TARGETS)
TARGET_NAME=all
CLEANTARGETS=$(wildcard *.fifo) $(wildcard */*.fifo)
rule rnaseq-sbatch:
	shell: "@echo Running rnaseq-sbatch recipe: submit data to queue manager@echo All sample targets: $(TARGETS)@echo All samples: $(SAMPLES)@$(eval NTARGETS=$(shell echo $(SAMPLES) | awk '{print NF}'))@i=0; sbatchcounter=0; SAMPLELIST=; JOBNAMES=""; for SMP in $(SAMPLES); do \let i=i+1; \SAMPLELIST="$$SAMPLELIST$$SMP "; \if [ $$(($$i % $(RNASEQ_BATCH_SIZE))) == 0 ] || [ $$i == $(NTARGETS) ]; then  \let sbatchcounter=sbatchcounter+1; \$(MAKE) rnaseq-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) $(TARGET_NAME) SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=rnaseq-$$sbatchcounter.out SLURM_ERR=rnaseq-$$sbatchcounter.err ;\$(MAKE) BROAD_rnaseq-settings SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) $(TARGET_NAME) SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=rnaseq-$$sbatchcounter.out SLURM_ERR=rnaseq-$$sbatchcounter.err > BROAD_rnaseq-$$sbatchcounter.settings; \$(MAKE) -pn rnaseq-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) $(TARGET_NAME) SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=rnaseq-$$sbatchcounter.out SLURM_ERR=rnaseq-$$sbatchcounter.err >> BROAD_rnaseq-$$sbatchcounter.settings; \JOBNAMES="$$JOBNAMES,rnaseq-$$sbatchcounter"; \SAMPLELIST=; \fi \done ;"
rule clean:
	shell: "rm -f $(CLEANTARGETS)"
rule .PHONY:
	input: " %.log rnaseq all BROAD_rnaseq-settings BROAD_rnaseq-header"
rule BROAD_rnaseq-header:
	shell: "@echo -e "\nMakefile.BROAD_rnaseq options"@echo "=====================""
rule BROAD_rnaseq-targets:
	shell: "@echo -e "\nrnaseq targets"@echo "---------------------""
rule BROAD_rnaseq-settings:
	input: " BROAD_rnaseq-header print-RNASEQ_BATCH_SIZE"
rule .PRECIOUS:
	input: " %$(READ1_LABEL)$(FASTQ_SUFFIX).gz %$(READ2_LABEL)$(FASTQ_SUFFIX).gz %.rsem %.tophat2 %.cufflinks %.cufflinks_quant %/$(TOPHAT2_OPTION_OUTPUT_DIR)/accepted_hits.bam"
rule rsem:
	input: " $(RSEM_TARGETS)"
rule tophat2:
	input: " $(TOPHAT2_TARGETS)"
rule cufflinks:
	input: " $(CUFFLINKS_TARGETS)"
rule cufflinks_quant:
	input: " $(CUFFLINKS_QUANT_TARGETS)"
rule rule_12:
	input: " $(TOPHAT2_TARGETS)"
	output: "{prefix}/$(TOPHAT2_OPTION_OUTPUT_DIR)/accepted_hits.bam"
	shell: "echo "Checkpoint that all tophat2 targets have been completed before running rnaseqc""
rule rnaseqc_targets:
	input: " $(RNASEQC_TARGETS)"
rule all:
	input: " rsem tophat2 cufflinks cufflinks_quant rnaseqc_targets"
