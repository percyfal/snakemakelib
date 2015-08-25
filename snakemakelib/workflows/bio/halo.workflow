halo-sbatch:
	@echo Running halo-sbatch recipe: submit data to queue manager
	@echo All sample targets: $(SAMPLE_TARGETS)
	@echo All samples: $(SAMPLES)
	@$(eval NTARGETS=$(shell echo $(SAMPLES) | awk '{print NF}'))
	@i=0; sbatchcounter=0; SAMPLELIST=; JOBNAMES=""; for SMP in $(SAMPLES); do \
	let i=i+1; \
	SAMPLELIST="$$SAMPLELIST$$SMP "; \
	if [ $$(($$i % $(HALO_BATCH_SIZE))) == 0 ] || [ $$i == $(NTARGETS) ]; then  \
	let sbatchcounter=sbatchcounter+1; \
	$(MAKE) halo-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) samples SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=halo-$$sbatchcounter.out SLURM_ERR=halo-$$sbatchcounter.err ;\
	$(MAKE) halo-settings SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) samples SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=halo-$$sbatchcounter.out SLURM_ERR=halo-$$sbatchcounter.err > halo-$$sbatchcounter.settings; \
	$(MAKE) -pn halo-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) samples SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=halo-$$sbatchcounter.out SLURM_ERR=halo-$$sbatchcounter.err >> halo-$$sbatchcounter.settings; \
	JOBNAMES="$$JOBNAMES,halo-$$sbatchcounter"; \
	SAMPLELIST=; \
	fi \
	done ;
%.raw.vcf: %.bam
	$(GATK_COMMAND) -T UnifiedGenotyper $(GATK_UNIFIEDGENOTYPER_RAW_OPTIONS) -I $< -o $@.tmp && mv $@.tmp $@  && mv $@.tmp.idx $@.idx
%.intervals: %.bam %.raw.vcf
	$(GATK_COMMAND) -T RealignerTargetCreator $(GATK_REALIGN_TARGET_CREATOR_OPTIONS) -I $< -o $@.tmp && mv $@.tmp $@
%.merge.realign.bam: %.merge.bam %.merge.intervals %.merge.raw.vcf
	$(GATK_COMMAND) -T IndelRealigner $(GATK_INDELREALIGNER_OPTIONS) -I $< -known $(word 3, $^) -o $@.tmp --targetIntervals $(word 2, $^) && mv $@.tmp $@  && mv $@.tmp.bai $@.bai
.PRECIOUS: %.bai %$(READ1_LABEL).trimmed.sync.fastq.gz %$(READ2_LABEL).trimmed.sync.fastq.gz %.vcf %.raw.vcf %.filtered.vcf %.sort.merge.realign.recal.clip.bam
CLEANTARGETS:=$(wildcard $(SAMPLE_PREFIX)*/*/*trimmed*) $(wildcard $(SAMPLE_PREFIX)*/*/*bam*) $(wildcard $(SAMPLE_PREFIX)*/*/*metrics) $(wildcard $(SAMPLE_PREFIX)*/*merge*) all.*
clean:
	rm -f $(CLEANTARGETS)
.PHONY: trimsync flowcells samples %.log halo all halo-settings halo-header
halo-header:
	@echo -e "\nMakefile.halo options"
	@echo "====================="
halo-targets:
	@echo -e "\nhalo targets"
	@echo "---------------------"
halo-settings: halo-header print-HALO_TARGET_PREFIX print-HALO_BATCH_SIZE print-GATK_UNIFIEDGENOTYPER_RAW_OPTIONS halo-targets print-TRIMSYNC_TARGETS print-FLOWCELL_TARGETS print-SAMPLE_TARGETS general-settings ngsvars-settings sequenceprocessing-settings samtools-settings bwa-settings picard-settings gatk-settings
flowcells: $(FLOWCELL_TARGETS)
samples: $(FLOWCELL_TARGETS) $(VCF_TARGETS) $(SAMPLE_TARGETS)
$(HALO_TARGET_PREFIX).bam: $(VCF_TARGETS) $(SAMPLE_TARGETS)
	@$(eval NTARGETS=$(shell echo $(SAMPLE_TARGETS) | awk '{print NF}'))
	@if [ $(NTARGETS) != 0 ]; then \
	$(PICARD_JAVA) $(PICARD_HOME)/MergeSamFiles.jar $(addprefix INPUT=,$(SAMPLE_TARGETS)) O=$@.tmp $(PICARD_OPTIONS) $(PICARD_MERGESAM_OPTIONS) && mv $@.tmp $@ && mv $@.tmp.bai $(@:.bam=).bai; \
	fi
halo: all
all: $(FLOWCELL_TARGETS) $(HALO_TARGET_PREFIX).filtered.eval_metrics metrics.txt
%.trimmed.sync.bam: %$(READ1_LABEL).trimmed.sync.fastq.gz %$(READ2_LABEL).trimmed.sync.fastq.gz
	$(BWA) mem $(BWA_OPTIONS) $(BWA_REF) $^ | $(SAMTOOLS) view -Sbh - > $@.tmp && mv $@.tmp $@
%.sort.rg.bam: %.sort.bam
	java -Xmx2g -jar $(PICARD_HOME)/AddOrReplaceReadGroups.jar INPUT=$< OUTPUT=$@.tmp SORT_ORDER=coordinate \
	RGID=$(firstword $(subst ., ,$*)) RGLB=lib RGPL=ILLUMINA RGPU=$(firstword $(subst ., ,$*)) \
	RGSM=$(firstword $(subst /, ,$(firstword $(subst ., ,$*)))) CREATE_INDEX=true && mv $@.tmp $@; mv $@.tmp.bai $(@:.bam=).bai
%.sort.merge.bam: $(FLOWCELL_TARGETS)
	@$(eval INPUTFILES=$(addprefix INPUT=,$(filter $(dir $*)%, $(FLOWCELL_TARGETS))))
	$(PICARD_JAVA) $(PICARD_HOME)/MergeSamFiles.jar $(INPUTFILES) O=$@.tmp $(PICARD_COMMON_OPTIONS) $(PICARD_MERGESAM_OPTIONS) && mv $@.tmp $@ && mv $@.tmp.bai $(@:.bam=).bai

# -*- snakemake -*-
import os
from snakemakelib.utils import update_config, sml_rules_path

# Start by including the general snakefile
include: os.path.join(sml_rules_path(), 'base_settings.rules')

config_default = { 
	"halo" : {
		"GATK_VARIANT_EVAL_OPTIONS" : "-ST Filter -l INFO --doNotUseAllStandardModules --evalModule CompOverlap --evalModule CountVariants --evalModule TiTvVariantEvaluator --evalModule ValidationReport --stratificationModule Filter",
		"HALO_TARGET_PREFIX" : "halo",
		"GATK_VARIANTFILTRATION_OPTIONS" : "--clusterWindowSize 10 --clusterSize 3 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 10" --filterName "LowCoverage" --filterExpression "QUAL < 30.0" --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0" --filterName "LowQual" --filterExpression "QD < 1.5" --filterName "LowQD"",
		"GATK_CLIPREADS_OPTIONS" : "--cyclesToTrim 1-5 --clipRepresentation WRITE_NS -R $(GATK_REF)",
		"GATK_UNIFIEDGENOTYPER_RAW_OPTIONS" : "-stand_call_conf 30.0 -stand_emit_conf 10.0  --downsample_to_coverage 200 --output_mode EMIT_VARIANTS_ONLY -glm BOTH -nt $(GATK_THREADS) -R $(GATK_REF)",
		"HALO_BATCH_SIZE" : "8",
		"GATK_UNIFIEDGENOTYPER_OPTIONS" : "-stand_call_conf 30.0 -stand_emit_conf 10.0  --downsample_to_coverage 200 --output_mode EMIT_VARIANTS_ONLY -glm BOTH -nt $(GATK_THREADS) -R $(GATK_REF) --dbsnp $(GATK_DBSNP)",
		},
	},
}

config = update_config(config_default, config)

MERGE_INPUT_SUFFIX=.trimmed.sync.sort.rg.bam
SAMPLE_TARGET_SUFFIX=.sort.merge.realign.recal.clip.bam
VCF_TARGET_SUFFIX=.sort.merge.raw.vcf
FLOWCELL_TARGETS=$(subst $(READ1_LABEL).fastq.gz,$(MERGE_INPUT_SUFFIX),$(FASTQFILES))
SAMPLE_TARGETS=$(foreach s,$(SAMPLES),$(s)/$(s)$(SAMPLE_TARGET_SUFFIX))
VCF_TARGETS=$(foreach s,$(SAMPLES),$(s)/$(s)$(VCF_TARGET_SUFFIX))
PICARD_DUPMETRICS_TARGETS=$(realpath $(subst .recal.clip.bam,.dup_metrics,$(SAMPLE_TARGETS)))
PICARD_HSMETRICS_TARGETS=$(realpath $(subst .recal.clip.bam,.hs_metrics,$(SAMPLE_TARGETS)) $(subst .bam,.hs_metrics,$(FLOWCELL_TARGETS)))
PICARD_INSERTMETRICS_TARGETS=$(realpath $(subst .recal.clip.bam,.insert_metrics,$(SAMPLE_TARGETS)) $(subst .bam,.insert_metrics,$(FLOWCELL_TARGETS)))
PICARD_ALIGNMETRICS_TARGETS=$(subst .recal.clip.bam,.align_metrics,$(SAMPLE_TARGETS)) $(subst .bam,.align_metrics,$(FLOWCELL_TARGETS))
rule halo-sbatch:
	shell: "@echo Running halo-sbatch recipe: submit data to queue manager@echo All sample targets: $(SAMPLE_TARGETS)@echo All samples: $(SAMPLES)@$(eval NTARGETS=$(shell echo $(SAMPLES) | awk '{print NF}'))@i=0; sbatchcounter=0; SAMPLELIST=; JOBNAMES=""; for SMP in $(SAMPLES); do \let i=i+1; \SAMPLELIST="$$SAMPLELIST$$SMP "; \if [ $$(($$i % $(HALO_BATCH_SIZE))) == 0 ] || [ $$i == $(NTARGETS) ]; then  \let sbatchcounter=sbatchcounter+1; \$(MAKE) halo-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) samples SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=halo-$$sbatchcounter.out SLURM_ERR=halo-$$sbatchcounter.err ;\$(MAKE) halo-settings SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) samples SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=halo-$$sbatchcounter.out SLURM_ERR=halo-$$sbatchcounter.err > halo-$$sbatchcounter.settings; \$(MAKE) -pn halo-$$sbatchcounter.sbatch SLURM_COMMAND="$(MAKE) -j $(SLURM_MAKE_J) $(SLURM_MAKE_OPTIONS) samples SAMPLES=\'$$SAMPLELIST\'" SLURM_OUT=halo-$$sbatchcounter.out SLURM_ERR=halo-$$sbatchcounter.err >> halo-$$sbatchcounter.settings; \JOBNAMES="$$JOBNAMES,halo-$$sbatchcounter"; \SAMPLELIST=; \fi \done ;"
rule rule_2:
	input: " {prefix}.bam"
	output: "{prefix}.raw.vcf"
	shell: "$(GATK_COMMAND) -T UnifiedGenotyper $(GATK_UNIFIEDGENOTYPER_RAW_OPTIONS) -I $< -o $@.tmp && mv $@.tmp $@  && mv $@.tmp.idx $@.idx"
rule rule_3:
	input: " {prefix}.bam {prefix}.raw.vcf"
	output: "{prefix}.intervals"
	shell: "$(GATK_COMMAND) -T RealignerTargetCreator $(GATK_REALIGN_TARGET_CREATOR_OPTIONS) -I $< -o $@.tmp && mv $@.tmp $@"
rule rule_4:
	input: " {prefix}.merge.bam {prefix}.merge.intervals {prefix}.merge.raw.vcf"
	output: "{prefix}.merge.realign.bam"
	shell: "$(GATK_COMMAND) -T IndelRealigner $(GATK_INDELREALIGNER_OPTIONS) -I $< -known $(word 3, $^) -o $@.tmp --targetIntervals $(word 2, $^) && mv $@.tmp $@  && mv $@.tmp.bai $@.bai"
rule .PRECIOUS:
	input: " %.bai %$(READ1_LABEL).trimmed.sync.fastq.gz %$(READ2_LABEL).trimmed.sync.fastq.gz %.vcf %.raw.vcf %.filtered.vcf %.sort.merge.realign.recal.clip.bam"
rule clean:
	shell: "rm -f $(CLEANTARGETS)"
rule .PHONY:
	input: " trimsync flowcells samples %.log halo all halo-settings halo-header"
rule halo-header:
	shell: "@echo -e "\nMakefile.halo options"@echo "=====================""
rule halo-targets:
	shell: "@echo -e "\nhalo targets"@echo "---------------------""
rule halo-settings:
	input: " halo-header print-HALO_TARGET_PREFIX print-HALO_BATCH_SIZE print-GATK_UNIFIEDGENOTYPER_RAW_OPTIONS halo-targets print-TRIMSYNC_TARGETS print-FLOWCELL_TARGETS print-SAMPLE_TARGETS general-settings ngsvars-settings sequenceprocessing-settings samtools-settings bwa-settings picard-settings gatk-settings"
rule samples:
	input: " $(FLOWCELL_TARGETS) $(VCF_TARGETS) $(SAMPLE_TARGETS)"
rule $(HALO_TARGET_PREFIX).bam:
	input: " $(VCF_TARGETS) $(SAMPLE_TARGETS)"
	shell: "@$(eval NTARGETS=$(shell echo $(SAMPLE_TARGETS) | awk '{print NF}'))@if [ $(NTARGETS) != 0 ]; then \$(PICARD_JAVA) $(PICARD_HOME)/MergeSamFiles.jar $(addprefix INPUT=,$(SAMPLE_TARGETS)) O=$@.tmp $(PICARD_OPTIONS) $(PICARD_MERGESAM_OPTIONS) && mv $@.tmp $@ && mv $@.tmp.bai $(@:.bam=).bai; \fi"
rule halo:
	input: " all"
rule all:
	input: " $(FLOWCELL_TARGETS) $(HALO_TARGET_PREFIX).filtered.eval_metrics metrics.txt"
rule rule_15:
	input: " {prefix}$(READ1_LABEL).trimmed.sync.fastq.gz {prefix}$(READ2_LABEL).trimmed.sync.fastq.gz"
	output: "{prefix}.trimmed.sync.bam"
	shell: "$(BWA) mem $(BWA_OPTIONS) $(BWA_REF) $^ | $(SAMTOOLS) view -Sbh - > $@.tmp && mv $@.tmp $@"
rule rule_16:
	input: " {prefix}.sort.bam"
	output: "{prefix}.sort.rg.bam"
	shell: "java -Xmx2g -jar $(PICARD_HOME)/AddOrReplaceReadGroups.jar INPUT=$< OUTPUT=$@.tmp SORT_ORDER=coordinate \RGID=$(firstword $(subst ., ,$*)) RGLB=lib RGPL=ILLUMINA RGPU=$(firstword $(subst ., ,$*)) \RGSM=$(firstword $(subst /, ,$(firstword $(subst ., ,$*)))) CREATE_INDEX=true && mv $@.tmp $@; mv $@.tmp.bai $(@:.bam=).bai"
