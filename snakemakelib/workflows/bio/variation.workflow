# -*- snakemake -*-
import os
from os.path import join
from snakemake.utils import update_config
from snakemakelib.config import SNAKEMAKELIB_RULES_PATH
from snakemakelib.bio.ngs.targets import generic_target_generator


def read_backed_phasing_create_input(wildcards):
    bamfile = wildcards.prefix.replace(".bp_variants", ".bam")
    return ["{prefix}.vcf".format(prefix=wildcards.prefix),
            bamfile, bamfile.replace(".bam", ".bai")]

variation_config = {
    'bio.ngs.tools.gatk': {
        'combine_variants': {
            'options': "--genotypemergeoption UNSORTED",
        },
        'read_backed_phasing': {
            'inputfun': read_backed_phasing_create_input,
        },
    },
}

update_config(variation_config, config)
config = variation_config

include: join(SNAKEMAKELIB_RULES_PATH, 'settings.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'utils.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'bio/ngs/variation', 'snpeff.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'bio/ngs/tools', 'gatk.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'bio/ngs/qc', 'picard.rules')
include: join(SNAKEMAKELIB_RULES_PATH, 'bio/ngs/align', 'bwa.rules')

ruleorder: gatk_print_reads > picard_build_bam_index

ruleorder: picard_build_bam_index > samtools_index
ruleorder: picard_add_or_replace_read_groups > picard_build_bam_index

path = config['bio.ngs.settings'].get('path') if not config['bio.ngs.settings'].get('path') is None else os.curdir

# Target suffices
TARGET_SUFFIX = ".sort.merge.rg.dup.realign.recal.bp_variants.phased.annotated.vcf"
DUP_METRICS_SUFFIX = ".sort.merge.rg.dup.dup_metrics"
ALIGN_METRICS_SUFFIX = ".sort.merge.rg.dup.align_metrics"
INSERT_METRICS_SUFFIX = ".sort.merge.rg.dup.insert_metrics"
HS_METRICS_SUFFIX = ".sort.merge.rg.dup.hs_metrics"

# Default targets
VCF_TARGETS = generic_target_generator(
    tgt_re=config['bio.ngs.settings']['sampleorg'].sample_re,
    target_suffix=TARGET_SUFFIX,
    src_re=config['bio.ngs.settings']['sampleorg'].raw_run_re,
    **config['bio.ngs.settings'])

VCF_TXT_TARGETS = generic_target_generator(
    tgt_re=config['bio.ngs.settings']['sampleorg'].sample_re,
    target_suffix=TARGET_SUFFIX.replace(".vcf", ".txt"),
    src_re=config['bio.ngs.settings']['sampleorg'].raw_run_re,
    **config['bio.ngs.settings'])

DUP_METRICS_TARGETS = generic_target_generator(
    tgt_re=config['bio.ngs.settings']['sampleorg'].sample_re,
    target_suffix=DUP_METRICS_SUFFIX,
    src_re=config['bio.ngs.settings']['sampleorg'].raw_run_re,
    **config['bio.ngs.settings'])

ALIGN_METRICS_TARGETS = generic_target_generator(
    tgt_re=config['bio.ngs.settings']['sampleorg'].sample_re,
    target_suffix=ALIGN_METRICS_SUFFIX,
    src_re=config['bio.ngs.settings']['sampleorg'].raw_run_re,
    **config['bio.ngs.settings'])

INSERT_METRICS_TARGETS = generic_target_generator(
    tgt_re=config['bio.ngs.settings']['sampleorg'].sample_re,
    target_suffix=INSERT_METRICS_SUFFIX,
    src_re=config['bio.ngs.settings']['sampleorg'].raw_run_re,
    **config['bio.ngs.settings'])

HS_METRICS_TARGETS = generic_target_generator(
    tgt_re=config['bio.ngs.settings']['sampleorg'].sample_re,
    target_suffix=HS_METRICS_SUFFIX,
    src_re=config['bio.ngs.settings']['sampleorg'].raw_run_re,
    **config['bio.ngs.settings'])

rule variation_all:
    input: VCF_TARGETS + VCF_TXT_TARGETS + DUP_METRICS_TARGETS + ALIGN_METRICS_TARGETS + INSERT_METRICS_TARGETS + HS_METRICS_TARGETS

# Run metrics only
rule variation_metrics:
    input: DUP_METRICS_TARGETS + ALIGN_METRICS_TARGETS + INSERT_METRICS_TARGETS + HS_METRICS_TARGETS

rule variation_snp_filtration:
    """Run variant filtration and variant recalibration


    The rule does snp variant filtration/variant recalibration. It
    sets up different filtration tasks depending on the *cov_interval*
    setting, which can be one of "regional", "exome", or None. The
    effects of the different choices are as follows:

    regional
      Do filtering based on JEXL-expressions. See `section 3, subtitle
      Recommendations for very small data sets
      <http://www.broadinstitute.org/gatk/guide/topic?name=best-practices>`_

    exome
      Use VQSR with modified argument settings (`--maxGaussians 4
      --percentBad 0.05`) as recommended in `3. Notes about small
      whole exome projects
      <http://www.broadinstitute.org/gatk/guide/topic?name=best-practices>`_

    None
      Perform "standard" VQSR

    """
    input: "{prefix}.snp.vcf"
    output: "{prefix}.snp.filtSNP.vcf"
    run:
        def _regional_JEXL_filter():
            cmd = config['bio.ngs.tools.gatk']['cmd'] + " -T VariantFiltration"
            options = " ".join([
                " ".join(["-R", config['bio.ngs.tools.gatk']['variant_snp_JEXL_filtration']['ref']]),
                " ".join(["--filterName GATKStandard{e} --filterExpression '{exp}'".format(e=exp.split()[0], exp=exp)
                          for exp in config['bio.ngs.tools.gatk']['variant_snp_JEXL_filtration']['expressions']])
            ])
            shell("{cmd} {opts} --variant {input} --out {out}".format(cmd=cmd, opts=options, input=input, out=output))

        def _exome_VQSR_filter():
            print("exome VQSR")

        def _standard_VQSR_filter():
            print("VQSR")

        if config['bio.ngs.tools.gatk']['cov_interval'] == "regional":
            _regional_JEXL_filter()
        elif config['bio.ngs.tools.gatk']['cov_interval'] == "exome":
            try:
                _exome_VQSR_filter()
            except:
                _regional_JEXL_filter()
        else:
            try:
                _standard_VQSR_filter()
            except:
                _regional_JEXL_filter()

rule variation_indel_filtration:
    """Run variant filtration and variant recalibration

    """
    input: "{prefix}.indel.vcf"
    output: "{prefix}.indel.filtINDEL.vcf"
    run:
        def _regional_JEXL_filter():
            cmd = config['bio.ngs.tools.gatk']['cmd'] + " -T VariantFiltration"
            options = " ".join([
                " ".join(["-R", config['bio.ngs.tools.gatk']['variant_indel_JEXL_filtration']['ref']]),
                " ".join(["--filterName GATKStandard{e} --filterExpression '{exp}'".format(e=exp.split()[0], exp=exp)
                          for exp in config['bio.ngs.tools.gatk']['variant_indel_JEXL_filtration']['expressions']])
            ])
            shell("{cmd} {opts} --variant {input} --out {out}".format(cmd=cmd, opts=options, input=input, out=output))

        def _exome_VQSR_filter():
            print("exome VQSR")

        def _standard_VQSR_filter():
            print("VQSR")

        if config['bio.ngs.tools.gatk']['cov_interval'] == "regional":
            _regional_JEXL_filter()
        elif config['bio.ngs.tools.gatk']['cov_interval'] == "exome":
            try:
                _exome_VQSR_filter()
            except:
                _regional_JEXL_filter()
        else:
            try:
                _standard_VQSR_filter()
            except:
                _regional_JEXL_filter()


rule variation_combine_variants:
    """Run GATK CombineVariants to combine variant files.

    The default rule combines files with suffixes filtSNP.vcf and
    filtINDEL.vcf.

    """
    params: cmd = config['bio.ngs.tools.gatk']['cmd'] + " -T " + config['bio.ngs.tools.gatk']['combine_variants']['cmd'],
            options = " ".join(["-R", config['bio.ngs.tools.gatk']['combine_variants']['ref'],
                                config['bio.ngs.tools.gatk']['combine_variants']['options']])

    input: "{prefix}.snp.filtSNP.vcf", "{prefix}.indel.filtINDEL.vcf"
    output: "{prefix}.bp_variants.vcf"
    run:
        inputstr = " ".join(["-V {}".format(x) for x in input])
        shell("{cmd} {ips} -o {out} {opt}".format(cmd=params.cmd,
                                                  ips=inputstr,
                                                  out=output,
                                                  opt=params.options))

rule clean:
    """Clean working directory. WARNING: will remove all files except
    (.fastq|.fastq.gz) and csv files
    """
    params: d = workflow._workdir
    shell: 'for f in `find  {params.d} -type f -name "*" | grep -v ".fastq$" | grep -v ".fastq.gz$" | grep -v ".csv$"`; do echo removing $f; rm -f $f; done'
