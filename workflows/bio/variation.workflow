# -*- snakemake -*-
import os
from snakemakelib.config import update_sml_config, sml_rules_path, get_sml_config, init_sml_config

def read_backed_phasing_create_input(wildcards):
    bamfile = wildcards.prefix.replace(".bp_variants", ".bam")
    return ["{prefix}.vcf".format(prefix=wildcards.prefix), bamfile, bamfile.replace(".bam", ".bai")]

variation_config = {
    'bio.ngs.tools.gatk' : {
        'combine_variants' : {
            'options' : "--genotypemergeoption UNSORTED",
        },
        'read_backed_phasing' : {
            'inputfun' : read_backed_phasing_create_input,
        },
    },
}

update_sml_config(variation_config)

include: os.path.join(sml_rules_path(), 'settings.rules')
include: os.path.join(sml_rules_path(), 'utils.rules')
include: os.path.join(sml_rules_path(), 'bio/ngs/variation', 'variation.rules')
include: os.path.join(sml_rules_path(), 'bio/ngs/tools', 'gatk.rules')
include: os.path.join(sml_rules_path(), 'bio/ngs/qc', 'picard.rules')
include: os.path.join(sml_rules_path(), 'bio/ngs/align', 'bwa.rules')

cfg = get_sml_config()

ruleorder: gatk_print_reads > picard_build_bam_index

# Default targets: expand samples and flowcells
TARGET_SUFFIX=".sort.merge.rg.dup.realign.recal.bp_variants.phased.annotated.vcf"

rule all:
    input: ["{sample}{sep}{sample}{sfx}".format(sep=os.sep, sample=x, sfx=TARGET_SUFFIX) for x in cfg['bio.ngs.settings']['samples']] + ["{sample}{sep}{sample}{sfx}".format(sep=os.sep, sample=x, sfx=TARGET_SUFFIX).replace(".vcf", ".txt") for x in cfg['bio.ngs.settings']['samples']]

rule variation_snp_filtration:
    """Run variant filtration and variant recalibration


    The FiltrationWrapper wraps snp and indel variant filtration
    tasks. It sets up different filtration tasks depending on the
    *cov_interval* setting, which can be one of "regional", "exome",
    or None. The effects of the different choices are as follows:

    regional
      Do filtering based on JEXL-expressions. See `section 3, subtitle Recommendations for very small data sets <http://www.broadinstitute.org/gatk/guide/topic?name=best-practices>`_

    exome
      Use VQSR with modified argument settings (`--maxGaussians 4 --percentBad 0.05`) as recommended in  `3. Notes about small whole exome projects <http://www.broadinstitute.org/gatk/guide/topic?name=best-practices>`_

    None
      Perform "standard" VQSR
    """
    input: "{prefix}.snp.vcf"
    output: "{prefix}.snp.filtSNP.vcf"
    run:
          def _regional_JEXL_filter():
              cmd = cfg['bio.ngs.tools.gatk']['cmd'] + " -T VariantFiltration"
              options = " ".join([
                  " ".join(["-R", cfg['bio.ngs.tools.gatk']['variant_snp_JEXL_filtration']['ref']]),
                  " ".join(["--filterName GATKStandard{e} --filterExpression '{exp}'".format(e=exp.split()[0], exp=exp) \
                            for exp in cfg['bio.ngs.tools.gatk']['variant_snp_JEXL_filtration']['expressions']])
              ])
              shell("{cmd} {opts} --variant {input} --out {out}".format(cmd=cmd, opts=options, input=input, out=output))
          def _exome_VQSR_filter():
              print ("exome VQSR")
          def _standard_VQSR_filter():
              print ("VQSR")
          
          if cfg['bio.ngs.tools.gatk']['cov_interval'] == "regional":
              _regional_JEXL_filter()
          elif cfg['bio.ngs.tools.gatk']['cov_interval'] == "exome":
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
              cmd = cfg['bio.ngs.tools.gatk']['cmd'] + " -T VariantFiltration"
              options = " ".join([
                  " ".join(["-R", cfg['bio.ngs.tools.gatk']['variant_indel_JEXL_filtration']['ref']]),
                  " ".join(["--filterName GATKStandard{e} --filterExpression '{exp}'".format(e=exp.split()[0], exp=exp) \
                            for exp in cfg['bio.ngs.tools.gatk']['variant_indel_JEXL_filtration']['expressions']])
              ])
              shell("{cmd} {opts} --variant {input} --out {out}".format(cmd=cmd, opts=options, input=input, out=output))
          def _exome_VQSR_filter():
              print ("exome VQSR")
          def _standard_VQSR_filter():
              print ("VQSR")
          
          if cfg['bio.ngs.tools.gatk']['cov_interval'] == "regional":
              _regional_JEXL_filter()
          elif cfg['bio.ngs.tools.gatk']['cov_interval'] == "exome":
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
    params: cmd = cfg['bio.ngs.tools.gatk']['cmd'] + " -T " + cfg['bio.ngs.tools.gatk']['combine_variants']['cmd'],
            options = " ".join(["-R", cfg['bio.ngs.tools.gatk']['combine_variants']['ref'],
                                cfg['bio.ngs.tools.gatk']['combine_variants']['options']])

    input: "{prefix}.snp.filtSNP.vcf", "{prefix}.indel.filtINDEL.vcf"
    output: "{prefix}.bp_variants.vcf"
    run: 
        inputstr = " ".join(["-V {}".format(x) for x in input])
        shell("{cmd} {ips} -o {out} {opt}".format(cmd=params.cmd, ips=inputstr, out=output, opt=params.options))

rule clean:
    """Clean working directory. WARNING: will remove all files except
    (.fastq|.fastq.gz) and csv files
    """
    params: d = workflow._workdir
    shell: 'for f in `find  {params.d} -type f -name "*" |grep -v ".snakemake" | grep -v ".fastq$" | grep -v ".fastq.gz$" | grep -v ".csv$"`; do echo removing $f; rm -f $f; done'
