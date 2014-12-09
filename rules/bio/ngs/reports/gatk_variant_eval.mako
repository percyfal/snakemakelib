GATK variant evaluation summary report
======================================

:Project: ${project_name}
:Application: ${application}
:Input: ${input}
% if region:
:Region: ${region}
% endif

Method overview
===============

This document contains plots of metrics that can be obtained by
running `GATK VariantEval
<https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php>`.
The output consists of a number of tables from which the plots have
been generated:

1. CompOverlap - comparison of called variants with a known set of
   variants, typically dbSNP. In addition, if the analysis is based on
   a sequence capture experiment, the variants have been stratified by
   the bait/target regions so as to assess the number of variants
   falling inside/outside the bait/target regions.

2. CountVariants - counts and classifies variants, e.g. insertions,
   deletions etc

3. TiTvVariantEvaluator - evaluation of statistics related to the
   number of transitions and transversions

4. ValidationReport - evaluation of specificity vs sensitivity.
   Currently no plots are reported based on this data.

% if region:

Region-based variant summary
-----------------------------

Table 1 shows a variant summary stratified for the region.

% else:

Sample-based variant summary
-----------------------------

Table 1 shows a variant summary stratified by sample. 

% end if

Results from CompOverlap
-------------------------

.. figure:: ${variants_per_sample}

   Number of variants per sample, stratified by dbSNP status (known/novel) and call status (called/filtered), where in the latter case called=variants that are kept for future analysis (PASS), and filtered=variants that fail filtering quality criteria.

.. figure:: ${known_site_freq}

   Frequency of known sites. This figure gives an overview of 1) the percentage of known sites calls for every sample and 2) how many of the known site calls are called versus filtered

.. figure:: ${dbsnp_concordance_known}

   dbSNP concordance for raw, called, and filtered calls.


Results from CountVariants
---------------------------

.. figure:: ${nSNP}

   The number of SNPs per sample, stratified by dbSNP status and call status

.. figure:: ${nIns}

   The number of insertions per sample, stratified by dbSNP status and call status

.. figure:: ${nDel}

   The number of deletions per sample, stratified by dbSNP status and call status

.. figure:: ${nComp}

   The number of complex variants per sample, stratified by dbSNP status and call status

.. figure:: ${nMNP}

   The number of multiple-nucleotide polymorphisms (MNPs) per sample,
   stratified by dbSNP status and call status

.. figure:: ${nHets}

   The number of heterozygotes per sample, stratified by dbSNP status and call status

.. figure:: ${nHomRef}

   The number of homozygous reference calls per sample, stratified by
   dbSNP status and call status

.. figure:: ${nHomVar}

   The number of homozygous variant calls per sample, stratified by
   dbSNP status and call status

.. figure:: ${nNoCalls}

   The number of non-calls per sample, stratified by dbSNP status and call status

.. figure:: ${nSingletons}

   The number of singletons per sample, stratified by dbSNP status and call status

Results from TiTvVariantEvaluator
-----------------------------------

.. figure:: ${TiTv}

   The transition/transversion ratio per sample, stratified by dbSNP
   status and call status


Attachments
------------
