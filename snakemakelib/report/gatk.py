# Copyright (c) 2014 Per Unneberg
import os
from snakemake.utils import R

def gatk_eval_report_plots_R(input, output):
    R("""
      library(gsalib)
      library(lattice)

      # lattice.options(trellis.par.set(simpleTheme(pch=19)))
      
      ev <- gsa.read.gatkreport('{input}')

      # Plot results from CompOverlap
      x <- ev$CompOverlap
      x <- x[x$Sample != "all", ]
      
      png('{output.variants_per_sample}')
      print(stripplot(nEvalVariants ~ Sample | Filter + Novelty, data=x,
      subset=(Filter != "raw" & Novelty != "all"),
      scales=list(x=list(draw=FALSE)),
      main="Number of variants per sample\n(variants are classified as known/novel and called/filtered)",
      ylab="Number of variants", xlab="Sample"
      ))
      dev.off()

      png('{output.known_site_freq}')
      print(stripplot(compRate ~ Sample | Filter, data=x,
      subset=(Novelty == "all"), index.cond=list(c(3,1,2)),
      scales=list(x=list(draw=FALSE)), layout=c(3,1),
      main="Frequency of known sites",
      ylab="Known sites (%)", xlab="Sample",
      ))
      dev.off()

      png('{output.dbsnp_concordance_known}')
      print(stripplot(concordantRate ~ Sample | Filter, data=x,
                subset=(Novelty == "known"), index.cond=list(c(3,1,2)),
                scales=list(x=list(draw=FALSE)), layout=c(3,1),
                main="Concordance with dbSNP for variants at known sites",
                ylab="Concordant variants (%)", xlab="Sample"
                ))
      dev.off()

      # Plot results from CountVariants
      x <- ev$CountVariants
      x <- x[x$Sample != "all", ]


      count.metrics <- list(nSNP = c("Number of SNPs", "nSNPs", '{output.nSNP}'),
      nIns = c("Number of insertions", "nInsertions", '{output.nIns}'),
      nDel = c("Number of deletions", "nDeletions", '{output.nDel}'),
      nComp = c("Number of complex variants", "nComplex", '{output.nComp}'),
      nMNP = c("Number of other variants (MNP, symbolic and mixed)", "I(nMNPs + nSymbolic + nMixed)", '{output.nMNP}'),
      nHets = c("Number of heterozygous sites", "nHets", '{output.nHets}'),
      nHomRef = c("Number of homozygous-reference sites", "nHomRef", '{output.nHomRef}'),
      nHomVar = c("Number of homozygous-alternative sites", "nHomVar", '{output.nHomVar}'),
      nNoCalls = c("Number of 'no calls'", "nNoCalls", '{output.nNoCalls}'),
      nSingletons = c("Number of private variants", "nSingletons", '{output.nSingletons}')
      )

      for(metric.name in names(count.metrics)) {{
      png(count.metrics[[metric.name]][3])
      print(stripplot(as.formula(paste(count.metrics[[metric.name]][2], "~ Sample | Filter + Novelty")),
      data=x, subset=(Filter != "raw" & Novelty != "all"),
      scales=list(x=list(draw=FALSE)),
      main=paste(metric.name, "per sample"), ylab=count.metrics[[metric.name]][1], xlab="Sample"
      ))
      dev.off()
      }}

      ## Plot results from TiTvVariantEvaluator

      x <- ev$TiTvVariantEvaluator
      x <- x[x$Sample != "all", ]
      
      png('{output.TiTv}')
      print(stripplot(tiTvRatio ~ Sample | Filter + Novelty, data=x,
      subset=(Filter != "raw" & Novelty != "all"),
      scales=list(x=list(draw=FALSE)),
      main="Transition / transversion ratio for each sample",
      ylab="Transition / transversion ratio", xlab="Sample"
      ))
      dev.off()
    """)
