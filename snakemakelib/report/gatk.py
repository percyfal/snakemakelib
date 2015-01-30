# Copyright (c) 2014 Per Unneberg
import os
from snakemake.utils import R

def gatk_eval_report_plots_R(input, output):
    R("""
      library(gsalib)
      library(lattice)
      library(RColorBrewer)
      cpal <- colorRampPalette(brewer.pal(9,"Paired"))(10)

      # lattice.options(trellis.par.set(simpleTheme(pch=19)))
      
      ev <- gsa.read.gatkreport('{input.metrics}')

      # Plot results from CompOverlap
      x <- ev$CompOverlap
      xs <- x

      # Check if IntervalStratification present; if so, add conditioning to plots
      stratstr <- ""
      if ("IntervalStratification" %in% colnames(x)) {{
          stratstr <- " + IntervalStratification"
          xs <- x[x$IntervalStratification == "all",]
      }}

      png('{output.variants_per_sample}')
      print(stripplot(as.formula(paste("nEvalVariants ~ Sample | Filter + Novelty")), data=xs,
      subset=(Filter != "raw" & Novelty != "all" & Sample != "all"),
      scales=list(x=list(draw=FALSE)),
      main="Number of variants per sample\n(variants are classified as known/novel and called/filtered)",
      ylab="Number of variants", xlab="Sample",
      par.settings=simpleTheme(pch=19, col=cpal), auto.key=TRUE
      ))
      dev.off()

      png('{output.variants_in_regions}')
      i <- c(3,1,2)
      if (stratstr != "") 
         i <- c(3,1,2,6,4,5,9,7,8)
      print(stripplot(as.formula(paste("nEvalVariants ~ Novelty | Filter", stratstr)), 
            data=x, subset = (Sample=="all"),
                      main="Number of variants by novelty stratified by filter and interval", par.settings=simpleTheme(pch=19, col=cpal), auto.key=TRUE))
      dev.off()

      png('{output.known_site_freq}')
      print(stripplot(as.formula(paste("compRate ~ Sample | Filter")), data=xs,
      subset=(Novelty == "all"), index.cond=list(c(3,1,2)),
      scales=list(x=list(draw=FALSE)), layout=c(3,1),
      main="Frequency of known sites",
      ylab="Known sites (%)", xlab="Sample", par.settings=simpleTheme(pch=19, col=cpal)
      ))
      dev.off()

      png('{output.dbsnp_concordance_known}')
      print(stripplot(as.formula(paste("concordantRate ~ Sample | Filter")), data=xs,
                subset=(Novelty == "known"), index.cond=list(c(3,1,2)),
                scales=list(x=list(draw=FALSE)), layout=c(3,1),
                main="Concordance with dbSNP for variants at known sites",
                ylab="Concordant variants (%)", xlab="Sample",
                par.settings=simpleTheme(pch=19, col=cpal)
                ))
      dev.off()
   
      # Plot results from CountVariants
      x <- ev$CountVariants
      x <- x[x$Sample == "all", ]
      xs <- x
      if ("IntervalStratification" %in% colnames(xs)) {{
          xs <- x[x$IntervalStratification == "all",]
      }}

      png('{output.nVariants}')
      print(stripplot(as.formula(paste("nVariantLoci + nSNPs + nInsertions + nDeletions ~ Novelty | Filter", stratstr)),
      data=xs, jitter.data=TRUE,
      main=paste("Total number of variants"), par.settings=simpleTheme(pch=19, col=cpal), auto.key=TRUE
      ))
      dev.off()

      x <- ev$CountVariants
      x <- x[x$Sample != "all", ]
      xs <- x
      if ("IntervalStratification" %in% colnames(xs)) {{
          xs <- x[x$IntervalStratification == "all",]
      }}

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
      print(stripplot(as.formula(paste(count.metrics[[metric.name]][2], "~ Sample | Novelty + Filter")),
      data=xs, subset=(Filter != "raw" & Novelty != "all"),
      scales=list(x=list(draw=FALSE)),
      main=paste(metric.name, "per sample"), ylab=count.metrics[[metric.name]][1], xlab="Sample", par.settings=simpleTheme(pch=19, col=cpal)
      ))
      dev.off()
      }}

      ## Plot results from TiTvVariantEvaluator

      x <- ev$TiTvVariantEvaluator
      x <- x[x$Sample != "all", ]
      xs <- x
      if ("IntervalStratification" %in% colnames(xs)) {{
          xs <- x[x$IntervalStratification == "all",]
      }}

      png('{output.TiTv}')
      print(stripplot(as.formula(paste("tiTvRatio ~ Sample | Novelty + Filter")), data=xs,
      subset=(Filter != "raw" & Novelty != "all"),
      scales=list(x=list(draw=FALSE)),
      main="Transition / transversion ratio for each sample",
      ylab="Transition / transversion ratio", xlab="Sample",
      par.settings=simpleTheme(pch=19, col=cpal)
      ))
      dev.off()

      # Finally make the output summary table
      if (stratstr != "") {{
         cv <- ev$CountVariants[ev$CountVariants$Sample!="all" & ev$CountVariants$Filter=="called" & ev$CountVariants$IntervalStratification=="all",]
         titv <- ev$TiTvVariantEvaluator[ev$TiTvVariantEvaluator$Sample!="all" & ev$TiTvVariantEvaluator$Filter=="called" & ev$TiTvVariantEvaluator$IntervalStratification=="all",]
    
         df.out <- cbind(cv[cv$Novelty=="all", c("Sample", "nVariantLoci")],
                         cv[cv$Novelty=="known", "nVariantLoci"],
                         titv[titv$Novelty=="all", "tiTvRatio"],
                         titv[titv$Novelty=="known", "tiTvRatio"],
                         titv[titv$Novelty=="novel", "tiTvRatio"])
         colnames(df.out) <- c("Sample", "Total", "Known", "Ti/Tv", "Ti/Tv known", "Ti/Tv novel")
         levels(df.out$Sample) <- c(levels(df.out$Sample), "all")

         cv.all <- ev$CountVariants[ev$CountVariants$Sample=="all" & ev$CountVariants$Filter=="called" & ev$CountVariants$IntervalStratification=="all",]
         titv.all <- ev$TiTvVariantEvaluator[ev$TiTvVariantEvaluator$Sample=="all" & ev$TiTvVariantEvaluator$Filter=="called" & ev$TiTvVariantEvaluator$IntervalStratification=="all",]
         df.all <- cbind(cv.all[cv.all$Novelty=="all", c("Sample", "nVariantLoci")],
                         cv.all[cv.all$Novelty=="known", "nVariantLoci"],
                         titv.all[titv.all$Novelty=="all", "tiTvRatio"],
                         titv.all[titv.all$Novelty=="known", "tiTvRatio"],
                         titv.all[titv.all$Novelty=="novel", "tiTvRatio"])
         colnames(df.all) <- c("Sample", "Total", "Known", "Ti/Tv", "Ti/Tv known", "Ti/Tv novel")
         df.out <- rbind(df.out, df.all)
         write.csv(format(df.out, digits=3), file='{output.varianttable}', row.names=FALSE)
      }} else {{
         cv <- ev$CountVariants[ev$CountVariants$Sample!="all" & ev$CountVariants$Filter=="called",]
         titv <- ev$TiTvVariantEvaluator[ev$TiTvVariantEvaluator$Sample!="all" & ev$TiTvVariantEvaluator$Filter=="called",]
    
         df.out <- cbind(cv[cv$Novelty=="all", c("Sample", "nVariantLoci")],
                         cv[cv$Novelty=="known", "nVariantLoci"],
                         titv[titv$Novelty=="all", "tiTvRatio"],
                         titv[titv$Novelty=="known", "tiTvRatio"],
                         titv[titv$Novelty=="novel", "tiTvRatio"])
         colnames(df.out) <- c("Sample", "Total", "Known", "Ti/Tv", "Ti/Tv known", "Ti/Tv novel")
         levels(df.out$Sample) <- c(levels(df.out$Sample), "all")

         cv.all <- ev$CountVariants[ev$CountVariants$Sample=="all" & ev$CountVariants$Filter=="called",]
         titv.all <- ev$TiTvVariantEvaluator[ev$TiTvVariantEvaluator$Sample=="all" & ev$TiTvVariantEvaluator$Filter=="called",]
         df.all <- cbind(cv.all[cv.all$Novelty=="all", c("Sample", "nVariantLoci")],
                         cv.all[cv.all$Novelty=="known", "nVariantLoci"],
                         titv.all[titv.all$Novelty=="all", "tiTvRatio"],
                         titv.all[titv.all$Novelty=="known", "tiTvRatio"],
                         titv.all[titv.all$Novelty=="novel", "tiTvRatio"])
         colnames(df.all) <- c("Sample", "Total", "Known", "Ti/Tv", "Ti/Tv known", "Ti/Tv novel")
         df.out <- rbind(df.out, df.all)
         write.csv(format(df.out, digits=3), file='{output.varianttable}', row.names=FALSE)
      }}
    """)
