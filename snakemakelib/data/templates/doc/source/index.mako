Project summary
=============================

:Project: ${project_name}
:Application: ${application}
:Date: ${date}

Samples
--------

.. toctree::
   :maxdepth: 2

   samples/index

QC Metrics
----------

Sequence statistics
^^^^^^^^^^^^^^^^^^^

.. plot::

   import os
   import csv
   from pylab import *
   import matplotlib.pyplot as plt
   
   pmccsv = [m[2]['.align_metrics'] for m in sorted(metrics.values())]
   for c in pmccsv:
       df = [row for row in csv.DictReader(c)]
       for row in df:
            nseq[row["CATEGORY"]].append(int(row["TOTAL_READS"]))
       pct_aligned[row["CATEGORY"]].append(100 * float(row["PCT_PF_READS_ALIGNED"]))
       n = len(metrics.keys())

   pmccsv = [m[2]['.dup_metrics'] for m in sorted(metrics.values())]
   dup = []
   for c in pmccsv:
       df = [row for row in csv.DictReader(c)]
       for row in df:
           dup.append(100 * float(df[0]["PERCENT_DUPLICATION"]))

   if len(dup) == 0:
       sdup = [100 for i in range(0, n)]
   else: 
       sdup = [int(100 + 10 * x) for x in dup]
   plt.scatter(pct_aligned['PAIR'], nseq['PAIR'], s=sdup, alpha=0.75)
   plt.xlabel(r'Percent aligned', fontsize=14)
   plt.yscale('log', **{'basey':10})
   plt.ylabel(r'Read count', fontsize=14)
   plt.title("Sequence summary.\nPoint sizes correspond to duplication levels.", fontsize=14)
   plt.tight_layout()
   plt.show()

