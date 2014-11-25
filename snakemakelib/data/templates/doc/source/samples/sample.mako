<%!
import os
%>

${sample_id}
=======================================

Metrics
--------

Test

.. plot::

   import os
   import csv
   from pylab import *
   import matplotlib.pyplot as plt

   pmccsv = ${pmcmetrics.get('.align_metrics', None)}
   if not pmccsv is None:
       nseq = {'FIRST_OF_PAIR':[], 'SECOND_OF_PAIR':[], 'PAIR':[]}
       df = [row for row in csv.DictReader(pmccsv)]
       for row in df:
           nseq[row["CATEGORY"]].append(int(row["TOTAL_READS"]))
           
       n = len(${pmcidlist})
       xticks(range(0,n), [x for x in ${pmcidlist}], rotation=45)
       xlim(-.1, (n-1)*1.1)
       plt.plot(range(0,n), nseq['PAIR'], "o")
       plt.tight_layout()
       plt.show()


Sample runs
------------

:No. of reads:


