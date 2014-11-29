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

   from pylab import *
   import matplotlib.pyplot as plt
   
   plt.scatter(${pct_aligned}, ${nseq}, s=${sdup}, alpha=0.75)
   plt.xlabel(r'Percent aligned', fontsize=14)
   plt.yscale('log', **{'basey':10})
   plt.ylabel(r'Read count', fontsize=14)
   plt.title("Sequence summary.\nPoint sizes correspond to duplication levels.", fontsize=14)
   plt.tight_layout()
   plt.show()

Alignment metrics
^^^^^^^^^^^^^^^^^

.. plot::

   from pylab import *
   import matplotlib.pyplot as plt
   ids = ${idlist}
   n = len(ids)
   xticks(range(0,n), [x for x in ids], rotation=45)
   xlim(-.1, (n-1)*1.1)
   plt.plot(range(0,n), ${pct_aligned}, "o")
   plt.tight_layout()
   plt.show()


Duplication metrics
^^^^^^^^^^^^^^^^^^^^

.. plot::

   from pylab import *
   import matplotlib.pyplot as plt
   
   ids = ${idlist}
   n = len(ids)
   xticks(range(0,n), [x for x in ids], rotation=45)
   xlim(-.1, (n-1)*1.1)
   plt.plot(range(0,n), ${sdup}, "o")
   plt.tight_layout()
   plt.show()


Hybridization metrics
^^^^^^^^^^^^^^^^^^^^^

.. plot::

   from pylab import *
   import matplotlib.pyplot as plt

   columns = ${columns}
   hticks = ${hticks}
   xticks(range(0,len(hticks)), [x for x in hticks])
   hsmetrics = ${hsmetrics}
   plt.boxplot(np.array(hsmetrics))
   plt.show()

.. plot::

   import math
   from pylab import *
   import matplotlib.pyplot as plt

   columns = ${columns}
   hticks = ${hticks}
   ids = ${idlist}
   n = len(ids)
   hsmetrics = ${hsmetrics}
   nsubplots = int(math.ceil(n/9))
   nrow = int(math.ceil(n/3))
   k = 0
   for i_subplot in range(0, nsubplots + 1):
      f, axarr = plt.subplots(3, 3, sharex='col', sharey='row')
      for i in range(0, 3):
      	  for j in range(0, 3):
       	      if k < n:
	      	  x = range(0, len(hticks))
               	  axarr[i,j].plot(x, hsmetrics[k], "o")
	       	  axarr[i,j].set_xticks(x)
	       	  axarr[i,j].set_title(ids[k])
	       	  axarr[i,j].set_xlim(-.1, (len(hticks)-1)*1.1)
	       	  axarr[i,j].set_ylim(-5, 105)
	       	  axarr[i,j].set_xticklabels(hticks)
   	      else:
		  axarr[i,j].axis('off')
              k += 1
   plt.show()
