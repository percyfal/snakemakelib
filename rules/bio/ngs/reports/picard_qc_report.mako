Picard QC report
=============================


+-----------------------+----------------------+
| Project               | ${"{:>20s}".format(project_name)} |
+-----------------------+----------------------+
| Application           | ${"{:>20s}".format(application)} |
+-----------------------+----------------------+
% if region:
| Region                | ${"{:>20s}".format(region)} |
+-----------------------+----------------------+
% endif
% if regionsummary:
| Number of samples     | ${"{:>20d}".format(len(samples))} |
+-----------------------+----------------------+



Region summary
----------------

.. csv-table:: Region QC summary. Columns show region name and for each region the *average* numbers over samples for total number of reads, percent aligned reads, percent duplication, mean insert size, mean coverage over target regions, percent sequenced bases that have aligned to target, followed by the percent bases in target regions covered at 10X and 30X
   :class: docutils
   :file: ${reportdir}${summarytable}
   :header-rows: 1


% else:

Sample QC summary
------------------


.. csv-table:: Sample QC summary. Columns show sample name, total number of reads, percent aligned reads, percent duplication, mean insert size, mean coverage over target regions, percent sequenced bases that have aligned to target, followed by the percent bases in target regions covered at 10X and 30X
   :class: docutils
   :file: ${reportdir}${summarytable}
   :header-rows: 1

% endif


% if not regionsummary:

QC Metrics
----------

Sequence statistics
^^^^^^^^^^^^^^^^^^^

.. figure:: ${seqstats}


Alignment metrics
^^^^^^^^^^^^^^^^^

.. figure:: ${alnmet}

Duplication metrics
^^^^^^^^^^^^^^^^^^^^

.. figure:: ${dupmet}

Insert metrics
^^^^^^^^^^^^^^^^^^^^

.. figure:: ${insmet}

Target metrics
^^^^^^^^^^^^^^^^^^^^

.. figure:: ${targetmet}


.. figure:: ${target2dup}


Hybridization metrics
^^^^^^^^^^^^^^^^^^^^^

.. figure:: ${hsmet}


Sample-based hybridization metrics
----------------------------------

% for f in hsmetsub:
.. figure:: ${f}
   :align: center

% endfor

% endif
