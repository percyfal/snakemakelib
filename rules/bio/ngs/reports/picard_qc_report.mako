Picard QC report
=============================

:Project: ${project_name}
:Application: ${application}
% if region:
:Region: ${region}
% endif


Sample QC summary
------------------

.. csv-table:: Sample QC summary. Columns show sample name, total number of reads, percent aligned reads, percent duplication, mean insert size, mean coverage over target regions, percent sequenced bases that have aligned to target, followed by the percent bases in target regions covered at 10X and 30X
   :class: docutils
   :file: report/picardmetricssummary.csv
   :header-rows: 1

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
