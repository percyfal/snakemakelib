:Project: ${project_name}
:Application: ${application}

Rule graphs and DAGs
---------------------

% for r in rulegraphs:

${r}
${"^" * (len(r) + 1)}

.. figure:: ${r}
   :align: center

   Rulegraph for ${r}

% endfor
