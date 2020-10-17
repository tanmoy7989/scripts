This is a loose collection of scripts that I find useful time to time in my research work. 

1) ```visualize_multiscale.py```:

This script creates a typically requested multiscale rendering of an integrative model by combining detailed structures from the PDB files of rigid bodies used in the model construction as well as flexible beads extracted from a representative RMF file (usually the centroid of the top cluster after structural clustering).

RMF is a specialized file format for storing multiscale, multi-representation coordinate data from integrative models generated with the Salilab software [IMP](https://integrativemodeling.org/). Further details about the RMF format can be found [here](https://github.com/salilab/rmf).

This script will create:
a) A PDB file constructed by stitching together the different rigid
body PDB files after aligning them with their corresponding particles in
the given RMF file.

b) If requested, will write out a chimerax script that displays the
multiscale (PDB + unstructured parts from RMF) nature of the model.

c) Additionally, if crosslink (XL) data is provided, then it will also
map the crosslinks on to the multiscale structure.

Requirements:
[IMP](https://integrativemodeling.org/)
[Biopython](https://biopython.org/)
[pandas](https://pandas.pydata.org/)

==== more to come ====



