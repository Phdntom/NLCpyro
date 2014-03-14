NLCpyro
=======

A collection numerical linked cluster routines for finding the thermodynamic properties of pyrochlore quantum spin ices.

clusters.py
===========




constants.py
============
Contains constants for a collection of Quantum Spice Ice Models.

ORDER specifies the expansion order in powers of number of clusters.
e.g. ORDER=3, uses all cluster embeddings using 3 tetrahedra

other constants include lattice basis vectors xi,yi,zi as well as
phase factors a,b in the zeta and gamma matrix




diagonalize.py
==============
Builds and diagonalizes clusters to 3rd or 4th order. It performs necessary
subgraph subraction and calculates final thermodynamic properties.
The output files are named as

[PROPERTY][ORDER]_J_[Jzz]_[Jpm]_[Jpmpm]_[Jzpm]_h_[fieldmag]_d_[fielddir]

[PROPERTY] = C, S, M ( specific heat, entropy, magnetization)

[ORDER] = 3 or 4

The file is output to files organized in blocks...

The file is organized with rows representing a temperature point.

Each point contains the different orders of the calculation, e.g.

In the case of 3rd ORDER the block is:
TEMP ORDER_0 ORDER_1 ORDER_2 ORDER_3 EULER_3

In the case of 4th ORDER the block is:
TEMP ORDER_0 ORDER_1 ORDER_2 ORDER_3 ORDER_4 EULER_3 EULER_4

The file contains as many lines as specified by [T_POINTS] in the file
"constants.py"

To see the parameters listed prior to running, you can type

$ python contants.py

To run the program with the parameters in "consants.py" type

$ python diagonalize.py

Python must be installed along with the numpy library.
