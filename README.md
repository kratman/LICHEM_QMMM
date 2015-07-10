***
<div align=center> <h2>
LICHEM: Layered Interacting CHEmical Models
</h2> </div>

<div align=center> <h4> By: Eric G. Kratz </h4> </div>

<div align=center> <h3> Symbiotic Computational Chemistry </h3> </div>
***

### LICHEM: A QMMM interface for polarizable force fields

### Introduction

This package is designed to be an open source (GPLv3) interface between QM
and MM software so that QMMM calculations can be performed with polarizable
and frozen density force fields. Functionality is also present for standard
point-charge based force fields, pure MM, and pure QM calculations.

Available calculations: single-point energies, optimized geometries, classical
molecular dynamics, classical Monte Carlo, path-integral Monte Carlo, and
reaction pathways

```
Available QM wrappers: Gaussian, NWChem, PSI4

Available MM wrappers: TINKER, AMBER, LAMMPS
```

Note: Some features will be held back until the research is published. Please
report all bugs so that the code can be updated quickly. Bug fixes are
released on a rolling basis, while major updates are developed and tested in a
private repository. A list of known issues and the development timeline are
available in the doc directory. 

### Installation

Currently, the binary and user's manual are not included in the repository.
The Makefile can be used to generate both files. Since LICHEM is designed to
be simple, only a small number of packages are required to compile the code.
An approximate list of packages is given below.
```
 LICHEM binary: OpenMP, Eigen3
 LICHEM manual: LaTeX, BibTeX
```

To install LICHEM, clone the git repository

```
user:$ mkdir LICHEM 
user:$ git clone https://github.com/kratman/LICHEM_QMMM.git ./LICHEM/
```

or unpack the zipped source code

```
user:$ mkdir LICHEM
user:$ cd LICHEM/
user:$ unzip LICHEM_QMMM-master.zip
user:$ mv LICHEM_QMMM-master/* .
user:$ rmdir LICHEM_QMMM-master
```

The LICHEM binary will eventually be included with the zipped source code,
however, modified or git source code can be compiled with the Makefile
provided with LICHEM. On Ubuntu boxes, the Makefile should function
without modifications. However, with other operating systems it may be
necessary to change the path to the Eigen3 package.

```
Default: -I/usr/include/eigen3/
```

The Makefile can produce both the documentation and the binary.

```
user:$ make install
```

For development, debugging, or testing an alternate binary can be compiled.

```
user:$ make Dev
```

Additional make rules can be found in the Makefile.
