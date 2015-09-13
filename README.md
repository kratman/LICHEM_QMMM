
[//]: # (Mixture of GitHub markdown and HTML. HTML is needed for formatting.)

***
<div align=center> <h2>
LICHEM: Layered Interacting CHEmical Models
</h2> </div>

<div align=center> <h4> By: Eric G. Kratz </h4> </div>

<div align=center> <h3> Symbiotic Computational Chemistry </h3> </div>
***

### LICHEM: A QMMM interface for polarizable force fields

<h4>
NOTICE: The multipole routines are currently only in the private development
repository. The multipole functionality will be uploaded when the LICHEM paper
has been accepted.
</h4>

### Introduction

This package is designed to be an open source (GPLv3) interface between QM
and MM software so that QMMM calculations can be performed with polarizable
and frozen electron density force fields. Functionality is also present for
standard point-charge based force fields, pure MM, and pure QM calculations.

Available calculations: single-point energies, optimized geometries, classical
molecular dynamics, classical Monte Carlo, path-integral Monte Carlo, and
reaction pathways
```
Available QM wrappers: Gaussian, NWChem, PSI4

Available MM wrappers: TINKER

MM wrappers in development: AMBER, LAMMPS
```

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

### Updates

LICHEM is still in the early stages of development. Please report all bugs so
that the code can be updated quickly. A list of known issues and the
development timeline are available in the doc directory.

Development versions of LICHEM are kept in a private repository. Large changes
to the code or the addition of new features (multipoles, diffuse charges,
reaction path ways, etc) are tested before they are merged into the public
repository. To ensure that LICHEM compiles properly, dummy function calls are
added for incomplete or untested features.

Relatively minor updates (documentation, printing, bug fixes, comments, etc)
are released ASAP. Updates from the private development branch are generally
merged into the public repository as a single commit around midnight.

### Testing

Test calculations can be performed with the runtests script in the tests
directory.

Tests are performed for pairs of QM and MM wrappers.
```
user:$ ./runtests Ncpus QMPackage MMPackage
```

A dry run can be performed to check packages without performing the
calculations.
```
user:$ ./runtests Ncpus QMPackage MMPackage Dry
```

Currently, only single-point energies and packages are tested. Optimizations
and pseudo-bond calculations will be added in the future.

### Development

The development of LICHEM was supported by funding from the NIH (Grant No.
R01GM108583) and Wayne State University. LICHEM is maintained by the Cisneros
research group at Wayne State University.

Developers: <br>
  Cisneros group, Wayne State University <br>
  Piquemal group, UPMC University
