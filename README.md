
[//]: # (Mixture of GitHub markdown and HTML. HTML is needed for formatting.)

***
<div align=center> <h2>
LICHEM: Layered Interacting CHEmical Models
</h2> </div>

<div align=center> <h4> By: Eric G. Kratz </h4> </div>

<div align=center> <h3> Symbiotic Computational Chemistry </h3> </div>
***

### LICHEM: A QM/MM interface for polarizable force fields

[![GPL license](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat)](https://github.com/kratman/LICHEM_QMMM/blob/master/src/GPL_LICENSE)

### Automated tests

**Travis CI build:** [![Build Status](https://travis-ci.org/kratman/LICHEM_QMMM.svg?branch=master)](https://travis-ci.org/kratman/LICHEM_QMMM)

### Introduction

This package is an open-source (GPLv3) interface between QM and MM software
so that QM/MM calculations can be performed with polarizable and frozen
electron density force fields. Functionality is also present for standard
point-charge based force fields, pure MM, and pure QM calculations.

Available calculations: single-point energies, geometry optimizations,
classical Monte Carlo, path-integral Monte Carlo, and reaction pathways.
```
Available QM wrappers: Gaussian, NWChem, PSI4

Available MM wrappers: TINKER

MM wrappers in development: LAMMPS
```

### Citing LICHEM

If you use LICHEM in your research, please cite the following paper:

```
Kratz, E.G.; Walker, A.R.; Lagardere, L; Lipparini, F; Piquemal, J.-P.; and
Cisneros, G.A.; "LICHEM: A QM/MM Program for Simulations with Multipolar and
Polarizable Force Fields", J. Comput. Chem., 37, 11, 1019, (2016)
```

### Installation

Currently, the binary and user's manual are not included in the repository.
However, the Makefile can be used to generate both files. Since LICHEM is
designed to be simple, only a small number of packages are required to compile
the code. An approximate list of packages is given below.
```
 LICHEM binary: OpenMP
 LICHEM test suite: Python
 LICHEM manual: LaTeX, BibTeX, TeXLive
```
A copy of the Eigen3 library is included with the LICHEM source code.
However, other versions of Eigen3 can be used to build LICHEM by modifying
the Makefile.

To install LICHEM, clone the git repository:
```
user:$ mkdir LICHEM
user:$ git clone https://github.com/kratman/LICHEM_QMMM.git ./LICHEM/
```

On Ubuntu boxes, the Makefile should function without modifications. However,
it may be necessary to install additional LaTeX packages. On OSX machines,
the SEDI, TEX, BIB, and CXXFLAGS variables will need to be modified.

The Makefile can produce both the documentation and the binary.
```
user:$ make install
```

For development, debugging, or testing an alternate binary can be compiled.
```
user:$ make Dev
```

Additional make rules can be found in the Makefile.

### QM and MM packages

LICHEM uses QM and MM packages that are pre-installed and located in the
user's path. The LICHEM wrappers are kept up-to-date with the development of
the QM and MM packages. It is best to use packages from conda, git, or other
package management systems.

### Updates

LICHEM is still in the early stages of development. Please report all bugs so
that the code can be updated quickly. Bugs should be reported through the
GitHub issues interface. A list of known issues and the development timeline
are available in the doc directory or on the GitHub page.

Development versions of LICHEM are kept in a private repository. Large changes
to the code or the addition of new features (diffuse charges, reaction path
ways, etc) are tested and published before they are merged into the public
repository. To ensure that LICHEM compiles properly, dummy function calls are
added for incomplete or untested features.

Relatively minor updates (documentation, output printing, bug fixes, comments,
etc) are released ASAP. Updates from the private development branch are
generally merged into the public repository as a single commit.

Further development details can be found in the doc directory, or the GitHub
issues section.

### Jokes

Jokes and Easter eggs can be included by changing
```
 const bool JOKES = 0; //Print humorous comments
```
in LICHEM_options.h to
```
 const bool JOKES = 1; //Print humorous comments
```
before compiling the code.

### Testing

Test calculations can be performed with the runtests script in the tests
directory.

Tests can be performed for pairs of QM and MM wrappers.
```
user:$ ./runtests Ncpus QMPackage MMPackage
```

Additionally, tests can be run for all wrappers at once.
```
user:$ ./runtests Ncpus All
```

A dry run can be performed to check packages without performing the
calculations.
```
user:$ ./runtests Ncpus QMPackage MMPackage Dry
```

If tests are consistently failing, please post details on the GitHub issues
section.

### Development

The development of LICHEM was supported by funding from the NIH (Grant No.
R01GM108583) and Wayne State University. LICHEM is maintained by the Cisneros
research group at University of North Texas.

Developers:
<ul>
  <li>Cisneros group, University of North Texas
</ul>

