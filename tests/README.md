
[//]: # (Mixture of GitHub markdown and HTML. HTML is needed for formatting.)

***
<div align=center> <h2>
LICHEM: Layered Interacting CHEmical Models
</h2> </div>

<div align=center> <h4> By: Eric G. Kratz </h4> </div>

<div align=center> <h3> Symbiotic Computational Chemistry </h3> </div>
***

### LICHEM: A QMMM interface for polarizable force fields

### Test suite

This directory contains calculations for testing the LICHEM wrappers,
efficiency, and accuracy.

Tests are performed for pairs of QM and MM wrappers.
```
user:$ ./runtests Ncpus QMPackage MMPackage
```

Additionally, tests can be run for all wrappers at once.
```
user:$ ./runtests Ncpus All
```

A dry run can be performed to check packages without perfoming the
calculations.
```
user:$ ./runtests Ncpus QMPackage MMPackage Dry
```

If tests are consistently failing, please post details on the GitHub issues
section.

### Notes

The test suite prints run times for the tests. Since LICHEM has different
settings for different wrappers, the times only represent the efficiency of
the wrappers. The run times are not for the comparison of the efficiency of
different packages.

The test suite is not compatable with python 3.0 or higher.

###Tests

[//]: # (Table entries cannot have newlines)

| Test | Description | QM | MM |
| :--- | :--- | :---: | :---: |
| HF energy | Hartree-Fock energy calculated using only the QM wrapper. | PSI4,Gaussian | N/A |
| PBE0 energy | Density functional theory energy calculated using only the QM wrapper. | PSI4,Gaussian,NWChem | N/A |
| CCSD energy | Coupled-cluster energy calculated using only the QM wrapper. | PSI4 | N/A |
| PM6 energy | Semi-empirical energy calculated using only the QM wrapper. | Gaussian | N/A |
| Frequencies | Harmonic frequencies using only the QM wrapper. | PSI4,Gaussian,NWChem | N/A |
| NEB TS energy | Nudged elastic band optimization using only the QM wrapper. | PSI4,Gaussian,NWChem | N/A |
| TIP3P energy | MM energy of the water dimer with the TIP3P model. | N/A | TINKER |
| AMOEBA/GK energy | MM energy of the water dimer in the generalized Kirkwood implicit solvent. | N/A | TINKER |
| PBE0/TIP3P energy | QMMM energy of a water dimer calculated with PBE0 and TIP3P. | PSI4,Gaussian,NWChem | TINKER |
| PBE0/AMOEBA energy | Polarizable QMMM energy of a water dimer calculated with PBE0 and AMOEBA. | PSI4,Gaussian,NWChem | TINKER |
| DFP/Pseudobonds | QMMM Davidon-Fletcher-Powell optimization of 2-Butyne with the two methyl groups replaced by pseudobond/boundary atoms. | Gaussian,NWChem | TINKER |
