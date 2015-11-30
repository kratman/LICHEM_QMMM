
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

### Notes

The test suite prints run times for the tests. Since LICHEM has different
settings for different wrappers, the times only represent the efficiency of
the wrappers. The run times are not for the comparison of the efficiency of
different packages.

The test suite is not compatable with python 3.0 or higher.

