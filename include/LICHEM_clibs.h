/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Standard C++ headers used in LICHEM.

*/

//Make including safe
#ifndef LICHEM_CLIBS
#define LICHEM_CLIBS

//Header files for parallelization
#ifdef _OPENMP
 //Use OpenMP
 #pragma message("OpenMP is enabled.")
 #include <omp.h>
#else
 #pragma message("OpenMP is disabled.")
#endif
#ifdef _OPENACC
 //Use OpenMP
 #pragma message("OpenACC is enabled.")
 #include <openacc.h>
#else
 #pragma message("OpenACC is disabled.")
#endif

//General header files
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <algorithm>

#endif

