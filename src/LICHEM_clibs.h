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

//Header files
#ifdef _OPENMP
 //Use OpenMP
 #pragma message("OpenMP is enabled.")
 #include <omp.h>
#endif
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <algorithm>

#endif

