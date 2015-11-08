/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for frozen density QMMM calculations using the
 Gaussian electrostatic model (GEM).

 NB: Frozen density and GEM routines are being held back until they are
 tested and the work is published.

 Reference for GEM:
 

 Reference for conversion to point-charges:
 

*/

//GEM utility functions
double GEMC6(double C6, Coord& POSi, Coord& POSj, double Rcut)
{
  //Function to calculate the LJ style dispersion
  double Eij = 0; //Dispersion energy
  
  //Return energy
  return Eij;
};

double GEMBuffC7()
{
  //Function to calculate buffered 14-7 style dispersion
  double Eij = 0; //Dispersion energy
  
  //Return energy
  return Eij;
};

//Functions to calculate GEM energy


//Functions for conversion to point-charges

