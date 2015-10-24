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
 

 References for integrals:
 

*/

//Gaussian function classes and structs
double GauDen1s::ChrgNuc(double qpc, Coord pos, double Rcut)
{
  //Function to calculate interactions between nuclei and charges
  double E = 0;
  
  return E;
};

double GauDen1s::NucNuc(GauDen1s gau2, double Rcut)
{
  //Function to calculate interactions between nuclei
  double E = 0;
  
  return E;
};

double GauDen1s::TwoOver(GauDen1s gau2)
{
  //Calculates the overlap of two Gaussian functions
  double Sij = 0.0;
  
  return Sij;
};

double GauDen1s::OneCoulPC(double qpc, Coord pos, double Rcut)
{
  //Calculates the electron-charge interactions for a point-charge
  double E = 0.0;
  
  return E;
};

double GauDen1s::OneCoulNuc(GauDen1s gau2, double Rcut)
{
  //Calculates the electron-charge interactions for a point-charge
  double E = 0.0;
  
  return E;
};

double GauDen1s::TwoCoul(GauDen1s gau2, double Rcut)
{
  //Calculates the electron-electron interactions
  double E = 0.0;
  
  return E;
};

//Functions for conversion to point-charges

