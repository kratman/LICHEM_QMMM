/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for AMBER.

 Reference for AMBER:
 Case et al., AMBER 14, (2015)

*/

//MM utility functions


//MM wrapper functions
double AMBERForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                   QMMMSettings& QMMMOpts, int Bead)
{
  //Runs AMBER energy calculations
  double E = 0.0;
  
  //Change units
  E *= kcal2eV;
  return E;
};

double AMBEREnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs AMBER energy calculations
  double E = 0.0;
  
  //Change units
  E *= kcal2eV;
  return E;
};

double AMBEROpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs AMBER optimizations
  double E = 0.0;
  
  //Change units
  E *= kcal2eV;
  return E;
};

