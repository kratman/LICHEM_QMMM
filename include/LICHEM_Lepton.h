/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Variables and options for the eFF semi-classical electron model.

 Reference for eFF:
 Su et al., Phys. Rev. Lett., 18, 99, 185002, (2007)

*/

//Make including safe
#ifndef LICHEM_LEPTON_HEADERS
#define LICHEM_LEPTON_HEADERS

namespace LICHEMLepton
{
  //Compile time eFF options
  const double Rad = 0.10; //Initial electron radius
  const double rho = 1.0; //Parameter for eFF VB mixing
  const double sbar = 1.0; //Parameter for eFF radius scaling
  const double rbar = 1.0; //Parameter for eFF distance scaling
  const double Charm = 0; //Force constant for eFF harmonic constraint
  const bool Scale_eFF = 0; //Scale the eFF kinetic energy by ElrtNbeads
  const double Scale_POW = 0.5; //ElrtNbeads = pow(P,Scale_POW)
  const bool qGrid = 0; //Restrict point-charge movements
  const double RadMin = 0.01; //Minimum electron radius
  const double RadMax = 25.0; //Maximum electron radius
  const double ElecCutoff = 15.0; //Cutoff for electron electrostatics

  //PIMC move probabilities
  double ElBeadProb = 0.25; //Probability to move an electron bead
  double ElCentProb = 0.25; //Probability to move an electron centroid
  double RadProb = 0; //Probability to change electron radius
  double SwapProb = 0.05; //Probability to swap spins
  double FlipProb = 0.05; //Probability to flip a single spin

  //Global variables
  double ElrtNbeads; //Scale factor for eFF kinetic energy
};

#endif

