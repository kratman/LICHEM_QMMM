/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Global variables used in LICHEM.

*/

//Make including safe
#ifndef LICHEM_GLOBALS
#define LICHEM_GLOBALS

//Namespace for global variables
namespace LICHEMGlobal
{
  //Global variables
  int globalSys = 0; //Global dummy return for all system calls
  string xyzFilename; //Saves a filename given in the arguments
  string conFilename; //Saves a filename given in the arguments
  string regFilename; //Saves a filename given in the arguments
  int Nthreads = 1; //Total number of threads available
  int Ncpus = 1; //Number of processors for QM calculations
  int Nfreeze = 0; //Number of frozen atoms
  int Npseudo = 0; //Number of pseudo-bonds
  int Nbound = 0; //Number of boundary-atoms
  int Natoms = 0; //Total number of atoms
  int Nqm = 0; //Number of QM atoms
  int Nmm = 0; //Number of MM atoms
  double mcStep = 2*stepMin; //Monte Carlo step size
  double Lx = 10000.0; //Box length
  double Ly = 10000.0; //Box length
  double Lz = 10000.0; //Box length

  //Flags for simulation options
  bool GEM = 0; //Flag for frozen density QMMM potential
  bool AMOEBA = 0; //Flag for polarizable QMMM potential
  bool CHRG = 0; //Flag for point-charge QMMM potential
  bool PSI4 = 0; //Wrapper flag
  bool NWChem = 0; //Wrapper flag
  bool Gaussian = 0; //Wrapper flag
  bool TINKER = 0; //Wrapper flag
  bool LAMMPS = 0; //Wrapper flag
  bool PBCon = 0; //Flag for the boundary conditions
  bool QMMM = 0; //Flag for the type of wrapper
  bool MMonly = 0; //Flag for the type of wrapper
  bool QMonly = 0; //Flag for the type of wrapper
  bool OptSim = 0; //Flag for energy minimization with QM packages
  bool SteepSim = 0; //Flag for steepest descent minimization in LICHEM
  bool DFPSim = 0; //Flag for DFP minimization in LICHEM
  bool NEBSim = 0; //Flag for NEB path optimization in LICHEM
  bool PIMCSim = 0; //Flag for Monte Carlo
  bool FBNEBSim = 0; //Flag for force-bias NEB Monte Carlo
  bool FreqCalc = 0; //Flag for a frequency calculation
  bool SinglePoint = 0; //Flag for energy calculation
  bool GauExternal = 0; //Runs Gaussian with External

  //Timers
  int startTime = 0; //Time the calculation starts
  int endTime = 0; //Time the calculation ends
  int QMTime = 0; //Sum of QM wrapper times
  int MMTime = 0; //Sum of MM wrapper times
};

#endif

