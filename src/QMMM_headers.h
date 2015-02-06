/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 Headers and globals for FLUKE. This must be the first file
 imported into main().

*/

//Header Files
#include <omp.h>
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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <Eigen/StdList>
#include <Eigen/Eigen>
#include <Eigen/StdVector>
using namespace Eigen;
using namespace std;

//Compile options
const bool Jokes = 1; //Print humorous comments
const bool Debug = 0; //Turn debugging on/off
const bool Isotrop = 1; //Force isotropic expansion
const double StepMin = 0.01; //Minimum step size
const double StepMax = 1.0; //Maximum step size
const double Centratio= 5.0; //Scales 'step' for centroids
const int Acc_Check = 5000; //Eq steps before checking accratio

//Move Probabilities for PIMC
//Note: These probabilities allow for multi-particle moves
double BeadProb = 0.55; //Probability to move a single bead
double CentProb = 0.55; //Probability to move a centroid
double VolProb = 0.10; //Volume change probability

//Global constants
const double k = 8.6173324e-5; //Boltzmann constant (eV)
const double hbar = 6.58211928e-16; //Reduced Planck Constant (eV)
const double hbarSI = 1.054571726e-34; //Reduced Planck Constant (SI)
const double kb = 0.69503476; //Boltzmann constant (cm-1)
const double kSI = 1.3806488e-23; //Boltzmann constant (SI)
const double m2Ang = 1.0e10; //Angstroms to meters
const double amu2kg = 1.660538921e-27; //Atomic mass units to kg
const double cs = 2.99792458e8; //Speed of light (m)
const double pi = 4*atan(1); //Pi
const double h = 2*pi*hbar; //Planck Constant (eV)
const double SI2eV = 1/(1.602176565e-19); //Convert SI to eV
const double ToeV = amu2kg*SI2eV/(m2Ang*m2Ang); //Convert to eV units
const double C2eV = m2Ang/(4*pi*SI2eV*8.854187817e-12); //Coulomb to eV
const double Masse = 9.10938291e-31; //Mass of an electron (kg)
const double ElecMass = 5.4857990943e-4; //Mass of an electron (amu)
const double BohrRad = 0.52917721092; //Bohr radius (Ang)
const double Har2eV = 27.21138386; //Hartrees to eV
const double atm2eV = SI2eV*1.01325e-25; //atmA^3 to eV
const double Na = 6.02214179e23; //Avogadro's number
const double kcal2eV = 4184*SI2eV/Na; //kcal/mol to eV
const double sqrt2 = pow(2,0.5); //Square root of 2
const double HugeNum = 100000.0; //Large number to reject step
const double fs2s = 1e-12; //Convert fs to s

//Globals
string xyzfilename; //Saves the filename given in the arguments
string confilename; //Saves the filename given in the arguments
string regfilename; //Saves the filename given in the arguments
int Ncpus = 1; //Number of processors for QM calculations
int Nfreeze = 0; //Number of frozen atoms
int Npseudo = 0; //Number of pseudo-atoms
int Nbound = 0; //Number of boundary-atoms
int Natoms = 0; //Total number of atoms
int Nqm = 0; //Number of QM atoms
int Nmm = 0; //Number of MM atoms
double step = StepMin; //PIMC step size
double Lx = 500.0; //Box length
double Ly = 500.0; //Box length
double Lz = 500.0; //Box length

//Flags for simulation options
int GEM = 0; //Flag for frozen density QMMM potential
int AMOEBA = 0; //Flag for polarizable QMMM potential
int CHRG = 0; //Flag for point-charge QMMM potential
int PBCon = 0; //Flag for the boundary conditions
int PSI4 = 0; //Wrapper flag
int Gaussian = 0; //Wrapper flag
int TINKER = 0; //Wrapper flag
int LAMMPS = 0; //Wrapper flag
int AMBER = 0; //Wrapper flag
bool QMMM = 0; //Flag for the type of wrapper
bool MMonly = 0; //Flag for the type of wrapper
bool QMonly = 0; //Flag for the type of wrapper
bool OptSim = 0; //Flag for energy minimization with QM packages
bool SteepSim = 0; //Flag for steepest descent minimization in FLUKE
bool DFPSim = 0; //Flag for DFP minimization in FLUKE
bool MDSim = 0; //Flag for a NVT MD simulation
bool PIMCSim = 0; //Flag for Monte Carlo
bool PathSim = 0; //Flag for reaction paths
bool SinglePoint = 0; //Flag for energy calculation
bool GauExternal = 0; //Runs Gaussian with External
bool TranState = 0; //Flag to search for a transition state

//Timers
int StartTime = 0; //Time the calculation starts
int EndTime = 0; //Time the calculation ends
int QMTime = 0; //Sum of QM wrapper times
int MMTime = 0; //Sum of MM wrapper times

//Custom data types
struct Coord
{
  double x; //x position
  double y; //y position
  double z; //z position
};

struct Mpole
{
  //Cartesian multipoles
  //Atoms for the local frame of reference
  bool ChiralFlip; //Flip y axis
  string Type; //Bisector, Z-then-X, Z-Only, 3-Fold, or Z-Bisect
  int Atom1; //Atom which defines the z axis
  int Atom2; //Atom which defines the x axis
  int Atom3; //Atom which defines the y axis (chiral only)
  //Monopole moment
  double q;
  //Cartesian dipole moments
  double Dx;
  double Dy;
  double Dz;
  //Cartesian induced dipole moments
  double IDx;
  double IDy;
  double IDz;
  //Cartesian quadrupole moments
  double Qxx;
  double Qyy;
  double Qzz;
  double Qxy;
  double Qxz;
  double Qyz;
};

struct RedMpole
{
  //Reduced multipole from sph. harm. and diagonalization
  //Monopole
  double Q00;
  //Dipole moments
  double Q10;
  double Q11c;
  double Q11s;
  //Quadrupole moments
  double Q20;
  double Q22c;
  //Spherical harmonic vectors
  Matrix3d Vec;
};

struct OctCharges
{
  //A grid of point-charges which replaces multipoles
  double q1; //Charge in the +x direction
  double q2; //Charge in the +y direction
  double q3; //Charge in the +z direction
  double q4; //Charge in the -x direction
  double q5; //Charge in the -y direction
  double q6; //Charge in the -z direction
  //Vectors for the quadrupole frame of reference
  Vector3d Vecx;
  Vector3d Vecy;
  Vector3d Vecz;
  //Positions of the charges in the global frame
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  double x3;
  double y3;
  double z3;
  double x4;
  double y4;
  double z4;
  double x5;
  double y5;
  double z5;
  double x6;
  double y6;
  double z6;
};

struct QMMMAtom
{
  //Data type for atomic information
  double m; //Mass of atom
  bool QMregion; //QM, MM, pseudo-atom, or boundary-atom
  bool MMregion; //QM, MM, pseudo-atom, or boundary-atom
  bool PAregion; //QM, MM, pseudo-atom, or boundary-atom
  bool BAregion; //QM, MM, pseudo-atom, or boundary-atom
  bool Frozen; //Part of a frozen shell
  string QMTyp; //Real atom type
  string MMTyp; //Force field atom type
  int NumTyp; //Numerical atom type (if used)
  int NumClass; //Numerical atom class (if used)
  int id; //Atom number, starts at zero
  vector<int> Bonds; //Connectivity
  double Ep; //Storage for PI energies
  vector<Coord> P; //Array of beads
  vector<Coord> Vel; //Array of velocities
  vector<Mpole> MP; //Multipoles
  vector<OctCharges> PC; //Point-charge multipoles
};

struct QMMMElec
{
  //Data type for electronic information (eFF model)
  string typ; //Lepton type
  double m; //mass (amu)
  double q; //Charge (au)
  int spin; //Spin
  double Ep; //Temporary energy for parallel
  vector<Coord> P; //Bead coordinates
  vector<double> rad; //Radius (Ang)
};

struct QMMMSettings
{
  //Input needed for QM wrappers
  string Func; //DFT functional
  string Basis; //Basis set for QM calculations
  string RAM; //Ram for QM calculations
  string Charge; //QM total charge
  string Spin; //QM total spin
  double Eqm; //QM total energy
  //Input needed for MM wrappers
  double Emm; //MM total energy
  //Input needed for MC and MD functions
  string Ensemble; //NVT or NPT
  double Temp; //Temperature
  double Beta; //Inverse temperature
  double Press; //External pressure
  int Neq; //Number of equilibration run steps
  int Nsteps; //Number of production run steps
  int Nbeads; //Number of time-slices or beads
  double accratio; //Target acceptance ratio
  int Nprint; //Number of steps before printing
  double dt; //MD timestep
  double tautemp; //Thermostat time constant
  //Input needed for optimizations
  int MaxOptSteps; //Maximum iterative optimization steps
  double MMOptTol; //Criteria to end the optimization
  double QMOptTol; //Criteria to end the optimization
  double StepScale; //Steepest descent step size (Ang)
  double MaxStep; //Maximum size of the optimization step
};

//Function declarations
void PrintFancyTitle();

bool CheckFile(const string&);

double Bohring(double);

double CoordDist2(Coord&,Coord&);

string Typing(int);

int RevTyping(string);

void ExtractTINKpoles(vector<QMMMAtom>&,int);

void RotateTINKCharges(vector<QMMMAtom>&,int);

void WriteTINKMpole(vector<QMMMAtom>&,fstream&,int,int);

RedMpole Cart2SphHarm(Mpole&);

OctCharges SphHarm2Charges(RedMpole);

void FindTINKERClasses(vector<QMMMAtom>&);

double TINKEREnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKEROpt(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKERForces(vector<QMMMAtom>&,vector<Coord>&,QMMMSettings&,int);

double TINKERMMForces(vector<QMMMAtom>&,vector<Coord>&,QMMMSettings&,int);

double TINKERPolEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double AMBEREnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double AMBEROpt(vector<QMMMAtom>&,QMMMSettings&,int);

double LAMMPSWrapper(string,vector<QMMMAtom>&,QMMMSettings&,int);

double LAMMPSEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double LAMMPSOpt(vector<QMMMAtom>&,QMMMSettings&,int);

double LAMMPSForces(vector<QMMMAtom>&,vector<Coord>&,QMMMSettings&,int);

double GaussianEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double GaussianOpt(vector<QMMMAtom>&,QMMMSettings&,int);

void ExternalGaussian(int&,char**&);

double GaussianForces(vector<QMMMAtom>&,vector<Coord>&,QMMMSettings&,int);

void GaussianCharges(vector<QMMMAtom>&,QMMMSettings&,int);

double PSIEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double PSIOpt(vector<QMMMAtom>&,QMMMSettings&,int);

double PSIForces(vector<QMMMAtom>&,vector<Coord>&,QMMMSettings&,int);

void PSICharges(vector<QMMMAtom>&,QMMMSettings&,int);

double EFFEnergy(QMMMAtom&,QMMMElec&,int);

double KineticE_eFF(vector<QMMMElec>&,QMMMSettings&);

double EFFCorr(QMMMElec&,QMMMElec&,int);

double Get_EeFF(vector<QMMMAtom>&,vector<QMMMElec>&,QMMMSettings&);

bool OptConverged(vector<QMMMAtom>&,vector<QMMMAtom>&,vector<Coord>&,
     int,QMMMSettings& QMMMOpts,int,bool);

void FLUKESteepest(vector<QMMMAtom>&,QMMMSettings&,int);

void FLUKEDFP(vector<QMMMAtom>&,QMMMSettings&,int);

double BerendsenThermo(vector<QMMMAtom>&,QMMMSettings&,int);

void VerletUpdate(vector<QMMMAtom>&,QMMMSettings&,fstream&,bool,int);

double SpringEnergy(double,double);

double Get_PI_Espring(vector<QMMMAtom>&,QMMMSettings&);

double Get_PI_Epot(vector<QMMMAtom>&,QMMMSettings&);

bool MCMove(vector<QMMMAtom>&,QMMMSettings&);

void Get_Centroid(QMMMAtom&,QMMMSettings&);

Coord Get_COM(vector<QMMMAtom>&,QMMMSettings&);

void TINK2FLUKE(int&,char**&);

void Remove_COM(vector<QMMMAtom>&,QMMMSettings&);

void Print_traj(vector<QMMMAtom>&,fstream&,QMMMSettings&);

void ReadArgs(int&,char**&,fstream&,fstream&,fstream&,fstream&);

void ReadFLUKEInput(fstream&,fstream&,fstream&,
     vector<QMMMAtom>&,QMMMSettings&);

void FLUKEErrorChecker(QMMMSettings&);

void FLUKEPrintSettings(QMMMSettings&);

void GetQuotes(vector<string>&);

void BurstTraj(vector<QMMMAtom>&,string&,QMMMSettings&);

//Function definitions
#include "Core_funcs.cpp"
#include "Input_Reader.cpp"
#include "TINK2FLUKE.cpp"
#ifdef DEVCOMP
#include "Real_Multipoles.cpp"
#endif
#ifndef DEVCOMP
#include "Multipoles.cpp"
#endif
#ifdef DEVCOMP
#include "Real_Frozen_density.cpp"
#endif
#ifndef DEVCOMP
#include "Frozen_density.cpp"
#endif
#include "PathIntegral.cpp"
#ifdef DEVCOMP
#include "Real_ReactionPath.cpp"
#endif
#ifndef DEVCOMP
#include "ReactionPath.cpp"
#endif
#include "Optimizers.cpp"
#include "Dynamics.cpp"
#include "Analysis.cpp"

//Wrapper definitions
#include "Lepton_eng.cpp"
#include "Gaussian.cpp"
#include "TINKER.cpp"
#include "LAMMPS.cpp"
#include "AMBER.cpp"
#include "PSI4.cpp"

