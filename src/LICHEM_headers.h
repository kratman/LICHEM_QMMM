/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Headers and globals for LICHEM.

*/

//Make including safe
#ifndef LICHEM_HEADERS
#define LICHEM_HEADERS

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

//Set namespaces for common libraries
using namespace Eigen;
using namespace std;

//Compile options
const bool Jokes = 0; //Print humorous comments
const bool Isotrop = 1; //Force isotropic expansion in NPT Monte Carlo
const double StepMin = 0.005; //Minimum Monte Carlo step size (Angstroms)
const double StepMax = 1.0; //Maximum Monte Carlo step size (Angstroms)
const double CentRatio= 5.0; //Scales step size for path-integral centroids
const int Acc_Check = 2000; //Eq Monte Carlo steps before checking accratio

//Move Probabilities for PIMC
/*

NB: The Monte Carlo engine allows for multi-particle moves:

If (BeadProb+CentProb) > 1) then there is a chance that multiple particles
move during a step.

If (BeadProb+CentProb) < 1) then there is a chance that no prarticles move
during a step.

If (VolProb > 0) and the simulation is not periodic, then VolProb does
nothing.

If the simulation is periodic, multi-particle moves will be common since
changing the box size moves all particles.

*/
double BeadProb = 0.55; //Probability to move all beads for an atom
double CentProb = 0.55; //Probability to move a centroid
double VolProb = 0.25; //Volume change probability

//Global exact constants
const double pi = 4*atan(1); //Pi
const double sqrt2 = pow(2,0.5); //Square root of 2
const double HugeNum = 1e50; //Large number to reject step
const double fs2s = 1e-12; //Convert fs to s
const double m2Ang = 1.0e10; //Angstroms to meters
const double atm2Pa = 1.01325e5; //Atmospheres to Pascal

//Global measured constants (NIST, CODATA 2010)
const double EpsZero = 8.54187817e-12; //Electric constant
const double hbar = 6.58211928e-16; //Reduced Planck Constant (eV)
const double k = 8.6173324e-5; //Boltzmann constant (eV)
const double kSI = 1.3806488e-23; //Boltzmann constant (SI)
const double amu2kg = 1.660538921e-27; //Atomic mass units to kg
const double SI2eV = 1/(1.602176565e-19); //Convert SI to eV
const double Masse = 9.10938291e-31; //Mass of an electron (kg)
const double BohrRad = 0.52917721092; //Bohr radius (Ang)
const double Har2eV = 27.21138505; //Hartrees to eV
const double Na = 6.02214129e23; //Avogadro's number
const double Debye2au = 0.393430307; //Convert from Debye to au

//Global derived constants
const double atm2eV = SI2eV*atm2Pa/(m2Ang*m2Ang*m2Ang); //atmA^3 to eV
const double C2eV = m2Ang/(4*pi*SI2eV*EpsZero); //Coulomb to eV
const double ElecMass = Masse/amu2kg; //Mass of an electron (amu)
const double ToeV = amu2kg*SI2eV/(m2Ang*m2Ang); //Convert to eV units (PIMC)
const double kcal2eV = 4184*SI2eV/Na; //kcal/mol to eV

//Globals
int GlobalSys = 0; //Global dummy return for all system calls
string xyzfilename; //Saves a filename given in the arguments
string confilename; //Saves a filename given in the arguments
string regfilename; //Saves a filename given in the arguments
int Nthreads = 1; //Total number of threads available
int Ncpus = 1; //Number of processors for QM calculations
int Nfreeze = 0; //Number of frozen atoms
int Npseudo = 0; //Number of pseudo-bonds
int Nbound = 0; //Number of boundary-atoms
int Natoms = 0; //Total number of atoms
int Nqm = 0; //Number of QM atoms
int Nmm = 0; //Number of MM atoms
double step = 2*StepMin; //Monte Carlo step size
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
bool AMBER = 0; //Wrapper flag
bool PBCon = 0; //Flag for the boundary conditions
bool QMMM = 0; //Flag for the type of wrapper
bool MMonly = 0; //Flag for the type of wrapper
bool QMonly = 0; //Flag for the type of wrapper
bool OptSim = 0; //Flag for energy minimization with QM packages
bool SteepSim = 0; //Flag for steepest descent minimization in LICHEM
bool QuickSim = 0; //Flag for QuickMin optimization in LICHEM
bool DFPSim = 0; //Flag for DFP minimization in LICHEM
bool NEBSim = 0; //Flag for NEB path optimization in LICHEM
bool ESDSim = 0; //Flag for ensemble steepest descent
bool PIMCSim = 0; //Flag for Monte Carlo
bool ENEBSim = 0; //Flag for ensemble NEB reaction paths
bool SinglePoint = 0; //Flag for energy calculation
bool GauExternal = 0; //Runs Gaussian with External

//Timers
int StartTime = 0; //Time the calculation starts
int EndTime = 0; //Time the calculation ends
int QMTime = 0; //Sum of QM wrapper times
int MMTime = 0; //Sum of MM wrapper times

//Custom classes
class PeriodicTable
{
  //Class for storing periodic table data
  private:
    //Atom types
    vector<string> Typs; //Atomic symbols
    //Approximate 1s Gaussian widths of the atoms
    vector<double> GauWids; //Diffuse charge widths
    //Bond distances
    vector<double> CovRadii; //Covalent radii
  public:
    //Set data (hard coded)
    PeriodicTable();
    //Retrieve data
    string Typing(int); //Atom type
    int RevTyping(string); //Atomic number
    double GetGauWid(string); //Gaussian (1s) width
    double GetCovRadius(string); //Covalent radius
};

class Coord
{
  public:
    double x; //x position
    double y; //y position
    double z; //z position
};

class Mpole
{
  //Cartesian multipoles
  public:
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
    //Cartesian induced dipole moments (global frame)
    double IDx;
    double IDy;
    double IDz;
    //Cartesian quadrupole moments (Q_ij = Q_ji)
    double Qxx;
    double Qxy;
    double Qxz;
    double Qyy;
    double Qyz;
    double Qzz;
};

class RedMpole
{
  //Reduced multipole from sph. harm. and diagonalization
  public:
    //Monopole
    double Q00;
    //Dipole moments
    double Q10; //Z component
    double Q11c; //X component
    double Q11s; //Y component
    //Quadrupole moments
    double Q20; //Z^2 component
    double Q22c; //X^2-Y^2 component
    //Spherical harmonic vectors
    Vector3d Vecx; //X direction in quadrupole frame
    Vector3d Vecy; //Y direction in quadrupole frame
    Vector3d Vecz; //Z direction in quadrupole frame
};

class OctCharges
{
  //A grid of point-charges which replaces multipoles
  public:
    double q1; //Charge in the +x direction
    double q2; //Charge in the +y direction
    double q3; //Charge in the +z direction
    double q4; //Charge in the -x direction
    double q5; //Charge in the -y direction
    double q6; //Charge in the -z direction
    //Positions of the charges in the global frame
    double x1; //Position of charge 1
    double y1; //Position of charge 1
    double z1; //Position of charge 1
    double x2; //Position of charge 2
    double y2; //Position of charge 2
    double z2; //Position of charge 2
    double x3; //Position of charge 3
    double y3; //Position of charge 3
    double z3; //Position of charge 3
    double x4; //Position of charge 4
    double y4; //Position of charge 4
    double z4; //Position of charge 4
    double x5; //Position of charge 5
    double y5; //Position of charge 5
    double z5; //Position of charge 5
    double x6; //Position of charge 6
    double y6; //Position of charge 6
    double z6; //Position of charge 6
};

class GauDen1s
{
  //Simple 1s Gaussian class
  private:
    //Properties
    double mag; //Magnitude/population (prefactor)
    double wid; //Width in a.u.
    double q; //Nuclear charge
    double x; //X position
    double y; //Y position
    double z; //Z position
  public:
    //Constructor
    GauDen1s(double magi,double widi,double qi,double xi,double yi,double zi)
    {
      //Save data given to the constructor
      mag = magi;
      wid = widi;
      q = qi;
      x = xi;
      y = yi;
      z = zi;
      //Convert to a.u.
      wid *= (BohrRad*BohrRad);
      return;
    }
    //Point-charge interactions
    double ChrgNuc(double,Coord,double); //Nuclei-charge electrostatic term
    double NucNuc(GauDen1s,double); //Nuclei-nuclei electrostatic term
    //Electron density integrals
    double TwoOver(GauDen1s); //Density-density overlap
    double OneCoulPC(double,Coord,double); //Density-charge (MM)
    double OneCoulNuc(GauDen1s,double); //Density-nucleus
    double TwoCoul(GauDen1s,double); //Density-density Coulomb repulsion
};

class HermGau
{
  //Class for Hermite Gaussians
  private:
    //Gaussian properties
    double mag; //Coefficient in front of the Gaussian
    double alpha; //Exponent (width)
    int powx; //Power in the x direction
    int powy; //Power in the y direction
    int powz; //Power in the z direction
    //Position
    double x; //X position in Angstroms
    double y; //Y position in Angstroms
    double z; //Z position in Angstroms
  public:
    //Constructor
    HermGau(double ci, double a, int ix, int iy, int iz,
            double xi, double yi, double zi)
    {
      //Save data given to the constructor
      mag = ci;
      alpha = a;
      powx = ix;
      powy = iy;
      powz = iz;
      x = xi;
      y = yi;
      z = zi;
      //Convert to a.u.
      alpha *= (BohrRad*BohrRad);
      return;
    }
    double Coeff(); //Return the coefficient
    double XPos(); //Return the x position
    double YPos(); //Return the y position
    double ZPos(); //Return the z position
    double Alpha(); //Return the Gaussian coefficient (width)
    int XPow(); //Return the Hermite power in the x direction
    int YPow(); //Return the Hermite power in the y direction
    int ZPow(); //Return the Hermite power in the z direction
};

class QMMMAtom
{
  //Data type for atomic information
  public:
    //Temporary storage
    double Ep; //Path-integral energies
    //Regions
    bool QMregion; //QM, MM, pseudo-bond, or boundary-atom
    bool MMregion; //QM, MM, pseudo-bond, or boundary-atom
    bool PBregion; //QM, MM, pseudo-bond, or boundary-atom
    bool BAregion; //QM, MM, pseudo-bond, or boundary-atom
    bool Frozen; //Part of a frozen shell
    //Force field information
    double m; //Mass of atom
    string QMTyp; //Real atom type
    string MMTyp; //Force field atom type
    int NumTyp; //Numerical atom type (if used)
    int NumClass; //Numerical atom class (if used)
    int id; //Atom number, starts at zero
    vector<int> Bonds; //Connectivity
    //Coordinates
    vector<Coord> P; //Array of beads
    //Multipoles
    vector<Mpole> MP; //Multipoles
    vector<OctCharges> PC; //Point-charge multipoles
};

class QMMMElec
{
  //Data type for electronic information (eFF model)
  public:
    //Temporary storage
    double Ep; //Path-integral energies
    //Particle properties
    string typ; //Lepton type
    double m; //mass (amu)
    double q; //Charge (au)
    //Coordinates
    vector<int> spin; //Spin
    vector<Coord> P; //Bead XYZ coordinates
    vector<double> rad; //Radius (Ang)
};

class QMMMSettings
{
  //Settings for the simulation and wrappers
  public:
    //Input needed for QM wrappers
    string Func; //QM method (functional, HF, etc)
    string Basis; //Basis set for QM calculations
    string RAM; //Ram for QM calculations
    bool MemMB; //Is the RAM in mb or gb
    string Charge; //QM total charge
    string Spin; //QM total spin
    string BackDir; //Directory for log file backups
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
    //Input needed for reaction paths
    double Kspring; //Elastic band spring constant
    int TSBead; //Current guess of the transition state
    bool Climb; //Flag to turn on climbing image NEB
    bool FrznEnds; //Flag to freeze the NEB end points
    //Storage of energies (NEB and PIMC)
    double Eold; //Temporary storage
    double Ereact; //Reactant energy
    double Eprod; //Product energy
    double Ets; //Transition state energy
};

//Set up periodic table
PeriodicTable PTable;

//Function declarations (alphabetical)
double AMBERForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

double AMBEREnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double AMBEROpt(vector<QMMMAtom>&,QMMMSettings&,int);

bool Angled(vector<QMMMAtom>&,int,int);

double Bohring(double);

bool Bonded(vector<QMMMAtom>&,int,int);

double BoysFunc(int,double);

void BurstTraj(vector<QMMMAtom>&,QMMMSettings&);

RedMpole Cart2SphHarm(Mpole&);

bool CheckFile(const string&);

double CoordDist2(Coord&,Coord&);

bool Dihedraled(vector<QMMMAtom>&,int,int);

double EFFCorr(QMMMElec&,QMMMElec&,int);

double EFFEnergy(QMMMAtom&,QMMMElec&,int);

void EnsembleNEB(vector<QMMMAtom>&,fstream&,QMMMSettings&);

void EnsembleSD(vector<QMMMAtom>&,fstream&,QMMMSettings&,int);

void ExternalGaussian(int&,char**&);

void ExtractGlobalPoles(int& argc, char**& argv);

void ExtractTINKpoles(vector<QMMMAtom>&,int);

int FindMaxThreads();

void FindTINKERClasses(vector<QMMMAtom>&);

void GaussianCharges(vector<QMMMAtom>&,QMMMSettings&,int);

double GaussianEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double GaussianForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

double GaussianOpt(vector<QMMMAtom>&,QMMMSettings&,int);

double Get_EeFF(vector<QMMMAtom>&,vector<QMMMElec>&,QMMMSettings&);

double Get_PI_Espring(vector<QMMMAtom>&,QMMMSettings&);

double Get_PI_Epot(vector<QMMMAtom>&,QMMMSettings&);

void GetQuotes(vector<string>&);

double HermCoul1e(HermGau&,double,Coord&);

double HermCoul2e(HermGau&,HermGau&);

double HermOverlap(HermGau&,HermGau&);

VectorXd KabschDisplacement(MatrixXd&,MatrixXd&,int);

void KabschRotation(MatrixXd&,MatrixXd&,int);

double KineticE_eFF(vector<QMMMElec>&,QMMMSettings&);

double LAMMPSEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double LAMMPSForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

double LAMMPSOpt(vector<QMMMAtom>&,QMMMSettings&,int);

void LAMMPSTopology(vector<QMMMAtom>&,stringstream&,int);

void LICHEMDFP(vector<QMMMAtom>&,QMMMSettings&,int);

void LICHEMErrorChecker(QMMMSettings&);

void LICHEMNEB(vector<QMMMAtom>&,QMMMSettings&);

void LICHEMPrintSettings(QMMMSettings&);

void LICHEMQuickMin(vector<QMMMAtom>&,QMMMSettings&,int);

void LICHEMSteepest(vector<QMMMAtom>&,QMMMSettings&,int);

void LICHEM2BASIS(int&,char**&);

void LICHEM2TINK(int&,char**&);

bool MCMove(vector<QMMMAtom>&,QMMMSettings&,double&);

VectorXd NEBTangent(VectorXd&,VectorXd&,QMMMSettings&,int);

void NWChemCharges(vector<QMMMAtom>&,QMMMSettings&,int);

double NWChemEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double NWChemForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

double NWChemOpt(vector<QMMMAtom>&,QMMMSettings&,int);

bool OptConverged(vector<QMMMAtom>&,vector<QMMMAtom>&,vector<Coord>&,
     int,QMMMSettings& QMMMOpts,int,bool);

bool PathConverged(vector<QMMMAtom>&,vector<QMMMAtom>&,MatrixXd&,
     int,QMMMSettings&,bool);

void PBCCenter(vector<QMMMAtom>&,QMMMSettings&);

void PrintFancyTitle();

void Print_traj(vector<QMMMAtom>&,fstream&,QMMMSettings&);

void PSICharges(vector<QMMMAtom>&,QMMMSettings&,int);

double PSIEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double PSIForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

double PSIOpt(vector<QMMMAtom>&,QMMMSettings&,int);

void ReadArgs(int&,char**&,fstream&,fstream&,fstream&,fstream&);

void ReadLICHEMInput(fstream&,fstream&,fstream&,
     vector<QMMMAtom>&,QMMMSettings&);

void RotateTINKCharges(vector<QMMMAtom>&,int);

OctCharges SphHarm2Charges(RedMpole);

double TINKEREnergy(vector<QMMMAtom>&,QMMMSettings&,int);

void TINKERDynamics(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKEROpt(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKERForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

double TINKERPolForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

void TINKERInduced(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKERPolEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

void TINK2LICHEM(int&,char**&);

vector<int> TraceBoundary(vector<QMMMAtom>&,int);

void WriteGauInput(vector<QMMMAtom>&,string,QMMMSettings&,int);

void WriteNWChemInput(vector<QMMMAtom>&,string,QMMMSettings&,int);

void WritePSIInput(vector<QMMMAtom>&,string,QMMMSettings&,int);

void WriteTINKMpole(vector<QMMMAtom>&,fstream&,int,int);

//Function definitions (alphabetical)
#include "Analysis.cpp"
#include "Basis.cpp"
#include "Core_funcs.cpp"
#include "Frozen_density.cpp"
#include "Hermite_eng.cpp"
#include "Input_Reader.cpp"
#include "Multipoles.cpp"
#include "Optimizers.cpp"
#include "PathIntegral.cpp"
#include "ReactionPath.cpp"
#include "StructWriter.cpp"
#include "TINK2LICHEM.cpp"

//Wrapper definitions (alphabetical)
#include "AMBER.cpp"
#include "Gaussian.cpp"
#include "LAMMPS.cpp"
#include "Lepton_eng.cpp"
#include "NWChem.cpp"
#include "PSI4.cpp"
#include "TINKER.cpp"

#endif

