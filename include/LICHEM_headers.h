/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Primary header for LICHEM.

*/

//Make including safe
#ifndef LICHEM_HEADERS
#define LICHEM_HEADERS

/*!
  \defgroup  Analysis Analysis
  \brief Routines for analyzing structures and trajectories
*/

/*!
  \defgroup Classes Classes
  \brief Class definitions
*/

/*!
  \defgroup FrznDen FrzDen
  \brief Routines for manipulating frozen electron density
*/

/*!
  \defgroup Gaussian Gaussian
  \brief Wrapper functions for Gaussian calculations
*/

/*!
  \defgroup LAMMPS LAMMPS
  \brief Wrapper functions for LAMMPS calculations
*/

/*!
  \defgroup Main Main
  \brief LICHEM
*/

/*!
  \defgroup Misc Misc
  \brief Misc. utility functions
*/

/*!
  \defgroup MonteCarlo MonteCarlo
  \brief Path-integral Monte Carlo routines
*/

/*!
  \defgroup Multipoles Multipoles
  \brief Routines for manipulating atomic multipole moments
*/

/*!
  \defgroup NWChem NWChem
  \brief Wrapper functions for NWChem calculations
*/

/*!
  \defgroup Optimizers Optimizers
  \brief Optimization routines
*/

/*!
  \defgroup PSI4 PSI4
  \brief Wrapper functions for psi4 calculations
*/

/*!
  \defgroup tests tests
  \brief LICHEM semi-automated test suite
*/

/*!
  \defgroup TINKER TINKER
  \brief Wrapper functions for TINKER calculations
*/

//C++ libraries
#include "LICHEM_clibs.h"
using namespace std;

//Matrix libraries
#include "LICHEM_matrix.h"
using namespace Eigen;

//Predefined constants, variables, and flags
#include "LICHEM_options.h"
using namespace LICHEMOpts;
#include "LICHEM_constants.h"
using namespace LICHEMConst;
#include "LICHEM_globals.h"
using namespace LICHEMGlobal;

//LICHEM headers and libraries
#include "LICHEM_base_classes.h"
#include "LICHEM_Hermite.h"
#include "LICHEM_classes.h"

//Set up periodic table
PeriodicTable chemTable;

//Function declarations (alphabetical)
bool Angled(vector<QMMMAtom>&,int,int);

double Bohring(double);

bool Bonded(vector<QMMMAtom>&,int,int);

double BoysFunc(int,double);

void BurstTraj(vector<QMMMAtom>&,QMMMSettings&);

RedMPole Cart2SphHarm(MPole&);

VectorXd CINEBTangent(VectorXd&,VectorXd&,QMMMSettings&,int);

bool CheckFile(const string&);

void CheckNEBTangent(VectorXd&);

Coord CoordDist2(Coord&,Coord&);

bool Dihedraled(vector<QMMMAtom>&,int,int);

void ExternalGaussian(int&,char**&);

void ExtractGlobalPoles(int& argc, char**& argv);

void ExtractTINKpoles(vector<QMMMAtom>&,int);

int FBNEBMCMove(vector<QMMMAtom>&,vector<VectorXd>&,QMMMSettings&,
                VectorXd&);

void FetchQuotes(vector<string>&);

int FindMaxThreads();

Coord FindQMCOM(vector<QMMMAtom>&,QMMMSettings&,int);

void FindTINKERClasses(vector<QMMMAtom>&);

void GaussianCharges(vector<QMMMAtom>&,QMMMSettings&,int);

double GaussianEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double GaussianExternOpt(vector<QMMMAtom>&,QMMMSettings&,int);

double GaussianForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

MatrixXd GaussianHessian(vector<QMMMAtom>&,QMMMSettings&,int);

double GaussianOpt(vector<QMMMAtom>&,QMMMSettings&,int);

double GEMBuffC7(double,double,Coord&,Coord&,double);

double GEMC6(double,Coord&,Coord&,double);

double Get_PI_Espring(vector<QMMMAtom>&,QMMMSettings&);

double Get_PI_Epot(vector<QMMMAtom>&,QMMMSettings&);

vector<HermGau> HermBasis(string,string);

double HermCoul1e(HermGau&,double,Coord&);

double HermCoul2e(HermGau&,HermGau&);

double HermOverlap(HermGau&,HermGau&);

VectorXd KabschDisplacement(MatrixXd&,MatrixXd&,int);

void KabschRotation(MatrixXd&,MatrixXd&,int);

double LAMMPSEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double LAMMPSForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

MatrixXd LAMMPSHessian(vector<QMMMAtom>&,QMMMSettings&,int);

double LAMMPSOpt(vector<QMMMAtom>&,QMMMSettings&,int);

void LAMMPSTopology(vector<QMMMAtom>&,stringstream&,int);

void LICHEM2BASIS(int&,char**&);

void LICHEM2TINK(int&,char**&);

template<typename T> int LICHEMCount(T);

double LICHEMDensity(vector<QMMMAtom>&,QMMMSettings&);

void LICHEMDFP(vector<QMMMAtom>&,QMMMSettings&,int);

void LICHEMErrorChecker(QMMMSettings&);

double LICHEMFactorial(int);

void LICHEMFixSciNot(string&);

template<typename T> string LICHEMFormFloat(T,int);

VectorXd LICHEMFreq(vector<QMMMAtom>&,MatrixXd&,QMMMSettings&,int,int&);

void LICHEMLowerText(string&);

void LICHEMNEB(vector<QMMMAtom>&,QMMMSettings&,int);

void LICHEMPrintSettings(vector<QMMMAtom>&,QMMMSettings&);

void LICHEMSteepest(vector<QMMMAtom>&,QMMMSettings&,int);

void LICHEMUpperText(string&);

double LRECFunction(Coord&,QMMMSettings&);

bool MCMove(vector<QMMMAtom>&,QMMMSettings&,double&);

VectorXd NEBTangent(VectorXd&,VectorXd&,QMMMSettings&,int);

void NWChemCharges(vector<QMMMAtom>&,QMMMSettings&,int);

double NWChemEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double NWChemForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

MatrixXd NWChemHessian(vector<QMMMAtom>&,QMMMSettings&,int);

double NWChemOpt(vector<QMMMAtom>&,QMMMSettings&,int);

bool OptConverged(vector<QMMMAtom>&,vector<QMMMAtom>&,VectorXd&,
                  int,QMMMSettings&,int,bool);

bool PathConverged(vector<QMMMAtom>&,vector<QMMMAtom>&,MatrixXd&,
                   int,QMMMSettings&,bool);

void PathLinInterpolate(int&,char**&);

void PBCCenter(vector<QMMMAtom>&,QMMMSettings&);

void PrintFancyTitle();

void PrintLapin();

void Print_traj(vector<QMMMAtom>&,fstream&,QMMMSettings&);

void PSI4Charges(vector<QMMMAtom>&,QMMMSettings&,int);

double PSI4Energy(vector<QMMMAtom>&,QMMMSettings&,int);

double PSI4Forces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

MatrixXd PSI4Hessian(vector<QMMMAtom>&,QMMMSettings&,int);

double PSI4Opt(vector<QMMMAtom>&,QMMMSettings&,int);

void ReadArgs(int&,char**&,fstream&,fstream&,fstream&,fstream&);

void ReadLICHEMInput(fstream&,fstream&,fstream&,
                     vector<QMMMAtom>&,QMMMSettings&);

void RotateTINKCharges(vector<QMMMAtom>&,int);

OctCharges SphHarm2Charges(RedMPole);

void SplitPathTraj(int&,char**&);

void TINK2LICHEM(int&,char**&);

double TINKEREnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKERForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

MatrixXd TINKERHessian(vector<QMMMAtom>&,QMMMSettings&,int);

void TINKERInduced(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKERMMForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

double TINKEROpt(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKERPolEnergy(vector<QMMMAtom>&,QMMMSettings&,int);

double TINKERPolForces(vector<QMMMAtom>&,VectorXd&,QMMMSettings&,int);

vector<int> TraceBoundary(vector<QMMMAtom>&,int);

void WriteChargeFile(vector<QMMMAtom>&,QMMMSettings&,int);

void WriteGauInput(vector<QMMMAtom>&,string,QMMMSettings&,int);

void WriteNWChemInput(vector<QMMMAtom>&,string,QMMMSettings&,int);

void WriteModes(vector<QMMMAtom>&,bool,VectorXd&,MatrixXd&,QMMMSettings&,int);

void WritePSI4Input(vector<QMMMAtom>&,string,QMMMSettings&,int);

void WriteTINKMPole(vector<QMMMAtom>&,fstream&,int,int);

void WriteQMConnect(int&,char**&);

//Function definitions (alphabetical)
#include "Analysis.cpp"
#include "Basis.cpp"
#include "Basis_sets.cpp"
#include "Core_funcs.cpp"
#include "Frozen_density.cpp"
#include "Hermite_eng.cpp"
#include "Input_reader.cpp"
#include "LAMMPS2LICHEM.cpp"
#include "LICHEM_classes.cpp"
#include "Multipoles.cpp"
#include "Optimizers.cpp"
#include "Path_integral.cpp"
#include "Reaction_path.cpp"
#include "Struct_writer.cpp"
#include "Text_format.cpp"
#include "TINK2LICHEM.cpp"

//Wrapper definitions (alphabetical)
#include "GauExternal.cpp"
#include "Gaussian.cpp"
#include "LAMMPS.cpp"
#include "NWChem.cpp"
#include "PSI4.cpp"
#include "TINKER.cpp"

#endif

