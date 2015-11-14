/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Primary data structures used in LICHEM.

*/

//Make including safe
#ifndef LICHEM_STRUCTURES
#define LICHEM_STRUCTURES

//GEM electron density
class GEMDen
{
  //Class for GEM diffuse charge density
  private:
    //Definition of the local frame of reference
    bool ChiralFlip; //Flip y axis
    string Type; //Bisector, Z-then-X, Z-Only, 3-Fold, or Z-Bisect
    int Atom1; //Atom which defines the z axis
    int Atom2; //Atom which defines the x axis
    int Atom3; //Atom which defines the y axis (chiral only)
    //Basis functions and density
    vector<HermGau> Dens; //GEM density
  public:
    //Constructor
    GEMDen();
    //Destructor
    ~GEMDen();
    //Functions to manipulate GEM density
    Mpole GEMDM(); //Function to generate multipoles from density
};

//LICHEM data structures
class QMMMAtom
{
  //Data type for atomic information
  public:
    //Constructor
    QMMMAtom();
    //Destructor
    ~QMMMAtom();
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
    //Electron density
    vector<GEMDen> GEM; //GEM frozen density
};

class QMMMElec
{
  //Data type for electronic information (eFF model)
  public:
    //Constructor
    QMMMElec();
    //Destructor
    ~QMMMElec();
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
    //Constructor
    QMMMSettings();
    //Destructor
    ~QMMMSettings();
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

//GEMDen class function definitions
GEMDen::GEMDen()
{
  //Generic constructor
  return;
};

GEMDen::~GEMDen()
{
  //Generic destructor
  return;
};

Mpole GEMDen::GEMDM()
{
  //Function to convert GEM density to distributed multipoles
  Mpole dmpole; //Blank set of multipoles
  //Save frame of reference
  dmpole.ChiralFlip = ChiralFlip;
  dmpole.Type = Type;
  dmpole.Atom1 = Atom1;
  dmpole.Atom2 = Atom2;
  dmpole.Atom3 = Atom3;
  //Initialize multipoles
  dmpole.q = 0;
  dmpole.Dx = 0;
  dmpole.Dy = 0;
  dmpole.Dz = 0;
  dmpole.IDx = 0;
  dmpole.IDy = 0;
  dmpole.IDz = 0;
  dmpole.Qxx = 0;
  dmpole.Qxy = 0;
  dmpole.Qxz = 0;
  dmpole.Qyy = 0;
  dmpole.Qxz = 0;
  dmpole.Qzz = 0;
  //Convert Hermite Gaussians to multipoles
  for (unsigned int i=0;i<Dens.size();i++)
  {
    //Check for a monopole
    if ((Dens[i].XPow() == 0) and (Dens[i].YPow() == 0) and
       (Dens[i].ZPow() == 0))
    {
      //Update monopole
      dmpole.q += Dens[i].Coeff();
    }
    //Check for a dipole
    if ((Dens[i].XPow() == 1) and (Dens[i].YPow() == 0) and
       (Dens[i].ZPow() == 0))
    {
      //Update x dipole
      dmpole.Dx += Dens[i].Coeff();
    }
    if ((Dens[i].XPow() == 0) and (Dens[i].YPow() == 1) and
       (Dens[i].ZPow() == 0))
    {
      //Update y dipole
      dmpole.Dy += Dens[i].Coeff();
    }
    if ((Dens[i].XPow() == 0) and (Dens[i].YPow() == 0) and
       (Dens[i].ZPow() == 1))
    {
      //Update z dipole
      dmpole.Dz += Dens[i].Coeff();
    }
    //Check for a quadrupole
    if ((Dens[i].XPow() == 2) and (Dens[i].YPow() == 0) and
       (Dens[i].ZPow() == 0))
    {
      //Update xx quadrupole
      dmpole.Qxx += Dens[i].Coeff();
    }
    if ((Dens[i].XPow() == 1) and (Dens[i].YPow() == 1) and
       (Dens[i].ZPow() == 0))
    {
      //Update xy quadrupole
      dmpole.Qxy += Dens[i].Coeff();
    }
    if ((Dens[i].XPow() == 1) and (Dens[i].YPow() == 0) and
       (Dens[i].ZPow() == 1))
    {
      //Update xz quadrupole
      dmpole.Qxz += Dens[i].Coeff();
    }
    if ((Dens[i].XPow() == 0) and (Dens[i].YPow() == 2) and
       (Dens[i].ZPow() == 0))
    {
      //Update xx quadrupole
      dmpole.Qyy += Dens[i].Coeff();
    }
    if ((Dens[i].XPow() == 0) and (Dens[i].YPow() == 1) and
       (Dens[i].ZPow() == 1))
    {
      //Update yz quadrupole
      dmpole.Qyz += Dens[i].Coeff();
    }
    if ((Dens[i].XPow() == 0) and (Dens[i].YPow() == 0) and
       (Dens[i].ZPow() == 2))
    {
      //Update yz quadrupole
      dmpole.Qzz += Dens[i].Coeff();
    }
  }
  return dmpole;
};

//QMMMAtom class function definitions
QMMMAtom::QMMMAtom()
{
  //Generic constructor
  return;
};

QMMMAtom::~QMMMAtom()
{
  //Generic destructor
  return;
};

//QMMMElec class function definitions
QMMMElec::QMMMElec()
{
  //Generic constructor
  return;
};

QMMMElec::~QMMMElec()
{
  //Generic destructor
  return;
};

//QMMMSettings class function definitions
QMMMSettings::QMMMSettings()
{
  //Generic constructor
  return;
};

QMMMSettings::~QMMMSettings()
{
  //Generic destructor
  return;
};

#endif

