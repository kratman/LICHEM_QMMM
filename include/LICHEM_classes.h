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
    //Constructors
    GEMDen();
    GEMDen(string,string);
    //Destructor
    ~GEMDen();
    //Functions to manipulate GEM density
    void SetBasis(string,string); //Sets the Hermite basis
    void SetFrame(bool,string,int,int,int); //Sets the frame of reference
    Mpole GEMDM(); //Function to generate multipoles from density
};

//LICHEM particle data structures
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

//LICHEM simulation data
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
    int RAM; //Ram for QM calculations
    bool MemMB; //Is the RAM in mb or gb
    int Charge; //QM total charge
    int Spin; //QM total spin
    string BackDir; //Directory for log file backups
    bool UseLREC; //Use a long-range correction
    double LRECCut; //Cutoff for the long-range correction
    //Input needed for MM wrappers
    bool UseMMCut; //Flag to turn the cutoff on or off
    double MMOptCut; //Electrostatic cutoff for MM optimzations (Ang)
    bool UseEwald; //Use Ewald summation for MM energy and optimizations
    bool UseImpSolv; //Use implicit solvents for MM energy and optimizations
    string SolvModel; //Type of implicit solvent
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
    bool NEBFreq; //Flag to calculate TS frequencies after NEB
    bool PrintNormModes; //Print normal modes for pure QM calculations
    bool StartPathChk; //Flag to initially use checkpoints from nearby beads
    //Storage of energies (NEB and PIMC)
    double Eold; //Temporary storage
    double Ereact; //Reactant energy
    double Eprod; //Product energy
    double Ets; //Transition state energy
};

#endif

