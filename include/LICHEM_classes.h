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
    bool chiralFlip_; //Flip y axis
    string type_; //Bisector, Z-then-X, Z-Only, 3-Fold, or Z-Bisect
    int atom1_; //Atom which defines the z axis
    int atom2_; //Atom which defines the x axis
    int atom3_; //Atom which defines the y axis (chiral only)
    //Basis functions and density
    vector<HermGau> dens_; //GEM density
  public:
    //Constructors
    GEMDen();
    GEMDen(string,string);
    //Destructor
    ~GEMDen();
    //Functions to manipulate GEM density
    void setBasis(string,string); //Sets the Hermite basis
    void setFrame(bool,string,int,int,int); //Sets the frame of reference
    MPole GEMDM(); //Function to generate multipoles from density
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
    bool NEBActive; //Included in NEB tangent calculations
    bool QMRegion; //QM, MM, pseudo-bond, or boundary-atom
    bool MMRegion; //QM, MM, pseudo-bond, or boundary-atom
    bool PBRegion; //QM, MM, pseudo-bond, or boundary-atom
    bool BARegion; //QM, MM, pseudo-bond, or boundary-atom
    bool frozen; //Part of a frozen shell
    //Force field information
    double m; //Mass of atom
    string QMTyp; //Real atom type
    string MMTyp; //Force field atom type
    int numTyp; //Numerical atom type (if used)
    int numClass; //Numerical atom class (if used)
    int id; //Atom number, starts at zero
    vector<int> bonds; //Connectivity
    //Coordinates
    vector<Coord> P; //Array of beads
    //Multipoles
    vector<MPole> MP; //Multipoles
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
    string func; //QM method (functional, HF, etc)
    string basis; //Basis set for QM calculations
    int RAM; //Ram for QM calculations
    bool memMB; //Is the RAM in mb or gb
    int charge; //QM total charge
    int spin; //QM total spin
    string unitsQM; //Specifies the units for the QM calculations
    string backDir; //Directory for log file backups
    //Input needed for QMMM long-range electrostatics
    bool useLREC; //Use a long-range correction
    double LRECCut; //Cutoff for the long-range correction
    int LRECPow; //Exponent for the LREC smoothing function
    //Input needed for MM wrappers
    bool useMMCut; //Flag to turn the cutoff on or off
    double MMOptCut; //Electrostatic cutoff for MM optimzations (Ang)
    bool useEwald; //Use Ewald summation for MM energy and optimizations
    bool useImpSolv; //Use implicit solvents for MM energy and optimizations
    string solvModel; //Type of implicit solvent
    //Input needed for MC and reaction path functions
    string ensemble; //NVT or NPT
    double temp; //Temperature
    double beta; //Inverse temperature
    double press; //External pressure
    int NEq; //Number of equilibration run steps
    int NSteps; //Number of production run steps
    int NBeads; //Number of time-slices or beads
    double accRatio; //Target acceptance ratio
    int NPrint; //Number of steps before printing
    //Input needed for optimizations
    int maxOptSteps; //Maximum iterative optimization steps
    double MMOptTol; //Criteria to end the optimization
    double QMOptTol; //Criteria to end the optimization
    double stepScale; //Steepest descent step size (Ang)
    double maxStep; //Maximum size of the optimization step
    //Input needed for reaction paths
    double kSpring; //Elastic band spring constant
    int TSBead; //Current guess of the transition state
    bool climb; //Flag to turn on climbing image NEB
    bool frznEnds; //Flag to freeze the NEB end points
    bool NEBFreq; //Flag to calculate TS frequencies after NEB
    bool printNormModes; //Print normal modes for pure QM calculations
    bool startPathChk; //Flag to initially use checkpoints from nearby beads
    //Storage of energies (NEB and PIMC)
    double EOld; //Temporary storage
    double EReact; //Reactant energy
    double EProd; //Product energy
    double ETrans; //Transition state energy
};

#endif

