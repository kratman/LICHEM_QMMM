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

//LICHEM data structures
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

#endif

