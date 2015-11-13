/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Base classes used in LICHEM.

*/

//Make including safe
#ifndef LICHEM_BASE_CLASSES
#define LICHEM_BASE_CLASSES

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
    //Set data (hard coded constructor)
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
    //Constructor
    Coord()
    {
      x = 0;
      y = 0;
      z = 0;
    };
    //Positions or displacements
    double x; //x position
    double y; //y position
    double z; //z position
    //Functions
    double VecMag(); //Return the squared vector magnitude
};

class Mpole
{
  //Cartesian multipoles
  public:
    //Definition of the local frame of reference
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
    //Octahedral charges
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

//Base class function definitions
double Coord::VecMag()
{
  //Return the squared vector magnitude
  double R2; //Squared distance
  R2 = x*x+y*y+z*z; //Calculate using local variables
  //Return value
  return R2;
};

//Periodic Table functions
PeriodicTable::PeriodicTable()
{
  //Set atomic properties
  //NB: The widths were calculated for the all-electron density
  Typs.push_back("H"); //Hydrogen, 1
  GauWids.push_back(0.33046422838582712); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("He"); //Helium, 2
  GauWids.push_back(0.44937774826454929); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Li"); //Lithium, 3
  GauWids.push_back(0.24891493951151511); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Be"); //Berylium, 4
  GauWids.push_back(0.28393617417828898); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("B"); //Boron, 5
  GauWids.push_back(0.31688766703304322); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("C"); //Carbon, 6
  GauWids.push_back(0.34748552120700499); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("N"); //Nitrogen, 7
  GauWids.push_back(0.38230385698895192); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("O"); //Oxygen, 8
  GauWids.push_back(0.40455289847810627); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("F"); //Fluorine, 9
  GauWids.push_back(0.43227871563688025); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ne"); //Neon, 10
  GauWids.push_back(0.45940792416326015); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Na"); //Sodium, 11
  GauWids.push_back(0.28603607723272795); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mg"); //Magnesium, 12
  GauWids.push_back(0.29359655766949572); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Al"); //Aluminum, 13
  GauWids.push_back(0.30283227512360922); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Si"); //Silicon, 14
  GauWids.push_back(0.31462038964458622); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("P"); //Phosphorus, 15
  GauWids.push_back(0.33301416291877084); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("S"); //Sulfur, 16
  GauWids.push_back(0.35056486830985728); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cl"); //Chlorine, 17
  GauWids.push_back(0.37229368179032268); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ar"); //Argon, 18
  GauWids.push_back(0.39498740756144213); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("K"); //Potassium, 19
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ca"); //Calcium, 20
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sc"); //Scandium, 21
  GauWids.push_back(0.29098058180517811); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ti"); //Titanium, 22
  GauWids.push_back(0.29703979129006947); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("V"); //Vanadium, 23
  GauWids.push_back(0.30365511638085596); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cr"); //Chromium, 24
  GauWids.push_back(0.30996489546441913); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mn"); //Manganese, 25
  GauWids.push_back(0.32461806746013033); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Fe"); //Iron, 26
  GauWids.push_back(0.33074715374383107); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Co"); //Cobalt, 27
  GauWids.push_back(0.33482470199473979); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ni"); //Nickel, 28
  GauWids.push_back(0.33983030212559556); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cu"); //Copper, 29
  GauWids.push_back(0.34432982070592233); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Zn"); //Zinc, 30
  GauWids.push_back(0.34411081636918234); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ga"); //Gallium, 31
  GauWids.push_back(0.33702898576385676); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ge"); //Germanium, 32
  GauWids.push_back(0.33935669473845503); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("As"); //Arsenic, 33
  GauWids.push_back(0.34818738929821641); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Se"); //Selenium, 34
  GauWids.push_back(0.35888462095130624); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Br"); //Bromine, 35
  GauWids.push_back(0.37243540751493021); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Kr"); //Krypton, 36
  GauWids.push_back(0.3880309254592651); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rb"); //Rubidium, 37
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sr"); //Strontium, 38
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Y"); //Yttrium, 39
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Zr"); //Zirconium, 40
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Nb"); //Niobium, 41
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mo"); //Molybdenum, 42
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tc"); //Technetium, 43
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ru"); //Ruthenium, 44
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rh"); //Rhodium, 45
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pd"); //Palladium, 46
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ag"); //Silver, 47
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cd"); //Cadmium, 48
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("In"); //Indium, 49
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sn"); //Tin, 50
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sb"); //Antimony, 51
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Te"); //Tellurium, 52
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("I"); //Iodine, 53
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Xe"); //Xenon, 54
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cs"); //Caesium, 55
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ba"); //Barium, 56
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("La"); //Lanthanum, 57
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ce"); //Cerium, 58
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pr"); //Praseodymium, 59
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Nd"); //Neodymium, 60
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pm"); //Promethium, 61
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sm"); //Samarium, 62
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Eu"); //Europium, 63
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Gd"); //Gadolinium, 64
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tb"); //Terbium, 65
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Dy"); //Dysprosium, 66
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ho"); //Holmium, 67
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Er"); //Erbium, 68
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tm"); //Thulium, 69
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Yb"); //Ytterbium, 70
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Lu"); //Lutetium, 71
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Hf"); //Hafnium, 72
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ta"); //Tantalum, 73
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("W"); //Tungsten, 74
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Re"); //Rhenium, 75
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Os"); //Osmium, 76
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ir"); //Iridium, 77
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pt"); //Platinum, 78
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Au"); //Gold, 79
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Hg"); //Mercury, 80
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tl"); //Thallium, 81
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pb"); //Lead, 82
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Bi"); //Bismuth, 83
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Po"); //Polonium, 84
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("At"); //Astatine, 85
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rn"); //Radon, 86
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Fr"); //Francium, 87
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ra"); //Radium, 88
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ac"); //Actinium, 89
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Th"); //Thorium, 90
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pa"); //Protactium, 91
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("U"); //Uranium, 92
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Np"); //Neptunium, 93
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pu"); //Plutonium, 94
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Am"); //Americium, 95
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cm"); //Curium, 96
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Bk"); //Berkelium, 97
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cf"); //Californium, 98
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Es"); //Einsteinium, 99
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Fm"); //Fermium, 100
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Md"); //Mendelevium, 101
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("No"); //Nobelium, 102
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Lr"); //Lawrencium, 103
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rf"); //Rutherfordium, 104
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Db"); //Dubnium, 105
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sg"); //Seaborgium, 106
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Bh"); //Bohrium, 107
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Hs"); //Hasium, 108
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mt"); //Meitnerium, 109
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ds"); //Darmstadtium, 110
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rg"); //Roentgenium, 111
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cn"); //Copernicium, 112
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back(""); //Num. 113
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Fl"); //Flerovium, 114
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back(""); //Num. 115
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Lv"); //Livermorium, 116
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back(""); //Num. 117
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back(""); //Num. 118
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
};

string PeriodicTable::Typing(int Z)
{
  //Function to convert nuclear charges to atom types
  return Typs[Z-1];
};

int PeriodicTable::RevTyping(string AtName)
{
  //Function to convert atom types to nuclear charges
  int Z = 0;
  #pragma omp parallel for num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      Z = i+1;
    }
  }
  #pragma omp barrier
  return Z;
};

double PeriodicTable::GetGauWid(string AtName)
{
  //Function to find the 1s Gaussian width of an atom
  double gauwid = 0;
  #pragma omp parallel for num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      gauwid = GauWids[i];
    }
  }
  #pragma omp barrier
  return gauwid;
};

double PeriodicTable::GetCovRadius(string AtName)
{
  //Function to find the 1s Gaussian width of an atom
  double radius = 0;
  #pragma omp parallel for num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      radius = CovRadii[i];
    }
  }
  #pragma omp barrier
  return radius;
};

#endif

