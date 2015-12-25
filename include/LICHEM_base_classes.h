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
    //Radii
    vector<double> vdWRadii; //Van der Waals radii
    //Masses
    vector<double> AtMasses; //Atomic masses (amu)
  public:
    //Set data (hard coded constructor)
    PeriodicTable();
    //Destructor
    ~PeriodicTable();
    //Retrieve data
    string Typing(int); //Atom type
    int RevTyping(string); //Atomic number
    double GetGauWid(string); //Gaussian (1s) width
    double GetCovRadius(string); //Covalent radius
    double GetRadius(string); //Van der Waals radius
    double GetAtMass(string); //Atomic mass
};

class Coord
{
  public:
    //Constructor
    Coord();
    //Destructor
    ~Coord();
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
    //Constructor
    Mpole();
    //Destructor
    ~Mpole();
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
    //Constructor
    RedMpole();
    //Destructor
    ~RedMpole();
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
    //Constructor
    OctCharges();
    //Destructor
    ~OctCharges();
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

#endif

