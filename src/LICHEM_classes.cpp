/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Functions for the data structures used in LICHEM.

*/

//Coord class function definitions
Coord::Coord()
{
  //Constructor
  x = 0;
  y = 0;
  z = 0;
  return;
};

Coord::~Coord()
{
  //Generic destructor
  return;
};

double Coord::VecMag()
{
  //Return the squared vector magnitude
  double R2; //Squared distance
  R2 = x*x+y*y+z*z; //Calculate using local variables
  //Return value
  return R2;
};

//Mpole class function definitions
Mpole::Mpole()
{
  //Generic constructor
  return;
};

Mpole::~Mpole()
{
  //Generic destructor
  return;
};

//RedMpole class function definitions
RedMpole::RedMpole()
{
  //Generic constructor
  return;
};

RedMpole::~RedMpole()
{
  //Generic destructor
  return;
};

//OctCharges class function definitions
OctCharges::OctCharges()
{
  //Generic constructor
  return;
};

OctCharges::~OctCharges()
{
  //Generic destructor
  return;
};

//GEMDen class function definitions
GEMDen::GEMDen()
{
  //Generic constructor
  return;
};

GEMDen::GEMDen(string Typ, string BasName)
{
  //Fancy constructor
  SetBasis(Typ,BasName);
  return;
};

GEMDen::~GEMDen()
{
  //Generic destructor
  return;
};

void GEMDen::SetBasis(string Typ, string BasName)
{
  //Set the basis set
  Dens = HermBasis(Typ,BasName);
  return;
};

void GEMDen::SetFrame(bool flip, string frame, int at1, int at2, int at3)
{
  //Set the local frame of reference
  ChiralFlip = flip;
  Type = frame;
  Atom1 = at1;
  Atom2 = at2;
  Atom3 = at3;
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
  //Add the nuclear charge from the periodic table
  dmpole.q += PTable.RevTyping(Type);
  //Convert Hermite Gaussians to multipoles
  for (unsigned int i=0;i<Dens.size();i++)
  {
    //Check for a monopole
    if ((Dens[i].XPow() == 0) and (Dens[i].YPow() == 0) and
       (Dens[i].ZPow() == 0))
    {
      //Update monopole
      dmpole.q += Dens[i].Coeff();
      //Update diagonal quadrupole moments
      dmpole.Qxx += Dens[i].Coeff()/(2*Dens[i].Alpha());
      dmpole.Qyy += Dens[i].Coeff()/(2*Dens[i].Alpha());
      dmpole.Qzz += Dens[i].Coeff()/(2*Dens[i].Alpha());
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
      dmpole.Qxx += 2*Dens[i].Coeff();
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
      dmpole.Qyy += 2*Dens[i].Coeff();
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
      dmpole.Qzz += 2*Dens[i].Coeff();
    }
  }
  //Convert to a traceless quadrupole
  double QTrace = dmpole.Qxx+dmpole.Qyy+dmpole.Qzz;
  dmpole.Qxx -= QTrace;
  dmpole.Qyy -= QTrace;
  dmpole.Qzz -= QTrace;
  //Return GEM multipole
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

//PeriodicTable class function definitions
PeriodicTable::PeriodicTable()
{
  //Set atomic properties
  //NB: The widths were calculated for the all-electron density
  //Hydrogen, 1
  Typs.push_back("H");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Helium, 2
  Typs.push_back("He");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lithium, 3
  Typs.push_back("Li");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Berylium, 4
  Typs.push_back("Be");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Boron, 5
  Typs.push_back("B");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Carbon, 6
  Typs.push_back("C");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Nitrogen, 7
  Typs.push_back("N");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Oxygen, 8
  Typs.push_back("O");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Fluorine, 9
  Typs.push_back("F");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Neon, 10
  Typs.push_back("Ne");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Sodium, 11
  Typs.push_back("Na");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Magnesium, 12
  Typs.push_back("Mg");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Aluminum, 13
  Typs.push_back("Al");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Silicon, 14
  Typs.push_back("Si");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Phosphorus, 15
  Typs.push_back("P");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Sulfur, 16
  Typs.push_back("S");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Chlorine, 17
  Typs.push_back("Cl");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Argon, 18
  Typs.push_back("Ar");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Potassium, 19
  Typs.push_back("K");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Calcium, 20
  Typs.push_back("Ca");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Scandium, 21
  Typs.push_back("Sc");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Titanium, 22
  Typs.push_back("Ti");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Vanadium, 23
  Typs.push_back("V");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Chromium, 24
  Typs.push_back("Cr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Manganese, 25
  Typs.push_back("Mn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Iron, 26
  Typs.push_back("Fe");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Cobalt, 27
  Typs.push_back("Co");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Nickel, 28
  Typs.push_back("Ni");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Copper, 29
  Typs.push_back("Cu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Zinc, 30
  Typs.push_back("Zn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Gallium, 31
  Typs.push_back("Ga");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Germanium, 32
  Typs.push_back("Ge");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Arsenic, 33
  Typs.push_back("As");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Selenium, 34
  Typs.push_back("Se");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Bromine, 35
  Typs.push_back("Br");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Krypton, 36
  Typs.push_back("Kr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rubidium, 37
  Typs.push_back("Rb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Strontium, 38
  Typs.push_back("Sr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Yttrium, 39
  Typs.push_back("Y");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Zirconium, 40
  Typs.push_back("Zr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Niobium, 41
  Typs.push_back("Nb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Molybdenum, 42
  Typs.push_back("Mo");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Technetium, 43
  Typs.push_back("Tc");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Ruthenium, 44
  Typs.push_back("Ru");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rhodium, 45
  Typs.push_back("Rh");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Palladium, 46
  Typs.push_back("Pd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Silver, 47
  Typs.push_back("Ag");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Cadmium, 48
  Typs.push_back("Cd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Indium, 49
  Typs.push_back("In");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tin, 50
  Typs.push_back("Sn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Antimony, 51
  Typs.push_back("Sb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tellurium, 52
  Typs.push_back("Te");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Iodine, 53
  Typs.push_back("I");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Xenon, 54
  Typs.push_back("Xe");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Caesium, 55
  Typs.push_back("Cs");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Barium, 56
  Typs.push_back("Ba");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lanthanum, 57
  Typs.push_back("La");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Cerium, 58
  Typs.push_back("Ce");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Praseodymium, 59
  Typs.push_back("Pr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Neodymium, 60
  Typs.push_back("Nd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Promethium, 61
  Typs.push_back("Pm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Samarium, 62
  Typs.push_back("Sm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Europium, 63
  Typs.push_back("Eu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Gadolinium, 64
  Typs.push_back("Gd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Terbium, 65
  Typs.push_back("Tb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Dysprosium, 66
  Typs.push_back("Dy");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Holmium, 67
  Typs.push_back("Ho");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Erbium, 68
  Typs.push_back("Er");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Thulium, 69
  Typs.push_back("Tm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Ytterbium, 70
  Typs.push_back("Yb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lutetium, 71
  Typs.push_back("Lu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Hafnium, 72
  Typs.push_back("Hf");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tantalum, 73
  Typs.push_back("Ta");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tungsten, 74
  Typs.push_back("W");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rhenium, 75
  Typs.push_back("Re");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Osmium, 76
  Typs.push_back("Os");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Iridium, 77
  Typs.push_back("Ir");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Platinum, 78
  Typs.push_back("Pt");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Gold, 79
  Typs.push_back("Au");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Mercury, 80
  Typs.push_back("Hg");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Thallium, 81
  Typs.push_back("Tl");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lead, 82
  Typs.push_back("Pb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Bismuth, 83
  Typs.push_back("Bi");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Polonium, 84
  Typs.push_back("Po");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Astatine, 85
  Typs.push_back("At");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Radon, 86
  Typs.push_back("Rn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Francium, 87
  Typs.push_back("Fr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Radium, 88
  Typs.push_back("Ra");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Actinium, 89
  Typs.push_back("Ac");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Thorium, 90
  Typs.push_back("Th");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Protactium, 91
  Typs.push_back("Pa");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Uranium, 92
  Typs.push_back("U");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Neptunium, 93
  Typs.push_back("Np");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Plutonium, 94
  Typs.push_back("Pu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Americium, 95
  Typs.push_back("Am");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Curium, 96
  Typs.push_back("Cm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Berkelium, 97
  Typs.push_back("Bk");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Californium, 98
  Typs.push_back("Cf");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Einsteinium, 99
  Typs.push_back("Es");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Fermium, 100
  Typs.push_back("Fm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Mendelevium, 101
  Typs.push_back("Md");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Nobelium, 102
  Typs.push_back("No");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lawrencium, 103
  Typs.push_back("Lr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rutherfordium, 104
  Typs.push_back("Rf");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Dubnium, 105
  Typs.push_back("Db");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Seaborgium, 106
  Typs.push_back("Sg");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Bohrium, 107
  Typs.push_back("Bh");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Hasium, 108
  Typs.push_back("Hs");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Meitnerium, 109
  Typs.push_back("Mt");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Darmstadtium, 110
  Typs.push_back("Ds");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Roentgenium, 111
  Typs.push_back("Rg");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Copernicium, 112
  Typs.push_back("Cn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 113
  Typs.push_back("");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Flerovium, 114
  Typs.push_back("Fl");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 115
  Typs.push_back("");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Livermorium, 116
  Typs.push_back("Lv");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 117
  Typs.push_back(""); //Unnamed
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 118
  Typs.push_back(""); //Unnamed
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  return;
};

PeriodicTable::~PeriodicTable()
{
  //Generic destructor
  return;
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
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      Z = i+1;
    }
  }
  return Z;
};

double PeriodicTable::GetCovRadius(string AtName)
{
  //Function to find the covalent radius of an atom
  double radius = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      radius = CovRadii[i];
    }
  }
  return radius;
};

double PeriodicTable::GetRadius(string AtName)
{
  //Function to find the vdW radius of an atom
  double radius = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      radius = vdWRadii[i];
    }
  }
  return radius;
};

double PeriodicTable::GetAtMass(string AtName)
{
  //Function to find the atomic mass of an atom
  double mass = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      mass = AtMasses[i];
    }
  }
  return mass;
};

