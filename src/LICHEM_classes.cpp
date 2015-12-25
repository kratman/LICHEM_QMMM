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
  GauWids.push_back(0.33046422838582712); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Helium, 2
  Typs.push_back("He");
  GauWids.push_back(0.44937774826454929); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lithium, 3
  Typs.push_back("Li");
  GauWids.push_back(0.24891493951151511); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Berylium, 4
  Typs.push_back("Be");
  GauWids.push_back(0.28393617417828898); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Boron, 5
  Typs.push_back("B");
  GauWids.push_back(0.31688766703304322); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Carbon, 6
  Typs.push_back("C");
  GauWids.push_back(0.34748552120700499); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Nitrogen, 7
  Typs.push_back("N");
  GauWids.push_back(0.38230385698895192); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Oxygen, 8
  Typs.push_back("O");
  GauWids.push_back(0.40455289847810627); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Fluorine, 9
  Typs.push_back("F");
  GauWids.push_back(0.43227871563688025); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Neon, 10
  Typs.push_back("Ne");
  GauWids.push_back(0.45940792416326015); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Sodium, 11
  Typs.push_back("Na");
  GauWids.push_back(0.28603607723272795); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Magnesium, 12
  Typs.push_back("Mg");
  GauWids.push_back(0.29359655766949572); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Aluminum, 13
  Typs.push_back("Al");
  GauWids.push_back(0.30283227512360922); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Silicon, 14
  Typs.push_back("Si");
  GauWids.push_back(0.31462038964458622); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Phosphorus, 15
  Typs.push_back("P");
  GauWids.push_back(0.33301416291877084); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Sulfur, 16
  Typs.push_back("S");
  GauWids.push_back(0.35056486830985728); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Chlorine, 17
  Typs.push_back("Cl");
  GauWids.push_back(0.37229368179032268); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Argon, 18
  Typs.push_back("Ar");
  GauWids.push_back(0.39498740756144213); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Potassium, 19
  Typs.push_back("K");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Calcium, 20
  Typs.push_back("Ca");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Scandium, 21
  Typs.push_back("Sc");
  GauWids.push_back(0.29098058180517811); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Titanium, 22
  Typs.push_back("Ti");
  GauWids.push_back(0.29703979129006947); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Vanadium, 23
  Typs.push_back("V");
  GauWids.push_back(0.30365511638085596); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Chromium, 24
  Typs.push_back("Cr");
  GauWids.push_back(0.30996489546441913); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Manganese, 25
  Typs.push_back("Mn");
  GauWids.push_back(0.32461806746013033); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Iron, 26
  Typs.push_back("Fe");
  GauWids.push_back(0.33074715374383107); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Cobalt, 27
  Typs.push_back("Co");
  GauWids.push_back(0.33482470199473979); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Nickel, 28
  Typs.push_back("Ni");
  GauWids.push_back(0.33983030212559556); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Copper, 29
  Typs.push_back("Cu");
  GauWids.push_back(0.34432982070592233); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Zinc, 30
  Typs.push_back("Zn");
  GauWids.push_back(0.34411081636918234); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Gallium, 31
  Typs.push_back("Ga");
  GauWids.push_back(0.33702898576385676); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Germanium, 32
  Typs.push_back("Ge");
  GauWids.push_back(0.33935669473845503); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Arsenic, 33
  Typs.push_back("As");
  GauWids.push_back(0.34818738929821641); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Selenium, 34
  Typs.push_back("Se");
  GauWids.push_back(0.35888462095130624); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Bromine, 35
  Typs.push_back("Br");
  GauWids.push_back(0.37243540751493021); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Krypton, 36
  Typs.push_back("Kr");
  GauWids.push_back(0.3880309254592651); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rubidium, 37
  Typs.push_back("Rb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Strontium, 38
  Typs.push_back("Sr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Yttrium, 39
  Typs.push_back("Y");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Zirconium, 40
  Typs.push_back("Zr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Niobium, 41
  Typs.push_back("Nb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Molybdenum, 42
  Typs.push_back("Mo");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Technetium, 43
  Typs.push_back("Tc");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Ruthenium, 44
  Typs.push_back("Ru");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rhodium, 45
  Typs.push_back("Rh");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Palladium, 46
  Typs.push_back("Pd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Silver, 47
  Typs.push_back("Ag");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Cadmium, 48
  Typs.push_back("Cd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Indium, 49
  Typs.push_back("In");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tin, 50
  Typs.push_back("Sn");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Antimony, 51
  Typs.push_back("Sb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tellurium, 52
  Typs.push_back("Te");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Iodine, 53
  Typs.push_back("I");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Xenon, 54
  Typs.push_back("Xe");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Caesium, 55
  Typs.push_back("Cs");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Barium, 56
  Typs.push_back("Ba");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lanthanum, 57
  Typs.push_back("La");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Cerium, 58
  Typs.push_back("Ce");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Praseodymium, 59
  Typs.push_back("Pr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Neodymium, 60
  Typs.push_back("Nd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Promethium, 61
  Typs.push_back("Pm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Samarium, 62
  Typs.push_back("Sm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Europium, 63
  Typs.push_back("Eu");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Gadolinium, 64
  Typs.push_back("Gd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Terbium, 65
  Typs.push_back("Tb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Dysprosium, 66
  Typs.push_back("Dy");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Holmium, 67
  Typs.push_back("Ho");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Erbium, 68
  Typs.push_back("Er");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Thulium, 69
  Typs.push_back("Tm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Ytterbium, 70
  Typs.push_back("Yb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lutetium, 71
  Typs.push_back("Lu");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Hafnium, 72
  Typs.push_back("Hf");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tantalum, 73
  Typs.push_back("Ta");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Tungsten, 74
  Typs.push_back("W");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rhenium, 75
  Typs.push_back("Re");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Osmium, 76
  Typs.push_back("Os");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Iridium, 77
  Typs.push_back("Ir");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Platinum, 78
  Typs.push_back("Pt");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Gold, 79
  Typs.push_back("Au");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Mercury, 80
  Typs.push_back("Hg");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Thallium, 81
  Typs.push_back("Tl");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lead, 82
  Typs.push_back("Pb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Bismuth, 83
  Typs.push_back("Bi");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Polonium, 84
  Typs.push_back("Po");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Astatine, 85
  Typs.push_back("At");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Radon, 86
  Typs.push_back("Rn");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Francium, 87
  Typs.push_back("Fr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Radium, 88
  Typs.push_back("Ra");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Actinium, 89
  Typs.push_back("Ac");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Thorium, 90
  Typs.push_back("Th");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Protactium, 91
  Typs.push_back("Pa");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Uranium, 92
  Typs.push_back("U");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Neptunium, 93
  Typs.push_back("Np");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Plutonium, 94
  Typs.push_back("Pu");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Americium, 95
  Typs.push_back("Am");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Curium, 96
  Typs.push_back("Cm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Berkelium, 97
  Typs.push_back("Bk");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Californium, 98
  Typs.push_back("Cf");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Einsteinium, 99
  Typs.push_back("Es");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Fermium, 100
  Typs.push_back("Fm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Mendelevium, 101
  Typs.push_back("Md");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Nobelium, 102
  Typs.push_back("No");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Lawrencium, 103
  Typs.push_back("Lr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Rutherfordium, 104
  Typs.push_back("Rf");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Dubnium, 105
  Typs.push_back("Db");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Seaborgium, 106
  Typs.push_back("Sg");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Bohrium, 107
  Typs.push_back("Bh");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Hasium, 108
  Typs.push_back("Hs");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Meitnerium, 109
  Typs.push_back("Mt");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Darmstadtium, 110
  Typs.push_back("Ds");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Roentgenium, 111
  Typs.push_back("Rg");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Copernicium, 112
  Typs.push_back("Cn");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 113
  Typs.push_back("");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Flerovium, 114
  Typs.push_back("Fl");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 115
  Typs.push_back("");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Livermorium, 116
  Typs.push_back("Lv");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 117
  Typs.push_back(""); //Unnamed
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.0); //Default value
  //Num. 118
  Typs.push_back(""); //Unnamed
  GauWids.push_back(0.3); //Default value
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

double PeriodicTable::GetGauWid(string AtName)
{
  //Function to find the 1s Gaussian width of an atom
  double gauwid = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      gauwid = GauWids[i];
    }
  }
  return gauwid;
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

