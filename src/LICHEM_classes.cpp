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
  //QM wrapper settings
  Func = "N/A";
  Basis = "N/A";
  RAM = 256;
  MemMB = 1;
  Charge = 0;
  Spin = 1;
  BackDir = "Old_files";
  UseLREC = 0;
  LRECCut = 1000.0; //Effectively infinite
  //MM wrapper settings
  UseMMCut = 0;
  MMOptCut = 1000.0; //Effectively infinite
  UseEwald = 0;
  UseImpSolv = 0;
  SolvModel = "N/A";
  //MC, MD, and RP settings
  Ensemble = "N/A";
  Temp = 300.0;
  Beta = 1/(300.0*k);
  Press = 0.0;
  Neq = 0;
  Nsteps = 0;
  Nbeads = 1; //Key for printing
  accratio = 0.5;
  Nprint = 5000;
  dt = 1.0;
  tautemp = 1000.0;
  //Optimization settings
  MaxOptSteps = 200;
  MMOptTol = 1e-2;
  QMOptTol = 5e-4;
  StepScale = 1.0;
  MaxStep = 0.1;
  //Additional RP settings
  Kspring = 1.0;
  TSBead = 0;
  Climb = 0;
  FrznEnds = 0;
  NEBFreq = 0;
  StartPathChk = 1; //Speeds up reaction pathways
  //Temporary energy storage
  Eold = 0.0;
  Ereact = 0.0;
  Eprod = 0.0;
  Ets = 0.0;
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
  //Hydrogen, 1
  Typs.push_back("H");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(1.00784); //NIST, 2015
  //Helium, 2
  Typs.push_back("He");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(4.002602); //NIST, 2015
  //Lithium, 3
  Typs.push_back("Li");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(6.938); //NIST, 2015
  //Berylium, 4
  Typs.push_back("Be");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(9.0121831); //NIST, 2015
  //Boron, 5
  Typs.push_back("B");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(10.806); //NIST, 2015
  //Carbon, 6
  Typs.push_back("C");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(12.0096); //NIST, 2015
  //Nitrogen, 7
  Typs.push_back("N");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(14.00643); //NIST, 2015
  //Oxygen, 8
  Typs.push_back("O");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(15.99903); //NIST, 2015
  //Fluorine, 9
  Typs.push_back("F");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(18.998403163); //NIST, 2015
  //Neon, 10
  Typs.push_back("Ne");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(20.1797); //NIST, 2015
  //Sodium, 11
  Typs.push_back("Na");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(22.98976928); //NIST, 2015
  //Magnesium, 12
  Typs.push_back("Mg");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(24.304); //NIST, 2015
  //Aluminum, 13
  Typs.push_back("Al");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(26.9815385); //NIST, 2015
  //Silicon, 14
  Typs.push_back("Si");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(28.084); //NIST, 2015
  //Phosphorus, 15
  Typs.push_back("P");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(30.973761998); //NIST, 2015
  //Sulfur, 16
  Typs.push_back("S");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(32.059); //NIST, 2015
  //Chlorine, 17
  Typs.push_back("Cl");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(35.446); //NIST, 2015
  //Argon, 18
  Typs.push_back("Ar");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(39.948); //NIST, 2015
  //Potassium, 19
  Typs.push_back("K");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(39.0983); //NIST, 2015
  //Calcium, 20
  Typs.push_back("Ca");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(40.078); //NIST, 2015
  //Scandium, 21
  Typs.push_back("Sc");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(44.955908); //NIST, 2015
  //Titanium, 22
  Typs.push_back("Ti");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(47.867); //NIST, 2015
  //Vanadium, 23
  Typs.push_back("V");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(50.9415); //NIST, 2015
  //Chromium, 24
  Typs.push_back("Cr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(51.9961); //NIST, 2015
  //Manganese, 25
  Typs.push_back("Mn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(54.938044); //NIST, 2015
  //Iron, 26
  Typs.push_back("Fe");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(55.845); //NIST, 2015
  //Cobalt, 27
  Typs.push_back("Co");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(58.933194); //NIST, 2015
  //Nickel, 28
  Typs.push_back("Ni");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(58.6934); //NIST, 2015
  //Copper, 29
  Typs.push_back("Cu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(63.546); //NIST, 2015
  //Zinc, 30
  Typs.push_back("Zn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(65.38); //NIST, 2015
  //Gallium, 31
  Typs.push_back("Ga");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(69.723); //NIST, 2015
  //Germanium, 32
  Typs.push_back("Ge");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(72.630); //NIST, 2015
  //Arsenic, 33
  Typs.push_back("As");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(74.921595); //NIST, 2015
  //Selenium, 34
  Typs.push_back("Se");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(78.971); //NIST, 2015
  //Bromine, 35
  Typs.push_back("Br");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(79.901); //NIST, 2015
  //Krypton, 36
  Typs.push_back("Kr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(83.798); //NIST, 2015
  //Rubidium, 37
  Typs.push_back("Rb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(85.4678); //NIST, 2015
  //Strontium, 38
  Typs.push_back("Sr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(87.62); //NIST, 2015
  //Yttrium, 39
  Typs.push_back("Y");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(88.90584); //NIST, 2015
  //Zirconium, 40
  Typs.push_back("Zr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(91.224); //NIST, 2015
  //Niobium, 41
  Typs.push_back("Nb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(92.90637); //NIST, 2015
  //Molybdenum, 42
  Typs.push_back("Mo");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(95.95); //NIST, 2015
  //Technetium, 43
  Typs.push_back("Tc");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(98); //NIST, 2015
  //Ruthenium, 44
  Typs.push_back("Ru");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(101.07); //NIST, 2015
  //Rhodium, 45
  Typs.push_back("Rh");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(102.90550); //NIST, 2015
  //Palladium, 46
  Typs.push_back("Pd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(106.42); //NIST, 2015
  //Silver, 47
  Typs.push_back("Ag");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(107.8682); //NIST, 2015
  //Cadmium, 48
  Typs.push_back("Cd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(112.414); //NIST, 2015
  //Indium, 49
  Typs.push_back("In");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(114.818); //NIST, 2015
  //Tin, 50
  Typs.push_back("Sn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(118.710); //NIST, 2015
  //Antimony, 51
  Typs.push_back("Sb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(121.760); //NIST, 2015
  //Tellurium, 52
  Typs.push_back("Te");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(127.60); //NIST, 2015
  //Iodine, 53
  Typs.push_back("I");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(126.90447); //NIST, 2015
  //Xenon, 54
  Typs.push_back("Xe");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(131.293); //NIST, 2015
  //Caesium, 55
  Typs.push_back("Cs");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(132.90545196); //NIST, 2015
  //Barium, 56
  Typs.push_back("Ba");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(137.327); //NIST, 2015
  //Lanthanum, 57
  Typs.push_back("La");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(138.90547); //NIST, 2015
  //Cerium, 58
  Typs.push_back("Ce");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(140.116); //NIST, 2015
  //Praseodymium, 59
  Typs.push_back("Pr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(140.90766); //NIST, 2015
  //Neodymium, 60
  Typs.push_back("Nd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(144.242); //NIST, 2015
  //Promethium, 61
  Typs.push_back("Pm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(145); //NIST, 2015
  //Samarium, 62
  Typs.push_back("Sm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(150.36); //NIST, 2015
  //Europium, 63
  Typs.push_back("Eu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(151.964); //NIST, 2015
  //Gadolinium, 64
  Typs.push_back("Gd");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(157.25); //NIST, 2015
  //Terbium, 65
  Typs.push_back("Tb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(158.92535); //NIST, 2015
  //Dysprosium, 66
  Typs.push_back("Dy");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(162.500); //NIST, 2015
  //Holmium, 67
  Typs.push_back("Ho");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(164.93033); //NIST, 2015
  //Erbium, 68
  Typs.push_back("Er");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(167.259); //NIST, 2015
  //Thulium, 69
  Typs.push_back("Tm");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(168.93422); //NIST, 2015
  //Ytterbium, 70
  Typs.push_back("Yb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(173.054); //NIST, 2015
  //Lutetium, 71
  Typs.push_back("Lu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(174.9668); //NIST, 2015
  //Hafnium, 72
  Typs.push_back("Hf");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(178.49); //NIST, 2015
  //Tantalum, 73
  Typs.push_back("Ta");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(180.94788); //NIST, 2015
  //Tungsten, 74
  Typs.push_back("W");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(183.84); //NIST, 2015
  //Rhenium, 75
  Typs.push_back("Re");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(186.207); //NIST, 2015
  //Osmium, 76
  Typs.push_back("Os");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(190.23); //NIST, 2015
  //Iridium, 77
  Typs.push_back("Ir");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(192.217); //NIST, 2015
  //Platinum, 78
  Typs.push_back("Pt");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(195.084); //NIST, 2015
  //Gold, 79
  Typs.push_back("Au");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(196.966569); //NIST, 2015
  //Mercury, 80
  Typs.push_back("Hg");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(200.592); //NIST, 2015
  //Thallium, 81
  Typs.push_back("Tl");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(204.382); //NIST, 2015
  //Lead, 82
  Typs.push_back("Pb");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(207.2); //NIST, 2015
  //Bismuth, 83
  Typs.push_back("Bi");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(208.98040); //NIST, 2015
  //Polonium, 84
  Typs.push_back("Po");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(209); //NIST, 2015
  //Astatine, 85
  Typs.push_back("At");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(210); //NIST, 2015
  //Radon, 86
  Typs.push_back("Rn");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(222); //NIST, 2015
  //Francium, 87
  Typs.push_back("Fr");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(223); //NIST, 2015
  //Radium, 88
  Typs.push_back("Ra");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(226); //NIST, 2015
  //Actinium, 89
  Typs.push_back("Ac");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(227); //NIST, 2015
  //Thorium, 90
  Typs.push_back("Th");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(232.0377); //NIST, 2015
  //Protactium, 91
  Typs.push_back("Pa");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(231.03588); //NIST, 2015
  //Uranium, 92
  Typs.push_back("U");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(238.02891); //NIST, 2015
  //Neptunium, 93
  Typs.push_back("Np");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(237); //NIST, 2015
  //Plutonium, 94
  Typs.push_back("Pu");
  CovRadii.push_back(1.0); //Default value
  vdWRadii.push_back(1.0); //Default value
  AtMasses.push_back(244); //NIST, 2015
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

