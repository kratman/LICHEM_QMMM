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

double Coord::vecMag()
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

GEMDen::GEMDen(string typ, string basName)
{
  //Fancy constructor
  setBasis(typ,basName);
  return;
};

GEMDen::~GEMDen()
{
  //Generic destructor
  return;
};

void GEMDen::setBasis(string typ, string basName)
{
  //Set the basis set
  dens_ = HermBasis(typ,basName);
  return;
};

void GEMDen::setFrame(bool flip, string frame, int at1, int at2, int at3)
{
  //Set the local frame of reference
  chiralFlip_ = flip;
  type_ = frame;
  atom1_ = at1;
  atom2_ = at2;
  atom3_ = at3;
  return;
};

Mpole GEMDen::GEMDM()
{
  //Function to convert GEM density to distributed multipoles
  Mpole dmpole; //Blank set of multipoles
  //Save frame of reference
  dmpole.chiralFlip = chiralFlip_;
  dmpole.type = type_;
  dmpole.atom1 = atom1_;
  dmpole.atom2 = atom2_;
  dmpole.atom3 = atom3_;
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
  dmpole.q += chemTable.revTyping(type_);
  //Convert Hermite Gaussians to multipoles
  for (unsigned int i=0;i<dens_.size();i++)
  {
    //Check for a monopole
    if ((dens_[i].xPow() == 0) and (dens_[i].yPow() == 0) and
       (dens_[i].zPow() == 0))
    {
      //Update monopole
      dmpole.q += dens_[i].coeff();
      //Update diagonal quadrupole moments
      dmpole.Qxx += dens_[i].coeff()/(2*dens_[i].getAlpha());
      dmpole.Qyy += dens_[i].coeff()/(2*dens_[i].getAlpha());
      dmpole.Qzz += dens_[i].coeff()/(2*dens_[i].getAlpha());
    }
    //Check for a dipole
    if ((dens_[i].xPow() == 1) and (dens_[i].yPow() == 0) and
       (dens_[i].zPow() == 0))
    {
      //Update x dipole
      dmpole.Dx += dens_[i].coeff();
    }
    if ((dens_[i].xPow() == 0) and (dens_[i].yPow() == 1) and
       (dens_[i].zPow() == 0))
    {
      //Update y dipole
      dmpole.Dy += dens_[i].coeff();
    }
    if ((dens_[i].xPow() == 0) and (dens_[i].yPow() == 0) and
       (dens_[i].zPow() == 1))
    {
      //Update z dipole
      dmpole.Dz += dens_[i].coeff();
    }
    //Check for a quadrupole
    if ((dens_[i].xPow() == 2) and (dens_[i].yPow() == 0) and
       (dens_[i].zPow() == 0))
    {
      //Update xx quadrupole
      dmpole.Qxx += 2*dens_[i].coeff();
    }
    if ((dens_[i].xPow() == 1) and (dens_[i].yPow() == 1) and
       (dens_[i].zPow() == 0))
    {
      //Update xy quadrupole
      dmpole.Qxy += dens_[i].coeff();
    }
    if ((dens_[i].xPow() == 1) and (dens_[i].yPow() == 0) and
       (dens_[i].zPow() == 1))
    {
      //Update xz quadrupole
      dmpole.Qxz += dens_[i].coeff();
    }
    if ((dens_[i].xPow() == 0) and (dens_[i].yPow() == 2) and
       (dens_[i].zPow() == 0))
    {
      //Update xx quadrupole
      dmpole.Qyy += 2*dens_[i].coeff();
    }
    if ((dens_[i].xPow() == 0) and (dens_[i].yPow() == 1) and
       (dens_[i].zPow() == 1))
    {
      //Update yz quadrupole
      dmpole.Qyz += dens_[i].coeff();
    }
    if ((dens_[i].xPow() == 0) and (dens_[i].yPow() == 0) and
       (dens_[i].zPow() == 2))
    {
      //Update yz quadrupole
      dmpole.Qzz += 2*dens_[i].coeff();
    }
  }
  //Convert to a traceless quadrupole
  double qTrace = dmpole.Qxx+dmpole.Qyy+dmpole.Qzz;
  dmpole.Qxx -= qTrace;
  dmpole.Qyy -= qTrace;
  dmpole.Qzz -= qTrace;
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
  func = "N/A";
  basis = "N/A";
  RAM = 256;
  memMB = 1;
  charge = 0;
  spin = 1;
  backDir = "Old_files";
  //QMMM long-range electrostatics settings
  useLREC = 0;
  LRECCut = 1000.0; //Effectively infinite
  LRECPow = 3;
  //MM wrapper settings
  useMMCut = 0;
  MMOptCut = 1000.0; //Effectively infinite
  useEwald = 0;
  useImpSolv = 0;
  solvModel = "N/A";
  //MC, MD, and RP settings
  ensemble = "N/A";
  temp = 300.0;
  beta = 1/(300.0*kBoltz);
  press = 0.0;
  NEq = 0;
  NSteps = 0;
  NBeads = 1; //Key for printing
  accRatio = 0.5;
  NPrint = 5000;
  dt = 1.0;
  tauTemp = 1000.0;
  //Optimization settings
  maxOptSteps = 200;
  MMOptTol = 1e-2;
  QMOptTol = 5e-4;
  stepScale = 1.0;
  maxStep = 0.1;
  //Additional RP settings
  kSpring = 1.0;
  TSBead = 0;
  climb = 0;
  frznEnds = 0;
  NEBFreq = 0;
  printNormModes = 0;
  startPathChk = 1; //Speeds up reaction pathways
  //Temporary energy storage
  EOld = 0.0;
  EReact = 0.0;
  EProd = 0.0;
  ETrans = 0.0;
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
  typs_.push_back("H");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.00784); //NIST, 2015
  //Helium, 2
  typs_.push_back("He");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(4.002602); //NIST, 2015
  //Lithium, 3
  typs_.push_back("Li");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(6.938); //NIST, 2015
  //Berylium, 4
  typs_.push_back("Be");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(9.0121831); //NIST, 2015
  //Boron, 5
  typs_.push_back("B");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(10.806); //NIST, 2015
  //Carbon, 6
  typs_.push_back("C");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(12.0096); //NIST, 2015
  //Nitrogen, 7
  typs_.push_back("N");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(14.00643); //NIST, 2015
  //Oxygen, 8
  typs_.push_back("O");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(15.99903); //NIST, 2015
  //Fluorine, 9
  typs_.push_back("F");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(18.998403163); //NIST, 2015
  //Neon, 10
  typs_.push_back("Ne");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(20.1797); //NIST, 2015
  //Sodium, 11
  typs_.push_back("Na");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(22.98976928); //NIST, 2015
  //Magnesium, 12
  typs_.push_back("Mg");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(24.304); //NIST, 2015
  //Aluminum, 13
  typs_.push_back("Al");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(26.9815385); //NIST, 2015
  //Silicon, 14
  typs_.push_back("Si");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(28.084); //NIST, 2015
  //Phosphorus, 15
  typs_.push_back("P");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(30.973761998); //NIST, 2015
  //Sulfur, 16
  typs_.push_back("S");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(32.059); //NIST, 2015
  //Chlorine, 17
  typs_.push_back("Cl");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(35.446); //NIST, 2015
  //Argon, 18
  typs_.push_back("Ar");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(39.948); //NIST, 2015
  //Potassium, 19
  typs_.push_back("K");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(39.0983); //NIST, 2015
  //Calcium, 20
  typs_.push_back("Ca");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(40.078); //NIST, 2015
  //Scandium, 21
  typs_.push_back("Sc");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(44.955908); //NIST, 2015
  //Titanium, 22
  typs_.push_back("Ti");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(47.867); //NIST, 2015
  //Vanadium, 23
  typs_.push_back("V");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(50.9415); //NIST, 2015
  //Chromium, 24
  typs_.push_back("Cr");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(51.9961); //NIST, 2015
  //Manganese, 25
  typs_.push_back("Mn");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(54.938044); //NIST, 2015
  //Iron, 26
  typs_.push_back("Fe");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(55.845); //NIST, 2015
  //Cobalt, 27
  typs_.push_back("Co");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(58.933194); //NIST, 2015
  //Nickel, 28
  typs_.push_back("Ni");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(58.6934); //NIST, 2015
  //Copper, 29
  typs_.push_back("Cu");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(63.546); //NIST, 2015
  //Zinc, 30
  typs_.push_back("Zn");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(65.38); //NIST, 2015
  //Gallium, 31
  typs_.push_back("Ga");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(69.723); //NIST, 2015
  //Germanium, 32
  typs_.push_back("Ge");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(72.630); //NIST, 2015
  //Arsenic, 33
  typs_.push_back("As");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(74.921595); //NIST, 2015
  //Selenium, 34
  typs_.push_back("Se");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(78.971); //NIST, 2015
  //Bromine, 35
  typs_.push_back("Br");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(79.901); //NIST, 2015
  //Krypton, 36
  typs_.push_back("Kr");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(83.798); //NIST, 2015
  //Rubidium, 37
  typs_.push_back("Rb");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(85.4678); //NIST, 2015
  //Strontium, 38
  typs_.push_back("Sr");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(87.62); //NIST, 2015
  //Yttrium, 39
  typs_.push_back("Y");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(88.90584); //NIST, 2015
  //Zirconium, 40
  typs_.push_back("Zr");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(91.224); //NIST, 2015
  //Niobium, 41
  typs_.push_back("Nb");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(92.90637); //NIST, 2015
  //Molybdenum, 42
  typs_.push_back("Mo");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(95.95); //NIST, 2015
  //Technetium, 43
  typs_.push_back("Tc");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(98); //NIST, 2015
  //Ruthenium, 44
  typs_.push_back("Ru");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(101.07); //NIST, 2015
  //Rhodium, 45
  typs_.push_back("Rh");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(102.90550); //NIST, 2015
  //Palladium, 46
  typs_.push_back("Pd");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(106.42); //NIST, 2015
  //Silver, 47
  typs_.push_back("Ag");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(107.8682); //NIST, 2015
  //Cadmium, 48
  typs_.push_back("Cd");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(112.414); //NIST, 2015
  //Indium, 49
  typs_.push_back("In");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(114.818); //NIST, 2015
  //Tin, 50
  typs_.push_back("Sn");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(118.710); //NIST, 2015
  //Antimony, 51
  typs_.push_back("Sb");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(121.760); //NIST, 2015
  //Tellurium, 52
  typs_.push_back("Te");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(127.60); //NIST, 2015
  //Iodine, 53
  typs_.push_back("I");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(126.90447); //NIST, 2015
  //Xenon, 54
  typs_.push_back("Xe");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(131.293); //NIST, 2015
  //Caesium, 55
  typs_.push_back("Cs");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(132.90545196); //NIST, 2015
  //Barium, 56
  typs_.push_back("Ba");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(137.327); //NIST, 2015
  //Lanthanum, 57
  typs_.push_back("La");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(138.90547); //NIST, 2015
  //Cerium, 58
  typs_.push_back("Ce");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(140.116); //NIST, 2015
  //Praseodymium, 59
  typs_.push_back("Pr");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(140.90766); //NIST, 2015
  //Neodymium, 60
  typs_.push_back("Nd");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(144.242); //NIST, 2015
  //Promethium, 61
  typs_.push_back("Pm");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(145); //NIST, 2015
  //Samarium, 62
  typs_.push_back("Sm");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(150.36); //NIST, 2015
  //Europium, 63
  typs_.push_back("Eu");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(151.964); //NIST, 2015
  //Gadolinium, 64
  typs_.push_back("Gd");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(157.25); //NIST, 2015
  //Terbium, 65
  typs_.push_back("Tb");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(158.92535); //NIST, 2015
  //Dysprosium, 66
  typs_.push_back("Dy");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(162.500); //NIST, 2015
  //Holmium, 67
  typs_.push_back("Ho");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(164.93033); //NIST, 2015
  //Erbium, 68
  typs_.push_back("Er");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(167.259); //NIST, 2015
  //Thulium, 69
  typs_.push_back("Tm");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(168.93422); //NIST, 2015
  //Ytterbium, 70
  typs_.push_back("Yb");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(173.054); //NIST, 2015
  //Lutetium, 71
  typs_.push_back("Lu");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(174.9668); //NIST, 2015
  //Hafnium, 72
  typs_.push_back("Hf");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(178.49); //NIST, 2015
  //Tantalum, 73
  typs_.push_back("Ta");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(180.94788); //NIST, 2015
  //Tungsten, 74
  typs_.push_back("W");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(183.84); //NIST, 2015
  //Rhenium, 75
  typs_.push_back("Re");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(186.207); //NIST, 2015
  //Osmium, 76
  typs_.push_back("Os");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(190.23); //NIST, 2015
  //Iridium, 77
  typs_.push_back("Ir");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(192.217); //NIST, 2015
  //Platinum, 78
  typs_.push_back("Pt");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(195.084); //NIST, 2015
  //Gold, 79
  typs_.push_back("Au");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(196.966569); //NIST, 2015
  //Mercury, 80
  typs_.push_back("Hg");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(200.592); //NIST, 2015
  //Thallium, 81
  typs_.push_back("Tl");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(204.382); //NIST, 2015
  //Lead, 82
  typs_.push_back("Pb");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(207.2); //NIST, 2015
  //Bismuth, 83
  typs_.push_back("Bi");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(208.98040); //NIST, 2015
  //Polonium, 84
  typs_.push_back("Po");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(209); //NIST, 2015
  //Astatine, 85
  typs_.push_back("At");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(210); //NIST, 2015
  //Radon, 86
  typs_.push_back("Rn");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(222); //NIST, 2015
  //Francium, 87
  typs_.push_back("Fr");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(223); //NIST, 2015
  //Radium, 88
  typs_.push_back("Ra");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(226); //NIST, 2015
  //Actinium, 89
  typs_.push_back("Ac");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(227); //NIST, 2015
  //Thorium, 90
  typs_.push_back("Th");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(232.0377); //NIST, 2015
  //Protactium, 91
  typs_.push_back("Pa");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(231.03588); //NIST, 2015
  //Uranium, 92
  typs_.push_back("U");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(238.02891); //NIST, 2015
  //Neptunium, 93
  typs_.push_back("Np");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(237); //NIST, 2015
  //Plutonium, 94
  typs_.push_back("Pu");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(244); //NIST, 2015
  //Americium, 95
  typs_.push_back("Am");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Curium, 96
  typs_.push_back("Cm");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Berkelium, 97
  typs_.push_back("Bk");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Californium, 98
  typs_.push_back("Cf");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Einsteinium, 99
  typs_.push_back("Es");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Fermium, 100
  typs_.push_back("Fm");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Mendelevium, 101
  typs_.push_back("Md");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Nobelium, 102
  typs_.push_back("No");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Lawrencium, 103
  typs_.push_back("Lr");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Rutherfordium, 104
  typs_.push_back("Rf");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Dubnium, 105
  typs_.push_back("Db");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Seaborgium, 106
  typs_.push_back("Sg");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Bohrium, 107
  typs_.push_back("Bh");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Hasium, 108
  typs_.push_back("Hs");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Meitnerium, 109
  typs_.push_back("Mt");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Darmstadtium, 110
  typs_.push_back("Ds");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Roentgenium, 111
  typs_.push_back("Rg");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Copernicium, 112
  typs_.push_back("Cn");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Num. 113
  typs_.push_back("");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Flerovium, 114
  typs_.push_back("Fl");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Num. 115
  typs_.push_back("");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Livermorium, 116
  typs_.push_back("Lv");
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Num. 117
  typs_.push_back(""); //Unnamed
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  //Num. 118
  typs_.push_back(""); //Unnamed
  covRadii_.push_back(1.0); //Default value
  vdWRadii_.push_back(1.0); //Default value
  atMasses_.push_back(1.0); //Default value
  return;
};

PeriodicTable::~PeriodicTable()
{
  //Generic destructor
  return;
};

string PeriodicTable::typing(int Z)
{
  //Function to convert nuclear charges to atom types
  return typs_[Z-1];
};

int PeriodicTable::revTyping(string atName)
{
  //Function to convert atom types to nuclear charges
  int Z = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<typs_.size();i++)
  {
    if (atName == typs_[i])
    {
      Z = i+1;
    }
  }
  return Z;
};

double PeriodicTable::getCovRadius(string atName)
{
  //Function to find the covalent radius of an atom
  double radius = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<typs_.size();i++)
  {
    if (atName == typs_[i])
    {
      //Save value
      radius = covRadii_[i];
    }
  }
  return radius;
};

double PeriodicTable::getRadius(string atName)
{
  //Function to find the vdW radius of an atom
  double radius = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<typs_.size();i++)
  {
    if (atName == typs_[i])
    {
      //Save value
      radius = vdWRadii_[i];
    }
  }
  return radius;
};

double PeriodicTable::getAtMass(string atName)
{
  //Function to find the atomic mass of an atom
  double mass = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (unsigned int i=0;i<typs_.size();i++)
  {
    if (atName == typs_[i])
    {
      //Save value
      mass = atMasses_[i];
    }
  }
  return mass;
};

