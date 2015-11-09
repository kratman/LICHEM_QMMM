/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Physical constants defined in LICHEM.

*/

//Make including safe
#ifndef LICHEM_CONST
#define LICHEM_CONST

//Namespace for constants
namespace LICHEMConst
{
  //Global exact constants
  const double pi = 4*atan(1); //Pi
  const double sqrt2 = pow(2,0.5); //Square root of 2
  const double HugeNum = 1e50; //Large number to reject step
  const double fs2s = 1e-12; //Convert fs to s
  const double m2Ang = 1.0e10; //Angstroms to meters
  const double atm2Pa = 1.01325e5; //Atmospheres to Pascal

  //Global measured constants (NIST, CODATA 2010)
  const double EpsZero = 8.54187817e-12; //Electric constant
  const double hbar = 6.58211928e-16; //Reduced Planck Constant (eV)
  const double k = 8.6173324e-5; //Boltzmann constant (eV)
  const double kSI = 1.3806488e-23; //Boltzmann constant (SI)
  const double amu2kg = 1.660538921e-27; //Atomic mass units to kg
  const double SI2eV = 1/(1.602176565e-19); //Convert SI to eV
  const double Masse = 9.10938291e-31; //Mass of an electron (kg)
  const double BohrRad = 0.52917721092; //Bohr radius (Ang)
  const double Har2eV = 27.21138505; //Hartrees to eV
  const double Na = 6.02214129e23; //Avogadro's number
  const double Debye2au = 0.393430307; //Convert from Debye to au

  //Global derived constants
  const double atm2eV = SI2eV*atm2Pa/(m2Ang*m2Ang*m2Ang); //atmA^3 to eV
  const double C2eV = m2Ang/(4*pi*SI2eV*EpsZero); //Coulomb to eV
  const double ElecMass = Masse/amu2kg; //Mass of an electron (amu)
  const double ToeV = amu2kg*SI2eV/(m2Ang*m2Ang); //Convert to eV units (PIMC)
  const double kcal2eV = 4184*SI2eV/Na; //kcal/mol to eV
};

#endif

