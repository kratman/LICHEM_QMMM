/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for NWChem.

 Reference for NWChem:
 

*/

//MM utility functions


//MM wrapper functions
double NWChemForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem force calculations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys; //Dummy return for system calls
  
  //Change units
  E *= Har2eV;
  return E;
};

void NWChemCharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Calculates atomic charges with NWChem
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys; //Dummy return for system calls
  
  return;
};

double NWChemEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Runs NWChem energy calculations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys; //Dummy return for system calls
  
  //Change units
  E *= Har2eV;
  return E;
};

double NWChemOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem optimizations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys; //Dummy return for system calls
  
  //Change units
  E *= Har2eV;
  return E;
};

