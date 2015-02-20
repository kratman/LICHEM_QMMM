/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for AMBER.

 Reference for AMBER:
 Case et al. AMBER 14, 2014

*/

//MM utility functions


//MM wrapper functions
double AMBERForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs AMBER energy calculations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys; //Dummy return for system calls
  
  //Change units
  E *= kcal2eV;
  return E;
};

double AMBEREnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs AMBER energy calculations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys; //Dummy return for system calls
  
  //Change units
  E *= kcal2eV;
  return E;
};

double AMBEROpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs AMBER optimizations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys; //Dummy return for system calls
  
  //Change units
  E *= kcal2eV;
  return E;
};

