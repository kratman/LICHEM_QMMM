/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for LAMMPS.

*/

//MM utility functions


//MM wrapper functions
double LammpsWrapper(string RunTyp, vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs LAMMPS
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys;
  //Change units
  E *= kcal2eV;
  return E;
};

