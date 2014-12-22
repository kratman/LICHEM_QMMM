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
double LAMMPSForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double Emm = 0.0;
  int ct;
  int sys;
  call.str("");
  //Construct LAMMPS input

  //Return
  return Emm;
};

//MM wrapper functions
double LAMMPSWrapper(string RunTyp, vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs LAMMPS
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys;
  //Construct LAMMPS input

  //Change units
  E *= kcal2eV;
  return E;
};

