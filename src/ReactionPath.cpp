/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                             and Alice Walker                               #
#                                                                            #
##############################################################################

 Reaction path and transition state search functions for FLUKE.

*/

//Convergence test functions
bool PathConverged(vector<QMMMAtom>& Struct, vector<QMMMAtom>& OldStruct,
     vector<vector<double> >& ForceStats, int stepct, QMMMSettings& QMMMOpts,
     bool QMregion)
{
  //Check convergence of QMMM optimizations
  bool PathDone = 0;
  
  return PathDone;
};

//Path optimization routines
void EnsembleNEB(vector<QMMMAtom>& Struct, fstream& traj,
     QMMMSettings& QMMMOpts)
{
  //Cartesian Ensemble SD/NEB optimizer which runs sequentially
  
  return;
};
