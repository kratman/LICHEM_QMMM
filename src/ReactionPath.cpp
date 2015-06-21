/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                              and Alice Walker                               #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Reaction path and transition state search functions for LICHEM.

 References for NEB:
  

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
