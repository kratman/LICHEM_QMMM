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

//Tangent functions
VectorXd SymmTangent(vector<QMMMAtom>& OldStruct,
         QMMMSettings& QMMMOpts, int Bead)
{
  //Calculate NEB tangents with the transition state in the center
  
  //Initialize tangent and structures
  VectorXd QMTangent(3*(Nqm+Npseudo));
  
  return QMTangent;
};

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
void LICHEMNEB(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Cartesian NEB optimizer which runs sequentially
  
  return;
};

void EnsembleNEB(vector<QMMMAtom>& Struct, fstream& traj,
     QMMMSettings& QMMMOpts)
{
  //Cartesian Ensemble SD/NEB optimizer which runs sequentially
  
  return;
};

