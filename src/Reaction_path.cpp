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

 NB: Reaction path routines are being held back until they are tested and
 the work is published.

 References for NEB:
 

 References for global DFP:
 

*/

//Tangent functions
VectorXd CINEBTangent(VectorXd& Distp1, VectorXd& Distm1,
                      QMMMSettings& QMMMOpts, int Bead)
{
  //Calculate climbing image nudged elastic band tangents
  
  //Initialize tangent and structures
  VectorXd QMTangent(3*(Nqm+Npseudo));
  
  return QMTangent;
};

VectorXd NEBTangent(VectorXd& Distp1, VectorXd& Distm1,
                    QMMMSettings& QMMMOpts, int Bead)
{
  //Calculate nudged elastic band tangents
  
  //Initialize tangent and structures
  VectorXd QMTangent(3*(Nqm+Npseudo));
  
  return QMTangent;
};

//Convergence test functions
bool PathConverged(vector<QMMMAtom>& Struct, vector<QMMMAtom>& OldStruct,
                   MatrixXd& ForceStats, int stepct, QMMMSettings& QMMMOpts,
                   bool QMregion)
{
  //Check convergence of QMMM optimizations
  bool PathDone = 1;
  
  return PathDone;
};

//Path optimization routines
void LICHEMNEB(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int optct)
{
  //Cartesian NEB optimizer which runs sequentially
  
  return;
};

//Path ensemble samping routines
bool FBNEBMCMove(vector<QMMMAtom>& Struct, vector<VectorXd>& AllForces,
                 QMMMSettings& QMMMOpts, VectorXd& Emc)
{
  //Function to try a force-bias Monte Carlo move

  return 0;
};

void EnsembleNEB(vector<QMMMAtom>& Struct, fstream& traj,
                 QMMMSettings& QMMMOpts)
{
  //Cartesian Ensemble SD/NEB optimizer which runs sequentially
  
  return;
};

