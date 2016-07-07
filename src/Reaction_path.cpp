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
void CheckNEBTangent(VectorXd& tangent)
{
  //Check if the tangent is reasonable
  
  return;
};

VectorXd CINEBTangent(VectorXd& distp1, VectorXd& distm1,
                      QMMMSettings& QMMMOpts, int bead)
{
  //Calculate climbing image nudged elastic band tangents
  
  //Initialize tangent and structures
  VectorXd QMTangent(3*(Nqm+Npseudo));
  
  return QMTangent;
};

VectorXd NEBTangent(VectorXd& distp1, VectorXd& distm1,
                    QMMMSettings& QMMMOpts, int bead)
{
  //Calculate nudged elastic band tangents
  
  //Initialize tangent and structures
  VectorXd QMTangent(3*(Nqm+Npseudo));
  
  return QMTangent;
};

//Convergence test functions
bool PathConverged(vector<QMMMAtom>& QMMMData, vector<QMMMAtom>& oldQMMMData,
                   MatrixXd& forceStats, int stepCt, QMMMSettings& QMMMOpts,
                   bool QMRegion)
{
  //Check convergence of QMMM optimizations
  bool pathDone = 1;
  
  return pathDone;
};

//Path optimization routines
void LICHEMNEB(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int optCt)
{
  //Cartesian NEB optimizer which runs sequentially
  
  return;
};

//Path ensemble samping routines
bool FBNEBMCMove(vector<QMMMAtom>& QMMMData, vector<VectorXd>& allForces,
                 QMMMSettings& QMMMOpts, VectorXd& Emc)
{
  //Function to try a force-bias Monte Carlo move
  
  return 0;
};

