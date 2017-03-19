/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                               Alice R. Walker                               #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Reaction path and transition state search functions for LICHEM.

 References for NEB:
 Henkelman et al., J. Chem. Phys., 113, 22, 9901, (2000)
 Henkelman et al., J. Chem. Phys., 113, 22, 9978, (2000)

 References for global DFP:
 Press et al., Numerical Recipes 3nd Edition, (2007)
 Sheppard et al., J. Chem. Phys., 128, 13, 134106, (2008)

 NB: The NEB and DFP algorithms used in LICHEM are based on the references
 above, but were modified to better handle QMMM systems.

*/

/*!
  \ingroup Optimizers
*/
///@{

//Tangent functions
void CheckNEBTangent(VectorXd& tangent)
{
  //Check if the tangent is reasonable
  double tanNorm = tangent.squaredNorm();
  if (tanNorm < (1e-12))
  {
    //Zero
    tangent.setZero();
  }
  else if (tanNorm != tanNorm)
  {
    //NaN
    tangent.setZero();
  }
  else if ((tanNorm >= hugeNum) || (tanNorm <= (-hugeNum)))
  {
    //Inf
    tangent.setZero();
  }
  return;
};

VectorXd CINEBTangent(VectorXd& distp1, VectorXd& distm1,
                      QMMMSettings& QMMMOpts, int bead)
{
  //Calculate climbing image nudged elastic band tangents
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Calculate tangent
  VectorXd QMTangent;
  if ((bead < QMMMOpts.TSBead) && (bead != 0))
  {
    //Reactant side of the TS
    QMTangent = distp1;
  }
  if ((bead > QMMMOpts.TSBead) && (bead != (QMMMOpts.NBeads-1)))
  {
    //Product side of the TS
    QMTangent = distm1;
  }
  if ((bead == QMMMOpts.TSBead) && (bead != 0) &&
     (bead != (QMMMOpts.NBeads-1)))
  {
    //Transition state
    VectorXd QMTangent2(Ndof); //Second tagent vector
    //First tangent
    QMTangent = distp1;
    //Second tangent
    QMTangent2 = distm1;
    //Normalize and add vectors
    QMTangent.normalize();
    QMTangent2.normalize();
    QMTangent += QMTangent2;
  }
  //Normalize tangent
  QMTangent.normalize();
  CheckNEBTangent(QMTangent);
  //Return final tangent
  return QMTangent;
};

VectorXd NEBTangent(VectorXd& distp1, VectorXd& distm1,
                    QMMMSettings& QMMMOpts, int bead)
{
  //Calculate nudged elastic band tangents
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize tangent and structures
  VectorXd QMTangent(Ndof);
  VectorXd QMTangent2(Ndof);
  //Calculate tangent
  QMTangent = distp1; //Copy p+1 direction
  QMTangent2 = distm1; //Copy p-1 direction
  QMTangent.normalize(); //Normalize p+1 direction
  QMTangent2.normalize(); //Normalize p-1 direction
  QMTangent += QMTangent2; //Combine tangents
  QMTangent.normalize(); //Normalize combined tangent
  CheckNEBTangent(QMTangent);
  //Return final tangent
  return QMTangent;
};

//Convergence test functions
bool PathConverged(vector<QMMMAtom>& QMMMData, vector<QMMMAtom>& oldQMMMData,
                   MatrixXd& forceStats, int stepCt, QMMMSettings& QMMMOpts,
                   bool QMRegion)
{
  //Check convergence of QMMM optimizations
  //Initialize stats variables
  bool pathDone = 0;
  double RMSDiff = 0;
  double RMSForce = 0;
  double maxForce = 0;
  double sumE = 0;
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Convergence criteria
  double maxFTol = 20*QMMMOpts.QMOptTol; //Opt. tolerance for max. force
  double RMSFTol = 10*QMMMOpts.QMOptTol; //Opt. tolerance for RMS force
  //Check progress of all beads
  if (QMRegion)
  {
    //Check if a QM calculation is converged
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Find max forces and RMSforce
      if (maxForce < forceStats(p,0))
      {
        maxForce = forceStats(p,0);
      }
      RMSForce += forceStats(p,1);
      //Find RMS deviation for the whole path
      #pragma omp parallel for schedule(dynamic) reduction(+:RMSDiff)
      for (int i=0;i<Natoms;i++)
      {
        //Calculate RMS displacement
        double RMSTemp = 0; //Store a local sum
        if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
        {
          for (int j=0;j<i;j++)
          {
            if (QMMMData[j].QMRegion || QMMMData[j].PBRegion)
            {
              double RNew = 0;
              double ROld = 0;
              RNew = CoordDist2(QMMMData[i].P[p],
                                QMMMData[j].P[p]).vecMag();
              ROld = CoordDist2(oldQMMMData[i].P[p],
                                oldQMMMData[j].P[p]).vecMag();
              RNew = sqrt(RNew);
              ROld = sqrt(ROld);
              //Update local sum
              RMSTemp += (RNew-ROld)*(RNew-ROld);
            }
          }
        }
        //Update sum
        RMSDiff += RMSTemp;
      }
    }
    int adjustedBeads; //Number of moving beads
    if (QMMMOpts.frznEnds)
    {
      //End points do not count
      adjustedBeads = QMMMOpts.NBeads-2;
    }
    else
    {
      //All beads
      adjustedBeads = QMMMOpts.NBeads;
    }
    RMSDiff /= (Nqm+Npseudo)*(Nqm+Npseudo-1)/2;
    RMSDiff /= adjustedBeads; //Adjust for multiple replicas
    RMSDiff = sqrt(RMSDiff);
    RMSForce /= Ndof;
    RMSForce /= adjustedBeads; //Adjust for multiple replicas
    RMSForce = sqrt(RMSForce);
    //Print progress
    cout << "    QM step: " << stepCt;
    cout << " | RMS dev: ";
    cout << LICHEMFormFloat(RMSDiff,12);
    cout << " \u212B" << '\n';
    cout << "    Max. force: ";
    cout << LICHEMFormFloat(maxForce,12);
    cout << " eV/\u212B | RMS force: ";
    cout << LICHEMFormFloat(RMSForce,12);
    cout << " eV/\u212B" << '\n';
    cout << "    React. energy: ";
    cout << LICHEMFormFloat(QMMMOpts.EReact,14);
    cout << " eV | Prod. energy: ";
    cout << LICHEMFormFloat(QMMMOpts.EProd,14);
    cout << " eV" << '\n';
    cout << "    TS bead: " << QMMMOpts.TSBead;
    cout << " | TS energy: ";
    cout << LICHEMFormFloat(QMMMOpts.ETrans,14);
    cout << " eV" << '\n';
    //Check convergence criteria
    if ((RMSDiff <= QMMMOpts.QMOptTol) && (RMSForce <= RMSFTol) &&
       (maxForce <= maxFTol))
    {
      //Make sure the TS is converged too...
      if (NEBSim && (!QMMMOpts.climb))
      {
        //Turn on climbing image forces
        QMMMOpts.climb = 1;
        cout << "    QM is nearly converged. Starting climbing image NEB...";
        cout << '\n';
      }
      else
      {
        //Finish the optimization
        pathDone = 1;
        cout << "    QM optimization complete." << '\n';
      }
    }
    cout << '\n';
    cout.flush();
  }
  if (!QMRegion)
  {
    //Check if the MM region changed and gather statistics
    int pathStart = 0;
    int pathEnd = QMMMOpts.NBeads;
    if (QMMMOpts.frznEnds)
    {
      //Skip energy calculations on the end points
      pathStart = 1;
      pathEnd = QMMMOpts.NBeads-1;
    }
    VectorXd Es(QMMMOpts.NBeads);
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      sumE = 0;
      //Calculate QM energy
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        sumE += GaussianEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        sumE += PSI4Energy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        sumE += NWChemEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        sumE += TINKEREnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      //Calculate RMS displacement
      #pragma omp parallel for schedule(dynamic) reduction(+:RMSDiff)
      for (int i=0;i<Natoms;i++)
      {
        double RMSTemp = 0; //Store a local sum
        for (int j=0;j<i;j++)
        {
          double RNew = 0;
          double ROld = 0;
          RNew = CoordDist2(QMMMData[i].P[p],
                            QMMMData[j].P[p]).vecMag();
          ROld = CoordDist2(oldQMMMData[i].P[p],
                            oldQMMMData[j].P[p]).vecMag();
          RNew = sqrt(RNew);
          ROld = sqrt(ROld);
          //Update local sum
          RMSTemp += (RNew-ROld)*(RNew-ROld);
        }
        //Update sum
        RMSDiff += RMSTemp;
      }
      Es(p) = sumE;
    }
    int adjustedBeads; //Number of moving beads
    if (QMMMOpts.frznEnds)
    {
      //End points do not count
      adjustedBeads = QMMMOpts.NBeads-2;
    }
    else
    {
      //All beads
      adjustedBeads = QMMMOpts.NBeads;
    }
    RMSDiff /= (Natoms-Nfreeze)*(Natoms-Nfreeze-1)/2;
    RMSDiff /= adjustedBeads; //Adjust for multiple replicas
    RMSDiff = sqrt(RMSDiff);
    //Update energies
    sumE = 0; //Reusing this variable to avoid making a new one
    sumE = Es.maxCoeff();
    for (int p=pathStart;p<pathEnd;p++)
    {
      if (p == 0)
      {
        //Update energy if end points are active
        QMMMOpts.EReact = Es(p);
      }
      if (Es(p) == sumE)
      {
        //Save new energy
        QMMMOpts.TSBead = p;
        QMMMOpts.ETrans = Es(p);
      }
      if (p == (QMMMOpts.NBeads-1))
      {
        //Update energy if end points are active
        QMMMOpts.EProd = Es(p);
      }
    }
    //Calculate reaction coordinate
    VectorXd reactCoord(QMMMOpts.NBeads); //Reaction coordinate
    reactCoord.setZero();
    for (int p=0;p<(QMMMOpts.NBeads-1);p++)
    {
      MatrixXd geom1((Nqm+Npseudo),3); //Current replica
      MatrixXd geom2((Nqm+Npseudo),3); //Next replica
      VectorXd disp; //Store the displacement
      //Save geometries
      int ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
        {
          //Save current replica
          geom1(ct,0) = QMMMData[i].P[p].x;
          geom1(ct,1) = QMMMData[i].P[p].y;
          geom1(ct,2) = QMMMData[i].P[p].z;
          //Save replica p+1
          geom2(ct,0) = QMMMData[i].P[p+1].x;
          geom2(ct,1) = QMMMData[i].P[p+1].y;
          geom2(ct,2) = QMMMData[i].P[p+1].z;
          ct += 1;
        }
      }
      //Calculate displacement
      disp = KabschDisplacement(geom1,geom2,(Nqm+Npseudo));
      if (NEBSim)
      {
        //Remove inactive atoms
        ct = 0; //Reset counter for the number of atoms
        for (int i=0;i<Natoms;i++)
        {
          //Only include QM and PB regions
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            //Only include active atoms in the tangent
            if (!QMMMData[i].NEBActive)
            {
              //Delete distance components
              disp(ct) = 0;
              disp(ct+1) = 0;
              disp(ct+2) = 0;
            }
            //Advance counter
            ct += 3;
          }
        }
      }
      //Update reaction coordinate
      reactCoord(p+1) = reactCoord(p); //Start from previous bead
      reactCoord(p+1) += disp.norm(); //Add magnitude of the displacement
    }
    reactCoord /= reactCoord.maxCoeff(); //Must be between 0 and 1
    //Print progress
    cout << " | Opt. step: " << stepCt;
    cout << " | TS bead: " << QMMMOpts.TSBead;
    cout << " | TS energy: ";
    cout << LICHEMFormFloat(Es(QMMMOpts.TSBead),16);
    cout << " eV" << '\n';
    cout << " | RMS dev: " << LICHEMFormFloat(RMSDiff,12);
    cout << " \u212B | Bead energies:" << '\n';
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      cout << "   Bead: " << setw(3) << p;
      cout << " | React coord: " << LICHEMFormFloat(reactCoord(p),5);
      cout << " | Energy: " << LICHEMFormFloat(Es(p),16);
      cout << " eV" << '\n';
    }
    //Check convergence
    if (RMSDiff <= QMMMOpts.MMOptTol)
    {
      pathDone = 1;
      if (QMMM && (stepCt > 1))
      {
        cout << "    QMMM relaxation satisfactory.";
        cout << '\n';
      }
    }
    //Flush output
    cout.flush();
  }
  return pathDone;
};

//Path optimization routines
void LICHEMNEB(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int optCt)
{
  //Cartesian NEB-DFP optimizer which globally optimizes the reaction path
  fstream qmFile; //Generic file stream
  stringstream call; //Stream for system calls and reading/writing files
  int stepCt = 0; //Counter for optimization steps
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  double maxFTol = 20*QMMMOpts.QMOptTol; //Opt. tolerance for max. force
  double RMSFTol = 10*QMMMOpts.QMOptTol; //Opt. tolerance for RMS force
  //Initialize charges
  if (Nmm > 0)
  {
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Write charges for bead p
      WriteChargeFile(QMMMData,QMMMOpts,p);
    }
  }
  //Initialize trajectory file
  call.str("");
  call << "LICHMNEBOpt.xyz";
  qmFile.open(call.str().c_str(),ios_base::out);
  //Set end points for the optimization
  int pathStart = 0;
  int pathEnd = QMMMOpts.NBeads;
  if (QMMMOpts.frznEnds)
  {
    //Change the start and end points
    pathStart = 1;
    pathEnd = QMMMOpts.NBeads-1;
  }
  //Initialize optimization variables
  int newTS; //Storage for new TS ID
  double newTSEnergy; //Saved TS energy for updating
  double stepScale; //Local copy
  double EOld = 0; //Previous total energy
  double sumE = 0; //Current total energy
  double sdScale = 0.01; //Scale factor for SD steps
  double vecMax = 0; //Maximum displacement
  bool pathDone = 0; //Flag to end the optimization
  double forcesDotTan; //Overlap of vectors
  double springDist; //Distance between the neighboring beads
  //Create DFP arrays
  int matrixSize; //Number of elements in the matrices and arrays
  matrixSize = Ndof; //Only include QM and PB atoms
  matrixSize *= QMMMOpts.NBeads-2; //All replicas execpt the end points
  VectorXd optVecG(matrixSize); //Gradient descent direction
  VectorXd optVecR(Ndof); //Gradient descent direction (reactant)
  VectorXd optVecP(Ndof); //Gradient descent direction (product)
  VectorXd gradDiffG(matrixSize); //Change in the gradient (path)
  VectorXd gradDiffR(Ndof); //Change in the gradient (reactant)
  VectorXd gradDiffP(Ndof); //Change in the gradient (product)
  VectorXd gForces(matrixSize); //Global forces
  VectorXd rForces(Ndof); //Reactant forces
  VectorXd pForces(Ndof); //Product forces
  MatrixXd iHessG(matrixSize,matrixSize); //Path inverse Hessian
  MatrixXd iHessR(Ndof,Ndof); //Reactant inverse Hessian
  MatrixXd iHessP(Ndof,Ndof); //Product inverse Hessian
  VectorXd forces(Ndof); //Local forces
  //Initialize arrays
  optVecG.setZero();
  optVecR.setZero();
  optVecP.setZero();
  gradDiffG.setZero();
  gradDiffR.setZero();
  gradDiffP.setZero();
  gForces.setZero();
  rForces.setZero();
  pForces.setZero();
  forces.setZero();
  //Create an identity matrix as the initial Hessian
  iHessG.setIdentity(); //Already an "inverse" Hessian
  iHessR.setIdentity(); //Already an "inverse" Hessian
  iHessP.setIdentity(); //Already an "inverse" Hessian
  //Create array to store stats and check convergence
  MatrixXd forceStats(QMMMOpts.NBeads,2);
  forceStats.setZero();
  //Run optimization
  newTS = 0; //Reactant
  newTSEnergy = -hugeNum; //All energies will be higher
  for (int p=pathStart;p<pathEnd;p++)
  {
    double E = 0;
    double Emm = 0;
    double Eqm = 0;
    //Create blank force array
    forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      Eqm += GaussianForces(QMMMData,forces,QMMMOpts,p);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      Eqm += PSI4Forces(QMMMData,forces,QMMMOpts,p);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      Eqm += NWChemForces(QMMMData,forces,QMMMOpts,p);
      QMTime += (unsigned)time(0)-tStart;
    }
    E += Eqm;
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,forces,QMMMOpts,p);
      if (AMOEBA || QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,forces,QMMMOpts,p);
      }
      Emm += TINKEREnergy(QMMMData,QMMMOpts,p);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      E += LAMMPSForces(QMMMData,forces,QMMMOpts,p);
      Emm += LAMMPSEnergy(QMMMData,QMMMOpts,p);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Update energies
    sumE += E;
    if (p == 0)
    {
      //Update reactant energy
      QMMMOpts.EReact = Eqm+Emm;
    }
    if (p == QMMMOpts.TSBead)
    {
      //Update old TS energy
      QMMMOpts.ETrans = Eqm+Emm;
    }
    if (p == (QMMMOpts.NBeads-1))
    {
      //Update product energy
      QMMMOpts.EProd = Eqm+Emm;
    }
    if ((Eqm+Emm) > newTSEnergy)
    {
      //Find the current TS
      newTSEnergy = Eqm+Emm; //New energy
      newTS = p; //New TS
    }
    //Modify forces along the tangent
    VectorXd QMTangent(Ndof); //Tangent vector
    VectorXd distp1(Ndof); //Displacement for p+1
    VectorXd distm1(Ndof); //Displacement for p-1
    if ((p != 0) && (p != (QMMMOpts.NBeads-1)))
    {
      //Calculate tangents for middle replicas
      MatrixXd geom1((Nqm+Npseudo),3); //Current replica
      MatrixXd geom2((Nqm+Npseudo),3); //Second replica
      MatrixXd geom3((Nqm+Npseudo),3); //Third replica
      //Save geometries
      int ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
        {
          //Save current replica
          geom1(ct,0) = QMMMData[i].P[p].x;
          geom1(ct,1) = QMMMData[i].P[p].y;
          geom1(ct,2) = QMMMData[i].P[p].z;
          //Save replica p+1
          geom2(ct,0) = QMMMData[i].P[p+1].x;
          geom2(ct,1) = QMMMData[i].P[p+1].y;
          geom2(ct,2) = QMMMData[i].P[p+1].z;
          //Save replica p-1
          geom3(ct,0) = QMMMData[i].P[p-1].x;
          geom3(ct,1) = QMMMData[i].P[p-1].y;
          geom3(ct,2) = QMMMData[i].P[p-1].z;
          //Advance counter
          ct += 1;
        }
      }
      //Calculate displacements
      distp1 = KabschDisplacement(geom1,geom2,(Nqm+Npseudo));
      distp1 *= -1; //Change direction
      distm1 = KabschDisplacement(geom1,geom3,(Nqm+Npseudo));
      //Remove inactive atoms
      ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
        {
          //Only include active atoms in the tangent
          if (!QMMMData[i].NEBActive)
          {
            //Delete distance components
            distp1(ct) = 0;
            distp1(ct+1) = 0;
            distp1(ct+2) = 0;
            distm1(ct) = 0;
            distm1(ct+1) = 0;
            distm1(ct+2) = 0;
          }
          //Advance counter
          ct += 3;
        }
      }
      //Calculate tangent
      QMTangent = CINEBTangent(distp1,distm1,QMMMOpts,p);
    }
    else
    {
      //No tangent for the reactant and product
      QMTangent.setZero();
      //Get rid of displacements for spring force calculation
      distp1.setZero();
      distm1.setZero();
    }
    //Project forces onto tangent vector
    forcesDotTan = forces.dot(QMTangent); //Overlap of vectors
    //Remove forces along the tangent
    forces -= forcesDotTan*QMTangent; //End points have no tangents
    if ((p == QMMMOpts.TSBead) && QMMMOpts.climb)
    {
      //Climbing image for TS
      forces -= forcesDotTan*QMTangent;
      //Get rid of displacements for spring force calculations
      distp1.setZero();
      distm1.setZero();
    }
    //Add p+1 spring forces (distp1.norm()=0 at end points and TS)
    springDist = distp1.norm(); //Distance from p+1
    forces += (QMMMOpts.kSpring*springDist*QMTangent);
    //Add p-1 spring forces (distm1.norm()=0 at end points and TS)
    springDist = distm1.norm(); //Distance from p-1
    forces -= (QMMMOpts.kSpring*springDist*QMTangent);
    //Add forces to global arrays
    if (p == 0)
    {
      //Reactant
      #pragma omp parallel for schedule(dynamic)
      for (int i=0;i<Ndof;i++)
      {
        rForces(i) += forces(i);
      }
    }
    else if (p == (QMMMOpts.NBeads-1))
    {
      //Product
      #pragma omp parallel for schedule(dynamic)
      for (int i=0;i<Ndof;i++)
      {
        pForces(i) += forces(i);
      }
    }
    else
    {
      //Path
      int gfID; //Location of the bead in the global array
      gfID = Ndof; //Adjust for the QM and PB degrees of freedom
      gfID *= p-1; //Adjust for the bead ID
      #pragma omp parallel for schedule(dynamic)
      for (int i=0;i<Ndof;i++)
      {
        gForces(gfID+i) += forces(i);
      }
    }
  }
  //Update TS properties
  if (newTSEnergy >= QMMMOpts.EProd)
  {
    //Favor product as the TS over mid point
    if (newTSEnergy == QMMMOpts.EProd)
    {
      //Make the product the transition state
      newTS = QMMMOpts.NBeads-1;
    }
    //Check if the TS bead has changed
    if (newTS != QMMMOpts.TSBead)
    {
      //Signal that the TS changed
      sumE = -hugeNum;
    }
    //Safely update TS, even if endpoints are frozen
    QMMMOpts.ETrans = newTSEnergy;
    QMMMOpts.TSBead = newTS;
  }
  //Output initial RMS force
  vecMax = 0; //Using this variable to avoid creating a new one
  vecMax += gForces.squaredNorm(); //Add path RMS forces
  vecMax += rForces.squaredNorm(); //Add reactant forces
  vecMax += pForces.squaredNorm(); //Add product forces
  vecMax /= Ndof; //Scale by the number of atoms
  vecMax /= QMMMOpts.NBeads; //Scale by the number of replicas
  vecMax = sqrt(vecMax); //RMS force
  EOld = sumE; //Save the energy
  cout << "    Performing a steepest descent step..." << '\n';
  cout << "    QM step: 0";
  cout << " | RMS force: ";
  cout << LICHEMFormFloat(vecMax,14);
  cout << " eV/\u212B" << '\n';
  cout << "    React. energy: ";
  cout << LICHEMFormFloat(QMMMOpts.EReact,14);
  cout << " eV | Prod. energy: ";
  cout << LICHEMFormFloat(QMMMOpts.EProd,14);
  cout << " eV" << '\n';
  cout << "    TS bead: " << QMMMOpts.TSBead;
  cout << " | TS energy: ";
  cout << LICHEMFormFloat(QMMMOpts.ETrans,14);
  cout << " eV" << '\n' << '\n';
  cout.flush();
  //Optimize path
  stepScale = QMMMOpts.stepScale;
  stepScale *= sdScale; //Take a small first step
  while ((!pathDone) && (stepCt < QMMMOpts.maxOptSteps))
  {
    sumE = 0; //Reinitialize energy
    //Copy old structure and forces
    vector<QMMMAtom> oldQMMMData = QMMMData; //Save structure
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<Ndof;i++)
    {
      //Reinitialize the change in the gradient
      gradDiffR(i) = rForces(i); //Reactant
      gradDiffP(i) = pForces(i); //Product
    }
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<matrixSize;i++)
    {
      //Reinitialize the change in the gradient
      gradDiffG(i) = gForces(i); //Path
    }
    //Determine new structure
    optVecG = iHessG*gForces;
    optVecG *= stepScale;
    optVecR = iHessR*rForces;
    optVecR *= stepScale;
    optVecP = iHessP*pForces;
    optVecP *= stepScale;
    //Check average step size
    vecMax = optVecG.norm()/QMMMOpts.NBeads;
    if (vecMax > QMMMOpts.maxStep)
    {
      //Scale step size
      optVecG *= (QMMMOpts.maxStep/vecMax);
    }
    vecMax = optVecR.norm()/QMMMOpts.NBeads;
    if (vecMax > QMMMOpts.maxStep)
    {
      //Scale step size
      optVecR *= (QMMMOpts.maxStep/vecMax);
    }
    vecMax = optVecP.norm()/QMMMOpts.NBeads;
    if (vecMax > QMMMOpts.maxStep)
    {
      //Scale step size
      optVecP *= (QMMMOpts.maxStep/vecMax);
    }
    //Update positions
    #pragma omp parallel for schedule(dynamic)
    for (int p=pathStart;p<pathEnd;p++)
    {
      //Loop over all beads
      if (p == 0)
      {
        //Reactant
        int ct = 0;
        for (int i=0;i<Natoms;i++)
        {
          //Move QM atoms
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            QMMMData[i].P[p].x += optVecR(ct);
            QMMMData[i].P[p].y += optVecR(ct+1);
            QMMMData[i].P[p].z += optVecR(ct+2);
            ct += 3;
          }
        }
      }
      else if (p == (QMMMOpts.NBeads-1))
      {
        //Product
        int ct = 0;
        for (int i=0;i<Natoms;i++)
        {
          //Move QM atoms
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            QMMMData[i].P[p].x += optVecP(ct);
            QMMMData[i].P[p].y += optVecP(ct+1);
            QMMMData[i].P[p].z += optVecP(ct+2);
            ct += 3;
          }
        }
      }
      else
      {
        //Path
        int gfID; //Location in the global forces array
        gfID = p-1; //Bead ID
        gfID *= Ndof; //Number of QM/PB atoms
        int ct = 0; //Counter
        for (int i=0;i<Natoms;i++)
        {
          //Move QM atoms
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            QMMMData[i].P[p].x += optVecG(gfID+ct);
            QMMMData[i].P[p].y += optVecG(gfID+ct+1);
            QMMMData[i].P[p].z += optVecG(gfID+ct+2);
            ct += 3;
          }
        }
      }
    }
    //Print structure
    Print_traj(QMMMData,qmFile,QMMMOpts);
    //Calculate new forces
    gForces.setZero(); //Remove old forces
    rForces.setZero(); //Remove old forces (reactant)
    pForces.setZero(); //Remove old forces (product)
    newTS = 0; //Storage for new TS ID
    newTSEnergy = -hugeNum;
    for (int p=pathStart;p<pathEnd;p++)
    {
      double E = 0;
      double Eqm = 0;
      double Emm = 0;
      //Erase old forces
      forces.setZero();
      //Calculate forces (QM part)
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        Eqm += GaussianForces(QMMMData,forces,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        Eqm += PSI4Forces(QMMMData,forces,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        Eqm += NWChemForces(QMMMData,forces,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      E += Eqm; //Save the partial energy
      //Calculate forces (MM part)
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        E += TINKERForces(QMMMData,forces,QMMMOpts,p);
        if (AMOEBA || QMMMOpts.useImpSolv)
        {
          //Forces from MM polarization
          E += TINKERPolForces(QMMMData,forces,QMMMOpts,p);
        }
        Emm += TINKEREnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        E += LAMMPSForces(QMMMData,forces,QMMMOpts,p);
        Emm += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      //Save total energy
      sumE += E;
      //Update energies
      if (p == 0)
      {
        //Reactant energy
        if ((Eqm+Emm) > QMMMOpts.EReact)
        {
          //Force the Hessian to be rebuilt
          EOld = -hugeNum;
        }
        QMMMOpts.EReact = Eqm+Emm;
      }
      else if (p == (QMMMOpts.NBeads-1))
      {
        //Product energy
        if ((Eqm+Emm) > QMMMOpts.EProd)
        {
          //Force the Hessian to be rebuilt
          EOld = -hugeNum;
        }
        QMMMOpts.EProd = Eqm+Emm;
      }
      //Modify forces along the tangent
      VectorXd QMTangent(Ndof); //Tangent vector
      VectorXd distp1(Ndof); //Displacement for p+1
      VectorXd distm1(Ndof); //Displacement for p-1
      if ((p != 0) && (p != (QMMMOpts.NBeads-1)))
      {
        //Calculate tangents for middle replicas
        MatrixXd geom1((Nqm+Npseudo),3); //Current replica
        MatrixXd geom2((Nqm+Npseudo),3); //Second replica
        MatrixXd geom3((Nqm+Npseudo),3); //Third replica
        //Save geometries
        int ct = 0; //Reset counter for the number of atoms
        for (int i=0;i<Natoms;i++)
        {
          //Only include QM and PB regions
          if (oldQMMMData[i].QMRegion || oldQMMMData[i].PBRegion)
          {
            //Save current replica
            geom1(ct,0) = oldQMMMData[i].P[p].x;
            geom1(ct,1) = oldQMMMData[i].P[p].y;
            geom1(ct,2) = oldQMMMData[i].P[p].z;
            //Save replica p+1
            geom2(ct,0) = oldQMMMData[i].P[p+1].x;
            geom2(ct,1) = oldQMMMData[i].P[p+1].y;
            geom2(ct,2) = oldQMMMData[i].P[p+1].z;
            //Save replica p-1
            geom3(ct,0) = oldQMMMData[i].P[p-1].x;
            geom3(ct,1) = oldQMMMData[i].P[p-1].y;
            geom3(ct,2) = oldQMMMData[i].P[p-1].z;
            //Advance counter
            ct += 1;
          }
        }
        //Calculate displacements
        distp1 = KabschDisplacement(geom1,geom2,(Nqm+Npseudo));
        distp1 *= -1; //Change direction
        distm1 = KabschDisplacement(geom1,geom3,(Nqm+Npseudo));
        //Remove inactive atoms
        ct = 0; //Reset counter for the number of atoms
        for (int i=0;i<Natoms;i++)
        {
          //Only include QM and PB regions
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            //Only include active atoms in the tangent
            if (!QMMMData[i].NEBActive)
            {
              //Delete distance components
              distp1(ct) = 0;
              distp1(ct+1) = 0;
              distp1(ct+2) = 0;
              distm1(ct) = 0;
              distm1(ct+1) = 0;
              distm1(ct+2) = 0;
            }
            //Advance counter
            ct += 3;
          }
        }
        //Calculate tangent
        QMTangent = CINEBTangent(distp1,distm1,QMMMOpts,p);
      }
      else
      {
        //No tangent for the reactant and product
        QMTangent.setZero();
        //Get rid of displacements for spring force calculation
        distp1.setZero();
        distm1.setZero();
      }
      //Project forces onto tangent vector
      double forcesDotTan = forces.dot(QMTangent); //Overlap of vectors
      //Remove forces along the tangent
      forces -= forcesDotTan*QMTangent; //End points have no tangents
      if ((p == QMMMOpts.TSBead) && QMMMOpts.climb)
      {
        //Climbing image for TS after 10 steps
        forces -= forcesDotTan*QMTangent;
        //Get rid of displacements for spring force calculations
        distp1.setZero();
        distm1.setZero();
      }
      if (p == pathStart)
      {
        //Smoothly transition to climbing image NEB
        if (stepCt == 8)
        {
          //Turn on climbing image NEB for the next iteration
          QMMMOpts.climb = 1;
        }
        else if ((stepCt == 9) && (optCt == 0))
        {
          //Notify users and change Hessian settings
          cout << "    Starting climbing image NEB...";
          cout << '\n';
        }
      }
      //Add spring forces (Disti1.norm()=0 at end points and TS)
      double springDist = distp1.norm()-distm1.norm();
      forces += (QMMMOpts.kSpring*springDist*QMTangent);
      //Add forces to global arrays
      if (p == 0)
      {
        //Reactant
        #pragma omp parallel for schedule(dynamic)
        for (int i=0;i<Ndof;i++)
        {
          rForces(i) += forces(i);
        }
      }
      else if (p == (QMMMOpts.NBeads-1))
      {
        //Product
        #pragma omp parallel for schedule(dynamic)
        for (int i=0;i<Ndof;i++)
        {
          pForces(i) += forces(i);
        }
      }
      else
      {
        //Path
        int gfID; //Location of the bead in the global array
        gfID = Ndof;
        gfID *= p-1;
        #pragma omp parallel for schedule(dynamic)
        for (int i=0;i<Ndof;i++)
        {
          gForces(gfID+i) += forces(i);
        }
      }
      //Update statistics for convergence testing
      double maxForce;
      maxForce = abs(forces.maxCoeff());
      if (abs(forces.minCoeff()) > maxForce)
      {
        //Update max
        maxForce = abs(forces.minCoeff());
      }
      //Save statistics
      if ((Eqm+Emm) > newTSEnergy)
      {
        //Assuming Reactant->Product; Puts the TS on the low energy side
        newTSEnergy = Eqm+Emm; //New energy
        newTS = p; //New TS
      }
      forceStats(p,0) = maxForce;
      forceStats(p,1) = forces.squaredNorm(); //RMS force
    }
    //Check for unstable optimization vectors
    int ct = 0; //Use a counter as a safe way to check all replicas
    #pragma omp parallel for schedule(dynamic) reduction(+:ct)
    for (int p=pathStart;p<pathEnd;p++)
    {
      //Check if the beads are moving in the correct directions
      double vecDotForce; //Dot product of the OptVec and Forces
      vecDotForce = 0; //Set dot product equal to zero
      double normForce; //Local norm of the forces
      normForce = 0; //Set norm equal to zero
      double localMaxForce; //Local maximum force
      localMaxForce = 0; //Set max equal to zero
      if (p == 0)
      {
        //Reactant
        for (int i=0;i<Ndof;i++)
        {
          //Update the dot products
          vecDotForce += optVecR(i)*rForces(i);
          normForce += rForces(i)*rForces(i);
          if (localMaxForce < abs(rForces(i)))
          {
            localMaxForce = abs(rForces(i));
          }
        }
      }
      else if (p == (QMMMOpts.NBeads-1))
      {
        //Product
        for (int i=0;i<Ndof;i++)
        {
          //Update the dot products
          vecDotForce += optVecP(i)*pForces(i);
          normForce += pForces(i)*pForces(i);
          if (localMaxForce < abs(pForces(i)))
          {
            localMaxForce = abs(pForces(i));
          }
        }
      }
      else
      {
        //Path
        int gfID; //Location in arrays
        gfID = Ndof; //Adjust for number of QM and PB atoms
        gfID *= p-1; //Adjust for the bead ID
        for (int i=0;i<Ndof;i++)
        {
          //Update the dot products
          vecDotForce += optVecG(gfID+i)*gForces(gfID+i);
          normForce += gForces(gfID+i)*gForces(gfID+i);
          if (localMaxForce < abs(gForces(gfID+i)))
          {
            localMaxForce = abs(gForces(gfID+i));
          }
        }
      }
      normForce /= Ndof; //Make the dot product a norm
      normForce = sqrt(normForce); //Make the norm an RMS value
      if (((vecDotForce < 0) || (localMaxForce >= 1.0)) &&
         (normForce > RMSFTol) && (localMaxForce > maxFTol))
      {
        //Bead is moving in the wrong direction and is not converged
        ct += 1; //Log the instability
      }
    }
    if (ct > 0)
    {
      //One or more beads are moving in the wrong direction
      EOld = -hugeNum; //Force the Hessian to be rebuilt
    }
    //Update TS ID and energy
    if (newTS != QMMMOpts.TSBead)
    {
      //The TS has moved
      EOld = -hugeNum; //Force the Hessian to be rebuilt
    }
    QMMMOpts.ETrans = newTSEnergy;
    QMMMOpts.TSBead = newTS;
    //Update Hessian
    gradDiffG -= gForces; //Path
    gradDiffR -= rForces; //Reactant
    gradDiffP -= pForces; //Product
    if (((stepCt%25) == 0) || (stepCt < 15))
    {
      //Build a new Hessian after 30 steps
      cout << "    Performing a steepest descent step...";
      cout << '\n';
      //Shrink step size
      if (stepCt < 15)
      {
        //Reduce stepsize
        stepScale = sdScale*QMMMOpts.stepScale; //Small step
      }
      else
      {
        //Reduce step size further
        stepScale *= 0.75;
      }
      //Create new Hessian as an identity matrix
      iHessG.setIdentity(); //Already an "inverse" Hessian
      iHessR.setIdentity(); //Already an "inverse" Hessian (reactant)
      iHessP.setIdentity(); //Already an "inverse" Hessian (product)
    }
    else if (((stepCt+1)%25) == 0)
    {
      //Prepare for the upcoming SD step
      cout << "    Reducing the step size...";
      cout << '\n';
      //Shrink step size
      if (stepScale > (sdScale*QMMMOpts.stepScale))
      {
        //Reduce step size
        stepScale = sdScale*QMMMOpts.stepScale;
      }
      else if (stepScale > (0.25*sdScale*QMMMOpts.stepScale))
      {
        //Reduce step size further
        stepScale *= 0.75;
      }
      else
      {
        //Minimum step size in case TS is climbing
        stepScale = 0.25*sdScale*QMMMOpts.stepScale;
      }
      //Create new Hessian as an identity matrix
      iHessG.setIdentity(); //Already an "inverse" Hessian
      iHessR.setIdentity(); //Already an "inverse" Hessian (reactant)
      iHessP.setIdentity(); //Already an "inverse" Hessian (product)
    }
    else if (EOld != (-hugeNum))
    {
      //Update Hessian
      cout << "    Updating inverse Hessian...";
      cout << '\n';
      //Start really long "line" (path)
      iHessG = iHessG+((optVecG*optVecG.transpose())/(optVecG.transpose()
      *gradDiffG))-((iHessG*gradDiffG*gradDiffG.transpose()*iHessG)
      /(gradDiffG.transpose()*iHessG*gradDiffG));
      //End really long "line" (path)
      //Start really long "line" (reactant)
      iHessR = iHessR+((optVecR*optVecR.transpose())/(optVecR.transpose()
      *gradDiffR))-((iHessR*gradDiffR*gradDiffR.transpose()*iHessR)
      /(gradDiffR.transpose()*iHessR*gradDiffR));
      //End really long "line" (reactant)
      //Start really long "line" (product)
      iHessP = iHessP+((optVecP*optVecP.transpose())/(optVecP.transpose()
      *gradDiffP))-((iHessP*gradDiffP*gradDiffP.transpose()*iHessP)
      /(gradDiffP.transpose()*iHessP*gradDiffP));
      //End really long "line" (product)
      //Increase stepsize for the next iteration
      stepScale *= 1.20; //Does not reach the full StepScale by 25 steps
      if (stepScale > QMMMOpts.stepScale)
      {
        //Prevent step size from getting too large
        stepScale = QMMMOpts.stepScale;
      }
    }
    else
    {
      //Take a small steepest descent step and rebuild Hessian
      cout << "    Potentially unstable path. Constructing a new Hessian...";
      cout << '\n';
      //Reduce step size
      if (stepScale > (sdScale*QMMMOpts.stepScale))
      {
        stepScale = sdScale*QMMMOpts.stepScale;
      }
      else if (stepScale > (0.25*sdScale*QMMMOpts.stepScale))
      {
        //Reduce step size further
        stepScale *= 0.75;
      }
      else
      {
        //Minimum step size in case TS is climbing
        stepScale = 0.25*sdScale*QMMMOpts.stepScale;
      }
      //Construct new Hessian
      iHessG.setIdentity(); //Already an "inverse" Hessian (path)
      iHessR.setIdentity(); //Already an "inverse" Hessian (reactant)
      iHessP.setIdentity(); //Already an "inverse" Hessian (product)
    }
    //Update old energy
    EOld = sumE;
    //Check convergence
    stepCt += 1;
    pathDone = PathConverged(QMMMData,oldQMMMData,forceStats,stepCt,QMMMOpts,1);
  }
  //Clean up files
  call.str("");
  call << "rm -f LICHMNEBOpt.xyz MMCharges_*.txt";
  globalSys = system(call.str().c_str());
  //Finish and return
  return;
};

//Path ensemble samping routines
int FBNEBMCMove(vector<QMMMAtom>& QMMMData, vector<VectorXd>& allForces,
                QMMMSettings& QMMMOpts, VectorXd& Emc)
{
  //Function to try a force-bias Monte Carlo move
  int acc = 0; //Accept or reject
  double randNum; //Random number
  int Ndof = 3*(Nqm+Npseudo); //Number of QM degrees of freedom
  VectorXd allEnergies(QMMMOpts.NBeads); //Energies of individual beads
  VectorXd randNums(QMMMOpts.NBeads); //Array to store random numbers
  //Save structure and forces for rejected moves
  VectorXd oldEnergies = Emc;
  vector<QMMMAtom> oldQMMMData = QMMMData;
  vector<VectorXd> oldForces = allForces;
  vector<VectorXd> tempForces = allForces; //Includes tangents
  //Add tangent forces
  if (QMMMOpts.NBeads > 1)
  {
    //Only calculate tangent forces for multireplica simulations
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Calculate QM tangent forces
      int ct; //Generic counter
      double forcesDotTan; //Overlap of vectors
      double springDist; //Distance between the neighboring beads
      //Read old forces in the QM region
      VectorXd forces(Ndof); //Temporary force array
      forces.setZero();
      ct = 0;
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
        {
          //X component
          forces(ct) = tempForces[p](3*i);
          ct += 1;
          //Y component
          forces(ct) = tempForces[p](3*i+1);
          ct += 1;
          //Z component
          forces(ct) = tempForces[p](3*i+2);
          ct += 1;
        }
      }
      VectorXd QMTangent(Ndof); //Tangent vector
      VectorXd distp1(Ndof); //Displacement for p+1
      VectorXd distm1(Ndof); //Displacement for p-1
      if ((p != 0) && (p != (QMMMOpts.NBeads-1)))
      {
        //Calculate tangents for middle replicas
        MatrixXd geom1((Nqm+Npseudo),3); //Current replica
        MatrixXd geom2((Nqm+Npseudo),3); //Second replica
        MatrixXd geom3((Nqm+Npseudo),3); //Third replica
        //Save geometries
        int ct = 0; //Reset counter for the number of atoms
        for (int i=0;i<Natoms;i++)
        {
          //Only include QM and PB regions
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            //Save current replica
            geom1(ct,0) = QMMMData[i].P[p].x;
            geom1(ct,1) = QMMMData[i].P[p].y;
            geom1(ct,2) = QMMMData[i].P[p].z;
            //Save replica p+1
            geom2(ct,0) = QMMMData[i].P[p+1].x;
            geom2(ct,1) = QMMMData[i].P[p+1].y;
            geom2(ct,2) = QMMMData[i].P[p+1].z;
            //Save replica p-1
            geom3(ct,0) = QMMMData[i].P[p-1].x;
            geom3(ct,1) = QMMMData[i].P[p-1].y;
            geom3(ct,2) = QMMMData[i].P[p-1].z;
            //Advance counter
            ct += 1;
          }
        }
        //Calculate displacements
        distp1 = KabschDisplacement(geom1,geom2,(Nqm+Npseudo));
        distp1 *= -1; //Change direction
        distm1 = KabschDisplacement(geom1,geom3,(Nqm+Npseudo));
        //Calculate tangent
        QMTangent = NEBTangent(distp1,distm1,QMMMOpts,p);
      }
      else
      {
        //No tangent for the reactant and product
        QMTangent.setZero();
        //Get rid of displacements for spring force calculation
        distp1.setZero();
        distm1.setZero();
      }
      //Project forces onto tangent vector
      forcesDotTan = forces.dot(QMTangent); //Overlap of vectors
      //Remove forces along the tangent
      forces -= forcesDotTan*QMTangent; //End points have no tangents
      //Add spring forces (Distp1.norm()=0 and Distm1.norm()=0 at end points)
      springDist = distp1.norm()-distm1.norm(); //Combined distance
      forces += (QMMMOpts.kSpring*springDist*QMTangent);
      //Save modified forces
      ct = 0;
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
        {
          //X component
          tempForces[p](3*i) = forces(ct);
          ct += 1;
          //Y component
          tempForces[p](3*i+1) = forces(ct);
          ct += 1;
          //Z component
          tempForces[p](3*i+2) = forces(ct);
          ct += 1;
        }
      }
    }
  }
  //Randomly move along forces forces
  randNums.setZero();
  vector<VectorXd> randomNoise; //Adds random noise to MC moves
  for (int p=0;p<QMMMOpts.NBeads;p++)
  {
    //Separate the serial random numbers
    randNum = (((double)rand())/((double)RAND_MAX)); //Number between 0 and 1
    randNum *= 4.0; //Now between 0 and 4; Average: 2.0
    randNum -= 2.0; //Now between -2 and 2; Abs. average: 1.0
    randNums(p) = randNum*mcStep; //Scale displacement
    //Create random noise
    VectorXd tempRandNoise(3*Natoms);
    tempRandNoise.setRandom();
    if (tempRandNoise.squaredNorm() != 0)
    {
      //Normalize and scale vector
      tempRandNoise.normalize();
      tempRandNoise *= 0.05; //Noise is 5% of the forces
    }
    randomNoise.push_back(tempRandNoise);
  }
  #pragma omp parallel for schedule(dynamic)
  for (int p=0;p<QMMMOpts.NBeads;p++)
  {
    //Normalize and avoid divide by zero errors
    if (tempForces[p].squaredNorm() != 0)
    {
      //Normalize the vector
      tempForces[p].normalize();
    }
    //Scale vector
    tempForces[p] += randomNoise[p]; //Add noise vector
    tempForces[p] *= randNums(p); //Scale forces
    //Update postions
    for (int i=0;i<Natoms;i++)
    {
      QMMMData[i].P[p].x += tempForces[p](3*i);
      QMMMData[i].P[p].y += tempForces[p](3*i+1);
      QMMMData[i].P[p].z += tempForces[p](3*i+2);
    }
    //Clear array
    allForces[p].setZero();
  }
  //Fix parallel for classical MC
  int mcThreads = Nthreads;
  if (QMMMOpts.NBeads == 1)
  {
    mcThreads = 1;
  }
  //Calculate new energies and forces
  allEnergies.setZero();
  #pragma omp parallel for schedule(dynamic) num_threads(mcThreads) \
          reduction(+:QMTime,MMTime)
  for (int p=0;p<QMMMOpts.NBeads;p++)
  {
    //Calculate the energy of bead p
    double E = 0;
    double Emm = 0;
    double Eqm = 0;
    //Create blank force array
    VectorXd forces(Ndof);
    forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      Eqm += GaussianForces(QMMMData,forces,QMMMOpts,p);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      Eqm += PSI4Forces(QMMMData,forces,QMMMOpts,p);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      Eqm += NWChemForces(QMMMData,forces,QMMMOpts,p);
      QMTime += (unsigned)time(0)-tStart;
    }
    E += Eqm;
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,forces,QMMMOpts,p);
      if (AMOEBA || QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,forces,QMMMOpts,p);
      }
      Emm += TINKEREnergy(QMMMData,QMMMOpts,p);
      MMTime += (unsigned)time(0)-tStart;
    }
    allEnergies(p) = Eqm+Emm; //Save the energy for calculating statistics
    //Save QM forces to global array
    int ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
      {
        //Add X component
        allForces[p](3*i) = forces(ct);
        ct += 1;
        //Add Y component
        allForces[p](3*i+1) = forces(ct);
        ct += 1;
        //Add Z component
        allForces[p](3*i+2) = forces(ct);
        ct += 1;
      }
    }
    //Add MM forces
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E = TINKERMMForces(QMMMData,allForces[p],QMMMOpts,p);
      MMTime += (unsigned)time(0)-tStart;
    }
  }
  //Accept or reject
  randNums.setZero();
  for (int p=0;p<QMMMOpts.NBeads;p++)
  {
    //Separate the serial random numbers
    randNum = (((double)rand())/((double)RAND_MAX)); //Number between 0 and 1
    randNums(p) = randNum;
  }
  #pragma omp parallel for schedule(dynamic) reduction(+:acc)
  for (int p=0;p<QMMMOpts.NBeads;p++)
  {
    double dE = QMMMOpts.beta*(allEnergies(p)-oldEnergies(p));
    double prob = exp(-dE);
    if (randNums(p) <= prob)
    {
      //Accept and keep structure/forces
      acc += 1;
    }
    else
    {
      //Reject and revert to old structure and forces
      allEnergies(p) = oldEnergies(p);
      allForces[p] = oldForces[p];
      for (int i=0;i<Natoms;i++)
      {
        //Revert to old structure and forces for bead p
        QMMMData[i].P[p].x = oldQMMMData[i].P[p].x;
        QMMMData[i].P[p].y = oldQMMMData[i].P[p].y;
        QMMMData[i].P[p].z = oldQMMMData[i].P[p].z;
      }
    }
  }
  //Update energies and return number of decisions
  Emc = allEnergies;
  return acc;
};

//End of file group
///@}

