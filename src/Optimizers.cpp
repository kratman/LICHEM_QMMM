/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 A set of optimization routines for the QM part of QMMM calculculations.
 MM optimizations are performed with the native optimizers in the MM
 wrappers.

 Reference for optimization routines:
 Press et al., Numerical Recipes 3nd Edition, (2007)

*/

//Convergence test functions
bool OptConverged(vector<QMMMAtom>& QMMMData, vector<QMMMAtom>& oldQMMMData,
                  VectorXd& forces, int stepCt, QMMMSettings& QMMMOpts,
                  int bead, bool QMRegion)
{
  //Check convergence of QMMM optimizations
  bool optDone = 0; //Ends the simulation
  double RMSDiff = 0; //RMS deviation
  double RMSForce = 0; //RMS force
  double maxForce = 0; //Maximum force
  double sumE = 0; //Storage for energies
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Convergence criteria
  double maxFTol = 20*QMMMOpts.QMOptTol; //Opt. tolerance for max. force
  double RMSFTol = 10*QMMMOpts.QMOptTol; //Opt. tolerance for RMS force
  //Check progress
  if (QMRegion)
  {
    //Check if a QM calculation is converged
    maxForce = abs(forces.maxCoeff());
    if (maxForce < abs(forces.minCoeff()))
    {
      maxForce = abs(forces.minCoeff());
    }
    RMSForce = sqrt(forces.squaredNorm()/Ndof);
    #pragma omp parallel for schedule(dynamic) reduction(+:RMSDiff)
    for (int i=0;i<Natoms;i++)
    {
      //Calculate QM-QM distance matrix
      double RMSTemp = 0; //Store a local sum
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        for (int j=0;j<i;j++)
        {
          if (QMMMData[j].QMRegion or QMMMData[j].PBRegion)
          {
            double RNew = 0;
            double ROld = 0;
            RNew = CoordDist2(QMMMData[i].P[bead],
                              QMMMData[j].P[bead]).vecMag();
            ROld = CoordDist2(oldQMMMData[i].P[bead],
                              oldQMMMData[j].P[bead]).vecMag();
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
    RMSDiff /= (Nqm+Npseudo)*(Nqm+Npseudo-1)/2;
    RMSDiff = sqrt(RMSDiff);
    //Print progress
    cout << "    QM step: " << stepCt;
    cout << " | RMS dev: " << LICHEMFormFloat(RMSDiff,12);
    cout << " \u212B" << '\n';
    cout << "    Max. force: " << LICHEMFormFloat(maxForce,12);
    cout << " eV/\u212B | RMS force: " << LICHEMFormFloat(RMSForce,12);
    cout << " eV/\u212B" << '\n';
    //Check convergence criteria
    if ((RMSDiff <= QMMMOpts.QMOptTol) and (RMSForce <= RMSFTol) and
       (maxForce <= maxFTol))
    {
      optDone = 1;
      cout << "    QM optimization complete." << '\n';
    }
    cout << '\n';
    cout.flush();
  }
  if (!QMRegion)
  {
    //Check energy and convergence of the whole system
    sumE = 0; //Reinitialize the energy
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      sumE += GaussianEnergy(QMMMData,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      sumE += PSI4Energy(QMMMData,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      sumE += NWChemEnergy(QMMMData,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      sumE += TINKEREnergy(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      sumE += LAMMPSEnergy(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Calculate RMS displacement (distance matrix)
    #pragma omp parallel for schedule(dynamic) reduction(+:RMSDiff)
    for (int i=0;i<Natoms;i++)
    {
      double RMSTemp = 0; //Store a local sum
      for (int j=0;j<i;j++)
      {
        double RNew = 0;
        double ROld = 0;
        RNew = CoordDist2(QMMMData[i].P[bead],
                          QMMMData[j].P[bead]).vecMag();
        ROld = CoordDist2(oldQMMMData[i].P[bead],
                          oldQMMMData[j].P[bead]).vecMag();
        RNew = sqrt(RNew);
        ROld = sqrt(ROld);
        //Update local sum
        RMSTemp += (RNew-ROld)*(RNew-ROld);
      }
      //Update sum
      RMSDiff += RMSTemp;
    }
    RMSDiff /= (Natoms-Nfreeze)*(Natoms-Nfreeze-1)/2;
    RMSDiff = sqrt(RMSDiff);
    //Print progress
    cout << " | Opt. step: ";
    cout << stepCt << " | Energy: ";
    cout << LICHEMFormFloat(sumE,16) << " eV ";
    cout << " | RMS dev: " << LICHEMFormFloat(RMSDiff,12);
    cout << " \u212B" << '\n';
    //Check convergence
    if (RMSDiff <= QMMMOpts.MMOptTol)
    {
      optDone = 1;
      if (QMMM and (stepCt > 1))
      {
        cout << "    QMMM relaxation satisfactory.";
        cout << '\n';
      }
    }
    //Flush output
    cout.flush();
  }
  return optDone;
};

//Optimizer functions
void LICHEMSteepest(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
{
  //Cartesian steepest descent optimizer
  stringstream call; //Stream for system calls and reading/writing files
  int stepCt = 0; //Counter for optimization steps
  fstream qmFile, inFile, outFile; //Generic file names
  double EOld = 0; //Old saved energy
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize charges
  if (Nmm > 0)
  {
    WriteChargeFile(QMMMData,QMMMOpts,bead);
  }
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << bead << ".xyz";
  qmFile.open(call.str().c_str(),ios_base::out);
  //Initialize optimization variables
  double stepSize = 1;
  double vecMax = 0;
  bool optDone = 0;
  vector<QMMMAtom> oldQMMMData = QMMMData; //Previous structure
  //Run optimization
  double stepScale = QMMMOpts.stepScale;
  stepScale *= 0.70; //Take a smaller first step
  while ((!optDone) and (stepCt < QMMMOpts.maxOptSteps))
  {
    double E = 0;
    //Create blank force array
    VectorXd forces(Ndof);
    forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      E += GaussianForces(QMMMData,forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      E += PSI4Forces(QMMMData,forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      E += NWChemForces(QMMMData,forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,forces,QMMMOpts,bead);
      if (AMOEBA or QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,forces,QMMMOpts,bead);
      }
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      E += LAMMPSForces(QMMMData,forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Check step size
    if (E > EOld)
    {
      //Take smaller steps if the energy does not improve
      cout << "    Energy did not decrease. Reducing the step size...";
      cout << '\n';
      stepScale *= 0.60; //Reduce step size
    }
    //Save structure and energy
    EOld = E;
    oldQMMMData = QMMMData;
    //Check optimization step size
    vecMax = abs(forces.maxCoeff());
    if (abs(forces.minCoeff()) > vecMax)
    {
      vecMax = stepScale*abs(forces.minCoeff());
    }
    else
    {
      vecMax *= stepScale;
    }
    stepSize = stepScale;
    if (vecMax > QMMMOpts.maxStep)
    {
      //Scale step size
      cout << "    Scaling step size to match the maximum...";
      cout << '\n';
      stepSize *= (QMMMOpts.maxStep/vecMax);
    }
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        QMMMData[i].P[bead].x += stepSize*forces(ct);
        QMMMData[i].P[bead].y += stepSize*forces(ct+1);
        QMMMData[i].P[bead].z += stepSize*forces(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(QMMMData,qmFile,QMMMOpts);
    //Check convergence
    optDone = OptConverged(QMMMData,oldQMMMData,forces,stepCt,QMMMOpts,bead,1);
    stepCt += 1;
    //Increase step size
    stepScale *= 1.05;
    if (stepScale > QMMMOpts.stepScale)
    {
      //Prevent step size from getting too large
      stepScale = QMMMOpts.stepScale;
    }
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << bead << ".xyz";
  call << " MMCharges_" << bead << ".txt";
  globalSys = system(call.str().c_str());
  //Finish and return
  return;
};

void LICHEMDFP(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
{
  //A simple Davidon-Fletcher-Powell optimizer
  //NB: This optimizer does not have a true line search, instead
  //a steepest descent step is performed if the optimizer is unstable
  stringstream call; //Stream for system calls and reading/writing files
  int stepCt = 0; //Counter for optimization steps
  fstream qmFile,inFile,outFile; //Generic file streams
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  double RMSFTol = 10*QMMMOpts.QMOptTol; //Opt. tolerance for RMS force
  //Initialize charges
  if (Nmm > 0)
  {
    WriteChargeFile(QMMMData,QMMMOpts,bead);
  }
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << bead << ".xyz";
  qmFile.open(call.str().c_str(),ios_base::out);
  //Create DFP arrays
  VectorXd optVec(Ndof); //Gradient descent direction
  VectorXd gradDiff(Ndof); //Change in the gradient
  VectorXd forces(Ndof); //Forces
  MatrixXd iHess(Ndof,Ndof); //Inverse Hessian
  //Initialize arrays
  optVec.setZero();
  gradDiff.setZero();
  forces.setZero();
  //Create an identity matrix as the initial Hessian
  iHess.setIdentity(); //Already an "inverse" Hessian
  //Initialize optimization variables
  double stepScale; //Local copy
  double E = 0; //Energy
  double EOld = 0; //Energy from previous step
  double sdScale = 0.01; //Scale factor for SD steps
  double vecMax = 0; //Maxium atomic displacement
  bool optDone = 0; //Flag to end the optimization
  //Calculate forces (QM part)
  if (Gaussian)
  {
    int tStart = (unsigned)time(0);
    E += GaussianForces(QMMMData,forces,QMMMOpts,bead);
    QMTime += (unsigned)time(0)-tStart;
  }
  if (PSI4)
  {
    int tStart = (unsigned)time(0);
    E += PSI4Forces(QMMMData,forces,QMMMOpts,bead);
    QMTime += (unsigned)time(0)-tStart;
    //Delete annoying useless files
    globalSys = system("rm -f psi.* timer.*");
  }
  if (NWChem)
  {
    int tStart = (unsigned)time(0);
    E += NWChemForces(QMMMData,forces,QMMMOpts,bead);
    QMTime += (unsigned)time(0)-tStart;
  }
  //Calculate forces (MM part)
  if (TINKER)
  {
    int tStart = (unsigned)time(0);
    E += TINKERForces(QMMMData,forces,QMMMOpts,bead);
    if (AMOEBA or QMMMOpts.useImpSolv)
    {
      //Forces from MM polarization
      E += TINKERPolForces(QMMMData,forces,QMMMOpts,bead);
    }
    MMTime += (unsigned)time(0)-tStart;
  }
  if (LAMMPS)
  {
    int tStart = (unsigned)time(0);
    E += LAMMPSForces(QMMMData,forces,QMMMOpts,bead);
    MMTime += (unsigned)time(0)-tStart;
  }
  //Output initial RMS force
  vecMax = 0; //Using this variable to avoid creating a new one
  vecMax = forces.squaredNorm(); //Calculate initial RMS force
  vecMax = sqrt(vecMax/Ndof);
  cout << "    Performing a steepest descent step..." << '\n';
  cout << "    QM step: 0";
  cout << " | RMS force: " << LICHEMFormFloat(vecMax,12);
  cout << " eV/\u212B";
  cout << '\n' << '\n';
  cout.flush();
  //Optimize structure
  EOld = E; //Save energy
  stepScale = QMMMOpts.stepScale;
  stepScale *= sdScale; //Take a very small first step
  while ((!optDone) and (stepCt < QMMMOpts.maxOptSteps))
  {
    E = 0; // Reinitialize energy
    //Copy old structure and old forces
    vector<QMMMAtom> oldQMMMData = QMMMData;
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<Ndof;i++)
    {
      //Reinitialize the change in the gradient
      gradDiff(i) = forces(i);
    }
    //Determine new structure
    optVec = iHess*forces;
    optVec *= stepScale;
    //Check step size
    vecMax = optVec.norm();
    if (vecMax > QMMMOpts.maxStep)
    {
      //Scale step size
      optVec *= (QMMMOpts.maxStep/vecMax);
    }
    //Update positions
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        QMMMData[i].P[bead].x += optVec(ct);
        QMMMData[i].P[bead].y += optVec(ct+1);
        QMMMData[i].P[bead].z += optVec(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(QMMMData,qmFile,QMMMOpts);
    //Calculate forces (QM part)
    forces.setZero();
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      E += GaussianForces(QMMMData,forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      E += PSI4Forces(QMMMData,forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      E += NWChemForces(QMMMData,forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,forces,QMMMOpts,bead);
      if (AMOEBA or QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,forces,QMMMOpts,bead);
      }
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      E += LAMMPSForces(QMMMData,forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Check stability
    double vecDotForces; //Dot product of the forces and optimization vector
    vecDotForces = optVec.dot(forces);
    double normForce; //Local norm of the forces
    normForce = forces.norm(); //Take the norm of the foces
    normForce /= sqrt(Ndof); //Make the norm an RMS value
    double localMaxForce; //Maximum value in the force array
    localMaxForce = abs(forces.maxCoeff());
    if (((vecDotForces < 0) or (localMaxForce >= 1.0)) and
       (normForce > RMSFTol))
    {
      //Optimizer is going the wrong direction and is not converged
      EOld = -1*hugeNum; //Force the Hessian to be rebuilt
    }
    //Update Hessian
    gradDiff -= forces;
    if (((stepCt%30) == 0) or (stepCt < 15))
    {
      //Build a new Hessian after 30 steps
      cout << "    Performing a steepest descent step...";
      cout << '\n';
      //Shrink step size
      if ((stepCt < 15) and (E < EOld))
      {
        //Reduce step size
        stepScale = sdScale*QMMMOpts.stepScale; //Small step
      }
      else
      {
        //Reduce step size further
        stepScale *= 0.75;
      }
      //Check for minimum step size
      if (stepScale < (0.25*sdScale*QMMMOpts.stepScale))
      {
        stepScale = 0.25*sdScale*QMMMOpts.stepScale;
      }
      //Create new Hessian as an identity matrix
      iHess.setIdentity(); //Already an "inverse" Hessian
    }
    else if (((stepCt+1)%30) == 0)
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
      else
      {
        //Reduce step size further
        stepScale *= 0.75;
      }
      //Check for minimum step size
      if (stepScale < (0.25*sdScale*QMMMOpts.stepScale))
      {
        stepScale = 0.25*sdScale*QMMMOpts.stepScale;
      }
      //Create new Hessian as an identity matrix
      iHess.setIdentity(); //Already an "inverse" Hessian
    }
    else if (E < EOld)
    {
      //Update Hessian
      cout << "    Updating inverse Hessian...";
      cout << '\n';
      //Start really long "line"
      iHess = iHess+((optVec*optVec.transpose())/(optVec.transpose()
      *gradDiff))-((iHess*gradDiff*gradDiff.transpose()*iHess)
      /(gradDiff.transpose()*iHess*gradDiff));
      //End really long "line"
      //Increase stepsize for the next iteration
      stepScale *= 1.20;
      if (stepScale > QMMMOpts.stepScale)
      {
        //Prevent step size from getting too large
        stepScale = QMMMOpts.stepScale;
      }
    }
    else
    {
      //Take a small steepest descent step and rebuild Hessian
      cout << "    Potentially unstable Hessian.";
      cout << " Constructing new Hessian...";
      cout << '\n';
      //Reduce step size
      if (stepScale > (sdScale*QMMMOpts.stepScale))
      {
        stepScale = sdScale*QMMMOpts.stepScale;
      }
      else
      {
        //Reduce step size further
        stepScale *= 0.75;
      }
      //Check for minimum step size
      if (stepScale < (0.25*sdScale*QMMMOpts.stepScale))
      {
        stepScale = 0.25*sdScale*QMMMOpts.stepScale;
      }
      iHess.setIdentity(); //Already an "inverse" Hessian
    }
    //Save energy
    EOld = E;
    //Check convergence
    stepCt += 1;
    optDone = OptConverged(QMMMData,oldQMMMData,forces,stepCt,QMMMOpts,bead,1);
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << bead << ".xyz";
  call << " MMCharges_" << bead << ".txt";
  globalSys = system(call.str().c_str());
  //Finish and return
  return;
};

