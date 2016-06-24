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
bool OptConverged(vector<QMMMAtom>& QMMMData, vector<QMMMAtom>& OldQMMMData,
                  VectorXd& Forces, int stepct, QMMMSettings& QMMMOpts,
                  int bead, bool QMregion)
{
  //Check convergence of QMMM optimizations
  bool OptDone = 0; //Ends the simulation
  double RMSdiff = 0; //RMS deviation
  double RMSforce = 0; //RMS force
  double MAXforce = 0; //Maximum force
  double SumE = 0; //Storage for energies
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Convergence criteria
  double MaxFTol = 20*QMMMOpts.QMOptTol; //Opt. tolerance for max. force
  double RMSFTol = 10*QMMMOpts.QMOptTol; //Opt. tolerance for RMS force
  //Check progress
  if (QMregion)
  {
    //Check if a QM calculation is converged
    MAXforce = abs(Forces.maxCoeff());
    if (MAXforce < abs(Forces.minCoeff()))
    {
      MAXforce = abs(Forces.minCoeff());
    }
    RMSforce = sqrt(Forces.squaredNorm()/Ndof);
    #pragma omp parallel for schedule(dynamic) reduction(+:RMSdiff)
    for (int i=0;i<Natoms;i++)
    {
      //Calculate QM-QM distance matrix
      double RMStmp = 0; //Store a local sum
      if (QMMMData[i].QMregion or QMMMData[i].PBregion)
      {
        for (int j=0;j<i;j++)
        {
          if (QMMMData[j].QMregion or QMMMData[j].PBregion)
          {
            double Rnew = 0;
            double Rold = 0;
            Rnew = CoordDist2(QMMMData[i].P[bead],
                              QMMMData[j].P[bead]).vecMag();
            Rold = CoordDist2(OldQMMMData[i].P[bead],
                              OldQMMMData[j].P[bead]).vecMag();
            Rnew = sqrt(Rnew);
            Rold = sqrt(Rold);
            //Update local sum
            RMStmp += (Rnew-Rold)*(Rnew-Rold);
          }
        }
      }
      //Update sum
      RMSdiff += RMStmp;
    }
    RMSdiff /= (Nqm+Npseudo)*(Nqm+Npseudo-1)/2;
    RMSdiff = sqrt(RMSdiff);
    //Print progress
    cout << "    QM step: " << stepct;
    cout << " | RMS dev: " << LICHEMFormFloat(RMSdiff,12);
    cout << " \u212B" << '\n';
    cout << "    Max. force: " << LICHEMFormFloat(MAXforce,12);
    cout << " eV/\u212B | RMS force: " << LICHEMFormFloat(RMSforce,12);
    cout << " eV/\u212B" << '\n';
    //Check convergence criteria
    if ((RMSdiff <= QMMMOpts.QMOptTol) and (RMSforce <= RMSFTol) and
       (MAXforce <= MaxFTol))
    {
      OptDone = 1;
      cout << "    QM optimization complete." << '\n';
    }
    cout << '\n';
    cout.flush();
  }
  if (!QMregion)
  {
    //Check energy and convergence of the whole system
    SumE = 0; //Reinitialize the energy
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      SumE += GaussianEnergy(QMMMData,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      SumE += PSI4Energy(QMMMData,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      SumE += NWChemEnergy(QMMMData,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      SumE += TINKEREnergy(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      SumE += AMBEREnergy(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      SumE += LAMMPSEnergy(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Calculate RMS displacement (distance matrix)
    #pragma omp parallel for schedule(dynamic) reduction(+:RMSdiff)
    for (int i=0;i<Natoms;i++)
    {
      double RMStmp = 0; //Store a local sum
      for (int j=0;j<i;j++)
      {
        double Rnew = 0;
        double Rold = 0;
        Rnew = CoordDist2(QMMMData[i].P[bead],
                          QMMMData[j].P[bead]).vecMag();
        Rold = CoordDist2(OldQMMMData[i].P[bead],
                          OldQMMMData[j].P[bead]).vecMag();
        Rnew = sqrt(Rnew);
        Rold = sqrt(Rold);
        //Update local sum
        RMStmp += (Rnew-Rold)*(Rnew-Rold);
      }
      //Update sum
      RMSdiff += RMStmp;
    }
    RMSdiff /= (Natoms-Nfreeze)*(Natoms-Nfreeze-1)/2;
    RMSdiff = sqrt(RMSdiff);
    //Print progress
    cout << " | Opt. step: ";
    cout << stepct << " | Energy: ";
    cout << LICHEMFormFloat(SumE,16) << " eV ";
    cout << " | RMS dev: " << LICHEMFormFloat(RMSdiff,12);
    cout << " \u212B" << '\n';
    //Check convergence
    if (RMSdiff <= QMMMOpts.MMOptTol)
    {
      OptDone = 1;
      if (QMMM and (stepct > 1))
      {
        cout << "    QMMM relaxation satisfactory.";
        cout << '\n';
      }
    }
    //Flush output
    cout.flush();
  }
  return OptDone;
};

//Optimizer functions
void LICHEMSteepest(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
{
  //Cartesian steepest descent optimizer
  stringstream call; //Stream for system calls and reading/writing files
  int stepct = 0; //Counter for optimization steps
  fstream qmfile, inFile, outFile; //Generic file names
  double Eold = 0; //Old saved energy
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize charges
  if (Nmm > 0)
  {
    WriteChargeFile(QMMMData,QMMMOpts,bead);
  }
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize optimization variables
  double stepsize = 1;
  double VecMax = 0;
  bool OptDone = 0;
  vector<QMMMAtom> OldQMMMData = QMMMData; //Previous structure
  //Run optimization
  double StepScale = QMMMOpts.stepScale;
  StepScale *= 0.70; //Take a smaller first step
  while ((!OptDone) and (stepct < QMMMOpts.maxOptSteps))
  {
    double E = 0;
    //Create blank force array
    VectorXd Forces(Ndof);
    Forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      E += GaussianForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      E += PSI4Forces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      E += NWChemForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,Forces,QMMMOpts,bead);
      if (AMOEBA or QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,Forces,QMMMOpts,bead);
      }
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      E += AMBERForces(QMMMData,Forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      E += LAMMPSForces(QMMMData,Forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Check step size
    if (E > Eold)
    {
      //Take smaller steps if the energy does not improve
      cout << "    Energy did not decrease. Reducing the step size...";
      cout << '\n';
      StepScale *= 0.60; //Reduce step size
    }
    //Save structure and energy
    Eold = E;
    OldQMMMData = QMMMData;
    //Check optimization step size
    VecMax = abs(Forces.maxCoeff());
    if (abs(Forces.minCoeff()) > VecMax)
    {
      VecMax = StepScale*abs(Forces.minCoeff());
    }
    else
    {
      VecMax *= StepScale;
    }
    stepsize = StepScale;
    if (VecMax > QMMMOpts.maxStep)
    {
      //Scale step size
      cout << "    Scaling step size to match the maximum...";
      cout << '\n';
      stepsize *= (QMMMOpts.maxStep/VecMax);
    }
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (QMMMData[i].QMregion or QMMMData[i].PBregion)
      {
        QMMMData[i].P[bead].x += stepsize*Forces(ct);
        QMMMData[i].P[bead].y += stepsize*Forces(ct+1);
        QMMMData[i].P[bead].z += stepsize*Forces(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(QMMMData,qmfile,QMMMOpts);
    //Check convergence
    OptDone = OptConverged(QMMMData,OldQMMMData,Forces,stepct,QMMMOpts,bead,1);
    stepct += 1;
    //Increase step size
    StepScale *= 1.05;
    if (StepScale > QMMMOpts.stepScale)
    {
      //Prevent step size from getting too large
      StepScale = QMMMOpts.stepScale;
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

void LICHEMQuickMin(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
{
  //Cartesian damped Verlet optimizer (aka QuickMin)
  stringstream call; //Stream for system calls and reading/writing files
  int stepct = 0; //Counter for optimization steps
  fstream qmfile, inFile, outFile; //Generic file names
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize charges
  if (Nmm > 0)
  {
    WriteChargeFile(QMMMData,QMMMOpts,bead);
  }
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize optimization variables
  double VecMax = 0;
  bool OptDone = 0;
  vector<QMMMAtom> OldQMMMData = QMMMData; //Previous structure
  VectorXd QMVel(Ndof); //Velocity vector
  QMVel.setZero(); //Start at zero Kelvin
  //Run optimization
  double sdscale = 0.01; //Scale factor for SD steps
  double TimeStep; //Timestep for the Verlet algorithm
  TimeStep = sdscale*QMMMOpts.stepScale; //Make a local copy
  while ((!OptDone) and (stepct < QMMMOpts.maxOptSteps))
  {
    double E = 0;
    OldQMMMData = QMMMData; //Save old structure
    //Create blank force array
    VectorXd Forces(Ndof);
    Forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      E += GaussianForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      E += PSI4Forces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      E += NWChemForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,Forces,QMMMOpts,bead);
      if (AMOEBA or QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,Forces,QMMMOpts,bead);
      }
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      E += AMBERForces(QMMMData,Forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      E += LAMMPSForces(QMMMData,Forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Project velocities
    double VdotF = QMVel.dot(Forces); //Overlap of forces and velocities
    bool NoVelScale = 0; //Do not increase the timestep
    if (VdotF <= 0)
    {
      //Delete velocities and take a steepest descent step
      cout << "    Taking a SD step...";
      TimeStep = sdscale*QMMMOpts.stepScale; //Reset timestep
      QMVel = TimeStep*Forces; //SD step
      NoVelScale = 1; //Skip velocity scaling
    }
    else
    {
      //Damp velocities
      cout << "    Taking a DV step...";
      QMVel = VdotF*Forces; //Scale forces based on curvature
      TimeStep *= 1.30; //Increase timestep
      if (TimeStep > QMMMOpts.stepScale)
      {
        //Set to the maximum value
        TimeStep = QMMMOpts.stepScale;
      }
    }
    //Check optimization step size
    VecMax = TimeStep*QMVel.norm(); //Displacement after the update
    if (VecMax > QMMMOpts.maxStep)
    {
      //Take a smaller step
      cout << " Reducing timestep and velocites...";
      VecMax = (QMMMOpts.maxStep/VecMax); //Save to scale forces
      QMVel *= VecMax; //Reduce velocities
      TimeStep *= VecMax; //Reduce timestep
      if (TimeStep < (sdscale*QMMMOpts.stepScale))
      {
        //Revert to the minimum value
        TimeStep = sdscale*QMMMOpts.stepScale;
      }
    }
    else if ((TimeStep < QMMMOpts.stepScale) and (!NoVelScale))
    {
      //Finish update method output
      cout << " Increasing timestep...";
    }
    cout << '\n'; //Print a black line after printing the update method
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (QMMMData[i].QMregion or QMMMData[i].PBregion)
      {
        QMMMData[i].P[bead].x += TimeStep*QMVel(ct);
        QMMMData[i].P[bead].y += TimeStep*QMVel(ct+1);
        QMMMData[i].P[bead].z += TimeStep*QMVel(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(QMMMData,qmfile,QMMMOpts);
    //Check convergence
    OptDone = OptConverged(QMMMData,OldQMMMData,Forces,stepct,QMMMOpts,bead,1);
    stepct += 1;
    //Update velocities
    QMVel += TimeStep*Forces;
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
  int stepct = 0; //Counter for optimization steps
  fstream qmfile,inFile,outFile; //Generic file streams
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize charges
  if (Nmm > 0)
  {
    WriteChargeFile(QMMMData,QMMMOpts,bead);
  }
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Create DFP arrays
  VectorXd OptVec(Ndof); //Gradient descent direction
  VectorXd GradDiff(Ndof); //Change in the gradient
  VectorXd Forces(Ndof); //Forces
  MatrixXd IHess(Ndof,Ndof); //Inverse Hessian
  //Initialize arrays
  OptVec.setZero();
  GradDiff.setZero();
  Forces.setZero();
  //Create an identity matrix as the initial Hessian
  IHess.setIdentity(); //Already an "inverse" Hessian
  //Initialize optimization variables
  double StepScale; //Local copy
  double E = 0; //Energy
  double Eold = 0; //Energy from previous step
  double sdscale = 0.01; //Scale factor for SD steps
  double VecMax = 0; //Maxium atomic displacement
  bool OptDone = 0; //Flag to end the optimization
  //Calculate forces (QM part)
  if (Gaussian)
  {
    int tStart = (unsigned)time(0);
    E += GaussianForces(QMMMData,Forces,QMMMOpts,bead);
    QMTime += (unsigned)time(0)-tStart;
  }
  if (PSI4)
  {
    int tStart = (unsigned)time(0);
    E += PSI4Forces(QMMMData,Forces,QMMMOpts,bead);
    QMTime += (unsigned)time(0)-tStart;
    //Delete annoying useless files
    globalSys = system("rm -f psi.* timer.*");
  }
  if (NWChem)
  {
    int tStart = (unsigned)time(0);
    E += NWChemForces(QMMMData,Forces,QMMMOpts,bead);
    QMTime += (unsigned)time(0)-tStart;
  }
  //Calculate forces (MM part)
  if (TINKER)
  {
    int tStart = (unsigned)time(0);
    E += TINKERForces(QMMMData,Forces,QMMMOpts,bead);
    if (AMOEBA or QMMMOpts.useImpSolv)
    {
      //Forces from MM polarization
      E += TINKERPolForces(QMMMData,Forces,QMMMOpts,bead);
    }
    MMTime += (unsigned)time(0)-tStart;
  }
  if (AMBER)
  {
    int tStart = (unsigned)time(0);
    E += AMBERForces(QMMMData,Forces,QMMMOpts,bead);
    MMTime += (unsigned)time(0)-tStart;
  }
  if (LAMMPS)
  {
    int tStart = (unsigned)time(0);
    E += LAMMPSForces(QMMMData,Forces,QMMMOpts,bead);
    MMTime += (unsigned)time(0)-tStart;
  }
  //Output initial RMS force
  VecMax = 0; //Using this variable to avoid creating a new one
  VecMax = Forces.squaredNorm(); //Calculate initial RMS force
  VecMax = sqrt(VecMax/Ndof);
  cout << "    Performing a steepest descent step..." << '\n';
  cout << "    QM step: 0";
  cout << " | RMS force: " << LICHEMFormFloat(VecMax,12);
  cout << " eV/\u212B";
  cout << '\n' << '\n';
  cout.flush();
  //Optimize structure
  Eold = E; //Save energy
  StepScale = QMMMOpts.stepScale;
  StepScale *= sdscale; //Take a very small first step
  while ((!OptDone) and (stepct < QMMMOpts.maxOptSteps))
  {
    E = 0; // Reinitialize energy
    //Copy old structure and old forces
    vector<QMMMAtom> OldQMMMData = QMMMData;
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<Ndof;i++)
    {
      //Reinitialize the change in the gradient
      GradDiff(i) = Forces(i);
    }
    //Determine new structure
    OptVec = IHess*Forces;
    OptVec *= StepScale;
    //Check step size
    VecMax = OptVec.norm();
    if (VecMax > QMMMOpts.maxStep)
    {
      //Scale step size
      OptVec *= (QMMMOpts.maxStep/VecMax);
    }
    //Update positions
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (QMMMData[i].QMregion or QMMMData[i].PBregion)
      {
        QMMMData[i].P[bead].x += OptVec(ct);
        QMMMData[i].P[bead].y += OptVec(ct+1);
        QMMMData[i].P[bead].z += OptVec(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(QMMMData,qmfile,QMMMOpts);
    //Calculate forces (QM part)
    Forces.setZero();
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      E += GaussianForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      E += PSI4Forces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      E += NWChemForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,Forces,QMMMOpts,bead);
      if (AMOEBA or QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,Forces,QMMMOpts,bead);
      }
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      E += AMBERForces(QMMMData,Forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      E += LAMMPSForces(QMMMData,Forces,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Check stability
    double VecDotForces; //Dot product of the forces and optimization vector
    VecDotForces = OptVec.dot(Forces);
    double NormForce; //Local norm of the forces
    NormForce = Forces.norm(); //Take the norm of the foces
    NormForce /= sqrt(Ndof); //Make the norm an RMS value
    double LocalMaxForce; //Maximum value in the force array
    LocalMaxForce = abs(Forces.maxCoeff());
    if (((VecDotForces < 0) or (LocalMaxForce >= 1.0)) and
       (NormForce > (10*QMMMOpts.QMOptTol)))
    {
      //Optimizer is going the wrong direction and is not converged
      Eold = -1*hugeNum; //Force the Hessian to be rebuilt
    }
    //Update Hessian
    GradDiff -= Forces;
    if (((stepct%30) == 0) or (stepct < 15))
    {
      //Build a new Hessian after 30 steps
      cout << "    Performing a steepest descent step...";
      cout << '\n';
      //Shrink step size
      if ((stepct < 15) and (E < Eold))
      {
        //Reduce step size
        StepScale = sdscale*QMMMOpts.stepScale; //Small step
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      //Check for minimum step size
      if (StepScale < (0.25*sdscale*QMMMOpts.stepScale))
      {
        StepScale = 0.25*sdscale*QMMMOpts.stepScale;
      }
      //Create new Hessian as an identity matrix
      IHess.setIdentity(); //Already an "inverse" Hessian
    }
    else if (((stepct+1)%30) == 0)
    {
      //Prepare for the upcoming SD step
      cout << "    Reducing the step size...";
      cout << '\n';
      //Shrink step size
      if (StepScale > (sdscale*QMMMOpts.stepScale))
      {
        //Reduce step size
        StepScale = sdscale*QMMMOpts.stepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      //Check for minimum step size
      if (StepScale < (0.25*sdscale*QMMMOpts.stepScale))
      {
        StepScale = 0.25*sdscale*QMMMOpts.stepScale;
      }
      //Create new Hessian as an identity matrix
      IHess.setIdentity(); //Already an "inverse" Hessian
    }
    else if (E < Eold)
    {
      //Update Hessian
      cout << "    Updating inverse Hessian...";
      cout << '\n';
      //Start really long "line"
      IHess = IHess+((OptVec*OptVec.transpose())/(OptVec.transpose()
      *GradDiff))-((IHess*GradDiff*GradDiff.transpose()*IHess)
      /(GradDiff.transpose()*IHess*GradDiff));
      //End really long "line"
      //Increase stepsize for the next iteration
      StepScale *= 1.20;
      if (StepScale > QMMMOpts.stepScale)
      {
        //Prevent step size from getting too large
        StepScale = QMMMOpts.stepScale;
      }
    }
    else
    {
      //Take a small steepest descent step and rebuild Hessian
      cout << "    Potentially unstable Hessian.";
      cout << " Constructing new Hessian...";
      cout << '\n';
      //Reduce step size
      if (StepScale > (sdscale*QMMMOpts.stepScale))
      {
        StepScale = sdscale*QMMMOpts.stepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      //Check for minimum step size
      if (StepScale < (0.25*sdscale*QMMMOpts.stepScale))
      {
        StepScale = 0.25*sdscale*QMMMOpts.stepScale;
      }
      IHess.setIdentity(); //Already an "inverse" Hessian
    }
    //Save energy
    Eold = E;
    //Check convergence
    stepct += 1;
    OptDone = OptConverged(QMMMData,OldQMMMData,Forces,stepct,QMMMOpts,bead,1);
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << bead << ".xyz";
  call << " MMCharges_" << bead << ".txt";
  globalSys = system(call.str().c_str());
  //Finish and return
  return;
};

//Ensemble optimizers
void EnsembleSD(vector<QMMMAtom>& QMMMData, fstream& traj,
                QMMMSettings& QMMMOpts, int bead)
{
  //Ensemble steepest descent optimizer
  stringstream call; //Stream for system calls and reading/writing files
  int stepct = 0; //Counter for optimization steps
  fstream inFile, outFile; //Generic file names
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize optimization variables
  double stepsize = 1;
  double StepScale = QMMMOpts.stepScale; //Saved copy
  //Run optimization
  while (stepct < QMMMOpts.maxOptSteps)
  {
    //Run MD
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      TINKERDynamics(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
      if (AMOEBA)
      {
        //Set up current multipoles
        RotateTINKCharges(QMMMData,bead);
      }
    }
    //Perform SD step
    double E = 0;
    double SumE = 0;
    //Create blank force array
    VectorXd Forces(Ndof);
    Forces.setZero();
    //Calculate forces and energy (QM part)
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      SumE += GaussianForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      SumE += PSI4Forces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      SumE += NWChemForces(QMMMData,Forces,QMMMOpts,bead);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate forces and energy (MM part)
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      E += TINKERForces(QMMMData,Forces,QMMMOpts,bead);
      SumE += TINKEREnergy(QMMMData,QMMMOpts,bead);
      if (AMOEBA or QMMMOpts.useImpSolv)
      {
        //Forces from MM polarization
        E += TINKERPolForces(QMMMData,Forces,QMMMOpts,bead);
      }
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      E += AMBERForces(QMMMData,Forces,QMMMOpts,bead);
      SumE += AMBEREnergy(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      E += LAMMPSForces(QMMMData,Forces,QMMMOpts,bead);
      SumE += LAMMPSEnergy(QMMMData,QMMMOpts,bead);
      MMTime += (unsigned)time(0)-tStart;
    }
    //Determine new structure
    int ct = 0; //Counter for QM and PB atoms
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (QMMMData[i].QMregion or QMMMData[i].PBregion)
      {
        //Check X step size
        stepsize = StepScale*Forces(ct);
        if (abs(stepsize) > QMMMOpts.maxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.maxStep/abs(stepsize);
        }
        QMMMData[i].P[bead].x += stepsize;
        //Check Y step size
        stepsize = StepScale*Forces(ct+1);
        if (abs(stepsize) > QMMMOpts.maxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.maxStep/abs(stepsize);
        }
        QMMMData[i].P[bead].y += stepsize;
        //Check Z step size
        stepsize = StepScale*Forces(ct+2);
        if (abs(stepsize) > QMMMOpts.maxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.maxStep/abs(stepsize);
        }
        QMMMData[i].P[bead].z += stepsize;
        ct += 3;
      }
    }
    //Print structure and energy
    stepct += 1;
    Print_traj(QMMMData,traj,QMMMOpts);
    cout << " | Step: " << stepct;
    cout << " | Simulation time: ";
    cout << (stepct*QMMMOpts.dt*QMMMOpts.NSteps/1000);
    cout << " ps | Energy: " << LICHEMFormFloat(SumE,16);
    cout << " eV" << '\n';
    cout.flush();
  }
  //Clean up files and return
  call.str("");
  call << "rm -f LICHM_" << bead << ".dyn";
  globalSys = system(call.str().c_str());
  return;
};

