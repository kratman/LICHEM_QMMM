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
bool OptConverged(vector<QMMMAtom>& Struct, vector<QMMMAtom>& OldStruct,
                  VectorXd& Forces, int stepct, QMMMSettings& QMMMOpts,
                  int Bead, bool QMregion)
{
  //Check convergence of QMMM optimizations
  bool OptDone = 0; //Ends the simulation
  double RMSdiff = 0; //RMS deviation
  double RMSforce = 0; //RMS force
  double MAXforce = 0; //Maximum force
  double SumE = 0; //Storage for energies
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
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
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        for (int j=0;j<i;j++)
        {
          if (Struct[j].QMregion or Struct[j].PBregion)
          {
            double Rnew = 0;
            double Rold = 0;
            Rnew = CoordDist2(Struct[i].P[Bead],Struct[j].P[Bead]).VecMag();
            Rold = CoordDist2(OldStruct[i].P[Bead],
                              OldStruct[j].P[Bead]).VecMag();
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
    if ((RMSdiff <= QMMMOpts.QMOptTol) and
       (RMSforce <= (10*QMMMOpts.QMOptTol)) and
       (MAXforce <= (20*QMMMOpts.QMOptTol)))
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
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSI4Energy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
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
        Rnew = CoordDist2(Struct[i].P[Bead],Struct[j].P[Bead]).VecMag();
        Rold = CoordDist2(OldStruct[i].P[Bead],OldStruct[j].P[Bead]).VecMag();
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
void LICHEMSteepest(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                    int Bead)
{
  //Cartesian steepest descent optimizer
  stringstream call; //Stream for system calls and reading/writing files
  int stepct = 0; //Counter for optimization steps
  fstream qmfile, ifile, ofile; //Generic file names
  double Eold = 0; //Old saved energy
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize charges
  WriteChargeFile(Struct,QMMMOpts,Bead);
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize optimization variables
  double stepsize = 1;
  double VecMax = 0;
  bool OptDone = 0;
  vector<QMMMAtom> OldStruct = Struct; //Previous structure
  //Run optimization
  double StepScale = QMMMOpts.StepScale;
  StepScale *= 0.70; //Take a smaller first step
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    double E = 0;
    //Create blank force array
    VectorXd Forces(Ndof);
    Forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      E += PSI4Forces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
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
    OldStruct = Struct;
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
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Scale step size
      cout << "    Scaling step size to match the maximum...";
      cout << '\n';
      stepsize *= (QMMMOpts.MaxStep/VecMax);
    }
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        Struct[i].P[Bead].x += stepsize*Forces(ct);
        Struct[i].P[Bead].y += stepsize*Forces(ct+1);
        Struct[i].P[Bead].z += stepsize*Forces(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Check convergence
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
    stepct += 1;
    //Increase step size
    StepScale *= 1.05;
    if (StepScale > QMMMOpts.StepScale)
    {
      //Prevent step size from getting too large
      StepScale = QMMMOpts.StepScale;
    }
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Finish and return
  return;
};

void LICHEMQuickMin(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                    int Bead)
{
  //Cartesian damped Verlet optimizer (aka QuickMin)
  stringstream call; //Stream for system calls and reading/writing files
  int stepct = 0; //Counter for optimization steps
  fstream qmfile, ifile, ofile; //Generic file names
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize charges
  WriteChargeFile(Struct,QMMMOpts,Bead);
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize optimization variables
  double VecMax = 0;
  bool OptDone = 0;
  vector<QMMMAtom> OldStruct = Struct; //Previous structure
  VectorXd QMVel(Ndof); //Velocity vector
  QMVel.setZero(); //Start at zero Kelvin
  //Run optimization
  double sdscale = 0.01; //Scale factor for SD steps
  double TimeStep; //Timestep for the Verlet algorithm
  TimeStep = sdscale*QMMMOpts.StepScale; //Make a local copy
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    double E = 0;
    OldStruct = Struct; //Save old structure
    //Create blank force array
    VectorXd Forces(Ndof);
    Forces.setZero();
    //Calculate forces (QM part)
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      E += PSI4Forces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Project velocities
    double VdotF = QMVel.dot(Forces); //Overlap of forces and velocities
    bool NoVelScale = 0; //Do not increase the timestep
    if (VdotF <= 0)
    {
      //Delete velocities and take a steepest descent step
      cout << "    Taking a SD step...";
      TimeStep = sdscale*QMMMOpts.StepScale; //Reset timestep
      QMVel = TimeStep*Forces; //SD step
      NoVelScale = 1; //Skip velocity scaling
    }
    else
    {
      //Damp velocities
      cout << "    Taking a DV step...";
      QMVel = VdotF*Forces; //Scale forces based on curvature
      TimeStep *= 1.30; //Increase timestep
      if (TimeStep > QMMMOpts.StepScale)
      {
        //Set to the maximum value
        TimeStep = QMMMOpts.StepScale;
      }
    }
    //Check optimization step size
    VecMax = TimeStep*QMVel.norm(); //Displacement after the update
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Take a smaller step
      cout << " Reducing timestep and velocites...";
      VecMax = (QMMMOpts.MaxStep/VecMax); //Save to scale forces
      QMVel *= VecMax; //Reduce velocities
      TimeStep *= VecMax; //Reduce timestep
      if (TimeStep < (sdscale*QMMMOpts.StepScale))
      {
        //Revert to the minimum value
        TimeStep = sdscale*QMMMOpts.StepScale;
      }
    }
    else if ((TimeStep < QMMMOpts.StepScale) and (!NoVelScale))
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
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        Struct[i].P[Bead].x += TimeStep*QMVel(ct);
        Struct[i].P[Bead].y += TimeStep*QMVel(ct+1);
        Struct[i].P[Bead].z += TimeStep*QMVel(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Check convergence
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
    stepct += 1;
    //Update velocities
    QMVel += TimeStep*Forces;
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Finish and return
  return;
};

void LICHEMDFP(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //A simple Davidon-Fletcher-Powell optimizer
  //NB: This optimizer does not have a true line search, instead
  //a steepest descent step is performed if the optimizer is unstable
  stringstream call; //Stream for system calls and reading/writing files
  int stepct = 0; //Counter for optimization steps
  fstream qmfile,ifile,ofile; //Generic file streams
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize charges
  WriteChargeFile(Struct,QMMMOpts,Bead);
  //Initialize QM trajectory file
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
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
    int tstart = (unsigned)time(0);
    E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
  }
  if (PSI4)
  {
    int tstart = (unsigned)time(0);
    E += PSI4Forces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
    //Delete annoying useless files
    GlobalSys = system("rm -f psi.* timer.*");
  }
  if (NWChem)
  {
    int tstart = (unsigned)time(0);
    E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
  }
  //Calculate forces (MM part)
  if (TINKER)
  {
    int tstart = (unsigned)time(0);
    E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
    if (AMOEBA)
    {
      E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
    }
    MMTime += (unsigned)time(0)-tstart;
  }
  if (AMBER)
  {
    int tstart = (unsigned)time(0);
    E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
    MMTime += (unsigned)time(0)-tstart;
  }
  if (LAMMPS)
  {
    int tstart = (unsigned)time(0);
    E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
    MMTime += (unsigned)time(0)-tstart;
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
  StepScale = QMMMOpts.StepScale;
  StepScale *= sdscale; //Take a very small first step
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    E = 0; // Reinitialize energy
    //Copy old structure and old forces
    vector<QMMMAtom> OldStruct = Struct;
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
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Scale step size
      OptVec *= (QMMMOpts.MaxStep/VecMax);
    }
    //Update positions
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        Struct[i].P[Bead].x += OptVec(ct);
        Struct[i].P[Bead].y += OptVec(ct+1);
        Struct[i].P[Bead].z += OptVec(ct+2);
        ct += 3;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Calculate forces (QM part)
    Forces.setZero();
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      E += PSI4Forces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Check stability
    double VecDotForces; //Dot product of the forces and optimization vector
    VecDotForces = OptVec.dot(Forces);
    double NormForce; //Local norm of the forces
    NormForce = Forces.norm(); //Take the norm of the foces
    NormForce /= sqrt(Ndof); //Make the norm an RMS value
    if ((VecDotForces < 0) and (NormForce > (10*QMMMOpts.QMOptTol)))
    {
      //Optimizer is going the wrong direction and is not converged
      Eold = -1*HugeNum; //Force the Hessian to be rebuilt
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
        StepScale = sdscale*QMMMOpts.StepScale; //Small step
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      //Check for minimum step size
      if (StepScale < (0.25*sdscale*QMMMOpts.StepScale))
      {
        StepScale = 0.25*sdscale*QMMMOpts.StepScale;
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
      if (StepScale > (sdscale*QMMMOpts.StepScale))
      {
        //Reduce step size
        StepScale = sdscale*QMMMOpts.StepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      //Check for minimum step size
      if (StepScale < (0.25*sdscale*QMMMOpts.StepScale))
      {
        StepScale = 0.25*sdscale*QMMMOpts.StepScale;
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
      if (StepScale > QMMMOpts.StepScale)
      {
        //Prevent step size from getting too large
        StepScale = QMMMOpts.StepScale;
      }
    }
    else
    {
      //Take a small steepest descent step and rebuild Hessian
      cout << "    Potentially unstable structure.";
      cout << " Constructing new Hessian...";
      cout << '\n';
      //Reduce step size
      if (StepScale > (sdscale*QMMMOpts.StepScale))
      {
        StepScale = sdscale*QMMMOpts.StepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      //Check for minimum step size
      if (StepScale < (0.25*sdscale*QMMMOpts.StepScale))
      {
        StepScale = 0.25*sdscale*QMMMOpts.StepScale;
      }
      IHess.setIdentity(); //Already an "inverse" Hessian
    }
    //Save energy
    Eold = E;
    //Check convergence
    stepct += 1;
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Finish and return
  return;
};

//Ensemble optimizers
void EnsembleSD(vector<QMMMAtom>& Struct, fstream& traj,
                QMMMSettings& QMMMOpts, int Bead)
{
  //Ensemble steepest descent optimizer
  stringstream call; //Stream for system calls and reading/writing files
  int stepct = 0; //Counter for optimization steps
  fstream ifile, ofile; //Generic file names
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Initialize optimization variables
  double stepsize = 1;
  double StepScale = QMMMOpts.StepScale; //Saved copy
  //Run optimization
  while (stepct < QMMMOpts.MaxOptSteps)
  {
    //Run MD
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      TINKERDynamics(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
      if (AMOEBA)
      {
        //Set up current multipoles
        RotateTINKCharges(Struct,Bead);
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
      int tstart = (unsigned)time(0);
      SumE += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSI4Forces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces and energy (MM part)
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      SumE += TINKEREnergy(Struct,QMMMOpts,Bead);
      if (AMOEBA)
      {
        //Force from MM polarization
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      SumE += AMBEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Determine new structure
    int ct = 0; //Counter for QM and PB atoms
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Check X step size
        stepsize = StepScale*Forces(ct);
        if (abs(stepsize) > QMMMOpts.MaxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.MaxStep/abs(stepsize);
        }
        Struct[i].P[Bead].x += stepsize;
        //Check Y step size
        stepsize = StepScale*Forces(ct+1);
        if (abs(stepsize) > QMMMOpts.MaxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.MaxStep/abs(stepsize);
        }
        Struct[i].P[Bead].y += stepsize;
        //Check Z step size
        stepsize = StepScale*Forces(ct+2);
        if (abs(stepsize) > QMMMOpts.MaxStep)
        {
          //Scale step
          stepsize *= QMMMOpts.MaxStep/abs(stepsize);
        }
        Struct[i].P[Bead].z += stepsize;
        ct += 3;
      }
    }
    //Print structure and energy
    stepct += 1;
    Print_traj(Struct,traj,QMMMOpts);
    cout << " | Step: " << stepct;
    cout << " | Simulation time: ";
    cout << (stepct*QMMMOpts.dt*QMMMOpts.Nsteps/1000);
    cout << " ps | Energy: " << LICHEMFormFloat(SumE,16);
    cout << " eV" << '\n';
    cout.flush();
  }
  //Clean up files and return
  call.str("");
  call << "rm -f LICHM_" << Bead << ".dyn";
  GlobalSys = system(call.str().c_str());
  return;
};

