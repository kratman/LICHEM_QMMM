/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 A set of optimization routines for the QM part of QMMM calculculations. These
 are inefficient, however, they are useful for testing and debugging.

 Reference for optimization routines:
 Press et al., Numerical Recipes 2nd Edition, 1997

*/

//Convergence test functions
bool OptConverged(vector<QMMMAtom>& Struct, vector<QMMMAtom>& OldStruct,
     vector<Coord>& Forces, int stepct, QMMMSettings& QMMMOpts, int Bead,
     bool QMregion)
{
  //Check convergence of QMMM optimizations
  stringstream call;
  string dummy; //Generic string
  call.str("");
  //Initialize stats variables
  bool OptDone = 0;
  double RMSdiff = 0;
  double RMSforce = 0;
  double MAXforce = 0;
  double SumE = 0;
  //Check progress
  if (QMregion)
  {
    //Check if a QM calculation is converged
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Calculate RMS forces
      double Fx = Forces[i].x;
      double Fy = Forces[i].y;
      double Fz = Forces[i].z;
      if (abs(Fx) > MAXforce)
      {
        MAXforce = abs(Fx);
      }
      if (abs(Fy) > MAXforce)
      {
        MAXforce = abs(Fy);
      }
      if (abs(Fz) > MAXforce)
      {
        MAXforce = abs(Fz);
      }
      RMSforce += Fx*Fx+Fy*Fy+Fz*Fz;
    }
    #pragma omp parallel for reduction(+:RMSdiff)
    for (int i=0;i<Natoms;i++)
    {
      //Calculate RMS displacement
      if (Struct[i].QMregion or Struct[i].PAregion)
      {
        for (int j=0;j<i;j++)
        {
          if (Struct[j].QMregion or Struct[j].PAregion)
          {
            double Rnew = 0;
            double Rold = 0;
            Rnew = CoordDist2(Struct[i].P[Bead],Struct[j].P[Bead]);
            Rold = CoordDist2(OldStruct[i].P[Bead],OldStruct[j].P[Bead]);
            Rnew = sqrt(Rnew);
            Rold = sqrt(Rold);
            RMSdiff += (Rnew-Rold)*(Rnew-Rold);
          }
        }
      }
    }
    #pragma omp barrier
    RMSdiff /= (Nqm+Npseudo)*(Nqm+Npseudo-1)/2;
    RMSdiff = sqrt(RMSdiff);
    RMSforce /= 3*(Nqm+Npseudo);
    RMSforce = sqrt(RMSforce);
    //Print progress
    call.copyfmt(cout); //Save settings
    cout << setprecision(12);
    cout << "    QM step: " << stepct;
    cout << " | RMS dev: " << RMSdiff;
    cout << " \u212B" << '\n';
    cout << "    Max. force: " << MAXforce;
    cout << " eV/\u212B | RMS force: " << RMSforce;
    cout << " eV/\u212B" << '\n';
    cout.copyfmt(call); //Return to previous settings
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
    //Check if the MM region changed and gather statistics
    SumE = 0;
    //Calculate QM energy
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += PSIEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Calculate RMS displacement
    #pragma omp parallel for reduction(+:RMSdiff)
    for (int i=0;i<Natoms;i++)
    {
      for (int j=0;j<i;j++)
      {
        double Rnew = 0;
        double Rold = 0;
        Rnew = CoordDist2(Struct[i].P[Bead],Struct[j].P[Bead]);
        Rold = CoordDist2(OldStruct[i].P[Bead],OldStruct[j].P[Bead]);
        Rnew = sqrt(Rnew);
        Rold = sqrt(Rold);
        RMSdiff += (Rnew-Rold)*(Rnew-Rold);
      }
    }
    #pragma omp barrier
    RMSdiff /= (Natoms-Nfreeze)*(Natoms-Nfreeze-1)/2;
    RMSdiff = sqrt(RMSdiff);
    //Print progress
    call.copyfmt(cout); //Save settings
    cout << setprecision(12);
    cout << " | Opt. step: ";
    cout << stepct << " | Energy: ";
    cout << SumE << " eV ";
    cout << " | RMS dev: " << RMSdiff;
    cout << " \u212B" << '\n';
    cout.copyfmt(call); //Replace settings
    //Check convergence
    if (RMSdiff <= QMMMOpts.MMOptTol)
    {
      OptDone = 1;
      if (QMMM)
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
void FLUKESteepest(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Cartesian steepest descent optimizer
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  int stepct = 0; //Counter for optimization steps
  fstream qmfile, ifile, ofile; //Generic file names
  double Eold = 0; //Old saved energy
  //Initialize files
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize charges for Gaussian
  if ((AMOEBA == 1) and ((Gaussian == 1) or (NWChem == 1)))
  {
    if (TINKER == 1)
    {
      //Set up current multipoles
      RotateTINKCharges(Struct,Bead);
    }
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile << '\n';
      }
    }
    ofile.copyfmt(cout);
    ofile.flush();
    ofile.close();
  }
  //Initialize optimization variables
  double stepsize = 1;
  double VecMax = 0;
  bool OptDone = 0;
  vector<QMMMAtom> OldStruct = Struct; //Previous structure
  //Run optimization
  double StepScale = QMMMOpts.StepScale;
  StepScale *= 0.75; //Take a smaller first step
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    double E = 0;
    //Create blank force array
    vector<Coord> Forces;
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Create arrays with zeros
      Coord tmp;
      tmp.x = 0;
      tmp.y = 0;
      tmp.z = 0;
      Forces.push_back(tmp);
    }
    //Calculate forces (QM part)
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      E += PSIForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA == 1)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER == 1)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS == 1)
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
    VecMax = 0;
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Check if the step size is too large
      if (abs(StepScale*Forces[i].x) > VecMax)
      {
        VecMax = abs(StepScale*Forces[i].x);
      }
      if (abs(StepScale*Forces[i].y) > VecMax)
      {
        VecMax = abs(StepScale*Forces[i].y);
      }
      if (abs(StepScale*Forces[i].z) > VecMax)
      {
        VecMax = abs(StepScale*Forces[i].z);
      }
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
      if (Struct[i].QMregion or Struct[i].PAregion)
      {
        Struct[i].P[Bead].x += stepsize*Forces[ct].x;
        Struct[i].P[Bead].y += stepsize*Forces[ct].y;
        Struct[i].P[Bead].z += stepsize*Forces[ct].z;
        ct += 1;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Check convergence
    stepct += 1;
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
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

void FLUKEDFP(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //A simple DFP optimizer, which is similar to BFGS updating
  //Note: This optimizer does not have a true line search, instead
  //a steepest descent step is performed if the energy rises
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  int stepct = 0; //Counter for optimization steps
  fstream qmfile,ifile,ofile; //Generic file names
  //Initialize files
  call.str("");
  call << "QMOpt_" << Bead << ".xyz";
  qmfile.open(call.str().c_str(),ios_base::out);
  //Initialize variables
  double E = 0; //Energy
  double Eold = 0; //Energy from previous step
  double VecMax = 0;
  //Initialize multipoles for Gaussian optimizations
  if ((AMOEBA == 1) and ((Gaussian == 1) or (NWChem == 1)))
  {
    if (TINKER == 1)
    {
      //Set up current multipoles
      RotateTINKCharges(Struct,Bead);
    }
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile << '\n';
      }
    }
    ofile.copyfmt(cout);
    ofile.flush();
    ofile.close();
  }
  //Create DFP arrays
  VectorXd OptVec(3*(Nqm+Npseudo)); //Gradient descent direction
  VectorXd GradDiff(3*(Nqm+Npseudo)); //Change in the gradient
  VectorXd NGrad(3*(Nqm+Npseudo)); //Negative of the gradient
  MatrixXd IHess(3*(Nqm+Npseudo),3*(Nqm+Npseudo)); //Inverse Hessian
  #pragma omp parallel for
  for (int i=0;i<(3*(Nqm+Npseudo));i++)
  {
    //Initialize arrays
    OptVec(i) = 0;
    GradDiff(i) = 0;
    NGrad(i) = 0;
    //Create an identity matrix as the initial Hessian
    for (int j=0;j<i;j++)
    {
      //Set off diagonal terms
      IHess(i,j) = 0.0;
      IHess(j,i) = 0.0;
    }
    IHess(i,i) = 1.0; //Already an "inverse Hessian"
  }
  #pragma omp barrier
  vector<Coord> Forces;
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Create arrays with zeros
    Coord tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    Forces.push_back(tmp);
  }
  //Calculate forces (QM part)
  if (Gaussian == 1)
  {
    int tstart = (unsigned)time(0);
    E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
  }
  if (PSI4 == 1)
  {
    int tstart = (unsigned)time(0);
    E += PSIForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
    //Clean up annoying useless files
    GlobalSys = system("rm -f psi.*");
  }
  if (NWChem == 1)
  {
    int tstart = (unsigned)time(0);
    E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
    QMTime += (unsigned)time(0)-tstart;
  }
  //Calculate forces (MM part)
  if (TINKER == 1)
  {
    int tstart = (unsigned)time(0);
    E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
    if (AMOEBA == 1)
    {
      E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
    }
    MMTime += (unsigned)time(0)-tstart;
  }
  if (AMBER == 1)
  {
    int tstart = (unsigned)time(0);
    E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
    MMTime += (unsigned)time(0)-tstart;
  }
  if (LAMMPS == 1)
  {
    int tstart = (unsigned)time(0);
    E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
    MMTime += (unsigned)time(0)-tstart;
  }
  //Save forces
  int ct = 0; //Counter
  VecMax = 0; //Using this variable to avoid creating a new one
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Change forces array
    NGrad(ct) = Forces[i].x;
    NGrad(ct+1) = Forces[i].y;
    NGrad(ct+2) = Forces[i].z;
    ct += 3;
    //Calculate initial RMS force
    VecMax += Forces[i].x*Forces[i].x;
    VecMax += Forces[i].y*Forces[i].y;
    VecMax += Forces[i].z*Forces[i].z;
  }
  //Output initial RMS force
  VecMax = sqrt(VecMax/(3*(Nqm+Npseudo)));
  call.copyfmt(cout); //Save settings
  cout << setprecision(12);
  cout << "    QM step: 0";
  cout << " | RMS force: " << VecMax;
  cout << " eV/\u212B";
  cout << '\n' << '\n';
  cout.flush();
  cout.copyfmt(call); //Return to previous settings
  //Optimize structure
  Eold = E;
  bool OptDone = 0;
  double StepScale = QMMMOpts.StepScale;
  StepScale *= 0.05; //Take a very small first step
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    E = 0; // Reinitialize energy
    //Copy old structure and delete old forces force array
    vector<QMMMAtom> OldStruct = Struct;
    ct = 0; //Counter
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Reinitialize the change in the gradient
      GradDiff(ct) = NGrad(ct);
      GradDiff(ct+1) = NGrad(ct+1);
      GradDiff(ct+2) = NGrad(ct+2);
      ct += 3;
      //Delete old forces
      Forces[i].x = 0;
      Forces[i].y = 0;
      Forces[i].z = 0;
    }
    //Determine new structure
    OptVec = IHess*NGrad;
    OptVec *= StepScale;
    VecMax = 0;
    for (int i=0;i<(3*(Nqm+Npseudo));i++)
    {
      //Check if the step size is too large
      if (abs(OptVec(i)) > VecMax)
      {
        VecMax = abs(OptVec(i));
      }
    }
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Scale step size
      OptVec *= (QMMMOpts.MaxStep/VecMax);
    }
    ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if (Struct[i].QMregion or Struct[i].PAregion)
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
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      E += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      E += PSIForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      E += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces (MM part)
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      if (AMOEBA == 1)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER == 1)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Update Hessian
    ct = 0;
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      GradDiff(ct) -= Forces[i].x;
      GradDiff(ct+1) -= Forces[i].y;
      GradDiff(ct+2) -= Forces[i].z;
      NGrad(ct) = Forces[i].x;
      NGrad(ct+1) = Forces[i].y;
      NGrad(ct+2) = Forces[i].z;
      ct += 3;
    }
    if (((stepct%30) == 0) and (stepct != 0))
    {
      //Build a new Hessian after 30 steps
      cout << "    Constructing new Hessian...";
      cout << '\n';
      //Shrink step size
      if (StepScale > (0.02*QMMMOpts.StepScale))
      {
        StepScale = 0.02*QMMMOpts.StepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      #pragma omp parallel for
      for (int i=0;i<(3*(Nqm+Npseudo));i++)
      {
        //Create identity matrix
        for (int j=0;j<i;j++)
        {
          IHess(i,j) = 0.0;
          IHess(j,i) = 0.0;
        }
        IHess(i,i) = 1.0; //Already an "inverse Hessian"
      }
      #pragma omp barrier
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
      //Increase stepsize
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
      cout << "    Energy did not decrease. Constructing new Hessian...";
      cout << '\n';
      //Shrink step size
      if (StepScale > (0.02*QMMMOpts.StepScale))
      {
        StepScale = 0.02*QMMMOpts.StepScale;
      }
      else
      {
        //Reduce step size further
        StepScale *= 0.75;
      }
      #pragma omp parallel for
      for (int i=0;i<(3*(Nqm+Npseudo));i++)
      {
        //Create identity matrix
        for (int j=0;j<i;j++)
        {
          IHess(i,j) = 0.0;
          IHess(j,i) = 0.0;
        }
        IHess(i,i) = 1.0; //Already an "inverse Hessian"
      }
      #pragma omp barrier
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
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  int stepct = 0; //Counter for optimization steps
  fstream ifile, ofile; //Generic file names
  //Initialize optimization variables
  double stepsize = 1;
  double VecMax = 0;
  //Run optimization
  double StepScale = QMMMOpts.StepScale;
  while (stepct < QMMMOpts.MaxOptSteps)
  {
    //Run MD
    if (TINKER == 1)
    {
      TINKERDynamics(Struct,QMMMOpts,Bead);
      if (AMOEBA == 1)
      {
        //Set up current multipoles
        RotateTINKCharges(Struct,Bead);
      }
    }
    //Perform SD step
    double E = 0;
    double SumE = 0;
    //Create blank force array
    vector<Coord> Forces;
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Create arrays with zeros
      Coord tmp;
      tmp.x = 0;
      tmp.y = 0;
      tmp.z = 0;
      Forces.push_back(tmp);
    }
    //Calculate forces and energy (QM part)
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += PSIForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate forces and energy (MM part)
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      SumE += TINKEREnergy(Struct,QMMMOpts,Bead);
      if (AMOEBA == 1)
      {
        E += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
      }
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER == 1)
    {
      int tstart = (unsigned)time(0);
      E += AMBERForces(Struct,Forces,QMMMOpts,Bead);
      SumE += AMBEREnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      E += LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Check optimization step size
    VecMax = 0;
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Check if the step size is too large
      if (abs(StepScale*Forces[i].x) > VecMax)
      {
        VecMax = abs(StepScale*Forces[i].x);
      }
      if (abs(StepScale*Forces[i].y) > VecMax)
      {
        VecMax = abs(StepScale*Forces[i].y);
      }
      if (abs(StepScale*Forces[i].z) > VecMax)
      {
        VecMax = abs(StepScale*Forces[i].z);
      }
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
      if (Struct[i].QMregion or Struct[i].PAregion)
      {
        Struct[i].P[Bead].x += stepsize*Forces[ct].x;
        Struct[i].P[Bead].y += stepsize*Forces[ct].y;
        Struct[i].P[Bead].z += stepsize*Forces[ct].z;
        ct += 1;
      }
    }
    //Print structure and energy
    stepct += 1;
    Print_traj(Struct,traj,QMMMOpts);
    call.copyfmt(cout); //Save settings
    cout << setprecision(16);
    cout << " | Opt. step: ";
    cout << stepct << " | Energy: ";
    cout << SumE << " eV" << '\n';
    cout.copyfmt(call); //Replace settings
    cout.flush();
  }
  //Clean up files and return
  call.str("");
  call << "rm -f QMMM_" << Bead << ".dyn";
  GlobalSys = system(call.str().c_str());
  return;
};

