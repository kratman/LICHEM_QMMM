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
    for (int i=0;i<Natoms;i++)
    {
      //Calculate RMS displacement
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        for (int j=0;j<i;j++)
        {
          if ((Struct[j].QMregion == 1) or (Struct[j].PAregion == 1))
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
    RMSdiff /= (Nqm+Npseudo)*(Nqm+Npseudo-1)/2;
    RMSdiff = sqrt(RMSdiff);
    RMSforce /= 3*(Nqm+Npseudo);
    RMSforce = sqrt(RMSforce);
    //Print progress
    call.copyfmt(cout); //Save settings
    cout << setprecision(8);
    cout << "    QM Step: " << stepct;
    cout << " | RMS dev: " << RMSdiff;
    cout << " \u212B" << '\n';
    cout << "    Max force: " << MAXforce;
    cout << " eV/\u212B | RMS force: " << RMSforce;
    cout << " eV/\u212B";
    cout << '\n' << '\n';
    cout.flush();
    cout.copyfmt(call); //Return to previous settings
    //Check convergence criteria
    if ((RMSdiff <= QMMMOpts.QMOptTol) and
       (RMSforce <= (50*QMMMOpts.QMOptTol)) and
       (MAXforce <= (100*QMMMOpts.QMOptTol)))
    {
      OptDone = 1;
    }
  }
  if (!QMregion)
  {
    SumE = 0;
    //Check if the MM region changed and gather statistics
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += PSIEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      int sys = system("rm -f psi.*");
    }
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Calculate RMS displacement
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
    RMSdiff /= Natoms*(Natoms-1)/2;
    RMSdiff = sqrt(RMSdiff);
    //Print progress
    cout << " | Opt. Step: ";
    cout << stepct << " | Energy: ";
    cout << SumE << " eV ";
    cout << " | RMS dev: " << RMSdiff;
    cout << " \u212B" << '\n';
    //Check convergence
    if (RMSdiff <= QMMMOpts.MMOptTol)
    {
      OptDone = 1;
      cout << "    MM relaxation satisfactory.";
      cout << '\n';
    }
  }
  return OptDone;
};

//Optimizer functions
void FLUKESteepest(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Steepest descent optimizer
  int sys; //Dummy return for system calls
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
  if ((AMOEBA == 1) and (Gaussian == 1))
  {
    if (TINKER == 1)
    {
      //Set up current charges
      RotateTINKCharges(Struct,Bead);
    }
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
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
  //Optimize structure
  double stepsize = 1;
  double VecMax = 0;
  bool OptDone = 0;
  while ((!OptDone) and (stepct < QMMMOpts.MaxOptSteps))
  {
    double E = 0;
    //Copy old structure and create blank force array
    vector<QMMMAtom> OldStruct = Struct;
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
      int sys = system("rm -f psi.*");
    }
    //Calculate forces (MM part)
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Check step size
    if (E >= Eold)
    {
      //Take smaller steps if the energy does not improve
      cout << "    Energy did not decrease. Reducing step size by 15%...";
      cout << '\n';
      QMMMOpts.StepScale *= 0.85;
    }
    else
    {
      //Take larger steps if the energy is still decreasing
      cout << "    Energy is decreasing. Increasing step size by 1%...";
      cout << '\n';
      QMMMOpts.StepScale *= 1.01;
    }
    //Check optimization step size
    VecMax = 0;
    for (int i=0;i<(Nqm+Npseudo);i++)
    {
      //Check if the step size is too large
      if (abs(QMMMOpts.StepScale*Forces[i].x) > VecMax)
      {
        VecMax = abs(QMMMOpts.StepScale*Forces[i].x);
      }
      if (abs(QMMMOpts.StepScale*Forces[i].y) > VecMax)
      {
        VecMax = abs(QMMMOpts.StepScale*Forces[i].y);
      }
      if (abs(QMMMOpts.StepScale*Forces[i].z) > VecMax)
      {
        VecMax = abs(QMMMOpts.StepScale*Forces[i].z);
      }
    }
    stepsize = QMMMOpts.StepScale;
    if (VecMax > QMMMOpts.MaxStep)
    {
      //Scale step size
      stepsize *= (QMMMOpts.MaxStep/VecMax);
    }
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        Struct[i].P[Bead].x += stepsize*Forces[ct].x;
        Struct[i].P[Bead].y += stepsize*Forces[ct].y;
        Struct[i].P[Bead].z += stepsize*Forces[ct].z;
        ct += 1;
      }
    }
    Eold = E; //Save energy
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Check convergence
    stepct += 1;
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  sys = system(call.str().c_str());
  //Return for MM optimization
  return;
};

void FLUKEDFP(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //A simple DFP optimizer, which is similar to BFGS updating
  //Note: This optimizer does not have a true line search, instead
  //a steepest descent step is performed if the energy rises
  int sys; //Dummy return for system calls
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
  if ((AMOEBA == 1) and (Gaussian == 1))
  {
    if (TINKER == 1)
    {
      //Set up current charges
      RotateTINKCharges(Struct,Bead);
    }
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
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
  for (int i=0;i<(3*(Nqm+Npseudo));i++)
  {
    //Initialize arrays
    OptVec(i) = 0;
    GradDiff(i) = 0;
    NGrad(i) = 0;
    //Create a scaled identity matrix as the initial Hessian
    for (int j=0;j<i;j++)
    {
      //Set off diagonal terms
      IHess(i,j) = 0.0;
      IHess(j,i) = 0.0;
    }
    IHess(i,i) = 0.1; //Already an "inverse Hessian"
  }
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
    int sys = system("rm -f psi.*");
  }
  //Calculate forces (MM part)
  if (TINKER == 1)
  {
    int tstart = (unsigned)time(0);
    E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
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
  cout << setprecision(8);
  cout << "    QM Step: 0";
  cout << " | RMS force: " << VecMax;
  cout << " eV/\u212B";
  cout << '\n' << '\n';
  cout.flush();
  cout.copyfmt(call); //Return to previous settings
  //Optimize structure
  Eold = E;
  bool OptDone = 0;
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
    OptVec *= QMMMOpts.StepScale;
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
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
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
      int sys = system("rm -f psi.*");
    }
    //Calculate forces (MM part)
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
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
    if (E < Eold)
    {
      //Start really long "line"
      IHess = IHess+((OptVec*OptVec.transpose())/(OptVec.transpose()
      *GradDiff))-((IHess*GradDiff*GradDiff.transpose()*IHess)
      /(GradDiff.transpose()*IHess*GradDiff));
      //End really long "line"
    }
    else
    {
      //Take a small steepest descent step and rebuild Hessian
      cout << "    Energy did not decrease. Constructing new Hessian...";
      cout << '\n';
      for (int i=0;i<(3*(Nqm+Npseudo));i++)
      {
        //Create scaled identity matrix
        for (int j=0;j<i;j++)
        {
          IHess(i,j) = 0.0;
          IHess(j,i) = 0.0;
        }
        IHess(i,i) = 0.10; //Already an "inverse Hessian"
      }
    }
    Eold = E; //Save energy
    //Check convergence
    stepct += 1;
    OptDone = OptConverged(Struct,OldStruct,Forces,stepct,QMMMOpts,Bead,1);
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt_" << Bead << ".xyz";
  call << " MMCharges_" << Bead << ".txt";
  sys = system(call.str().c_str());
  //Return for MM optimization
  return;
};

