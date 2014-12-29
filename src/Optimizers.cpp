/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 A set of optimization routines for the QM part of QMMM calculculations. These
 are inefficient, however, they are  useful for testing.

*/

void FLUKESteepest(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Steepest descent optimizer
  int sys;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  double RMSdiff = 10000;
  double RMSforce = 10000;
  double MAXforce = 10000;
  int stepct = 0;
  fstream qmfile;
  qmfile.open("QMOpt.xyz",ios_base::out);
  //Optimize structure
  while (((RMSdiff >= QMMMOpts.QMOptTol) or
        (RMSforce >= (100*QMMMOpts.QMOptTol)) or
        (MAXforce >= (200*QMMMOpts.QMOptTol))) and
        (stepct <= QMMMOpts.MaxOptSteps))
  {
    double E = 0;
    RMSforce = 0;
    RMSdiff = 0;
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
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Move QM atoms
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        Struct[i].P[Bead].x += QMMMOpts.SteepStep*Forces[ct].x;
        Struct[i].P[Bead].y += QMMMOpts.SteepStep*Forces[ct].y;
        Struct[i].P[Bead].z += QMMMOpts.SteepStep*Forces[ct].z;
        ct += 1;
      }
    }
    //Print structure
    Print_traj(Struct,qmfile,QMMMOpts);
    //Check convergence
    MAXforce = 0;
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
        double dx = Struct[i].P[Bead].x-OldStruct[i].P[Bead].x;
        double dy = Struct[i].P[Bead].y-OldStruct[i].P[Bead].y;
        double dz = Struct[i].P[Bead].z-OldStruct[i].P[Bead].z;
        RMSdiff += dx*dx+dy*dy+dz*dz;
      }
    }
    RMSdiff /= 3*(Nqm+Npseudo);
    RMSdiff = sqrt(RMSdiff);
    RMSforce /= 3*(Nqm+Npseudo);
    RMSforce = sqrt(RMSforce);
    stepct += 1;
    //Print progress
    call.copyfmt(cout); //Save settings
    cout << setprecision(8);
    cout << "    QM Step: " << (stepct-1);
    cout << " | RMS Disp: " << RMSdiff;
    cout << '\n';
    cout << "    Max force: " << MAXforce;
    cout << " | RMS force: " << RMSforce;
    cout << '\n' << endl;
    cout.copyfmt(call); //Return to previous settings
  }
  //Calculate new point charges
  if (QMMM == 1)
  {
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      GaussianCharges(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      PSICharges(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      int sys = system("rm -f psi.*");
    }
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt.xyz";
  sys = system(call.str().c_str());
  //Return for MM optimization
  return;
};

void FLUKEBFGS(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //BFGS optimizer for QM atoms
  double E = 0; //Energy
  int sys;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  double RMSdiff = 10000;
  double RMSforce = 10000;
  double MAXforce = 10000;
  int stepct = 0;
  fstream qmfile;
  qmfile.open("QMOpt.xyz",ios_base::out);
  //Create BFGS arrays
  VectorXd OptVec(3*(Nqm+Npseudo));
  VectorXd GradDiff(3*(Nqm+Npseudo));
  VectorXd NGrad(3*(Nqm+Npseudo));
  MatrixXd Hess(3*(Nqm+Npseudo),3*(Nqm+Npseudo));
  for (int i=0;i<(3*(Nqm+Npseudo));i++)
  {
    //Initialize arrays
    OptVec(i) = 0;
    GradDiff(i) = 0;
    NGrad(i) = 0;
    //Create matrix for the initial Hessian
    for (int j=0;j<(3*(Nqm+Npseudo));j++)
    {
      //Set off diagonal terms
      Hess(i,j) = 0.0;
    }
    Hess(i,i) = 10.0; //Scale diagonal identity matrix for small initial steps
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
  //Determine new structure
  int ct = 0; //Counter
  int ct2 = 0; //Secondary counter
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Change forces array
    NGrad(ct2) = Forces[i].x;
    NGrad(ct2+1) = Forces[i].y;
    NGrad(ct2+2) = Forces[i].z;
    ct2 += 3;
  }
  //Optimize structure
  while (((RMSdiff >= QMMMOpts.QMOptTol) or
        (RMSforce >= (100*QMMMOpts.QMOptTol)) or
        (MAXforce >= (200*QMMMOpts.QMOptTol))) and
        (stepct <= QMMMOpts.MaxOptSteps))
  {
    RMSforce = 0;
    RMSdiff = 0;
    MAXforce = 0;
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
    OptVec = Hess.inverse()*NGrad;
    OptVec *= QMMMOpts.SteepStep;
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
    //Start really long "line"
    Hess = Hess+((GradDiff*GradDiff.transpose())/(GradDiff.transpose()*
    OptVec))-((Hess*OptVec*OptVec.transpose()*Hess)/(OptVec.transpose()*
    Hess*OptVec)); //End really long line
    //Check convergence
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
        double dx = Struct[i].P[Bead].x-OldStruct[i].P[Bead].x;
        double dy = Struct[i].P[Bead].y-OldStruct[i].P[Bead].y;
        double dz = Struct[i].P[Bead].z-OldStruct[i].P[Bead].z;
        RMSdiff += dx*dx+dy*dy+dz*dz;
      }
    }
    RMSdiff /= 3*(Nqm+Npseudo);
    RMSdiff = sqrt(RMSdiff);
    RMSforce /= 3*(Nqm+Npseudo);
    RMSforce = sqrt(RMSforce);
    stepct += 1;
    //Print progress
    call.copyfmt(cout); //Save settings
    cout << setprecision(8);
    cout << "    QM Step: " << (stepct-1);
    cout << " | RMS Disp: " << RMSdiff;
    cout << '\n';
    cout << "    Max force: " << MAXforce;
    cout << " | RMS force: " << RMSforce;
    cout << '\n' << endl;
    cout.copyfmt(call); //Return to previous settings
  }
  //Calculate new point charges
  if (QMMM == 1)
  {
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      GaussianCharges(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      PSICharges(Struct,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      int sys = system("rm -f psi.*");
    }
  }
  //Clean up files
  call.str("");
  call << "rm -f QMOpt.xyz";
  sys = system(call.str().c_str());
  //Return for MM optimization
  return;
};

