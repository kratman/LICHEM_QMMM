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
  double RMSdiff = 1;
  double RMSforce = 1;
  double MAXforce = 1;
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
    if (Psi4 == 1)
    {
      int tstart = (unsigned)time(0);
      E += PsiForces(Struct,Forces,QMMMOpts,Bead);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      int sys = system("rm -f psi.*");
    }
    //Calculate forces (MM part)
    if (Tinker == 1)
    {
      int tstart = (unsigned)time(0);
      E += TinkerForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Determine new structure
    int ct = 0; //Counter
    for (int i=0;i<Natoms;i++)
    {
      //Calculate RMS displacement
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        Struct[i].x += QMMMOpts.SteepStep*Forces[ct].x;
        Struct[i].y += QMMMOpts.SteepStep*Forces[ct].y;
        Struct[i].z += QMMMOpts.SteepStep*Forces[ct].z;
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
        double dx = Struct[i].x-OldStruct[i].x;
        double dy = Struct[i].y-OldStruct[i].y;
        double dz = Struct[i].z-OldStruct[i].z;
        RMSdiff += dx*dx+dy*dy+dz*dz;
      }
    }
    RMSdiff /= 3*(Nqm+Npseudo);
    RMSdiff = sqrt(RMSdiff);
    RMSforce /= 3*(Nqm+Npseudo);
    RMSforce = sqrt(RMSforce);
    stepct += 1;
    //Print progress
    cout << "    QM Step: " << (stepct-1);
    cout << " | RMS Disp: " << RMSdiff;
    cout << '\n';
    cout << "    Max force: " << MAXforce;
    cout << " | RMS force: " << RMSforce;
    cout << '\n' << endl;
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
    if (Psi4 == 1)
    {
      int tstart = (unsigned)time(0);
      PsiCharges(Struct,QMMMOpts,Bead);
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
  //BFGS optimizer
  
  return;
};

