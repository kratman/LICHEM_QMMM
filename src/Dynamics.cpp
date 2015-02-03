/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 Molecular dynamics, thermostat, and barostat functions for FLUKE.

 Reference for the Berendsen thermostat:
 Berendsen et al. J. Chem. Phys. 81, 8, 3684, 1984

 Reference for the update algorithm:
 

*/

//Thermostats
double BerendsenThermo(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Berendsen thermostat to maintain a constant temperature and remove
  //center of mass velocities
  int sys; //Dummy return for system calls
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  double T = 0; //Temperature
  double Ek = 0; //Kinetic energy
  double VelScale = 0; //Scale factor for velocities
  //Calculate center of mass translation
  double vxcom = 0;
  double vycom = 0;
  double vzcom = 0;
  for (int i=0;i<Natoms;i++)
  {
    vxcom += Struct[i].Vel[Bead].x;
    vycom += Struct[i].Vel[Bead].y;
    vzcom += Struct[i].Vel[Bead].z;
  }
  vxcom /= Natoms;
  vxcom /= Natoms;
  vxcom /= Natoms;
  #pragma omp parallel for
  for (int i=0;i<Natoms;i++)
  {
    //Remove center of mass translation
    Struct[i].Vel[Bead].x -= vxcom;
    Struct[i].Vel[Bead].y -= vycom;
    Struct[i].Vel[Bead].z -= vzcom;
  }
  #pragma omp barrier
  //Calculate temperature
  for (int i=0;i<Natoms;i++)
  {
    //Calculate the kinetic energy
    double v2 = 0;
    v2 += Struct[i].Vel[Bead].x*Struct[i].Vel[Bead].x;
    v2 += Struct[i].Vel[Bead].y*Struct[i].Vel[Bead].y;
    v2 += Struct[i].Vel[Bead].z*Struct[i].Vel[Bead].z;
    Ek += Struct[i].m*v2; //Two left out below
  }
  T = (Ek*amu2kg)/(3*Natoms*kSI*m2Ang*m2Ang*fs2s*fs2s); //Two left out above
  //Calculate the scale factor for the velocities
  VelScale = (T/QMMMOpts.Temp);
  VelScale -= 1;
  VelScale *= (QMMMOpts.dt/QMMMOpts.tautemp);
  VelScale += 1;
  VelScale = sqrt(VelScale);
  //Scale velocities
  #pragma omp parallel for
  for (int i=0;i<Natoms;i++)
  {
    Struct[i].Vel[Bead].x *= VelScale;
    Struct[i].Vel[Bead].y *= VelScale;
    Struct[i].Vel[Bead].z *= VelScale;
  }
  #pragma omp barrier
  //Return
  return T;
};

//Update algorithms
void VerletUpdate(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     fstream& outfile, bool ProdRun, int Bead)
{
  //Runs the velocity Verlet algorithm
  double E = 0; //Energy
  int sys; //Dummy return for system calls
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  int ct; //Generic counter
  int MDSteps; //Number of steps for the MD simulation
  double T = 0; //Instantaneous temperature
  double Eavg = 0; //Average energy
  double Tavg = 0; //Average temperature
  int avgct = 0; //Another counter
  vector<Coord> Forces; //QM forces
  vector<Coord> MMForces; //MM forces
  double AccelConst; //Conversion constant for forces
  AccelConst = ((m2Ang*m2Ang*fs2s*fs2s)/(amu2kg*SI2eV));
  double VelConst; //Conversion constant for velocity
  VelConst = ((m2Ang*m2Ang*fs2s)/(amu2kg*SI2eV));
  //Set up the run
  if (ProdRun)
  {
    MDSteps = QMMMOpts.Nsteps;
    cout << "Starting production run:" << '\n';
    cout << '\n';
  }
  else
  {
    MDSteps = QMMMOpts.Neq;
    cout << "Starting equilibration:" << '\n';
    cout << '\n';
  }
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Create forces array for QM and PA regions
    Coord tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    Forces.push_back(tmp);
  }
  for (int i=0;i<Natoms;i++)
  {
    //Create forces array for MM, BA, and PA regions
    Coord tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    MMForces.push_back(tmp);
  }
  //Run MD
  for (int n=0;n<(MDSteps+1);n++)
  {
    E = 0;
    //Update QM forces
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
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      E += TINKERForces(Struct,Forces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Update MM forces
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      E += TINKERMMForces(Struct,MMForces,QMMMOpts,Bead);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Sum forces and delete old QM forces
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion == 1)
      {
        //Use the QM forces
        MMForces[i].x = Forces[ct].x;
        MMForces[i].y = Forces[ct].y;
        MMForces[i].z = Forces[ct].z;
        Forces[ct].x = 0;
        Forces[ct].y = 0;
        Forces[ct].z = 0;
        ct += 1;
      }
      if (Struct[i].PAregion == 1)
      {
        //Take the average of the MM and QM forces
        MMForces[i].x += Forces[ct].x;
        MMForces[i].y += Forces[ct].y;
        MMForces[i].z += Forces[ct].z;
        MMForces[i].x *= 0.5;
        MMForces[i].y *= 0.5;
        MMForces[i].z *= 0.5;
        Forces[ct].x = 0;
        Forces[ct].y = 0;
        Forces[ct].z = 0;
        ct += 1;
      }
    }
    //Update postions and delete old MM forces
    #pragma omp parallel for
    for (int i=0;i<Natoms;i++)
    {
      //Update from velocity
      Struct[i].P[Bead].x += Struct[i].Vel[Bead].x*QMMMOpts.dt;
      Struct[i].P[Bead].y += Struct[i].Vel[Bead].y*QMMMOpts.dt;
      Struct[i].P[Bead].z += Struct[i].Vel[Bead].z*QMMMOpts.dt;
      //Update from acceleration (multiline)
      Struct[i].P[Bead].x += 0.5*MMForces[i].x*QMMMOpts.dt*QMMMOpts.dt
      *AccelConst/Struct[i].m;
      Struct[i].P[Bead].y += 0.5*MMForces[i].y*QMMMOpts.dt*QMMMOpts.dt
      *AccelConst/Struct[i].m;
      Struct[i].P[Bead].z += 0.5*MMForces[i].z*QMMMOpts.dt*QMMMOpts.dt
      *AccelConst/Struct[i].m;
    }
    #pragma omp barrier
    //Update velocities and delete old forces
    #pragma omp parallel for
    for (int i=0;i<Natoms;i++)
    {
      Struct[i].Vel[Bead].x += 0.5*MMForces[i].x*QMMMOpts.dt*VelConst
      /Struct[i].m;
      Struct[i].Vel[Bead].y += 0.5*MMForces[i].y*QMMMOpts.dt*VelConst
      /Struct[i].m;
      Struct[i].Vel[Bead].z += 0.5*MMForces[i].z*QMMMOpts.dt*VelConst
      /Struct[i].m;
      //Delete old MM forces
      MMForces[i].x = 0;
      MMForces[i].y = 0;
      MMForces[i].z = 0;
    }
    #pragma omp barrier
    //Correct temperature
    T = BerendsenThermo(Struct,QMMMOpts,Bead);
    //Print trajectory
    if ((n == 0) or ((n%QMMMOpts.Nprint) == 0))
    {
      E = 0;
      if (Gaussian == 1)
      {
        int tstart = (unsigned)time(0);
        E += GaussianEnergy(Struct,QMMMOpts,Bead);
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4 == 1)
      {
        int tstart = (unsigned)time(0);
        E += PSIEnergy(Struct,QMMMOpts,Bead);
        QMTime += (unsigned)time(0)-tstart;
        //Clean up annoying useless files
        int sys = system("rm -f psi.*");
      }
      if (TINKER == 1)
      {
        int tstart = (unsigned)time(0);
        E += TINKEREnergy(Struct,QMMMOpts,Bead);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER == 1)
      {
        int tstart = (unsigned)time(0);
        E += AMBEREnergy(Struct,QMMMOpts,Bead);
        MMTime += (unsigned)time(0)-tstart;
      }
      Tavg += T;
      Eavg += E;
      cout << " | MD Step: ";
      cout << n << " | Temperature: ";
      cout << T << " K | Energy: ";
      cout << E << " eV";
      cout << endl; //Print progress
      if (ProdRun == 1)
      {
        Print_traj(Struct,outfile,QMMMOpts);
      }
      avgct += 1;
    }
  }
  Eavg /= avgct;
  Tavg /= avgct;
  cout << '\n';
  cout << "MD simulation complete.";
  cout << '\n' << '\n';
  cout << "Average energy: ";
  cout << Eavg << " eV | Average temperature: ";
  cout << Tavg << " K";
  cout << '\n' << '\n';
  return;
};

