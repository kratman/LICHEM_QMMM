/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 Molecular dynamics, thermostat, and barostat functions for FLUKE.

 Reference for the Berendsen thermostat
 

 Reference for the update algorithm
 

*/

//Thermostats
double BerendsenThermo(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Berendsen thermostat to maintain a constant temperature
  int sys;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  double T = 0; //Temperature

  return T;
};

//Update algorithms
void VerletUpdate(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     fstream& outfile, bool ProdRun, int Bead)
{
  //Runs the velocity Verlet algorithm
  double E = 0; //Energy
  int sys;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  int MDSteps;
  double T = 0; //Instantaneous temperature
  double Eavg = 0; //Average energy
  double Tavg = 0; //Average temperature
  int avgct = 0; //Counter
  vector<Coord> Forces; //QM forces
  vector<Coord> MMForces; //MM forces
  //Set up the run
  if (ProdRun == 1)
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
    Coord tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    Forces.push_back(tmp);
  }
  for (int i=0;i<Natoms;i++)
  {
    Coord tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    MMForces.push_back(tmp);
  }
  //Run MD
  for (int n=0;n<MDSteps;n++)
  {
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
      
    }
    //Calculate velocities and delete old forces
    
    //Correct temperature
    T = BerendsenThermo(Struct,QMMMOpts,Bead);
    //Update postions
    
    //Print trajectory
    if ((n == 0) or ((n%QMMMOpts.Nprint) == 0))
    {
      E = 0;
      if (Gaussian == 1)
      {
        int tstart = (unsigned)time(0);
        E += GaussianEnergy(Struct,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4 == 1)
      {
        int tstart = (unsigned)time(0);
        E += PSIEnergy(Struct,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tstart;
        //Clean up annoying useless files
        int sys = system("rm -f psi.*");
      }
      if (TINKER == 1)
      {
        int tstart = (unsigned)time(0);
        E += TINKEREnergy(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER == 1)
      {
        int tstart = (unsigned)time(0);
        E += AMBEREnergy(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      Tavg += T;
      Eavg += E;
      cout << " | MD Step: ";
      cout << (n+1) << " | Temperature: ";
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
