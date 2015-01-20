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
void BerendsenThermo(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Berendsen thermostat to maintain a constant temperature

  return;
};

//Update algorithms
void VerletUpdate(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     bool ProdRun, int Bead)
{
  //Runs the velocity Verlet algorithm
  double E = 0; //Energy
  int sys;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  int MDSteps;
  vector<Coord> Forces; //QM forces
  vector<Coord> MMForces; //MM forces
  //Set up the run
  if (ProdRun == 1)
  {
    MDSteps = QMMMOpts.Nsteps;
  }
  else
  {
    MDSteps = QMMMOpts.Neq;
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
    BerendsenThermo(Struct,QMMMOpts,Bead);
    //Update postions
    
    //Print trajectory
    
  }
  return;
};
