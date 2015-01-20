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

  int MDSteps;
  //Set up the run
  if (ProdRun == 1)
  {
    MDSteps = QMMMOpts.Nsteps;
  }
  else
  {
    MDSteps = QMMMOpts.Neq;
  }
  //Run MD
  for (int n=0;n<MDSteps;n++)
  {
    
  }
  return;
};
