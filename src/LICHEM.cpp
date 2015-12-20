/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM is licensed under GPLv3, for more information see GPL_LICENSE

 References for the LICHEM package:
 Kratz et al., J. Comp. Chem., X, X, X, (2016)

*/

//Primary LICHEM header
#include "LICHEM_headers.h"

int main(int argc, char* argv[])
{
  //Misc. initialization
  StartTime = (unsigned)time(0); //Time the program starts
  srand((unsigned)time(0)); //Serial only random numbers
  //End of section

  //Output stream settings
  //NB: The streams should always be returned to these settings
  cout.precision(16);
  cerr.precision(16);
  //End of section

  //Initialize local variables
  string dummy; //Generic string
  double SumE,SumE2,VolAvg,Ek;
  fstream xyzfile,connectfile,regionfile,outfile; //Input and output files
  vector<QMMMAtom> Struct; //Atom list
  vector<QMMMAtom> OldStruct; //A copy of the atoms list
  vector<QMMMElec> Elecs; //Semi-classical electrons (eFF model)
  QMMMSettings QMMMOpts; //QM and MM wrapper settings
  int randnum; //Random integer
  //End of section

  //Print title and compile date
  PrintFancyTitle();
  cout << '\n';
  cout << "Last modification: ";
  cout << __TIME__ << " on ";
  cout << __DATE__ << '\n';
  cout << '\n';
  cout.flush();
  //End of section

  //Print early messages
  cout << "Reading input..." << '\n';
  cout << '\n';
  cout.flush();
  //End of section

  //Read arguments and look for errors
  ReadArgs(argc,argv,xyzfile,connectfile,regionfile,outfile);
  //End of section

  //Read input and check for errors
  InitializeVariables(QMMMOpts);
  ReadLICHEMInput(xyzfile,connectfile,regionfile,Struct,QMMMOpts);
  //End of section

  //Check input for even more errors
  LICHEMErrorChecker(QMMMOpts);
  LICHEMPrintSettings(QMMMOpts);
  //End of section

  //Center the system
  if (PBCon)
  {
    //PBC corrections that are compatible with all packages
    PBCCenter(Struct,QMMMOpts); //Center the atoms in the box
  }
  //End of section

  /*
    NB: All optional simulation types should be wrapped in comments and
    else-if statements. The first comment should define what calculation is
    going to be performed, then the simulation should be enclosed in an
    else-if statement. After the else-if, an "//End of section" comment
    should be added to mark where the next simulation type begins.
  */

  //Calculate single-point energy
  if (SinglePoint)
  {
    double Eqm = 0; //QM energy
    double Emm = 0; //MM energy
    cout << '\n'; //Print blank line
    cout << "Single-point energy:" << '\n';
    cout.flush(); //Print progress
    //Calculate QM energy
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      Eqm += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      Eqm += PSI4Energy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      Eqm += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (QMMM or QMonly)
    {
      //Print QM partial energy
      cout << "QM energy: " << LICHEMFormDouble(Eqm,16) << " eV";
      cout << '\n';
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      Emm += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      Emm += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      Emm += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    //Print the rest of the energies
    if (QMMM or MMonly)
    {
      //Print MM partial energy
      cout << "MM energy: " << LICHEMFormDouble(Emm,16) << " eV";
      cout << '\n';
    }
    SumE = Eqm+Emm; //Total energy
    if (QMMM)
    {
      //Print total energy
      cout << "QMMM energy: ";
      cout << LICHEMFormDouble(SumE,16) << " eV";
      cout << " ";
      cout << LICHEMFormDouble(SumE/Har2eV,16) << " a.u.";
      cout << '\n';
    }
    cout << '\n';
    cout.flush(); //Print output
  }
  //End of section

  //Optimize structure (native QM and MM package optimizers)
  else if (OptSim)
  {
    //NB: Currently only Gaussian works with this option
    VectorXd Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    cout << "Optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSI4Energy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    cout << " | Opt. step: ";
    cout << optct << " | Energy: ";
    cout << LICHEMFormDouble(SumE,16) << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    //Run optimization
    bool OptDone = 0;
    if (QMMMOpts.MaxOptSteps == 0)
    {
      OptDone = 1;
    }
    while (!OptDone)
    {
      //Copy structure
      OldStruct = Struct;
      //Run MM optimization
      if (TINKER)
      {
        int tstart = (unsigned)time(0);
        SumE = TINKEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER)
      {
        int tstart = (unsigned)time(0);
        SumE = AMBEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS)
      {
        int tstart = (unsigned)time(0);
        SumE = LAMMPSOpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n';
        cout.flush();
      }
      cout << '\n';
      //Run QM optimization
      if (Gaussian)
      {
        int tstart = (unsigned)time(0);
        SumE = GaussianOpt(Struct,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4)
      {
        int tstart = (unsigned)time(0);
        SumE = PSI4Opt(Struct,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tstart;
        //Delete annoying useless files
        GlobalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tstart = (unsigned)time(0);
        SumE = NWChemOpt(Struct,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tstart;
      }
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Steepest descent optimization
  else if (SteepSim)
  {
    VectorXd Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    cout << "Steepest descent optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSI4Energy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.* ");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    cout << " | Opt. step: ";
    cout << optct << " | Energy: ";
    cout << LICHEMFormDouble(SumE,16) << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    //Run optimization
    bool OptDone = 0;
    while (!OptDone)
    {
      //Copy structure
      OldStruct = Struct;
      //Run MM optimization
      if (TINKER)
      {
        int tstart = (unsigned)time(0);
        SumE = TINKEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER)
      {
        int tstart = (unsigned)time(0);
        SumE = AMBEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS)
      {
        int tstart = (unsigned)time(0);
        SumE = LAMMPSOpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n';
        cout.flush();
      }
      cout << '\n';
      //Run QM optimization
      LICHEMSteepest(Struct,QMMMOpts,0);
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Damped Verlet (QuickMin) optimization
  else if (QuickSim)
  {
    VectorXd Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    cout << "Damped Verlet optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSI4Energy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.* ");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    cout << " | Opt. step: ";
    cout << optct << " | Energy: ";
    cout << LICHEMFormDouble(SumE,16) << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    //Run optimization
    bool OptDone = 0;
    while (!OptDone)
    {
      //Copy structure
      OldStruct = Struct;
      //Run MM optimization
      if (TINKER)
      {
        int tstart = (unsigned)time(0);
        SumE = TINKEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER)
      {
        int tstart = (unsigned)time(0);
        SumE = AMBEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS)
      {
        int tstart = (unsigned)time(0);
        SumE = LAMMPSOpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n';
        cout.flush();
      }
      cout << '\n';
      //Run QM optimization
      LICHEMQuickMin(Struct,QMMMOpts,0);
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //DFP minimization
  else if (DFPSim)
  {
    VectorXd Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Change optimization tolerance for the first step
    double SavedQMOptTol = QMMMOpts.QMOptTol; //Save value from input
    double SavedMMOptTol = QMMMOpts.MMOptTol; //Save value from input
    if (QMMMOpts.QMOptTol < 0.005)
    {
      QMMMOpts.QMOptTol = 0.005; //Speedy convergance on the first step
    }
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; //Speedy convergance on the first step
    }
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    cout << "DFP optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSI4Energy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    cout << " | Opt. step: ";
    cout << optct << " | Energy: ";
    cout << LICHEMFormDouble(SumE,16) << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    //Run optimization
    bool OptDone = 0;
    while (!OptDone)
    {
      //Copy structure
      OldStruct = Struct;
      //Run MM optimization
      if (TINKER)
      {
        int tstart = (unsigned)time(0);
        SumE = TINKEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER)
      {
        int tstart = (unsigned)time(0);
        SumE = AMBEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS)
      {
        int tstart = (unsigned)time(0);
        SumE = LAMMPSOpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n';
        cout.flush();
      }
      cout << '\n';
      //Run QM optimization
      LICHEMDFP(Struct,QMMMOpts,0);
      //Reset tolerance before optimization check
      QMMMOpts.QMOptTol = SavedQMOptTol;
      QMMMOpts.MMOptTol = SavedMMOptTol;
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
      if (optct == 1)
      {
        //Avoid terminating restarts on the loose tolerance step
        OptDone = 0; //Not converged
      }
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Run Monte Carlo
  else if (PIMCSim)
  {
    //Adjust probabilities
    if (Natoms == 1)
    {
      //Remove atom centroid moves
      CentProb = 0.0;
      BeadProb = 1.0;
    }
    if (QMMMOpts.Ensemble == "NVT")
    {
      //Remove volume changes
      VolProb = 0.0;
    }
    //Run simulations
    cout << '\n';
    SumE = 0; //Average energy
    SumE2 = 0; //Average squared energy
    VolAvg = 0; //Average volume
    Ek = 0; //PIMC kinietic energy
    if (QMMMOpts.Nbeads > 1)
    {
      //Set kinetic energy
      Ek = 3*Natoms*QMMMOpts.Nbeads/(2*QMMMOpts.Beta);
    }
    //Initialize local variables
    int Nct = 0; //Step counter
    int ct = 0; //Secondary counter
    double Nacc = 0; //Number of accepted moves
    double Nrej = 0; //Number of rejected moves
    double Emc = 0; //Monte Carlo energy
    double Et = 0; //Total energy for printing
    bool acc; //Flag for accepting a step
    //Find the number of characters to print for the step counter
    int SimCharLen;
    SimCharLen = QMMMOpts.Neq+QMMMOpts.Nsteps;
    SimCharLen = LICHEMCountInt(SimCharLen);
    //Start equilibration run and calculate initial energy
    cout << "Monte Carlo equilibration:" << '\n';
    cout.flush();
    QMMMOpts.Eold = 0;
    QMMMOpts.Eold += Get_PI_Epot(Struct,QMMMOpts);
    QMMMOpts.Eold += Get_PI_Espring(Struct,QMMMOpts);
    if (VolProb > 0)
    {
      //Add PV term
      QMMMOpts.Eold += QMMMOpts.Press*Lx*Ly*Lz*atm2eV;
    }
    Emc = QMMMOpts.Eold; //Needed if equilibration is skipped
    Nct = 0;
    while (Nct < QMMMOpts.Neq)
    {
      Emc = 0;
      //Check step size
      if(ct == Acc_Check)
      {
        if ((Nacc/(Nrej+Nacc)) > QMMMOpts.accratio)
        {
          //Increase step size
          double randval;
          randval = (((double)rand())/((double)RAND_MAX));
          //Use random values to keep from cycling up and down
          if (randval >= 0.5)
          {
            step *= 1.10;
          }
          else
          {
            step *= 1.09;
          }
        }
        if ((Nacc/(Nrej+Nacc)) < QMMMOpts.accratio)
        {
          //Decrease step size
          double randval;
          randval = (((double)rand())/((double)RAND_MAX));
          //Use random values to keep from cycling up and down
          if (randval >= 0.5)
          {
            step *= 0.90;
          }
          else
          {
            step *= 0.91;
          }
        }
        if (step < StepMin)
        {
          //Set to minimum
          step = StepMin;
        }
        if (step > StepMax)
        {
          //Set to maximum
          step = StepMax;
        }
        //Statistics
        cout << " | Step: " << setw(SimCharLen) << Nct;
        cout << " | Step size: ";
        cout << LICHEMFormDouble(step,6);
        cout << " | Accept ratio: ";
        cout << LICHEMFormDouble((Nacc/(Nrej+Nacc)),6);
        cout << '\n';
        cout.flush(); //Print stats
        //Reset counters
        ct = 0;
        Nacc = 0;
        Nrej = 0;
      }
      //Continue simulation
      ct += 1;
      acc = MCMove(Struct,QMMMOpts,Emc);
      if (acc)
      {
        Nct += 1;
        Nacc += 1;
      }
      else
      {
        Nrej += 1;
      }
    }
    cout << " Equilibration complete." << '\n';
    //Start production run
    Nct = 0;
    Nacc = 0;
    Nrej = 0;
    cout << '\n';
    cout << "Monte Carlo production:" << '\n';
    cout.flush();
    //Print starting conditions
    Print_traj(Struct,outfile,QMMMOpts);
    Et = Ek+Emc; //Calculate total energy using previous saved energy
    Et -= 2*Get_PI_Espring(Struct,QMMMOpts);
    cout << " | Step: " << setw(SimCharLen) << 0;
    cout << " | Energy: " << LICHEMFormDouble(Et,12);
    cout << " eV";
    if (QMMMOpts.Ensemble == "NPT")
    {
      cout << " | Volume: ";
      cout << LICHEMFormDouble(Lx*Ly*Lz,10);
      cout << " \u212B^3";
    }
    cout << '\n';
    cout.flush(); //Print results
    //Continue simulation
    while (Nct < QMMMOpts.Nsteps)
    {
      Emc = 0; //Set energy to zero
      acc = MCMove(Struct,QMMMOpts,Emc);
      if (acc)
      {
        //Increase counters
        Nct += 1;
        Nacc += 1;
        //Calculate energy
        Et = 0;
        Et += Ek+Emc;
        Et -= 2*Get_PI_Espring(Struct,QMMMOpts);
        VolAvg += Lx*Ly*Lz;
        SumE += Et;
        SumE2 += Et*Et;
        if ((Nct%QMMMOpts.Nprint) == 0)
        {
          //Print progress
          Print_traj(Struct,outfile,QMMMOpts);
          cout << " | Step: " << setw(SimCharLen) << Nct;
          cout << " | Energy: " << LICHEMFormDouble(Et,12);
          cout << " eV";
          if (QMMMOpts.Ensemble == "NPT")
          {
            cout << " | Volume: ";
            cout << LICHEMFormDouble(Lx*Ly*Lz,10);
            cout << " \u212B^3";
          }
          cout << '\n';
          cout.flush(); //Print results
        }
      }
      else
      {
        Nrej += 1;
      }
    }
    if ((Nct%QMMMOpts.Nprint) != 0)
    {
      //Print final geometry if it was not already written
      Print_traj(Struct,outfile,QMMMOpts);
    }
    SumE /= QMMMOpts.Nsteps; //Average energy
    SumE2 /= QMMMOpts.Nsteps; //Variance of the energy
    VolAvg /= QMMMOpts.Nsteps; //Average volume
    //Print simulation details and statistics
    cout << '\n';
    if (QMMMOpts.Nbeads > 1)
    {
      cout << "PI";
    }
    cout << "MC statistics:" << '\n';
    cout << " | Temperature: ";
    cout << QMMMOpts.Temp;
    cout << " K ";
    if (QMMMOpts.Ensemble == "NPT")
    {
      cout << " | Volume: ";
      cout << LICHEMFormDouble(VolAvg,12);
      cout << " \u212B^3";
    }
    cout << '\n';
    cout << " | Average energy: ";
    cout << LICHEMFormDouble(SumE,16);
    cout << " eV | Variance: ";
    cout << LICHEMFormDouble((SumE2-(SumE*SumE)),12);
    cout << " eV^2";
    cout << '\n';
    cout << " | Acceptance ratio: ";
    cout << LICHEMFormDouble((Nacc/(Nrej+Nacc)),6);
    cout << " | Optimum step size: ";
    cout << LICHEMFormDouble(step,6);
    cout << " \u212B";
    cout << '\n';
    cout << '\n';
    cout.flush();
  }
  //End of section

  //NEB optimization
  else if (NEBSim)
  {
    MatrixXd ForceStats; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Check number of beads for Climbing image
    if (QMMMOpts.Nbeads < 4)
    {
      //If the system is well behaved, then the TS bead is known.
      QMMMOpts.Climb = 1; //Turn on climbing image NEB
    }
    //Change optimization tolerance for the first step
    double SavedQMOptTol = QMMMOpts.QMOptTol; //Save value from input
    double SavedMMOptTol = QMMMOpts.MMOptTol; //Save value from input
    if (QMMMOpts.QMOptTol < 0.005)
    {
      QMMMOpts.QMOptTol = 0.005; //Speedy convergance on the first step
    }
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; //Speedy convergance on the first step
    }
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    cout << "Nudged elastic band optimization:" << '\n';
    if (QMMMOpts.Climb)
    {
      cout << " | Short path detected. Starting climbing image NEB.";
      cout << '\n' << '\n';
    }
    cout << " | Opt. step: 0 | Bead energies:";
    cout << '\n';
    cout.flush(); //Print progress
    //Calculate reaction coordinate positions
    VectorXd ReactCoord(QMMMOpts.Nbeads); //Reaction coordinate
    ReactCoord.setZero();
    for (int p=0;p<(QMMMOpts.Nbeads-1);p++)
    {
      MatrixXd Geom1((Nqm+Npseudo),3); //Current replica
      MatrixXd Geom2((Nqm+Npseudo),3); //Next replica
      VectorXd Disp; //Store the displacement
      //Save geometries
      int ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (Struct[i].QMregion or Struct[i].PBregion)
        {
          //Save current replica
          Geom1(ct,0) = Struct[i].P[p].x;
          Geom1(ct,1) = Struct[i].P[p].y;
          Geom1(ct,2) = Struct[i].P[p].z;
          //Save replica p+1
          Geom2(ct,0) = Struct[i].P[p+1].x;
          Geom2(ct,1) = Struct[i].P[p+1].y;
          Geom2(ct,2) = Struct[i].P[p+1].z;
          ct += 1;
        }
      }
      //Calculate displacement
      Disp = KabschDisplacement(Geom1,Geom2,(Nqm+Npseudo));
      ReactCoord(p+1) = ReactCoord(p); //Start from previous bead
      ReactCoord(p+1) += Disp.norm(); //Add magnitude of the displacement
    }
    ReactCoord /= ReactCoord.maxCoeff(); //Must be between 0 and 1
    //Calculate initial energies
    QMMMOpts.Ets = -1*HugeNum; //Locate the initial transition state
    for (int p=0;p<QMMMOpts.Nbeads;p++)
    {
      SumE = 0; //Clear old energies
      //Calculate QM energy
      if (Gaussian)
      {
        int tstart = (unsigned)time(0);
        SumE += GaussianEnergy(Struct,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4)
      {
        int tstart = (unsigned)time(0);
        SumE += PSI4Energy(Struct,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tstart;
        //Delete annoying useless files
        GlobalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tstart = (unsigned)time(0);
        SumE += NWChemEnergy(Struct,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tstart;
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tstart = (unsigned)time(0);
        SumE += TINKEREnergy(Struct,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER)
      {
        int tstart = (unsigned)time(0);
        SumE += AMBEREnergy(Struct,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS)
      {
        int tstart = (unsigned)time(0);
        SumE += LAMMPSEnergy(Struct,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (p == 0)
      {
        //Save reactant energy
        QMMMOpts.Ereact = SumE;
      }
      else if (p == (QMMMOpts.Nbeads-1))
      {
        //Save product energy
        QMMMOpts.Eprod = SumE;
      }
      cout << "   Bead: ";
      cout << setw(3) << p << " | React. coord: ";
      cout << LICHEMFormDouble(ReactCoord(p),5);
      cout << " | Energy: ";
      cout << LICHEMFormDouble(SumE,16) << " eV";
      cout << '\n';
      cout.flush(); //Print progress
      //Update transition state
      if (SumE > QMMMOpts.Ets)
      {
        //Save new properties
        QMMMOpts.TSBead = p;
        QMMMOpts.Ets = SumE;
      }
      //Copy checkpoint data to speed up first step
      if (p != (QMMMOpts.Nbeads-1))
      {
        stringstream call;
        if (Gaussian and (QMMMOpts.Func != "SemiEmp"))
        {
          call.str("");
          call << "cp LICHM_" << (p);
          call << ".chk";
          call << " LICHM_" << (p+1);
          call << ".chk";
          GlobalSys = system(call.str().c_str());
        }
        if (PSI4)
        {
          call.str("");
          call << "cp LICHM_" << (p);
          call << ".32";
          call << " LICHM_" << (p+1);
          call << ".32; ";
          call << "cp LICHM_" << (p);
          call << ".180";
          call << " LICHM_" << (p+1);
          call << ".180";
          GlobalSys = system(call.str().c_str());
        }
      }
    }
    //Run optimization
    bool PathDone = 0;
    int PathStart = 0; //First bead to optimize
    int PathEnd = QMMMOpts.Nbeads; //Last bead to optimize
    if (QMMMOpts.FrznEnds)
    {
      //Change the start and end points
      PathStart = 1;
      PathEnd = QMMMOpts.Nbeads-1;
    }
    while (!PathDone)
    {
      //Copy structure
      OldStruct = Struct;
      //Run MM optimization
      for (int p=PathStart;p<PathEnd;p++)
      {
        if (TINKER)
        {
          int tstart = (unsigned)time(0);
          SumE = TINKEROpt(Struct,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tstart;
        }
        if (AMBER)
        {
          int tstart = (unsigned)time(0);
          SumE = AMBEROpt(Struct,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tstart;
        }
        if (LAMMPS)
        {
          int tstart = (unsigned)time(0);
          SumE = LAMMPSOpt(Struct,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tstart;
        }
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n';
        cout.flush();
      }
      cout << '\n';
      //Run QM optimization
      LICHEMNEB(Struct,QMMMOpts,optct);
      //Reset tolerance before optimization check
      QMMMOpts.QMOptTol = SavedQMOptTol;
      QMMMOpts.MMOptTol = SavedMMOptTol;
      //Print optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      PathDone = PathConverged(Struct,OldStruct,ForceStats,optct,QMMMOpts,0);
      if (optct == 1)
      {
        //Avoid terminating restarts on the loose tolerance step
        PathDone = 0; //Not converged
      }
    }
    BurstTraj(Struct,QMMMOpts);
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    //Print the reaction barriers
    double dEfor = QMMMOpts.Ets-QMMMOpts.Ereact; //Forward barrier
    double dErev = QMMMOpts.Ets-QMMMOpts.Eprod; //Reverse barrier
    cout << "NEB Results:" << '\n';
    cout << " | Transition state bead: " << QMMMOpts.TSBead;
    cout << '\n' << '\n';
    cout << " | Forward barrier: ";
    cout << LICHEMFormDouble(dEfor,16) << " eV ," << '\n';
    cout << " |   ";
    cout << LICHEMFormDouble(dEfor/Har2eV,16) << " a.u.";
    cout << " , ";
    cout << LICHEMFormDouble(dEfor/kcal2eV,16) << " kcal/mol";
    cout << '\n' << '\n';
    cout << " | Reverse barrier: ";
    cout << LICHEMFormDouble(dErev,16) << " eV ," << '\n';
    cout << " |   ";
    cout << LICHEMFormDouble(dErev/Har2eV,16) << " a.u.";
    cout << " , ";
    cout << LICHEMFormDouble(dErev/kcal2eV,16) << " kcal/mol";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Ensemble minimization
  else if (ESDSim)
  {
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    cout << "Ensemble optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tstart = (unsigned)time(0);
      SumE += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4)
    {
      int tstart = (unsigned)time(0);
      SumE += PSI4Energy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tstart = (unsigned)time(0);
      SumE += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER)
    {
      int tstart = (unsigned)time(0);
      SumE += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    cout << " | Opt. step: 0";
    cout << " | Energy: ";
    cout << LICHEMFormDouble(SumE,16);
    cout << " eV" << '\n';
    cout.flush(); //Print progress
    //Run optimization
    EnsembleSD(Struct,outfile,QMMMOpts,0);
    //Finish output
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Ensemble NEB simulation
  else if (ENEBSim)
  {
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Optimize path
    cout << "Ensemble NEB path optimization:" << '\n';
    cout.flush(); //Print progress
    EnsembleNEB(Struct,outfile,QMMMOpts);
    //Finish output
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Inform the user if no simulations were performed
  else
  {
    cout << "Nothing was done..." << '\n';
    cout << "Check the simulation type in " << regfilename;
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Clean up
  xyzfile.close();
  outfile.close();
  regionfile.close();
  connectfile.close();
  if (Gaussian)
  {
    //Clear any remaining Gaussian files
    stringstream call; //Stream for system calls and reading/writing files
    call.str("");
    call << "rm -f Gau-*"; //Produced if there is a crash
    GlobalSys = system(call.str().c_str());
  }
  if (PSI4)
  {
    //Clear any remaining PSI4 files
    stringstream call; //Stream for system calls and reading/writing files
    call.str("");
    call << "rm -f psi*";
    GlobalSys = system(call.str().c_str());
  }
  if (SinglePoint)
  {
    //Clear worthless output xyz file
    stringstream call; //Stream for system calls and reading/writing files
    call.str("");
    call << "rm -f ";
    for (int i=0;i<argc;i++)
    {
      //Find filename
      dummy = string(argv[i]);
      if (dummy == "-o")
      {
        call << argv[i+1];
      }
    }
    GlobalSys = system(call.str().c_str());
  }
  //End of section

  //Print usage statistics
  EndTime = (unsigned)time(0); //Time the program completes
  double TotalHours = (double(EndTime)-double(StartTime));
  double TotalQM = double(QMTime);
  if (PIMCSim and (QMMMOpts.Nbeads > 1))
  {
    //Average over the number of running simulations
    TotalQM /= Nthreads;
  }
  double TotalMM = double(MMTime);
  if (PIMCSim and (QMMMOpts.Nbeads > 1))
  {
    //Average over the number of running simulations
    TotalMM /= Nthreads;
  }
  double OtherTime = TotalHours-TotalQM-TotalMM;
  TotalHours /= 3600.0; //Convert from seconds to hours
  TotalQM /= 3600.0; //Convert from seconds to hours
  TotalMM /= 3600.0; //Convert from seconds to hours
  OtherTime /= 3600.0; //Convert from seconds to hours
  cout << "################# Usage Statistics #################";
  cout << '\n';
  cout << "  Total wall time:                     ";
  cout << LICHEMFormDouble(TotalHours,6) << " hours";
  cout << '\n';
  cout << "  Wall time for QM Wrappers:           ";
  cout << LICHEMFormDouble(TotalQM,6) << " hours";
  cout << '\n';
  cout << "  Wall time for MM Wrappers:           ";
  cout << LICHEMFormDouble(TotalMM,6) << " hours";
  cout << '\n';
  cout << "  Wall time for LICHEM:                ";
  cout << LICHEMFormDouble(OtherTime,6) << " hours";
  cout << '\n';
  cout << "####################################################";
  cout << '\n';
  cout.flush();
  //End of section

  //Print a quote
  if (Jokes)
  {
    cout << '\n';
    cout << "Random quote:";
    cout << '\n';
    string quote; //Random quote
    vector<string> Quotes; //Stores all possible quotes
    GetQuotes(Quotes); //Fetch list of quotes
    randnum = rand() % 1000; //Randomly pick 1 of 1000 quotes
    cout << Quotes[randnum]; //Print quote
    cout << '\n';
  }
  //End of section

  //Finish output
  cout << '\n';
  cout << "Done.";
  cout << '\n';
  cout << '\n';
  cout.flush();
  //End of section

  //Useless but supresses unused return errors for system calls
  int RetValue = GlobalSys;
  RetValue = 0; //This can be changed to error messages later
  //End of section

  //Quit
  return RetValue;
};

