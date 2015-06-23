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

*/

//Main header
#include "QMMM_headers.h"

int main(int argc, char* argv[])
{
  //Misc. initialization
  StartTime = (unsigned)time(0); //Time the program starts
  cout.precision(12);
  //End of section

  //Initialize local variables
  srand((unsigned)time(0)); //Serial only random numbers
  string dummy; //Generic string
  double SumE,SumE2,VolAvg,Ek;
  fstream xyzfile,connectfile,regionfile,outfile; //Input and output files
  vector<QMMMAtom> Struct; //Atom list
  vector<QMMMAtom> OldStruct; //A copy of the atoms list
  vector<QMMMElec> Elecs; //Semi-classical electrons
  QMMMSettings QMMMOpts; //QM wrapper settings
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
    double Eqm = 0;
    double Emm = 0;
    cout << fixed;
    //Calculate QM energy
    if (Gaussian == 1)
    {
      int tstart = (unsigned)time(0);
      Eqm += GaussianEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (PSI4 == 1)
    {
      int tstart = (unsigned)time(0);
      Eqm += PSIEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
      //Clean up annoying useless files
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      Eqm += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    if (QMMM or QMonly)
    {
      //Print QM partial energy
      cout << "QM energy: " << Eqm << " eV";
      cout << '\n';
    }
    //Calculate MM energy
    if (TINKER == 1)
    {
      int tstart = (unsigned)time(0);
      Emm += TINKEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (AMBER == 1)
    {
      int tstart = (unsigned)time(0);
      Emm += AMBEREnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      Emm += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    if (QMMM or MMonly)
    {
      //Print MM partial energy
      cout << "MM energy: " << Emm << " eV";
      cout << '\n';
    }
    SumE = Eqm+Emm;
    if (QMMM)
    {
      //Print total energy
      cout << "QMMM energy: ";
      cout << SumE << " eV";
      cout << " ";
      cout << SumE/Har2eV << " a.u.";
    }
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
      VolProb = 0.0;
    }
    //Calculate initial energy
    QMMMOpts.Eold = 0;
    QMMMOpts.Eold += Get_PI_Epot(Struct,QMMMOpts);
    QMMMOpts.Eold += Get_PI_Espring(Struct,QMMMOpts);
    if (VolProb > 0)
    {
      QMMMOpts.Eold += QMMMOpts.Press*Lx*Ly*Lz*atm2eV;
    }
    //Run simulations
    cout << '\n';
    SumE = 0;
    SumE2 = 0;
    VolAvg = 0;
    Ek = 0;
    if (QMMMOpts.Nbeads > 1)
    {
      //PIMC kinetic energy
      Ek = 3*Natoms*QMMMOpts.Nbeads/(2*QMMMOpts.Beta);
    }
    int Nct = 0; //Step counter
    int ct = 0; //Secondary counter
    double Nacc = 0; //Number of accepted moves
    double Nrej = 0; //Number of rejected moves
    double Emc = 0; //Monte Carlo energy
    bool acc; //Flag for accepting a step
    cout << "Starting equilibration..." << '\n';
    cout.flush();
    Nct = 0;
    while (Nct < QMMMOpts.Neq) //Equilibration
    {
      Emc = 0;
      if(ct == Acc_Check)
      {
        if ((Nacc/(Nrej+Nacc)) > QMMMOpts.accratio)
        {
          step *= 1.10;
        }
        if ((Nacc/(Nrej+Nacc)) < QMMMOpts.accratio)
        {
          step *= 0.91;
        }
        if (step < StepMin)
        {
          step = StepMin;
        }
        if (step > StepMax)
        {
          step = StepMax;
        }
        ct = 0;
        Nacc = 0;
        Nrej = 0;
      }
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
    Nct = 0;
    Nacc = 0;
    Nrej = 0;
    cout << "Starting production run..." << '\n';
    cout.flush();
    Print_traj(Struct,outfile,QMMMOpts);
    while (Nct < QMMMOpts.Nsteps)
    {
      Emc = 0; //Set energy to zero
      acc = MCMove(Struct,QMMMOpts,Emc);
      if (acc)
      {
        Nct += 1;
        Nacc += 1;
        double Et = Ek+Emc;
        Et -= 2*Get_PI_Espring(Struct,QMMMOpts);
        VolAvg += Lx*Ly*Lz;
        SumE += Et;
        SumE2 += Et*Et;
        if ((Nct%QMMMOpts.Nprint) == 0)
        {
          //Print progress
          Print_traj(Struct,outfile,QMMMOpts);
          cout << " | Step: " << Nct;
          cout << " | Energy: " << Et << " eV";
          if (QMMMOpts.Ensemble == "NPT")
          {
            cout << " | Volume: " << Lx*Ly*Lz << " \u212B^3";
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
    //Print output
    cout << '\n';
    cout << "Temperature: ";
    cout << QMMMOpts.Temp;
    cout << " K    ";
    if (QMMMOpts.Ensemble == "NPT")
    {
      cout << "Volume: ";
      cout << VolAvg;
      cout << " \u212B^3";
    }
    cout << '\n';
    cout << "Average energy: ";
    cout << SumE;
    cout << " eV    ";
    cout << "Variance: ";
    cout << (SumE2-SumE*SumE);
    cout << " eV^2";
    cout << '\n';
    cout << "Acceptance ratio: ";
    cout << (Nacc/(Nrej+Nacc));
    cout << "    ";
    cout << "Optimum step size: ";
    cout << step;
    cout << " \u212B";
    cout << '\n';
    cout << '\n';
    cout.flush();
  }
  //End of section

  //Optimize structure
  else if (OptSim)
  {
    vector<Coord> Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
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
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
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
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    stringstream call;
    call.copyfmt(cout); //Save settings
    cout << " | Opt. step: ";
    cout << optct << " | Energy: ";
    cout << setprecision(12) << SumE << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    cout.copyfmt(call); //Replace settings
    //Run optimization
    bool OptDone = 0;
    while (!OptDone)
    {
      //Copy structure
      OldStruct = Struct;
      //Run MM optimization
      if (TINKER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = TINKEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = AMBEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = LAMMPSOpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n' << '\n';
        cout.flush();
      }
      //Run QM optimization
      if (Gaussian == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = GaussianOpt(Struct,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4 == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = PSIOpt(Struct,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tstart;
        //Clean up annoying useless files
        GlobalSys = system("rm -f psi.*");
      }
      if (NWChem == 1)
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
    vector<Coord> Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
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
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
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
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    stringstream call;
    call.copyfmt(cout); //Save settings
    cout << " | Opt. step: ";
    cout << optct << " | Energy: ";
    cout << setprecision(12) << SumE << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    cout.copyfmt(call); //Replace settings
    //Run optimization
    bool OptDone = 0;
    while (!OptDone)
    {
      //Copy structure
      double SavedStepSize = QMMMOpts.StepScale; //Save old step size
      OldStruct = Struct;
      //Run MM optimization
      if (TINKER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = TINKEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = AMBEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = LAMMPSOpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n' << '\n';
        cout.flush();
      }
      //Run QM optimization
      LICHEMSteepest(Struct,QMMMOpts,0);
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
      QMMMOpts.StepScale = SavedStepSize;
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
    vector<Coord> Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
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
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
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
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    stringstream call;
    call.copyfmt(cout); //Save settings
    cout << " | Opt. step: ";
    cout << optct << " | Energy: ";
    cout << setprecision(12) << SumE << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    cout.copyfmt(call); //Replace settings
    //Run optimization
    bool OptDone = 0;
    while (!OptDone)
    {
      //Copy structure
      double SavedStepSize = QMMMOpts.StepScale; //Save old step size
      OldStruct = Struct;
      //Run MM optimization
      if (TINKER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = TINKEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = AMBEROpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS == 1)
      {
        int tstart = (unsigned)time(0);
        SumE = LAMMPSOpt(Struct,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n' << '\n';
        cout.flush();
      }
      //Run QM optimization
      LICHEMDFP(Struct,QMMMOpts,0);
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
      QMMMOpts.StepScale = SavedStepSize;
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //NEB optimization
  else if (NEBSim)
  {
    vector<Coord> Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energies
    for (int p=0;p<QMMMOpts.Nbeads;p++)
    {
      SumE = 0; //Clear old energies
      //Calculate QM energy
      if (Gaussian == 1)
      {
        int tstart = (unsigned)time(0);
        SumE += GaussianEnergy(Struct,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tstart;
      }
      if (PSI4 == 1)
      {
        int tstart = (unsigned)time(0);
        SumE += PSIEnergy(Struct,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tstart;
        //Clean up annoying useless files
        GlobalSys = system("rm -f psi.*");
      }
      if (NWChem == 1)
      {
        int tstart = (unsigned)time(0);
        SumE += NWChemEnergy(Struct,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tstart;
      }
      //Calculate MM energy
      if (TINKER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE += TINKEREnergy(Struct,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (AMBER == 1)
      {
        int tstart = (unsigned)time(0);
        SumE += AMBEREnergy(Struct,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tstart;
      }
      if (LAMMPS == 1)
      {
        int tstart = (unsigned)time(0);
        SumE += LAMMPSEnergy(Struct,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tstart;
      }
      stringstream call;
      call.copyfmt(cout); //Save settings
      cout << "   Bead: ";
      cout << p << " | Energy: ";
      cout << setprecision(12) << SumE << " eV";
      cout << '\n';
      cout.flush(); //Print progress
      cout.copyfmt(call); //Replace settings
      //Copy checkpoint data to speed up first step
      if (p != (QMMMOpts.Nbeads-1))
      {
        if (Gaussian == 1)
        {
          call.str("");
          call << "cp QMMM_" << (p);
          call << ".chk";
          call << " QMMM_" << (p+1);
          call << ".chk";
          GlobalSys = system(call.str().c_str());
        }
      }
    }
    //Run optimization
    bool OptDone = 0;
    while (!OptDone)
    {
      //Copy structure
      double SavedStepSize = QMMMOpts.StepScale; //Save old step size
      OldStruct = Struct;
      //Run MM optimization
      for (int p=0;p<QMMMOpts.Nbeads;p++)
      {
        if (TINKER == 1)
        {
          int tstart = (unsigned)time(0);
          SumE = TINKEROpt(Struct,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tstart;
        }
        if (AMBER == 1)
        {
          int tstart = (unsigned)time(0);
          SumE = AMBEROpt(Struct,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tstart;
        }
        if (LAMMPS == 1)
        {
          int tstart = (unsigned)time(0);
          SumE = LAMMPSOpt(Struct,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tstart;
        }
        if (QMMM)
        {
          cout << "    MM optimization complete.";
          cout << '\n' << '\n';
          cout.flush();
        }
      }
      //Run QM optimization
      LICHEMNEB(Struct,QMMMOpts);
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
      QMMMOpts.StepScale = SavedStepSize;
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Ensemble minimization
  else if (ESDSim)
  {
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energy
    SumE = 0; //Clear old energies
    //Calculate QM energy
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
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += NWChemEnergy(Struct,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tstart;
    }
    //Calculate MM energy
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
    if (LAMMPS == 1)
    {
      int tstart = (unsigned)time(0);
      SumE += LAMMPSEnergy(Struct,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tstart;
    }
    stringstream call;
    call.copyfmt(cout); //Save settings
    cout << " | Opt. step: 0";
    cout << " | Energy: ";
    cout << setprecision(16) << SumE;
    cout << " eV" << '\n';
    cout.flush(); //Print progress
    cout.copyfmt(call); //Replace settings
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
  if (Gaussian == 1)
  {
    //Clear any remaining Gaussian files
    stringstream call;
    call.copyfmt(cout);
    call.str("");
    call << "rm -f Gau-*"; //Produced if there is a crash
    GlobalSys = system(call.str().c_str());
  }
  if (PSI4 == 1)
  {
    //Clear any remaining PSI4 files
    stringstream call;
    call.copyfmt(cout);
    call.str("");
    call << "rm -f psi*";
    GlobalSys = system(call.str().c_str());
  }
  if (SinglePoint)
  {
    //Clear worthless output xyz file
    stringstream call;
    call.copyfmt(cout);
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
  if (QMMMOpts.Nbeads > 1)
  {
    //Average over the number of beads
    TotalQM /= QMMMOpts.Nbeads;
  }
  double TotalMM = double(MMTime);
  if (QMMMOpts.Nbeads > 1)
  {
    //Average over the number of beads
    TotalMM /= QMMMOpts.Nbeads;
  }
  double OtherTime = TotalHours-TotalQM-TotalMM;
  TotalHours /= 3600.0; //Convert from seconds to hours
  TotalQM /= 3600.0; //Convert from seconds to hours
  TotalMM /= 3600.0; //Convert from seconds to hours
  OtherTime /= 3600.0; //Convert from seconds to hours
  cout << "################# Usage Statistics #################";
  cout << '\n';
  cout << "  Total wall time: ";
  cout << TotalHours << " hours";
  cout << '\n';
  cout << "  Wall time for QM Wrappers: ";
  cout << TotalQM << " hours";
  cout << '\n';
  cout << "  Wall time for MM Wrappers: ";
  cout << TotalMM << " hours";
  cout << '\n';
  cout << "  Wall time for LICHEM: ";
  cout << OtherTime << " hours";
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

  //Quit
  cout << '\n';
  cout << "Done.";
  cout << '\n';
  cout << '\n';
  cout.flush();

  //Useless but supresses unused return errors for system calls
  int RetValue = GlobalSys;
  RetValue = 0; //This can be changed to error messages later
  return RetValue;
};
