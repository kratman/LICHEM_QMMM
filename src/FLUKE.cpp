/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE is licensed under GPLv3, for more information see GPL_LICENSE

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

  //Print title
  PrintFancyTitle();
  cout << '\n';
  cout << "Last modification: ";
  cout << __TIME__ << " on ";
  cout << __DATE__ << '\n';
  cout << endl;
  //End of section

  //Read arguments and look for errors
  ReadArgs(argc,argv,xyzfile,connectfile,regionfile,outfile);
  //End of section

  //Read input and check for errors
  ReadFLUKEInput(xyzfile,connectfile,regionfile,Struct,QMMMOpts);
  //End of section

  //Check input for even more errors
  FLUKEErrorChecker(QMMMOpts);
  FLUKEPrintSettings(QMMMOpts);
  //End of section

  //Find induced dipoles for polarizable simulations
  if ((AMOEBA == 1) and (QMMM) and (!PIMCSim))
  {
    if (Gaussian == 1)
    {
      #pragma omp parallel for
      for (int i=0;i<QMMMOpts.Nbeads;i++)
      {
        GaussianCharges(Struct,QMMMOpts,i);
      }
      #pragma omp barrier
    }
    if (PSI4 == 1)
    {
      #pragma omp parallel for
      for (int i=0;i<QMMMOpts.Nbeads;i++)
      {
        PSICharges(Struct,QMMMOpts,i);
      }
      #pragma omp barrier
    }
    if (TINKER == 1)
    {
      #pragma omp parallel for
      for (int i=0;i<QMMMOpts.Nbeads;i++)
      {
        TINKERInduced(Struct,QMMMOpts,i);
      }
      #pragma omp barrier
    }
  }
  //End of section

  //Calculate single-point energy (optional)
  if (SinglePoint == 1)
  {
    double Eqm = 0;
    double Emm = 0;
    cout << fixed;
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
      int sys = system("rm -f psi.*");
    }
    if ((QMMM == 1) or (QMonly == 1))
    {
      //Print QM partial energy
      cout << "QM energy: " << Eqm << " eV";
      cout << endl;
    }
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
    if ((QMMM == 1) or (MMonly == 1))
    {
      //Print MM partial energy
      cout << "MM energy: " << Emm << " eV";
      cout << endl;
    }
    SumE = Eqm+Emm;
    if (QMMM == 1)
    {
      //Print total energy
      cout << "QMMM total energy: ";
      cout << SumE << " eV";
      cout << " ";
      cout << SumE/Har2eV << " a.u.";
    }
    cout << '\n' << endl;
  }
  //End of section

  //Run Monte Carlo (optional)
  if (PIMCSim == 1)
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
    //Run simulations
    cout << '\n';
    SumE = 0;
    SumE2 = 0;
    VolAvg = 0;
    Ek = 3*Natoms*QMMMOpts.Nbeads/(2*QMMMOpts.Beta);
    int Nct = 0; //Step counter
    int ct = 0; //Secondary counter
    double Nacc = 0; //Number of accepted moves
    double Nrej = 0; //Number of rejected moves
    double Emc = 0; //Monte Carlo energy
    bool acc; //Flag for accepting a step
    cout << "Starting equilibration..." << endl;
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
    cout << "Starting production run..." << endl;
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
    cout << "Average Energy: ";
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
    cout << endl;
  }
  //End of section

  //Run MD simulation
  if (MDSim == 1)
  {
    //Equilibration
    VerletUpdate(Struct,QMMMOpts,outfile,0,0);
    //Production
    VerletUpdate(Struct,QMMMOpts,outfile,1,0);
  }
  //End of section

  //Optimize structure (optional)
  if (OptSim == 1)
  {
    vector<Coord> Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energy
    SumE = 0; //Clear old energies
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
    cout << " | Opt. Step: ";
    cout << optct << " | Energy: ";
    cout << SumE << " eV";
    cout << endl; //Print progress
    //Run optimization
    bool OptDone = 0;
    while (OptDone == 0)
    {
      //Copy structure
      OldStruct = Struct;
      //Run optimization
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
        int sys = system("rm -f psi.*");
      }
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
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << endl;
  }
  //End of section

  //Steepest descent optimization
  if (SteepSim == 1)
  {
    vector<Coord> Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energy
    SumE = 0; //Clear old energies
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
    cout << " | Opt. Step: ";
    cout << optct << " | Energy: ";
    cout << SumE << " eV";
    cout << endl; //Print progress
    //Run optimization
    bool OptDone = 0;
    while (OptDone == 0)
    {
      //Copy structure
      OldStruct = Struct;
      //Run optimization
      FLUKESteepest(Struct,QMMMOpts,0);
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
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << endl;
  }
  //End of sections

  //DFP minimization
  if (DFPSim == 1)
  {
    vector<Coord> Forces; //Dummy array needed for convergence tests
    int optct = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(Struct,outfile,QMMMOpts);
    //Calculate initial energy
    SumE = 0; //Clear old energies
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
    cout << " | Opt. Step: ";
    cout << optct << " | Energy: ";
    cout << SumE << " eV";
    cout << endl; //Print progress
    //Run optimization
    bool OptDone = 0;
    while (OptDone == 0)
    {
      //Copy structure
      OldStruct = Struct;
      //Run optimization
      FLUKEDFP(Struct,QMMMOpts,0);
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
      //Print Optimized geometry
      Print_traj(Struct,outfile,QMMMOpts);
      //Check convergence
      optct += 1;
      OptDone = OptConverged(Struct,OldStruct,Forces,optct,QMMMOpts,0,0);
    }
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << endl;
  }
  //End of sections

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
    int sys = system(call.str().c_str());
  }
  if (PSI4 == 1)
  {
    //Clear any remaining PSI4 files
    stringstream call;
    call.copyfmt(cout);
    call.str("");
    call << "rm -f psi*";
    int sys = system(call.str().c_str());
  }
  if (SinglePoint == 1)
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
    int sys = system(call.str().c_str());
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
  cout << TotalHours << " hour";
  cout << '\n';
  cout << "  Wall time for QM Wrappers: ";
  cout << TotalQM << " hours";
  cout << '\n';
  cout << "  Wall time for MM Wrappers: ";
  cout << TotalMM << " hours";
  cout << '\n';
  cout << "  Wall time for FLUKE: ";
  cout << OtherTime << " hours";
  cout << '\n';
  if (Jokes == 1)
  {
    randnum = rand() % 10; //Randomly injure animals
    cout << "  Animals injured in the making of this science: ";
    cout << randnum;
    cout << '\n';
  }
  cout << "####################################################";
  cout << endl;
  //End of section

  //Print a quote
  if (Jokes == 1)
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
  cout << endl;
  return 0;
};
