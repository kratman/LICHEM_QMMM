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
 Kratz et al., J. Comput. Chem., 37, 11, 1019, (2016)

*/

//Primary LICHEM header
#include "LICHEM_headers.h"

int main(int argc, char* argv[])
{
  //Misc. initialization
  startTime = (unsigned)time(0); //Time the program starts
  srand((unsigned)time(0)); //Serial only random numbers
  //End of section

  //Output stream settings
  //NB: The streams should always be returned to these settings
  cout.precision(16);
  cerr.precision(16);
  //End of section

  //Initialize local variables
  string dummy; //Generic string
  double sumE,sumE2,denAvg,LxAvg,LyAvg,LzAvg,Ek; //Energies and properties
  fstream xyzFile,connectFile,regionFile,outFile; //Input and output files
  vector<QMMMAtom> QMMMData; //Atom list
  vector<QMMMAtom> OldQMMMData; //A copy of the atoms list
  vector<QMMMElec> Elecs; //Semi-classical electrons (eFF model)
  QMMMSettings QMMMOpts; //QM and MM wrapper settings
  int randNum; //Random integer
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
  ReadArgs(argc,argv,xyzFile,connectFile,regionFile,outFile);
  //End of section

  //Read input and check for errors
  ReadLICHEMInput(xyzFile,connectFile,regionFile,QMMMData,QMMMOpts);
  LICHEMErrorChecker(QMMMOpts);
  LICHEMPrintSettings(QMMMData,QMMMOpts);
  //End of section

  //Fix PBC
  if (PBCon)
  {
    //Relatively safe PBC correction
    if (!TINKER)
    {
      PBCCenter(QMMMData,QMMMOpts); //Center the atoms in the box
    }
  }
  //End of section

  //Create backup directories
  if (CheckFile("BACKUPQM"))
  {
    stringstream call;
    call.str("");
    //Delete old files
    call << "rm -rf " << QMMMOpts.backDir << "; ";
    //Create new directory
    call << "mkdir " << QMMMOpts.backDir;
    globalSys = system(call.str().c_str());
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
    double Eqm; //QM energy
    double Emm; //MM energy
    if (QMMMOpts.NBeads == 1)
    {
      cout << "Single-point energy:";
    }
    if (QMMMOpts.NBeads > 1)
    {
      cout << "Multi-point energies:";
    }
    cout << '\n' << '\n';
    cout.flush(); //Print progress
    //Loop over all beads
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Calculate QMMM energy
      Eqm = 0; //Reset QM energy
      Emm = 0; //Reset MM energy
      //Calculate QM energy
      if (QMMMOpts.NBeads > 1)
      {
        cout << " Energy for bead: " << p << '\n';
        cout.flush();
      }
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        Eqm += GaussianEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        Eqm += PSI4Energy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        Eqm += NWChemEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (QMMM or QMonly)
      {
        //Print QM partial energy
        cout << "  QM energy: " << LICHEMFormFloat(Eqm,16) << " eV";
        cout << '\n';
        //Print progress
        cout.flush();
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        Emm += TINKEREnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (AMBER)
      {
        int tStart = (unsigned)time(0);
        Emm += AMBEREnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        Emm += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      //Print the rest of the energies
      if (QMMM or MMonly)
      {
        //Print MM partial energy
        cout << "  MM energy: " << LICHEMFormFloat(Emm,16) << " eV";
        cout << '\n';
      }
      sumE = Eqm+Emm; //Total energy
      if (QMMM)
      {
        //Print total energy
        cout << "  QMMM energy: ";
        cout << LICHEMFormFloat(sumE,16) << " eV";
        cout << " ";
        cout << LICHEMFormFloat(sumE/har2eV,16) << " a.u.";
        cout << '\n';
      }
      cout << '\n';
      cout.flush(); //Print output
    }
  }
  //End of section

  //LICHEM frequency calculation
  else if (FreqCalc)
  {
    int remCt = 0; //Number of deleted translation and rotation modes
    int Ndof = 3*(Nqm+Npseudo); //Number of degrees of freedom
    MatrixXd QMMMHess(Ndof,Ndof);
    VectorXd QMMMFreqs(Ndof);
    if (QMMMOpts.NBeads == 1)
    {
      cout << "Single-point frequencies:";
    }
    if (QMMMOpts.NBeads > 1)
    {
      cout << "Multi-point frequencies:";
    }
    cout << '\n';
    cout.flush(); //Print progress
    //Loop over all beads
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Calculate QMMM frequencies
      QMMMHess.setZero(); //Reset frequencies
      QMMMFreqs.setZero();
      //Calculate QM energy
      if (QMMMOpts.NBeads > 1)
      {
        cout << '\n';
        cout << " Frequencies for bead: " << p << '\n';
        cout.flush();
      }
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += GaussianHessian(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += PSI4Hessian(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += NWChemHessian(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += TINKERHessian(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (AMBER)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += AMBERHessian(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += LAMMPSHessian(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      //Calculate frequencies
      QMMMFreqs = LICHEMFreq(QMMMData,QMMMHess,QMMMOpts,p,remCt);
      //Print the frequencies
      if (remCt > 0)
      {
        cout << "  | Identified " << remCt;
        cout << " translation/rotation modes";
        cout << '\n';
      }
      cout << "  | Frequencies:" << '\n' << '\n';
      cout << "   ";
      remCt = 0; //Reuse as a counter
      for (int i=0;i<Ndof;i++)
      {
        if (abs(QMMMFreqs(i)) > 0)
        {
          cout << " ";
          cout << LICHEMFormFloat(QMMMFreqs(i),10);
          remCt += 1;
          if (remCt == 3)
          {
            //Start a new line
            cout << '\n';
            cout << "   "; //Add extra space
            remCt = 0;
          }
        }
      }
      if (remCt != 0)
      {
        //Terminate trailing line
        cout << '\n';
      }
      cout << '\n';
      cout.flush(); //Print output
    }
  }
  //End of section

  //Optimize structure (native QM and MM package optimizers)
  else if (OptSim)
  {
    //NB: Currently only Gaussian works with this option
    VectorXd forces; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    cout << "Optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    sumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      sumE += GaussianEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      sumE += PSI4Energy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      sumE += NWChemEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      sumE += TINKEREnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      sumE += AMBEREnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      sumE += LAMMPSEnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    cout << " | Opt. step: ";
    cout << optCt << " | Energy: ";
    cout << LICHEMFormFloat(sumE,16) << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    //Run optimization
    bool optDone = 0;
    if (QMMMOpts.maxOptSteps == 0)
    {
      optDone = 1;
    }
    while (!optDone)
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Run MM optimization
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        sumE = TINKEROpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (AMBER)
      {
        int tStart = (unsigned)time(0);
        sumE = AMBEROpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE = LAMMPSOpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
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
        int tStart = (unsigned)time(0);
        sumE = GaussianOpt(QMMMData,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        sumE = PSI4Opt(QMMMData,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        sumE = NWChemOpt(QMMMData,QMMMOpts,0);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Print Optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      optDone = OptConverged(QMMMData,OldQMMMData,forces,optCt,QMMMOpts,0,0);
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
    VectorXd forces; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    cout << "Steepest descent optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    sumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      sumE += GaussianEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      sumE += PSI4Energy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.* ");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      sumE += NWChemEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      sumE += TINKEREnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      sumE += AMBEREnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      sumE += LAMMPSEnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    cout << " | Opt. step: ";
    cout << optCt << " | Energy: ";
    cout << LICHEMFormFloat(sumE,16) << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    //Run optimization
    bool optDone = 0;
    while (!optDone)
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Run MM optimization
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        sumE = TINKEROpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (AMBER)
      {
        int tStart = (unsigned)time(0);
        sumE = AMBEROpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE = LAMMPSOpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n';
        cout.flush();
      }
      cout << '\n';
      //Run QM optimization
      LICHEMSteepest(QMMMData,QMMMOpts,0);
      //Print Optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      optDone = OptConverged(QMMMData,OldQMMMData,forces,optCt,QMMMOpts,0,0);
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
    VectorXd forces; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Change optimization tolerance for the first step
    double savedQMOptTol = QMMMOpts.QMOptTol; //Save value from input
    double savedMMOptTol = QMMMOpts.MMOptTol; //Save value from input
    if (QMMMOpts.QMOptTol < 0.005)
    {
      QMMMOpts.QMOptTol = 0.005; //Speedy convergance on the first step
    }
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; //Speedy convergance on the first step
    }
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    cout << "DFP optimization:" << '\n';
    cout.flush(); //Print progress
    //Calculate initial energy
    sumE = 0; //Clear old energies
    //Calculate QM energy
    if (Gaussian)
    {
      int tStart = (unsigned)time(0);
      sumE += GaussianEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    if (PSI4)
    {
      int tStart = (unsigned)time(0);
      sumE += PSI4Energy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
      //Delete annoying useless files
      globalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      int tStart = (unsigned)time(0);
      sumE += NWChemEnergy(QMMMData,QMMMOpts,0);
      QMTime += (unsigned)time(0)-tStart;
    }
    //Calculate MM energy
    if (TINKER)
    {
      int tStart = (unsigned)time(0);
      sumE += TINKEREnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (AMBER)
    {
      int tStart = (unsigned)time(0);
      sumE += AMBEREnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    if (LAMMPS)
    {
      int tStart = (unsigned)time(0);
      sumE += LAMMPSEnergy(QMMMData,QMMMOpts,0);
      MMTime += (unsigned)time(0)-tStart;
    }
    cout << " | Opt. step: ";
    cout << optCt << " | Energy: ";
    cout << LICHEMFormFloat(sumE,16) << " eV";
    cout << '\n';
    cout.flush(); //Print progress
    //Run optimization
    bool optDone = 0;
    while (!optDone)
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Run MM optimization
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        sumE = TINKEROpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (AMBER)
      {
        int tStart = (unsigned)time(0);
        sumE = AMBEROpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE = LAMMPSOpt(QMMMData,QMMMOpts,0);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (QMMM)
      {
        cout << "    MM optimization complete.";
        cout << '\n';
        cout.flush();
      }
      cout << '\n';
      //Run QM optimization
      LICHEMDFP(QMMMData,QMMMOpts,0);
      //Reset tolerance before optimization check
      QMMMOpts.QMOptTol = savedQMOptTol;
      QMMMOpts.MMOptTol = savedMMOptTol;
      //Print Optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      optDone = OptConverged(QMMMData,OldQMMMData,forces,optCt,QMMMOpts,0,0);
      if (optCt == 1)
      {
        //Avoid terminating restarts on the loose tolerance step
        optDone = 0; //Not converged
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
    //Change units
    QMMMOpts.press *= atm2eV; //Pressure in eV/Ang^3
    //Adjust probabilities
    if (Natoms == 1)
    {
      //Remove atom centroid moves
      centProb = 0.0;
      beadProb = 1.0;
    }
    if (QMMMOpts.ensemble == "NVT")
    {
      //Remove volume changes
      volProb = 0.0;
    }
    //Initialize local variables
    sumE = 0; //Average energy
    sumE2 = 0; //Average squared energy
    denAvg = 0; //Average density
    LxAvg = 0; //Average box length
    LyAvg = 0; //Average box length
    LzAvg = 0; //Average box length
    Ek = 0; //PIMC kinietic energy
    if (QMMMOpts.NBeads > 1)
    {
      //Set kinetic energy
      Ek = 3*Natoms*QMMMOpts.NBeads/(2*QMMMOpts.beta);
    }
    int Nct = 0; //Step counter
    int ct = 0; //Secondary counter
    double Nacc = 0; //Number of accepted moves
    double Nrej = 0; //Number of rejected moves
    double Emc = 0; //Monte Carlo energy
    double Et = 0; //Total energy for printing
    bool acc; //Flag for accepting a step
    //Find the number of characters to print for the step counter
    int simCharLen;
    simCharLen = QMMMOpts.NEq+QMMMOpts.NSteps;
    simCharLen = LICHEMCount(simCharLen);
    //Start equilibration run and calculate initial energy
    cout << "Monte Carlo equilibration:" << '\n';
    cout.flush();
    QMMMOpts.EOld = 0;
    QMMMOpts.EOld += Get_PI_Epot(QMMMData,QMMMOpts);
    QMMMOpts.EOld += Get_PI_Espring(QMMMData,QMMMOpts);
    if (volProb > 0)
    {
      //Add PV term
      QMMMOpts.EOld += QMMMOpts.press*Lx*Ly*Lz;
    }
    Emc = QMMMOpts.EOld; //Needed if equilibration is skipped
    Nct = 0; //Reset counter to zero
    while (Nct < QMMMOpts.NEq)
    {
      Emc = 0;
      //Check step size
      if(ct == acc_Check)
      {
        if ((Nacc/(Nrej+Nacc)) > QMMMOpts.accRatio)
        {
          //Increase step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 1.001+randVal;
        }
        if ((Nacc/(Nrej+Nacc)) < QMMMOpts.accRatio)
        {
          //Decrease step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 0.999-randVal;
        }
        if (mcStep < stepMin)
        {
          //Set to minimum
          mcStep = stepMin;
        }
        if (mcStep > stepMax)
        {
          //Set to maximum
          mcStep = stepMax;
        }
        //Statistics
        cout << " | Step: " << setw(simCharLen) << Nct;
        cout << " | Step size: ";
        cout << LICHEMFormFloat(mcStep,6);
        cout << " | Accept ratio: ";
        cout << LICHEMFormFloat((Nacc/(Nrej+Nacc)),6);
        cout << '\n';
        cout.flush(); //Print stats
        //Reset counters
        ct = 0;
        Nacc = 0;
        Nrej = 0;
      }
      //Continue simulation
      ct += 1;
      acc = MCMove(QMMMData,QMMMOpts,Emc);
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
    Nct = 0; //Reset counter to zero
    Nacc = 0; //Reset counter to zero
    Nrej = 0; //Reset counter to zero
    cout << '\n';
    cout << "Monte Carlo production:" << '\n';
    cout.flush();
    //Print starting conditions
    Print_traj(QMMMData,outFile,QMMMOpts);
    Et = Ek+Emc; //Calculate total energy using previous saved energy
    Et -= 2*Get_PI_Espring(QMMMData,QMMMOpts);
    cout << " | Step: " << setw(simCharLen) << 0;
    cout << " | Energy: " << LICHEMFormFloat(Et,12);
    cout << " eV";
    if (QMMMOpts.ensemble == "NPT")
    {
      double rho;
      rho = LICHEMDensity(QMMMData,QMMMOpts);
      cout << " | Density: ";
      cout << LICHEMFormFloat(rho,8);
      cout << " g/cm\u00B3";
    }
    cout << '\n';
    cout.flush(); //Print results
    //Continue simulation
    while (Nct < QMMMOpts.NSteps)
    {
      Emc = 0; //Set energy to zero
      acc = MCMove(QMMMData,QMMMOpts,Emc);
      //Update averages
      Et = 0;
      Et += Ek+Emc;
      Et -= 2*Get_PI_Espring(QMMMData,QMMMOpts);
      denAvg += LICHEMDensity(QMMMData,QMMMOpts);
      LxAvg += Lx;
      LyAvg += Ly;
      LzAvg += Lz;
      sumE += Et;
      sumE2 += Et*Et;
      //Update counters and print output
      if (acc)
      {
        //Increase counters
        Nct += 1;
        Nacc += 1;
        //Print trajectory and instantaneous energies
        if ((Nct%QMMMOpts.NPrint) == 0)
        {
          //Print progress
          Print_traj(QMMMData,outFile,QMMMOpts);
          cout << " | Step: " << setw(simCharLen) << Nct;
          cout << " | Energy: " << LICHEMFormFloat(Et,12);
          cout << " eV";
          if (QMMMOpts.ensemble == "NPT")
          {
            double rho;
            rho = LICHEMDensity(QMMMData,QMMMOpts);
            cout << " | Density: ";
            cout << LICHEMFormFloat(rho,8);
            cout << " g/cm\u00B3";
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
    if ((Nct%QMMMOpts.NPrint) != 0)
    {
      //Print final geometry if it was not already written
      Print_traj(QMMMData,outFile,QMMMOpts);
    }
    sumE /= Nrej+Nacc; //Average energy
    sumE2 /= Nrej+Nacc; //Variance of the energy
    denAvg /= Nrej+Nacc; //Average density
    LxAvg /= Nrej+Nacc; //Average box size
    LyAvg /= Nrej+Nacc; //Average box size
    LzAvg /= Nrej+Nacc; //Average box size
    //Print simulation details and statistics
    cout << '\n';
    if (QMMMOpts.NBeads > 1)
    {
      cout << "PI";
    }
    cout << "MC statistics:" << '\n';
    if (QMMMOpts.ensemble == "NPT")
    {
      //Print simulation box information
      cout << " | Density: ";
      cout << LICHEMFormFloat(denAvg,8);
      cout << " g/cm\u00B3" << '\n';
      cout << " | Average box size (\u212B): " << '\n';
      cout << "  "; //Indent
      cout << " Lx = " << LICHEMFormFloat(LxAvg,12);
      cout << " Ly = " << LICHEMFormFloat(LyAvg,12);
      cout << " Lz = " << LICHEMFormFloat(LzAvg,12);
      cout << '\n';
    }
    cout << " | Average energy: ";
    cout << LICHEMFormFloat(sumE,16);
    cout << " eV | Variance: ";
    cout << LICHEMFormFloat((sumE2-(sumE*sumE)),12);
    cout << " eV\u00B2";
    cout << '\n';
    cout << " | Acceptance ratio: ";
    cout << LICHEMFormFloat((Nacc/(Nrej+Nacc)),6);
    cout << " | Optimized step size: ";
    cout << LICHEMFormFloat(mcStep,6);
    cout << " \u212B";
    cout << '\n';
    cout << '\n';
    cout.flush();
  }
  //End of section

  //Force-bias NEB Monte Carlo
  else if (FBNEBSim)
  {
    //Initialize local variables
    int Nct = 0; //Step counter
    int ct = 0; //Secondary counter
    double Nacc = 0; //Number of accepted moves
    double Nrej = 0; //Number of rejected moves
    int acc; //Number of steps accepted along the path
    vector<VectorXd> allForces; //Stores forces between MC steps
    VectorXd sumE(QMMMOpts.NBeads); //Average energy array
    VectorXd sumE2(QMMMOpts.NBeads); //Average squared energy array
    VectorXd Emc(QMMMOpts.NBeads); //Current MC energy
    sumE.setZero();
    sumE2.setZero();
    Emc.setZero();
    //Initialize force arrays
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Zero force vectors make the first move an energy calculation
      VectorXd tmp(3*Natoms);
      tmp.setZero();
      allForces.push_back(tmp);
    }
    //Find the number of characters to print for the step counter
    int simCharLen;
    simCharLen = QMMMOpts.NEq+QMMMOpts.NSteps;
    simCharLen = LICHEMCount(simCharLen);
    //Start equilibration run
    cout << "Monte Carlo equilibration:" << '\n';
    cout.flush();
    int savedNPrint = QMMMOpts.NPrint;
    if (QMMMOpts.NPrint < 100)
    {
      //Prevent the print rate from breaking the tuning
      QMMMOpts.NPrint = 100; //Minimum value
    }
    QMMMOpts.EOld = hugeNum; //Forces the first step to be accepted
    Nct = 0; //Reset counter to zero
    while (Nct < QMMMOpts.NEq)
    {
      //Check step size
      if (ct == QMMMOpts.NPrint)
      {
        if ((Nacc/(Nrej+Nacc)) > QMMMOpts.accRatio)
        {
          //Increase step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 1.001+randVal;
        }
        if ((Nacc/(Nrej+Nacc)) < QMMMOpts.accRatio)
        {
          //Decrease step size
          double randVal; //Use random values to keep from cycling up and down
          randVal = (((double)rand())/((double)RAND_MAX));
          randVal /= 10.0;
          mcStep *= 0.999-randVal;
        }
        if (mcStep < stepMin)
        {
          //Set to minimum
          mcStep = stepMin;
        }
        if (mcStep > stepMax)
        {
          //Set to maximum
          mcStep = stepMax;
        }
        //Statistics
        cout << " | Accepted: " << setw(simCharLen) << Nct;
        cout << " | Step size: ";
        cout << LICHEMFormFloat(mcStep,6);
        cout << " | Accept ratio: ";
        cout << LICHEMFormFloat((Nacc/(Nrej+Nacc)),6);
        cout << '\n';
        cout.flush(); //Print stats
        //Reset counters for the next round of tuning
        ct = 0;
        Nacc = 0;
        Nrej = 0;
      }
      //Continue simulation
      ct += 1; //Equilibration counts cycles instead of steps
      acc = FBNEBMCMove(QMMMData,allForces,QMMMOpts,Emc);
      Nct += acc; //Equilibration counts acceptances instead of steps
      Nacc += acc;
      Nrej += QMMMOpts.NBeads-acc;
    }
    QMMMOpts.NPrint = savedNPrint; //Restore user defined sample rate
    cout << " Equilibration complete." << '\n';
    //Start production run
    Nct = 0; //Reset counter to zero
    Nacc = 0; //Reset counter to zero
    Nrej = 0; //Reset counter to zero
    cout << '\n';
    cout << "Monte Carlo production:" << '\n';
    cout.flush();
    //Print starting conditions
    Print_traj(QMMMData,outFile,QMMMOpts);
    cout << " | Steps: " << setw(simCharLen) << 0;
    cout << " | Accepted: " << setw(simCharLen) << 0;
    cout << '\n';
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      cout << "    Bead: ";
      cout << setw(3) << p << " | Energy: ";
      cout << LICHEMFormFloat(Emc(p),16) << " eV" << '\n';
    }
    cout.flush(); //Print results
    //Continue simulation
    while (Nacc < QMMMOpts.NSteps)
    {
      acc = FBNEBMCMove(QMMMData,allForces,QMMMOpts,Emc);
      //Update statistics
      #pragma omp parallel for schedule(dynamic)
      for (int p=0;p<QMMMOpts.NBeads;p++)
      {
        sumE(p) += Emc(p);
        sumE2(p) += Emc(p)*Emc(p);
      }
      //Update counters
      Nct += QMMMOpts.NBeads;
      Nacc += acc;
      Nrej += QMMMOpts.NBeads;
      //Print output
      if ((((Nct/QMMMOpts.NBeads)%QMMMOpts.NPrint) == 0) or
         (Nacc == QMMMOpts.NSteps))
      {
        //Print progress
        Print_traj(QMMMData,outFile,QMMMOpts);
        cout << " | Steps: " << setw(simCharLen) << Nct;
        cout << " | Accepted: " << setw(simCharLen) << Nacc;
        cout << '\n';
        for (int p=0;p<QMMMOpts.NBeads;p++)
        {
          cout << "    Bead: ";
          cout << setw(3) << p << " | Energy: ";
          cout << LICHEMFormFloat(Emc(p),16) << " eV" << '\n';
        }
        cout.flush(); //Print results
      }
    }
    if (((Nct/QMMMOpts.NBeads)%QMMMOpts.NPrint) != 0)
    {
      //Print final geometry if it was not already written
      Print_traj(QMMMData,outFile,QMMMOpts);
    }
    //Calculate statistics
    double scaleBy = 0.0; //Avoid integer division
    scaleBy += QMMMOpts.NBeads; //Adjust for the number of replicas
    scaleBy /= Nct; //Adjust for the total number of samples
    sumE *= scaleBy; //Scale the sum to calculate the average
    sumE2 *= scaleBy; //Scale the sum to calculate the average
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      //Standard deviation
      sumE2(p) = sqrt(sumE2(p)-sumE(p)*sumE(p));
    }
    //Print simulation details and statistics
    cout << '\n';
    cout << "Monte Carlo statistics:" << '\n';
    cout << " | Acceptance ratio: ";
    cout << LICHEMFormFloat((Nacc/Nct),6);
    cout << " | Optimized step size: ";
    cout << LICHEMFormFloat(mcStep,6);
    cout << " \u212B";
    cout << '\n';
    cout << " | Average energies:" << '\n';
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      cout << "    Bead: ";
      cout << setw(3) << p << " | Energy: ";
      cout << LICHEMFormFloat(sumE(p),16);
      cout << " +/- " << LICHEMFormFloat(sumE2(p),16);
      cout << " eV" << '\n';
    }
    cout << '\n';
    cout.flush();
  }
  //End of section

  //NEB optimization
  else if (NEBSim)
  {
    MatrixXd forceStats; //Dummy array needed for convergence tests
    int optCt = 0; //Counter for optimization steps
    //Check number of beads for Climbing image
    if (QMMMOpts.NBeads < 4)
    {
      //If the system is well behaved, then the TS bead is known.
      QMMMOpts.climb = 1; //Turn on climbing image NEB
    }
    //Change optimization tolerance for the first step
    double savedQMOptTol = QMMMOpts.QMOptTol; //Save value from input
    double savedMMOptTol = QMMMOpts.MMOptTol; //Save value from input
    if (QMMMOpts.QMOptTol < 0.005)
    {
      QMMMOpts.QMOptTol = 0.005; //Speedy convergance on the first step
    }
    if (QMMMOpts.MMOptTol < 0.25)
    {
      QMMMOpts.MMOptTol = 0.25; //Speedy convergance on the first step
    }
    //Print initial structure
    Print_traj(QMMMData,outFile,QMMMOpts);
    cout << "Nudged elastic band optimization:" << '\n';
    if (QMMMOpts.climb)
    {
      cout << " | Short path detected. Starting climbing image NEB.";
      cout << '\n' << '\n';
    }
    cout << " | Opt. step: 0 | Bead energies:";
    cout << '\n';
    cout.flush(); //Print progress
    //Calculate reaction coordinate positions
    VectorXd reactCoord(QMMMOpts.NBeads); //Reaction coordinate
    reactCoord.setZero();
    for (int p=0;p<(QMMMOpts.NBeads-1);p++)
    {
      MatrixXd geom1((Nqm+Npseudo),3); //Current replica
      MatrixXd geom2((Nqm+Npseudo),3); //Next replica
      VectorXd disp; //Store the displacement
      //Save geometries
      int ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
        {
          //Save current replica
          geom1(ct,0) = QMMMData[i].P[p].x;
          geom1(ct,1) = QMMMData[i].P[p].y;
          geom1(ct,2) = QMMMData[i].P[p].z;
          //Save replica p+1
          geom2(ct,0) = QMMMData[i].P[p+1].x;
          geom2(ct,1) = QMMMData[i].P[p+1].y;
          geom2(ct,2) = QMMMData[i].P[p+1].z;
          ct += 1;
        }
      }
      //Calculate displacement
      disp = KabschDisplacement(geom1,geom2,(Nqm+Npseudo));
      //Remove inactive atoms
      ct = 0; //Reset counter for the number of atoms
      for (int i=0;i<Natoms;i++)
      {
        //Only include QM and PB regions
        if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
        {
          //Only include active atoms in the tangent
          if (!QMMMData[i].NEBActive)
          {
            //Delete distance components
            disp(ct) = 0;
            disp(ct+1) = 0;
            disp(ct+2) = 0;
          }
          //Advance counter
          ct += 3;
        }
      }
      //Update reaction coordinate
      reactCoord(p+1) = reactCoord(p); //Start from previous bead
      reactCoord(p+1) += disp.norm(); //Add magnitude of the displacement
    }
    reactCoord /= reactCoord.maxCoeff(); //Must be between 0 and 1
    //Calculate initial energies
    QMMMOpts.ETrans = -1*hugeNum; //Locate the initial transition state
    for (int p=0;p<QMMMOpts.NBeads;p++)
    {
      sumE = 0; //Clear old energies
      //Calculate QM energy
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        sumE += GaussianEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        sumE += PSI4Energy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        sumE += NWChemEnergy(QMMMData,QMMMOpts,p);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Calculate MM energy
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        sumE += TINKEREnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (AMBER)
      {
        int tStart = (unsigned)time(0);
        sumE += AMBEREnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        sumE += LAMMPSEnergy(QMMMData,QMMMOpts,p);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (p == 0)
      {
        //Save reactant energy
        QMMMOpts.EReact = sumE;
      }
      else if (p == (QMMMOpts.NBeads-1))
      {
        //Save product energy
        QMMMOpts.EProd = sumE;
      }
      cout << "   Bead: ";
      cout << setw(LICHEMCount(QMMMOpts.NBeads)) << p;
      cout << " | React. coord: ";
      cout << LICHEMFormFloat(reactCoord(p),5);
      cout << " | Energy: ";
      cout << LICHEMFormFloat(sumE,16) << " eV";
      cout << '\n';
      cout.flush(); //Print progress
      //Update transition state
      if (sumE > QMMMOpts.ETrans)
      {
        //Save new properties
        QMMMOpts.TSBead = p;
        QMMMOpts.ETrans = sumE;
      }
      //Copy checkpoint data to speed up first step
      if ((p != (QMMMOpts.NBeads-1)) and QMMMOpts.startPathChk)
      {
        stringstream call;
        if (Gaussian and (QMMMOpts.func != "SemiEmp"))
        {
          call.str("");
          call << "cp LICHM_" << p << ".chk ";
          call << "LICHM_" << (p+1) << ".chk";
          call << " 2> LICHM_" << (p+1) << ".trash; ";
          call << "rm -f LICHM_" << (p+1) << ".trash";
          globalSys = system(call.str().c_str());
        }
        if (PSI4)
        {
          call.str("");
          call << "cp LICHM_" << p << ".180 ";
          call << "LICHM_" << (p+1) << ".180";
          call << " 2> LICHM_" << (p+1) << ".trash; ";
          call << "rm -f LICHM_" << (p+1) << ".trash";
          globalSys = system(call.str().c_str());
        }
      }
    }
    //Run optimization
    bool pathDone = 0;
    int pathStart = 0; //First bead to optimize
    int pathEnd = QMMMOpts.NBeads; //Last bead to optimize
    if (QMMMOpts.frznEnds)
    {
      //Change the start and end points
      pathStart = 1;
      pathEnd = QMMMOpts.NBeads-1;
    }
    while (!pathDone)
    {
      //Copy structure
      OldQMMMData = QMMMData;
      //Run MM optimization
      for (int p=pathStart;p<pathEnd;p++)
      {
        if (TINKER)
        {
          int tStart = (unsigned)time(0);
          sumE = TINKEROpt(QMMMData,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tStart;
        }
        if (AMBER)
        {
          int tStart = (unsigned)time(0);
          sumE = AMBEROpt(QMMMData,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tStart;
        }
        if (LAMMPS)
        {
          int tStart = (unsigned)time(0);
          sumE = LAMMPSOpt(QMMMData,QMMMOpts,p);
          MMTime += (unsigned)time(0)-tStart;
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
      LICHEMNEB(QMMMData,QMMMOpts,optCt);
      //Reset tolerance before optimization check
      QMMMOpts.QMOptTol = savedQMOptTol;
      QMMMOpts.MMOptTol = savedMMOptTol;
      //Print optimized geometry
      Print_traj(QMMMData,outFile,QMMMOpts);
      //Check convergence
      optCt += 1;
      pathDone = PathConverged(QMMMData,OldQMMMData,forceStats,optCt,
                               QMMMOpts,0);
      if (optCt == 1)
      {
        //Avoid terminating restarts on the loose tolerance step
        pathDone = 0; //Not converged
      }
    }
    BurstTraj(QMMMData,QMMMOpts);
    cout << '\n';
    cout << "Optimization complete.";
    cout << '\n' << '\n';
    //Print the reaction barriers
    double dEfor = QMMMOpts.ETrans-QMMMOpts.EReact; //Forward barrier
    double dErev = QMMMOpts.ETrans-QMMMOpts.EProd; //Reverse barrier
    cout << "NEB Results:" << '\n';
    cout << " | Transition state bead: " << QMMMOpts.TSBead;
    cout << '\n' << '\n';
    cout << " | Forward barrier: ";
    cout << LICHEMFormFloat(dEfor,16) << " eV ," << '\n';
    cout << " |   ";
    cout << LICHEMFormFloat(dEfor/har2eV,16) << " a.u.";
    cout << " , ";
    cout << LICHEMFormFloat(dEfor/kcal2eV,16) << " kcal/mol";
    cout << '\n' << '\n';
    cout << " | Reverse barrier: ";
    cout << LICHEMFormFloat(dErev,16) << " eV ," << '\n';
    cout << " |   ";
    cout << LICHEMFormFloat(dErev/har2eV,16) << " a.u.";
    cout << " , ";
    cout << LICHEMFormFloat(dErev/kcal2eV,16) << " kcal/mol";
    cout << '\n' << '\n';
    //Check stability
    if (QMMMOpts.NEBFreq)
    {
      //Calculate TS frequencies
      int remCt = 0; //Number of deleted translation and rotation modes
      int Ndof = 3*(Nqm+Npseudo); //Number of degrees of freedom
      MatrixXd QMMMHess(Ndof,Ndof);
      VectorXd QMMMFreqs(Ndof);
      cout << '\n'; //Print blank line
      cout << "TS frequencies:";
      cout << '\n';
      cout.flush(); //Print progress
      //Calculate QMMM frequencies
      QMMMHess.setZero(); //Reset Hessian
      QMMMFreqs.setZero(); //Reset frequencies
      //Calculate QM Hessian
      if (Gaussian)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += GaussianHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        QMTime += (unsigned)time(0)-tStart;
      }
      if (PSI4)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += PSI4Hessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        QMTime += (unsigned)time(0)-tStart;
        //Delete annoying useless files
        globalSys = system("rm -f psi.* timer.*");
      }
      if (NWChem)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += NWChemHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        QMTime += (unsigned)time(0)-tStart;
      }
      //Calculate MM Hessian
      if (TINKER)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += TINKERHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (AMBER)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += AMBERHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        MMTime += (unsigned)time(0)-tStart;
      }
      if (LAMMPS)
      {
        int tStart = (unsigned)time(0);
        QMMMHess += LAMMPSHessian(QMMMData,QMMMOpts,QMMMOpts.TSBead);
        MMTime += (unsigned)time(0)-tStart;
      }
      //Calculate frequencies
      QMMMFreqs = LICHEMFreq(QMMMData,QMMMHess,QMMMOpts,QMMMOpts.TSBead,remCt);
      //Print the frequencies
      if (remCt > 0)
      {
        cout << "  | Identified " << remCt;
        cout << " translation/rotation modes";
        cout << '\n';
      }
      cout << "  | Frequencies:" << '\n' << '\n';
      cout << "   ";
      remCt = 0; //Reuse as a counter
      for (int i=0;i<Ndof;i++)
      {
        if (abs(QMMMFreqs(i)) > 0)
        {
          cout << " ";
          cout << LICHEMFormFloat(QMMMFreqs(i),10);
          remCt += 1;
          if (remCt == 3)
          {
            //Start a new line
            cout << '\n';
            cout << "   "; //Add extra space
            remCt = 0;
          }
        }
      }
      if (remCt != 0)
      {
        //Terminate trailing line
        cout << '\n';
      }
      cout << '\n';
    }
    cout.flush();
  }
  //End of section

  //Inform the user if no simulations were performed
  else
  {
    cout << "Nothing was done..." << '\n';
    cout << "Check the simulation type in " << regFilename;
    cout << '\n' << '\n';
    cout.flush();
  }
  //End of section

  //Clean up
  if (Gaussian)
  {
    //Clear any remaining Gaussian files
    stringstream call; //Stream for system calls and reading/writing files
    call.str("");
    call << "rm -f Gau-*"; //Produced if there is a crash
    globalSys = system(call.str().c_str());
  }
  if (PSI4)
  {
    //Clear any remaining PSI4 files
    stringstream call; //Stream for system calls and reading/writing files
    call.str("");
    call << "rm -f psi*";
    globalSys = system(call.str().c_str());
  }
  if (SinglePoint or FreqCalc)
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
    globalSys = system(call.str().c_str());
  }
  //End of section

  //Print usage statistics
  endTime = (unsigned)time(0); //Time the program completes
  double totalHours = (double(endTime)-double(startTime));
  double totalQM = double(QMTime);
  if ((QMMMOpts.NBeads > 1) and (PIMCSim or FBNEBSim))
  {
    //Average over the number of running simulations
    totalQM /= Nthreads;
  }
  double totalMM = double(MMTime);
  if ((QMMMOpts.NBeads > 1) and (PIMCSim or FBNEBSim))
  {
    //Average over the number of running simulations
    totalMM /= Nthreads;
  }
  double otherTime = totalHours-totalQM-totalMM;
  totalHours /= 3600.0; //Convert from seconds to hours
  totalQM /= 3600.0; //Convert from seconds to hours
  totalMM /= 3600.0; //Convert from seconds to hours
  otherTime /= 3600.0; //Convert from seconds to hours
  cout << "################# Usage Statistics #################";
  cout << '\n';
  cout << "  Total wall time:                     ";
  cout << LICHEMFormFloat(totalHours,6) << " hours";
  cout << '\n';
  cout << "  Wall time for QM Wrappers:           ";
  cout << LICHEMFormFloat(totalQM,6) << " hours";
  cout << '\n';
  cout << "  Wall time for MM Wrappers:           ";
  cout << LICHEMFormFloat(totalMM,6) << " hours";
  cout << '\n';
  cout << "  Wall time for LICHEM:                ";
  cout << LICHEMFormFloat(otherTime,6) << " hours";
  cout << '\n';
  cout << "####################################################";
  cout << '\n';
  cout.flush();
  //End of section

  //Print a quote
  if (JOKES)
  {
    cout << '\n';
    cout << "Random quote:";
    cout << '\n';
    string quote; //Random quote
    vector<string> Quotes; //Stores all possible quotes
    FetchQuotes(Quotes); //Fetch list of quotes
    randNum = rand() % 1000; //Randomly pick 1 of 1000 quotes
    cout << Quotes[randNum]; //Print quote
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
  int retValue = globalSys;
  retValue = 0; //This can be changed to error messages later
  //End of section

  //Quit
  return retValue;
};

