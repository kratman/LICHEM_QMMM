/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 Routines for reading and checking the input for FLUKE.

*/

//Misc. functions
void ReadArgs(int& argc, char**& argv, fstream& xyzfile,
     fstream& connectfile, fstream& regionfile, fstream& outfile)
{
  int sys; //Dummy return for system calls
  string dummy; //Generic string
  stringstream call;
  //Read command line arguments
  if (argc == 1)
  {
    //Escape if there are no arguments
    cout << '\n';
    cout << "Missing arguments..." << '\n' << '\n';
    cout << "Usage: FLUKE -n Ncpus -x Input.xyz -c Connectivity.inp ";
    cout << "-r Regions.inp -o Output.xyz" << '\n';
    cout << '\n';
    cout << "Use -h or --help for detailed instructions.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  dummy = string(argv[1]);
  if (dummy == "-GauExtern")
  {
    //Escape to GauExternal
    ExternalGaussian(argc,argv);
  }
  if (dummy == "-convert")
  {
    //Attempt to create FLUKE input from other formats
    dummy = string(argv[2]);
    if (dummy == "-t")
    {
      TINK2FLUKE(argc,argv);
    }
    else
    {
      cout << '\n';
      cout << "Unrecognized file format.";
      cout << '\n';
      cout << '\n';
      cout.flush();
    }
  }
  if (dummy == "-tinker")
  {
    //Attempt to create a TINKER XYZ file from FLUKE input
    FLUKE2TINK(argc,argv);
  }
  if (dummy == "-GlobalPoles")
  {
    ExtractGlobalPoles(argc,argv);
  }
  if ((argc % 2) != 1)
  {
    dummy = string(argv[1]);
    if ((dummy != "-h") and (dummy != "--help"))
    {
      //Escape if there are missing arguments
      cout << '\n';
      cout << "Odd number of arguments..." << '\n' << '\n';
      cout << "Usage: FLUKE -n Ncpus -x Input.xyz -c Connectivity.inp ";
      cout << "-r Regions.inp -o Output.xyz" << '\n';
      cout << '\n';
      cout << "Use -h or --help for detailed instructions.";
      cout << '\n';
      cout.flush();
      exit(0);
    }
  }
  for (int i=0;i<argc;i++)
  {
    //Read file names and CPUs
    dummy = string(argv[i]);
    if ((dummy == "-h") or (dummy == "--help"))
    {
      //Print helpful information and exit
      cout << '\n';
      cout << "Usage: FLUKE -n Ncpus -x Input.xyz -c Connectivity.inp ";
      cout << "-r Regions.inp -o Output.xyz" << '\n';
      cout << '\n';
      cout << "Command line arguments:" << '\n' << '\n';
      cout << "  -n    Number of CPUs used for the QM calculation." << '\n';
      cout << '\n';
      cout << "  -x    Input xyz file." << '\n' << '\n';
      cout << "  -c    Connectivity and force field input file." << '\n';
      cout << '\n';
      cout << "  -r    Information about how the system is subdivided" << '\n';
      cout << "        into QM, MM, and psuedo-atom regions." << '\n' << '\n';
      cout << "  -o    Output xyz file for the optimized structures.";
      cout << '\n' << '\n';
      cout.flush();
      exit(0);
    }
    if (dummy == "-n")
    {
      Ncpus = atoi(argv[i+1]);
      call.str("");
      call << "export OMP_NUM_THREADS=" << Ncpus;
      sys = system(call.str().c_str());
    }
    if (dummy == "-x")
    {
      xyzfilename = string(argv[i+1]);
      xyzfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-c")
    {
      confilename = string(argv[i+1]);
      connectfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-r")
    {
      regfilename = string(argv[i+1]);
      regionfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-o")
    {
      outfile.open(argv[i+1],ios_base::out);
    }
  }
  for (int i=0;i<argc;i++)
  {
    //Detect bad arguments
    dummy = string(argv[i]);
    if (dummy[0] == '-')
    {
      bool BadArgs = 0; //Bad argument found
      if ((dummy != "-n") and (dummy != "-x") and
      (dummy != "-c") and (dummy != "-r") and
      (dummy != "-o"))
      {
        BadArgs = 1;
      }
      if (BadArgs == 1)
      {
        cout << '\n';
        cout << "Unrecognized flag..." << '\n' << '\n';
        cout << "Usage: FLUKE -n Ncpus -x Input.xyz -c Connectivity.inp ";
        cout << "-r Regions.inp -o Output.xyz" << '\n';
        cout << '\n';
        cout << "Use -h or --help for detailed instructions.";
        cout << '\n';
        cout.flush();
        exit(0);
      }
    }
  }
  if (argc != 11)
  {
    //Escape if there are too few arguments
    cout << '\n';
    cout << "Missing arguments..." << '\n' << '\n';
    cout << "Usage: FLUKE -n Ncpus -x Input.xyz -c Connectivity.inp ";
    cout << "-r Regions.inp -o Output.xyz" << '\n';
    cout << '\n';
    cout << "Use -h or --help for detailed instructions.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  //Make sure input files can be read
  bool DoQuit = 0;
  if (!xyzfile.good())
  {
    cout << "Error: Could not open xyz file.";
    cout << '\n';
    cout.flush();
    DoQuit = 1;
  }
  if (!connectfile.good())
  {
    cout << "Error: Could not open connectivity file.";
    cout << '\n';
    cout.flush();
    DoQuit = 1;
  }
  if (!regionfile.good())
  {
    cout << "Error: Could not open region file.";
    cout << '\n';
    cout.flush();
    DoQuit = 1;
  }
  if (!outfile.good())
  {
    cout << "Error: Could not create output file.";
    cout << '\n';
    cout.flush();
    DoQuit = 1;
  }
  if (DoQuit == 1)
  {
    exit(0);
  }
  return;
};

void ReadFLUKEInput(fstream& xyzfile, fstream& connectfile,
     fstream& regionfile, vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Read input
  string dummy; //Generic string
  if (GauExternal == 0)
  {
    xyzfile >> Natoms;
    for (int i=0;i<Natoms;i++)
    {
      //Save atom information
      QMMMAtom tmp;
      xyzfile >> tmp.QMTyp;
      Coord tmp2;
      xyzfile >> tmp2.x >> tmp2.y >> tmp2.z;
      tmp.P.push_back(tmp2); //Set up zeroth replica
      tmp.id = i;
      tmp.QMregion = 0;
      tmp.MMregion = 1;
      tmp.PAregion = 0;
      tmp.BAregion = 0;
      tmp.Frozen = 0;
      Mpole tmp3; //Initialize charges and multipoles
      OctCharges tmp4; //Initialize charges and multipoles
      tmp.MP.push_back(tmp3);
      tmp.PC.push_back(tmp4);
      Struct.push_back(tmp);
    }
  }
  for (int i=0;i<Natoms;i++)
  {
    //Save connectivity information
    int tmp;
    //id MMTyp NumTyp q Nbonds [connectivity]
    connectfile >> tmp; //Atom ID
    if (tmp != Struct[i].id)
    {
      //Escape if connectivity errors are found
      cout << "Error: Atoms in the connectivity file are out of order.";
      cout << '\n';
      cout.flush();
      exit(0); //Escape
    }
    connectfile >> Struct[i].MMTyp >> Struct[i].NumTyp;
    connectfile >> Struct[i].m >> Struct[i].MP[0].q;
    connectfile >> tmp; //Number of bonds
    for (int j=0;j<tmp;j++)
    {
      //Save each bond to the atom's connectivity list
      int AtomID;
      connectfile >> AtomID;
      if (AtomID >= Natoms)
      {
        //Search for more connectivity errors
        cout << "Error: Atom index out of range in connectivity.";
        cout << '\n';
        cout << "Atom " << i << " bonded to non-existant atom ";
        cout << AtomID << '\n';
        cout.flush();
        exit(0); //Escape
      }
      Struct[i].Bonds.push_back(AtomID); //Add bond
    }
  }
  //Collect misc. simulation options
  regionfile >> dummy >> dummy; //Potential type
  if ((dummy == "QM") or (dummy == "qm"))
  {
    //Read only QM options
    QMonly = 1;
    MMonly = 0;
    QMMM = 0;
    regionfile >> dummy >> dummy;
    if ((dummy == "psi4") or (dummy == "Psi4")
       or (dummy == "PSI4"))
    {
      PSI4 = 1;
    }
    if ((dummy == "gaussian") or (dummy == "Gaussian"))
    {
      Gaussian = 1;
    }
    regionfile >> dummy >> QMMMOpts.Func;
    regionfile >> dummy >> QMMMOpts.Basis;
    regionfile >> dummy >> QMMMOpts.RAM;
    regionfile >> dummy >> QMMMOpts.Charge;
    regionfile >> dummy >> QMMMOpts.Spin;
    //Place all atoms in the QM region
    for (int i=0;i<Natoms;i++)
    {
      Struct[i].QMregion = 1;
      Struct[i].MMregion = 0;
      Struct[i].PAregion = 0;
      Struct[i].BAregion = 0;
    }
  }
  if ((dummy == "QMMM") or (dummy == "qmmm"))
  {
    //Read QM and MM options
    QMMM = 1;
    QMonly = 0;
    MMonly = 0;
    regionfile >> dummy >> dummy;
    if ((dummy == "psi4") or (dummy == "Psi4")
       or (dummy == "PSI4"))
    {
      PSI4 = 1;
    }
    if ((dummy == "gaussian") or (dummy == "Gaussian"))
    {
      Gaussian = 1;
    }
    regionfile >> dummy >> QMMMOpts.Func;
    regionfile >> dummy >> QMMMOpts.Basis;
    regionfile >> dummy >> QMMMOpts.RAM;
    regionfile >> dummy >> QMMMOpts.Charge;
    regionfile >> dummy >> QMMMOpts.Spin;
    regionfile >> dummy >> dummy; //MM wrapper
    if ((dummy == "Tinker") or (dummy == "TINKER")
       or (dummy == "tinker"))
    {
      TINKER = 1;
    }
    if ((dummy == "AMBER") or (dummy == "Amber") or
       (dummy == "amber"))
    {
      AMBER = 1;
    }
    if ((dummy == "LAMMPS") or (dummy == "lammps") or
       (dummy == "Lammps"))
    {
      LAMMPS = 1;
    }
    regionfile >> dummy >> dummy; //Potential type
    if ((dummy == "AMOEBA") or (dummy == "amoeba"))
    {
      //AMOEBA polarizable force field
      AMOEBA = 1;
      regionfile >> dummy >> QMMMOpts.Nind;
      if (TINKER == 1)
      {
        ExtractTINKpoles(Struct,0);
      }
    }
    if ((dummy == "Charges") or (dummy == "charges") or
       (dummy == "Charge") or (dummy == "charge") or
       (dummy == "point-charge"))
    {
      //Point-charge force fields
      CHRG = 1;
    }
    if ((dummy == "GEM") or (dummy == "gem") or (dummy == "Gem"))
    {
      //Frozen density
      GEM = 1;
      regionfile >> dummy >> QMMMOpts.Nind;
    }
  }
  if ((dummy == "MM") or (dummy == "mm"))
  {
    //Read only MM options
    MMonly = 1;
    QMonly = 0;
    QMMM = 0;
    regionfile >> dummy >> dummy; //MM wrapper
    if ((dummy == "Tinker") or (dummy == "TINKER")
       or (dummy == "tinker"))
    {
      TINKER = 1;
    }
    if ((dummy == "AMBER") or (dummy == "Amber") or
       (dummy == "amber"))
    {
      AMBER = 1;
    }
    if ((dummy == "LAMMPS") or (dummy == "lammps") or
       (dummy == "Lammps"))
    {
      LAMMPS = 1;
    }
    regionfile >> dummy >> dummy; //Potential type
    if ((dummy == "AMOEBA") or (dummy == "amoeba"))
    {
      //AMOEBA polarizable force field
      AMOEBA = 1;
      if (TINKER == 1)
      {
        ExtractTINKpoles(Struct,0);
      }
    }
    if ((dummy == "Charges") or (dummy == "charges") or
       (dummy == "Charge") or (dummy == "charge") or
       (dummy == "point-charge"))
    {
      //Point-charge force fields
      CHRG = 1;
    }
    if ((dummy == "GEM") or (dummy == "gem") or (dummy == "Gem"))
    {
      //Frozen density
      GEM = 1;
    }
  }
  regionfile >> dummy >> dummy; //Calculation type
  if ((dummy == "PIMC") or (dummy == "pimc"))
  {
    //Read MC and PIMC options
    MDSim = 0;
    OptSim = 0;
    DFPSim = 0;
    PIMCSim = 1;
    SteepSim = 0;
    SinglePoint = 0;
    regionfile >> dummy >> dummy; //Ensemble
    if ((dummy == "NVT") or (dummy == "nvt"))
    {
      QMMMOpts.Ensemble = "NVT";
    }
    if ((dummy == "NPT") or (dummy == "npt"))
    {
      QMMMOpts.Ensemble = "NPT";
    }
    regionfile >> dummy >> QMMMOpts.Temp;
    QMMMOpts.Beta = 1/(k*QMMMOpts.Temp);
    regionfile >> dummy >> QMMMOpts.Press;
    regionfile >> dummy >> QMMMOpts.Neq;
    regionfile >> dummy >> QMMMOpts.Nsteps;
    regionfile >> dummy >> QMMMOpts.Nbeads;
    regionfile >> dummy >> QMMMOpts.accratio;
    regionfile >> dummy >> QMMMOpts.Nprint;
    for (int i=0;i<Natoms;i++)
    {
      //Create path-integral beads
      for (int j=0;j<(QMMMOpts.Nbeads-1);j++)
      {
        //Pick random displacements
        double randx = (((double)rand())/((double)RAND_MAX));
        double randy = (((double)rand())/((double)RAND_MAX));
        double randz = (((double)rand())/((double)RAND_MAX));
        //Place the first bead at the initial position
        if (j == 0)
        {
          randx = 0.5;
          randy = 0.5;
          randz = 0.5;
        }
        Coord temp; //Bead coordinates
        //Set random bead displacements
        temp.x = Struct[i].P[0].x+(randx-0.5)*StepMin*Centratio;
        temp.y = Struct[i].P[0].y+(randy-0.5)*StepMin*Centratio;
        temp.z = Struct[i].P[0].z+(randz-0.5)*StepMin*Centratio;
        Struct[i].P.push_back(temp);
        Mpole temp2 = Struct[i].MP[0];
        Struct[i].MP.push_back(temp2);
        OctCharges temp3 = Struct[i].PC[0];
        Struct[i].PC.push_back(temp3);
      }
    }
  }
  if ((dummy == "MD") or (dummy == "md"))
  {
    //Read MC and PIMC options
    MDSim = 1;
    OptSim = 0;
    DFPSim = 0;
    PIMCSim = 0;
    SteepSim = 0;
    SinglePoint = 0;
    QMMMOpts.Ensemble = "NVT";
    QMMMOpts.Nbeads = 1;
    regionfile >> dummy >> QMMMOpts.dt; //Timestep
    regionfile >> dummy >> QMMMOpts.Temp; //Temperature
    QMMMOpts.Beta = 1/(k*QMMMOpts.Temp); //Inverse temperature
    regionfile >> dummy >> QMMMOpts.tautemp; //Thermostat time constant
    regionfile >> dummy >> QMMMOpts.Neq; //Number of equil. steps
    regionfile >> dummy >> QMMMOpts.Nsteps; //Number of prod. steps
    regionfile >> dummy >> QMMMOpts.Nprint; //Print frequency
    //Initialize velocity array
    for (int i=0;i<Natoms;i++)
    {
      //Only 1 bead in array
      Coord tmp;
      //Initialize with a uniform velocity distribution
      double randnum;
      randnum = (((double)rand())/((double)RAND_MAX));
      tmp.x = sqrt(kSI*QMMMOpts.Temp/(Struct[i].m*amu2kg));
      tmp.x *= m2Ang*fs2s;
      if (randnum < 0.5)
      {
        tmp.x *= -1;
      }
      randnum = (((double)rand())/((double)RAND_MAX));
      tmp.y = sqrt(kSI*QMMMOpts.Temp/(Struct[i].m*amu2kg));
      tmp.y *= m2Ang*fs2s;
      if (randnum < 0.5)
      {
        tmp.y *= -1;
      }
      randnum = (((double)rand())/((double)RAND_MAX));
      tmp.z = sqrt(kSI*QMMMOpts.Temp/(Struct[i].m*amu2kg));
      tmp.z *= m2Ang*fs2s;
      if (randnum < 0.5)
      {
        tmp.z *= -1;
      }
      Struct[i].Vel.push_back(tmp);
    }
  }
  if ((dummy == "OPT") or (dummy == "Opt") or (dummy == "opt")
     or (dummy == "TS") or (dummy == "ts"))
  {
    //Read energy minimization options
    if ((dummy == "TS") or (dummy == "ts"))
    {
      //Search for a transition state
      TranState = 1;
    }
    MDSim = 0;
    OptSim = 1;
    DFPSim = 0;
    PIMCSim = 0;
    SteepSim = 0;
    SinglePoint = 0;
    QMMMOpts.Temp = 0;
    QMMMOpts.Beta = 0;
    QMMMOpts.Press = 0;
    QMMMOpts.Neq = 0;
    QMMMOpts.Nsteps = 0;
    QMMMOpts.Nbeads = 1;
    QMMMOpts.accratio = 0;
    QMMMOpts.Nprint = 0;
    QMMMOpts.Ensemble = "N/A";
    regionfile >> dummy >> QMMMOpts.MaxStep;
    regionfile >> dummy >> QMMMOpts.MMOptTol;
    regionfile >> dummy >> QMMMOpts.MaxOptSteps;
  }
  if ((dummy == "Steep") or (dummy == "steep") or
  (dummy == "SD") or (dummy == "sd"))
  {
    //Read energy minimization options
    MDSim = 0;
    OptSim = 0;
    DFPSim = 0;
    PIMCSim = 0;
    SteepSim = 1;
    SinglePoint = 0;
    QMMMOpts.Temp = 0;
    QMMMOpts.Beta = 0;
    QMMMOpts.Press = 0;
    QMMMOpts.Neq = 0;
    QMMMOpts.Nsteps = 0;
    QMMMOpts.Nbeads = 1;
    QMMMOpts.accratio = 0;
    QMMMOpts.Nprint = 0;
    QMMMOpts.Ensemble = "N/A";
    regionfile >> dummy >> QMMMOpts.StepScale;
    regionfile >> dummy >> QMMMOpts.MaxStep;
    regionfile >> dummy >> QMMMOpts.QMOptTol;
    regionfile >> dummy >> QMMMOpts.MMOptTol;
    regionfile >> dummy >> QMMMOpts.MaxOptSteps;
  }
  if ((dummy == "bfgs") or (dummy == "BFGS") or
     (dummy == "dfp") or (dummy == "DFP"))
  {
    //Read energy minimization options for the DFP optimizer
    MDSim = 0;
    OptSim = 0;
    DFPSim = 1;
    PIMCSim = 0;
    SteepSim = 0;
    SinglePoint = 0;
    QMMMOpts.Temp = 0;
    QMMMOpts.Beta = 0;
    QMMMOpts.Press = 0;
    QMMMOpts.Neq = 0;
    QMMMOpts.Nsteps = 0;
    QMMMOpts.Nbeads = 1;
    QMMMOpts.accratio = 0;
    QMMMOpts.Nprint = 0;
    QMMMOpts.Ensemble = "N/A";
    regionfile >> dummy >> QMMMOpts.StepScale;
    regionfile >> dummy >> QMMMOpts.MaxStep;
    regionfile >> dummy >> QMMMOpts.QMOptTol;
    regionfile >> dummy >> QMMMOpts.MMOptTol;
    regionfile >> dummy >> QMMMOpts.MaxOptSteps;
  }
  if ((dummy == "SP") or (dummy == "sp") or
  (dummy == "energy") or (dummy == "Energy"))
  {
    //Read energy minimization options
    MDSim = 0;
    OptSim = 0;
    DFPSim = 0;
    PIMCSim = 0;
    SteepSim = 0;
    SinglePoint = 1;
    QMMMOpts.Temp = 0;
    QMMMOpts.Beta = 0;
    QMMMOpts.Press = 0;
    QMMMOpts.Neq = 0;
    QMMMOpts.Nsteps = 0;
    QMMMOpts.Nbeads = 1;
    QMMMOpts.accratio = 0;
    QMMMOpts.Nprint = 0;
    QMMMOpts.Ensemble = "N/A";
  }
  regionfile >> dummy >> dummy; //PBC options
  if ((dummy == "Yes") or (dummy == "yes") or
  (dummy == "YES") or (dummy == "true") or
  (dummy == "True") or (dummy == "TRUE"))
  {
    //Read box sizes
    PBCon = 1;
    regionfile >> dummy;
    regionfile >> Lx >> Ly >> Lz;
  }
  regionfile >> dummy >> Nqm; //Number of QM atoms
  for (int i=0;i<Nqm;i++)
  {
    int AtomID;
    regionfile >> AtomID;
    Struct[AtomID].QMregion = 1;
    Struct[AtomID].MMregion = 0;
  }
  regionfile >> dummy >> Npseudo; //Number of pseudo-atoms
  for (int i=0;i<Npseudo;i++)
  {
    int AtomID;
    regionfile >> AtomID;
    Struct[AtomID].PAregion = 1;
    Struct[AtomID].MMregion = 0;
  }
  regionfile >> dummy >> Nbound; //Number of boundary-atoms
  for (int i=0;i<Nbound;i++)
  {
    int AtomID;
    regionfile >> AtomID;
    Struct[AtomID].BAregion = 1;
    Struct[AtomID].MMregion = 0;
  }
  Nmm = Natoms-Nqm-Npseudo-Nbound; //Number of MM atoms
  regionfile >> dummy >> Nfreeze; //Number of frozen atoms
  for (int i=0;i<Nfreeze;i++)
  {
    int AtomID;
    regionfile >> AtomID;
    Struct[AtomID].Frozen = 1;
    if (PIMCSim == 1)
    {
      //Frozen atoms must be purely classical
      #pragma omp parallel for
      for (int j=0;j<QMMMOpts.Nbeads;j++)
      {
        Struct[AtomID].P[j].x = Struct[AtomID].P[0].x;
        Struct[AtomID].P[j].y = Struct[AtomID].P[0].y;
        Struct[AtomID].P[j].z = Struct[AtomID].P[0].z;
      }
      #pragma omp barrier
    }
  }
  //Read initial structures for all beads
  if ((CheckFile("BeadStartStruct.xyz")) and (!GauExternal))
  {
    //Print output
    cout << "Reading restart information...";
    cout << '\n' << '\n';;
    //Open file
    fstream beadfile;
    beadfile.open("BeadStartStruct.xyz",ios_base::in);
    //Read and discard number of atoms
    beadfile >> dummy;
    //Read atom/bead positions
    for (int i=0;i<Natoms;i++)
    {
      for (int j=0;j<QMMMOpts.Nbeads;j++)
      {
        //Read atom type and discard
        beadfile >> dummy;
        //Read XYZ coordinates
        beadfile >> Struct[i].P[j].x;
        beadfile >> Struct[i].P[j].y;
        beadfile >> Struct[i].P[j].z;
      }
    }
  }
  //Collect additonal TINKER input
  if ((TINKER == 1) and (!GauExternal))
  {
    //Classes are not used in the QMMM, but looking for them can spot errors
    FindTINKERClasses(Struct);
  }
  return;
};

void FLUKEErrorChecker(QMMMSettings& QMMMOpts)
{
  //Checks for basic errors and conflicts
  bool DoQuit = 0; //Bool, quit with error
  if (((TINKER+AMBER+LAMMPS) == 0) and (QMonly != 1))
  {
    //Check the MM wrappers
    cout << " Error: No valid MM wrapper selected.";
    cout << '\n';
    cout << "  Select a wrapper if you want to run this type ";
    cout << "of calculation.";
    cout << '\n';
    DoQuit = 1;
  }
  if (((PSI4+Gaussian) == 0) and (MMonly != 1))
  {
    //Check the QM wrappers
    cout << " Error: No valid QM wrapper selected.";
    cout << '\n';
    cout << "  Select a wrapper if you want to run this type ";
    cout << "of calculation.";
    cout << '\n';
    DoQuit = 1;
  }
  if ((QMMMOpts.Ensemble == "NPT") and (PBCon == 0))
  {
    //Check the PBC options
    cout << " Error: NPT simulation without PBC.";
    cout << '\n';
    cout << "  Turn PBC on if you want to run this type ";
    cout << "of calculation.";
    cout << '\n';
    DoQuit = 1;
  }
  if (Ncpus < 1)
  {
    //Checks the number of threads and continue
    cout << " Warning: Calculations cannot run with ";
    cout << Ncpus << " CPUs.";
    cout << '\n';
    if (Jokes == 1)
    {
      cout << " Do you know how computers work?";
    }
    cout << " Ncpus set to 1";
    cout << '\n';
    Ncpus = 1;
  }
  if ((PSI4 == 1) and (QMMM == 1))
  {
    if (OptSim == 1)
    {
      cout << " Error: QMMM PSI4 optimizations can only be performed with";
      cout << '\n';
      cout << " the native steepest descent or DFP.";
      cout << '\n';
      DoQuit = 1;
    }
    if ((Npseudo != 0) or (Nbound != 0))
    {
      cout << " Error: The PSI4 wrapper can only use QM and MM atoms.";
      cout << '\n';
      cout << " Remove the pseudo-atoms and boundary-atoms.";
      cout << '\n';
      DoQuit = 1;
    }
  }
  if (DoQuit == 1)
  {
    //Quits
    cout << '\n';
    cout.flush();
    exit(0);
  }
  if (DoQuit == 0)
  {
    //Sarcastically continues
    cout << "No fatal errors detected.";
    cout << '\n';
    if (Jokes == 1)
    {
      cout << " And there was much rejoicing. Yay...";
      cout << '\n';
      cout << '\n';
      cout.flush();
    }
  }
  return;
};

void FLUKEPrintSettings(QMMMSettings& QMMMOpts)
{
  //Prints out the simulation details
  cout << "Setting up simulation..." << '\n';
  cout << '\n';
  cout << "Atoms: " << Natoms << '\n';
  if (QMMM == 1)
  {
    cout << " QM atoms: " << Nqm << '\n';
    cout << " MM atoms: " << Nmm << '\n';
    cout << " Pseudo-atoms: " << Npseudo << '\n';
    cout << " Boundary-atoms: " << Nbound << '\n';
  }
  if (Nfreeze > 0)
  {
    cout << " Frozen atoms: " << Nfreeze << '\n';
  }
  if (PIMCSim == 1)
  {
    //Print input for error checking
    if (QMMMOpts.Nbeads > 1)
    {
      cout << " PI Beads: " << QMMMOpts.Nbeads << '\n';
    }
    cout << '\n';
    cout << "Equilibration steps: " << QMMMOpts.Neq << '\n';
    cout << "Steps for production run: " << QMMMOpts.Nsteps << '\n';
    cout << "Simulation mode: ";
    if (QMMM == 1)
    {
      cout << "QMMM";
    }
    if (QMonly == 1)
    {
      cout << "Pure QM";
    }
    if (MMonly == 1)
    {
      cout << "Pure MM";
    }
    cout << " " << QMMMOpts.Ensemble;
    if (QMMMOpts.Nbeads > 1)
    {
      cout << " path-integral";
    }
    cout << " Monte Carlo" << '\n';
  }
  if (MDSim == 1)
  {
    cout << '\n';
    cout << "Simulation mode: ";
    if (QMMM == 1)
    {
      cout << "QMMM";
    }
    if (QMonly == 1)
    {
      cout << "Pure QM";
    }
    if (MMonly == 1)
    {
      cout << "Pure MM";
    }
    cout << " molecular dynamics" << '\n';
    cout << " Time step: " << QMMMOpts.dt << " fs" << '\n';
    cout << " Temperature: " << QMMMOpts.Temp << " K" << '\n';
    cout << " Berendsen time-constant: " << QMMMOpts.tautemp << " fs" << '\n';
    cout << " Equilibration MD steps: " << QMMMOpts.Neq << '\n';
    cout << " Production MD steps: " << QMMMOpts.Nsteps << '\n';
  }
  if ((OptSim == 1) or (SteepSim == 1) or (DFPSim == 1))
  {
    cout << '\n';
    cout << "Simulation mode: ";
    if (QMMM == 1)
    {
      cout << "QMMM";
    }
    if (QMonly == 1)
    {
      cout << "Pure QM";
    }
    if (MMonly == 1)
    {
      cout << "Pure MM";
    }
    cout << " energy minimization" << '\n';
    if ((QMMM == 1) or (QMonly == 1))
    {
      cout << " QM";
      if (QMMM == 1)
      {
        cout << "MM";
      }
      cout << " minimizer: ";
      if (OptSim == 1)
      {
        cout << "Native QM optimizer" << '\n';
      }
      if (SteepSim == 1)
      {
        cout << "FLUKE steepest descent" << '\n';
      }
      if (DFPSim == 1)
      {
        cout << "FLUKE DFP" << '\n';
      }
    }
  }
  if (SinglePoint == 1)
  {
    cout << '\n';
    cout << "Simulation mode: ";
    if (QMMM == 1)
    {
      cout << "QMMM";
    }
    if (QMonly == 1)
    {
      cout << "Pure QM";
    }
    if (MMonly == 1)
    {
      cout << "Pure MM";
    }
    cout << " single-point energy" << '\n';
  }
  if ((QMonly == 1) or (QMMM == 1))
  {
    cout << " QM wrapper: ";
    if (PSI4 == 1)
    {
      cout << "PSI4" << '\n';
    }
    if (Gaussian == 1)
    {
      cout << "Gaussian" << '\n';
    }
    cout << " QM method: ";
    cout << QMMMOpts.Func << "/";
    cout << QMMMOpts.Basis << '\n';
  }
  if ((MMonly == 1) or (QMMM == 1))
  {
    cout << " MM wrapper: ";
    if (TINKER == 1)
    {
      cout << "TINKER" << '\n';
    }
    if (AMBER == 1)
    {
      cout << "AMBER" << '\n';
    }
    if (LAMMPS == 1)
    {
      cout << "LAMMPS" << '\n';
    }
    if (QMMM == 1)
    {
      cout << " QMMM potential: ";
      if (CHRG == 1)
      {
        cout << "Point-charge force field" << '\n';
      }
      if (AMOEBA == 1)
      {
        cout << "Polarizable force field" << '\n';
      }
    }
  }
  cout << '\n';
  //Print convergence criteria for optimizations
  if ((OptSim == 1) or (SteepSim == 1) or (DFPSim == 1))
  {
    cout << "Optimization settings:" << '\n';
    if ((SteepSim == 1) or (DFPSim == 1))
    {
      cout << " Step scale factor: " << QMMMOpts.MaxStep;
      cout << '\n';
    }
    cout << " Max step size: " << QMMMOpts.MaxStep;
    cout << " \u212B" << '\n';
    cout << " Max steps: " << QMMMOpts.MaxOptSteps;
    cout << '\n' << '\n';
    if ((SteepSim == 1) or (DFPSim == 1))
    {
      cout << "QM convergence criteria:" << '\n';
      cout << "  RMS deviation: " << QMMMOpts.QMOptTol;
      cout << " \u212B" << '\n';
      cout << "  Max force: " << (100*QMMMOpts.QMOptTol);
      cout << " eV/\u212B" << '\n';
      cout << "  RMS force: " << (50*QMMMOpts.QMOptTol);
      cout << " eV/\u212B" << '\n';
      cout << '\n';
    }
    cout << "MM convergence criteria:" << '\n';
    cout << "  RMS deviation: " << QMMMOpts.MMOptTol;
    cout << " \u212B" << '\n';
    cout << "  RMS force: " << QMMMOpts.MMOptTol*kcal2eV;
    cout << " eV/\u212B" << '\n';
    cout << '\n';
  }
  cout.flush(); //Flush for output being redirected to a file
  return;
};

void GetQuotes(vector<string>& Quotes)
{
  //Generate random quotes
  string dummy; //Generic string
  dummy = "\'It is difficult to prove that this quote is not random.\'";
  dummy += '\n';
  dummy += "                                           -Eric G. Kratz";
  for (int i=0;i<1000;i++)
  {
    //Add quotes to the list
    Quotes.push_back(dummy);
  }
  return;
};
