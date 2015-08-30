/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for Gaussian. These routines are written for g09
 Rev D. Note that the external function needs to call MM codes.

 Reference for Gaussian:
 Frisch et al., Gaussian 09 Rev D.01, (2009)

*/

//QM utility functions
void ExternalGaussian(int& argc, char**& argv)
{
  //This function is an "external script" that can be called by
  //Gaussian's external interface
  double Eqm = 0; //Stores partial energies
  double Emm = 0; //Stores partial energies
  vector<QMMMAtom> Struct; //Atomic data
  QMMMSettings QMMMOpts; //Simulation settings
  int DerType,Bead;
  stringstream call;
  call.copyfmt(cout); //Copy print settings
  string dummy,Stub; //Generic strings
  //Declare lots of file streams
  fstream xyzfile,connectfile,regionfile; //LICHEM streams
  fstream GauInput,GauOutput,GauMsg,GauFchk,GauMatrix; //Gaussian streams
  fstream ofile,ifile; //Generic streams
  //Read arguments
  for (int i=0;i<argc;i++)
  {
    //Read file names and CPUs
    dummy = string(argv[i]);
    if (dummy == "-n")
    {
      Ncpus = atoi(argv[i+1]);
      //Set OpenMP threads for the external routine
      omp_set_num_threads(Ncpus);
    }
    if (dummy == "-GauExtern")
    {
      //Get the QMMM filename and xyzfile
      Stub = string(argv[i+1]);
      call.str("");
      call << Stub << ".xyz";
      xyzfile.open(call.str().c_str(),ios_base::in);
    }
    if (dummy == "-c")
    {
      connectfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-r")
    {
      regionfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-b")
    {
      Bead = atoi(argv[i+1]);
    }
  }
  GauInput.open(argv[12],ios_base::in);
  GauOutput.open(argv[13],ios_base::out);
  GauMsg.open(argv[14],ios_base::out);
  //Read LICHEM input
  ReadLICHEMInput(xyzfile,connectfile,regionfile,Struct,QMMMOpts);
  //Read g09 input for new QM atom positions
  getline(GauInput,dummy);
  stringstream line(dummy);
  line >> dummy >> DerType;
  if (DerType == 2)
  {
    cerr << "Error: Second derivatives of the energy were requested!!!";
    cerr << '\n';
    cerr << "Something is wrong.";
    cerr << '\n';
    cerr.flush();
    exit(0);
  }
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion or Struct[i].PBregion)
    {
      //Save atom information
      getline(GauInput,dummy);
      stringstream line(dummy);
      line >> dummy;
      line >> Struct[i].P[Bead].x;
      line >> Struct[i].P[Bead].y;
      line >> Struct[i].P[Bead].z;
      //Change units
      Struct[i].P[Bead].x *= BohrRad;
      Struct[i].P[Bead].y *= BohrRad;
      Struct[i].P[Bead].z *= BohrRad;
    }
  }
  GauInput.close();
  //Calculate the QMMM forces
  VectorXd Forces(3*(Nqm+Npseudo)); //Forces for QM and PB
  Forces.setZero();
  fstream MMgrad,QMlog; //QMMM output
  //QM forces
  Eqm = GaussianForces(Struct,Forces,QMMMOpts,Bead);
  //MM forces
  if (TINKER == 1)
  {
    Emm = TINKERForces(Struct,Forces,QMMMOpts,Bead);
    if (AMOEBA)
    {
      Emm += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
    }
  }
  if (AMBER == 1)
  {
    Emm = AMBERForces(Struct,Forces,QMMMOpts,Bead);
  }
  if (LAMMPS == 1)
  {
    Emm = LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
  }
  //Write formatted output for g09
  double E = (Eqm+Emm)/Har2eV; //Calculate
  GauOutput << fixed; //Formatting
  GauOutput << left; //More formatting
  GauOutput.precision(12);
  GauOutput << setw(20) << E; //QM+MM partial energy
  GauOutput << setw(20) << 0.0; //Dipole moment
  GauOutput << setw(20) << 0.0; //Dipole moment
  GauOutput << setw(20) << 0.0; //Dipole moment
  GauOutput << '\n';
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Write forces
    GauOutput << setw(20) << (-1*Forces(3*i)*BohrRad/Har2eV);
    GauOutput << setw(20) << (-1*Forces(3*i+1)*BohrRad/Har2eV);
    GauOutput << setw(20) << (-1*Forces(3*i+2)*BohrRad/Har2eV);
    GauOutput << '\n';
  }
  GauOutput << setw(20) << 0.0; //Polarizability
  GauOutput << setw(20) << 0.0; //Polarizability
  GauOutput << setw(20) << 0.0; //Polarizability
  GauOutput << '\n';
  GauOutput << setw(20) << 0.0; //Polarizability
  GauOutput << setw(20) << 0.0; //Polarizability
  GauOutput << setw(20) << 0.0; //Polarizability
  GauOutput << '\n';
  for (int i=0;i<(3*(Nqm+Npseudo));i++)
  {
    //Dipole derivatives
    GauOutput << setw(20) << 0.0;
    GauOutput << setw(20) << 0.0;
    GauOutput << setw(20) << 0.0;
    GauOutput << '\n';
  }
  GauOutput.flush();
  GauOutput.close();
  //Write new XYZ for recovery of failed optimizations
  call.str("");
  call << Stub << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Write XYZ coordinates
    ofile << Struct[i].QMTyp << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].x << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].y << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].z << '\n';
  }
  ofile.flush();
  ofile.close();
  //Return to Gaussian
  cout << "Forces were returned to Gaussian..." << '\n';
  cout.flush();
  exit(0);
  return;
};

//QM wrapper functions
double GaussianForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces on a set of atoms
  stringstream call;
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream ofile,ifile,QMlog; //Generic input files
  double Eqm = 0; //QM energy
  double Eself = 0; //External field self-energy
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  bool UseCheckPoint = CheckFile(call.str());
  //Construct Gaussian input
  call.str("");
  call << "#P " << QMMMOpts.Func << "/";
  call << QMMMOpts.Basis << " Force=NoStep Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)"; //Line ended further below
  if (UseCheckPoint)
  {
    //Restart if possible
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if (Npseudo > 0)
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    call << "Charge=angstroms "; //Read charges
    call << "Population(MK,ReadRadii)";
    call << '\n';
  }
  WriteGauInput(Struct,call.str(),QMMMOpts,Bead);
  //Run Gaussian
  call.str("");
  call << "g09 " << "LICHM_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Extract forces
  call.str("");
  call << "LICHM_" << Bead << ".log";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool GradDone = 0;
  while ((!QMlog.eof()) and (!GradDone))
  {
    //Parse file line by line
    getline(QMlog,dummy);
    stringstream line(dummy);
    line >> dummy;
    //This only works with #P
    if (dummy == "Center")
    {
      line >> dummy >> dummy;
      if (dummy == "Forces")
      {
        GradDone = 1; //Not grad school, that lasts forever
        getline(QMlog,dummy); //Clear junk
        getline(QMlog,dummy); //Clear more junk
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Extract forces; Convoluted, but "easy"
          getline(QMlog,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> Fx;
          line >> Fy;
          line >> Fz;
          //Save forces
          if (abs(Fx) >= 1e-8)
          {
            Forces(3*i) += Fx*Har2eV/BohrRad;
          }
          if (abs(Fy) >= 1e-8)
          {
            Forces(3*i+1) += Fy*Har2eV/BohrRad;
          }
          if (abs(Fz) >= 1e-8)
          {
            Forces(3*i+2) += Fz*Har2eV/BohrRad;
          }
        }
      }
    }
    if (dummy == "Self")
    {
      line >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> Eself;
      }
    }
    //Check for partial QMMM energy
    if (dummy == "SCF")
    {
      line >> dummy;
      if (dummy == "Done:")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> Eqm;
      }
    }
    //Check for charges
    if (dummy == "ESP")
    {
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
  }
  QMlog.close();
  if (!GradDone)
  {
    cerr << "Warning: No forces recovered!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to recover...";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".com";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  Eqm -= Eself;
  Eqm *= Har2eV;
  return Eqm;
};

void GaussianCharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Function to update QM point-charges
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout); //Copy print settings
  //Remove multipoles file
  call.str("");
  call << "rm -f MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  bool UseCheckPoint = CheckFile(call.str());
  //Construct Gaussian input
  call.str("");
  call << "#P " << QMMMOpts.Func << "/";
  call << QMMMOpts.Basis << " SP Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)" << '\n';
  if (QMMM)
  {
    if (Npseudo > 0)
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (UseCheckPoint)
    {
      //Read pseudo potential
      call << "Guess=TCheck ";
    }
    call << "Charge=angstroms "; //Read charges
    call << "Population(MK,ReadRadii)";
    call << '\n';
  }
  WriteGauInput(Struct,call.str(),QMMMOpts,Bead);
  //Run QM calculation
  call.str("");
  call << "g09 LICHM_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "ESP")
    {
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
  }
  ifile.close();
  //Clean up files and save checkpoint file
  call.str("");
  call << "mv LICHM_" << Bead;
  call << ".chk tmp_" << Bead;
  call << ".chk && ";
  call << "rm -f ";
  call << "LICHM_" << Bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".com";
  call << " && mv tmp_" << Bead;
  call << ".chk LICHM_" << Bead;
  call << ".chk";
  GlobalSys = system(call.str().c_str());
  return;
};

double GaussianEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Calculates the QM energy with Gaussian
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  double Eself = 0.0; //External field self-energy
  //Remove multipole file
  call.str("");
  call << "rm -f MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  bool UseCheckPoint = CheckFile(call.str());
  //Construct Gaussian input
  call.str("");
  call << "#P " << QMMMOpts.Func << "/";
  call << QMMMOpts.Basis << " SP Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)";
  if (UseCheckPoint)
  {
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if (Npseudo > 0)
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    call << "Charge=angstroms "; //Read charges
    call << "Population(MK,ReadRadii)";
    call << '\n';
  }
  WriteGauInput(Struct,call.str(),QMMMOpts,Bead);
  //Calculate energy
  call.str("");
  call << "g09 ";
  call << "LICHM_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Read output
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  bool QMfinished = 0;
  while (!ifile.eof())
  {
    stringstream line;
    getline(ifile,dummy);
    line.str(dummy);
    line >> dummy;
    //Search for field self-energy
    if (dummy == "Self")
    {
      line >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> Eself;
      }
    }
    //Search for energy
    if (dummy == "SCF")
    {
      line >> dummy;
      if (dummy == "Done:")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> E;
        QMfinished = 1;
      }
    }
    //Check for charges
    if (dummy == "ESP")
    {
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
  }
  if (!QMfinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  ifile.close();
  //Clean up files and save checkpoint file
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".com";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E -= Eself;
  E *= Har2eV;
  return E;
};

double GaussianOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Runs Gaussian optimizations with GauExternal
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  int ExtCPUs = 1; //Number of CPUs for GauExternal
  if (AMOEBA and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Write a new XYZ
  //Note: GauExternal needs different input than the rest of the wrappers
  call.str("");
  call << "LICHMExt_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Print XYZ coordinates
    ofile << setprecision(12) << Struct[i].QMTyp << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].x << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].y << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].z << '\n';
  }
  ofile.flush();
  ofile.close();
  //Write multipole point-charges
  if (AMOEBA and QMMM)
  {
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile.copyfmt(cout);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile << '\n';
      }
    }
    ofile.copyfmt(cout);
    ofile.flush();
    ofile.close();
  }
  //Write Gaussian input
  call.str("");
  call << "LICHMExt_" << Bead << ".com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=LICHMExt_" << Bead << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM;
  if (QMMMOpts.MemMB)
  {
    call << "MB";
  }
  else
  {
    call << "GB";
  }
  call << '\n';
  call << "%NprocShared=";
  if (Ncpus <= 2)
  {
    call << 1;
  }
  else
  {
    call << 2;
    ExtCPUs = Ncpus-2;
  }
  call << '\n';
  call << "#P " << "external=\"lichem -GauExtern ";
  call << "LICHMExt_" << Bead;  //Just the stub
  call << " -n " << ExtCPUs;
  call << " -c " << confilename;
  call << " -r " << regfilename;
  call << " -b " << Bead;
  call << "\"" << '\n';
  call << "Symmetry=None Opt=(";
  call << "MaxCycles=" << QMMMOpts.MaxOptSteps;
  call << ",MaxStep=" << int(round((QMMMOpts.MaxStep/(0.01*BohrRad))));
  call << ")" << '\n';
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.Charge << " " << QMMMOpts.Spin << '\n';
  //Add atoms
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion)
    {
      call << Struct[i].QMTyp;
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
    if (Struct[i].PBregion)
    {
      call << "F";
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
  }
  call << '\n'; //Blank line needed
  //Write Gaussian input
  ofile << call.str();
  ofile.close();
  //Run Optimization
  call.str("");
  call << "g09 ";
  call << "LICHMExt_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHMExt_";
  call << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  bool Optfinished = 0;
  while (!ifile.eof())
  {
    //This loop will find the last geometry even if the calculation
    //did not complete
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Center")
    {
      line >> dummy >> dummy;
      line >> dummy;
      if (dummy == "Coordinates")
      {
        //Clear junk
        getline(ifile,dummy);
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Get new coordinates
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy;
            line >> Struct[i].P[Bead].x;
            line >> Struct[i].P[Bead].y;
            line >> Struct[i].P[Bead].z;
          }
        }
      }
    }
    if (dummy == "--")
    {
      line >> dummy;
      if (dummy == "Stationary")
      {
        Optfinished = 1;
      }
    }
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f LICHMExt_";
  call << Bead << ".*";
  call << " MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Print warnings and errors
  if (!Optfinished)
  {
    cerr << "Warning: Optimization did not converge!!!";
    cerr << '\n';
    cerr << "An older geometry will be recovered...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Calculate new point-charges and return
  GaussianCharges(Struct,QMMMOpts,Bead);
  return E;
};

