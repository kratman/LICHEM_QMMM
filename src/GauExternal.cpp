/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper for Gaussian's external interface. These routines are written
 for g09 Rev D. Note that the external function needs to call MM codes.

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
  vector<QMMMAtom> QMMMData; //Atomic data
  QMMMSettings QMMMOpts; //Simulation settings
  int derType = 0; //Type of derivatives
  int bead = 0; //Which replica
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy,Stub; //Generic strings
  //Declare lots of file streams
  fstream xyzFile,connectFile,regionFile; //LICHEM streams
  fstream gauInput,gauOutput,gauMsg,gauFchk,gauMatrix; //Gaussian streams
  fstream outFile; //Generic streams
  //Read arguments
  for (int i=0;i<argc;i++)
  {
    //Read file names and CPUs
    dummy = string(argv[i]);
    if (dummy == "-n")
    {
      Ncpus = atoi(argv[i+1]);
      //Set OpenMP threads for the external routine
      #ifdef _OPENMP
        omp_set_num_threads(Ncpus);
      #endif
    }
    if (dummy == "-GauExtern")
    {
      //Get the QMMM filename and xyzfile
      Stub = string(argv[i+1]);
      call.str("");
      call << Stub << ".xyz";
      xyzFile.open(call.str().c_str(),ios_base::in);
    }
    if (dummy == "-c")
    {
      //Open the connectivity file
      connectFile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-r")
    {
      //Open the region file
      regionFile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-b")
    {
      //Read the current bead
      bead = atoi(argv[i+1]);
    }
  }
  //Open files passed by Gaussian
  gauInput.open(argv[12],ios_base::in);
  gauOutput.open(argv[13],ios_base::out);
  gauMsg.open(argv[14],ios_base::out);
  //Read LICHEM input
  ReadLICHEMInput(xyzFile,connectFile,regionFile,QMMMData,QMMMOpts);
  //Set degrees of freedom
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Read g09 input for new QM atom positions
  getline(gauInput,dummy);
  stringstream line(dummy);
  line >> dummy >> derType;
  if (derType == 2)
  {
    //Make sure Gaussian is only requesting energies or forces
    cerr << "Error: Second derivatives of the energy were requested!!!";
    cerr << '\n';
    cerr << "Something is wrong.";
    cerr << '\n';
    cerr.flush();
    exit(0);
  }
  //Read updated positions from Gaussian files
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
    {
      //Save atom information
      getline(gauInput,dummy);
      stringstream line(dummy);
      line >> dummy;
      line >> QMMMData[i].P[bead].x;
      line >> QMMMData[i].P[bead].y;
      line >> QMMMData[i].P[bead].z;
      //Change units
      QMMMData[i].P[bead].x *= bohrRad;
      QMMMData[i].P[bead].y *= bohrRad;
      QMMMData[i].P[bead].z *= bohrRad;
    }
  }
  gauInput.close();
  //Calculate the QMMM forces
  VectorXd forces(Ndof); //Forces for QM and PB
  forces.setZero();
  fstream MMgrad,QMlog; //QMMM output
  //QM forces
  Eqm = GaussianForces(QMMMData,forces,QMMMOpts,bead);
  //MM forces
  if (TINKER)
  {
    Emm = TINKERForces(QMMMData,forces,QMMMOpts,bead);
    if (AMOEBA or QMMMOpts.useImpSolv)
    {
      //Calculate polarization forces for AMOEBA
      Emm += TINKERPolForces(QMMMData,forces,QMMMOpts,bead);
    }
  }
  if (LAMMPS)
  {
    Emm = LAMMPSForces(QMMMData,forces,QMMMOpts,bead);
  }
  //Write formatted output for g09
  double E = (Eqm+Emm)/har2eV; //Calculate
  gauOutput << left; //More formatting
  gauOutput << LICHEMFormFloat(E,20); //QM+MM partial energy
  gauOutput << LICHEMFormFloat(0.0,20); //Dipole moment
  gauOutput << LICHEMFormFloat(0.0,20); //Dipole moment
  gauOutput << LICHEMFormFloat(0.0,20); //Dipole moment
  gauOutput << '\n';
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Write forces
    gauOutput << LICHEMFormFloat(-1*forces(3*i)*bohrRad/har2eV,20);
    gauOutput << LICHEMFormFloat(-1*forces(3*i+1)*bohrRad/har2eV,20);
    gauOutput << LICHEMFormFloat(-1*forces(3*i+2)*bohrRad/har2eV,20);
    gauOutput << '\n';
  }
  gauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  gauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  gauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  gauOutput << '\n';
  gauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  gauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  gauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  gauOutput << '\n';
  for (int i=0;i<Ndof;i++)
  {
    //Dipole derivatives
    gauOutput << LICHEMFormFloat(0.0,20);
    gauOutput << LICHEMFormFloat(0.0,20);
    gauOutput << LICHEMFormFloat(0.0,20);
    gauOutput << '\n';
  }
  //Write output and close the file
  gauOutput.flush();
  gauOutput.close();
  //Write new XYZ for recovery of failed optimizations
  call.str("");
  call << Stub << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Write XYZ coordinates
    outFile << QMMMData[i].QMTyp << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16) << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16) << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16) << '\n';
  }
  outFile.flush();
  outFile.close();
  //Return to Gaussian
  cout << "Forces were returned to Gaussian..." << '\n';
  cout.flush();
  exit(0);
  return;
};

double GaussianExternOpt(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                         int bead)
{
  //Runs Gaussian optimizations with GauExternal
  fstream inFile,QMlog; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  int extCPUs = 1; //Number of CPUs for GauExternal
  if (AMOEBA and TINKER)
  {
    RotateTINKCharges(QMMMData,bead);
  }
  //Write a new XYZ
  //NB: GauExternal needs different input than the rest of the wrappers
  call.str("");
  call << "LICHMExt_" << bead << ".xyz";
  inFile.open(call.str().c_str(),ios_base::out);
  inFile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Print XYZ coordinates
    inFile << QMMMData[i].QMTyp << " ";
    inFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16) << " ";
    inFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16) << " ";
    inFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16) << '\n';
  }
  inFile.flush();
  inFile.close();
  //Write Gaussian input
  call.str("");
  call << "LICHMExt_" << bead << ".com";
  inFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=LICHMExt_" << bead << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM;
  if (QMMMOpts.memMB)
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
    extCPUs = Ncpus-2;
  }
  call << '\n';
  call << "#P " << "external=\"lichem -GauExtern ";
  call << "LICHMExt_" << bead;  //Just the stub
  call << " -n " << extCPUs;
  call << " -c " << conFilename;
  call << " -r " << regFilename;
  call << " -b " << bead;
  call << "\"" << '\n';
  call << "Symmetry=None Opt=(";
  call << "MaxCycles=" << QMMMOpts.maxOptSteps;
  call << ",MaxStep=" << int(round((QMMMOpts.maxStep/(0.01*bohrRad))));
  call << ")" << '\n';
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.charge << " " << QMMMOpts.spin << '\n';
  //Add atoms
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].QMRegion)
    {
      call << QMMMData[i].QMTyp;
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
      call << '\n';
    }
    if (QMMMData[i].PBRegion)
    {
      call << "F";
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
      call << '\n';
    }
  }
  call << '\n'; //Blank line needed
  //Write Gaussian input
  inFile << call.str();
  inFile.close();
  //Run Optimization
  call.str("");
  call << "g09 ";
  call << "LICHMExt_" << bead;
  globalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHMExt_";
  call << bead << ".log";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool optFinished = 0;
  while (!QMlog.eof())
  {
    //This loop will find the last geometry even if the calculation
    //did not complete
    getline(QMlog,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Center")
    {
      line >> dummy >> dummy;
      line >> dummy;
      if (dummy == "Coordinates")
      {
        //Clear junk
        getline(QMlog,dummy);
        getline(QMlog,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Get new coordinates
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy;
            line >> QMMMData[i].P[bead].x;
            line >> QMMMData[i].P[bead].y;
            line >> QMMMData[i].P[bead].z;
          }
        }
      }
    }
    if (dummy == "--")
    {
      line >> dummy;
      if (dummy == "Stationary")
      {
        optFinished = 1;
      }
    }
  }
  QMlog.close();
  //Clean up files
  call.str("");
  call << "rm -f LICHMExt_";
  call << bead << ".*";
  globalSys = system(call.str().c_str());
  //Print warnings and errors
  if (!optFinished)
  {
    cerr << "Warning: Optimization did not converge!!!";
    cerr << '\n';
    cerr << "An older geometry will be recovered...";
    cerr << '\n';
    E = hugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Calculate new point-charges and return
  GaussianCharges(QMMMData,QMMMOpts,bead);
  return E;
};

