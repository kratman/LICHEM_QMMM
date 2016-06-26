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
  if (AMBER)
  {
    Emm = AMBERForces(QMMMData,forces,QMMMOpts,bead);
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

