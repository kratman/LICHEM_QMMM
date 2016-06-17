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
  vector<QMMMAtom> Struct; //Atomic data
  QMMMSettings QMMMOpts; //Simulation settings
  int DerType = 0; //Type of derivatives
  int Bead = 0; //Which replica
  stringstream call; //Stream for system calls and reading/writing files
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
      xyzfile.open(call.str().c_str(),ios_base::in);
    }
    if (dummy == "-c")
    {
      //Open the connectivity file
      connectfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-r")
    {
      //Open the region file
      regionfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-b")
    {
      //Read the current bead
      Bead = atoi(argv[i+1]);
    }
  }
  //Open files passed by Gaussian
  GauInput.open(argv[12],ios_base::in);
  GauOutput.open(argv[13],ios_base::out);
  GauMsg.open(argv[14],ios_base::out);
  //Read LICHEM input
  ReadLICHEMInput(xyzfile,connectfile,regionfile,Struct,QMMMOpts);
  //Set degrees of freedom
  int Ndof = 3*(Nqm+Npseudo); //Number of QM and PB degrees of freedom
  //Read g09 input for new QM atom positions
  getline(GauInput,dummy);
  stringstream line(dummy);
  line >> dummy >> DerType;
  if (DerType == 2)
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
  VectorXd Forces(Ndof); //Forces for QM and PB
  Forces.setZero();
  fstream MMgrad,QMlog; //QMMM output
  //QM forces
  Eqm = GaussianForces(Struct,Forces,QMMMOpts,Bead);
  //MM forces
  if (TINKER)
  {
    Emm = TINKERForces(Struct,Forces,QMMMOpts,Bead);
    if (AMOEBA or QMMMOpts.UseImpSolv)
    {
      //Calculate polarization forces for AMOEBA
      Emm += TINKERPolForces(Struct,Forces,QMMMOpts,Bead);
    }
  }
  if (AMBER)
  {
    Emm = AMBERForces(Struct,Forces,QMMMOpts,Bead);
  }
  if (LAMMPS)
  {
    Emm = LAMMPSForces(Struct,Forces,QMMMOpts,Bead);
  }
  //Write formatted output for g09
  double E = (Eqm+Emm)/Har2eV; //Calculate
  GauOutput << left; //More formatting
  GauOutput << LICHEMFormFloat(E,20); //QM+MM partial energy
  GauOutput << LICHEMFormFloat(0.0,20); //Dipole moment
  GauOutput << LICHEMFormFloat(0.0,20); //Dipole moment
  GauOutput << LICHEMFormFloat(0.0,20); //Dipole moment
  GauOutput << '\n';
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Write forces
    GauOutput << LICHEMFormFloat(-1*Forces(3*i)*BohrRad/Har2eV,20);
    GauOutput << LICHEMFormFloat(-1*Forces(3*i+1)*BohrRad/Har2eV,20);
    GauOutput << LICHEMFormFloat(-1*Forces(3*i+2)*BohrRad/Har2eV,20);
    GauOutput << '\n';
  }
  GauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  GauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  GauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  GauOutput << '\n';
  GauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  GauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  GauOutput << LICHEMFormFloat(0.0,20); //Polarizability
  GauOutput << '\n';
  for (int i=0;i<Ndof;i++)
  {
    //Dipole derivatives
    GauOutput << LICHEMFormFloat(0.0,20);
    GauOutput << LICHEMFormFloat(0.0,20);
    GauOutput << LICHEMFormFloat(0.0,20);
    GauOutput << '\n';
  }
  //Write output and close the file
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
    ofile << LICHEMFormFloat(Struct[i].P[Bead].x,16) << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].y,16) << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].z,16) << '\n';
  }
  ofile.flush();
  ofile.close();
  //Return to Gaussian
  cout << "Forces were returned to Gaussian..." << '\n';
  cout.flush();
  exit(0);
  return;
};

