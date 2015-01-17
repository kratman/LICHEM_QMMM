/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for Gaussian. These routines are written for g09.
 Note that the external function needs to interface with all MM codes.

 Citation for Gaussian:
 Frisch et al. Gaussian 09 Rev D.01 2009

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
  int sys,DerType,ct;
  stringstream call;
  call.copyfmt(cout);
  string dummy,Stub;
  fstream xyzfile,connectfile,regionfile;
  fstream GauInput,GauOutput,GauMsg,GauFchk,GauMatrix;
  fstream ofile,ifile;
  string TINKKeyFile = "tinker.key";
  int MaxTINKERNum = 3500;
  int MaxTINKERClass = 100;
  //Read arguments
  for (int i=0;i<argc;i++)
  {
    //Read file names and CPUs
    dummy = string(argv[i]);
    if (dummy == "-n")
    {
      Ncpus = atoi(argv[i+1]);
      call.str("");
      call << "export OMP_NUM_THREADS=" << Ncpus;
      sys = system(call.str().c_str());
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
  }
  GauInput.open(argv[10],ios_base::in);
  GauOutput.open(argv[11],ios_base::out);
  GauMsg.open(argv[12],ios_base::out);
  //Read FLUKE input
  ReadFLUKEInput(xyzfile,connectfile,regionfile,Struct,QMMMOpts);
  //Read g09 input for new QM atom positions
  getline(GauInput,dummy);
  stringstream line(dummy);
  line >> dummy >> DerType;
  if (DerType == 2)
  {
    cout << "Error: Second derivatives of the energy were requested!!!";
    cout << '\n';
    cout << "Something is wrong.";
    cout << endl;
    exit(0);
  }
  for (int i=0;i<Natoms;i++)
  {
    if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
    {
      //Save atom information
      getline(GauInput,dummy);
      stringstream line(dummy);
      line >> dummy;
      line >> Struct[i].P[0].x;
      line >> Struct[i].P[0].y;
      line >> Struct[i].P[0].z;
      //Change units
      Struct[i].P[0].x *= BohrRad;
      Struct[i].P[0].y *= BohrRad;
      Struct[i].P[0].z *= BohrRad;
    }
  }
  GauInput.close();
  //Calculate the QMMM forces
  vector<Coord> Forces; //Forces for QM and PA
  fstream MMgrad,QMlog; //QMMM output
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Create arrays with zeros
    Coord tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    Forces.push_back(tmp);
  }
  //MM forces
  if (TINKER == 1)
  {
    Emm = TINKERForces(Struct,Forces,QMMMOpts,0);
  }
  if (AMBER == 1)
  {
    
  }
  if (LAMMPS == 1)
  {
    
  }
  //QM forces
  Eqm = GaussianForces(Struct,Forces,QMMMOpts,0);
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
    GauOutput << setw(20) << (-1*Forces[i].x*BohrRad/Har2eV);
    GauOutput << setw(20) << (-1*Forces[i].y*BohrRad/Har2eV);
    GauOutput << setw(20) << (-1*Forces[i].z*BohrRad/Har2eV);
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
    ofile << Struct[i].QMTyp << " ";
    ofile << setprecision(12) << Struct[i].P[0].x << " ";
    ofile << setprecision(12) << Struct[i].P[0].y << " ";
    ofile << setprecision(12) << Struct[i].P[0].z << '\n';
  }
  ofile.flush();
  ofile.close();
  //Return to Gaussian
  exit(0);
  return;
};

//QM wrapper functions
double GaussianForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces on a set of atoms
  int sys;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  fstream ofile,ifile,QMlog;
  double Eqm = 0;
  double Eself = 0;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Construct g09 input
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=QMMM";
  call << "_" << Bead;
  call << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM << "GB" << '\n';
  call << "%NprocShared=" << Ncpus << '\n';
  call << "%NoSave" << '\n'; //Deletes files
  call << "#P " << QMMMOpts.Func << "/";
  call << QMMMOpts.Basis << " Force=NoStep Symmetry=None" << '\n';
  if (QMMM == 1)
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
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.Charge << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion == 1)
    {
      call << Struct[i].QMTyp;
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
    if (Struct[i].PAregion == 1)
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
  //Add the MM field
  if (CHRG == 1)
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].P[Bead].x;
        call << " " << setprecision(12) << Struct[i].P[Bead].y;
        call << " " << setprecision(12) << Struct[i].P[Bead].z;
        call << " " << setprecision(12) << Struct[i].MP[Bead].q;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  if (AMOEBA == 1)
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  //Add basis set information from the BASIS file
  ifile.open("BASIS",ios_base::in);
  if (ifile.good())
  {
    while (!ifile.eof())
    {
       //Copy BASIS line by line, if BASIS exists
       getline(ifile,dummy);
       call << dummy << '\n';
    }
    ifile.close();
    call << '\n'; //Blank line needed
  }
  ofile << call.str();
  ofile.flush();
  ofile.close();
  //Run Gaussian
  call.str("");
  call << "g09 " << "QMMM";
  call << "_" << Bead;
  sys = system(call.str().c_str());
  //Extract forces
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Create arrays with zeros
    Coord tmp;
    tmp.x = 0;
    tmp.y = 0;
    tmp.z = 0;
    Forces.push_back(tmp);
  }
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".log";
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
          Forces[i].x += Fx*Har2eV/BohrRad;
          Forces[i].y += Fy*Har2eV/BohrRad;
          Forces[i].z += Fz*Har2eV/BohrRad;
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
          if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            if (Struct[i].QMregion == 1)
            {
              //Only collect charges for QM atoms
              stringstream line(dummy);
              line >> dummy >> dummy;
              line >> Struct[i].MP[Bead].q;
            }
          }
        }
      }
    }
  }
  QMlog.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM";
  call << "_" << Bead;
  call << ".*";
  sys = system(call.str().c_str());
  //Return
  Eqm -= Eself;
  Eqm *= Har2eV;
  return Eqm;
};

void GaussianCharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Function to update QM point charges
  int sys;
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Construct input
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=QMMM";
  call << "_" << Bead;
  call << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM << "GB" << '\n';
  call << "%NprocShared=" << Ncpus << '\n';
  call << "%NoSave" << '\n'; //Deletes files
  call << "#P " << QMMMOpts.Func << "/";
  call << QMMMOpts.Basis << " SP Symmetry=None" << '\n';
  if (QMMM == 1)
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
  call << '\n';
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.Charge << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion == 1)
    {
      call << Struct[i].QMTyp;
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
    if (Struct[i].PAregion == 1)
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
  //Add the MM field
  if ((CHRG == 1) and (QMMM == 1))
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].P[Bead].x;
        call << " " << setprecision(12) << Struct[i].P[Bead].y;
        call << " " << setprecision(12) << Struct[i].P[Bead].z;
        call << " " << setprecision(12) << Struct[i].MP[Bead].q;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  if ((AMOEBA == 1) and (QMMM == 1))
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  //Add basis set information from the BASIS file
  ifile.open("BASIS",ios_base::in);
  if (ifile.good())
  {
    while (!ifile.eof())
    {
      //Copy BASIS line by line, if BASIS exists
      getline(ifile,dummy);
      call << dummy << '\n';
    }
    ifile.close();
    call << '\n'; //Blank line needed
  }
  //Write Gaussian input
  ofile << call.str();
  ofile.close();
  //Run QM calculation
  call.str("");
  call << "g09 QMMM";
  call << "_" << Bead;
  sys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".log";
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
          if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            if (Struct[i].QMregion == 1)
            {
              //Only collect charges for QM atoms
              stringstream line(dummy);
              line >> dummy >> dummy;
              line >> Struct[i].MP[Bead].q;
            }
          }
        }
      }
    }
  }
  //Clean up files
  call.str("");
  call << "rm -f QMMM";
  call << "_" << Bead;
  call << ".*";
  sys = system(call.str().c_str());
  return;
};

double GaussianEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Calculates the QM energy with Gaussian
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0; //QM energy
  double Eself = 0.0; //Field self-energy
  int sys;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  call.str("");
  call << "QMMM_" << Bead << ".com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=QMMM";
  call << "_" << Bead << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM << "GB" << '\n';
  call << "%NprocShared=" << Ncpus << '\n';
  call << "%NoSave" << '\n'; //Deletes files
  call << "#P " << QMMMOpts.Func << "/";
  call << QMMMOpts.Basis << " SP Symmetry=None" << '\n';
  if (QMMM == 1)
  {
    if (Npseudo > 0)
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    call << "Charge=angstroms"; //Read charges
    call << '\n';
  }
  call << '\n';
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.Charge << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion == 1)
    {
      call << Struct[i].QMTyp;
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
    if (Struct[i].PAregion == 1)
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
  //Add the MM field
  if ((CHRG == 1) and (QMMM == 1))
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].P[Bead].x;
        call << " " << setprecision(12) << Struct[i].P[Bead].y;
        call << " " << setprecision(12) << Struct[i].P[Bead].z;
        call << " " << setprecision(12) << Struct[i].MP[Bead].q;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  if ((AMOEBA == 1) and (QMMM == 1))
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  //Add basis set information from the BASIS file
  ifile.open("BASIS",ios_base::in);
  if (ifile.good())
  {
    while (!ifile.eof())
    {
      //Copy BASIS line by line, if BASIS exists
      getline(ifile,dummy);
      call << dummy << '\n';
    }
    ifile.close();
    call << '\n'; //Blank line needed
  }
  //Write Gaussian input
  ofile << call.str();
  ofile.close();
  //Calculate energy
  call.str("");
  call << "g09 ";
  call << "QMMM_" << Bead;
  sys = system(call.str().c_str());
  //Read output
  call.str("");
  call << "QMMM_" << Bead << ".log";
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
  }
  if (QMfinished != 1)
  {
    cout << "Warning: SCF did not converge!!!";
    cout << '\n';
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
  }
  //Remove files
  call.str("");
  call << "rm -f QMMM";
  call << "_" << Bead;
  call << ".*";
  sys = system(call.str().c_str());
  //Change units
  E -= Eself;
  E *= Har2eV;
  return E;
};

double GaussianOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Runs Gaussian optimizations with GauExternal
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0; //QM energy
  double Eself = 0.0; //Field self-energy
  int sys;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Write a new XYZ
  call.str("");
  call << "QMMMExt_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    ofile << setprecision(12) << Struct[i].QMTyp << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].x << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].y << " ";
    ofile << setprecision(12) << Struct[i].P[Bead].z << '\n';
  }
  ofile.flush();
  ofile.close();
  //Write multipole point-charges
  if ((AMOEBA == 1) and (QMMM == 1))
  {
    ofile.open("CHARGES",ios_base::out);
    ofile.copyfmt(cout);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
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
  call << "QMMMExt_" << Bead << ".com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=QMMMExt";
  call << "_" << Bead << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM << "GB" << '\n';
  call << "%NprocShared=" << Ncpus << '\n';
  call << "%NoSave" << '\n'; //Deletes files
  call << "#P " << "external=\"FLUKE -GauExtern ";
  call << "QMMMExt"; //Just the stub
  call << "_" << Bead;
  call << " -n " << Ncpus;
  call << " -c " << confilename;
  call << " -r " << regfilename;
  call << "\"" << '\n';
  call << "Symmetry=None Opt=(MaxCycles=";
  call << QMMMOpts.MaxOptSteps;
  call << ",MaxStep=" << int(round((QMMMOpts.MaxStep/(0.01*BohrRad))));
  call << ")" << '\n';
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.Charge << " " << QMMMOpts.Spin << '\n';
  //Add atoms
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion == 1)
    {
      call << Struct[i].QMTyp;
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
    if (Struct[i].PAregion == 1)
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
  call << "QMMMExt_" << Bead;
  sys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "QMMMExt_";
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
          if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
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
  //Clean up files
  call.str("");
  call << "rm -f QMMMExt_";
  call << Bead << ".*";
  sys = system(call.str().c_str());
  //Calculate new point-charges
  GaussianCharges(Struct,QMMMOpts,Bead);
  if (Optfinished == 0)
  {
    cout << "Warning: Optimization did not converge!!!";
    cout << '\n';
    cout << "An older geometry will be recovered...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
  }
  return E;
};

