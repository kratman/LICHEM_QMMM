/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for Gaussian. These routines are written for g09.
 Note that the external function needs to interface with all MM codes.

*/

//QM utility functions
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
    call << "Charge=angstroms"; //Read charges
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
  }
  Eqm -= Eself;
  QMlog.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM";
  call << "_" << Bead;
  call << ".*";
  sys = system(call.str().c_str());
  //Return
  return Eqm;
};

void ExternalGaussian(int& argc, char**& argv)
{
  //This function is an "external script" that can be called by
  //Gaussian's external interface
  double Eqm = 0; //Stores partial energies
  double Emm = 0; //Stores partial energies
  double Eself = 0; //Stores partial energies
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
  //Construct g09 input
  call.str("");
  call << Stub << "_extern.com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=" << Stub;
  call << "_extern.chk";
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
    call << "Charge=angstroms"; //Read charges
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
      call << " " << setprecision(12) << Struct[i].P[0].x;
      call << " " << setprecision(12) << Struct[i].P[0].y;
      call << " " << setprecision(12) << Struct[i].P[0].z;
      call.copyfmt(cout);
      call << '\n';
    }
    if (Struct[i].PAregion == 1)
    {
      call << "F";
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[0].x;
      call << " " << setprecision(12) << Struct[i].P[0].y;
      call << " " << setprecision(12) << Struct[i].P[0].z;
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
        call << " " << setprecision(12) << Struct[i].P[0].x;
        call << " " << setprecision(12) << Struct[i].P[0].y;
        call << " " << setprecision(12) << Struct[i].P[0].z;
        call << " " << setprecision(12) << Struct[i].MP[0].q;
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
  //Construct MM input
  if (TINKER == 1)
  {
    //Forces input for TINKER
    call.str("");
    call << "cp " << TINKKeyFile << " ";
    call << Stub << "_extern.key";
    sys = system(call.str().c_str());
    //Save new keyfile name
    call.str("");
    call << Stub << "_extern.key";
    TINKKeyFile = call.str(); //Save the new name
    //Add QM atoms to force field parameters list
    ofile.open(TINKKeyFile.c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        if (ct == 0)
        {
          //Start a new active line
          ofile << "active ";
        }
        else
        {
          //Place a space to separate values
          ofile << " ";
        }
        ofile << (Struct[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          ofile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      ofile << '\n';
    }
    ofile << "group-inter" << '\n'; //Modify interactions
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add group 1 atoms
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        if (ct == 0)
        {
          //Start a new group line
          ofile << "group 1 ";
        }
        else
        {
          //Place a space to separate values
          ofile << " ";
        }
        ofile << (Struct[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate a group line
          ct = 0;
          ofile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing group line
      ofile << '\n';
    }
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add group 1 atoms
      if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1))
      {
        if (ct == 0)
        {
          //Start a new group line
          ofile << "group 2 ";
        }
        else
        {
          //Place a space to separate values
          ofile << " ";
        }
        ofile << (Struct[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate a group line
          ct = 0;
          ofile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing group line
      ofile << '\n';
    }
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      //Add atom types
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        ofile << "atom " << (MaxTINKERNum+ct) << " ";
        ofile << Struct[i].NumClass << " ";
        ofile << Struct[i].MMTyp << " ";
        ofile << "\"Dummy QM atom type\" ";
        ofile << RevTyping(Struct[i].QMTyp) << " ";
        ofile << Struct[i].m << " ";
        ofile << Struct[i].Bonds.size();
        ofile << '\n';
        ct += 1;
      }
    }
    if (CHRG == 1)
    {
      ct = 0;
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
        {
          ofile << "charge " << (MaxTINKERNum+ct) << " ";
          ofile << 0.0; //Delete charges
          ofile << '\n';
          ct += 1;
        }
      }
    }
    ofile.flush();
    ofile.close();
    //Create TINKER xyz file from the structure
    call.str("");
    call << Stub << "_extern.xyz";
    ofile.open(call.str().c_str(),ios_base::out);
    //Write atoms to the xyz file
    ofile << Natoms << '\n';
    if (PBCon == 1)
    {
      //Write box size
      ofile << Lx << " " << Ly << " " << Lz;
      ofile << " 90.0 90.0 90.0";
      ofile << '\n';
    }
    ct = 0; //Counter for QM atoms
    for (int i=0;i<Natoms;i++)
    {
      ofile << setw(6) << (Struct[i].id+1);
      ofile << " ";
      ofile << setw(3) << Struct[i].MMTyp;
      ofile << " ";
      ofile << setw(12) << Struct[i].P[0].x;
      ofile << " ";
      ofile << setw(12) << Struct[i].P[0].y;
      ofile << " ";
      ofile << setw(12) << Struct[i].P[0].z;
      ofile << " ";
      if ((Struct[i].QMregion != 1) and (Struct[i].PAregion != 1))
      {
        ofile << setw(4) << Struct[i].NumTyp;
      }
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        ofile << setw(4) << (MaxTINKERNum+ct);
        ct += 1; //Count number of qm atoms
      }
      for (int j=0;j<Struct[i].Bonds.size();j++)
      {
        ofile << " "; //Avoids trailing spaces
        ofile << setw(6) << (Struct[i].Bonds[j]+1);
      }
      ofile.copyfmt(cout);
      ofile << '\n';
    }
    ofile.flush();
    ofile.close();
  }
  //Run g09 and MM
  if (TINKER == 1)
  {
    call.str("");
    call << "testgrad " << Stub;
    call << "_extern.xyz Y N > QMMM_extern.grad";
    sys = system(call.str().c_str());
  }
  call.str("");
  call << "g09 " << Stub;
  call << "_extern";
  sys = system(call.str().c_str());
  //Read the QMMM forces
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
    //Open files
    call.str("");
    call << Stub << "_extern.grad";
    MMgrad.open(call.str().c_str(),ios_base::in);
    //Read derivatives
    bool GradDone = 0;
    while ((!MMgrad.eof()) and (!GradDone))
    {
      getline(MMgrad,dummy);
      stringstream line(dummy);
      line >> dummy;
      if (dummy == "Type")
      {
        line >> dummy >> dummy;
        if (dummy == "dE/dX")
        {
          GradDone = 1; //Not grad school, that lasts forever
          getline(MMgrad,dummy);
          for (int i=0;i<(Nqm+Npseudo);i++)
          {
            //Convoluted, but "easy"
            getline(MMgrad,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy; //Clear junk
            line >> Forces[i].x;
            line >> Forces[i].y;
            line >> Forces[i].z;
            //Switch to a.u. and change sign
            Forces[i].x *= -1*(kcal2eV*BohrRad)/Har2eV;
            Forces[i].y *= -1*(kcal2eV*BohrRad)/Har2eV;
            Forces[i].z *= -1*(kcal2eV*BohrRad)/Har2eV;
          }
        }
      }
    }
  }
  MMgrad.close();
  //Save MM forces to log file
  call.str("");
  call << Stub << "_extern.grad";
  MMgrad.open(call.str().c_str(),ios_base::in);
  if (MMgrad.good())
  {
    call.str("");
    while (!MMgrad.eof())
    {
      //Copy log line by line
      getline(MMgrad,dummy);
      call << dummy << '\n';
    }
    MMgrad.close();
    GauMsg << call.str();
  }
  //QM forces
  call.str("");
  call << Stub << "_extern.log";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool GradDone = 0; //Check if the gradient finished
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
          Forces[i].x += Fx;
          Forces[i].y += Fy;
          Forces[i].z += Fz;
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
  }
  Eqm -= Eself;
  QMlog.close();
  //Write formatted output for g09
  double E = Eqm+Emm; //Calculate
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
    GauOutput << setw(20) << (-1*Forces[i].x);
    GauOutput << setw(20) << (-1*Forces[i].y);
    GauOutput << setw(20) << (-1*Forces[i].z);
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
  //Save QM force calculations to log file
  call.str("");
  call << Stub << "_extern.log";
  QMlog.open(call.str().c_str(),ios_base::in);
  if (QMlog.good())
  {
    call.str("");
    while (!QMlog.eof())
    {
      //Copy log line by line
      getline(QMlog,dummy);
      call << dummy << '\n';
    }
    QMlog.close();
    GauMsg << call.str();
  }
  GauMsg.flush();
  GauMsg.close();
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
  //Clean up output external
  call.str("");
  call << "rm -f ";
  call << Stub << "_extern.*";
  sys = system(call.str().c_str());
  //Return to Gaussian
  exit(0);
  return;
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
  call << ".com QMMM";
  call << "_" << Bead;
  call << ".log";
  sys = system(call.str().c_str());
  return;
};

//QM wrappers
double GaussianWrapper(string RunTyp, vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs Gaussian
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0; //QM energy
  double Eself = 0.0; //Field self-energy
  int sys;
  call.str("");
  if (RunTyp == "Opt")
  {
    //Write a new XYZ
    call.str("");
    call << "QMMM_" << Bead << ".xyz";
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
    //Write Gaussian input
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
    call << "#P " << "external=\"FLUKE -GauExtern ";
    call << "QMMM"; //Just the stub
    call << "_" << Bead;
    call << " -n " << Ncpus;
    call << " -c " << confilename;
    call << " -r " << regfilename;
    call << "\"" << '\n';
    call << "Symmetry=None Opt=(MaxCycles=";
    call << QMMMOpts.MaxOptSteps;
    call << ",MaxStep=" << (30*QMMMOpts.SteepStep);
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
    call << "QMMM_" << Bead;
    sys = system(call.str().c_str());
    //Read new structure
    call.str("");
    call << "QMMM";
    call << "_" << Bead;
    call << ".log";
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
    call << "rm -f QMMM";
    call << "_" << Bead;
    call << ".com ";
    call << "QMMM";
    call << "_" << Bead;
    call << ".log";
    sys = system(call.str().c_str());
    //Calculate new point-charges
    GaussianCharges(Struct,QMMMOpts,Bead);
    if (Optfinished == 0)
    {
      cout << "Warning: Optimization did not converge!!!";
      cout << '\n';
      E = 10000.0; //Large number to reject step
    }
  }
  if (RunTyp == "Enrg")
  {
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
      cout << " The calculation attempt will continue...";
      cout << '\n';
      E = 10000.0; //Large number to reject step
    }
  }
  //Remove files
  call.str("");
  call << "rm -f QMMM";
  call << "_" << Bead;
  call << ".com QMMM";
  call << "_" << Bead;
  call << ".log";
  sys = system(call.str().c_str());
  //Change units
  E -= Eself;
  E *= Har2eV;
  return E;
};

