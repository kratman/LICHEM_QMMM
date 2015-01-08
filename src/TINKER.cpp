/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for TINKER.

*/

//MM utility functions
double TINKERPolEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Function to extract the polarization energy
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  string TINKKeyFile = "tinker.key";
  int MaxTINKERNum = 3500;
  int MaxTINKERClass = 100;
  double Epol = 0;
  double E = 0;
  int sys;
  int ct;
  //Create TINKER xyz file
  call.str("");
  call << "QMMM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
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
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    if (Struct[i].QMregion != 1)
    {
      ofile << setw(4) << Struct[i].NumTyp;
    }
    if (Struct[i].QMregion == 1)
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
  //Create new TINKER key file
  call.str("");
  call << "cp " << TINKKeyFile << " QMMM";
  call << "_" << Bead;
  call << ".key";
  sys = system(call.str().c_str());
  //Save new keyfile name
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".key";
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
    if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1))
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
  ct = 0;
  for (int i=0;i<Natoms;i++)
  {
    //Add atom types
    if (Struct[i].QMregion == 1)
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
  ct = 0;
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if ((Struct[i].PAregion == 1) or (Struct[i].QMregion == 1))
    {
      //Write new multipole definition for the atom ID
      WriteTINKMpole(Struct,ofile,i,Bead);
    }
  }
  ofile.flush();
  ofile.close();
  //Calculate QMMM energy
  call.str("");
  call << "analyze QMMM_";
  call << Bead << ".xyz E > QMMM_";
  call << Bead << ".log";
  sys = system(call.str().c_str());
  //Extract polarization energy
  call.str("");
  call << "QMMM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  bool Efound = 0;
  while (!ifile.eof())
  {
    ifile >> dummy;
    if (dummy == "Total")
    {
      ifile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        ifile >> dummy >> E;
        Efound = 1;
      }
    }
    if (dummy == "Polarization")
    {
      ifile >> Epol;
    }
  }
  if (Efound == 0)
  {
    //Warn user if no energy was found
    cout << "Warning: No MM energy found after a calculation!!!";
    cout << '\n';
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  //Return polarization energy in kcal/mol
  return Epol;
};

double TINKERForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  string TINKKeyFile = "tinker.key";
  int MaxTINKERNum = 3500;
  int MaxTINKERClass = 100;
  double Emm = 0.0;
  int ct;
  int sys;
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp " << TINKKeyFile << " ";
  call << "QMMM";
  call << "_" << Bead;
  call << ".key";
  sys = system(call.str().c_str());
  //Save new keyfile name
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".key";
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
  if (AMOEBA == 1)
  {
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        WriteTINKMpole(Struct,ofile,i,Bead);
      }
    }
  }
  ofile.flush();
  ofile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".xyz";
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
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
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
  //Run MM
  call.str("");
  call << "testgrad " << "QMMM";
  call << "_" << Bead;
  call << ".xyz Y N > QMMM";
  call << "_" << Bead;
  call << ".grad";
  sys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".grad";
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
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Convoluted, but "easy"
          getline(MMgrad,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> Fx;
          line >> Fy;
          line >> Fz;
          //Switch to a.u. and change sign
          Forces[i].x += -1*Fx*kcal2eV;
          Forces[i].y += -1*Fy*kcal2eV;
          Forces[i].z += -1*Fz*kcal2eV;
        }
      }
    }
  }
  MMgrad.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM";
  call << "_" << Bead;
  call << ".*";
  sys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

void FindTINKERClasses(vector<QMMMAtom>& Struct)
{
  //Parses TINKER parameter files to find atom classes
  fstream ifile;
  string dummy;
  string TINKKeyFile = "tinker.key";
  ifile.open(TINKKeyFile.c_str(),ios_base::in);
  if (!ifile.good())
  {
    //Exit if files do not exist
    cout << "Error: Missing tinker.key file.";
    cout << endl;
    exit(0);
  }
  bool FileFound = 0; //Bool to break loops
  while ((!ifile.eof()) and (!FileFound))
  {
    //Detect the name of the force field file
    ifile >> dummy;
    if (dummy == "parameters")
    {
      ifile >> dummy;
      FileFound = 1;
    }
  }
  ifile.close();
  ifile.open(dummy.c_str(),ios_base::in);
  if (!FileFound)
  {
    //Exit if parameter file is not found
    cout << "Error: Cannot find TINKER parameter file.";
    cout << endl;
    exit(0);
  }
  if (!ifile.good())
  {
    //Exit if parameter file does not exist
    cout << "Error: Cannot read TINKER ";
    cout << dummy;
    cout << " parameter file.";
    cout << endl;
    exit(0);
  }
  int ct = 0; //Generic counter
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream FullLine(dummy);
    FullLine >> dummy;
    if (dummy == "atom")
    {
      int AtType,AtClass;
      FullLine >> AtType;
      FullLine >> AtClass;
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].NumTyp == AtType)
        {
          Struct[i].NumClass = AtClass;
          ct += 1;
        }
      }
    }
  }
  if (ct < Natoms)
  {
    cout << "Error: Atom type not found in TINKER parameters.";
    cout << '\n';
    cout << " Please check the input.";
    cout << endl;
    exit(0);
  }
  return;
}

//MM wrappers
double TINKERWrapper(string RunTyp, vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs TINKER MM
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy;
  string TINKKeyFile = "tinker.key";
  int MaxTINKERNum = 3500;
  int MaxTINKERClass = 100;
  double Epol = 0;
  double E = 0;
  int sys;
  int ct;
  call.str("");
  //Copy the original key file and make changes
  if (QMMM == 1)
  {
    call.str("");
    call << "cp " << TINKKeyFile << " QMMM";
    call << "_" << Bead;
    call << ".key";
    sys = system(call.str().c_str());
    //Save new keyfile name
    call.str("");
    call << "QMMM";
    call << "_" << Bead;
    call << ".key";
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
      if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1)
         or ((Struct[i].PAregion == 1) and (RunTyp == "Opt")))
      {
        if ((Struct[i].Frozen == 0) or (RunTyp == "Enrg"))
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
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      ofile << '\n';
    }
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      //Add atom types
      if (Struct[i].QMregion == 1)
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
        if (Struct[i].QMregion == 1)
        {
          //New atom types are only needed for QM atoms
          ofile << "charge " << (MaxTINKERNum+ct) << " ";
          if (RunTyp == "Opt")
          {
            ofile << Struct[i].MP[Bead].q;
          }
          if (RunTyp == "Enrg")
          {
            ofile << 0.0;
          }
          ofile << '\n';
          ct += 1;
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].PAregion == 1)
        {
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
        }
        if (Struct[i].QMregion == 1)
        {
          double qi = 0;
          if (RunTyp == "Enrg")
          {
            //remove charge
            qi = Struct[i].MP[Bead].q;
            Struct[i].MP[Bead].q = 0;
          }
          WriteTINKMpole(Struct,ofile,i,Bead);
          Struct[i].MP[Bead].q += qi; //Restore charge
        }
      }
    }
    ofile.flush();
    ofile.close();
  }
  //Create TINKER xyz file from the structure
  call.str("");
  call << "QMMM_" << Bead << ".xyz";
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
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    if (Struct[i].QMregion != 1)
    {
      ofile << setw(4) << Struct[i].NumTyp;
    }
    if (Struct[i].QMregion == 1)
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
  //Run optimization
  if (RunTyp == "Opt")
  {
    //Run TINKER
    call.str("");
    call << "newton QMMM_";
    call << Bead << ".xyz A A 0.01 > QMMM_";
    call << Bead << ".log";
    sys = system(call.str().c_str());
    //Read new structure
    call.str("");
    call << "QMMM_" << Bead << ".xyz_2";
    ifile.open(call.str().c_str(),ios_base::in);
    getline(ifile,dummy);
    if (PBCon == 1)
    {
      getline(ifile,dummy);
    }
    for (int i=0;i<Natoms;i++)
    {
      getline(ifile,dummy);
      stringstream line(dummy);
      //Read new positions
      line >> dummy >> dummy;
      line >> Struct[i].P[Bead].x;
      line >> Struct[i].P[Bead].y;
      line >> Struct[i].P[Bead].z;
    }
    ifile.close();
  }
  //Calculate MM potential energy
  if (RunTyp == "Enrg")
  {
    //Run TINKER
    call.str("");
    call << "analyze QMMM_";
    call << Bead << ".xyz E > QMMM_";
    call << Bead << ".log";
    sys = system(call.str().c_str());
    call.str("");
    call << "QMMM_" << Bead << ".log";
    ifile.open(call.str().c_str(),ios_base::in);
    //Read MM potential energy
    bool Efound = 0;
    while (!ifile.eof())
    {
      ifile >> dummy;
      if (dummy == "Total")
      {
        ifile >> dummy >> dummy;
        if (dummy == "Energy")
        {
          ifile >> dummy >> E;
          Efound = 1;
        }
      }
      if (dummy == "Polarization")
      {
        ifile >> Epol;
      }
    }
    if (Efound == 0)
    {
      //Warn user if no energy was found
      cout << "Warning: No MM energy found after a calculation!!!";
      cout << '\n';
      cout << " FLUKE will attempt to continue...";
      cout << '\n';
      E = HugeNum; //Large number to reject step
    }
    ifile.close();
  }
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".xyz_*";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  //Calculate polarization energy
  if ((AMOEBA == 1) and (QMMM == 1))
  {
    //Correct polarization energy for QMMM simulations
    E -= Epol; //Incorrect polarization energy
    E += TINKERPolEnergy(Struct,QMMMOpts,Bead);
  }
  //Change units
  E *= kcal2eV;
  return E;
};

