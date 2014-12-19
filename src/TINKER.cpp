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
double TinkerForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  string TinkKeyFile = "tinker.key";
  int MaxTinkerNum = 3500;
  int MaxTinkerClass = 100;
  double Emm = 0.0;
  int ct;
  int sys;
  call.str("");
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp " << TinkKeyFile << " ";
  call << "QMMM";
  if (Bead != -1)
  {
    call << "_" << Bead;
  }
  call << ".key";
  sys = system(call.str().c_str());
  //Save new keyfile name
  call.str("");
  call << "QMMM";
  if (Bead != -1)
  {
    call << "_" << Bead;
  }
  call << ".key";
  TinkKeyFile = call.str(); //Save the new name
  //Add QM atoms to force field parameters list
  ofile.open(TinkKeyFile.c_str(),ios_base::app|ios_base::out);
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
      ofile << "atom " << (MaxTinkerNum+ct) << " ";
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
        ofile << "charge " << (MaxTinkerNum+ct) << " ";
        ofile << 0.0; //Delete charges
        ofile << '\n';
        ct += 1;
      }
    }
  }
  ofile.flush();
  ofile.close();
  //Create Tinker xyz file from the structure
  call.str("");
  call << "QMMM";
  if (Bead != -1)
  {
    call << "_" << Bead;
  }
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
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    if (Bead == -1)
    {
      ofile << setw(12) << Struct[i].x;
      ofile << " ";
      ofile << setw(12) << Struct[i].y;
      ofile << " ";
      ofile << setw(12) << Struct[i].z;
    }
    if (Bead != -1)
    {
      ofile << setw(12) << Struct[i].P[Bead].x;
      ofile << " ";
      ofile << setw(12) << Struct[i].P[Bead].y;
      ofile << " ";
      ofile << setw(12) << Struct[i].P[Bead].z;
    }
    ofile << " ";
    if ((Struct[i].QMregion != 1) and (Struct[i].PAregion != 1))
    {
      ofile << setw(4) << Struct[i].NumTyp;
    }
    if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
    {
      ofile << setw(4) << (MaxTinkerNum+ct);
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
  if (Bead != -1)
  {
    call << "_" << Bead;
  }
  call << ".xyz Y N > QMMM";
  if (Bead != -1)
  {
    call << "_" << Bead;
  }
  call << ".grad";
  sys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "QMMM";
  if (Bead != -1)
  {
    call << "_" << Bead;
  }
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
  if (Bead != -1)
  {
    call << "_" << Bead;
  }
  call << ".*";
  sys = system(call.str().c_str());
  //Return
  return Emm;
};

void FindTinkerClasses(vector<QMMMAtom>& Struct)
{
  //Parses tinker parameter files to find atom classes
  fstream ifile;
  string dummy;
  string TinkKeyFile = "tinker.key";
  ifile.open(TinkKeyFile.c_str(),ios_base::in);
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
    cout << "Error: Cannot find tinker parameter file.";
    cout << endl;
    exit(0);
  }
  if (!ifile.good())
  {
    //Exit if parameter file does not exist
    cout << "Error: Cannot read tinker ";
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
    cout << "Error: Atom type not found in Tinker parameters.";
    cout << '\n';
    cout << " Please check the input.";
    cout << endl;
    exit(0);
  }
  return;
}

//MM wrappers
double TinkerWrapper(string RunTyp, vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs Tinker MM
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  string TinkKeyFile = "tinker.key";
  int MaxTinkerNum = 3500;
  int MaxTinkerClass = 100;
  double E = 0.0;
  int ct;
  int sys;
  call.str("");
  //Copy the original key file and make changes
  if (QMMM == 1)
  {
    call.str("");
    call << "cp " << TinkKeyFile << " QMMM";
    if (Bead != -1)
    {
      call << "_" << Bead;
    }
    call << ".key";
    sys = system(call.str().c_str());
    //Save new keyfile name
    call.str("");
    call << "QMMM";
    if (Bead != -1)
    {
      call << "_" << Bead;
    }
    call << ".key";
    TinkKeyFile = call.str(); //Save the new name
    //Add QM atoms to force field parameters list
    ofile.open(TinkKeyFile.c_str(),ios_base::app|ios_base::out);
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
        ofile << "atom " << (MaxTinkerNum+ct) << " ";
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
          ofile << "charge " << (MaxTinkerNum+ct) << " ";
          if (RunTyp == "Opt")
          {
            ofile << Struct[i].q;
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
    ofile.flush();
    ofile.close();
  }
  //Create Tinker xyz file from the structure
  if (Bead == -1)
  {
    ofile.open("QMMM.xyz",ios_base::out);
  }
  if (Bead != -1)
  {
    call.str("");
    call << "QMMM_" << Bead << ".xyz";
    ofile.open(call.str().c_str(),ios_base::out);
  }
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
    if (Bead == -1)
    {
      ofile << setw(12) << Struct[i].x;
      ofile << " ";
      ofile << setw(12) << Struct[i].y;
      ofile << " ";
      ofile << setw(12) << Struct[i].z;
    }
    if (Bead != -1)
    {
      ofile << setw(12) << Struct[i].P[Bead].x;
      ofile << " ";
      ofile << setw(12) << Struct[i].P[Bead].y;
      ofile << " ";
      ofile << setw(12) << Struct[i].P[Bead].z;
    }
    ofile << " ";
    if (Struct[i].QMregion != 1)
    {
      ofile << setw(4) << Struct[i].NumTyp;
    }
    if (Struct[i].QMregion == 1)
    {
      ofile << setw(4) << (MaxTinkerNum+ct);
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
    //Run tinker
    call.str("");
    call << "newton QMMM.xyz A A 0.01 > QMMM.log";
    sys = system(call.str().c_str());
    //Read new structure
    fstream ifile;
    ifile.open("QMMM.xyz_2",ios_base::in);
    getline(ifile,dummy);
    if (PBCon == 1)
    {
      getline(ifile,dummy);
    }
    for (int i=0;i<Natoms;i++)
    {
      getline(ifile,dummy);
      stringstream line(dummy);
      if (Bead == -1)
      {
        //Read new positions
        line >> dummy >> dummy;
        line >> Struct[i].x;
        line >> Struct[i].y;
        line >> Struct[i].z;
      }
      if (Bead != -1)
      {
        //Read new positions
        line >> dummy >> dummy;
        line >> Struct[i].P[Bead].x;
        line >> Struct[i].P[Bead].y;
        line >> Struct[i].P[Bead].z;
      }
    }
    ifile.close();
  }
  //Calculate MM potential energy
  if (RunTyp == "Enrg")
  {
    //Run tinker
    if (Bead != -1)
    {
      call.str("");
      call << "analyze QMMM_";
      call << Bead << ".xyz E > QMMM_";
      call << Bead << ".log";
      sys = system(call.str().c_str());
      call.str("");
      call << "QMMM_" << Bead << ".log";
      ifile.open(call.str().c_str(),ios_base::in);
    }
    if (Bead == -1)
    {
      sys = system("analyze QMMM.xyz E > QMMM.log");
      ifile.open("QMMM.log",ios_base::in);
    }
    //Read MM potential energy
    bool contread = 1;
    while ((!ifile.eof()) and (contread == 1))
    {
      ifile >> dummy;
      if (dummy == "Total")
      {
        ifile >> dummy >> dummy;
        if (dummy == "Energy")
        {
          ifile >> dummy >> E;
          contread = 0;
        }
      }
    }
    if (contread == 1)
    {
      //Warn user if no energy was found
      cout << "Warning: No MM energy found after calculation!!!";
      cout << '\n';
      cout << " The calculation attempt will continue...";
      cout << '\n';
      E = 10000.0; //Large number to reject step
    }
    ifile.close();
  }
  //Clean up files
  call.str("");
  call << "rm -f";
  if (Bead != -1)
  {
    call << " QMMM_" << Bead << ".xyz";
    call << " QMMM_" << Bead << ".log";
    call << " QMMM_" << Bead << ".xyz_*";
    call << " QMMM_" << Bead << ".key";
  }
  if (Bead == -1)
  {
    call << " QMMM.xyz";
    call << " QMMM.log";
    call << " QMMM.xyz_*";
    call << " QMMM.key";
  }
  sys = system(call.str().c_str());
  //Change units
  E *= kcal2eV;
  return E;
};

