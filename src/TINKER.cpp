/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for TINKER.

 Reference for TINKER:
 Ponder, TINKER - Software Tools for Molecular Design

*/

//MM utility functions
void FindTINKERClasses(vector<QMMMAtom>& Struct)
{
  //Parses TINKER parameter files to find atom classes
  fstream ifile;
  string dummy; //Generic string
  ifile.open("tinker.key",ios_base::in);
  if (!ifile.good())
  {
    //Exit if files do not exist
    cout << "Error: Missing tinker.key file.";
    cout << '\n';
    cout.flush();
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
    cout << '\n';
    cout.flush();
    exit(0);
  }
  if (!ifile.good())
  {
    //Exit if parameter file does not exist
    cout << "Error: Cannot read TINKER ";
    cout << dummy;
    cout << " parameter file.";
    cout << '\n';
    cout.flush();
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
    cout << '\n';
    cout.flush();
    exit(0);
  }
  return;
}

//MM wrapper functions
void TINKERInduced(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Function to extract induced dipoles
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  int sys; //Dummy return for system calls
  int ct; //Generic counter
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  call << "cp tinker.key QMMM_";
  call << Bead << ".key";
  sys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "QMMM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n'; //Make sure current line is empty
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "save-induced" << '\n'; //Save induced dipoles
  ofile << "thermostat berendsen" << '\n';
  ofile << "tau-temperature 0.1" << '\n';
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (Struct[i].MMregion or Struct[i].BAregion)
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
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if (Struct[i].QMregion or Struct[i].PAregion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMpole(Struct,ofile,i,Bead);
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
  }
  ofile.flush();
  ofile.close();
  //Calculate induced dipoles using dynamic
  call.str("");
  call << "dynamic QMMM_" << Bead << ".xyz ";
  call << "1 1e-4 1e-7 2 0 > QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  //Extract induced dipoles from the MD cycle file
  call.str("");
  call << "QMMM_" << Bead << ".001u";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy); //Clear number of atoms
  while (ifile.good())
  {
    int AtNum; //Identifies which atom was polarized
    //Parse file line by line
    getline(ifile,dummy);
    stringstream line(dummy);
    //Save dipoles for later
    line >> AtNum >> dummy; //Collect atom number and clear junk
    if (line.good())
    {
      AtNum -= 1; //Fixes array indexing
      line >> Struct[AtNum].MP[Bead].IDx;
      line >> Struct[AtNum].MP[Bead].IDy;
      line >> Struct[AtNum].MP[Bead].IDz;
      //Change units from Debye to a.u.
      Struct[AtNum].MP[Bead].IDx *= Debye2au;
      Struct[AtNum].MP[Bead].IDy *= Debye2au;
      Struct[AtNum].MP[Bead].IDz *= Debye2au;
    }
  }
  ifile.close();
  //Delete junk files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".xyz ";
  call << "QMMM_" << Bead << ".key ";
  call << "QMMM_" << Bead << ".0* ";
  call << "QMMM_" << Bead << ".dyn ";
  call << "QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  return;
};

double TINKERPolEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Function to extract the polarization energy
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  double Epol = 0;
  double E = 0;
  int sys; //Dummy return for system calls
  int ct; //Generic counter
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  call << "cp tinker.key QMMM_";
  call << Bead << ".key";
  sys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "QMMM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "polarizeterm only" << '\n'; //Get rid of other interactions
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (Struct[i].MMregion or Struct[i].BAregion)
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
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if (Struct[i].PAregion or Struct[i].QMregion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMpole(Struct,ofile,i,Bead);
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
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
  if (!Efound)
  {
    //Warn user if no energy was found
    cout << "Warning: No MM energy found after a calculation!!!";
    cout << '\n';
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
    cout.flush(); //Print warning immediately
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
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double Emm = 0.0;
  int ct; //Generic counter
  int sys; //Dummy return for system calls
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key QMMM_";
  call << Bead << ".key";
  sys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "QMMM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (Struct[i].QMregion or Struct[i].PAregion)
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
    if (Struct[i].QMregion or Struct[i].PAregion)
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
  if (CHRG == 1)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion or Struct[i].PAregion)
      {
        //New charges are needed for QM and PA atoms
        ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
        ofile << "0.0"; //Delete charges
        ofile << '\n';
      }
    }
  }
  if (AMOEBA == 1)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion or Struct[i].PAregion)
      {
        double qi = 0;
        //remove charge
        qi = Struct[i].MP[Bead].q;
        Struct[i].MP[Bead].q = 0;
        WriteTINKMpole(Struct,ofile,i,Bead);
        Struct[i].MP[Bead].q += qi; //Restore charge
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  call << ".xyz Y N N > QMMM";
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
          //Switch to eV/A and change sign
          if (abs(Fx) >= 1e-4)
          {
            Forces[i].x += -1*Fx*kcal2eV;
          }
          if (abs(Fy) >= 1e-4)
          {
            Forces[i].y += -1*Fy*kcal2eV;
          }
          if (abs(Fz) >= 1e-4)
          {
            Forces[i].z += -1*Fz*kcal2eV;
          }
        }
      }
    }
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "Energy")
      {
        //Collect partial MM energy
        line >> dummy >> Emm;
      }
    }
  }
  MMgrad.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".key";
  call << " QMMM_" << Bead << ".grad";
  sys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKERPolForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double Emm = 0.0;
  int ct; //Generic counter
  int sys; //Dummy return for system calls
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key QMMM_";
  call << Bead << ".key";
  sys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "QMMM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "polarizeterm only" << '\n';
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (Struct[i].QMregion or Struct[i].PAregion)
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
  if (AMOEBA == 1)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion or Struct[i].PAregion)
      {
        WriteTINKMpole(Struct,ofile,i,Bead);
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  call << ".xyz Y N N > QMMM";
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
          //Switch to eV/A and change sign
          if (abs(Fx) >= 1e-4)
          {
            Forces[i].x += -1*Fx*kcal2eV;
          }
          if (abs(Fy) >= 1e-4)
          {
            Forces[i].y += -1*Fy*kcal2eV;
          }
          if (abs(Fz) >= 1e-4)
          {
            Forces[i].z += -1*Fz*kcal2eV;
          }
        }
      }
    }
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "Energy")
      {
        //Collect partial MM energy
        line >> dummy >> Emm;
      }
    }
  }
  MMgrad.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".key";
  call << " QMMM_" << Bead << ".grad";
  sys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKERMMForces(vector<QMMMAtom>& Struct, vector<Coord>& MMForces,
     QMMMSettings& QMMMOpts, int Bead)
{
  //A routine to extract MM forces from TINKER
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double Emm = 0.0;
  int ct; //Generic counter
  int sys; //Dummy return for system calls
  //Copy the original key file and make changes
  if (QMMM)
  {
    call.str("");
    call << "cp tinker.key QMMM_";
    call << Bead << ".key";
    sys = system(call.str().c_str());
    //Update key file
    call.str("");
    call << "QMMM_";
    call << Bead << ".key";
    ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (Struct[i].MMregion or Struct[i].BAregion or Struct[i].PAregion)
      {
        if (!Struct[i].Frozen)
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
    if (CHRG == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion)
        {
          //New charges are only needed for QM atoms
          ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
          ofile << Struct[i].MP[Bead].q;
          ofile << '\n';
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion or Struct[i].PAregion)
        {
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  call << ".xyz Y N N > QMMM";
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
        for (int i=0;i<Natoms;i++)
        {
          if ((!Struct[i].Frozen) and (!Struct[i].QMregion))
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
            //Switch to eV/A and change sign
            if (abs(Fx) >= 1e-4)
            {
              MMForces[i].x += -1*Fx*kcal2eV;
            }
            if (abs(Fy) >= 1e-4)
            {
              MMForces[i].y += -1*Fy*kcal2eV;
            }
            if (abs(Fz) >= 1e-4)
            {
              MMForces[i].z += -1*Fz*kcal2eV;
            }
          }
        }
      }
    }
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "Energy")
      {
        //Collect partial MM energy
        line >> dummy >> Emm;
      }
    }
  }
  MMgrad.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".key";
  call << " QMMM_" << Bead << ".grad";
  sys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKEREnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs TINKER MM energy calculations
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  double E = 0;
  int sys; //Dummy return for system calls
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  if (QMMM)
  {
    call.str("");
    call << "cp tinker.key QMMM_";
    call << Bead << ".key";
    sys = system(call.str().c_str());
    //Update key file
    call.str("");
    call << "QMMM_";
    call << Bead << ".key";
    ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ofile << "polarizeterm none" << '\n'; //Remove polarization energy
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (Struct[i].MMregion or Struct[i].BAregion)
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
    if (CHRG == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion or Struct[i].PAregion or Struct[i].BAregion)
        {
          //New charges are only needed for QM atoms
          ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
          ofile << "0.0" << '\n';
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add multipoles
        if (Struct[i].QMregion or Struct[i].PAregion or Struct[i].BAregion)
        {
          double qi = 0;
          //remove charge
          qi = Struct[i].MP[Bead].q;
          Struct[i].MP[Bead].q = 0;
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
          //Restore charge
          Struct[i].MP[Bead].q = qi;
          //Remove polarization
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  //Calculate MM potential energy
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
  }
  if (!Efound)
  {
    //Warn user if no energy was found
    cout << "Warning: No MM energy found after a calculation!!!";
    cout << '\n';
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
    cout.flush(); //Print warning immediately
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  //Calculate polarization energy
  if ((AMOEBA == 1) and QMMM)
  {
    //Correct polarization energy for QMMM simulations
    E += TINKERPolEnergy(Struct,QMMMOpts,Bead);
  }
  //Change units
  E *= kcal2eV;
  return E;
};

void TINKERDynamics(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Runs TINKER MD (for MM atoms only)
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  double Epol = 0;
  double E = 0;
  int sys; //Dummy return for system calls
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  if (QMMM)
  {
    call.str("");
    call << "cp tinker.key QMMM_";
    call << Bead << ".key";
    sys = system(call.str().c_str());
    //Update key file
    call.str("");
    call << "QMMM_";
    call << Bead << ".key";
    ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ofile << "thermostat berendsen";
    ofile << '\n';
    ofile << "tau-temperature ";
    ofile << QMMMOpts.tautemp << '\n';
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (Struct[i].MMregion or Struct[i].BAregion)
      {
        if (!Struct[i].Frozen)
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
    if (CHRG == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion)
        {
          //New charges are only needed for QM atoms
          ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
          ofile << Struct[i].MP[Bead].q;
          ofile << '\n';
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion or Struct[i].PAregion)
        {
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  call.str("");
  call << "dynamic ";
  call << "QMMM_" << Bead << ".xyz ";
  call << QMMMOpts.Nsteps << " ";
  call << QMMMOpts.dt << " ";
  call << (QMMMOpts.Nsteps*QMMMOpts.dt/1000) << " ";
  call << "2 " << QMMMOpts.Temp;
  call << " > QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "QMMM_" << Bead << ".001";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy); //Discard number of atoms
  if (PBCon == 1)
  {
    //Discard PBC information
    getline(ifile,dummy);
  }
  for (int i=0;i<Natoms;i++)
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    //Read new positions
    line >> dummy >> dummy; //Discard atom ID and type
    line >> Struct[i].P[Bead].x;
    line >> Struct[i].P[Bead].y;
    line >> Struct[i].P[Bead].z;
  }
  ifile.close();
  //Clean up all files except the .dyn files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".0*";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  return;
};

double TINKEROpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs TINKER MM optimization
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  double Epol = 0;
  double E = 0;
  int sys; //Dummy return for system calls
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  if (QMMM)
  {
    call.str("");
    call << "cp tinker.key QMMM_";
    call << Bead << ".key";
    sys = system(call.str().c_str());
    //Update key file
    call.str("");
    call << "QMMM_";
    call << Bead << ".key";
    ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (Struct[i].MMregion or Struct[i].BAregion)
      {
        if (!Struct[i].Frozen)
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
    if (CHRG == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion)
        {
          //New charges are only needed for QM atoms
          ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
          ofile << Struct[i].MP[Bead].q;
          ofile << '\n';
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion or Struct[i].PAregion)
        {
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
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
    ofile << setw(4) << Struct[i].NumTyp;
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
  call.str("");
  call << "minimize QMMM_";
  call << Bead << ".xyz ";
  call << QMMMOpts.MMOptTol << " > QMMM_";
  call << Bead << ".log";
  sys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "QMMM_" << Bead << ".xyz_2";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy); //Discard number of atoms
  if (PBCon == 1)
  {
    //Discard PBC information
    getline(ifile,dummy);
  }
  for (int i=0;i<Natoms;i++)
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    //Read new positions
    line >> dummy >> dummy; //Discard atom ID and type
    line >> Struct[i].P[Bead].x;
    line >> Struct[i].P[Bead].y;
    line >> Struct[i].P[Bead].z;
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".xyz_*";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  //Change units
  E *= kcal2eV;
  return E;
};
