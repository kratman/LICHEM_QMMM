/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for TINKER.

 Reference for TINKER:
 Ponder, TINKER - Software Tools for Molecular Design

*/

//MM utility functions
void FindTINKERClasses(vector<QMMMAtom>& Struct)
{
  //Parses TINKER parameter files to find atom classes
  fstream ifile; //Generic file stream
  string dummy; //Generic string
  int ct; //Generic counter
  //Open generic key file
  ifile.open("tinker.key",ios_base::in);
  if (!ifile.good())
  {
    //Exit if files do not exist
    cout << "Error: Missing tinker.key file.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  //Find the parameter file
  bool FileFound = 0; //Bool to break loops
  while ((!ifile.eof()) and (!FileFound))
  {
    //Detect the name of the force field file
    ifile >> dummy;
    LICHEMLowerText(dummy);
    if (dummy == "parameters")
    {
      ifile >> dummy;
      FileFound = 1;
    }
  }
  //Open the parameters
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
  //Find parameters for all atoms
  ct = 0; //Generic counter
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream FullLine(dummy);
    FullLine >> dummy;
    LICHEMLowerText(dummy);
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
  //Check for errors
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
  fstream ofile,ifile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  int ct; //Generic counter
  //Create TINKER xyz file
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    ofile << LICHEMFormDouble(Lx,12) << " ";
    ofile << LICHEMFormDouble(Ly,12) << " ";
    ofile << LICHEMFormDouble(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Create new TINKER key file
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n'; //Make sure current line is empty
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    ofile << "a-axis " << LICHEMFormDouble(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormDouble(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormDouble(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }
  ofile << "save-induced" << '\n'; //Save induced dipoles
  ofile << "thermostat berendsen" << '\n';
  ofile << "tau-temperature 0.1" << '\n';
  ct = 0; //Generic counter
  if (QMMM)
  {
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
  }
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if (Struct[i].QMregion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMpole(Struct,ofile,i,Bead);
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
    if (Struct[i].PBregion)
    {
      //Modify the charge to force charge balance with the boundaries
      double qi = Struct[i].MP[Bead].q; //Save a copy
      vector<int> Boundaries;
      Boundaries = TraceBoundary(Struct,i);
      double qnew = qi;
      for (unsigned int j=0;j<Boundaries.size();j++)
      {
        //Subtract boundary atom charge
        qnew -= Struct[Boundaries[j]].MP[Bead].q;
      }
      Struct[i].MP[Bead].q = qnew; //Save modified charge
      WriteTINKMpole(Struct,ofile,i,Bead);
      Struct[i].MP[Bead].q = qi; //Return to unmodified charge
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
    if (Struct[i].BAregion)
    {
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
  }
  ofile.flush();
  ofile.close();
  //Calculate induced dipoles using dynamic
  call.str("");
  call << "dynamic LICHM_" << Bead << ".xyz ";
  call << "1 1e-4 1e-7 2 0 > LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Extract induced dipoles from the MD cycle file
  call.str("");
  call << "LICHM_" << Bead << ".001u";
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
  call << "LICHM_" << Bead << ".xyz ";
  call << "LICHM_" << Bead << ".key ";
  call << "LICHM_" << Bead << ".0* ";
  call << "LICHM_" << Bead << ".dyn ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  return;
};

double TINKERPolEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                       int Bead)
{
  //Function to extract the polarization energy
  fstream ofile,ifile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double Epol = 0;
  double E = 0;
  int ct; //Generic counter
  //Create TINKER xyz file
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    ofile << LICHEMFormDouble(Lx,12) << " ";
    ofile << LICHEMFormDouble(Ly,12) << " ";
    ofile << LICHEMFormDouble(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Create new TINKER key file
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    ofile << "a-axis " << LICHEMFormDouble(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormDouble(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormDouble(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }
  ofile << "polarizeterm only" << '\n'; //Get rid of other interactions
  ct = 0; //Generic counter
  if (QMMM)
  {
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
  }
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if (Struct[i].QMregion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMpole(Struct,ofile,i,Bead);
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
    if (Struct[i].PBregion)
    {
      //Modify the charge to force charge balance with the boundaries
      double qi = Struct[i].MP[Bead].q; //Save a copy
      vector<int> Boundaries;
      Boundaries = TraceBoundary(Struct,i);
      double qnew = qi;
      for (unsigned int j=0;j<Boundaries.size();j++)
      {
        //Subtract boundary atom charge
        qnew -= Struct[Boundaries[j]].MP[Bead].q;
      }
      Struct[i].MP[Bead].q = qnew; //Save modified charge
      WriteTINKMpole(Struct,ofile,i,Bead);
      Struct[i].MP[Bead].q = qi; //Return to unmodified charge
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
    if (Struct[i].BAregion)
    {
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
  }
  ofile.flush();
  ofile.close();
  //Calculate QMMM energy
  call.str("");
  call << "analyze LICHM_";
  call << Bead << ".xyz E > LICHM_";
  call << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Extract polarization energy
  call.str("");
  call << "LICHM_" << Bead << ".log";
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
    cerr << "Warning: No MM energy found after a calculation!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    Epol = 0; //Prevents errors when polarization is off
    cerr.flush(); //Print warning immediately
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << Bead << ".xyz";
  call << " LICHM_" << Bead << ".log";
  call << " LICHM_" << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Return polarization energy in kcal/mol
  return Epol;
};

double TINKERForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                    QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double Emm = 0.0;
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    ofile << "a-axis " << LICHEMFormDouble(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormDouble(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormDouble(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (Struct[i].QMregion or Struct[i].PBregion)
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
    if (Struct[i].QMregion or Struct[i].PBregion)
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
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion or Struct[i].PBregion or Struct[i].BAregion)
      {
        //New charges are needed for QM and PB atoms
        ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
        ofile << "0.0"; //Delete charges
        ofile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion or Struct[i].PBregion or Struct[i].BAregion)
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
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    ofile << LICHEMFormDouble(Lx,12) << " ";
    ofile << LICHEMFormDouble(Ly,12) << " ";
    ofile << LICHEMFormDouble(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Run MM
  call.str("");
  call << "testgrad ";
  call << "LICHM_" << Bead << ".xyz";
  call << " Y N N > ";
  call << "LICHM_" << Bead << ".grad";
  GlobalSys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << Bead << ".grad";
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
          //Change from gradient to force
          Fx *= -1;
          Fy *= -1;
          Fz *= -1;
          //Switch to eV/A and save forces
          Forces(3*i) += Fx*kcal2eV;
          Forces(3*i+1) += Fy*kcal2eV;
          Forces(3*i+2) += Fz*kcal2eV;
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
  call << " LICHM_" << Bead << ".xyz";
  call << " LICHM_" << Bead << ".key";
  call << " LICHM_" << Bead << ".grad";
  GlobalSys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKERPolForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double Emm = 0.0;
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    ofile << "a-axis " << LICHEMFormDouble(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormDouble(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormDouble(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }
  ofile << "polarizeterm only" << '\n'; //Get rid of other interactions
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (Struct[i].QMregion or Struct[i].PBregion)
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
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion)
      {
        WriteTINKMpole(Struct,ofile,i,Bead);
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (Struct[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = Struct[i].MP[Bead].q; //Save a copy
        vector<int> Boundaries;
        Boundaries = TraceBoundary(Struct,i);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= Struct[Boundaries[j]].MP[Bead].q;
        }
        Struct[i].MP[Bead].q = qnew; //Save modified charge
        WriteTINKMpole(Struct,ofile,i,Bead);
        Struct[i].MP[Bead].q = qi; //Return to unmodified charge
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (Struct[i].BAregion)
      {
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
    }
  }
  ofile.flush();
  ofile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    ofile << LICHEMFormDouble(Lx,12) << " ";
    ofile << LICHEMFormDouble(Ly,12) << " ";
    ofile << LICHEMFormDouble(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Run MM
  call.str("");
  call << "testgrad ";
  call << "LICHM_" << Bead << ".xyz";
  call << " Y N N > ";
  call << "LICHM_" << Bead << ".grad";
  GlobalSys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << Bead << ".grad";
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
          //Change from gradient to force
          Fx *= -1;
          Fy *= -1;
          Fz *= -1;
          //Switch to eV/A and change sign
          Forces(3*i) += Fx*kcal2eV;
          Forces(3*i+1) += Fy*kcal2eV;
          Forces(3*i+2) += Fz*kcal2eV;
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
  call << " LICHM_" << Bead << ".xyz";
  call << " LICHM_" << Bead << ".key";
  call << " LICHM_" << Bead << ".grad";
  GlobalSys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKEREnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs TINKER MM energy calculations
  fstream ofile,ifile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0;
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    ofile << "a-axis " << LICHEMFormDouble(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormDouble(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormDouble(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }
  if (QMMM)
  {
    ofile << "polarizeterm none" << '\n'; //Remove polarization energy
  }
  ct = 0; //Generic counter
  if (QMMM)
  {
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
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion or Struct[i].PBregion or Struct[i].BAregion)
      {
        //New charges are only needed for QM atoms
        ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
        ofile << "0.0" << '\n';
      }
    }
  }
  if (AMOEBA or GEM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add multipoles
      if (Struct[i].QMregion or Struct[i].PBregion or Struct[i].BAregion)
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
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    ofile << LICHEMFormDouble(Lx,12) << " ";
    ofile << LICHEMFormDouble(Ly,12) << " ";
    ofile << LICHEMFormDouble(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Calculate MM potential energy
  call.str("");
  call << "analyze LICHM_";
  call << Bead << ".xyz E > LICHM_";
  call << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  call.str("");
  call << "LICHM_" << Bead << ".log";
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
    cerr << "Warning: No MM energy found after a calculation!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << Bead << ".xyz";
  call << " LICHM_" << Bead << ".log";
  call << " LICHM_" << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Calculate polarization energy
  if ((AMOEBA or GEM) and QMMM)
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
  fstream ofile,ifile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    ofile << "a-axis " << LICHEMFormDouble(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormDouble(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormDouble(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }
  ofile << "thermostat berendsen";
  ofile << '\n';
  ofile << "tau-temperature ";
  ofile << QMMMOpts.tautemp << '\n';
  ct = 0; //Generic counter
  if (QMMM or (Nfreeze > 0))
  {
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
  }
  if (CHRG)
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
      if (Struct[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        vector<int> Boundaries;
        Boundaries = TraceBoundary(Struct,i);
        double qnew = Struct[i].MP[Bead].q;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= Struct[Boundaries[j]].MP[Bead].q;
        }
        ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
        ofile << qnew;
        ofile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion)
      {
        //Write new multipole definition for the atom ID
        WriteTINKMpole(Struct,ofile,i,Bead);
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (Struct[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = Struct[i].MP[Bead].q; //Save a copy
        vector<int> Boundaries;
        Boundaries = TraceBoundary(Struct,i);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= Struct[Boundaries[j]].MP[Bead].q;
        }
        Struct[i].MP[Bead].q = qnew; //Save modified charge
        WriteTINKMpole(Struct,ofile,i,Bead);
        Struct[i].MP[Bead].q = qi; //Return to unmodified charge
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (Struct[i].BAregion)
      {
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
    }
  }
  ofile.flush();
  ofile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    ofile << LICHEMFormDouble(Lx,12) << " ";
    ofile << LICHEMFormDouble(Ly,12) << " ";
    ofile << LICHEMFormDouble(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Run optimization
  call.str("");
  call << "dynamic ";
  call << "LICHM_" << Bead << ".xyz ";
  call << QMMMOpts.Nsteps << " ";
  call << QMMMOpts.dt << " ";
  call << (QMMMOpts.Nsteps*QMMMOpts.dt/1000) << " ";
  call << "2 " << QMMMOpts.Temp;
  call << " > LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHM_" << Bead << ".001";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy); //Discard number of atoms
  if (PBCon)
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
  call << " LICHM_" << Bead << ".xyz";
  call << " LICHM_" << Bead << ".log";
  call << " LICHM_" << Bead << ".0*";
  call << " LICHM_" << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  return;
};

double TINKEROpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs TINKER MM optimization
  fstream ofile,ifile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0;
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_";
  call << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << Bead << ".key";
  ofile.open(call.str().c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "openmp-threads " << Ncpus << '\n';
  ofile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    ofile << "a-axis " << LICHEMFormDouble(Lx,12) << '\n';
    ofile << "b-axis " << LICHEMFormDouble(Ly,12) << '\n';
    ofile << "c-axis " << LICHEMFormDouble(Lz,12) << '\n';
    ofile << "alpha 90.0" << '\n';
    ofile << "beta 90.0" << '\n';
    ofile << "gamma 90.0" << '\n';
  }
  ct = 0; //Generic counter
  if (QMMM or (Nfreeze > 0))
  {
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
  }
  if (CHRG)
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
      if (Struct[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        vector<int> Boundaries;
        Boundaries = TraceBoundary(Struct,i);
        double qnew = Struct[i].MP[Bead].q;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= Struct[Boundaries[j]].MP[Bead].q;
        }
        ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
        ofile << qnew;
        ofile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (Struct[i].QMregion)
      {
        //Write new multipole definition for the atom ID
        WriteTINKMpole(Struct,ofile,i,Bead);
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (Struct[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = Struct[i].MP[Bead].q; //Save a copy
        vector<int> Boundaries;
        Boundaries = TraceBoundary(Struct,i);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= Struct[Boundaries[j]].MP[Bead].q;
        }
        Struct[i].MP[Bead].q = qnew; //Save modified charge
        WriteTINKMpole(Struct,ofile,i,Bead);
        Struct[i].MP[Bead].q = qi; //Return to unmodified charge
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
      if (Struct[i].BAregion)
      {
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
    }
  }
  ofile.flush();
  ofile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    ofile << LICHEMFormDouble(Lx,12) << " ";
    ofile << LICHEMFormDouble(Ly,12) << " ";
    ofile << LICHEMFormDouble(Lz,12) << " ";
    ofile << "90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormDouble(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Run optimization
  call.str("");
  call << "minimize LICHM_";
  call << Bead << ".xyz ";
  call << QMMMOpts.MMOptTol << " > LICHM_";
  call << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz_2";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy); //Discard number of atoms
  if (PBCon)
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
  call << " LICHM_" << Bead << ".xyz";
  call << " LICHM_" << Bead << ".log";
  call << " LICHM_" << Bead << ".xyz_*";
  call << " LICHM_" << Bead << ".key";
  GlobalSys = system(call.str().c_str());
  //Change units
  E *= kcal2eV;
  return E;
};

