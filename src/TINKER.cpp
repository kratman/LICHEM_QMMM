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
void FindTINKERClasses(vector<QMMMAtom>& QMMMData)
{
  //Parses TINKER parameter files to find atom classes
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  int ct; //Generic counter
  //Open generic key file
  inFile.open("tinker.key",ios_base::in);
  if (!inFile.good())
  {
    //Exit if files do not exist
    cout << "Error: Missing tinker.key file.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  //Find the parameter file
  bool FileFound = 0; //Bool to break loops
  while ((!inFile.eof()) and (!FileFound))
  {
    //Detect the name of the force field file
    inFile >> dummy;
    LICHEMLowerText(dummy);
    if (dummy == "parameters")
    {
      inFile >> dummy;
      FileFound = 1;
    }
  }
  //Open the parameters
  inFile.close();
  inFile.open(dummy.c_str(),ios_base::in);
  if (!FileFound)
  {
    //Exit if parameter file is not found
    cout << "Error: Cannot find TINKER parameter file.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  if (!inFile.good())
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
  while (!inFile.eof())
  {
    getline(inFile,dummy);
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
        if (QMMMData[i].numTyp == AtType)
        {
          QMMMData[i].numClass = AtClass;
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
void TINKERInduced(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                   int bead)
{
  //Function to extract induced dipoles
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  int ct; //Generic counter
  //Create TINKER xyz file
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Create new TINKER key file
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n'; //Make sure current line is empty
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  outFile << "save-induced" << '\n'; //Save induced dipoles
  outFile << "thermostat berendsen" << '\n';
  outFile << "tau-temperature 0.1" << '\n';
  ct = 0; //Generic counter
  if (QMMM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMregion or QMMMData[i].BAregion)
      {
        if (ct == 0)
        {
          //Start a new active line
          outFile << "active ";
        }
        else
        {
          //Place a space to separate values
          outFile << " ";
        }
        outFile << (QMMMData[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          outFile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if (QMMMData[i].QMregion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMPole(QMMMData,outFile,i,bead);
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].PBregion)
    {
      //Modify the charge to force charge balance with the boundaries
      double qi = QMMMData[i].MP[bead].q; //Save a copy
      vector<int> Boundaries;
      Boundaries = TraceBoundary(QMMMData,i);
      double qnew = qi;
      for (unsigned int j=0;j<Boundaries.size();j++)
      {
        //Subtract boundary atom charge
        qnew -= QMMMData[Boundaries[j]].MP[bead].q;
      }
      QMMMData[i].MP[bead].q = qnew; //Save modified charge
      WriteTINKMPole(QMMMData,outFile,i,bead);
      QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].BAregion)
    {
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
  }
  outFile.flush();
  outFile.close();
  //Calculate induced dipoles using dynamic
  call.str("");
  call << "dynamic LICHM_" << bead << ".xyz ";
  call << "1 1e-4 1e-7 2 0 > LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Extract induced dipoles from the MD cycle file
  call.str("");
  call << "LICHM_" << bead << ".001u";
  inFile.open(call.str().c_str(),ios_base::in);
  getline(inFile,dummy); //Clear number of atoms
  while (inFile.good())
  {
    int AtNum; //Identifies which atom was polarized
    //Parse file line by line
    getline(inFile,dummy);
    stringstream line(dummy);
    //Save dipoles for later
    line >> AtNum >> dummy; //Collect atom number and clear junk
    if (line.good())
    {
      AtNum -= 1; //Fixes array indexing
      line >> QMMMData[AtNum].MP[bead].IDx;
      line >> QMMMData[AtNum].MP[bead].IDy;
      line >> QMMMData[AtNum].MP[bead].IDz;
      //Change units from Debye to a.u.
      QMMMData[AtNum].MP[bead].IDx *= debye2au;
      QMMMData[AtNum].MP[bead].IDy *= debye2au;
      QMMMData[AtNum].MP[bead].IDz *= debye2au;
    }
  }
  inFile.close();
  //Delete junk files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".0*";
  call << " LICHM_" << bead << ".dyn";
  call << " LICHM_" << bead << ".log";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  return;
};

double TINKERPolEnergy(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                       int bead)
{
  //Function to extract the polarization energy
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double Epol = 0; //Polarization energy
  double Esolv = 0; //Solvation energy
  double E = 0; //Total energy for error checking
  int ct; //Generic counter
  //Create TINKER xyz file
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Create new TINKER key file
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  if (AMOEBA)
  {
    //Get rid of non-polarization interactions
    outFile << "polarizeterm only" << '\n';
  }
  if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model
    outFile << "solvateterm" << '\n';
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  ct = 0; //Generic counter
  if (QMMM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMregion or QMMMData[i].BAregion)
      {
        if (ct == 0)
        {
          //Start a new active line
          outFile << "active ";
        }
        else
        {
          //Place a space to separate values
          outFile << " ";
        }
        outFile << (QMMMData[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          outFile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if (QMMMData[i].QMregion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMPole(QMMMData,outFile,i,bead);
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].PBregion)
    {
      //Modify the charge to force charge balance with the boundaries
      double qi = QMMMData[i].MP[bead].q; //Save a copy
      vector<int> Boundaries;
      Boundaries = TraceBoundary(QMMMData,i);
      double qnew = qi;
      for (unsigned int j=0;j<Boundaries.size();j++)
      {
        //Subtract boundary atom charge
        qnew -= QMMMData[Boundaries[j]].MP[bead].q;
      }
      QMMMData[i].MP[bead].q = qnew; //Save modified charge
      WriteTINKMPole(QMMMData,outFile,i,bead);
      QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].BAregion)
    {
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
  }
  outFile.flush();
  outFile.close();
  //Calculate QMMM energy
  call.str("");
  call << "analyze LICHM_";
  call << bead << ".xyz E > LICHM_";
  call << bead << ".log";
  globalSys = system(call.str().c_str());
  //Extract polarization energy
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool Efound = 0;
  while (!inFile.eof())
  {
    inFile >> dummy;
    if (dummy == "Total")
    {
      inFile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        inFile >> dummy >> E;
        Efound = 1;
      }
    }
    if (dummy == "Polarization")
    {
      inFile >> Epol;
    }
    if (dummy == "Implicit")
    {
      inFile >> dummy;
      if (dummy == "Solvation")
      {
        inFile >> Esolv;
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
    cerr.flush(); //Print warning immediately
    Epol = 0; //Prevents errors when polarization is off
    Esolv = 0; //Prevents errors when implicit solvation is off
  }
  inFile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".log";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  //Return polarization and solvation energy in kcal/mol
  return Epol+Esolv;
};

double TINKERForces(vector<QMMMAtom>& QMMMData, VectorXd& Forces,
                    QMMMSettings& QMMMOpts, int bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream outFile,inFile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double Emm = 0.0;
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (QMMMData[i].QMregion or QMMMData[i].PBregion)
    {
      if (ct == 0)
      {
        //Start a new active line
        outFile << "active ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    outFile << '\n';
  }
  outFile << "group-inter" << '\n'; //Modify interactions
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add group 1 atoms
    if (QMMMData[i].QMregion or QMMMData[i].PBregion)
    {
      if (ct == 0)
      {
        //Start a new group line
        outFile << "group 1 ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate a group line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing group line
    outFile << '\n';
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion or QMMMData[i].PBregion or QMMMData[i].BAregion)
      {
        //New charges are needed for QM and PB atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << "0.0"; //Delete charges
        outFile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion or QMMMData[i].PBregion or QMMMData[i].BAregion)
      {
        double qi = 0;
        //remove charge
        qi = QMMMData[i].MP[bead].q;
        QMMMData[i].MP[bead].q = 0;
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q += qi; //Restore charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Run MM
  call.str("");
  call << "testgrad ";
  call << "LICHM_" << bead << ".xyz";
  call << " Y N N > ";
  call << "LICHM_" << bead << ".grad";
  globalSys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << bead << ".grad";
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
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".grad";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKERMMForces(vector<QMMMAtom>& QMMMData, VectorXd& Forces,
                      QMMMSettings& QMMMOpts, int bead)
{
  //Function to calculate the forces on MM atoms
  //NB: QM atoms are included in the array, but their forces are not updated
  fstream outFile,inFile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double Emm = 0.0; //Energy returned by TINKER
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model
    outFile << "solvateterm" << '\n';
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  ct = 0; //Generic counter
  if (QMMM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add inactive atoms
      if (QMMMData[i].QMregion or QMMMData[i].PBregion or QMMMData[i].frozen)
      {
        if (ct == 0)
        {
          //Start a new inactive line
          outFile << "inactive ";
        }
        else
        {
          //Place a space to separate values
          outFile << " ";
        }
        outFile << (QMMMData[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          outFile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion)
      {
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qnew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BAregion)
      {
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Run MM
  call.str("");
  call << "testgrad ";
  call << "LICHM_" << bead << ".xyz";
  call << " Y N N > ";
  call << "LICHM_" << bead << ".grad";
  globalSys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << bead << ".grad";
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
          if ((!QMMMData[i].QMregion) or (!QMMMData[i].PBregion) or
             (!QMMMData[i].frozen))
          {
            //Only add MM forces to the array
            Forces(3*i) += Fx*kcal2eV;
            Forces(3*i+1) += Fy*kcal2eV;
            Forces(3*i+2) += Fz*kcal2eV;
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
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".grad";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  //Return energy for error checking purposes
  return Emm;
};

double TINKERPolForces(vector<QMMMAtom>& QMMMData, VectorXd& Forces,
                       QMMMSettings& QMMMOpts, int bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream outFile,inFile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double Emm = 0.0;
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  if (AMOEBA)
  {
    //Get rid of non-polarization interactions
    outFile << "polarizeterm only" << '\n';
  }
  if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model
    outFile << "solvateterm" << '\n';
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (QMMMData[i].QMregion or QMMMData[i].PBregion)
    {
      if (ct == 0)
      {
        //Start a new active line
        outFile << "active ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    outFile << '\n';
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion)
      {
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qnew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BAregion)
      {
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Run MM
  call.str("");
  call << "testgrad ";
  call << "LICHM_" << bead << ".xyz";
  call << " Y N N > ";
  call << "LICHM_" << bead << ".grad";
  globalSys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << bead << ".grad";
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
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".grad";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKEREnergy(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
{
  //Runs TINKER MM energy calculations
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0;
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  if (QMMM)
  {
    outFile << "polarizeterm none" << '\n'; //Remove polarization energy
  }
  else if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model for pure MM calculations
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  ct = 0; //Generic counter
  if (QMMM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMregion or QMMMData[i].BAregion)
      {
        if (ct == 0)
        {
          //Start a new active line
          outFile << "active ";
        }
        else
        {
          //Place a space to separate values
          outFile << " ";
        }
        outFile << (QMMMData[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          outFile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion or QMMMData[i].PBregion or QMMMData[i].BAregion)
      {
        //New charges are only needed for QM atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << "0.0" << '\n';
      }
    }
  }
  if (AMOEBA or GEM)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add multipoles
      if (QMMMData[i].QMregion or QMMMData[i].PBregion or QMMMData[i].BAregion)
      {
        double qi = 0;
        //remove charge
        qi = QMMMData[i].MP[bead].q;
        QMMMData[i].MP[bead].q = 0;
        //Write new multipole definition for the atom ID
        WriteTINKMPole(QMMMData,outFile,i,bead);
        //Restore charge
        QMMMData[i].MP[bead].q = qi;
        //Remove polarization
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Calculate MM potential energy
  call.str("");
  call << "analyze LICHM_";
  call << bead << ".xyz E > LICHM_";
  call << bead << ".log";
  globalSys = system(call.str().c_str());
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  //Read MM potential energy
  bool Efound = 0;
  while (!inFile.eof())
  {
    inFile >> dummy;
    if (dummy == "Total")
    {
      inFile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        inFile >> dummy >> E;
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
    cerr.flush(); //Print warning immediately
    E = hugeNum; //Large number to reject step
  }
  inFile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".log";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  //Calculate polarization energy
  if ((AMOEBA or GEM or QMMMOpts.useImpSolv) and QMMM)
  {
    //Correct polarization energy for QMMM simulations
    E += TINKERPolEnergy(QMMMData,QMMMOpts,bead);
  }
  //Change units
  E *= kcal2eV;
  return E;
};

void TINKERDynamics(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
{
  //Runs TINKER MD (for MM atoms only)
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  outFile << "thermostat berendsen";
  outFile << '\n';
  outFile << "tau-temperature ";
  outFile << (QMMMOpts.tauTemp/1000); //Converted from fs to ps
  outFile << '\n';
  ct = 0; //Generic counter
  if (QMMM or (Nfreeze > 0))
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMregion or QMMMData[i].BAregion)
      {
        if (!QMMMData[i].frozen)
        {
          if (ct == 0)
          {
            //Start a new active line
            outFile << "active ";
          }
          else
          {
            //Place a space to separate values
            outFile << " ";
          }
          outFile << (QMMMData[i].id+1);
          ct += 1;
          if (ct == 10)
          {
            //terminate an active line
            ct = 0;
            outFile << '\n';
          }
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion)
      {
        //New charges are only needed for QM atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << QMMMData[i].MP[bead].q;
        outFile << '\n';
      }
      if (QMMMData[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i);
        double qnew = QMMMData[i].MP[bead].q;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
        }
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << qnew;
        outFile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion)
      {
        //Write new multipole definition for the atom ID
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qnew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BAregion)
      {
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Run optimization
  call.str("");
  call << "dynamic ";
  call << "LICHM_" << bead << ".xyz ";
  call << QMMMOpts.NSteps << " ";
  call << QMMMOpts.dt << " ";
  call << (QMMMOpts.NSteps*QMMMOpts.dt/1000) << " ";
  call << "2 " << QMMMOpts.temp;
  call << " > LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHM_" << bead << ".001";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy); //Discard number of atoms
    if (PBCon)
    {
      //Discard PBC information
      getline(inFile,dummy);
    }
    for (int i=0;i<Natoms;i++)
    {
      getline(inFile,dummy);
      stringstream line(dummy);
      //Read new positions
      line >> dummy >> dummy; //Discard atom ID and type
      line >> QMMMData[i].P[bead].x;
      line >> QMMMData[i].P[bead].y;
      line >> QMMMData[i].P[bead].z;
    }
  }
  else
  {
    //Print error message
    cerr << "Warning: No structure found after the dynamics!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    cerr.flush();
    //Remove restart file
    call.str("");
    call << "rm -f";
    call << " LICHM_" << bead << ".dyn";
    globalSys = system(call.str().c_str());
  }
  inFile.close();
  //Clean up all files except the .dyn files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".log";
  call << " LICHM_" << bead << ".0*";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  return;
};

MatrixXd TINKERHessian(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                       int bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream outFile,inFile,MMlog; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  int Ndof = 3*(Nqm+Npseudo); //Number of degrees of freedom
  MatrixXd MMHess(Ndof,Ndof);
  MMHess.setZero();
  int ct; //Generic counter
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if (QMMMData[i].QMregion or QMMMData[i].PBregion)
    {
      if (ct == 0)
      {
        //Start a new active line
        outFile << "active ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    outFile << '\n';
  }
  outFile << "group-inter" << '\n'; //Modify interactions
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add group 1 atoms
    if (QMMMData[i].QMregion or QMMMData[i].PBregion)
    {
      if (ct == 0)
      {
        //Start a new group line
        outFile << "group 1 ";
      }
      else
      {
        //Place a space to separate values
        outFile << " ";
      }
      outFile << (QMMMData[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate a group line
        ct = 0;
        outFile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing group line
    outFile << '\n';
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion or QMMMData[i].PBregion or QMMMData[i].BAregion)
      {
        //New charges are needed for QM and PB atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << "0.0"; //Delete charges
        outFile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion or QMMMData[i].PBregion or QMMMData[i].BAregion)
      {
        double qi = 0;
        //remove charge
        qi = QMMMData[i].MP[bead].q;
        QMMMData[i].MP[bead].q = 0;
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q += qi; //Restore charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Run MM
  call.str("");
  call << "testhess ";
  call << "LICHM_" << bead << ".xyz";
  call << " Y N > ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Collect MM forces
  call.str("");
  call << "LICHM_" << bead << ".hes";
  MMlog.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool HessDone = 0;
  if (MMlog.good() and CheckFile(call.str()))
  {
    HessDone = 1;
    //Clear junk
    getline(MMlog,dummy);
    getline(MMlog,dummy);
    getline(MMlog,dummy);
    //Read diagonal elements
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMregion or QMMMData[i].PBregion)
      {
        //Read QM and PB diagonal elements
        MMlog >> MMHess(ct,ct);
        ct += 1;
        MMlog >> MMHess(ct,ct);
        ct += 1;
        MMlog >> MMHess(ct,ct);
        ct += 1;
      }
      else
      {
        //Read zeros
        MMlog >> dummy;
        MMlog >> dummy;
        MMlog >> dummy;
      }
    }
    //Read off-diagonal elements
    for (int i=0;i<Ndof;i++)
    {
      //Clear junk
      getline(MMlog,dummy);
      getline(MMlog,dummy);
      getline(MMlog,dummy);
      //Read elements
      for (int j=(i+1);j<Ndof;j++)
      {
        //Read value
        MMlog >> MMHess(i,j);
        //Apply symmetry
        MMHess(j,i) = MMHess(i,j);
      }
    }
  }
  //Change units
  MMHess *= (kcal2eV*bohrRad*bohrRad/har2eV); //Switch to a.u.
  //Check for errors
  if (!HessDone)
  {
    //Calculation did not finish
    cerr << "Error: No force constants recovered!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
  }
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".hes";
  call << " LICHM_" << bead << ".log";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  //Return
  return MMHess;
};

double TINKEROpt(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
{
  //Runs TINKER MM optimization
  fstream outFile,inFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy; //Generic string
  double E = 0;
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  call.str("");
  call << "cp tinker.key LICHM_";
  call << bead << ".key";
  globalSys = system(call.str().c_str());
  //Update key file
  call.str("");
  call << "LICHM_";
  call << bead << ".key";
  outFile.open(call.str().c_str(),ios_base::app|ios_base::out);
  outFile << '\n';
  if (QMMM)
  {
    outFile << "#LICHEM QMMM keywords"; //Marks the changes
  }
  else
  {
    outFile << "#LICHEM MM keywords"; //Marks the changes
  }
  outFile << '\n';
  if (QMMMOpts.useMMCut)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald and truncate vdW forces
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.MMOptCut,12);
      outFile << '\n';
      outFile << "ewald" << '\n';
    }
    else
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.MMOptCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.MMOptCut,12);
      outFile << '\n';
    }
  }
  else if (QMMMOpts.useLREC)
  {
    //Apply cutoff
    if (QMMMOpts.useEwald and PBCon)
    {
      //Use Ewald or PME
      outFile << "ewald" << '\n';
    }
    else if (!QMMMOpts.useImpSolv)
    {
      //Use smoothing functions
      outFile << "cutoff " << LICHEMFormFloat(QMMMOpts.LRECCut,12);
      outFile << '\n';
      outFile << "taper " << LICHEMFormFloat(0.90*QMMMOpts.LRECCut,12);
      outFile << '\n';
    }
  }
  if (QMMMOpts.useImpSolv)
  {
    //Add the implicit solvation model
    outFile << "solvate " << QMMMOpts.solvModel;
    outFile << '\n';
  }
  outFile << "openmp-threads " << Ncpus << '\n';
  outFile << "digits 12" << '\n'; //Increase precision
  if (PBCon)
  {
    //PBC defined twice for safety
    outFile << "a-axis " << LICHEMFormFloat(Lx,12) << '\n';
    outFile << "b-axis " << LICHEMFormFloat(Ly,12) << '\n';
    outFile << "c-axis " << LICHEMFormFloat(Lz,12) << '\n';
    outFile << "alpha 90.0" << '\n';
    outFile << "beta 90.0" << '\n';
    outFile << "gamma 90.0" << '\n';
  }
  ct = 0; //Generic counter
  if (QMMM or (Nfreeze > 0))
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if (QMMMData[i].MMregion or QMMMData[i].BAregion)
      {
        if (!QMMMData[i].frozen)
        {
          if (ct == 0)
          {
            //Start a new active line
            outFile << "active ";
          }
          else
          {
            //Place a space to separate values
            outFile << " ";
          }
          outFile << (QMMMData[i].id+1);
          ct += 1;
          if (ct == 10)
          {
            //terminate an active line
            ct = 0;
            outFile << '\n';
          }
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      outFile << '\n';
    }
  }
  if (CHRG)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion)
      {
        //New charges are only needed for QM atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << QMMMData[i].MP[bead].q;
        outFile << '\n';
      }
      if (QMMMData[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i);
        double qnew = QMMMData[i].MP[bead].q;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
        }
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << qnew;
        outFile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMregion)
      {
        //Write new multipole definition for the atom ID
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBregion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> Boundaries;
        Boundaries = TraceBoundary(QMMMData,i);
        double qnew = qi;
        for (unsigned int j=0;j<Boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qnew -= QMMMData[Boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qnew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BAregion)
      {
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
    }
  }
  outFile.flush();
  outFile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  if (PBCon)
  {
    //Write box size
    outFile << LICHEMFormFloat(Lx,12) << " ";
    outFile << LICHEMFormFloat(Ly,12) << " ";
    outFile << LICHEMFormFloat(Lz,12) << " ";
    outFile << "90.0 90.0 90.0";
    outFile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Run optimization
  call.str("");
  call << "minimize LICHM_";
  call << bead << ".xyz ";
  call << QMMMOpts.MMOptTol << " > LICHM_";
  call << bead << ".log";
  globalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHM_" << bead << ".xyz_2";
  inFile.open(call.str().c_str(),ios_base::in);
  getline(inFile,dummy); //Discard number of atoms
  if (PBCon)
  {
    //Discard PBC information
    getline(inFile,dummy);
  }
  for (int i=0;i<Natoms;i++)
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    //Read new positions
    line >> dummy >> dummy; //Discard atom ID and type
    line >> QMMMData[i].P[bead].x;
    line >> QMMMData[i].P[bead].y;
    line >> QMMMData[i].P[bead].z;
  }
  inFile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " LICHM_" << bead << ".xyz";
  call << " LICHM_" << bead << ".log";
  call << " LICHM_" << bead << ".xyz_*";
  call << " LICHM_" << bead << ".key";
  call << " LICHM_" << bead << ".err";
  globalSys = system(call.str().c_str());
  //Change units
  E *= kcal2eV;
  return E;
};

