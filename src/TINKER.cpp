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
  bool fileFound = 0; //Bool to break loops
  while ((!inFile.eof()) and (!fileFound))
  {
    //Detect the name of the force field file
    inFile >> dummy;
    LICHEMLowerText(dummy);
    if (dummy == "parameters")
    {
      inFile >> dummy;
      fileFound = 1;
    }
  }
  //Open the parameters
  inFile.close();
  inFile.open(dummy.c_str(),ios_base::in);
  if (!fileFound)
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
  while ((!inFile.eof()) and inFile.good())
  {
    getline(inFile,dummy);
    stringstream fullLine(dummy);
    fullLine >> dummy;
    LICHEMLowerText(dummy);
    if (dummy == "atom")
    {
      int atType,atClass;
      fullLine >> atType;
      fullLine >> atClass;
      for (int i=0;i<Natoms;i++)
      {
        if (QMMMData[i].numTyp == atType)
        {
          QMMMData[i].numClass = atClass;
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
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
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
    if (QMMMData[i].QMRegion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMPole(QMMMData,outFile,i,bead);
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].PBRegion)
    {
      //Modify the charge to force charge balance with the boundaries
      double qi = QMMMData[i].MP[bead].q; //Save a copy
      vector<int> boundaries;
      boundaries = TraceBoundary(QMMMData,i);
      double qNew = qi;
      for (unsigned int j=0;j<boundaries.size();j++)
      {
        //Subtract boundary atom charge
        qNew -= QMMMData[boundaries[j]].MP[bead].q;
      }
      QMMMData[i].MP[bead].q = qNew; //Save modified charge
      WriteTINKMPole(QMMMData,outFile,i,bead);
      QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].BARegion)
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
    int atNum; //Identifies which atom was polarized
    //Parse file line by line
    getline(inFile,dummy);
    stringstream line(dummy);
    //Save dipoles for later
    line >> atNum >> dummy; //Collect atom number and clear junk
    if (line.good())
    {
      atNum -= 1; //Fixes array indexing
      line >> QMMMData[atNum].MP[bead].IDx;
      line >> QMMMData[atNum].MP[bead].IDy;
      line >> QMMMData[atNum].MP[bead].IDz;
      //Change units from Debye to a.u.
      QMMMData[atNum].MP[bead].IDx *= debye2au;
      QMMMData[atNum].MP[bead].IDy *= debye2au;
      QMMMData[atNum].MP[bead].IDz *= debye2au;
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
  double EPol = 0; //Polarization energy
  double ESolv = 0; //Solvation energy
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
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
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
    if (QMMMData[i].QMRegion)
    {
      //Write new multipole definition for the atom ID
      WriteTINKMPole(QMMMData,outFile,i,bead);
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].PBRegion)
    {
      //Modify the charge to force charge balance with the boundaries
      double qi = QMMMData[i].MP[bead].q; //Save a copy
      vector<int> boundaries;
      boundaries = TraceBoundary(QMMMData,i);
      double qNew = qi;
      for (unsigned int j=0;j<boundaries.size();j++)
      {
        //Subtract boundary atom charge
        qNew -= QMMMData[boundaries[j]].MP[bead].q;
      }
      QMMMData[i].MP[bead].q = qNew; //Save modified charge
      WriteTINKMPole(QMMMData,outFile,i,bead);
      QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
      outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
      outFile << '\n';
    }
    if (QMMMData[i].BARegion)
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
  bool EFound = 0;
  while ((!inFile.eof()) and inFile.good())
  {
    inFile >> dummy;
    if (dummy == "Total")
    {
      inFile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        inFile >> dummy >> E;
        EFound = 1;
      }
    }
    if (dummy == "Polarization")
    {
      inFile >> EPol;
    }
    if (dummy == "Implicit")
    {
      inFile >> dummy;
      if (dummy == "Solvation")
      {
        inFile >> ESolv;
      }
    }
  }
  if (!EFound)
  {
    //Warn user if no energy was found
    cerr << "Warning: No MM energy found after a calculation!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    EPol = 0; //Prevents errors when polarization is off
    ESolv = 0; //Prevents errors when implicit solvation is off
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
  return EPol+ESolv;
};

double TINKERForces(vector<QMMMAtom>& QMMMData, VectorXd& forces,
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
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
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
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
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
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
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
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
      {
        double qi = 0;
        //Remove charge
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
  fstream MMGrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << bead << ".grad";
  MMGrad.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool gradDone = 0;
  while ((!MMGrad.eof()) and MMGrad.good() and (!gradDone))
  {
    getline(MMGrad,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Type")
    {
      line >> dummy >> dummy;
      if (dummy == "dE/dX")
      {
        gradDone = 1; //Not grad school, that lasts forever
        getline(MMGrad,dummy);
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          //Convoluted, but "easy"
          getline(MMGrad,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> fX;
          line >> fY;
          line >> fZ;
          //Change from gradient to force
          fX *= -1;
          fY *= -1;
          fZ *= -1;
          //Switch to eV/A and save forces
          forces(3*i) += fX*kcal2eV;
          forces(3*i+1) += fY*kcal2eV;
          forces(3*i+2) += fZ*kcal2eV;
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
  MMGrad.close();
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

double TINKERMMForces(vector<QMMMAtom>& QMMMData, VectorXd& forces,
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
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].frozen)
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
      if (QMMMData[i].QMRegion)
      {
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> boundaries;
        boundaries = TraceBoundary(QMMMData,i);
        double qNew = qi;
        for (unsigned int j=0;j<boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qNew -= QMMMData[boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qNew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BARegion)
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
  fstream MMGrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << bead << ".grad";
  MMGrad.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool gradDone = 0;
  while ((!MMGrad.eof()) and MMGrad.good() and (!gradDone))
  {
    getline(MMGrad,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Type")
    {
      line >> dummy >> dummy;
      if (dummy == "dE/dX")
      {
        gradDone = 1; //Not grad school, that lasts forever
        getline(MMGrad,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if ((QMMMData[i].MMRegion or QMMMData[i].BARegion) and
             (!QMMMData[i].frozen))
          {
            //Only update MM and BA forces in the array
            double fX = 0;
            double fY = 0;
            double fZ = 0;
            //Convoluted, but "easy"
            getline(MMGrad,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy; //Clear junk
            line >> fX;
            line >> fY;
            line >> fZ;
            //Change from gradient to force
            fX *= -1;
            fY *= -1;
            fZ *= -1;
            //Switch to eV/A and save the forces
            forces(3*i) = fX*kcal2eV;
            forces(3*i+1) = fY*kcal2eV;
            forces(3*i+2) = fZ*kcal2eV;
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
  MMGrad.close();
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

double TINKERPolForces(vector<QMMMAtom>& QMMMData, VectorXd& forces,
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
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
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
      if (QMMMData[i].QMRegion)
      {
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> boundaries;
        boundaries = TraceBoundary(QMMMData,i);
        double qNew = qi;
        for (unsigned int j=0;j<boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qNew -= QMMMData[boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qNew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BARegion)
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
  fstream MMGrad; //QMMM output
  //Open files
  call.str("");
  call << "LICHM_" << bead << ".grad";
  MMGrad.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool gradDone = 0;
  while ((!MMGrad.eof()) and MMGrad.good() and (!gradDone))
  {
    getline(MMGrad,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Type")
    {
      line >> dummy >> dummy;
      if (dummy == "dE/dX")
      {
        gradDone = 1; //Not grad school, that lasts forever
        getline(MMGrad,dummy);
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          //Convoluted, but "easy"
          getline(MMGrad,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> fX;
          line >> fY;
          line >> fZ;
          //Change from gradient to force
          fX *= -1;
          fY *= -1;
          fZ *= -1;
          //Switch to eV/A and change sign
          forces(3*i) += fX*kcal2eV;
          forces(3*i+1) += fY*kcal2eV;
          forces(3*i+2) += fZ*kcal2eV;
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
  MMGrad.close();
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

double TINKEREnergy(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
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
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
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
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
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
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
      {
        double qi = 0;
        //Remove charge
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
  bool EFound = 0;
  while ((!inFile.eof()) and inFile.good())
  {
    inFile >> dummy;
    if (dummy == "Total")
    {
      inFile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        inFile >> dummy >> E;
        EFound = 1;
      }
    }
  }
  if (!EFound)
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

MatrixXd TINKERHessian(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                       int bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream outFile,inFile,MMLog; //Generic file streams
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
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
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
    if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
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
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
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
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion or QMMMData[i].BARegion)
      {
        double qi = 0;
        //Remove charge
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
  MMLog.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool hessDone = 0;
  if (MMLog.good() and CheckFile(call.str()))
  {
    hessDone = 1;
    //Clear junk
    getline(MMLog,dummy);
    getline(MMLog,dummy);
    getline(MMLog,dummy);
    //Read diagonal elements
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
      {
        //Read QM and PB diagonal elements
        MMLog >> MMHess(ct,ct);
        ct += 1;
        MMLog >> MMHess(ct,ct);
        ct += 1;
        MMLog >> MMHess(ct,ct);
        ct += 1;
      }
      else
      {
        //Read zeros
        MMLog >> dummy;
        MMLog >> dummy;
        MMLog >> dummy;
      }
    }
    //Read off-diagonal elements
    for (int i=0;i<Ndof;i++)
    {
      //Clear junk
      getline(MMLog,dummy);
      getline(MMLog,dummy);
      getline(MMLog,dummy);
      //Read elements
      for (int j=(i+1);j<Ndof;j++)
      {
        //Read value
        MMLog >> MMHess(i,j);
        //Apply symmetry
        MMHess(j,i) = MMHess(i,j);
      }
    }
  }
  //Change units
  MMHess *= (kcal2eV*bohrRad*bohrRad/har2eV); //Switch to a.u.
  //Check for errors
  if (!hessDone)
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
      if (QMMMData[i].MMRegion or QMMMData[i].BARegion)
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
      if (QMMMData[i].QMRegion)
      {
        //New charges are only needed for QM atoms
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << QMMMData[i].MP[bead].q;
        outFile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        vector<int> boundaries;
        boundaries = TraceBoundary(QMMMData,i);
        double qNew = QMMMData[i].MP[bead].q;
        for (unsigned int j=0;j<boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qNew -= QMMMData[boundaries[j]].MP[bead].q;
        }
        outFile << "charge " << (-1*(QMMMData[i].id+1)) << " ";
        outFile << qNew;
        outFile << '\n';
      }
    }
  }
  if (AMOEBA)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if (QMMMData[i].QMRegion)
      {
        //Write new multipole definition for the atom ID
        WriteTINKMPole(QMMMData,outFile,i,bead);
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].PBRegion)
      {
        //Modify the charge to force charge balance with the boundaries
        double qi = QMMMData[i].MP[bead].q; //Save a copy
        vector<int> boundaries;
        boundaries = TraceBoundary(QMMMData,i);
        double qNew = qi;
        for (unsigned int j=0;j<boundaries.size();j++)
        {
          //Subtract boundary atom charge
          qNew -= QMMMData[boundaries[j]].MP[bead].q;
        }
        QMMMData[i].MP[bead].q = qNew; //Save modified charge
        WriteTINKMPole(QMMMData,outFile,i,bead);
        QMMMData[i].MP[bead].q = qi; //Return to unmodified charge
        outFile << "polarize -" << (QMMMData[i].id+1) << " 0.0 0.0";
        outFile << '\n';
      }
      if (QMMMData[i].BARegion)
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

