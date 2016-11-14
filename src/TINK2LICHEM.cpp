/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Functions to convert TINKER files to LICHEM format

 Reference for TINKER:
 Ponder, TINKER - Software Tools for Molecular Design

*/

/*!
  \ingroup TINKER
*/
///@{

//Parsers
void TINK2LICHEM(int& argc, char**& argv)
{
  //Local variables
  fstream tinkXYZ,tinkKey,paramFile; //Input
  fstream posFile,conFile,regFile; //Output
  stringstream line;
  string dummy; //Generic string
  int Ninact = 0;
  bool someFroz = 0;
  bool someAct = 0;
  int ct; //Generic counter
  //Open files
  posFile.open("xyzfile.xyz",ios_base::out);
  conFile.open("connect.inp",ios_base::out);
  regFile.open("regions.inp",ios_base::out);
  //Read arguments
  PBCon = 0;
  cout << '\n';
  cout << "Reading files: ";
  for (int i=0;i<argc;i++)
  {
    dummy = string(argv[i]);
    if (dummy == "-t")
    {
      //TINKER XYZ file
      tinkXYZ.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
      cout << " ";
    }
    if (dummy == "-k")
    {
      //TINKER key file
      tinkKey.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
      cout << " ";
    }
    if (dummy == "-p")
    {
      //Flag for PBC
      dummy = string(argv[i+1]);
      LICHEMLowerText(dummy);
      if ((dummy == "yes") or (dummy == "true"))
      {
        PBCon = 1;
      }
    }
  }
  //Start writing the files
  if (PBCon)
  {
    //Remind forgetful users that they said the system was periodic
    cout << '\n';
    cout << " with PBC";
  }
  cout << "...";
  cout << '\n';
  getline(tinkXYZ,dummy);
  line.str(dummy);
  line >> Natoms;
  line >> Nqm;
  if (!line.eof())
  {
    //Read QMMM information
    //NB: This only works if the atoms are in the correct order:
    // QM->PB->BA->MM
    line >> Npseudo;
    line >> Nbound;
  }
  else
  {
    //Not needed for a standard MM XYZ file
    Nqm = 0;
    Npseudo = 0;
    Nbound = 0;
  }
  Nmm = Natoms-Nqm-Npseudo-Nbound;
  if (Natoms != Nmm)
  {
    //QMMM input
    cout << '\n';
    cout << "Auto-detected file type: ";
    cout << "TINKER QMMM xyz";
    cout << '\n';
  }
  if (Natoms == Nmm)
  {
    //MM input
    cout << '\n';
    cout << "Auto-detected file type: ";
    cout << "TINKER xyz";
    cout << '\n';
  }
  //Print settings
  cout << "  Total atoms: " << Natoms << '\n';
  cout << "  QM atoms: " << Nqm << '\n';
  cout << "  Pseudo-atoms: " << Npseudo << '\n';
  cout << "  Boundary atoms: " << Nbound << '\n';
  cout << "  MM atoms: " << Nmm << '\n';
  cout << '\n';
  cout << "Creating xyzfile.xyz, connect.inp, and regions.inp";
  cout << "...";
  cout << '\n' << '\n';
  //Create lists for grouping frozen atoms
  vector<bool> actives;
  vector<bool> froz;
  for (int i=0;i<Natoms;i++)
  {
    //Seems contradictory, however, this lets the script
    //check for both active and inactive commands
    actives.push_back(0); //Define all as inactive
    froz.push_back(0); //Define all as active
  }
  //Find parameter file, active atoms, and inactive atoms
  if (!tinkKey.good())
  {
    cout << '\n';
    cout << "Error: Could not read TINKER key file!";
    cout << '\n' << '\n';
    cout.flush();
    exit(0);
  }
  while (!tinkKey.eof())
  {
    //Read key file line by line
    getline(tinkKey,dummy);
    stringstream line(dummy);
    //Read string item by item
    line >> dummy;
    LICHEMLowerText(dummy);
    if (dummy == "parameters")
    {
      line >> dummy;
      paramFile.open(dummy.c_str(),ios_base::in);
    }
    if (dummy == "active")
    {
      someAct = 1;
      int atNum;
      while (line >> atNum)
      {
        actives[atNum-1] = 1;
      }
    }
    if (dummy == "inactive")
    {
      someFroz = 1;
      int atNum;
      while (line >> atNum)
      {
        froz[atNum-1] = 1;
      }
    }
  }
  if (!paramFile.good())
  {
    cout << '\n';
    cout << "Error: Could not open TINKER parameter file!";
    cout << '\n' << '\n';
    cout.flush();
    exit(0);
  }
  if (someAct or someFroz)
  {
    for (int i=0;i<Natoms;i++)
    {
      if ((!actives[i]) and (!froz[i]))
      {
        if (someFroz)
        {
          froz[i] = 0;
        }
        if (someAct)
        {
          froz[i] = 1;
        }
      }
      if (actives[i] and (!froz[i]))
      {
        froz[i] = 0;
      }
      if ((!actives[i]) and froz[i])
      {
        froz[i] = 1;
      }
      if (actives[i] and froz[i])
      {
        froz[i] = 0;
        cout << '\n';
        cout << "Warning: Atom " << i;
        cout << " is listed as both active and inactive.";
        cout << '\n';
        cout << " The atom will be assume to be active.";
        cout << '\n';
      }
    }
    for (int i=0;i<Natoms;i++)
    {
      //Count frozen atoms
      if (i >= (Nqm+Npseudo+Nbound))
      {
        if (froz[i])
        {
          Ninact += 1;
        }
      }
      else
      {
        //QM, PB, and BA cannot be frozen
        froz[i] = 0;
      }
    }
  }
  if (Ninact > 0)
  {
    cout << '\n';
    cout << "Detected " << Ninact;
    cout << " frozen atoms in the key file...";
    cout << '\n';
  }
  //Collect PBC information
  if (PBCon)
  {
    //Grab whole line
    getline(tinkXYZ,dummy);
    stringstream line(dummy);
    //Read box lengths
    line >> Lx >> Ly >> Lz;
  }
  //Write region file
  regFile << fixed;
  if (PBCon)
  {
    //Print PBC settings
    regFile << "PBC: Yes" << '\n';
    regFile << "Box_size: ";
    regFile << Lx << " ";
    regFile << Ly << " ";
    regFile << Lz << '\n';
  }
  regFile << "QM_atoms: ";
  regFile << Nqm << '\n';
  ct = 0;
  for (int i=0;i<Nqm;i++)
  {
    regFile << i;
    ct += 1;
    if ((ct == 10) or (i == Nqm-1))
    {
      regFile << '\n';
      ct = 0;
    }
    else
    {
      regFile << " ";
    }
  }
  regFile << "Pseudobond_atoms: ";
  regFile << Npseudo << '\n';
  ct = 0;
  for (int i=0;i<Npseudo;i++)
  {
    regFile << i+Nqm;
    ct += 1;
    if ((ct == 10) or (i == Npseudo-1))
    {
      regFile << '\n';
      ct = 0;
    }
    else
    {
      regFile << " ";
    }
  }
  regFile << "Boundary_atoms: ";
  regFile << Nbound << '\n';
  ct = 0;
  for (int i=0;i<Nbound;i++)
  {
    regFile << i+Nqm+Npseudo;
    ct += 1;
    if ((ct == 10) or (i == Nbound-1))
    {
      regFile << '\n';
      ct = 0;
    }
    else
    {
      regFile << " ";
    }
  }
  regFile << "Frozen_atoms: ";
  regFile << Ninact << '\n';
  ct = 0;
  for (int i=0;i<Natoms;i++)
  {
    if (froz[i])
    {
      regFile << i;
      ct += 1;
      if (ct == 10)
      {
        regFile << '\n';
        ct = 0;
      }
      else
      {
        regFile << " ";
      }
    }
  }
  regFile << '\n';
  //Read force field parameters
  vector<vector<string> > masses;
  vector<vector<string> > nucCharges;
  vector<vector<string> > charges;
  while (!paramFile.eof())
  {
    //Read line by line for atom type info
    getline(paramFile,dummy);
    stringstream line(dummy);
    vector<string> fullLine;
    while (line >> dummy)
    {
      fullLine.push_back(dummy);
    }
    if (fullLine.size() > 2)
    {
      if (fullLine[0] == "atom")
      {
        vector<string> tmp;
        tmp.push_back(fullLine[1]);
        tmp.push_back(fullLine[fullLine.size()-2]);
        masses.push_back(tmp);
      }
      if (fullLine[0] == "atom")
      {
        vector<string> tmp;
        tmp.push_back(fullLine[1]);
        tmp.push_back(fullLine[fullLine.size()-3]);
        nucCharges.push_back(tmp);
      }
      if (fullLine[0] == "charge")
      {
        vector<string> tmp;
        tmp.push_back(fullLine[1]);
        tmp.push_back(fullLine[2]);
        charges.push_back(tmp);
      }
    }
  }
  //Create connectivity and xyx files
  posFile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    getline(tinkXYZ,dummy);
    stringstream line(dummy);
    vector<string> fullLine;
    while (line >> dummy)
    {
      fullLine.push_back(dummy);
    }
    string numTyp = fullLine[5];
    bool ZFound = 0;
    for (unsigned int j=0;j<nucCharges.size();j++)
    {
      if (nucCharges[j][0] == numTyp)
      {
        //Find QM atom type
        stringstream line(nucCharges[j][1]);
        int nuChg;
        line >> nuChg;
        posFile << chemTable.typing(nuChg);
        posFile << " ";
        ZFound = 1;
      }
    }
    if (!ZFound)
    {
      //Print error
      cout << "Warning: Missing nuclear charge!";
      cout << " The .prm file may be incomplete.";
      cout << '\n';
      cout << "LICHEM cannot continue...";
      cout << '\n';
      //Dump data to the files and exit
      posFile.flush();
      conFile.flush();
      regFile.flush();
      cout.flush();
      exit(0);
    }
    //Write positions
    posFile << fullLine[2] << " ";
    posFile << fullLine[3] << " ";
    posFile << fullLine[4] << '\n';
    //Write connectivity
    conFile << i << " ";
    conFile << fullLine[1] << " ";
    conFile << numTyp << " ";
    bool massFound = 0;
    for (unsigned int j=0;j<masses.size();j++)
    {
      if (masses[j][0] == numTyp)
      {
        conFile << masses[j][1] << " ";
        massFound = 1;
      }
    }
    if (!massFound)
    {
      //Print error
      cout << "Error: Missing mass! The .prm file may be incomplete.";
      cout << '\n';
      cout << "LICHEM cannot continue...";
      cout << '\n';
      //Dump output and quit
      posFile.flush();
      conFile.flush();
      regFile.flush();
      cout.flush();
      exit(0);
    }
    bool chargeFound = 0;
    for (unsigned int j=0;j<charges.size();j++)
    {
      if (charges[j][0] == numTyp)
      {
        conFile << charges[j][1] << " ";
        chargeFound = 1;
      }
    }
    if (!chargeFound)
    {
      //Print error
      cout << "Warning: Missing charge! The .prm file may be incomplete.";
      cout << '\n';
      conFile << "0.00" << " "; //Set charge to zero
    }
    int nBonds = fullLine.size()-6;
    conFile << nBonds;
    for (int j=0;j<nBonds;j++)
    {
      conFile << " ";
      int bondID;
      stringstream line(fullLine[j+6]);
      line >> bondID;
      conFile << (bondID-1);
    }
    conFile << '\n';
  }
  //Quit LICHEM
  cout << '\n';
  cout << "Do not forget to manually add the simulation input to the";
  cout << " regions.inp file.";
  cout << '\n' << '\n';
  cout << "Conversion complete.";
  cout << '\n';
  cout << '\n';
  cout.flush();
  tinkXYZ.close();
  tinkKey.close();
  paramFile.close();
  posFile.close();
  conFile.close();
  regFile.close();
  exit(0); //Quit
  return;
};

void LICHEM2TINK(int& argc, char**& argv)
{
  //Creates TINKER input from LICHEM input
  fstream posFile,conFile,outFile; //File streams
  string dummy; //Generic string
  //Read arguments
  bool doQuit = 0; //Exit with an error
  cout << "Reading LICHEM input: ";
  for (int i=0;i<argc;i++)
  {
    dummy = string(argv[i]);
    if (dummy == "-x")
    {
      //Open XYZ file
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open XYZ file!!!";
        cout << '\n';
        doQuit = 1;
      }
      xyzFilename = file.str();
      posFile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
      cout << " ";
    }
    if (dummy == "-c")
    {
      //Open connectivity file
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open connectivity file!!!";
        cout << '\n';
        doQuit = 1;
      }
      conFilename = file.str();
      conFile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
      cout << " ";
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if ((!CheckFile(xyzFilename)) or (!CheckFile(conFilename)))
  {
    cout << "Error: Missing files!!!";
    cout << '\n' << '\n';
    doQuit = 1;
  }
  //Parse input if files exist
  if (!doQuit)
  {
    //Parse XYZ file and connect file
    outFile.open("tinkxyz.xyz",ios_base::out);
    posFile >> Natoms; //Collect the number of atoms
    outFile << Natoms << '\n'; //Write the number of atoms
    for (int i=0;i<Natoms;i++)
    {
      //Collect positions
      double x,y,z;
      posFile >> dummy >> x >> y >> z;
      //Collect atom types
      string atType; //MM character atom type
      int atNum; //MM numerical atom type
      conFile >> dummy >> atType >> atNum;
      //Collect bonds
      int NBonds;
      vector<int> bonds;
      conFile >> dummy >> dummy >> NBonds;
      for (int j=0;j<NBonds;j++)
      {
        int bondID;
        conFile >> bondID;
        bondID += 1; //Fix index for TINKER
        bonds.push_back(bondID);
      }
      //Write line of TINKER XYZ
      outFile << (i+1) << " ";
      outFile << atType << " ";
      outFile << LICHEMFormFloat(x,12) << " ";
      outFile << LICHEMFormFloat(y,12) << " ";
      outFile << LICHEMFormFloat(z,12) << " ";
      outFile << atNum;
      for (unsigned int j=0;j<bonds.size();j++)
      {
        outFile << " "; //Prevents trailing spaces
        outFile << bonds[j];
      }
      outFile << '\n'; //End the line
    }
    //Finish up and exit
    cout << "TINKER XYZ data written to tinkxyz.xyz";
    cout << '\n' << '\n';
    cout.flush();
    outFile.close();
    posFile.close();
    conFile.close();
  }
  exit(0); //Quit
  return;
};

//End of file group
///@}

