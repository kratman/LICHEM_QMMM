/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Functions for creating and modifying basis sets for specific QM packages

 Reference for Gaussian:
 Frisch et al., Gaussian 09 Rev D.01, (2009)

 Reference for PSI4:
 Turney et al., WIREs Comp. Mol. Sci., 2, 4, 556, (2012)

 Reference for NWChem:
 Valiev et al., Comput. Phys. Commun., 181, 1477, (2010)

*/

/*!
  \ingroup Misc
*/
//! \{

//BASIS file creation functions

//! \brief Creates a Gaussian or NWChem basis set file.
//! \param argc - Number of arguments passed to LICHEM
//! \param argv - All arguments passed to LICHEM
void LICHEM2BASIS(int& argc,char**& argv)
{
  //Writes a BASIS file based on a LICHEM regions file
  fstream regFile,outFile; //File streams
  string dummy; //Generic string
  regFilename = "NOFILE"; //Global regions filename
  //Read arguments
  Nqm = 0; //For safety
  Npseudo = 0; //For safety
  bool doQuit = 0; //Exit with an error
  cout << "Reading LICHEM input: ";
  for (int i=0;i<argc;i++)
  {
    dummy = string(argv[i]);
    //Check regions file
    if (dummy == "-b")
    {
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open regions file!!!";
        cout << '\n';
        doQuit = 1;
      }
      regFilename = file.str();
      regFile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if (!CheckFile(regFilename))
  {
    //Missing flag
    cout << "Error: Missing region file!!!";
    cout << '\n' << '\n';
    doQuit = 1;
  }
  //Parse input
  if (!doQuit)
  {
    //Parse region file
    bool foundWrapper = 0;
    bool foundBasis = 0;
    bool foundQM = 0;
    vector<int> QMatoms;
    vector<int> PBatoms;
    string basisSetName = "???";
    string wrapperName = "N/A";
    outFile.open("BASIS",ios_base::out);
    if (regFile.good())
    {
      //Loop over the lines
      while (!regFile.eof())
      {
        //Look for key words
        stringstream line; //Generic stream
        getline(regFile,dummy); //Read the line
        line << dummy; //Copy to a better container
        line >> dummy; //Read first item in the string
        LICHEMLowerText(dummy);
        if (dummy == "qm_type:")
        {
          //Read basis set information
          line >> wrapperName;
          LICHEMLowerText(wrapperName);
          foundWrapper = 1;
        }
        if (dummy == "qm_basis:")
        {
          //Read basis set information
          line >> basisSetName;
          foundBasis = 1;
        }
        if (dummy == "qm_atoms:")
        {
          //Read the list QM atoms
          line >> Nqm;
          for (int i=0;i<Nqm;i++)
          {
            //Save ID to the array
            int atNum = 0;
            regFile >> atNum;
            QMatoms.push_back(atNum);
          }
          foundQM = 1;
        }
        if (dummy == "pseudobond_atoms:")
        {
          //Read the list of PB atoms
          line >> Npseudo;
          for (int i=0;i<Npseudo;i++)
          {
            //Save ID to the array
            int atNum = 0;
            regFile >> atNum;
            PBatoms.push_back(atNum);
          }
        }
      }
    }
    //Check for read errors
    if ((!foundWrapper) || (!foundBasis) || (!foundQM))
    {
      cout << "Error: Missing data in region file!!!";
      cout << '\n';
      cout << " Check keywords and try again...";
      cout << '\n' << '\n';
      exit(0);
    }
    //Modify atom order
    int Ntotal = Nqm+Npseudo; //Total number of QM and PB atoms
    vector<bool> atomList; //Region flags for all QM and PB atoms
    vector<int> QMPBatoms; //List of all QM and PB atoms
    for (int i=0;i<Nqm;i++)
    {
      //Add QM atoms
      QMPBatoms.push_back(QMatoms[i]);
    }
    for (int i=0;i<Npseudo;i++)
    {
      //Add PB atoms
      QMPBatoms.push_back(PBatoms[i]);
    }
    sort(QMPBatoms.begin(),QMPBatoms.end()); //Put atoms in ascending order
    for (int i=0;i<Ntotal;i++)
    {
      //Check regions
      bool inQMregion = 1;
      for (int j=0;j<Npseudo;j++)
      {
        if (PBatoms[j] == QMPBatoms[i])
        {
          //Atom is in the PB region
          inQMregion = 0;
        }
      }
      atomList.push_back(inQMregion);
    }
    //Write BASIS file
    if ((wrapperName == "gaussian") || (wrapperName == "g09"))
    {
      //Write Gaussian BASIS file for GEN input
      int ct; //Generic counter
      //Write QM atoms
      ct = 0; //Reset counter
      for (int i=0;i<Ntotal;i++)
      {
        if (atomList[i])
        {
          ct += 1; //Increase counter
          outFile << (i+1) << " ";
          if (ct == 8)
          {
            //Print basis info
            outFile << " 0" << '\n';
            outFile << basisSetName << '\n';
            outFile << "****" << '\n';
            //Start a new line
            ct = 0;
          }
        }
      }
      //Terminate QM basis set information
      if (ct != 0)
      {
        //Print basis info
        outFile << " 0" << '\n';
        outFile << basisSetName << '\n';
        outFile << "****" << '\n';
      }
      //Write PB atoms
      ct = 0; //Reset counter
      if (Npseudo > 0)
      {
        for (int i=0;i<Ntotal;i++)
        {
          if (!atomList[i])
          {
            ct += 1; //Increase counter
            outFile << (i+1) << " ";
            if (ct == 8)
            {
              //Print basis info
              outFile << " 0" << '\n';
              outFile << "[PB basis set]" << '\n';
              outFile << "****" << '\n';
              //Start a new line
              ct = 0;
            }
          }
        }
        //Terminate QM basis set information
        if (ct != 0)
        {
          //Print basis info
          outFile << " 0" << '\n';
          outFile << "[PB basis set]" << '\n';
          outFile << "****" << '\n';
        }
      }
      //Write pseudopotential information
      outFile << '\n';
      if (Npseudo > 0)
      {
        ct = 0; //Reset counter
        for (int i=0;i<Ntotal;i++)
        {
          if (!atomList[i])
          {
            ct += 1; //Increase counter
            outFile << (i+1) << " ";
            if (ct == 8)
            {
              //Print basis info
              outFile << " 0" << '\n';
              outFile << "[PB pseudopotentials]" << '\n';
              //Start a new line
              ct = 0;
            }
          }
        }
        //Terminate pseudopotential information
        if (ct != 0)
        {
          //Print basis info
          outFile << " 0" << '\n';
          outFile << "[PB pseudopotentials]" << '\n';
        }
        outFile << '\n';
      }
    }
    else if (wrapperName == "nwchem")
    {
      //Write NWChem BASIS file
      outFile << "basis" << '\n';
      outFile << " * library ";
      outFile << basisSetName;
      outFile << " except F2pb" << '\n';
      if (Npseudo > 0)
      {
        outFile << " F2pb [PB basis set]" << '\n';
      }
      outFile << "end" << '\n';
      if (Npseudo > 0)
      {
        outFile << "ecp" << '\n';
        outFile << " F2pb nelec 2" << '\n';
        outFile << " F2pb [Pseudopotential type]" << '\n';
        outFile << "  [PB pseudopotential]" << '\n';
        outFile << "end" << '\n';
      }
    }
    else if (wrapperName == "psi4")
    {
      cout << "Error: Pseudopotentials are not yet implemented in PSI4.";
      cout << '\n' << '\n';
      exit(0);
    }
    else
    {
      cout << "Error: Unrecognized QM package!!!";
      cout << '\n' << '\n';
      exit(0);
    }
    //Finish up and exit
    cout << "Basis set data written to BASIS";
    cout << '\n' << '\n';
    cout.flush();
    outFile.flush();
    outFile.close();
  }
  //Quit
  regFile.close();
  exit(0);
  return;
};

//End of file group
//! \}

