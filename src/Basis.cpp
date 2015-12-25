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

//BASIS file creation functions
void LICHEM2BASIS(int& argc,char**& argv)
{
  //Writes a BASIS file based on a LICHEM regions file
  fstream regfile,ofile; //File streams
  string dummy; //Generic string
  regfilename = "NOFILE"; //Global regions filename
  //Read arguments
  Nqm = 0; //For safety
  Npseudo = 0; //For safety
  bool DoQuit = 0; //Exit with an error
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
        DoQuit = 1;
      }
      regfilename = file.str();
      regfile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if (!CheckFile(regfilename))
  {
    cout << "Error: Missing region file!!!";
    cout << '\n' << '\n';
    DoQuit = 1;
  }
  //Parse input
  if (!DoQuit)
  {
    //Parse region file
    bool FoundWrapper = 0;
    bool FoundBasis = 0;
    bool FoundQM = 0;
    bool FoundPB = 0;
    vector<int> QMatoms;
    vector<int> PBatoms;
    string BasisSetName = "???";
    string WrapperName = "N/A";
    ofile.open("BASIS",ios_base::out);
    if (regfile.good())
    {
      //Loop over the lines
      while (!regfile.eof())
      {
        //Look for key words
        stringstream line; //Generic stream
        getline(regfile,dummy); //Read the line
        line << dummy; //Copy to a better container
        line >> dummy; //Read first item in the string
        LICHEMLowerText(dummy);
        if (dummy == "qm_type:")
        {
          //Read basis set information
          line >> WrapperName;
          LICHEMLowerText(WrapperName);
          FoundWrapper = 1;
        }
        if (dummy == "qm_basis:")
        {
          //Read basis set information
          line >> BasisSetName;
          FoundBasis = 1;
        }
        if (dummy == "qm_atoms:")
        {
          //Read the list QM atoms
          line >> Nqm;
          for (int i=0;i<Nqm;i++)
          {
            //Save ID to the array
            int AtNum = 0;
            regfile >> AtNum;
            QMatoms.push_back(AtNum);
          }
          FoundQM = 1;
        }
        if (dummy == "pseudobond_atoms:")
        {
          //Read the list of PB atoms
          line >> Npseudo;
          for (int i=0;i<Npseudo;i++)
          {
            //Save ID to the array
            int AtNum = 0;
            regfile >> AtNum;
            PBatoms.push_back(AtNum);
          }
          FoundPB = 1;
        }
      }
    }
    //Check for read errors
    if ((!FoundWrapper) or (!FoundBasis) or (!FoundQM) or (!FoundPB))
    {
      cout << "Error: Missing data in region file!!!";
      cout << '\n';
      cout << " Check keywords and try again...";
      cout << '\n' << '\n';
      exit(0);
    }
    //Modify atom order
    int Ntotal = Nqm+Npseudo; //Total number of QM and PB atoms
    vector<bool> AtomList; //Region flags for all QM and PB atoms
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
      bool InQMregion = 1;
      for (int j=0;j<Npseudo;j++)
      {
        if (PBatoms[j] == QMPBatoms[i])
        {
          //Atom is in the PB region
          InQMregion = 0;
        }
      }
      AtomList.push_back(InQMregion);
    }
    //Write BASIS file
    if ((WrapperName == "gaussian") or (WrapperName == "g09"))
    {
      //Write Gaussian BASIS file for GEN input
      int ct; //Generic counter
      //Write QM atoms
      ct = 0; //Reset counter
      for (int i=0;i<Ntotal;i++)
      {
        if (AtomList[i])
        {
          ct += 1; //Increase counter
          ofile << (i+1) << " ";
          if (ct == 8)
          {
            //Print basis info
            ofile << " 0" << '\n';
            ofile << BasisSetName << '\n';
            ofile << "****" << '\n';
            //Start a new line
            ct = 0;
          }
        }
      }
      //Terminate QM basis set information
      if (ct != 0)
      {
        //Print basis info
        ofile << " 0" << '\n';
        ofile << BasisSetName << '\n';
        ofile << "****" << '\n';
      }
      //Write PB atoms
      ct = 0; //Reset counter
      for (int i=0;i<Ntotal;i++)
      {
        if (!AtomList[i])
        {
          ct += 1; //Increase counter
          ofile << (i+1) << " ";
          if (ct == 8)
          {
            //Print basis info
            ofile << " 0" << '\n';
            ofile << "[PB basis set]" << '\n';
            ofile << "****" << '\n';
            //Start a new line
            ct = 0;
          }
        }
      }
      //Terminate QM basis set information
      if (ct != 0)
      {
        //Print basis info
        ofile << " 0" << '\n';
        ofile << "[PB basis set]" << '\n';
        ofile << "****" << '\n';
      }
      //Write pseudopotential information
      ofile << '\n';
      if (Npseudo > 0)
      {
        ct = 0; //Reset counter
        for (int i=0;i<Ntotal;i++)
        {
          if (!AtomList[i])
          {
            ct += 1; //Increase counter
            ofile << (i+1) << " ";
            if (ct == 8)
            {
              //Print basis info
              ofile << " 0" << '\n';
              ofile << "[PB pseudopotentials]" << '\n';
              //Start a new line
              ct = 0;
            }
          }
        }
        //Terminate pseudopotential information
        if (ct != 0)
        {
          //Print basis info
          ofile << " 0" << '\n';
          ofile << "[PB pseudopotentials]" << '\n';
        }
        ofile << '\n';
      }
    }
    else if (WrapperName == "nwchem")
    {
      //Write NWChem BASIS file
      ofile << "basis" << '\n';
      ofile << " * library ";
      ofile << BasisSetName;
      ofile << " except F2pb" << '\n';
      if (Npseudo > 0)
      {
        ofile << " F2pb [PB basis set]" << '\n';
      }
      ofile << "end" << '\n';
      if (Npseudo > 0)
      {
        ofile << "ecp" << '\n';
        ofile << " F2pb nelec 2" << '\n';
        ofile << " F2pb [Pseudopotential type]" << '\n';
        ofile << "  [PB pseudopotential]" << '\n';
        ofile << "end" << '\n';
      }
    }
    else if (WrapperName == "psi4")
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
    ofile.flush();
    ofile.close();
  }
  //Quit
  regfile.close();
  exit(0);
  return;
};

