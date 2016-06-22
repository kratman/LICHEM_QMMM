/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for NWChem.

 Reference for NWChem:
 Valiev et al., Comput. Phys. Commun., 181, 1477, (2010)

*/

//QM utility functions


//QM wrapper functions
void NWChemCharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                   int Bead)
{
  //Calculates atomic charges with NWChem
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  //Write NWChem input
  call.str("");
  call << "task dft energy" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    //Run in parallel
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << Bead << ".nw";
  call << " > LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for energy
  call.str("");
  call << "LICHM_" << Bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while (!inFile.eof())
  {
    stringstream line;
    getline(inFile,dummy);
    line.str(dummy);
    line >> dummy;
    //Search for energy
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> E; //Read energy
        QMFinished = 1;
      }
    }
  }
  inFile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Read charges
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> Struct[i].MP[Bead].q;
      }
    }
  }
  inFile.close();
  //Check for errors
  if (!QMFinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = hugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files and return
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".d*" << " ";
  call << "LICHM_" << Bead << ".f*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".h*" << " ";
  call << "LICHM_" << Bead << ".l*" << " ";
  call << "LICHM_" << Bead << ".n*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".x*" << " ";
  call << "LICHM_" << Bead << ".z*";
  globalSys = system(call.str().c_str());
  return;
};

double NWChemEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                    int Bead)
{
  //Runs NWChem energy calculations
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  //Write NWChem input
  call.str("");
  call << "task dft energy" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << Bead << ".nw";
  call << " > LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for energy
  call.str("");
  call << "LICHM_" << Bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while (!inFile.eof())
  {
    stringstream line;
    getline(inFile,dummy);
    line.str(dummy);
    line >> dummy;
    //Search for energy
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> E; //Read energy
        QMFinished = 1;
      }
    }
  }
  inFile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Read charges
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> Struct[i].MP[Bead].q;
      }
    }
  }
  inFile.close();
  //Check for errors
  if (!QMFinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = hugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_";
    call << Bead << ".nw " << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << "; ";
    call << "cp LICHM_" << Bead << ".movecs ";
    call << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << "; ";
    call << "cp LICHM_" << Bead << ".log ";
    call << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << " "; //Extra blank space before the next command
  }
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".d*" << " ";
  call << "LICHM_" << Bead << ".f*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".h*" << " ";
  call << "LICHM_" << Bead << ".l*" << " ";
  call << "LICHM_" << Bead << ".n*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".x*" << " ";
  call << "LICHM_" << Bead << ".z*";
  globalSys = system(call.str().c_str());
  //Change units and return
  E *= har2eV;
  return E;
};

double NWChemForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                    QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem force calculations
  fstream inFile; //Generic file stream
  string dummy; //Genric string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0;
  //Set up force calculation
  call.str("");
  call << "task dft gradient" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << Bead << ".nw";
  call << " > LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for forces and energies
  call.str("");
  call << "LICHM_" << Bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  bool gradDone = 0;
  while (!inFile.eof())
  {
    stringstream line;
    getline(inFile,dummy);
    line.str(dummy);
    line >> dummy;
    //Search for energy
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> E; //Read energy
        QMFinished = 1;
      }
    }
    if (dummy == "atom")
    {
      line >> dummy >> dummy;
      if (dummy == "gradient")
      {
        gradDone = 1; //Not grad school, that lasts forever
        getline(inFile,dummy); //Clear junk
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          //Initialize temporary force variables
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Extract forces; Convoluted, but "easy"
          getline(inFile,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> dummy >> dummy >> dummy; //Clear coordinates
          line >> Fx;
          line >> Fy;
          line >> Fz;
          //Change from gradient to force
          Fx *= -1;
          Fy *= -1;
          Fz *= -1;
          //Switch to eV/A and save forces
          Forces(3*i) += Fx*har2eV/bohrRad;
          Forces(3*i+1) += Fy*har2eV/bohrRad;
          Forces(3*i+2) += Fz*har2eV/bohrRad;
        }
      }
    }
  }
  inFile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        line >> Struct[i].MP[Bead].q;
      }
    }
  }
  inFile.close();
  //Check for errors
  if (!QMFinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = hugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  if (!gradDone)
  {
    cerr << "Warning: No forces recovered!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to recover...";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".d*" << " ";
  call << "LICHM_" << Bead << ".f*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".h*" << " ";
  call << "LICHM_" << Bead << ".l*" << " ";
  call << "LICHM_" << Bead << ".n*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".x*" << " ";
  call << "LICHM_" << Bead << ".z*";
  globalSys = system(call.str().c_str());
  //Change units and return
  E *= har2eV;
  return E;
};

MatrixXd NWChemHessian(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                       int Bead)
{
  //Function to calculate the QM Hessian
  fstream QMlog; //Generic file stream
  string dummy; //Genric string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  int Ndof = 3*(Nqm+Npseudo);
  MatrixXd QMHess(Ndof,Ndof);
  QMHess.setZero();
  //Set up Hessian calculation
  call.str("");
  call << "task dft hessian" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << Bead << ".nw";
  call << " > LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for Hessian
  call.str("");
  call << "LICHM_" << Bead << ".hess";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool hessDone = 0;
  if (QMlog.good())
  {
    hessDone = 1;
    for (int i=0;i<Ndof;i++)
    {
      for (int j=0;j<(i+1);j++)
      {
        //Read matrix element in scientific notation
        string matelmt;
        QMlog >> matelmt; //Save as a string
        //Change D notation to E notation
        LICHEMFixSciNot(matelmt);
        //Convert to a double and save matrix element
        QMHess(i,j) = atof(matelmt.c_str());
        //Apply symmetry
        QMHess(j,i) = QMHess(i,j);
      }
    }
  }
  QMlog.close();
  //Check for errors
  if (!hessDone)
  {
    //Calculation did not finish
    cerr << "Error: No force constants recovered!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".d*" << " ";
  call << "LICHM_" << Bead << ".f*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".h*" << " ";
  call << "LICHM_" << Bead << ".l*" << " ";
  call << "LICHM_" << Bead << ".n*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".x*" << " ";
  call << "LICHM_" << Bead << ".z*";
  globalSys = system(call.str().c_str());
  //Return Hessian
  return QMHess;
};

double NWChemOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem optimizations
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Save print settings
  double E = 0.0; //QM energy
  //Write NWChem input
  call.str("");
  call << "task dft optimize" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << Bead << ".nw";
  call << " > LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output
  call.str("");
  call << "LICHM_" << Bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while (!inFile.eof())
  {
    stringstream line;
    getline(inFile,dummy);
    line.str(dummy);
    line >> dummy;
    //Search for energy
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> E; //Read energy
        QMFinished = 1;
      }
    }
  }
  inFile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Read charges
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> Struct[i].MP[Bead].q;
      }
    }
  }
  inFile.close();
  //Check for errors
  if (!QMFinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = hugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".d*" << " ";
  call << "LICHM_" << Bead << ".f*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".h*" << " ";
  call << "LICHM_" << Bead << ".l*" << " ";
  call << "LICHM_" << Bead << ".n*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".x*" << " ";
  call << "LICHM_" << Bead << ".z*";
  globalSys = system(call.str().c_str());
  //Change units and return
  E *= har2eV;
  return E;
};

