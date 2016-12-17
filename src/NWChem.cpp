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

/*!
  \ingroup NWChem
*/
///@{

//QM wrapper functions
void NWChemCharges(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                   int bead)
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
  WriteNWChemInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    //Run in parallel
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << bead << ".nw";
  call << " > LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for energy
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while ((!inFile.eof()) && inFile.good())
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
  call << "LICHM_" << bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
      {
        //Read charges
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> QMMMData[i].MP[bead].q;
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
    call << "rm -f LICHM_" << bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files and return
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead << ".b*" << " ";
  call << "LICHM_" << bead << ".c*" << " ";
  call << "LICHM_" << bead << ".d*" << " ";
  call << "LICHM_" << bead << ".e*" << " ";
  call << "LICHM_" << bead << ".f*" << " ";
  call << "LICHM_" << bead << ".g*" << " ";
  call << "LICHM_" << bead << ".h*" << " ";
  call << "LICHM_" << bead << ".l*" << " ";
  call << "LICHM_" << bead << ".n*" << " ";
  call << "LICHM_" << bead << ".p*" << " ";
  call << "LICHM_" << bead << ".q*" << " ";
  call << "LICHM_" << bead << ".x*" << " ";
  call << "LICHM_" << bead << ".z*";
  globalSys = system(call.str().c_str());
  return;
};

double NWChemEnergy(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
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
  WriteNWChemInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << bead << ".nw";
  call << " > LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for energy
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while ((!inFile.eof()) && inFile.good())
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
  call << "LICHM_" << bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
      {
        //Read charges
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> QMMMData[i].MP[bead].q;
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
    call << "rm -f LICHM_" << bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_";
    call << bead << ".nw " << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << bead << ".trash; ";
    call << "rm -f LICHM_" << bead << ".trash";
    call << "; ";
    call << "cp LICHM_" << bead << ".movecs ";
    call << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << bead << ".trash; ";
    call << "rm -f LICHM_" << bead << ".trash";
    call << "; ";
    call << "cp LICHM_" << bead << ".log ";
    call << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << bead << ".trash; ";
    call << "rm -f LICHM_" << bead << ".trash";
    call << " "; //Extra blank space before the next command
  }
  call << "rm -f ";
  call << "LICHM_" << bead << ".b*" << " ";
  call << "LICHM_" << bead << ".c*" << " ";
  call << "LICHM_" << bead << ".d*" << " ";
  call << "LICHM_" << bead << ".e*" << " ";
  call << "LICHM_" << bead << ".f*" << " ";
  call << "LICHM_" << bead << ".g*" << " ";
  call << "LICHM_" << bead << ".h*" << " ";
  call << "LICHM_" << bead << ".l*" << " ";
  call << "LICHM_" << bead << ".n*" << " ";
  call << "LICHM_" << bead << ".p*" << " ";
  call << "LICHM_" << bead << ".q*" << " ";
  call << "LICHM_" << bead << ".x*" << " ";
  call << "LICHM_" << bead << ".z*";
  globalSys = system(call.str().c_str());
  //Change units and return
  E *= har2eV;
  return E;
};

double NWChemForces(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                    QMMMSettings& QMMMOpts, int bead)
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
  WriteNWChemInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << bead << ".nw";
  call << " > LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for forces and energies
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  bool gradDone = 0;
  while ((!inFile.eof()) && inFile.good())
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
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          //Extract forces; Convoluted, but "easy"
          getline(inFile,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> dummy >> dummy >> dummy; //Clear coordinates
          line >> fX;
          line >> fY;
          line >> fZ;
          //Change from gradient to force
          fX *= -1;
          fY *= -1;
          fZ *= -1;
          //Switch to eV/A and save forces
          forces(3*i) += fX*har2eV/bohrRad;
          forces(3*i+1) += fY*har2eV/bohrRad;
          forces(3*i+2) += fZ*har2eV/bohrRad;
        }
      }
    }
  }
  inFile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
      {
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        line >> QMMMData[i].MP[bead].q;
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
    call << "rm -f LICHM_" << bead << ".movecs";
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
  call << "LICHM_" << bead << ".b*" << " ";
  call << "LICHM_" << bead << ".c*" << " ";
  call << "LICHM_" << bead << ".d*" << " ";
  call << "LICHM_" << bead << ".e*" << " ";
  call << "LICHM_" << bead << ".f*" << " ";
  call << "LICHM_" << bead << ".g*" << " ";
  call << "LICHM_" << bead << ".h*" << " ";
  call << "LICHM_" << bead << ".l*" << " ";
  call << "LICHM_" << bead << ".n*" << " ";
  call << "LICHM_" << bead << ".p*" << " ";
  call << "LICHM_" << bead << ".q*" << " ";
  call << "LICHM_" << bead << ".x*" << " ";
  call << "LICHM_" << bead << ".z*";
  globalSys = system(call.str().c_str());
  //Change units and return
  E *= har2eV;
  return E;
};

MatrixXd NWChemHessian(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                       int bead)
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
  WriteNWChemInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << bead << ".nw";
  call << " > LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output for Hessian
  call.str("");
  call << "LICHM_" << bead << ".hess";
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
    call << "rm -f LICHM_" << bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead << ".b*" << " ";
  call << "LICHM_" << bead << ".c*" << " ";
  call << "LICHM_" << bead << ".d*" << " ";
  call << "LICHM_" << bead << ".e*" << " ";
  call << "LICHM_" << bead << ".f*" << " ";
  call << "LICHM_" << bead << ".g*" << " ";
  call << "LICHM_" << bead << ".h*" << " ";
  call << "LICHM_" << bead << ".l*" << " ";
  call << "LICHM_" << bead << ".n*" << " ";
  call << "LICHM_" << bead << ".p*" << " ";
  call << "LICHM_" << bead << ".q*" << " ";
  call << "LICHM_" << bead << ".x*" << " ";
  call << "LICHM_" << bead << ".z*";
  globalSys = system(call.str().c_str());
  //Return Hessian
  return QMHess;
};

double NWChemOpt(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
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
  WriteNWChemInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run calculation
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem LICHM_" << bead << ".nw";
  call << " > LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Parse output
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while ((!inFile.eof()) && inFile.good())
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
  call << "LICHM_" << bead << ".q";
  inFile.open(call.str().c_str(),ios_base::in);
  if (inFile.good())
  {
    getline(inFile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
      {
        //Read charges
        stringstream line;
        getline(inFile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> QMMMData[i].MP[bead].q;
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
    call << "rm -f LICHM_" << bead << ".movecs";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead << ".b*" << " ";
  call << "LICHM_" << bead << ".c*" << " ";
  call << "LICHM_" << bead << ".d*" << " ";
  call << "LICHM_" << bead << ".e*" << " ";
  call << "LICHM_" << bead << ".f*" << " ";
  call << "LICHM_" << bead << ".g*" << " ";
  call << "LICHM_" << bead << ".h*" << " ";
  call << "LICHM_" << bead << ".l*" << " ";
  call << "LICHM_" << bead << ".n*" << " ";
  call << "LICHM_" << bead << ".p*" << " ";
  call << "LICHM_" << bead << ".q*" << " ";
  call << "LICHM_" << bead << ".x*" << " ";
  call << "LICHM_" << bead << ".z*";
  globalSys = system(call.str().c_str());
  //Change units and return
  E *= har2eV;
  return E;
};

//End of file group
///@}

