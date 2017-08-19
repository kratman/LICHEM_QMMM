/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for PSI4.

 Reference for PSI4:
 Turney et al., WIREs Comp. Mol. Sci., 2, 4, 556, (2012)

*/

/*!
  \ingroup PSI4
*/
///@{

//QM wrapper functions
void PSI4Charges(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
{
  //Function to update QM point-charges
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  //Check if there is a checkpoint file
  bool useCheckPoint;
  call.str("");
  call << "LICHM_" << bead << ".180";
  useCheckPoint = CheckFile(call.str());
  //Set up charge calculation
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (useCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  if (QMMM)
  {
    call << "oeprop(qmwfn,'MULLIKEN_CHARGES')" << '\n';
  }
  WritePSI4Input(QMMMData,call.str(),QMMMOpts,bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << " -i ";
  call << "LICHM_" << bead << ".dat -o ";
  call << "LICHM_" << bead << ".out > ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << bead << ".180 ";
  call << "LICHM_" << bead << ".180 ";
  call << "2> LICHM_" << bead << ".trash; ";
  call << "rm -f LICHM_" << bead << ".trash";
  globalSys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "LICHM_" << bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(inFile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
  }
  inFile.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead << ".dat ";
  call << "LICHM_" << bead << ".out ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  return;
};

double PSI4Energy(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts, int bead)
{
  //Runs PSI4 for energy calculations
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0.0;
  //Check if there is a checkpoint file
  bool useCheckPoint;
  call.str("");
  call << "LICHM_" << bead << ".180";
  useCheckPoint = CheckFile(call.str());
  //Set up energy calculation
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (useCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  if (QMMM)
  {
    call << "oeprop(qmwfn,'MULLIKEN_CHARGES')" << '\n';
  }
  WritePSI4Input(QMMMData,call.str(),QMMMOpts,bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << " -i ";
  call << "LICHM_" << bead << ".dat -o ";
  call << "LICHM_" << bead << ".out > ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << bead << ".180 ";
  call << "LICHM_" << bead << ".180 ";
  call << "2> LICHM_" << bead << ".trash; ";
  call << "rm -f LICHM_" << bead << ".trash";
  globalSys = system(call.str().c_str());
  //Read energy
  call.str("");
  call << "LICHM_" << bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(inFile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    line >> dummy; //Get rid of junk
    if (dummy == "Final")
    {
      line >> dummy; //Check property
      if (dummy == "Energy:")
      {
        //Read energy
        line >> E; //Read energy
        QMFinished = 1;
      }
    }
  }
  inFile.close();
  //Collect energy (post-SCF)
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Energy:")
    {
      line >> E; //Read post-SCF energy
      QMFinished = 1;
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
    //Delete checkpoint
    call.str("");
    call << "rm -f *LICHM_" << bead << ".180";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_" << bead << ".* ";
    call << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << bead << ".trash; ";
    call << "rm -f LICHM_" << bead << ".trash";
    call << " "; //Extra blank space before the next command
  }
  call << "rm -f ";
  call << "LICHM_" << bead << ".dat ";
  call << "LICHM_" << bead << ".out ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Change units
  E *= har2eV;
  return E;
};

double PSI4Forces(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                  QMMMSettings& QMMMOpts, int bead)
{
  //Function for calculating the forces and charges on a set of atoms
  fstream inFile; //Generic file name
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0;
  //Check if there is a checkpoint file
  bool useCheckPoint;
  call.str("");
  call << "LICHM_" << bead << ".180";
  useCheckPoint = CheckFile(call.str());
  //Set up force calculation
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (useCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  call << "gradient('" << QMMMOpts.func << "'";
  call << ",bypass_scf=True)"; //Skip the extra SCF cycle
  call << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  if (QMMM)
  {
    call << "oeprop(qmwfn,'MULLIKEN_CHARGES')" << '\n';
  }
  WritePSI4Input(QMMMData,call.str(),QMMMOpts,bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << " -i ";
  call << "LICHM_" << bead << ".dat -o ";
  call << "LICHM_" << bead << ".out > ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << bead << ".180 ";
  call << "LICHM_" << bead << ".180 ";
  call << "2> LICHM_" << bead << ".trash; ";
  call << "rm -f LICHM_" << bead << ".trash";
  globalSys = system(call.str().c_str());
  //Extract forces
  call.str("");
  call << "LICHM_" << bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(inFile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    if (dummy == "-Total")
    {
      line >> dummy;
      if (dummy == "Gradient:")
      {
        getline(inFile,dummy);
        getline(inFile,dummy);
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          //Extract forces; Convoluted, but "easy"
          getline(inFile,dummy);
          stringstream line(dummy);
          line >> dummy; //Clear junk
          line >> fX; //Read value
          line >> fY; //Read value
          line >> fZ; //Read value
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
    line >> dummy; //Get rid of junk
    if (dummy == "Final")
    {
      line >> dummy; //Check property
      if (dummy == "Energy:")
      {
        line >> E; //Read energy
      }
    }
  }
  inFile.close();
  //Collect energy (post-SCF)
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy; //Check property
    if (dummy == "Energy:")
    {
      line >> E; //Read post-SCF energy
    }
  }
  inFile.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead << ".dat ";
  call << "LICHM_" << bead << ".out ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Change units
  E *= har2eV;
  return E;
};

MatrixXd PSI4Hessian(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                     int bead)
{
  //Runs PSI4 to calculate a Hessian
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  int Ndof = 3*(Nqm+Npseudo);
  MatrixXd QMHess(Ndof,Ndof);
  QMHess.setZero();
  //Check if there is a checkpoint file
  bool useCheckPoint;
  call.str("");
  call << "LICHM_" << bead << ".180";
  useCheckPoint = CheckFile(call.str());
  //Calculate Hessian
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (useCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  call << "QMHess = hessian('" << QMMMOpts.func << "'";
  call << ",bypass_scf=True)"; //Skip the extra SCF cycle
  call << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  call << "QMHess.print_out()" << '\n';
  if (QMMM)
  {
    call << "oeprop(qmwfn,'MULLIKEN_CHARGES')" << '\n';
  }
  WritePSI4Input(QMMMData,call.str(),QMMMOpts,bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << " -i ";
  call << "LICHM_" << bead << ".dat -o ";
  call << "LICHM_" << bead << ".out > ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << bead << ".180 ";
  call << "LICHM_" << bead << ".180 ";
  call << "2> LICHM_" << bead << ".trash; ";
  call << "rm -f LICHM_" << bead << ".trash";
  globalSys = system(call.str().c_str());
  //Extract Hessian
  call.str("");
  call << "LICHM_" << bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  bool hessDone = 0;
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(inFile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    line >> dummy; //Get rid of junk
    if (dummy == "Hessian")
    {
      hessDone = 1;
      getline(inFile,dummy); //Clear junk
      //Read Hessian in groups of Ndofx5
      int rowCt = 0; //Current row ID
      while (rowCt < Ndof)
      {
        //Clear junk
        getline(inFile,dummy); //Clear junk
        getline(inFile,dummy); //Clear junk
        getline(inFile,dummy); //Clear junk
        //Read elements
        for (int i=0;i<Ndof;i++)
        {
          inFile >> dummy; //Clear junk
          for (int j=0;j<5;j++)
          {
            if ((rowCt+j) < Ndof)
            {
              inFile >> QMHess(i,rowCt+j);
            }
          }
        }
        //Update row location
        rowCt += 5;
      }
    }
  }
  inFile.close();
  //Check for errors
  if (!hessDone)
  {
    //Calculation did not finish
    cerr << "Error: No force constants recovered!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << bead << ".180";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead << ".dat ";
  call << "LICHM_" << bead << ".out ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  return QMHess;
};

double PSI4Opt(vector<QMMMAtom>& QMMMData,
               QMMMSettings& QMMMOpts, int bead)
{
  //Runs PSI4 for pure QM optimizations
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0.0;
  //Set up QM only optimization
  call.str("");
  call << "optimize('" << QMMMOpts.func << "'";
  call << ")" << '\n';
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  call << ",return_wfn=True)" << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  if (QMMM)
  {
    call << "oeprop(qmwfn,'MULLIKEN_CHARGES')" << '\n';
  }
  WritePSI4Input(QMMMData,call.str(),QMMMOpts,bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << " -i ";
  call << "LICHM_" << bead << ".dat -o ";
  call << "LICHM_" << bead << ".out > ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << bead << ".180 ";
  call << "LICHM_" << bead << ".180 ";
  call << "2> LICHM_" << bead << ".trash; ";
  call << "rm -f LICHM_" << bead << ".trash";
  globalSys = system(call.str().c_str());
  //Read energy and structure
  call.str("");
  call << "LICHM_" << bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  bool optFinished = 0;
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Final")
    {
      line >> dummy;
      if (dummy == "optimized")
      {
        //Read new geometry
        optFinished = 1;
        for (int i=0;i<5;i++)
        {
          //Clear junk
          getline(inFile,dummy);
        }
        for (int i=0;i<Natoms;i++)
        {
          getline(inFile,dummy);
          stringstream line(dummy);
          line >> dummy; //Clear atom type
          double x,y,z;
          line >> x >> y >> z;
          QMMMData[i].P[bead].x = x;
          QMMMData[i].P[bead].y = y;
          QMMMData[i].P[bead].z = z;
        }
      }
    }
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(inFile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    line >> dummy; //Get rid of junk
    if (dummy == "Final")
    {
      line >> dummy; //Check value
      if (dummy == "Energy:")
      {
        //Read energy
        line >> E; //Read energy
        QMFinished = 1;
      }
    }
  }
  inFile.close();
  //Collect energy (post-SCF)
  call.str("");
  call << "LICHM_" << bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  while ((!inFile.eof()) && inFile.good())
  {
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy; //Check value
    if (dummy == "Energy:")
    {
      line >> E; //Read post-SCF energy
      QMFinished = 1;
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
    //Delete checkpoint
    call.str("");
    call << "rm -f *LICHM_" << bead << ".180";
    globalSys = system(call.str().c_str());
  }
  if (!optFinished)
  {
    cerr << "Warning: Optimization did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue using the";
    cerr << " old structure...";
    cerr << '\n';
    E = hugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f *LICHM_" << bead << ".180";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead << ".dat ";
  call << "LICHM_" << bead << ".out ";
  call << "LICHM_" << bead << ".log";
  globalSys = system(call.str().c_str());
  //Change units
  E *= har2eV;
  return E;
};

//End of file group
///@}

