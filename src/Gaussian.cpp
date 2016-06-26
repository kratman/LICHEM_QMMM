/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for Gaussian. These routines are written for g09
 Rev D. Note that the external function needs to call MM codes.

 Reference for Gaussian:
 Frisch et al., Gaussian 09 Rev D.01, (2009)

*/

//QM utility functions


//QM wrapper functions
void GaussianCharges(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                     int bead)
{
  //Function to update QM point-charges
  fstream QMlog; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".chk";
  bool useCheckPoint = CheckFile(call.str());
  if (QMMMOpts.func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    useCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.func << "/"; //Print the method
  }
  call << QMMMOpts.basis << " SP Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)" << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (useCheckPoint)
    {
      //Read pseudo potential
      call << "Guess=TCheck ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run QM calculation
  call.str("");
  call << "g09 LICHM_" << bead;
  globalSys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "LICHM_" << bead << ".log";
  QMlog.open(call.str().c_str(),ios_base::in);
  while (!QMlog.eof())
  {
    getline(QMlog,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      //Mulliken charges (fallback)
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMlog,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Count through all atoms in the QM calculations
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    if (dummy == "ESP")
    {
      //ESP (MK) charges
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMlog,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Count through all atoms in the QM calculations
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
  }
  QMlog.close();
  //Clean up files and save checkpoint file
  call.str("");
  call << "mv LICHM_" << bead;
  call << ".chk tmp_" << bead;
  call << ".chk; ";
  call << "rm -f ";
  call << "LICHM_" << bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << bead;
  call << ".com";
  call << "; mv tmp_" << bead;
  call << ".chk LICHM_" << bead;
  call << ".chk";
  globalSys = system(call.str().c_str());
  return;
};

double GaussianEnergy(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                      int bead)
{
  //Calculates the QM energy with Gaussian
  fstream QMlog; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  double Eself = 0.0; //External field self-energy
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".chk";
  bool useCheckPoint = CheckFile(call.str());
  if (QMMMOpts.func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    useCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.func << "/"; //Print the method
  }
  call << QMMMOpts.basis << " SP Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)";
  if (useCheckPoint)
  {
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(QMMMData,call.str(),QMMMOpts,bead);
  //Calculate energy
  call.str("");
  call << "g09 ";
  call << "LICHM_" << bead;
  globalSys = system(call.str().c_str());
  //Read output
  call.str("");
  call << "LICHM_" << bead << ".log";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while (!QMlog.eof())
  {
    stringstream line;
    getline(QMlog,dummy);
    line.str(dummy);
    line >> dummy;
    //Search for field self-energy
    if (dummy == "Self")
    {
      line >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> Eself; //Actual self-energy of the charges
      }
    }
    //Search for energy
    if (dummy == "SCF")
    {
      line >> dummy;
      if (dummy == "Done:")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> E; //QM energy
        QMFinished = 1;
      }
    }
    //Check for charges
    if (dummy == "Mulliken")
    {
      //Mulliken charges (fallback)
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMlog,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Count through all atoms in the QM calculations
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    if (dummy == "ESP")
    {
      //ESP (MK) charges
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMlog,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Count through all atoms in the QM calculations
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
  }
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
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  QMlog.close();
  //Clean up files and save checkpoint file
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_";
    call << bead << ".* ";
    call << QMMMOpts.backDir;
    call << "/.";
    call << " 2> LICHM_" << bead << ".trash; ";
    call << "rm -f LICHM_" << bead << ".trash";
    call << " "; //Extra blank space before the next command
  }
  call << "rm -f ";
  call << "LICHM_" << bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << bead;
  call << ".com";
  globalSys = system(call.str().c_str());
  //Change units and return
  E -= Eself;
  E *= har2eV;
  return E;
};

double GaussianForces(vector<QMMMAtom>& QMMMData, VectorXd& forces,
                      QMMMSettings& QMMMOpts, int bead)
{
  //Function for calculating the forces on a set of atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream QMlog; //Generic input files
  double Eqm = 0; //QM energy
  double Eself = 0; //External field self-energy
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".chk";
  bool useCheckPoint = CheckFile(call.str());
  if (QMMMOpts.func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    useCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.func << "/"; //Print the method
  }
  call << QMMMOpts.basis << " Force=NoStep Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)"; //Line ended further below
  if (useCheckPoint)
  {
    //Restart if possible
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run Gaussian
  call.str("");
  call << "g09 " << "LICHM_" << bead;
  globalSys = system(call.str().c_str());
  //Extract forces
  call.str("");
  call << "LICHM_" << bead << ".log";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool gradDone = 0;
  while ((!QMlog.eof()) and (!gradDone))
  {
    //Parse file line by line
    getline(QMlog,dummy);
    stringstream line(dummy);
    line >> dummy;
    //This only works with #P
    if (dummy == "Center")
    {
      line >> dummy >> dummy;
      if (dummy == "Forces")
      {
        gradDone = 1; //Not grad school, that lasts forever
        getline(QMlog,dummy); //Clear junk
        getline(QMlog,dummy); //Clear more junk
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double fX = 0;
          double fY = 0;
          double fZ = 0;
          //Extract forces; Convoluted, but "easy"
          getline(QMlog,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> fX;
          line >> fY;
          line >> fZ;
          //Save forces
          forces(3*i) += fX*har2eV/bohrRad;
          forces(3*i+1) += fY*har2eV/bohrRad;
          forces(3*i+2) += fZ*har2eV/bohrRad;
        }
      }
    }
    if (dummy == "Self")
    {
      line >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> dummy; //Ditto
        line >> Eself; //Actual self-energy of the charges
      }
    }
    //Check for partial QMMM energy
    if (dummy == "SCF")
    {
      line >> dummy;
      if (dummy == "Done:")
      {
        line >> dummy; //Clear junk
        line >> dummy; //Ditto
        line >> Eqm; //QM energy
      }
    }
    //Check for charges
    if (dummy == "Mulliken")
    {
      //Mulliken charges (fallback)
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMlog,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Count through all atoms in the QM calculations
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
    if (dummy == "ESP")
    {
      //ESP (MK) charges
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(QMlog,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Count through all atoms in the QM calculations
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> QMMMData[i].MP[bead].q;
          }
        }
      }
    }
  }
  QMlog.close();
  //Check for errors
  if (!gradDone)
  {
    cerr << "Warning: No forces recovered!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to recover...";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << bead;
  call << ".com";
  globalSys = system(call.str().c_str());
  //Change units and return
  Eqm -= Eself;
  Eqm *= har2eV;
  return Eqm;
};

MatrixXd GaussianHessian(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                         int bead)
{
  //Function for calculating the Hessian for a set of QM atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream QMlog; //Generic input files
  int Ndof = 3*(Nqm+Npseudo);
  MatrixXd QMHess(Ndof,Ndof);
  QMHess.setZero();
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".chk";
  bool useCheckPoint = CheckFile(call.str());
  if (QMMMOpts.func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    useCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.func << "/"; //Print the method
  }
  call << QMMMOpts.basis << " Freq Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)"; //Line ended further below
  if (useCheckPoint)
  {
    //Restart if possible
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(QMMMData,call.str(),QMMMOpts,bead);
  //Run Gaussian
  call.str("");
  call << "g09 " << "LICHM_" << bead;
  globalSys = system(call.str().c_str());
  //Generate formatted checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".chk";
  if (CheckFile(call.str()))
  {
    //Run formchk
    call.str("");
    call << "formchk ";
    call << "LICHM_" << bead << ".chk";
    call << " > LICHM_" << bead << ".trash; ";
    call << " rm -f LICHM_" << bead << ".trash";
    globalSys = system(call.str().c_str());
  }
  else
  {
    //Calculation did not finish
    cerr << "Error: No checkpoint file found!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Open checkpoint file
  call.str("");
  call << "LICHM_" << bead << ".fchk";
  QMlog.open(call.str().c_str(),ios_base::in);
  //Extract Hessian
  bool hessDone = 0;
  while (QMlog.good() and (!QMlog.eof()) and (!hessDone))
  {
    //Parse checkpoint file line by line
    getline(QMlog,dummy);
    stringstream line(dummy);
    line >> dummy >> dummy;
    if (dummy == "Force")
    {
      line >> dummy;
      if (dummy == "Constants")
      {
        hessDone = 1; //Recovered the Hessian
        //Read force constants (lower triangular)
        for (int i=0;i<Ndof;i++)
        {
          for (int j=0;j<(i+1);j++)
          {
            //Read matrix element in scientific notation
            string matElmt;
            QMlog >> matElmt; //Save as a string
            //Change D notation to E notation
            LICHEMFixSciNot(matElmt);
            //Convert to a double and save matrix element
            QMHess(i,j) = atof(matElmt.c_str());
            //Apply symmetry
            QMHess(j,i) = QMHess(i,j);
          }
        }
      }
    }
  }
  //Check for errors
  if (!hessDone)
  {
    //Calculation did not finish
    cerr << "Error: No force constants recovered!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << bead;
  call << ".com";
  call << " ";
  call << "LICHM_" << bead;
  call << ".fchk";
  globalSys = system(call.str().c_str());
  //Return
  return QMHess;
};

double GaussianOpt(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                   int bead)
{
  //Runs Gaussian optimizations with GauExternal
  fstream inFile,QMlog; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  int extCPUs = 1; //Number of CPUs for GauExternal
  if (AMOEBA and TINKER)
  {
    RotateTINKCharges(QMMMData,bead);
  }
  //Write a new XYZ
  //NB: GauExternal needs different input than the rest of the wrappers
  call.str("");
  call << "LICHMExt_" << bead << ".xyz";
  inFile.open(call.str().c_str(),ios_base::out);
  inFile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Print XYZ coordinates
    inFile << QMMMData[i].QMTyp << " ";
    inFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16) << " ";
    inFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16) << " ";
    inFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16) << '\n';
  }
  inFile.flush();
  inFile.close();
  //Write Gaussian input
  call.str("");
  call << "LICHMExt_" << bead << ".com";
  inFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=LICHMExt_" << bead << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM;
  if (QMMMOpts.memMB)
  {
    call << "MB";
  }
  else
  {
    call << "GB";
  }
  call << '\n';
  call << "%NprocShared=";
  if (Ncpus <= 2)
  {
    call << 1;
  }
  else
  {
    call << 2;
    extCPUs = Ncpus-2;
  }
  call << '\n';
  call << "#P " << "external=\"lichem -GauExtern ";
  call << "LICHMExt_" << bead;  //Just the stub
  call << " -n " << extCPUs;
  call << " -c " << conFilename;
  call << " -r " << regFilename;
  call << " -b " << bead;
  call << "\"" << '\n';
  call << "Symmetry=None Opt=(";
  call << "MaxCycles=" << QMMMOpts.maxOptSteps;
  call << ",MaxStep=" << int(round((QMMMOpts.maxStep/(0.01*bohrRad))));
  call << ")" << '\n';
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.charge << " " << QMMMOpts.spin << '\n';
  //Add atoms
  for (int i=0;i<Natoms;i++)
  {
    if (QMMMData[i].QMRegion)
    {
      call << QMMMData[i].QMTyp;
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
      call << '\n';
    }
    if (QMMMData[i].PBRegion)
    {
      call << "F";
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
      call << " " << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
      call << '\n';
    }
  }
  call << '\n'; //Blank line needed
  //Write Gaussian input
  inFile << call.str();
  inFile.close();
  //Run Optimization
  call.str("");
  call << "g09 ";
  call << "LICHMExt_" << bead;
  globalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHMExt_";
  call << bead << ".log";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool optFinished = 0;
  while (!QMlog.eof())
  {
    //This loop will find the last geometry even if the calculation
    //did not complete
    getline(QMlog,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Center")
    {
      line >> dummy >> dummy;
      line >> dummy;
      if (dummy == "Coordinates")
      {
        //Clear junk
        getline(QMlog,dummy);
        getline(QMlog,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (QMMMData[i].QMRegion or QMMMData[i].PBRegion)
          {
            //Get new coordinates
            getline(QMlog,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy;
            line >> QMMMData[i].P[bead].x;
            line >> QMMMData[i].P[bead].y;
            line >> QMMMData[i].P[bead].z;
          }
        }
      }
    }
    if (dummy == "--")
    {
      line >> dummy;
      if (dummy == "Stationary")
      {
        optFinished = 1;
      }
    }
  }
  QMlog.close();
  //Clean up files
  call.str("");
  call << "rm -f LICHMExt_";
  call << bead << ".*";
  globalSys = system(call.str().c_str());
  //Print warnings and errors
  if (!optFinished)
  {
    cerr << "Warning: Optimization did not converge!!!";
    cerr << '\n';
    cerr << "An older geometry will be recovered...";
    cerr << '\n';
    E = hugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << bead << ".chk";
    globalSys = system(call.str().c_str());
  }
  //Calculate new point-charges and return
  GaussianCharges(QMMMData,QMMMOpts,bead);
  return E;
};

