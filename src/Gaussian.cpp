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
void GaussianCharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                     int Bead)
{
  //Function to update QM point-charges
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  bool UseCheckPoint = CheckFile(call.str());
  if (QMMMOpts.Func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    UseCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.Func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.Func << "/"; //Print the method
  }
  call << QMMMOpts.Basis << " SP Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)" << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.Func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (UseCheckPoint)
    {
      //Read pseudo potential
      call << "Guess=TCheck ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.Func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(Struct,call.str(),QMMMOpts,Bead);
  //Run QM calculation
  call.str("");
  call << "g09 LICHM_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      //Mulliken charges (fallback)
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
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
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
  }
  ifile.close();
  //Clean up files and save checkpoint file
  call.str("");
  call << "mv LICHM_" << Bead;
  call << ".chk tmp_" << Bead;
  call << ".chk; ";
  call << "rm -f ";
  call << "LICHM_" << Bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".com";
  call << "; mv tmp_" << Bead;
  call << ".chk LICHM_" << Bead;
  call << ".chk";
  GlobalSys = system(call.str().c_str());
  return;
};

double GaussianEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                      int Bead)
{
  //Calculates the QM energy with Gaussian
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  double Eself = 0.0; //External field self-energy
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  bool UseCheckPoint = CheckFile(call.str());
  if (QMMMOpts.Func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    UseCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.Func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.Func << "/"; //Print the method
  }
  call << QMMMOpts.Basis << " SP Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)";
  if (UseCheckPoint)
  {
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.Func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.Func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(Struct,call.str(),QMMMOpts,Bead);
  //Calculate energy
  call.str("");
  call << "g09 ";
  call << "LICHM_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Read output
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  bool QMfinished = 0;
  while (!ifile.eof())
  {
    stringstream line;
    getline(ifile,dummy);
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
        QMfinished = 1;
      }
    }
    //Check for charges
    if (dummy == "Mulliken")
    {
      //Mulliken charges (fallback)
      line >> dummy;
      if (dummy == "charges:")
      {
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
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
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
  }
  //Check for errors
  if (!QMfinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  ifile.close();
  //Clean up files and save checkpoint file
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_";
    call << Bead << ".* ";
    call << QMMMOpts.BackDir;
    call << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << " "; //Extra blank space before the next command
  }
  call << "rm -f ";
  call << "LICHM_" << Bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".com";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E -= Eself;
  E *= Har2eV;
  return E;
};

double GaussianForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                      QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces on a set of atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream ofile,ifile,QMlog; //Generic input files
  double Eqm = 0; //QM energy
  double Eself = 0; //External field self-energy
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  bool UseCheckPoint = CheckFile(call.str());
  if (QMMMOpts.Func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    UseCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.Func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.Func << "/"; //Print the method
  }
  call << QMMMOpts.Basis << " Force=NoStep Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)"; //Line ended further below
  if (UseCheckPoint)
  {
    //Restart if possible
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.Func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.Func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(Struct,call.str(),QMMMOpts,Bead);
  //Run Gaussian
  call.str("");
  call << "g09 " << "LICHM_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Extract forces
  call.str("");
  call << "LICHM_" << Bead << ".log";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool GradDone = 0;
  while ((!QMlog.eof()) and (!GradDone))
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
        GradDone = 1; //Not grad school, that lasts forever
        getline(QMlog,dummy); //Clear junk
        getline(QMlog,dummy); //Clear more junk
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Extract forces; Convoluted, but "easy"
          getline(QMlog,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> Fx;
          line >> Fy;
          line >> Fz;
          //Save forces
          Forces(3*i) += Fx*Har2eV/BohrRad;
          Forces(3*i+1) += Fy*Har2eV/BohrRad;
          Forces(3*i+2) += Fz*Har2eV/BohrRad;
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
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
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
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Count through all atoms in the QM calculations
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
  }
  QMlog.close();
  //Check for errors
  if (!GradDone)
  {
    cerr << "Warning: No forces recovered!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to recover...";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".com";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  Eqm -= Eself;
  Eqm *= Har2eV;
  return Eqm;
};

MatrixXd GaussianHessian(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                         int Bead)
{
  //Function for calculating the Hessian for a set of QM atoms
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  string dummy; //Generic string
  fstream ofile,ifile,QMlog; //Generic input files
  int Ndof = 3*(Nqm+Npseudo);
  MatrixXd QMHess(Ndof,Ndof);
  QMHess.setZero();
  //Check if there is a checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  bool UseCheckPoint = CheckFile(call.str());
  if (QMMMOpts.Func == "SemiEmp")
  {
    //Disable checkpoints for the SemiEmp force calculations
    UseCheckPoint = 0;
    //Remove SemiEmp checkpoints to avoid errors
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Construct Gaussian input
  call.str("");
  call << "#P ";
  if (QMMMOpts.Func != "SemiEmp")
  {
    //Avoids defining both a basis set and method for semi-empirical
    call << QMMMOpts.Func << "/"; //Print the method
  }
  call << QMMMOpts.Basis << " Freq Symmetry=None" << '\n';
  call << "Int=UltraFine SCF=(YQC,Big,Direct)"; //Line ended further below
  if (UseCheckPoint)
  {
    //Restart if possible
    call << " Guess=TCheck";
  }
  call << '\n';
  if (QMMM)
  {
    if ((Npseudo > 0) and (QMMMOpts.Func != "SemiEmp"))
    {
      //Read pseudo potential
      call << "Pseudo=Read ";
    }
    if (Nmm > 0)
    {
      //Read charges
      call << "Charge=angstroms ";
    }
    if (QMMMOpts.Func != "SemiEmp")
    {
      //Avoids calculating ESP charges for semi-empirical
      call << "Population=(MK,ReadRadii)";
    }
    call << '\n';
  }
  WriteGauInput(Struct,call.str(),QMMMOpts,Bead);
  //Run Gaussian
  call.str("");
  call << "g09 " << "LICHM_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Generate formatted checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".chk";
  if (CheckFile(call.str()))
  {
    //Run formchk
    call.str("");
    call << "formchk ";
    call << "LICHM_" << Bead << ".chk";
    call << " > LICHM_" << Bead << ".trash; ";
    call << " rm -f LICHM_" << Bead << ".trash";
    GlobalSys = system(call.str().c_str());
  }
  else
  {
    //Calculation did not finish
    cerr << "Error: No checkpoint file found!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Open checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".fchk";
  QMlog.open(call.str().c_str(),ios_base::in);
  //Extract Hessian
  bool HessDone = 0;
  while (QMlog.good() and (!QMlog.eof()) and (!HessDone))
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
        HessDone = 1; //Recovered the Hessian
        //Read force constants (lower triangualr)
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
    }
  }
  //Check for errors
  if (!HessDone)
  {
    //Calculation did not finish
    cerr << "Error: No force constants recovered!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead;
  call << ".log";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".com";
  call << " ";
  call << "LICHM_" << Bead;
  call << ".fchk";
  GlobalSys = system(call.str().c_str());
  //Return
  return QMHess;
};

double GaussianOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                   int Bead)
{
  //Runs Gaussian optimizations with GauExternal
  fstream ofile,ifile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  int ExtCPUs = 1; //Number of CPUs for GauExternal
  if (AMOEBA and TINKER)
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Write a new XYZ
  //NB: GauExternal needs different input than the rest of the wrappers
  call.str("");
  call << "LICHMExt_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Print XYZ coordinates
    ofile << Struct[i].QMTyp << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].x,16) << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].y,16) << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].z,16) << '\n';
  }
  ofile.flush();
  ofile.close();
  //Write Gaussian input
  call.str("");
  call << "LICHMExt_" << Bead << ".com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=LICHMExt_" << Bead << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM;
  if (QMMMOpts.MemMB)
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
    ExtCPUs = Ncpus-2;
  }
  call << '\n';
  call << "#P " << "external=\"lichem -GauExtern ";
  call << "LICHMExt_" << Bead;  //Just the stub
  call << " -n " << ExtCPUs;
  call << " -c " << confilename;
  call << " -r " << regfilename;
  call << " -b " << Bead;
  call << "\"" << '\n';
  call << "Symmetry=None Opt=(";
  call << "MaxCycles=" << QMMMOpts.MaxOptSteps;
  call << ",MaxStep=" << int(round((QMMMOpts.MaxStep/(0.01*BohrRad))));
  call << ")" << '\n';
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.Charge << " " << QMMMOpts.Spin << '\n';
  //Add atoms
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion)
    {
      call << Struct[i].QMTyp;
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].z,16);
      call << '\n';
    }
    if (Struct[i].PBregion)
    {
      call << "F";
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].z,16);
      call << '\n';
    }
  }
  call << '\n'; //Blank line needed
  //Write Gaussian input
  ofile << call.str();
  ofile.close();
  //Run Optimization
  call.str("");
  call << "g09 ";
  call << "LICHMExt_" << Bead;
  GlobalSys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "LICHMExt_";
  call << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  bool Optfinished = 0;
  while (!ifile.eof())
  {
    //This loop will find the last geometry even if the calculation
    //did not complete
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Center")
    {
      line >> dummy >> dummy;
      line >> dummy;
      if (dummy == "Coordinates")
      {
        //Clear junk
        getline(ifile,dummy);
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            //Get new coordinates
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy;
            line >> Struct[i].P[Bead].x;
            line >> Struct[i].P[Bead].y;
            line >> Struct[i].P[Bead].z;
          }
        }
      }
    }
    if (dummy == "--")
    {
      line >> dummy;
      if (dummy == "Stationary")
      {
        Optfinished = 1;
      }
    }
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f LICHMExt_";
  call << Bead << ".*";
  GlobalSys = system(call.str().c_str());
  //Print warnings and errors
  if (!Optfinished)
  {
    cerr << "Warning: Optimization did not converge!!!";
    cerr << '\n';
    cerr << "An older geometry will be recovered...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".chk";
    GlobalSys = system(call.str().c_str());
  }
  //Calculate new point-charges and return
  GaussianCharges(Struct,QMMMOpts,Bead);
  return E;
};

