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

//QM utility functions


//QM wrapper functions
void PSI4Charges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Function to update QM point-charges
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  //Check if there is a checkpoint file
  bool UseCheckPoint;
  call.str("");
  call << "LICHM_" << Bead << ".180";
  UseCheckPoint = CheckFile(call.str());
  //Set up charge calculation
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (UseCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << Bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  if (QMMM)
  {
    call << "oeprop(qmwfn,'MULLIKEN_CHARGES')" << '\n';
  }
  WritePSI4Input(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  globalSys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "LICHM_" << Bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  while (!inFile.eof())
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
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
  }
  inFile.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  return;
};

double PSI4Energy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs PSI4 for energy calculations
  fstream inFile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0.0;
  //Check if there is a checkpoint file
  bool UseCheckPoint;
  call.str("");
  call << "LICHM_" << Bead << ".180";
  UseCheckPoint = CheckFile(call.str());
  //Set up energy calculation
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (UseCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << Bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  if (QMMM)
  {
    call << "oeprop(qmwfn,'MULLIKEN_CHARGES')" << '\n';
  }
  WritePSI4Input(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  globalSys = system(call.str().c_str());
  //Read energy
  call.str("");
  call << "LICHM_" << Bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  while (!inFile.eof())
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
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> Struct[i].MP[Bead].q;
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
  call << "LICHM_" << Bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  while (!inFile.eof())
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
    call << "rm -f *LICHM_" << Bead << ".180";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_" << Bead << ".* ";
    call << QMMMOpts.backDir << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << " "; //Extra blank space before the next command
  }
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Change units
  E *= har2eV;
  return E;
};

double PSI4Forces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                  QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces and charges on a set of atoms
  fstream inFile; //Generic file name
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0;
  //Check if there is a checkpoint file
  bool UseCheckPoint;
  call.str("");
  call << "LICHM_" << Bead << ".180";
  UseCheckPoint = CheckFile(call.str());
  //Set up force calculation
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (UseCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << Bead << ".180']";
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
  WritePSI4Input(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  globalSys = system(call.str().c_str());
  //Extract forces
  call.str("");
  call << "LICHM_" << Bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  while (!inFile.eof())
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
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> Struct[i].MP[Bead].q;
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
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Extract forces; Convoluted, but "easy"
          getline(inFile,dummy);
          stringstream line(dummy);
          line >> dummy; //Clear junk
          line >> Fx; //Read value
          line >> Fy; //Read value
          line >> Fz; //Read value
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
  call << "LICHM_" << Bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  while (!inFile.eof())
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
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Change units
  E *= har2eV;
  return E;
};

MatrixXd PSI4Hessian(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
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
  bool UseCheckPoint;
  call.str("");
  call << "LICHM_" << Bead << ".180";
  UseCheckPoint = CheckFile(call.str());
  //Calculate Hessian
  call.str("");
  call << "Eqm,qmwfn = energy('" << QMMMOpts.func << "'";
  if (UseCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << Bead << ".180']";
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
  WritePSI4Input(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  globalSys = system(call.str().c_str());
  //Extract Hessian
  call.str("");
  call << "LICHM_" << Bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  bool HessDone = 0;
  while (!inFile.eof())
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
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> Struct[i].MP[Bead].q;
          }
        }
      }
    }
    line >> dummy; //Get rid of junk
    if (dummy == "Hessian")
    {
      HessDone = 1;
      getline(inFile,dummy); //Clear junk
      //Read Hessian in groups of Ndofx5
      int rowct = 0; //Current row ID
      while (rowct < Ndof)
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
            if ((rowct+j) < Ndof)
            {
              inFile >> QMHess(i,rowct+j);
            }
          }
        }
        //Update row location
        rowct += 5;
      }
    }
  }
  inFile.close();
  //Check for errors
  if (!HessDone)
  {
    //Calculation did not finish
    cerr << "Error: No force constants recovered!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".180";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  return QMHess;
};

double PSI4Opt(vector<QMMMAtom>& Struct,
               QMMMSettings& QMMMOpts, int Bead)
{
  //Runs PSI4 for energy calculations
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
  WritePSI4Input(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  globalSys = system(call.str().c_str());
  //Read energy and structure
  call.str("");
  call << "LICHM_" << Bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  bool QMFinished = 0;
  bool OptFinished = 0;
  while (!inFile.eof())
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
        OptFinished = 1;
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
          Struct[i].P[Bead].x = x;
          Struct[i].P[Bead].y = y;
          Struct[i].P[Bead].z = z;
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
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy;
            line >> dummy >> dummy;
            line >> dummy;
            line >> Struct[i].MP[Bead].q;
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
  call << "LICHM_" << Bead << ".log";
  inFile.open(call.str().c_str(),ios_base::in);
  while (!inFile.eof())
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
    call << "rm -f *LICHM_" << Bead << ".180";
    globalSys = system(call.str().c_str());
  }
  if (!OptFinished)
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
    call << "rm -f *LICHM_" << Bead << ".180";
    globalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  globalSys = system(call.str().c_str());
  //Change units
  E *= har2eV;
  return E;
};

