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
  fstream ifile; //Generic file stream
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
  call << "Eqm,qmwfn = energy('" << QMMMOpts.Func << "'";
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
  GlobalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  GlobalSys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "LICHM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(ifile,dummy);
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
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  return;
};

double PSI4Energy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs PSI4 for energy calculations
  fstream ifile; //Generic file stream
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
  call << "Eqm,qmwfn = energy('" << QMMMOpts.Func << "'";
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
  GlobalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  GlobalSys = system(call.str().c_str());
  //Read energy
  call.str("");
  call << "LICHM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  bool QMfinished = 0;
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(ifile,dummy);
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
        QMfinished = 1;
      }
    }
  }
  ifile.close();
  //Collect energy (post-SCF)
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Energy:")
    {
      line >> E; //Read post-SCF energy
      QMfinished = 1;
    }
  }
  ifile.close();
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
    call << "rm -f *LICHM_" << Bead << ".180";
    GlobalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_" << Bead << ".* ";
    call << QMMMOpts.BackDir << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << " "; //Extra blank space before the next command
  }
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units
  E *= Har2eV;
  return E;
};

double PSI4Forces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                  QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces and charges on a set of atoms
  fstream ifile; //Generic file name
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
  call << "Eqm,qmwfn = energy('" << QMMMOpts.Func << "'";
  if (UseCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << Bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  call << "gradient('" << QMMMOpts.Func << "'";
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
  GlobalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  GlobalSys = system(call.str().c_str());
  //Extract forces
  call.str("");
  call << "LICHM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(ifile,dummy);
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
        getline(ifile,dummy);
        getline(ifile,dummy);
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Extract forces; Convoluted, but "easy"
          getline(ifile,dummy);
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
          Forces(3*i) += Fx*Har2eV/BohrRad;
          Forces(3*i+1) += Fy*Har2eV/BohrRad;
          Forces(3*i+2) += Fz*Har2eV/BohrRad;
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
  ifile.close();
  //Collect energy (post-SCF)
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy; //Check property
    if (dummy == "Energy:")
    {
      line >> E; //Read post-SCF energy
    }
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units
  E *= Har2eV;
  return E;
};

MatrixXd PSI4Hessian(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs PSI4 to calculate a Hessian
  fstream ifile; //Generic file stream
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
  call << "Eqm,qmwfn = energy('" << QMMMOpts.Func << "'";
  if (UseCheckPoint)
  {
    //Collect old wavefunction from restart file
    call << ",restart_file=[";
    call << "'./LICHM_" << Bead << ".180']";
  }
  call << ",return_wfn=True)" << '\n';
  call << "QMHess = hessian('" << QMMMOpts.Func << "'";
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
  GlobalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  GlobalSys = system(call.str().c_str());
  //Extract Hessian
  call.str("");
  call << "LICHM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  bool HessDone = 0;
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Mulliken")
    {
      line >> dummy;
      if (dummy == "Charges:")
      {
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(ifile,dummy);
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
      getline(ifile,dummy); //Clear junk
      //Read Hessian in groups of Ndofx5
      int rowct = 0; //Current row ID
      while (rowct < Ndof)
      {
        //Clear junk
        getline(ifile,dummy); //Clear junk
        getline(ifile,dummy); //Clear junk
        getline(ifile,dummy); //Clear junk
        //Read elements
        for (int i=0;i<Ndof;i++)
        {
          ifile >> dummy; //Clear junk
          for (int j=0;j<5;j++)
          {
            if ((rowct+j) < Ndof)
            {
              ifile >> QMHess(i,rowct+j);
            }
          }
        }
        //Update row location
        rowct += 5;
      }
    }
  }
  ifile.close();
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
    GlobalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  return QMHess;
};

double PSI4Opt(vector<QMMMAtom>& Struct,
               QMMMSettings& QMMMOpts, int Bead)
{
  //Runs PSI4 for energy calculations
  fstream ifile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0.0;
  //Set up QM only optimization
  call.str("");
  call << "optimize('" << QMMMOpts.Func << "'";
  call << ")" << '\n';
  call << "Eqm,qmwfn = energy('" << QMMMOpts.Func << "'";
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
  GlobalSys = system(call.str().c_str());
  //Save checkpoint file for the next calculation
  call.str("");
  call << "mv *.LICHM_" << Bead << ".180 ";
  call << "LICHM_" << Bead << ".180 ";
  call << "2> LICHM_" << Bead << ".trash; ";
  call << "rm -f LICHM_" << Bead << ".trash";
  GlobalSys = system(call.str().c_str());
  //Read energy and structure
  call.str("");
  call << "LICHM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  bool QMfinished = 0;
  bool Optfinished = 0;
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Final")
    {
      line >> dummy;
      if (dummy == "optimized")
      {
        //Read new geometry
        Optfinished = 1;
        for (int i=0;i<5;i++)
        {
          //Clear junk
          getline(ifile,dummy);
        }
        for (int i=0;i<Natoms;i++)
        {
          getline(ifile,dummy);
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
        getline(ifile,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if (Struct[i].QMregion or Struct[i].PBregion)
          {
            getline(ifile,dummy);
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
        QMfinished = 1;
      }
    }
  }
  ifile.close();
  //Collect energy (post-SCF)
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy; //Check value
    if (dummy == "Energy:")
    {
      line >> E; //Read post-SCF energy
      QMfinished = 1;
    }
  }
  ifile.close();
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
    call << "rm -f *LICHM_" << Bead << ".180";
    GlobalSys = system(call.str().c_str());
  }
  if (!Optfinished)
  {
    cerr << "Warning: Optimization did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue using the";
    cerr << " old structure...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f *LICHM_" << Bead << ".180";
    GlobalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units
  E *= Har2eV;
  return E;
};

