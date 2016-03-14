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
  fstream ifile; //Generic file stream
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
  GlobalSys = system(call.str().c_str());
  //Parse output for energy
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
    //Search for energy
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> E; //Read energy
        QMfinished = 1;
      }
    }
  }
  ifile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Read charges
        stringstream line;
        getline(ifile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> Struct[i].MP[Bead].q;
      }
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
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    GlobalSys = system(call.str().c_str());
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
  GlobalSys = system(call.str().c_str());
  return;
};

double NWChemEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                    int Bead)
{
  //Runs NWChem energy calculations
  fstream ifile; //Generic file stream
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
  GlobalSys = system(call.str().c_str());
  //Parse output for energy
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
    //Search for energy
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> E; //Read energy
        QMfinished = 1;
      }
    }
  }
  ifile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Read charges
        stringstream line;
        getline(ifile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> Struct[i].MP[Bead].q;
      }
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
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    GlobalSys = system(call.str().c_str());
  }
  //Clean up files
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "cp LICHM_";
    call << Bead << ".nw " << QMMMOpts.BackDir << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << "; ";
    call << "cp LICHM_" << Bead << ".movecs ";
    call << QMMMOpts.BackDir << "/.";
    call << " 2> LICHM_" << Bead << ".trash; ";
    call << "rm -f LICHM_" << Bead << ".trash";
    call << "; ";
    call << "cp LICHM_" << Bead << ".log ";
    call << QMMMOpts.BackDir << "/.";
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
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E *= Har2eV;
  return E;
};

double NWChemForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
                    QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem force calculations
  fstream ifile; //Generic file stream
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
  GlobalSys = system(call.str().c_str());
  //Parse output for forces and energies
  call.str("");
  call << "LICHM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  bool QMfinished = 0;
  bool GradDone = 0;
  while (!ifile.eof())
  {
    stringstream line;
    getline(ifile,dummy);
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
        QMfinished = 1;
      }
    }
    if (dummy == "atom")
    {
      line >> dummy >> dummy;
      if (dummy == "gradient")
      {
        GradDone = 1; //Not grad school, that lasts forever
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          //Initialize temporary force variables
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Extract forces; Convoluted, but "easy"
          getline(ifile,dummy);
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
          Forces(3*i) += Fx*Har2eV/BohrRad;
          Forces(3*i+1) += Fy*Har2eV/BohrRad;
          Forces(3*i+2) += Fz*Har2eV/BohrRad;
        }
      }
    }
  }
  ifile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        stringstream line;
        getline(ifile,dummy);
        line.str(dummy);
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        line >> Struct[i].MP[Bead].q;
      }
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
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    GlobalSys = system(call.str().c_str());
  }
  if (!GradDone)
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
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E *= Har2eV;
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
  GlobalSys = system(call.str().c_str());
  //Parse output for Hessian
  call.str("");
  call << "LICHM_" << Bead << ".hess";
  QMlog.open(call.str().c_str(),ios_base::in);
  bool HessDone = 0;
  if (QMlog.good())
  {
    HessDone = 1;
    for (int i=0;i<Ndof;i++)
    {
      for (int j=0;j<(i+1);j++)
      {
        //Read matrix element
        QMlog >> QMHess(i,j);
        //Apply symmetry
        QMHess(j,i) = QMHess(i,j);
      }
    }
  }
  QMlog.close();
  //Check for errors
  if (!HessDone)
  {
    //Calculation did not finish
    cerr << "Error: No force constants recovered!!!";
    cerr << '\n';
    cerr.flush(); //Print warning immediately
    //Delete checkpoint
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    GlobalSys = system(call.str().c_str());
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
  GlobalSys = system(call.str().c_str());
  //Return Hessian
  return QMHess;
};

double NWChemOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem optimizations
  fstream ifile; //Generic file stream
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
  GlobalSys = system(call.str().c_str());
  //Parse output
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
    //Search for energy
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "energy")
      {
        line >> dummy; //Clear junk
        line >> E; //Read energy
        QMfinished = 1;
      }
    }
  }
  ifile.close();
  //Parse output for charges
  call.str("");
  call << "LICHM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        //Read charges
        stringstream line;
        getline(ifile,dummy);
        line.str(dummy);
        //Clear junk
        line >> dummy >> dummy;
        line >> dummy >> dummy;
        //Save charge
        line >> Struct[i].MP[Bead].q;
      }
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
    //Remove checkpoint file
    call.str("");
    call << "rm -f LICHM_" << Bead << ".movecs";
    GlobalSys = system(call.str().c_str());
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
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E *= Har2eV;
  return E;
};

