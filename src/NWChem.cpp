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

//MM utility functions


//MM wrapper functions
double NWChemForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem force calculations
  fstream ifile; //Generic file stream
  string dummy; //Genric string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0;
  //Set up force calculations
  call.str("");
  call << "task dft gradient" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Calculate energy
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
        line >> E;
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
          //Change units and check precision
          if (abs(Fx) >= 1e-6)
          {
            Forces(3*i) -= Fx*Har2eV/BohrRad;
          }
          if (abs(Fy) >= 1e-6)
          {
            Forces(3*i+1) -= Fy*Har2eV/BohrRad;
          }
          if (abs(Fz) >= 1e-6)
          {
            Forces(3*i+2) -= Fz*Har2eV/BohrRad;
          }
        }
      }
    }
  }
  ifile.close();
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
  if (!QMfinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
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
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".m*" << " ";
  call << "LICHM_" << Bead << ".z*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".nw" << " ";
  call << "LICHM_" << Bead << ".db" << " ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E *= Har2eV;
  return E;
};

void NWChemCharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
     int Bead)
{
  //Calculates atomic charges with NWChem
  fstream ifile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy print settings
  double E = 0.0; //QM energy
  //Remove multipoles file
  call.str("");
  call << "rm -f MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Write NWChem input
  call.str("");
  call << "task dft energy" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Calculate energy
  call.str("");
  if (Ncpus > 1)
  {
    //Run in parallel
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
        line >> E;
        QMfinished = 1;
      }
    }
  }
  ifile.close();
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
  if (!QMfinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  //Clean up files and return
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".m*" << " ";
  call << "LICHM_" << Bead << ".z*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".nw" << " ";
  call << "LICHM_" << Bead << ".db" << " ";
  call << "LICHM_" << Bead << ".log";
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
  //Remove multipoles file
  call.str("");
  call << "rm -f MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Write NWChem input
  call.str("");
  call << "task dft energy" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Calculate energy
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
        line >> E;
        QMfinished = 1;
      }
    }
  }
  ifile.close();
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
  if (!QMfinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  //Clean up files
  call.str("");
  if (CheckFile("BACKUPQM"))
  {
    //Save old files
    call << "rm -rf ";
    call << QMMMOpts.BackDir;
    call << " && mkdir ";
    call << QMMMOpts.BackDir;
    call << " && ";
    call << "cp LICHM_";
    call << Bead << ".nw ";
    call << QMMMOpts.BackDir;
    call << "/. && ";
    call << "cp LICHM_";
    call << Bead << ".db ";
    call << QMMMOpts.BackDir;
    call << "/. && ";
    call << "cp LICHM_";
    call << Bead << ".log ";
    call << QMMMOpts.BackDir;
    call << "/. && ";
  }
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".m*" << " ";
  call << "LICHM_" << Bead << ".z*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".nw" << " ";
  call << "LICHM_" << Bead << ".db" << " ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E *= Har2eV;
  return E;
};

double NWChemOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem optimizations
  fstream ifile; //Generic file stream
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Save print settings
  double E = 0.0; //QM energy
  //Remove multipoles file
  call.str("");
  call << "rm -f MMCharges_" << Bead << ".txt";
  GlobalSys = system(call.str().c_str());
  //Write NWChem input
  call.str("");
  call << "task dft optimize" << '\n';
  call << "task esp" << '\n';
  WriteNWChemInput(Struct,call.str(),QMMMOpts,Bead);
  //Calculate energy
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
        line >> E;
        QMfinished = 1;
      }
    }
  }
  ifile.close();
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
  if (!QMfinished)
  {
    cerr << "Warning: SCF did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "LICHM_" << Bead << ".b*" << " ";
  call << "LICHM_" << Bead << ".c*" << " ";
  call << "LICHM_" << Bead << ".g*" << " ";
  call << "LICHM_" << Bead << ".m*" << " ";
  call << "LICHM_" << Bead << ".z*" << " ";
  call << "LICHM_" << Bead << ".p*" << " ";
  call << "LICHM_" << Bead << ".q*" << " ";
  call << "LICHM_" << Bead << ".nw" << " ";
  call << "LICHM_" << Bead << ".db" << " ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E *= Har2eV;
  return E;
};

