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
double PSIForces(vector<QMMMAtom>& Struct, VectorXd& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces and charges on a set of atoms
  fstream ifile; //Generic file name
  string dummy; //Generic string
  stringstream call; //Steam for system calls and reading/writing files
  call.copyfmt(cout);
  double E = 0;
  //Set up force calculation
  call.str("");
  call << "Eqm = gradient('" << QMMMOpts.Func << "')" << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  call << "oeprop('MULLIKEN_CHARGES')" << '\n';
  WritePSIInput(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
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
          line >> Fx;
          line >> Fy;
          line >> Fz;
          //Switch to eV/A and save forces
          if (abs(Fx) >= 1e-12)
          {
            Forces(3*i) += Fx*Har2eV/BohrRad;
          }
          if (abs(Fy) >= 1e-12)
          {
            Forces(3*i+1) += Fy*Har2eV/BohrRad;
          }
          if (abs(Fz) >= 1e-12)
          {
            Forces(3*i+2) += Fz*Har2eV/BohrRad;
          }
        }
      }
    }
    line >> dummy; //Get rid of junk
    if (dummy == "Final")
    {
      line >> dummy;
      if (dummy == "Energy:")
      {
        line >> E;
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
      line >> E;
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

void PSICharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Function to update QM point-charges
  fstream ifile;
  string dummy; //Generic string
  stringstream call; //Steam for system calls and reading/writing files
  call.copyfmt(cout);
  //Set up charge calculation
  call.str("");
  call << "energy('" << QMMMOpts.Func << "')" << '\n';
  call << "oeprop('MULLIKEN_CHARGES')" << '\n';
  WritePSIInput(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
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

double PSIEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs PSI4 for energy calculations
  fstream ifile;
  string dummy; //Generic string
  stringstream call; //Steam for system calls and reading/writing files
  call.copyfmt(cout);
  double E = 0.0;
  //Set up energy calculation
  call.str("");
  call << "Eqm = energy('" << QMMMOpts.Func << "')" << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  call << "oeprop('MULLIKEN_CHARGES')" << '\n';
  WritePSIInput(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
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
      line >> dummy;
      if (dummy == "Energy:")
      {
        //Read energy
        line >> E;
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
      line >> E;
      QMfinished = 1;
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
  call << "LICHM_" << Bead << ".dat ";
  call << "LICHM_" << Bead << ".out ";
  call << "LICHM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units
  E *= Har2eV;
  return E;
};

double PSIOpt(vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs psi4 for energy calculations
  fstream ifile;
  string dummy; //Generic string
  stringstream call; //Steam for system calls and reading/writing files
  call.copyfmt(cout);
  double E = 0.0;
  //Set up QM only optimization
  call.str("");
  call << "Eqm = optimize('" << QMMMOpts.Func << "')" << '\n';
  call << "print('Energy: '+`Eqm`)" << '\n';
  call << "oeprop('MULLIKEN_CHARGES')" << '\n';
  WritePSIInput(Struct,call.str(),QMMMOpts,Bead);
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus << "-i ";
  call << "LICHM_" << Bead << ".dat -o ";
  call << "LICHM_" << Bead << ".out > ";
  call << "LICHM_" << Bead << ".log";
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
      line >> dummy;
      if (dummy == "Energy:")
      {
        //Read energy
        line >> E;
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
      line >> E;
      QMfinished = 1;
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
  if (!Optfinished)
  {
    cerr << "Warning: Optimization did not converge!!!";
    cerr << '\n';
    cerr << " LICHEM will attempt to continue using the";
    cerr << " old structure...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
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

