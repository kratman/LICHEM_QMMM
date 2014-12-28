/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for PSI4.

*/

//QM utility functions
double PSIForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces on a set of atoms
  int sys;
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  call.str("");
  double Eqm = 0;
  //Set up memory
  call << "memory " << QMMMOpts.RAM;
  call << " gb" << '\n' << '\n';
  //Set up globals
  call << "set globals {" << '\n';
  call << "  basis ";
  call << QMMMOpts.Basis << '\n';
  call << "  guess sad" << '\n';
  call << "  ints_tolerance 1.0E-10" << '\n';
  call << "  scf_type df" << '\n';
  call << "}" << '\n' << '\n';
  //Set up molecules
  call << "molecule QMregion {" << '\n';
  call << "    " << QMMMOpts.Charge;
  call << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion == 1)
    {
      call << "    " << Struct[i].QMTyp;
      call << "    " << Struct[i].P[Bead].x;
      call << "    " << Struct[i].P[Bead].y;
      call << "    " << Struct[i].P[Bead].z;
      call << '\n';
    }
  }
  call << "    symmetry c1" << '\n';
  call << "    no_reorient" << '\n';
  call << "    no_com" << '\n';
  call << "}" << '\n' << '\n';
  //Set up MM field
  if ((QMMM == 1) and (CHRG == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].MP[Bead].q << ",";
        call << Struct[i].P[Bead].x << ",";
        call << Struct[i].P[Bead].y << ",";
        call << Struct[i].P[Bead].z;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  //Set up charge calculation
  call << "energy('" << QMMMOpts.Func << "')" << '\n';
  call << "gradient('" << QMMMOpts.Func << "')" << '\n';
  //Print file
  dummy = call.str(); //Store file as a temporary variable
  call.str("");
  call << "QMMM_" << Bead << ".dat";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << dummy << endl;
  ofile.close();
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus;
  call << "-i QMMM";
  call << "_" << Bead;
  call << ".dat -o QMMM";
  call << "_" << Bead;
  call << ".out > QMMM";
  call << "_" << Bead;
  call << ".log";
  sys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "QMMM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy;
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
          //Save forces
          Forces[i].x += Fx*Har2eV/BohrRad;
          Forces[i].y += Fy*Har2eV/BohrRad;
          Forces[i].z += Fz*Har2eV/BohrRad;
        }
      }
    }
    line >> dummy;
    if (dummy == "Final")
    {
      line >> dummy;
      if (dummy == "Energy:")
      {
        line >> Eqm;
      }
    }
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".dat ";
  call << "QMMM_" << Bead << ".out ";
  call << "QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  //Change units
  Eqm *= Har2eV;
  return Eqm;
};

void PSICharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Function to update QM point charges
  int sys;
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  call.str("");
  //Set up memory
  call << "memory " << QMMMOpts.RAM;
  call << " gb" << '\n' << '\n';
  //Set up globals
  call << "set globals {" << '\n';
  call << "  basis ";
  call << QMMMOpts.Basis << '\n';
  call << "  guess sad" << '\n';
  call << "  ints_tolerance 1.0E-10" << '\n';
  call << "  scf_type df" << '\n';
  call << "}" << '\n' << '\n';
  //Set up molecules
  call << "molecule QMregion {" << '\n';
  call << "    " << QMMMOpts.Charge;
  call << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion == 1)
    {
      call << "    " << Struct[i].QMTyp;
      call << "    " << Struct[i].P[Bead].x;
      call << "    " << Struct[i].P[Bead].y;
      call << "    " << Struct[i].P[Bead].z;
      call << '\n';
    }
  }
  call << "    symmetry c1" << '\n';
  call << "    no_reorient" << '\n';
  call << "    no_com" << '\n';
  call << "}" << '\n' << '\n';
  //Set up MM field
  if ((QMMM == 1) and (CHRG == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].MP[Bead].q << ",";
        call << Struct[i].P[Bead].x << ",";
        call << Struct[i].P[Bead].y << ",";
        call << Struct[i].P[Bead].z;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  //Set up charge calculation
  call << "energy('" << QMMMOpts.Func << "')" << '\n';
  call << "oeprop('MULLIKEN_CHARGES')" << '\n';
  //Print file
  dummy = call.str(); //Store as a temporary variable
  call.str("");
  call << "QMMM_" << Bead << ".dat";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << dummy << endl;
  ofile.close();
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus;
  call << "-i QMMM";
  call << "_" << Bead;
  call << ".dat -o QMMM";
  call << "_" << Bead;
  call << ".out > QMMM";
  call << "_" << Bead;
  call << ".log";
  sys = system(call.str().c_str());
  //Extract charges
  call.str("");
  call << "QMMM_" << Bead << ".out";
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
          if ((Struct[i].QMregion == 1) or (Struct[i].PAregion ==1))
          {
            getline(ifile,dummy);
            if (Struct[i].QMregion == 1)
            {
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
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".dat ";
  call << "QMMM_" << Bead << ".out ";
  call << "QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  return;
};

//QM wrapper functions
double PSIWrapper(string RunTyp, vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs psi4 for energy calculations
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys;
  call.str("");
  //Set up memory
  call << "memory " << QMMMOpts.RAM;
  call << " gb" << '\n' << '\n';
  //Set up globals
  call << "set globals {" << '\n';
  call << "  basis ";
  call << QMMMOpts.Basis << '\n';
  call << "  guess sad" << '\n';
  call << "  ints_tolerance 1.0E-10" << '\n';
  call << "  scf_type df" << '\n';
  call << "}" << '\n' << '\n';
  //Set up molecules
  call << "molecule QMregion {" << '\n';
  call << "    " << QMMMOpts.Charge;
  call << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion == 1)
    {
      call << "    " << Struct[i].QMTyp;
      call << "    " << Struct[i].P[Bead].x;
      call << "    " << Struct[i].P[Bead].y;
      call << "    " << Struct[i].P[Bead].z;
      call << '\n';
    }
  }
  call << "    symmetry c1" << '\n';
  call << "    no_reorient" << '\n';
  call << "    no_com" << '\n';
  call << "}" << '\n' << '\n';
  //Set up MM field
  if ((QMMM == 1) and (CHRG == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].MP[Bead].q << ",";
        call << Struct[i].P[Bead].x << ",";
        call << Struct[i].P[Bead].y << ",";
        call << Struct[i].P[Bead].z;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  //Set up energy calculation
  if (RunTyp == "Enrg")
  {
    call << "energy('" << QMMMOpts.Func << "')" << '\n';
  }
  //Set up QM only optimization
  if (RunTyp == "Opt")
  {
    call << "optimize('" << QMMMOpts.Func << "')" << '\n';
  }
  //Print file
  dummy = call.str(); //Store as a temporary variable
  call.str("");
  call << "QMMM_" << Bead << ".dat";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << dummy << endl;
  ofile.close();
  //Call PSI4
  call.str("");
  call << "psi4 -n " << Ncpus;
  call << "-u QMMM";
  call << "_" << Bead;
  call << ".dat -o QMMM";
  call << "_" << Bead;
  call << ".out > QMMM";
  call << "_" << Bead;
  call << ".log";
  sys = system(call.str().c_str());
  //Read energy and/or structure
  call.str("");
  call << "QMMM_" << Bead << ".out";
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
      if ((dummy == "optimized") and (RunTyp == "Opt"))
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
    line >> dummy;
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
  if (QMfinished == 0)
  {
    cout << "Warning: SCF did not converge!!!";
    cout << '\n';
    cout << " The calculation attempt will continue...";
    cout << '\n';
    E = 10000.0; //Large number to reject step
  }
  if ((Optfinished == 0) and (RunTyp == "Opt"))
  {
    cout << "Warning: Optimization did not converge!!!";
    cout << '\n';
    cout << " The calculation attempt will continue with the";
    cout << " old structure...";
    cout << '\n';
    E = 10000.0; //Large number to reject step
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".dat ";
  call << "QMMM_" << Bead << ".out ";
  call << "QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  //Change units
  E *= Har2eV;
  return E;
};

