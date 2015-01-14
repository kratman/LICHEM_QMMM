/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for PSI4.

 Citation for PSI4:
 Turney et al. WIREs Comp. Mol. Sci. 2, 4, 556, 2012

*/

//QM utility functions

//QM wrapper functions
double PSIForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces on a set of atoms
  int sys;
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double Eqm = 0;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Set up memory
  call.str("");
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
        call << Struct[i].P[Bead].x/BohrRad << ",";
        call << Struct[i].P[Bead].y/BohrRad << ",";
        call << Struct[i].P[Bead].z/BohrRad;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  if ((QMMM == 1) and (AMOEBA == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q1 << ",";
        call << Struct[i].PC[Bead].x1/BohrRad << ",";
        call << Struct[i].PC[Bead].y1/BohrRad << ",";
        call << Struct[i].PC[Bead].z1/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q2 << ",";
        call << Struct[i].PC[Bead].x2/BohrRad << ",";
        call << Struct[i].PC[Bead].y2/BohrRad << ",";
        call << Struct[i].PC[Bead].z2/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q3 << ",";
        call << Struct[i].PC[Bead].x3/BohrRad << ",";
        call << Struct[i].PC[Bead].y3/BohrRad << ",";
        call << Struct[i].PC[Bead].z3/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q4 << ",";
        call << Struct[i].PC[Bead].x4/BohrRad << ",";
        call << Struct[i].PC[Bead].y4/BohrRad << ",";
        call << Struct[i].PC[Bead].z4/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q5 << ",";
        call << Struct[i].PC[Bead].x5/BohrRad << ",";
        call << Struct[i].PC[Bead].y5/BohrRad << ",";
        call << Struct[i].PC[Bead].z5/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q6 << ",";
        call << Struct[i].PC[Bead].x6/BohrRad << ",";
        call << Struct[i].PC[Bead].y6/BohrRad << ",";
        call << Struct[i].PC[Bead].z6/BohrRad;
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
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Set up memory
  call.str("");
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
        call << Struct[i].P[Bead].x/BohrRad << ",";
        call << Struct[i].P[Bead].y/BohrRad << ",";
        call << Struct[i].P[Bead].z/BohrRad;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  if ((QMMM == 1) and (AMOEBA == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q1 << ",";
        call << Struct[i].PC[Bead].x1/BohrRad << ",";
        call << Struct[i].PC[Bead].y1/BohrRad << ",";
        call << Struct[i].PC[Bead].z1/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q2 << ",";
        call << Struct[i].PC[Bead].x2/BohrRad << ",";
        call << Struct[i].PC[Bead].y2/BohrRad << ",";
        call << Struct[i].PC[Bead].z2/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q3 << ",";
        call << Struct[i].PC[Bead].x3/BohrRad << ",";
        call << Struct[i].PC[Bead].y3/BohrRad << ",";
        call << Struct[i].PC[Bead].z3/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q4 << ",";
        call << Struct[i].PC[Bead].x4/BohrRad << ",";
        call << Struct[i].PC[Bead].y4/BohrRad << ",";
        call << Struct[i].PC[Bead].z4/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q5 << ",";
        call << Struct[i].PC[Bead].x5/BohrRad << ",";
        call << Struct[i].PC[Bead].y5/BohrRad << ",";
        call << Struct[i].PC[Bead].z5/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q6 << ",";
        call << Struct[i].PC[Bead].x6/BohrRad << ",";
        call << Struct[i].PC[Bead].y6/BohrRad << ",";
        call << Struct[i].PC[Bead].z6/BohrRad;
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

double PSIEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs psi4 for energy calculations
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Set up memory
  call.str("");
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
        call << Struct[i].P[Bead].x/BohrRad << ",";
        call << Struct[i].P[Bead].y/BohrRad << ",";
        call << Struct[i].P[Bead].z/BohrRad;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  if ((QMMM == 1) and (AMOEBA == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q1 << ",";
        call << Struct[i].PC[Bead].x1/BohrRad << ",";
        call << Struct[i].PC[Bead].y1/BohrRad << ",";
        call << Struct[i].PC[Bead].z1/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q2 << ",";
        call << Struct[i].PC[Bead].x2/BohrRad << ",";
        call << Struct[i].PC[Bead].y2/BohrRad << ",";
        call << Struct[i].PC[Bead].z2/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q3 << ",";
        call << Struct[i].PC[Bead].x3/BohrRad << ",";
        call << Struct[i].PC[Bead].y3/BohrRad << ",";
        call << Struct[i].PC[Bead].z3/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q4 << ",";
        call << Struct[i].PC[Bead].x4/BohrRad << ",";
        call << Struct[i].PC[Bead].y4/BohrRad << ",";
        call << Struct[i].PC[Bead].z4/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q5 << ",";
        call << Struct[i].PC[Bead].x5/BohrRad << ",";
        call << Struct[i].PC[Bead].y5/BohrRad << ",";
        call << Struct[i].PC[Bead].z5/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q6 << ",";
        call << Struct[i].PC[Bead].x6/BohrRad << ",";
        call << Struct[i].PC[Bead].y6/BohrRad << ",";
        call << Struct[i].PC[Bead].z6/BohrRad;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  //Set up energy calculation
  call << "energy('" << QMMMOpts.Func << "')" << '\n';
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
  //Read energy
  call.str("");
  call << "QMMM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  bool QMfinished = 0;
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    line >> dummy >> dummy;
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
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".dat ";
  call << "QMMM_" << Bead << ".out ";
  call << "QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  if ((AMOEBA == 1) and (QMMM == 1))
  {
    //Update point charges for the polarization energy calculation
    PSICharges(Struct,QMMMOpts,Bead);
  }
  //Change units
  E *= Har2eV;
  return E;
};

double PSIOpt(vector<QMMMAtom>& Struct,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs psi4 for energy calculations
  fstream ofile,ifile;
  string dummy;
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int sys;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    RotateTINKCharges(Struct,Bead);
  }
  //Set up memory
  call.str("");
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
        call << Struct[i].P[Bead].x/BohrRad << ",";
        call << Struct[i].P[Bead].y/BohrRad << ",";
        call << Struct[i].P[Bead].z/BohrRad;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  if ((QMMM == 1) and (AMOEBA == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion == 1)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q1 << ",";
        call << Struct[i].PC[Bead].x1/BohrRad << ",";
        call << Struct[i].PC[Bead].y1/BohrRad << ",";
        call << Struct[i].PC[Bead].z1/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q2 << ",";
        call << Struct[i].PC[Bead].x2/BohrRad << ",";
        call << Struct[i].PC[Bead].y2/BohrRad << ",";
        call << Struct[i].PC[Bead].z2/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q3 << ",";
        call << Struct[i].PC[Bead].x3/BohrRad << ",";
        call << Struct[i].PC[Bead].y3/BohrRad << ",";
        call << Struct[i].PC[Bead].z3/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q4 << ",";
        call << Struct[i].PC[Bead].x4/BohrRad << ",";
        call << Struct[i].PC[Bead].y4/BohrRad << ",";
        call << Struct[i].PC[Bead].z4/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q5 << ",";
        call << Struct[i].PC[Bead].x5/BohrRad << ",";
        call << Struct[i].PC[Bead].y5/BohrRad << ",";
        call << Struct[i].PC[Bead].z5/BohrRad;
        call << ")" << '\n';
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].PC[Bead].q6 << ",";
        call << Struct[i].PC[Bead].x6/BohrRad << ",";
        call << Struct[i].PC[Bead].y6/BohrRad << ",";
        call << Struct[i].PC[Bead].z6/BohrRad;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python(\'EXTERN\',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  //Set up QM only optimization
  call << "optimize('" << QMMMOpts.Func << "')" << '\n';
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
  //Read energy and structure
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
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
  }
  if (Optfinished == 0)
  {
    cout << "Warning: Optimization did not converge!!!";
    cout << '\n';
    cout << " FLUKE will attempt to continue using the";
    cout << " old structure...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".dat ";
  call << "QMMM_" << Bead << ".out ";
  call << "QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  if ((AMOEBA == 1) and (QMMM == 1))
  {
    //Update point charges for the polarization energy calculation
    PSICharges(Struct,QMMMOpts,Bead);
  }
  //Change units
  E *= Har2eV;
  return E;
};

