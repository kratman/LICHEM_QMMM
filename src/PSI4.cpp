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
double PsiForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the forces on a set of atoms
  double Eqm;
  
  return Eqm;
};

//QM wrapper functions
double PsiWrapper(string RunTyp, vector<QMMMAtom>& Struct,
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
      if (QMMM == 1)
      {
        Struct[i].q = 0; //Hack to make the QMMM wrapper function
      }
      call << "    " << Struct[i].QMTyp;
      if (Bead == -1)
      {
        call << "    " << Struct[i].x;
        call << "    " << Struct[i].y;
        call << "    " << Struct[i].z;
        call << '\n';
      }
      if (Bead != -1)
      {
        call << "    " << Struct[i].P[Bead].x;
        call << "    " << Struct[i].P[Bead].y;
        call << "    " << Struct[i].P[Bead].z;
        call << '\n';
      }
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
        call << Struct[i].q << ",";
        if (Bead == -1)
        {
          call << Struct[i].x << ",";
          call << Struct[i].y << ",";
          call << Struct[i].z;
        }
        if (Bead != -1)
        {
          call << Struct[i].P[Bead].x << ",";
          call << Struct[i].P[Bead].y << ",";
          call << Struct[i].P[Bead].z;
        }
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
  //Print file and call psi4
  if (Bead == -1)
  {
    ofile.open("QMMM.dat",ios_base::out);
    ofile << call.str() << endl;
    ofile.close();
  }
  if (Bead != -1)
  {
    dummy = call.str();
    call.str("");
    call << "QMMM_" << Bead << ".dat";
    ofile.open(call.str().c_str(),ios_base::out);
    ofile << dummy << endl;
    ofile.close();
  }
  call.str("");
  call << "psi4 -n " << Ncpus;
  call << "-u QMMM";
  if (Bead == -1)
  {
    call << ".dat -o QMMM.out > QMMM.log";
  }
  if (Bead != -1)
  {
    call << "_" << Bead;
    call << ".dat -o QMMM";
    call << "_" << Bead;
    call << ".out > QMMM";
    call << "_" << Bead;
    call << ".log";
  }
  sys = system(call.str().c_str());
  //Read energy and/or structure
  if (Bead == -1)
  {
    ifile.open("QMMM.out",ios_base::in);
  }
  if (Bead != -1)
  {
    call.str("");
    call << "QMMM_" << Bead << ".out";
    ifile.open(call.str().c_str(),ios_base::in);
  }
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
          if (Bead == -1)
          {
            Struct[i].x = x;
            Struct[i].y = y;
            Struct[i].z = z;
          }
          if (Bead != -1)
          {
            Struct[i].P[Bead].x = x;
            Struct[i].P[Bead].y = y;
            Struct[i].P[Bead].z = z;
          }
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
  if (Bead == -1)
  {
    call << "QMMM.dat QMMM.out QMMM.log";
  }
  if (Bead != -1)
  {
    call << "QMMM_" << Bead << ".dat ";
    call << "QMMM_" << Bead << ".out ";
    call << "QMMM_" << Bead << ".log";
  }
  sys = system(call.str().c_str());
  //Change units
  E *= Har2eV;
  return E;
};

