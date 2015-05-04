/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for NWChem.

 Reference for NWChem:
 Valiev et al. Comput. Phys. Commun. 181, 1477, (2010)

*/

//MM utility functions


//MM wrapper functions
double NWChemForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem force calculations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    //Set up multipoles
    RotateTINKCharges(Struct,Bead);
  }
  //Create NWChem input
  call.str("");
  call << "QMMM_" << Bead << ".nw";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << "start QMMM_" << Bead << '\n';
  ofile << "charge " << QMMMOpts.Charge << '\n';
  ofile << "geometry nocenter ";
  ofile << "noautoz noautosym" << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion or Struct[i].PAregion)
    {
      ofile << " " << Struct[i].QMTyp;
      ofile << " " << Struct[i].P[Bead].x;
      ofile << " " << Struct[i].P[Bead].y;
      ofile << " " << Struct[i].P[Bead].z;
      ofile << '\n';
    }
  }
  ofile << "end" << '\n';
  if (CheckFile("BASIS"))
  {
    //Add basis set and ecp info
    ifile.open("BASIS",ios_base::in);
    if (ifile.good());
    {
      while (!ifile.eof())
      {
        //Copy BASIS line by line, if BASIS exists
        getline(ifile,dummy);
        ofile << dummy << '\n';
      }
    }
    ifile.close();
  }
  else
  {
    ofile << "basis" << '\n';
    ofile << " * library " << QMMMOpts.Basis;
    ofile << '\n';
    ofile << "end" << '\n';
  }
  if (CHRG == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].P[Bead].x;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].y;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].z;
        ofile << " " << setprecision(12) << Struct[i].MP[Bead].q;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  if (AMOEBA == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  ofile << "dft" << '\n';
  ofile << " mult " << QMMMOpts.Spin << '\n';
  ofile << " direct" << '\n';
  ofile << " grid xfine" << '\n';
  ofile << " noio" << '\n';
  ofile << " tolerances tight" << '\n';
  ofile << " xc " << QMMMOpts.Func << '\n';
  ofile << "end" << '\n';
  ofile << "task dft gradient" << '\n';
  ofile << "task esp" << '\n';
  ofile.flush();
  ofile.close();
  //Calculate energy
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem QMMM_" << Bead << ".nw";
  call << " > QMMM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Parse output
  call.str("");
  call << "QMMM_" << Bead << ".log";
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
    if (dummy == "atom")
    {
      line >> dummy >> dummy;
      if (dummy == "gradient")
      {
        getline(ifile,dummy); //Clear junk
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
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
          if (abs(Fx) >= 1e-8)
          {
            Forces[i].x -= Fx*Har2eV/BohrRad;
          }
          if (abs(Fy) >= 1e-8)
          {
            Forces[i].y -= Fy*Har2eV/BohrRad;
          }
          if (abs(Fz) >= 1e-8)
          {
            Forces[i].z -= Fz*Har2eV/BohrRad;
          }
        }
      }
    }
  }
  ifile.close();
  call.str("");
  call << "QMMM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PAregion)
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
    cerr << " FLUKE will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  //Change units
  E *= Har2eV;
  return E;
};

void NWChemCharges(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Calculates atomic charges with NWChem
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    //Set up multipoles
    RotateTINKCharges(Struct,Bead);
  }
  //Create NWChem input
  call.str("");
  call << "QMMM_" << Bead << ".nw";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << "start QMMM_" << Bead << '\n';
  ofile << "charge " << QMMMOpts.Charge << '\n';
  ofile << "geometry nocenter ";
  ofile << "noautoz noautosym" << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion or Struct[i].PAregion)
    {
      ofile << " " << Struct[i].QMTyp;
      ofile << " " << Struct[i].P[Bead].x;
      ofile << " " << Struct[i].P[Bead].y;
      ofile << " " << Struct[i].P[Bead].z;
      ofile << '\n';
    }
  }
  ofile << "end" << '\n';
  if (CheckFile("BASIS"))
  {
    //Add basis set and ecp info
    ifile.open("BASIS",ios_base::in);
    if (ifile.good());
    {
      while (!ifile.eof())
      {
        //Copy BASIS line by line, if BASIS exists
        getline(ifile,dummy);
        ofile << dummy << '\n';
      }
    }
    ifile.close();
  }
  else
  {
    ofile << "basis" << '\n';
    ofile << " * library " << QMMMOpts.Basis;
    ofile << '\n';
    ofile << "end" << '\n';
  }
  if (CHRG == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].P[Bead].x;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].y;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].z;
        ofile << " " << setprecision(12) << Struct[i].MP[Bead].q;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  if (AMOEBA == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  ofile << "dft" << '\n';
  ofile << " mult " << QMMMOpts.Spin << '\n';
  ofile << " direct" << '\n';
  ofile << " grid xfine" << '\n';
  ofile << " noio" << '\n';
  ofile << " tolerances tight" << '\n';
  ofile << " xc " << QMMMOpts.Func << '\n';
  ofile << "end" << '\n';
  ofile << "task dft energy" << '\n';
  ofile << "task esp" << '\n';
  ofile.flush();
  ofile.close();
  //Calculate energy
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem QMMM_" << Bead << ".nw";
  call << " > QMMM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Parse output
  call.str("");
  call << "QMMM_" << Bead << ".log";
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
  call << "QMMM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PAregion)
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
    cerr << " FLUKE will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  //Clean up files and return
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".b*" << " ";
  call << "QMMM_" << Bead << ".c*" << " ";
  call << "QMMM_" << Bead << ".grid*" << " ";
  call << "QMMM_" << Bead << ".mo*" << " ";
  call << "QMMM_" << Bead << ".zmat" << " ";
  call << "QMMM_" << Bead << ".p*" << " ";
  call << "QMMM_" << Bead << ".db" << " ";
  call << "QMMM_" << Bead << ".q*" << " ";
  call << "QMMM_" << Bead << ".nw" << " ";
  call << "QMMM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  return;
};

double NWChemEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Runs NWChem energy calculations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    //Set up multipoles
    RotateTINKCharges(Struct,Bead);
  }
  //Create NWChem input
  call.str("");
  call << "QMMM_" << Bead << ".nw";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << "start QMMM_" << Bead << '\n';
  ofile << "charge " << QMMMOpts.Charge << '\n';
  ofile << "geometry nocenter ";
  ofile << "noautoz noautosym" << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion or Struct[i].PAregion)
    {
      ofile << " " << Struct[i].QMTyp;
      ofile << " " << Struct[i].P[Bead].x;
      ofile << " " << Struct[i].P[Bead].y;
      ofile << " " << Struct[i].P[Bead].z;
      ofile << '\n';
    }
  }
  ofile << "end" << '\n';
  if (CheckFile("BASIS"))
  {
    //Add basis set and ecp info
    ifile.open("BASIS",ios_base::in);
    if (ifile.good());
    {
      while (!ifile.eof())
      {
        //Copy BASIS line by line, if BASIS exists
        getline(ifile,dummy);
        ofile << dummy << '\n';
      }
    }
    ifile.close();
  }
  else
  {
    ofile << "basis" << '\n';
    ofile << " * library " << QMMMOpts.Basis;
    ofile << '\n';
    ofile << "end" << '\n';
  }
  if (CHRG == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].P[Bead].x;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].y;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].z;
        ofile << " " << setprecision(12) << Struct[i].MP[Bead].q;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  if (AMOEBA == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  ofile << "dft" << '\n';
  ofile << " mult " << QMMMOpts.Spin << '\n';
  ofile << " direct" << '\n';
  ofile << " grid xfine" << '\n';
  ofile << " noio" << '\n';
  ofile << " tolerances tight" << '\n';
  ofile << " xc " << QMMMOpts.Func << '\n';
  ofile << "end" << '\n';
  ofile << "task dft energy" << '\n';
  ofile << "task esp" << '\n';
  ofile.flush();
  ofile.close();
  //Calculate energy
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem QMMM_" << Bead << ".nw";
  call << " > QMMM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Parse output
  call.str("");
  call << "QMMM_" << Bead << ".log";
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
  call << "QMMM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PAregion)
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
    cerr << " FLUKE will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  //Clean up files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".b*" << " ";
  call << "QMMM_" << Bead << ".c*" << " ";
  call << "QMMM_" << Bead << ".grid*" << " ";
  call << "QMMM_" << Bead << ".mo*" << " ";
  call << "QMMM_" << Bead << ".zmat" << " ";
  call << "QMMM_" << Bead << ".p*" << " ";
  call << "QMMM_" << Bead << ".db" << " ";
  call << "QMMM_" << Bead << ".q*" << " ";
  call << "QMMM_" << Bead << ".nw" << " ";
  call << "QMMM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Change units and return
  E *= Har2eV;
  return E;
};

double NWChemOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs NWChem optimizations
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    //Set up multipoles
    RotateTINKCharges(Struct,Bead);
  }
  //Create NWChem input
  call.str("");
  call << "QMMM_" << Bead << ".nw";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << "start QMMM_" << Bead << '\n';
  ofile << "charge " << QMMMOpts.Charge << '\n';
  ofile << "geometry nocenter ";
  ofile << "noautoz noautosym" << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion or Struct[i].PAregion)
    {
      ofile << " " << Struct[i].QMTyp;
      ofile << " " << Struct[i].P[Bead].x;
      ofile << " " << Struct[i].P[Bead].y;
      ofile << " " << Struct[i].P[Bead].z;
      ofile << '\n';
    }
  }
  ofile << "end" << '\n';
  if (CheckFile("BASIS"))
  {
    //Add basis set and ecp info
    ifile.open("BASIS",ios_base::in);
    if (ifile.good());
    {
      while (!ifile.eof())
      {
        //Copy BASIS line by line, if BASIS exists
        getline(ifile,dummy);
        ofile << dummy << '\n';
      }
    }
    ifile.close();
  }
  else
  {
    ofile << "basis" << '\n';
    ofile << " * library " << QMMMOpts.Basis;
    ofile << '\n';
    ofile << "end" << '\n';
  }
  if (CHRG == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].P[Bead].x;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].y;
        ofile << " " << setprecision(12) << Struct[i].P[Bead].z;
        ofile << " " << setprecision(12) << Struct[i].MP[Bead].q;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  if (AMOEBA == 1)
  {
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        ofile << '\n';
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        ofile.copyfmt(cout);
        ofile << '\n';
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  ofile << "dft" << '\n';
  ofile << " mult " << QMMMOpts.Spin << '\n';
  ofile << " direct" << '\n';
  ofile << " grid xfine" << '\n';
  ofile << " noio" << '\n';
  ofile << " tolerances tight" << '\n';
  ofile << " xc " << QMMMOpts.Func << '\n';
  ofile << "end" << '\n';
  ofile << "task dft optimize" << '\n';
  ofile << "task esp" << '\n';
  ofile.flush();
  ofile.close();
  //Calculate energy
  call.str("");
  if (Ncpus > 1)
  {
    call << "mpirun -n " << Ncpus << " ";
  }
  call << "nwchem QMMM_" << Bead << ".nw";
  call << " > QMMM_" << Bead << ".log";
  GlobalSys = system(call.str().c_str());
  //Parse output
  call.str("");
  call << "QMMM_" << Bead << ".log";
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
  call << "QMMM_" << Bead << ".q";
  ifile.open(call.str().c_str(),ios_base::in);
  if (ifile.good())
  {
    getline(ifile,dummy);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PAregion)
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
    cerr << " FLUKE will attempt to continue...";
    cerr << '\n';
    E = HugeNum; //Large number to reject step
    cerr.flush(); //Print warning immediately
  }
  
  //Change units
  E *= Har2eV;
  return E;
};

