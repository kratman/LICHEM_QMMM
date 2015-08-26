/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Functions write input files for the wrappers.

*/

void WriteGauInput(vector<QMMMAtom>& Struct, string CalcTyp,
     QMMMSettings& QMMMOpts, int Bead)
{
  //Write Gaussian input files
  stringstream call;
  call.copyfmt(cout);
  string dummy,chrgfilename; //Generic strings
  fstream ifile,ofile; //Generic file names
  //Check if a list of point-charges exists
  call.str("");
  call << "MMCharges_" << Bead << ".txt";
  bool UseChrgFile = CheckFile(call.str());
  if (UseChrgFile)
  {
    chrgfilename = call.str();
  }
  if ((AMOEBA == 1) and (TINKER == 1) and (!UseChrgFile))
  {
    //Set up multipoles
    RotateTINKCharges(Struct,Bead);
  }
  //Construct g09 input
  call.str("");
  call << "LICHM_" << Bead << ".com";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=LICHM_" << Bead << ".chk";
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
  call << "%NprocShared=" << Ncpus << '\n';
  //Add ROUTE section
  call << CalcTyp;
  //Add structure
  call << '\n'; //Blank line
  call << "QMMM" << '\n' << '\n'; //Dummy title
  call << QMMMOpts.Charge << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion)
    {
      call << Struct[i].QMTyp;
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
    if (Struct[i].PBregion)
    {
      call << "F";
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout);
      call << '\n';
    }
  }
  call << '\n'; //Blank line needed
  //Add the MM field
  if ((CHRG == 1) and (!UseChrgFile))
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].P[Bead].x;
        call << " " << setprecision(12) << Struct[i].P[Bead].y;
        call << " " << setprecision(12) << Struct[i].P[Bead].z;
        call << " " << setprecision(12) << Struct[i].MP[Bead].q;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  if ((AMOEBA == 1) and (!UseChrgFile))
  {
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        call << fixed; //Forces numbers to be floats
        call << " " << setprecision(12) << Struct[i].PC[Bead].x1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z1;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q1;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z2;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q2;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z3;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q3;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z4;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q4;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z5;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q5;
        call << '\n';
        call << " " << setprecision(12) << Struct[i].PC[Bead].x6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].y6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].z6;
        call << " " << setprecision(12) << Struct[i].PC[Bead].q6;
        call.copyfmt(cout);
        call << '\n';
      }
    }
    call << '\n'; //Blank line needed
  }
  if (UseChrgFile)
  {
    //Add charges to g09 input
    ifile.open(chrgfilename.c_str(),ios_base::in);
    while (!ifile.eof())
    {
      //Copy charges line by line
      getline(ifile,dummy);
      call << dummy << '\n';
    }
    ifile.close();
  }
  //Add basis set information from the BASIS file
  ifile.open("BASIS",ios_base::in);
  if (ifile.good())
  {
    while (!ifile.eof())
    {
      //Copy BASIS line by line, if BASIS exists
      getline(ifile,dummy);
      call << dummy << '\n';
    }
    ifile.close();
  }
  ofile << call.str();
  ofile.flush();
  ofile.close();
  return;
};

void WriteNWChemInput(vector<QMMMAtom>& Struct, string CalcTyp,
     QMMMSettings& QMMMOpts, int Bead)
{
  //Write NWChem input files
  fstream ofile,ifile;
  string dummy,chrgfilename;
  stringstream call;
  call.copyfmt(cout);
  bool UseChrgFile = 0;
  //Calculate inverse box lengths for PBC
  double ix,iy,iz; //Inverse x,y,z
  ix = 1;
  iy = 1;
  iz = 1;
  if (PBCon)
  {
    ix /= Lx;
    iy /= Ly;
    iz /= Lz;
  }
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    //Check if charges are saved
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    UseChrgFile = CheckFile(call.str());
    if (UseChrgFile)
    {
      chrgfilename = call.str();
    }
    else
    {
      //Set up multipoles
      RotateTINKCharges(Struct,Bead);
    }
  }
  //Create NWChem input
  call.str("");
  call << "LICHM_" << Bead << ".nw";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "LICHM_" << Bead << ".db";
  if (CheckFile(call.str()))
  {
    ofile << "restart";
  }
  else
  {
    ofile << "start";
  }
  ofile << " LICHM_" << Bead << '\n';
  ofile << "memory " << QMMMOpts.RAM;
  if (QMMMOpts.MemMB)
  {
    ofile << " mb";
  }
  else
  {
    ofile << " gb";
  }
  ofile << '\n';
  ofile << "charge " << QMMMOpts.Charge << '\n';
  ofile << "geometry nocenter ";
  ofile << "noautoz noautosym" << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion or Struct[i].PBregion)
    {
      ofile << " " << Struct[i].QMTyp;
      ofile << " " << (Struct[i].P[Bead].x*ix);
      ofile << " " << (Struct[i].P[Bead].y*iy);
      ofile << " " << (Struct[i].P[Bead].z*iz);
      ofile << '\n';
    }
  }
  if (PBCon)
  {
    ofile << " system crystal" << '\n';
    ofile << "  lat_a " << Lx << '\n';
    ofile << "  lat_b " << Ly << '\n';
    ofile << "  lat_c " << Lz << '\n';
    ofile << "  alpha 90.0" << '\n';
    ofile << "  beta 90.0" << '\n';
    ofile << "  gamma 90.0" << '\n';
    ofile << " end" << '\n';
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
        stringstream line(dummy);
        if (line.str() != "")
        {
          //Avoid copying extra blank lines
          ofile << dummy << '\n';
        }
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
    ofile << "set bq:max_nbq " << Nmm << '\n';
    ofile << "bq mmchrg" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        ofile << fixed; //Forces numbers to be floats
        ofile << " " << setprecision(12) << (Struct[i].P[Bead].x*ix);
        ofile << " " << setprecision(12) << (Struct[i].P[Bead].y*iy);
        ofile << " " << setprecision(12) << (Struct[i].P[Bead].z*iz);
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
    ofile << "set bq:max_nbq " << (Nmm*6) << '\n';
    ofile << "bq mmchrg" << '\n';
    if (UseChrgFile)
    {
      //Add charges to NWChem input
      ifile.open(chrgfilename.c_str(),ios_base::in);
      while (!ifile.eof())
      {
        //Copy charges line by line
        getline(ifile,dummy);
        stringstream line(dummy);
        if (line.str() != "")
        {
          //Avoid copying extra blank lines
          ofile << dummy << '\n';
        }
      }
      ifile.close();
    }
    else
    {
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          ofile << fixed; //Forces numbers to be floats
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x1*ix);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y1*iy);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z1*iz);
          ofile << " " << setprecision(12) << Struct[i].PC[Bead].q1;
          ofile << '\n';
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x2*ix);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y2*iy);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z2*iz);
          ofile << " " << setprecision(12) << Struct[i].PC[Bead].q2;
          ofile << '\n';
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x3*ix);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y3*iy);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z3*iz);
          ofile << " " << setprecision(12) << Struct[i].PC[Bead].q3;
          ofile << '\n';
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x4*ix);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y4*iy);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z4*iz);
          ofile << " " << setprecision(12) << Struct[i].PC[Bead].q4;
          ofile << '\n';
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x5*ix);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y5*iy);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z5*iz);
          ofile << " " << setprecision(12) << Struct[i].PC[Bead].q5;
          ofile << '\n';
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].x6*ix);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].y6*iy);
          ofile << " " << setprecision(12) << (Struct[i].PC[Bead].z6*iz);
          ofile << " " << setprecision(12) << Struct[i].PC[Bead].q6;
          ofile.copyfmt(cout);
          ofile << '\n';
        }
      }
    }
    ofile << "end" << '\n';
    ofile << "set bq mmchrg" << '\n';
  }
  ofile << "dft" << '\n';
  ofile << " mult " << QMMMOpts.Spin << '\n';
  ofile << " direct" << '\n';
  ofile << " grid xfine nodisk" << '\n';
  ofile << " noio" << '\n';
  ofile << " cgmin" << '\n';
  ofile << " tolerances tight" << '\n';
  ofile << " xc " << QMMMOpts.Func << '\n';
  ofile << "end" << '\n';
  //Set calculation type
  ofile << CalcTyp;
  //Print file
  ofile.flush();
  ofile.close();
  return;
};

void WritePSIInput(vector<QMMMAtom>& Struct, string CalcTyp,
     QMMMSettings& QMMMOpts, int Bead)
{
  //Write PSI4 input files
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  fstream ifile,ofile; //Generic file names
  if ((AMOEBA == 1) and (TINKER == 1))
  {
    //Set up multipoles
    RotateTINKCharges(Struct,Bead);
  }
  //Set up memory
  call.str("");
  call << "memory " << QMMMOpts.RAM;
  if (QMMMOpts.MemMB)
  {
    call << " mb";
  }
  else
  {
    call << " gb";
  }
  call << '\n';
  //Set options
  if (QMMMOpts.Spin != "1")
  {
    if ((QMMMOpts.Func == "HF") or (QMMMOpts.Func == "hf")
       or (QMMMOpts.Func == "SCF") or (QMMMOpts.Func == "scf"))
    {
      //Hartree-Fock only setting
      call << "set reference uhf" << '\n';
    }
    else
    {
      //Assume it is a DFT method
      call << "set reference uks" << '\n';
    }
  }
  else
  {
    if ((QMMMOpts.Func == "HF") or (QMMMOpts.Func == "hf")
       or (QMMMOpts.Func == "SCF") or (QMMMOpts.Func == "scf"))
    {
      //Hartree-Fock only setting
      call << "set reference rhf" << '\n';
    }
    else
    {
      //Assume it is a DFT method
      call << "set reference rks" << '\n';
    }
  }
  call << "set basis ";
  call << QMMMOpts.Basis << '\n';
  call << "set guess sad" << '\n';
  call << "set scf_type df" << '\n';
  call << '\n';
  //Set up molecules
  call << "molecule QMregion {" << '\n';
  call << "  " << QMMMOpts.Charge;
  call << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion)
    {
      call << "  " << Struct[i].QMTyp;
      call << "  " << Struct[i].P[Bead].x;
      call << "  " << Struct[i].P[Bead].y;
      call << "  " << Struct[i].P[Bead].z;
      call << '\n';
    }
  }
  call << "  symmetry c1" << '\n';
  call << "  no_reorient" << '\n';
  call << "  no_com" << '\n';
  call << "}" << '\n' << '\n';
  //Set up MM field
  if (QMMM and (CHRG == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        call << "Chrgfield.extern.addCharge(";
        call << Struct[i].MP[Bead].q << ",";
        call << Struct[i].P[Bead].x/BohrRad << ",";
        call << Struct[i].P[Bead].y/BohrRad << ",";
        call << Struct[i].P[Bead].z/BohrRad;
        call << ")" << '\n';
      }
    }
    call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  if (QMMM and (AMOEBA == 1))
  {
    call << "Chrgfield = QMMM()" << '\n';
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
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
    call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
    call << '\n';
    call << '\n';
  }
  //Add calculation type
  call << CalcTyp;
  //Create file
  dummy = call.str(); //Store file as a temporary variable
  call.str("");
  call << "LICHM_" << Bead << ".dat";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << dummy << '\n';
  ofile.flush();
  ofile.close();
  return;
};
