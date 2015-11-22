/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Functions write input files for the QM wrappers.

*/

void WriteGauInput(vector<QMMMAtom>& Struct, string CalcTyp,
     QMMMSettings& QMMMOpts, int Bead)
{
  //Write Gaussian input files
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy,chrgfilename; //Generic strings
  fstream ifile,ofile; //Generic file names
  //Check for a charge file
  bool UseChargeFile = 0;
  call.str("");
  call << "MMCharges_" << Bead << ".txt";
  chrgfilename = call.str();
  UseChargeFile = CheckFile(call.str());
  //Initialize multipoles
  if (!UseChargeFile)
  {
    if (AMOEBA)
    {
      if (TINKER)
      {
        //Set up multipoles
        RotateTINKCharges(Struct,Bead);
      }
    }
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
      call.copyfmt(cout); //Copy settings from cout
      call << '\n';
    }
    if (Struct[i].PBregion)
    {
      call << "F";
      call << fixed; //Forces numbers to be floats
      call << " " << setprecision(12) << Struct[i].P[Bead].x;
      call << " " << setprecision(12) << Struct[i].P[Bead].y;
      call << " " << setprecision(12) << Struct[i].P[Bead].z;
      call.copyfmt(cout); //Copy settings from cout
      call << '\n';
    }
  }
  call << '\n'; //Blank line needed
  //Add the MM field
  if (QMMM and UseChargeFile)
  {
    ifile.open(chrgfilename.c_str(),ios_base::in);
    if (ifile.good())
    {
      while (!ifile.eof())
      {
        //Copy charge file line by line
        getline(ifile,dummy);
        call << dummy << '\n';
      }
      ifile.close();
    }
  }
  else if (QMMM)
  {
    if (CHRG)
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
          call.copyfmt(cout); //Copy settings from cout
          call << '\n';
        }
      }
      call << '\n'; //Blank line needed
    }
    if (AMOEBA)
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
          call.copyfmt(cout); //Copy settings from cout
          call << '\n';
        }
      }
      call << '\n'; //Blank line needed
    }
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
  fstream ofile,ifile; //Generic file streams
  string dummy,chrgfilename; //Generic strings
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  //Check for a charge file
  bool UseChargeFile = 0;
  call.str("");
  call << "MMCharges_" << Bead << ".txt";
  chrgfilename = call.str();
  UseChargeFile = CheckFile(call.str());
  //Initialize multipoles
  if (!UseChargeFile)
  {
    if (AMOEBA)
    {
      if (TINKER)
      {
        //Set up multipoles
        RotateTINKCharges(Struct,Bead);
      }
    }
  }
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
  if (QMMM and UseChargeFile)
  {
    ifile.open(chrgfilename.c_str(),ios_base::in);
    if (ifile.good())
    {
      ofile << "set bq:max_nbq " << (6*(Nmm+Nbound)) << '\n';
      ofile << "bq mmchrg";
      while (!ifile.eof())
      {
        //Copy charge file line by line
        ofile << '\n'; //Avoid adding an extra blank line
        getline(ifile,dummy);
        ofile << dummy; //Print the line
      }
      ifile.close();
      ofile << "end" << '\n';
      ofile << "set bq mmchrg" << '\n';
    }
  }
  else if (QMMM)
  {
    if (CHRG)
    {
      ofile << "set bq:max_nbq " << (Nmm+Nbound) << '\n';
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
          ofile.copyfmt(cout); //Copy settings from cout
          ofile << '\n';
        }
      }
      ofile << "end" << '\n';
      ofile << "set bq mmchrg" << '\n';
    }
    if (AMOEBA)
    {
      ofile << "set bq:max_nbq " << (6*(Nmm+Nbound)) << '\n';
      ofile << "bq mmchrg" << '\n';
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
          ofile.copyfmt(cout); //Copy settings from cout
          ofile << '\n';
        }
      }
      ofile << "end" << '\n';
      ofile << "set bq mmchrg" << '\n';
    }
  }
  //Add DFT settings
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

void WritePSI4Input(vector<QMMMAtom>& Struct, string CalcTyp,
     QMMMSettings& QMMMOpts, int Bead)
{
  //Write PSI4 input files
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy,chrgfilename; //Generic string
  fstream ifile,ofile; //Generic file names
  //Check for a charge file
  bool UseChargeFile = 0;
  call.str("");
  call << "MMCharges_" << Bead << ".txt";
  chrgfilename = call.str();
  UseChargeFile = CheckFile(call.str());
  //Initialize multipoles
  if (!UseChargeFile)
  {
    if (AMOEBA)
    {
      if (TINKER)
      {
        //Set up multipoles
        RotateTINKCharges(Struct,Bead);
      }
    }
  }
  //Check if there is a checkpoint file
  bool UseCheckPoint;
  call.str("");
  call << "LICHM_" << Bead << ".32";
  UseCheckPoint = CheckFile(call.str());
  if (UseCheckPoint)
  {
    //Check second file
    UseCheckPoint = 0;
    call.str("");
    call << "LICHM_" << Bead << ".180";
    UseCheckPoint = CheckFile(call.str());
  }
  //Set up memory
  call.str("");
  call << "set_num_threads(" << Ncpus << ")" << '\n';
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
  if (QMMMOpts.Spin == "1")
  {
    //Closed shell reference
    call << "set reference rhf";
  }
  else
  {
    //Open shell reference
    call << "set reference uhf";
  }
  call << '\n';
  call << "set basis ";
  call << QMMMOpts.Basis << '\n';
  if (UseCheckPoint)
  {
    call << "set guess read";
  }
  else
  {
    call << "set guess sad";
  }
  call << '\n';
  call << "set scf_type df" << '\n';
  call << '\n';
  //Keep the checkpoint files
  //NB: checkpoint->32, MOs->180
  call << "psi4_io.set_specific_path(32,'./')" << '\n';
  call << "psi4_io.set_specific_retention(32,True)" << '\n';
  call << "psi4_io.set_specific_path(180,'./')" << '\n';
  call << "psi4_io.set_specific_retention(180,True)" << '\n';
  call << '\n';
  //Set up molecules
  call << "molecule LICHM_";
  call << Bead << " {" << '\n';
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
  if (QMMM and UseChargeFile)
  {
    ifile.open(chrgfilename.c_str(),ios_base::in);
    if (ifile.good())
    {
      call << "Chrgfield = QMMM()" << '\n';
      while (!ifile.eof())
      {
        //Copy charge file line by line
        getline(ifile,dummy);
        call << dummy << '\n';
      }
      ifile.close();
      call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
      call << '\n' << '\n';
    }
  }
  else if (QMMM)
  {
    if (CHRG)
    {
      call << "Chrgfield = QMMM()" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          call << "Chrgfield.extern.addCharge(";
          call << Struct[i].MP[Bead].q << ",";
          call << Struct[i].P[Bead].x << ",";
          call << Struct[i].P[Bead].y << ",";
          call << Struct[i].P[Bead].z;
          call << ")" << '\n';
        }
      }
      call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
      call << '\n' << '\n';
    }
    if (AMOEBA)
    {
      call << "Chrgfield = QMMM()" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          call << "Chrgfield.extern.addCharge(";
          call << Struct[i].PC[Bead].q1 << ",";
          call << Struct[i].PC[Bead].x1 << ",";
          call << Struct[i].PC[Bead].y1 << ",";
          call << Struct[i].PC[Bead].z1;
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << Struct[i].PC[Bead].q2 << ",";
          call << Struct[i].PC[Bead].x2 << ",";
          call << Struct[i].PC[Bead].y2 << ",";
          call << Struct[i].PC[Bead].z2;
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << Struct[i].PC[Bead].q3 << ",";
          call << Struct[i].PC[Bead].x3 << ",";
          call << Struct[i].PC[Bead].y3 << ",";
          call << Struct[i].PC[Bead].z3;
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << Struct[i].PC[Bead].q4 << ",";
          call << Struct[i].PC[Bead].x4 << ",";
          call << Struct[i].PC[Bead].y4 << ",";
          call << Struct[i].PC[Bead].z4;
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << Struct[i].PC[Bead].q5 << ",";
          call << Struct[i].PC[Bead].x5 << ",";
          call << Struct[i].PC[Bead].y5 << ",";
          call << Struct[i].PC[Bead].z5;
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << Struct[i].PC[Bead].q6 << ",";
          call << Struct[i].PC[Bead].x6 << ",";
          call << Struct[i].PC[Bead].y6 << ",";
          call << Struct[i].PC[Bead].z6;
          call << ")" << '\n';
        }
      }
      call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
      call << '\n';
      call << '\n';
    }
    if (GEM)
    {
      //Add generic field field from a file (psithon)
      if (CheckFile("FIELD"))
      {
        //Read a block of psithon code
        ifile.open("FIELD",ios_base::in);
        while ((!ifile.eof()) and ifile.good())
        {
          getline(ifile,dummy);
          call << dummy << '\n';
        }
        //If the file was opened, save the field
        call << "activate(QMregion)" << '\n';
        call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
        call << '\n';
        call << '\n';
        ifile.close();
      }
    }
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

