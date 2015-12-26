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

//QM input writers
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
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].z,16);
      call << '\n';
    }
    if (Struct[i].PBregion)
    {
      call << "F";
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].z,16);
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
          call << " " << LICHEMFormDouble(Struct[i].P[Bead].x,16);
          call << " " << LICHEMFormDouble(Struct[i].P[Bead].y,16);
          call << " " << LICHEMFormDouble(Struct[i].P[Bead].z,16);
          call << " " << LICHEMFormDouble(Struct[i].MP[Bead].q,16);
          call << '\n';
        }
      }
      if (Nmm > 0)
      {
        call << '\n'; //Blank line needed
      }
    }
    if (AMOEBA)
    {
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].x1,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].y1,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].z1,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].q1,16);
          call << '\n';
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].x2,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].y2,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].z2,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].q2,16);
          call << '\n';
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].x3,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].y3,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].z3,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].q3,16);
          call << '\n';
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].x4,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].y4,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].z4,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].q4,16);
          call << '\n';
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].x5,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].y5,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].z5,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].q5,16);
          call << '\n';
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].x6,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].y6,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].z6,16);
          call << " " << LICHEMFormDouble(Struct[i].PC[Bead].q6,16);
          call << '\n';
        }
      }
      if (Nmm > 0)
      {
        call << '\n'; //Blank line needed
      }
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
    //NWChem uses fractional coordinates
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
    if (Struct[i].QMregion)
    {
      ofile << " " << Struct[i].QMTyp;
      ofile << " " << (Struct[i].P[Bead].x*ix);
      ofile << " " << (Struct[i].P[Bead].y*iy);
      ofile << " " << (Struct[i].P[Bead].z*iz);
      ofile << '\n';
    }
    if (Struct[i].PBregion)
    {
      ofile << " " << "F2pb";
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
  if (QMMM and UseChargeFile and (Nmm > 0))
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
  else if (QMMM and (Nmm > 0))
  {
    if (CHRG)
    {
      ofile << "set bq:max_nbq " << (Nmm+Nbound) << '\n';
      ofile << "bq mmchrg" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          ofile << " " << LICHEMFormDouble(Struct[i].P[Bead].x*ix,16);
          ofile << " " << LICHEMFormDouble(Struct[i].P[Bead].y*iy,16);
          ofile << " " << LICHEMFormDouble(Struct[i].P[Bead].z*iz,16);
          ofile << " " << LICHEMFormDouble(Struct[i].MP[Bead].q,16);
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
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].x1*ix,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].y1*iy,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].z1*iz,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].q1,16);
          ofile << '\n';
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].x2*ix,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].y2*iy,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].z2*iz,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].q2,16);
          ofile << '\n';
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].x3*ix,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].y3*iy,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].z3*iz,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].q3,16);
          ofile << '\n';
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].x4*ix,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].y4*iy,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].z4*iz,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].q4,16);
          ofile << '\n';
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].x5*ix,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].y5*iy,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].z5*iz,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].q5,16);
          ofile << '\n';
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].x6*ix,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].y6*iy,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].z6*iz,16);
          ofile << " " << LICHEMFormDouble(Struct[i].PC[Bead].q6,16);
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
  ofile << " tolerances tight" << '\n';
  ofile << " xc " << QMMMOpts.Func << '\n';
  //Use the checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".movecs";
  if (CheckFile(call.str()))
  {
    //Tell the DFT module to read the initial vectors
    ofile << " vectors input ";
    ofile << call.str(); //Defined above
    ofile << '\n';
  }
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
  call << " " << QMMMOpts.Charge;
  call << " " << QMMMOpts.Spin << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion)
    {
      call << " " << Struct[i].QMTyp;
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormDouble(Struct[i].P[Bead].z,16);
      call << '\n';
    }
  }
  call << " symmetry c1" << '\n';
  call << " no_reorient" << '\n';
  call << " no_com" << '\n';
  call << "}" << '\n' << '\n';
  //Set up MM field
  if (QMMM and UseChargeFile and (Nmm > 0))
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
  else if (QMMM and (Nmm > 0))
  {
    if (CHRG)
    {
      call << "Chrgfield = QMMM()" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          call << "Chrgfield.extern.addCharge(";
          call << LICHEMFormDouble(Struct[i].MP[Bead].q,16) << ",";
          call << LICHEMFormDouble(Struct[i].P[Bead].x,16) << ",";
          call << LICHEMFormDouble(Struct[i].P[Bead].y,16) << ",";
          call << LICHEMFormDouble(Struct[i].P[Bead].z,16);
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
          call << LICHEMFormDouble(Struct[i].PC[Bead].q1,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].x1,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].y1,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].z1,16);
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << LICHEMFormDouble(Struct[i].PC[Bead].q2,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].x2,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].y2,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].z2,16);
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << LICHEMFormDouble(Struct[i].PC[Bead].q3,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].x3,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].y3,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].z3,16);
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << LICHEMFormDouble(Struct[i].PC[Bead].q4,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].x4,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].y4,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].z4,16);
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << LICHEMFormDouble(Struct[i].PC[Bead].q5,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].x5,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].y5,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].z5,16);
          call << ")" << '\n';
          call << "Chrgfield.extern.addCharge(";
          call << LICHEMFormDouble(Struct[i].PC[Bead].q6,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].x6,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].y6,16) << ",";
          call << LICHEMFormDouble(Struct[i].PC[Bead].z6,16);
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

//Other QM files
void WriteQMConnect(int& argc,char**& argv)
{
  //Write the connectivity input for pure QM calculations
  fstream posfile,ofile; //File streams
  string dummy; //Generic string
  xyzfilename = "NOFILE"; //Global XYZ filename
  //Read arguments
  Nqm = 0; //For safety
  Npseudo = 0; //For safety
  bool DoQuit = 0; //Exit with an error
  cout << "Reading LICHEM input: ";
  for (int i=0;i<argc;i++)
  {
    dummy = string(argv[i]);
    //Check regions file
    if (dummy == "-q")
    {
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open XYZ file!!!";
        cout << '\n';
        DoQuit = 1;
      }
      xyzfilename = file.str();
      posfile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if (!CheckFile(xyzfilename))
  {
    cout << "Error: Missing XYZ file!!!";
    cout << '\n' << '\n';
    DoQuit = 1;
  }
  if (!DoQuit)
  {
    //Write connectivity information
    ofile.open("connect.inp",ios_base::out);
    posfile >> Natoms; //Number of atoms
    for (int i=0;i<Natoms;i++)
    {
      //Read the atom type
      string AtTyp;
      posfile >> AtTyp; //Read element
      //Clear junk position data
      posfile >> dummy >> dummy >> dummy;
      //Write connectivity line
      ofile << i << " "; //Index
      ofile << AtTyp << " "; //Element
      ofile << PTable.RevTyping(AtTyp) << " "; //Atomic number
      ofile << PTable.GetAtMass(AtTyp) << " "; //Mass
      ofile << "0.00 0" << '\n'; //Charge and bonds
    }
    cout << "Connectivity data written to connect.inp";
    cout << '\n' << '\n';
    cout.flush();
    ofile.flush();
    ofile.close();
  }
  //Quit
  posfile.close();
  exit(0);
  return;
};

