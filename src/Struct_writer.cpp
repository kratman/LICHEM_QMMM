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
  if (Nmm == 0)
  {
    //Skip blank charge files
    UseChargeFile = 0;
  }
  //Initialize multipoles and center of mass
  bool FirstCharge = 1; //Always write the first charge
  Coord QMCOM;
  if (!UseChargeFile)
  {
    if (PBCon or QMMMOpts.UseLREC)
    {
      QMCOM = FindQMCOM(Struct,QMMMOpts,Bead);
    }
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
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].z,16);
      call << '\n';
    }
    if (Struct[i].PBregion)
    {
      call << "F";
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].z,16);
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
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.UseLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = Struct[i].P[Bead].x-QMCOM.x;
            dy = Struct[i].P[Bead].y-QMCOM.y;
            dz = Struct[i].P[Bead].z-QMCOM.z;
            DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xshft = DistCent.x-dx;
              yshft = DistCent.y-dy;
              zshft = DistCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.UseLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.VecMag();
            if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
            {
              //Scale the charge
              rcom = sqrt(rcom);
              double scrqA,scrqB; //Temporary variables
              //Calculate temp. variables
              scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
              scrqB = -3*scrqA*scrqA;
              scrqA *= 2*scrqA*scrqA;
              //Combine temp. variables
              scrqA += scrqB+1;
              //Set the scale factor
              scrq -= pow(scrqA,QMMMOpts.LRECPow);
            }
            else
            {
              //Delete the charge
              scrq = 0;
            }
          }
          if ((scrq > 0) or FirstCharge)
          {
            FirstCharge = 0; //Skips writing the remaining zeros
            call << " ";
            call << LICHEMFormFloat(Struct[i].P[Bead].x+xshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].P[Bead].y+yshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].P[Bead].z+zshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16);
            call << '\n';
          }
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
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.UseLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = Struct[i].P[Bead].x-QMCOM.x;
            dy = Struct[i].P[Bead].y-QMCOM.y;
            dz = Struct[i].P[Bead].z-QMCOM.z;
            DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xshft = DistCent.x-dx;
              yshft = DistCent.y-dy;
              zshft = DistCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.UseLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.VecMag();
            if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
            {
              //Scale the charge
              rcom = sqrt(rcom);
              double scrqA,scrqB; //Temporary variables
              //Calculate temp. variables
              scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
              scrqB = -3*scrqA*scrqA;
              scrqA *= 2*scrqA*scrqA;
              //Combine temp. variables
              scrqA += scrqB+1;
              //Set the scale factor
              scrq -= pow(scrqA,QMMMOpts.LRECPow);
            }
            else
            {
              //Delete the charge
              scrq = 0;
            }
          }
          if ((scrq > 0) or FirstCharge)
          {
            FirstCharge = 0; //Skips writing the remaining zeros
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x1+xshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y1+yshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z1+zshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x2+xshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y2+yshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z2+zshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x3+xshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y3+yshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z3+zshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x4+xshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y4+yshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z4+zshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x5+xshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y5+yshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z5+zshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16);
            call << '\n';
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x6+xshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y6+yshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z6+zshft,16);
            call << " ";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16);
            call << '\n';
          }
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
  if (Nmm == 0)
  {
    //Skip blank charge files
    UseChargeFile = 0;
  }
  //Initialize multipoles and center of mass
  bool FirstCharge = 1; //Always write the first charge
  Coord QMCOM;
  if (!UseChargeFile)
  {
    if (PBCon or QMMMOpts.UseLREC)
    {
      QMCOM = FindQMCOM(Struct,QMMMOpts,Bead);
    }
    if (AMOEBA)
    {
      if (TINKER)
      {
        //Set up multipoles
        RotateTINKCharges(Struct,Bead);
      }
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
    if (Struct[i].QMregion)
    {
      ofile << " " << Struct[i].QMTyp;
      ofile << " " << (Struct[i].P[Bead].x);
      ofile << " " << (Struct[i].P[Bead].y);
      ofile << " " << (Struct[i].P[Bead].z);
      ofile << '\n';
    }
    if (Struct[i].PBregion)
    {
      ofile << " " << "F2pb";
      ofile << " " << (Struct[i].P[Bead].x);
      ofile << " " << (Struct[i].P[Bead].y);
      ofile << " " << (Struct[i].P[Bead].z);
      ofile << '\n';
    }
  }
  ofile << "end" << '\n';
  if (CheckFile("BASIS"))
  {
    //Add basis set and ecp info
    ifile.open("BASIS",ios_base::in);
    if (ifile.good())
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
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.UseLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = Struct[i].P[Bead].x-QMCOM.x;
            dy = Struct[i].P[Bead].y-QMCOM.y;
            dz = Struct[i].P[Bead].z-QMCOM.z;
            DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xshft = DistCent.x-dx;
              yshft = DistCent.y-dy;
              zshft = DistCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.UseLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.VecMag();
            if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
            {
              //Scale the charge
              rcom = sqrt(rcom);
              double scrqA,scrqB; //Temporary variables
              //Calculate temp. variables
              scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
              scrqB = -3*scrqA*scrqA;
              scrqA *= 2*scrqA*scrqA;
              //Combine temp. variables
              scrqA += scrqB+1;
              //Set the scale factor
              scrq -= pow(scrqA,QMMMOpts.LRECPow);
            }
            else
            {
              //Delete the charge
              scrq = 0;
            }
          }
          if ((scrq > 0) or FirstCharge)
          {
            FirstCharge = 0; //Skips writing the remaining zeros
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].x+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].y+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].z+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16);
            ofile << '\n';
          }
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
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.UseLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = Struct[i].P[Bead].x-QMCOM.x;
            dy = Struct[i].P[Bead].y-QMCOM.y;
            dz = Struct[i].P[Bead].z-QMCOM.z;
            DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xshft = DistCent.x-dx;
              yshft = DistCent.y-dy;
              zshft = DistCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.UseLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.VecMag();
            if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
            {
              //Scale the charge
              rcom = sqrt(rcom);
              double scrqA,scrqB; //Temporary variables
              //Calculate temp. variables
              scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
              scrqB = -3*scrqA*scrqA;
              scrqA *= 2*scrqA*scrqA;
              //Combine temp. variables
              scrqA += scrqB+1;
              //Set the scale factor
              scrq -= pow(scrqA,QMMMOpts.LRECPow);
            }
            else
            {
              //Delete the charge
              scrq = 0;
            }
          }
          if ((scrq > 0) or FirstCharge)
          {
            FirstCharge = 0; //Skips writing the remaining zeros
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x1+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y1+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z1+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x2+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y2+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z2+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x3+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y3+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z3+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x4+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y4+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z4+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x5+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y5+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z5+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x6+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y6+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z6+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16);
            ofile << '\n';
          }
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
  if (Nmm == 0)
  {
    //Skip blank charge files
    UseChargeFile = 0;
  }
  //Initialize multipoles and center of mass
  Coord QMCOM;
  if (!UseChargeFile)
  {
    if (PBCon or QMMMOpts.UseLREC)
    {
      QMCOM = FindQMCOM(Struct,QMMMOpts,Bead);
    }
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
  call << "LICHM_" << Bead << ".180";
  UseCheckPoint = CheckFile(call.str());
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
  if (QMMMOpts.Spin == 1)
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
  //NB: MOs->180
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
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].x,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].y,16);
      call << " " << LICHEMFormFloat(Struct[i].P[Bead].z,16);
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
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.UseLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = Struct[i].P[Bead].x-QMCOM.x;
            dy = Struct[i].P[Bead].y-QMCOM.y;
            dz = Struct[i].P[Bead].z-QMCOM.z;
            DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xshft = DistCent.x-dx;
              yshft = DistCent.y-dy;
              zshft = DistCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.UseLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.VecMag();
            if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
            {
              //Scale the charge
              rcom = sqrt(rcom);
              double scrqA,scrqB; //Temporary variables
              //Calculate temp. variables
              scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
              scrqB = -3*scrqA*scrqA;
              scrqA *= 2*scrqA*scrqA;
              //Combine temp. variables
              scrqA += scrqB+1;
              //Set the scale factor
              scrq -= pow(scrqA,QMMMOpts.LRECPow);
            }
            else
            {
              //Delete the charge
              scrq = 0;
            }
          }
          if (scrq > 0)
          {
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16) << ",";
            call << LICHEMFormFloat(Struct[i].P[Bead].x+xshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].P[Bead].y+yshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].P[Bead].z+zshft,16);
            call << ")" << '\n';
          }
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
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.UseLREC)
          {
            //Initialize displacements
            double dx,dy,dz; //Starting displacements
            dx = Struct[i].P[Bead].x-QMCOM.x;
            dy = Struct[i].P[Bead].y-QMCOM.y;
            dz = Struct[i].P[Bead].z-QMCOM.z;
            DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
            //Calculate the shift in positions
            //NB: Generally this work out to be +/- {Lx,Ly,Lz}
            if (PBCon)
            {
              xshft = DistCent.x-dx;
              yshft = DistCent.y-dy;
              zshft = DistCent.z-dz;
            }
          }
          //Check for long-range corrections
          double scrq = 1;
          if (QMMMOpts.UseLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.VecMag();
            if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
            {
              //Scale the charge
              rcom = sqrt(rcom);
              double scrqA,scrqB; //Temporary variables
              //Calculate temp. variables
              scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
              scrqB = -3*scrqA*scrqA;
              scrqA *= 2*scrqA*scrqA;
              //Combine temp. variables
              scrqA += scrqB+1;
              //Set the scale factor
              scrq -= pow(scrqA,QMMMOpts.LRECPow);
            }
            else
            {
              //Delete the charge
              scrq = 0;
            }
          }
          if (scrq > 0)
          {
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x1+xshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y1+yshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z1+zshft,16);
            call << ")" << '\n';
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x2+xshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y2+yshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z2+zshft,16);
            call << ")" << '\n';
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x3+xshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y3+yshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z3+zshft,16);
            call << ")" << '\n';
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x4+xshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y4+yshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z4+zshft,16);
            call << ")" << '\n';
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x5+xshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y5+yshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z5+zshft,16);
            call << ")" << '\n';
            call << "Chrgfield.extern.addCharge(";
            call << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].x6+xshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].y6+yshft,16) << ",";
            call << LICHEMFormFloat(Struct[i].PC[Bead].z6+zshft,16);
            call << ")" << '\n';
          }
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
  xyzFilename = "NOFILE"; //Global XYZ filename
  //Read arguments
  Natoms = 0; //For safety
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
      xyzFilename = file.str();
      posfile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if (!CheckFile(xyzFilename))
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

