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
  fstream inFile,outFile; //Generic file names
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
    if (PBCon or QMMMOpts.useLREC)
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
  outFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "%chk=LICHM_" << Bead << ".chk";
  call << '\n';
  call << "%Mem=" << QMMMOpts.RAM;
  if (QMMMOpts.memMB)
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
  call << QMMMOpts.charge << " " << QMMMOpts.spin << '\n';
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
    inFile.open(chrgfilename.c_str(),ios_base::in);
    if (inFile.good())
    {
      while (!inFile.eof())
      {
        //Copy charge file line by line
        getline(inFile,dummy);
        call << dummy << '\n';
      }
      inFile.close();
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
          if (PBCon or QMMMOpts.useLREC)
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
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.vecMag();
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
          if (PBCon or QMMMOpts.useLREC)
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
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.vecMag();
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
  inFile.open("BASIS",ios_base::in);
  if (inFile.good())
  {
    while (!inFile.eof())
    {
      //Copy BASIS line by line, if BASIS exists
      getline(inFile,dummy);
      call << dummy << '\n';
    }
    inFile.close();
  }
  outFile << call.str();
  outFile.flush();
  outFile.close();
  return;
};

void WriteNWChemInput(vector<QMMMAtom>& Struct, string CalcTyp,
                      QMMMSettings& QMMMOpts, int Bead)
{
  //Write NWChem input files
  fstream outFile,inFile; //Generic file streams
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
    if (PBCon or QMMMOpts.useLREC)
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
  outFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "LICHM_" << Bead << ".db";
  if (CheckFile(call.str()))
  {
    outFile << "restart";
  }
  else
  {
    outFile << "start";
  }
  outFile << " LICHM_" << Bead << '\n';
  outFile << "memory " << QMMMOpts.RAM;
  if (QMMMOpts.memMB)
  {
    outFile << " mb";
  }
  else
  {
    outFile << " gb";
  }
  outFile << '\n';
  outFile << "charge " << QMMMOpts.charge << '\n';
  outFile << "geometry nocenter ";
  outFile << "noautoz noautosym" << '\n';
  for (int i=0;i<Natoms;i++)
  {
    if (Struct[i].QMregion)
    {
      outFile << " " << Struct[i].QMTyp;
      outFile << " " << (Struct[i].P[Bead].x);
      outFile << " " << (Struct[i].P[Bead].y);
      outFile << " " << (Struct[i].P[Bead].z);
      outFile << '\n';
    }
    if (Struct[i].PBregion)
    {
      outFile << " " << "F2pb";
      outFile << " " << (Struct[i].P[Bead].x);
      outFile << " " << (Struct[i].P[Bead].y);
      outFile << " " << (Struct[i].P[Bead].z);
      outFile << '\n';
    }
  }
  outFile << "end" << '\n';
  if (CheckFile("BASIS"))
  {
    //Add basis set and ecp info
    inFile.open("BASIS",ios_base::in);
    if (inFile.good())
    {
      while (!inFile.eof())
      {
        //Copy BASIS line by line, if BASIS exists
        getline(inFile,dummy);
        stringstream line(dummy);
        if (line.str() != "")
        {
          //Avoid copying extra blank lines
          outFile << dummy << '\n';
        }
      }
    }
    inFile.close();
  }
  else
  {
    outFile << "basis" << '\n';
    outFile << " * library " << QMMMOpts.basis;
    outFile << '\n';
    outFile << "end" << '\n';
  }
  if (QMMM and UseChargeFile and (Nmm > 0))
  {
    inFile.open(chrgfilename.c_str(),ios_base::in);
    if (inFile.good())
    {
      outFile << "set bq:max_nbq " << (6*(Nmm+Nbound)) << '\n';
      outFile << "bq mmchrg";
      while (!inFile.eof())
      {
        //Copy charge file line by line
        outFile << '\n'; //Avoid adding an extra blank line
        getline(inFile,dummy);
        outFile << dummy; //Print the line
      }
      inFile.close();
      outFile << "end" << '\n';
      outFile << "set bq mmchrg" << '\n';
    }
  }
  else if (QMMM and (Nmm > 0))
  {
    if (CHRG)
    {
      outFile << "set bq:max_nbq " << (Nmm+Nbound) << '\n';
      outFile << "bq mmchrg" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.useLREC)
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
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.vecMag();
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
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].x+xshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].y+yshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].z+zshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16);
            outFile << '\n';
          }
        }
      }
      outFile << "end" << '\n';
      outFile << "set bq mmchrg" << '\n';
    }
    if (AMOEBA)
    {
      outFile << "set bq:max_nbq " << (6*(Nmm+Nbound)) << '\n';
      outFile << "bq mmchrg" << '\n';
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].MMregion)
        {
          //Check PBC (minimum image convention)
          Coord DistCent; //Distance from QM COM
          double xshft = 0;
          double yshft = 0;
          double zshft = 0;
          if (PBCon or QMMMOpts.useLREC)
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
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.vecMag();
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
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x1+xshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y1+yshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z1+zshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x2+xshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y2+yshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z2+zshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x3+xshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y3+yshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z3+zshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x4+xshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y4+yshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z4+zshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x5+xshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y5+yshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z5+zshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x6+xshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y6+yshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z6+zshft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16);
            outFile << '\n';
          }
        }
      }
      outFile << "end" << '\n';
      outFile << "set bq mmchrg" << '\n';
    }
  }
  //Add DFT settings
  outFile << "dft" << '\n';
  outFile << " mult " << QMMMOpts.spin << '\n';
  outFile << " direct" << '\n';
  outFile << " grid xfine nodisk" << '\n';
  outFile << " noio" << '\n';
  outFile << " tolerances tight" << '\n';
  outFile << " xc " << QMMMOpts.func << '\n';
  //Use the checkpoint file
  call.str("");
  call << "LICHM_" << Bead << ".movecs";
  if (CheckFile(call.str()))
  {
    //Tell the DFT module to read the initial vectors
    outFile << " vectors input ";
    outFile << call.str(); //Defined above
    outFile << '\n';
  }
  outFile << "end" << '\n';
  //Set calculation type
  outFile << CalcTyp;
  //Print file
  outFile.flush();
  outFile.close();
  return;
};

void WritePSI4Input(vector<QMMMAtom>& Struct, string CalcTyp,
                    QMMMSettings& QMMMOpts, int Bead)
{
  //Write PSI4 input files
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  string dummy,chrgfilename; //Generic string
  fstream inFile,outFile; //Generic file names
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
    if (PBCon or QMMMOpts.useLREC)
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
  if (QMMMOpts.memMB)
  {
    call << " mb";
  }
  else
  {
    call << " gb";
  }
  call << '\n';
  //Set options
  if (QMMMOpts.spin == 1)
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
  call << QMMMOpts.basis << '\n';
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
  call << " " << QMMMOpts.charge;
  call << " " << QMMMOpts.spin << '\n';
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
    inFile.open(chrgfilename.c_str(),ios_base::in);
    if (inFile.good())
    {
      call << "Chrgfield = QMMM()" << '\n';
      while (!inFile.eof())
      {
        //Copy charge file line by line
        getline(inFile,dummy);
        call << dummy << '\n';
      }
      inFile.close();
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
          if (PBCon or QMMMOpts.useLREC)
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
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.vecMag();
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
          if (PBCon or QMMMOpts.useLREC)
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
          if (QMMMOpts.useLREC)
          {
            //Use the long-range correction
            double rcom = 0; //Distance from center of mass
            //Calculate the distance from the center of mass
            rcom = DistCent.vecMag();
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
        inFile.open("FIELD",ios_base::in);
        while ((!inFile.eof()) and inFile.good())
        {
          getline(inFile,dummy);
          call << dummy << '\n';
        }
        //If the file was opened, save the field
        call << "activate(QMregion)" << '\n';
        call << "psi4.set_global_option_python('EXTERN',Chrgfield.extern)";
        call << '\n';
        call << '\n';
        inFile.close();
      }
    }
  }
  //Add calculation type
  call << CalcTyp;
  //Create file
  dummy = call.str(); //Store file as a temporary variable
  call.str("");
  call << "LICHM_" << Bead << ".dat";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << dummy << '\n';
  outFile.flush();
  outFile.close();
  return;
};

//Other QM files
void WriteQMConnect(int& argc,char**& argv)
{
  //Write the connectivity input for pure QM calculations
  fstream posfile,outFile; //File streams
  string dummy; //Generic string
  xyzFilename = "NOFILE"; //Global XYZ filename
  //Read arguments
  Natoms = 0; //For safety
  bool doQuit = 0; //Exit with an error
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
        doQuit = 1;
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
    doQuit = 1;
  }
  if (!doQuit)
  {
    //Write connectivity information
    outFile.open("connect.inp",ios_base::out);
    posfile >> Natoms; //Number of atoms
    for (int i=0;i<Natoms;i++)
    {
      //Read the atom type
      string AtTyp;
      posfile >> AtTyp; //Read element
      //Clear junk position data
      posfile >> dummy >> dummy >> dummy;
      //Write connectivity line
      outFile << i << " "; //Index
      outFile << AtTyp << " "; //Element
      outFile << PTable.revTyping(AtTyp) << " "; //Atomic number
      outFile << PTable.getAtMass(AtTyp) << " "; //Mass
      outFile << "0.00 0" << '\n'; //Charge and bonds
    }
    cout << "Connectivity data written to connect.inp";
    cout << '\n' << '\n';
    cout.flush();
    outFile.flush();
    outFile.close();
  }
  //Quit
  posfile.close();
  exit(0);
  return;
};

