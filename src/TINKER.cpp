/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE wrapper functions for TINKER.

 Reference for TINKER:
 Ponder, TINKER - Software Tools for Molecular Design

*/

//MM utility functions
void FindTINKERClasses(vector<QMMMAtom>& Struct)
{
  //Parses TINKER parameter files to find atom classes
  fstream ifile;
  string dummy; //Generic string
  string TINKKeyFile = "tinker.key";
  ifile.open(TINKKeyFile.c_str(),ios_base::in);
  if (!ifile.good())
  {
    //Exit if files do not exist
    cout << "Error: Missing tinker.key file.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  bool FileFound = 0; //Bool to break loops
  while ((!ifile.eof()) and (!FileFound))
  {
    //Detect the name of the force field file
    ifile >> dummy;
    if (dummy == "parameters")
    {
      ifile >> dummy;
      FileFound = 1;
    }
  }
  ifile.close();
  ifile.open(dummy.c_str(),ios_base::in);
  if (!FileFound)
  {
    //Exit if parameter file is not found
    cout << "Error: Cannot find TINKER parameter file.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  if (!ifile.good())
  {
    //Exit if parameter file does not exist
    cout << "Error: Cannot read TINKER ";
    cout << dummy;
    cout << " parameter file.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  int ct = 0; //Generic counter
  while (!ifile.eof())
  {
    getline(ifile,dummy);
    stringstream FullLine(dummy);
    FullLine >> dummy;
    if (dummy == "atom")
    {
      int AtType,AtClass;
      FullLine >> AtType;
      FullLine >> AtClass;
      for (int i=0;i<Natoms;i++)
      {
        if (Struct[i].NumTyp == AtType)
        {
          Struct[i].NumClass = AtClass;
          ct += 1;
        }
      }
    }
  }
  if (ct < Natoms)
  {
    cout << "Error: Atom type not found in TINKER parameters.";
    cout << '\n';
    cout << " Please check the input.";
    cout << '\n';
    cout.flush();
    exit(0);
  }
  return;
}

//MM wrapper functions
void TINKERInduced(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Function to extract induced dipoles
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  string TINKKeyFile = "tinker.key";
  int sys; //Dummy return for system calls
  int ct; //Generic counter
  //Create TINKER xyz file
  call.str("");
  call << "QMMM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n';
  if (PBCon == 1)
  {
    //Write box size
    ofile << Lx << " " << Ly << " " << Lz;
    ofile << " 90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile.copyfmt(cout);
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Create new TINKER key file
  call.str("");
  call << "cp " << TINKKeyFile << " QMMM";
  call << "_" << Bead;
  call << ".key";
  sys = system(call.str().c_str());
  //Save new keyfile name
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".key";
  TINKKeyFile = call.str(); //Save the new name
  //Add QM atoms to force field parameters list
  ofile.open(TINKKeyFile.c_str(),ios_base::app|ios_base::out);
  ofile << '\n'; //Make sure current line is empty
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ofile << "save-induced" << '\n'; //Save induced dipoles
  ofile << "thermostat berendsen" << '\n';
  ofile << "tau-temperature 0.1" << '\n';
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1))
    {
      if (ct == 0)
      {
        //Start a new active line
        ofile << "active ";
      }
      else
      {
        //Place a space to separate values
        ofile << " ";
      }
      ofile << (Struct[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        ofile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    ofile << '\n';
  }
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
    {
      //Write new multipole definition for the atom ID
      WriteTINKMpole(Struct,ofile,i,Bead);
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
  }
  ofile.flush();
  ofile.close();
  //Calculate induced dipoles using dynamic
  call.str("");
  call << "dynamic QMMM_" << Bead << ".xyz ";
  call << "1 1e-4 1e-7 2 0 > QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  //Extract induced dipoles from the MD cycle file
  call.str("");
  call << "QMMM_" << Bead << ".001u";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy); //Clear number of atoms
  while (ifile.good())
  {
    int AtNum; //Identifies which atom was polarized
    //Parse file line by line
    getline(ifile,dummy);
    stringstream line(dummy);
    //Save dipoles for later
    line >> AtNum >> dummy; //Collect atom number and clear junk
    if (line.good())
    {
      AtNum -= 1; //Fixes array indexing
      line >> Struct[AtNum].MP[Bead].IDx;
      line >> Struct[AtNum].MP[Bead].IDy;
      line >> Struct[AtNum].MP[Bead].IDz;
      //Change units from Debye to a.u.
      Struct[AtNum].MP[Bead].IDx *= Debye2au;
      Struct[AtNum].MP[Bead].IDy *= Debye2au;
      Struct[AtNum].MP[Bead].IDz *= Debye2au;
    }
  }
  ifile.close();
  //Delete junk files
  call.str("");
  call << "rm -f ";
  call << "QMMM_" << Bead << ".xyz ";
  call << "QMMM_" << Bead << ".key ";
  call << "QMMM_" << Bead << ".0* ";
  call << "QMMM_" << Bead << ".dyn ";
  call << "QMMM_" << Bead << ".log";
  sys = system(call.str().c_str());
  return;
};

double TINKERPolEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Function to extract the polarization energy
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  string TINKKeyFile = "tinker.key";
  double Epol = 0;
  double E = 0;
  int sys; //Dummy return for system calls
  int ct; //Generic counter
  //Create TINKER xyz file
  call.str("");
  call << "QMMM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << Natoms << '\n';
  if (PBCon == 1)
  {
    //Write box size
    ofile << Lx << " " << Ly << " " << Lz;
    ofile << " 90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile.copyfmt(cout);
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Create new TINKER key file
  call.str("");
  call << "cp " << TINKKeyFile << " QMMM";
  call << "_" << Bead;
  call << ".key";
  sys = system(call.str().c_str());
  //Save new keyfile name
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".key";
  TINKKeyFile = call.str(); //Save the new name
  //Add QM atoms to force field parameters list
  ofile.open(TINKKeyFile.c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1))
    {
      if (ct == 0)
      {
        //Start a new active line
        ofile << "active ";
      }
      else
      {
        //Place a space to separate values
        ofile << " ";
      }
      ofile << (Struct[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        ofile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    ofile << '\n';
  }
  for (int i=0;i<Natoms;i++)
  {
    //Add nuclear charges
    if ((Struct[i].PAregion == 1) or (Struct[i].QMregion == 1))
    {
      //Write new multipole definition for the atom ID
      WriteTINKMpole(Struct,ofile,i,Bead);
      ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
      ofile << '\n';
    }
  }
  ofile.flush();
  ofile.close();
  //Calculate QMMM energy
  call.str("");
  call << "analyze QMMM_";
  call << Bead << ".xyz E > QMMM_";
  call << Bead << ".log";
  sys = system(call.str().c_str());
  //Extract polarization energy
  call.str("");
  call << "QMMM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  bool Efound = 0;
  while (!ifile.eof())
  {
    ifile >> dummy;
    if (dummy == "Total")
    {
      ifile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        ifile >> dummy >> E;
        Efound = 1;
      }
    }
    if (dummy == "Polarization")
    {
      ifile >> Epol;
    }
  }
  if (Efound == 0)
  {
    //Warn user if no energy was found
    cout << "Warning: No MM energy found after a calculation!!!";
    cout << '\n';
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
    cout.flush(); //Print warning immediately
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  //Fix QM-MM double counting
  if (QMMMOpts.Nind > 0)
  {
    vector<double> Es; //List of energies
    vector<int> QMatoms; //List of the QM atoms
    for (int i=0;i<Natoms;i++)
    {
      //Create temporary storage arrays
      Es.push_back(0.0);
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        QMatoms.push_back(i);
      }
    }
    #pragma omp parallel for
    for (int i=0;i<Natoms;i++)
    {
      if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1))
      {
        //Store coordinates for dipoles
        Coord PoleXp = Struct[i].P[Bead];
        Coord PoleXm = Struct[i].P[Bead];
        Coord PoleYp = Struct[i].P[Bead];
        Coord PoleYm = Struct[i].P[Bead];
        Coord PoleZp = Struct[i].P[Bead];
        Coord PoleZm = Struct[i].P[Bead];
        //Adjust dipole positions
        PoleXp.x += 0.25*BohrRad;
        PoleXm.x -= 0.25*BohrRad;
        PoleYp.y += 0.25*BohrRad;
        PoleYm.y -= 0.25*BohrRad;
        PoleZp.z += 0.25*BohrRad;
        PoleZm.z -= 0.25*BohrRad;
        //Calculate energies
        for (int j=0;j<(Nqm+Npseudo);j++)
        {
          //Create variables
          double Etmp; //Energy of a dipole-charge interaction
          double r; //Distance between charges
          //Calculate X dipole interactions
          r = sqrt(CoordDist2(PoleXp,Struct[QMatoms[j]].P[Bead]));
          r /= BohrRad; //Change to a.u.
          Etmp = Struct[i].MP[Bead].IDx/(2*0.25);
          Etmp *= Struct[QMatoms[j]].MP[Bead].q;
          Etmp /= r;
          Es[i] += Etmp;
          r = sqrt(CoordDist2(PoleXm,Struct[QMatoms[j]].P[Bead]));
          r /= BohrRad; //Change to a.u.
          Etmp = Struct[i].MP[Bead].IDx/(2*0.25);
          Etmp *= Struct[QMatoms[j]].MP[Bead].q;
          Etmp /= r;
          Etmp *= -1; //Change sign for negative displacement
          Es[i] += Etmp;
          //Calculate X dipole interactions
          r = sqrt(CoordDist2(PoleYp,Struct[QMatoms[j]].P[Bead]));
          r /= BohrRad; //Change to a.u.
          Etmp = Struct[i].MP[Bead].IDy/(2*0.25);
          Etmp *= Struct[QMatoms[j]].MP[Bead].q;
          Etmp /= r;
          Es[i] += Etmp;
          r = sqrt(CoordDist2(PoleYm,Struct[QMatoms[j]].P[Bead]));
          r /= BohrRad; //Change to a.u.
          Etmp = Struct[i].MP[Bead].IDy/(2*0.25);
          Etmp *= Struct[QMatoms[j]].MP[Bead].q;
          Etmp /= r;
          Etmp *= -1; //Change sign for negative displacement
          Es[i] += Etmp;
          //Calculate Z dipole interactions
          r = sqrt(CoordDist2(PoleZp,Struct[QMatoms[j]].P[Bead]));
          r /= BohrRad; //Change to a.u.
          Etmp = Struct[i].MP[Bead].IDz/(2*0.25);
          Etmp *= Struct[QMatoms[j]].MP[Bead].q;
          Etmp /= r;
          Es[i] += Etmp;
          r = sqrt(CoordDist2(PoleZm,Struct[QMatoms[j]].P[Bead]));
          r /= BohrRad; //Change to a.u.
          Etmp = Struct[i].MP[Bead].IDz/(2*0.25);
          Etmp *= Struct[QMatoms[j]].MP[Bead].q;
          Etmp /= r;
          Etmp *= -1; //Change sign for negative displacement
          Es[i] += Etmp;
        }
      }
    }
    #pragma omp barrier
    for (int i=0;i<Natoms;i++)
    {
      //Sum up the energy in kcal/mol
      Epol -= Es[i]*Har2eV/kcal2eV;
    }
  }
  //Return polarization energy in kcal/mol
  return Epol;
};

double TINKERForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  string TINKKeyFile = "tinker.key";
  double Emm = 0.0;
  int ct; //Generic counter
  int sys; //Dummy return for system calls
  //Construct MM forces input for TINKER
  call.str("");
  call << "cp " << TINKKeyFile << " ";
  call << "QMMM";
  call << "_" << Bead;
  call << ".key";
  sys = system(call.str().c_str());
  //Save new keyfile name
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".key";
  TINKKeyFile = call.str(); //Save the new name
  //Add QM atoms to force field parameters list
  ofile.open(TINKKeyFile.c_str(),ios_base::app|ios_base::out);
  ofile << '\n';
  ofile << "#QM force field parameters"; //Marks the changes
  ofile << '\n';
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add active atoms
    if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
    {
      if (ct == 0)
      {
        //Start a new active line
        ofile << "active ";
      }
      else
      {
        //Place a space to separate values
        ofile << " ";
      }
      ofile << (Struct[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate an active line
        ct = 0;
        ofile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing actives line
    ofile << '\n';
  }
  ofile << "group-inter" << '\n'; //Modify interactions
  ct = 0; //Generic counter
  for (int i=0;i<Natoms;i++)
  {
    //Add group 1 atoms
    if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
    {
      if (ct == 0)
      {
        //Start a new group line
        ofile << "group 1 ";
      }
      else
      {
        //Place a space to separate values
        ofile << " ";
      }
      ofile << (Struct[i].id+1);
      ct += 1;
      if (ct == 10)
      {
        //terminate a group line
        ct = 0;
        ofile << '\n';
      }
    }
  }
  if (ct != 0)
  {
    //Terminate trailing group line
    ofile << '\n';
  }
  if (ct != 0)
  {
    //Terminate trailing group line
    ofile << '\n';
  }
  if (CHRG == 1)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        //New charges are needed for QM and PA atoms
        ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
        ofile << "0.0"; //Delete charges
        ofile << '\n';
      }
    }
  }
  if (AMOEBA == 1)
  {
    for (int i=0;i<Natoms;i++)
    {
      //Add nuclear charges
      if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
      {
        double qi = 0;
        //remove charge
        qi = Struct[i].MP[Bead].q;
        Struct[i].MP[Bead].q = 0;
        WriteTINKMpole(Struct,ofile,i,Bead);
        Struct[i].MP[Bead].q += qi; //Restore charge
        ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
        ofile << '\n';
      }
    }
  }
  ofile.flush();
  ofile.close();
  //Create TINKER xyz file from the structure
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon == 1)
  {
    //Write box size
    ofile << Lx << " " << Ly << " " << Lz;
    ofile << " 90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile.copyfmt(cout);
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Run MM
  call.str("");
  call << "testgrad " << "QMMM";
  call << "_" << Bead;
  call << ".xyz Y N N > QMMM";
  call << "_" << Bead;
  call << ".grad";
  sys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".grad";
  MMgrad.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool GradDone = 0;
  while ((!MMgrad.eof()) and (!GradDone))
  {
    getline(MMgrad,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Type")
    {
      line >> dummy >> dummy;
      if (dummy == "dE/dX")
      {
        GradDone = 1; //Not grad school, that lasts forever
        getline(MMgrad,dummy);
        for (int i=0;i<(Nqm+Npseudo);i++)
        {
          double Fx = 0;
          double Fy = 0;
          double Fz = 0;
          //Convoluted, but "easy"
          getline(MMgrad,dummy);
          stringstream line(dummy);
          line >> dummy >> dummy; //Clear junk
          line >> Fx;
          line >> Fy;
          line >> Fz;
          //Switch to eV/A and change sign
          if (abs(Fx) >= 1e-4)
          {
            Forces[i].x += -1*Fx*kcal2eV;
          }
          if (abs(Fy) >= 1e-4)
          {
            Forces[i].y += -1*Fy*kcal2eV;
          }
          if (abs(Fz) >= 1e-4)
          {
            Forces[i].z += -1*Fz*kcal2eV;
          }
        }
      }
    }
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "Energy")
      {
        //Collect partial MM energy
        line >> dummy >> Emm;
      }
    }
  }
  MMgrad.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".key";
  call << " QMMM_" << Bead << ".grad";
  sys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKERMMForces(vector<QMMMAtom>& Struct, vector<Coord>& MMForces,
     QMMMSettings& QMMMOpts, int Bead)
{
  //A routine to extract MM forces from TINKER
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  string TINKKeyFile = "tinker.key";
  double Emm = 0.0;
  int ct; //Generic counter
  int sys; //Dummy return for system calls
  //Copy the original key file and make changes
  if (QMMM == 1)
  {
    call.str("");
    call << "cp " << TINKKeyFile << " QMMM";
    call << "_" << Bead;
    call << ".key";
    sys = system(call.str().c_str());
    //Save new keyfile name
    call.str("");
    call << "QMMM";
    call << "_" << Bead;
    call << ".key";
    TINKKeyFile = call.str(); //Save the new name
    //Add QM atoms to force field parameters list
    ofile.open(TINKKeyFile.c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1)
         or (Struct[i].PAregion == 1))
      {
        if (Struct[i].Frozen == 0)
        {
          if (ct == 0)
          {
            //Start a new active line
            ofile << "active ";
          }
          else
          {
            //Place a space to separate values
            ofile << " ";
          }
          ofile << (Struct[i].id+1);
          ct += 1;
          if (ct == 10)
          {
            //terminate an active line
            ct = 0;
            ofile << '\n';
          }
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      ofile << '\n';
    }
    if (CHRG == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion == 1)
        {
          //New charges are only needed for QM atoms
          ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
          ofile << Struct[i].MP[Bead].q;
          ofile << '\n';
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
        {
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
        }
      }
    }
    ofile.flush();
    ofile.close();
  }
  //Create TINKER xyz file from the structure
  call.str("");
  call << "QMMM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon == 1)
  {
    //Write box size
    ofile << Lx << " " << Ly << " " << Lz;
    ofile << " 90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile.copyfmt(cout);
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Run MM
  call.str("");
  call << "testgrad " << "QMMM";
  call << "_" << Bead;
  call << ".xyz Y N N > QMMM";
  call << "_" << Bead;
  call << ".grad";
  sys = system(call.str().c_str());
  //Collect MM forces
  fstream MMgrad; //QMMM output
  //Open files
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".grad";
  MMgrad.open(call.str().c_str(),ios_base::in);
  //Read derivatives
  bool GradDone = 0;
  while ((!MMgrad.eof()) and (!GradDone))
  {
    getline(MMgrad,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Type")
    {
      line >> dummy >> dummy;
      if (dummy == "dE/dX")
      {
        GradDone = 1; //Not grad school, that lasts forever
        getline(MMgrad,dummy);
        for (int i=0;i<Natoms;i++)
        {
          if ((Struct[i].Frozen == 0) and (Struct[i].QMregion == 0))
          {
            double Fx = 0;
            double Fy = 0;
            double Fz = 0;
            //Convoluted, but "easy"
            getline(MMgrad,dummy);
            stringstream line(dummy);
            line >> dummy >> dummy; //Clear junk
            line >> Fx;
            line >> Fy;
            line >> Fz;
            //Switch to eV/A and change sign
            if (abs(Fx) >= 1e-4)
            {
              MMForces[i].x += -1*Fx*kcal2eV;
            }
            if (abs(Fy) >= 1e-4)
            {
              MMForces[i].y += -1*Fy*kcal2eV;
            }
            if (abs(Fz) >= 1e-4)
            {
              MMForces[i].z += -1*Fz*kcal2eV;
            }
          }
        }
      }
    }
    if (dummy == "Total")
    {
      line >> dummy >> dummy;
      if (dummy == "Energy")
      {
        //Collect partial MM energy
        line >> dummy >> Emm;
      }
    }
  }
  MMgrad.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".key";
  call << " QMMM_" << Bead << ".grad";
  sys = system(call.str().c_str());
  //Return
  Emm *= kcal2eV;
  return Emm;
};

double TINKEREnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs TINKER MM
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  string TINKKeyFile = "tinker.key";
  double Epol = 0;
  double E = 0;
  int sys; //Dummy return for system calls
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  if (QMMM == 1)
  {
    call.str("");
    call << "cp " << TINKKeyFile << " QMMM";
    call << "_" << Bead;
    call << ".key";
    sys = system(call.str().c_str());
    //Save new keyfile name
    call.str("");
    call << "QMMM";
    call << "_" << Bead;
    call << ".key";
    TINKKeyFile = call.str(); //Save the new name
    //Add QM atoms to force field parameters list
    ofile.open(TINKKeyFile.c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1))
      {
        if (ct == 0)
        {
          //Start a new active line
          ofile << "active ";
        }
        else
        {
          //Place a space to separate values
          ofile << " ";
        }
        ofile << (Struct[i].id+1);
        ct += 1;
        if (ct == 10)
        {
          //terminate an active line
          ct = 0;
          ofile << '\n';
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      ofile << '\n';
    }
    if (CHRG == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion == 1)
        {
          //New charges are only needed for QM atoms
          ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
          ofile << "0.0" << '\n';
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].PAregion == 1)
        {
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
        }
        if (Struct[i].QMregion == 1)
        {
          double qi = 0;
          //remove charge
          qi = Struct[i].MP[Bead].q;
          Struct[i].MP[Bead].q = 0;
          WriteTINKMpole(Struct,ofile,i,Bead);
          Struct[i].MP[Bead].q += qi; //Restore charge
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
        }
      }
    }
    ofile.flush();
    ofile.close();
  }
  //Create TINKER xyz file from the structure
  call.str("");
  call << "QMMM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon == 1)
  {
    //Write box size
    ofile << Lx << " " << Ly << " " << Lz;
    ofile << " 90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile.copyfmt(cout);
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Calculate MM potential energy
  call.str("");
  call << "analyze QMMM_";
  call << Bead << ".xyz E > QMMM_";
  call << Bead << ".log";
  sys = system(call.str().c_str());
  call.str("");
  call << "QMMM_" << Bead << ".log";
  ifile.open(call.str().c_str(),ios_base::in);
  //Read MM potential energy
  bool Efound = 0;
  while (!ifile.eof())
  {
    ifile >> dummy;
    if (dummy == "Total")
    {
      ifile >> dummy >> dummy;
      if (dummy == "Energy")
      {
        ifile >> dummy >> E;
        Efound = 1;
      }
    }
    if (dummy == "Polarization")
    {
      ifile >> Epol;
    }
  }
  if (Efound == 0)
  {
    //Warn user if no energy was found
    cout << "Warning: No MM energy found after a calculation!!!";
    cout << '\n';
    cout << " FLUKE will attempt to continue...";
    cout << '\n';
    E = HugeNum; //Large number to reject step
    cout.flush(); //Print warning immediately
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  //Calculate polarization energy
  if ((AMOEBA == 1) and (QMMM == 1))
  {
    //Correct polarization energy for QMMM simulations
    E -= Epol; //Incorrect polarization energy
    E += TINKERPolEnergy(Struct,QMMMOpts,Bead);
  }
  //Change units
  E *= kcal2eV;
  return E;
};

double TINKEROpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, int Bead)
{
  //Runs TINKER MM
  fstream ofile,ifile;
  stringstream call;
  call.copyfmt(cout);
  string dummy; //Generic string
  string TINKKeyFile = "tinker.key";
  double Epol = 0;
  double E = 0;
  int sys; //Dummy return for system calls
  int ct; //Generic counter
  call.str("");
  //Copy the original key file and make changes
  if (QMMM == 1)
  {
    call.str("");
    call << "cp " << TINKKeyFile << " QMMM";
    call << "_" << Bead;
    call << ".key";
    sys = system(call.str().c_str());
    //Save new keyfile name
    call.str("");
    call << "QMMM";
    call << "_" << Bead;
    call << ".key";
    TINKKeyFile = call.str(); //Save the new name
    //Add QM atoms to force field parameters list
    ofile.open(TINKKeyFile.c_str(),ios_base::app|ios_base::out);
    ofile << '\n';
    ofile << "#QM force field parameters"; //Marks the changes
    ofile << '\n';
    ct = 0; //Generic counter
    for (int i=0;i<Natoms;i++)
    {
      //Add active atoms
      if ((Struct[i].MMregion == 1) or (Struct[i].BAregion == 1)
         or (Struct[i].PAregion == 1))
      {
        if (Struct[i].Frozen == 0)
        {
          if (ct == 0)
          {
            //Start a new active line
            ofile << "active ";
          }
          else
          {
            //Place a space to separate values
            ofile << " ";
          }
          ofile << (Struct[i].id+1);
          ct += 1;
          if (ct == 10)
          {
            //terminate an active line
            ct = 0;
            ofile << '\n';
          }
        }
      }
    }
    if (ct != 0)
    {
      //Terminate trailing actives line
      ofile << '\n';
    }
    if (CHRG == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if (Struct[i].QMregion == 1)
        {
          //New charges are only needed for QM atoms
          ofile << "charge " << (-1*(Struct[i].id+1)) << " ";
          ofile << Struct[i].MP[Bead].q;
          ofile << '\n';
        }
      }
    }
    if (AMOEBA == 1)
    {
      for (int i=0;i<Natoms;i++)
      {
        //Add nuclear charges
        if ((Struct[i].QMregion == 1) or (Struct[i].PAregion == 1))
        {
          //Write new multipole definition for the atom ID
          WriteTINKMpole(Struct,ofile,i,Bead);
          ofile << "polarize -" << (Struct[i].id+1) << " 0.0 0.0";
          ofile << '\n';
        }
      }
    }
    ofile.flush();
    ofile.close();
  }
  //Create TINKER xyz file from the structure
  call.str("");
  call << "QMMM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  if (PBCon == 1)
  {
    //Write box size
    ofile << Lx << " " << Ly << " " << Lz;
    ofile << " 90.0 90.0 90.0";
    ofile << '\n';
  }
  ct = 0; //Counter for QM atoms
  for (int i=0;i<Natoms;i++)
  {
    ofile.precision(8);
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].x;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].y;
    ofile << " ";
    ofile << setw(10) << Struct[i].P[Bead].z;
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile.copyfmt(cout);
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Run optimization
  call.str("");
  call << "newton QMMM_";
  call << Bead << ".xyz A A ";
  call << QMMMOpts.MMOptTol << " > QMMM_";
  call << Bead << ".log";
  sys = system(call.str().c_str());
  //Read new structure
  call.str("");
  call << "QMMM_" << Bead << ".xyz_2";
  ifile.open(call.str().c_str(),ios_base::in);
  getline(ifile,dummy);
  if (PBCon == 1)
  {
    getline(ifile,dummy);
  }
  for (int i=0;i<Natoms;i++)
  {
    getline(ifile,dummy);
    stringstream line(dummy);
    //Read new positions
    line >> dummy >> dummy;
    line >> Struct[i].P[Bead].x;
    line >> Struct[i].P[Bead].y;
    line >> Struct[i].P[Bead].z;
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f";
  call << " QMMM_" << Bead << ".xyz";
  call << " QMMM_" << Bead << ".log";
  call << " QMMM_" << Bead << ".xyz_*";
  call << " QMMM_" << Bead << ".key";
  sys = system(call.str().c_str());
  //Calculate new induced dipoles
  if ((AMOEBA == 1) and (QMMM))
  {
    TINKERInduced(Struct,QMMMOpts,Bead);
  }
  //Change units
  E *= kcal2eV;
  return E;
};
