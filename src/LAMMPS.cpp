/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for LAMMPS.

 Reference for LAMMPS:
 Plimpton, Steve, J. Comp. Phys., 1, 117, 1, (1995)

*/

/*!
  \ingroup LAMMPS
*/
///@{

//MM utility functions
void LAMMPSTopology(vector<QMMMAtom>& QMMMData, stringstream& topology,
                    int bead)
{
  //Function to write bond and angle information for LAMMPS
  
  return;
};

//MM wrapper functions
double LAMMPSEnergy(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                    int bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream outFile,inFile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0.0;
  int ct; //Generic counter
  //Construct LAMMPS data file
  call.str("");
  call << "LICHM_" << bead << ".data";
  outFile.open(call.str().c_str(),ios_base::out);
  inFile.open("DATA",ios_base::in);
  call.str("");
  call << '\n';
  if (PBCon)
  {
    call << "0.0 " << Lx << " xlo xhi" << '\n';
    call << "0.0 " << Ly << " ylo yhi" << '\n';
    call << "0.0 " << Lz << " zlo zhi" << '\n';
  }
  else
  {
    call << (-1*Lx) << " " << Lx << " xlo xhi" << '\n';
    call << (-1*Ly) << " " << Ly << " ylo yhi" << '\n';
    call << (-1*Lz) << " " << Lz << " zlo zhi" << '\n';
  }
  while (!inFile.eof())
  {
     //Copy the potential line by line
     getline(inFile,dummy);
     call << dummy << '\n';
  }
  inFile.close();
  //Add PBC info
  call << "Atoms";
  call << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    call << (QMMMData[i].id+1);
    call << " 1 "; //Dummy molecule ID
    call << QMMMData[i].numTyp << " ";
    if (QMMMData[i].QMRegion || QMMMData[i].PBRegion || QMMMData[i].BARegion)
    {
      //Add zero charge
      call << 0.0;
    }
    else
    {
      call << QMMMData[i].MP[bead].q;
    }
    call << " ";
    call << QMMMData[i].P[bead].x;
    call << " ";
    call << QMMMData[i].P[bead].y;
    call << " ";
    call << QMMMData[i].P[bead].z;
    call << '\n';
  }
  call << '\n';
  inFile.open("TOPO",ios_base::in);
  while (!inFile.eof())
  {
     //Copy the potential line by line
     getline(inFile,dummy);
     call << dummy << '\n';
  }
  inFile.close();
  outFile << call.str();
  outFile.flush();
  outFile.close();
  //Construct input file
  call.str("");
  call << "LICHM_" << bead << ".in";
  outFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "atom_style full" << '\n';
  call << "units metal"; //eV,Ang,ps,bar,K
  call << '\n';
  if (PBCon)
  {
    call << "boundary p p p" << '\n';
  }
  else
  {
    call << "boundary s s s" << '\n';
  }
  call << '\n';
  call << "read_data LICHM_";
  call << bead << ".data";
  call << '\n';
  inFile.open("POTENTIAL",ios_base::in);
  while (!inFile.eof())
  {
     //Copy the potential line by line
     getline(inFile,dummy);
     call << dummy << '\n';
  }
  inFile.close();
  if (Nqm > 0)
  {
    //Partition atoms into groups
    call << "group qm id "; //QM and PB
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
      {
        call << (QMMMData[i].id+1); //LAMMPS id
        ct += 1;
        if (ct == 10)
        {
          //Break up lines
          call << " \\" << '\n';
          ct = 0;
        }
        else
        {
          //Add spaces
          call << " ";
        }
      }
    }
    if (ct != 0)
    {
      call << '\n';
    }
    //Partition atoms into MM group
    call << "group mm id "; //MM and BA
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].MMRegion || QMMMData[i].BARegion)
      {
        call << (QMMMData[i].id+1); //LAMMPS id
        ct += 1;
        if (ct == 10)
        {
          //Break up lines
          call << " \\" << '\n';
          ct = 0;
        }
        else
        {
          //Add spaces
          call << " ";
        }
      }
    }
    if (ct != 0)
    {
      call << '\n';
    }
  }
  call << "thermo 1" << '\n';
  if (MMonly)
  {
    call << "thermo_style custom step etotal" << '\n';
  }
  if (QMMM)
  {
    call << "compute mme mm group/group mm" << '\n';
    call << "compute qmmme mm group/group qm" << '\n';
    call << "thermo_style custom step c_mme c_qmmme" << '\n';
  }
  call << "run 0" << '\n'; //Only uses the initial energy
  outFile << call.str();
  outFile.flush();
  outFile.close();
  //Run calculation
  call.str("");
  call << "lammps -suffix omp -log LICHM_";
  call << bead;
  call << ".log < LICHM_";
  call << bead;
  call << ".in > LICHMlog_";
  call << bead;
  call << ".txt";
  globalSys = system(call.str().c_str());
  //Extract energy
  exit(0);
  
  //Clean up files
  call.str("");
  call << "rm -f LICHM";
  call << "_" << bead;
  call << ".in LICHM";
  call << "_" << bead;
  call << ".data LICHM";
  call << "_" << bead;
  call << ".log LICHMlog_";
  call << bead << ".txt";
  return E;
};

double LAMMPSForces(vector<QMMMAtom>& QMMMData, VectorXd& Forces,
                    QMMMSettings& QMMMOpts, int bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  double E = 0.0;
  
  return E;
};

MatrixXd LAMMPSHessian(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                       int bead)
{
  //Function for calculating the MM Hessian for the QM atoms
  MatrixXd MMHess((3*(Nqm+Npseudo)),(3*(Nqm+Npseudo)));
  MMHess.setZero();
  
  return MMHess;
};

double LAMMPSOpt(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                 int bead)
{
  //Function for optimizing with LAMMPS
  fstream outFile,inFile; //Generic file streams
  string dummy; //Generic string
  stringstream call; //Stream for system calls and reading/writing files
  call.copyfmt(cout); //Copy settings from cout
  double E = 0.0;
  int ct; //Generic counter
  //Construct LAMMPS data file
  call.str("");
  call << "LICHM_" << bead << ".data";
  inFile.open("DATA",ios_base::in);
  while (!inFile.eof())
  {
     //Copy the potential line by line
     getline(inFile,dummy);
     call << dummy << '\n';
  }
  inFile.close();
  call << '\n';
  call << "Atoms";
  call << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    call << (QMMMData[i].id+1);
    call << " 1 "; //Dummy molecule ID
    call << QMMMData[i].numTyp << " ";
    call << QMMMData[i].MP[bead].q << " ";
    call << QMMMData[i].P[bead].x;
    call << " ";
    call << QMMMData[i].P[bead].y;
    call << " ";
    call << QMMMData[i].P[bead].z;
    call << '\n';
  }
  outFile << call.str();
  outFile.flush();
  outFile.close();
  //Construct input file
  call.str("");
  call << "LICHM_" << bead << ".in";
  outFile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "atom_style full" << '\n';
  call << "units metal"; //eV,Ang,ps,bar,K
  call << '\n';
  call << "read_data ";
  call << "LICHM_" << bead << ".data";
  call << '\n' << '\n';
  inFile.open("POTENTIAL",ios_base::in);
  while (!inFile.eof())
  {
     //Copy the potential line by line
     getline(inFile,dummy);
     call << dummy << '\n';
  }
  inFile.close();
  call << '\n';
  if (Nqm > 0)
  {
    //Partition atoms into groups
    call << "group qm id "; //QM and PB
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
      {
        call << (QMMMData[i].id+1); //LAMMPS id
        ct += 1;
        if (ct == 10)
        {
          //Break up lines
          call << " \\" << '\n';
          ct = 0;
        }
        else
        {
          //Add spaces
          call << " ";
        }
      }
    }
    if (ct != 0)
    {
      call << '\n';
    }
    //Partition atoms into MM group
    call << "group mm id "; //MM and BA
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].MMRegion || QMMMData[i].BARegion)
      {
        call << (QMMMData[i].id+1); //LAMMPS id
        ct += 1;
        if (ct == 10)
        {
          //Break up lines
          call << " \\" << '\n';
          ct = 0;
        }
        else
        {
          //Add spaces
          call << " ";
        }
      }
    }
    if (ct != 0)
    {
      call << '\n';
    }
  }
  call << "thermo 1" << '\n';
  call << "thermo_style step etotal" << '\n';
  call << "compute mme mm pe" << '\n';
  call << "compute qmmme mm group/group qm" << '\n';
  //FIX THIS!!!!
  call << "min_style cg" << '\n';
  call << "minimize ";
  //Energy tolerance
  call << QMMMOpts.MMOptTol << " ";
  //Force tolerance
  call << QMMMOpts.MMOptTol << " ";
  //Max steps
  call << QMMMOpts.maxOptSteps << " ";
  //Max force iterations
  call << (3*Natoms*QMMMOpts.maxOptSteps);
  call << '\n';
  outFile << call.str();
  outFile.flush();
  outFile.close();
  //Run calculation
  call.str("");
  call << "lammps -suffix omp -log ";
  call << "LICHM_" << bead;
  call << "< LICHM_" << bead;
  call << ".in > LICHMlog_" << bead;
  call << ".txt";
  globalSys = system(call.str().c_str());
  //Extract new geometry
  
  //Clean up files
  call.str("");
  call << "rm -f LICHM";
  call << "_" << bead;
  call << ".in LICHM";
  call << "_" << bead;
  call << ".data LICHM";
  call << "_" << bead;
  call << ".log LICHMlog_";
  call << bead << ".txt";
  return E;
};

//End of file group
///@}

