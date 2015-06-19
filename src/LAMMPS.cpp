/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
###############################################################################

 LICHEM wrapper functions for LAMMPS.

 Reference for LAMMPS:
 Plimpton, Steve, J. Comp. Phys., 1, 117, 1, (1995)

*/

//MM utility functions

//MM wrapper functions
double LAMMPSForces(vector<QMMMAtom>& Struct, vector<Coord>& Forces,
       QMMMSettings& QMMMOpts, int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int ct;

  return E;
};

double LAMMPSEnergy(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Function for calculating the MM forces on a set of QM atoms
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int ct;
  //Construct LAMMPS data file
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".data";
  ifile.open("DATA",ios_base::in);
  while (!ifile.eof())
  {
     //Copy the potential line by line
     getline(ifile,dummy);
     call << dummy << '\n';
  }
  ifile.close();
  call << '\n';
  call << "Atoms";
  call << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    call << (Struct[i].id+1);
    call << " 1 "; //Dummy molecule ID
    call << Struct[i].NumTyp << " ";
    if (Struct[i].QMregion or Struct[i].PBregion or Struct[i].BAregion)
    {
      //Add zero charge
      call << 0.0;
    }
    else
    {
      call << Struct[i].MP[Bead].q;
    }
    call << " ";
    call << Struct[i].P[Bead].x;
    call << " ";
    call << Struct[i].P[Bead].y;
    call << " ";
    call << Struct[i].P[Bead].z;
    call << '\n';
  }
  ofile << call.str();
  ofile.flush();
  ofile.close();
  //Construct input file
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".in";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "atom_style full" << '\n';
  call << "units metal"; //eV,Ang,ps,bar,K
  call << '\n';
  call << "read_data QMMM_";
  call << Bead << ".data";
  call << '\n' << '\n';
  ifile.open("POTENTIAL",ios_base::in);
  while (!ifile.eof())
  {
     //Copy the potential line by line
     getline(ifile,dummy);
     call << dummy << '\n';
  }
  ifile.close();
  call << '\n';
  if (Nqm > 0)
  {
    //Partition atoms into groups
    call << "group qm id "; //QM and PB
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        call << (Struct[i].id+1); //LAMMPS id
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
      if (Struct[i].MMregion or Struct[i].BAregion)
      {
        call << (Struct[i].id+1); //LAMMPS id
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
  call << "run 1" << '\n';
  ofile << call.str();
  ofile.flush();
  ofile.close();
  //Run calculation
  call.str("");
  call << "lammps -suffix omp -log QMMM_";
  call << Bead;
  call << "< QMMM_";
  call << Bead;
  call << ".in > QMMMlog_";
  call << Bead;
  call << ".txt";
  GlobalSys = system(call.str().c_str());
  //Extract energy
  
  //Clean up files
  call.str("");
  call << "rm -f QMMM";
  call << "_" << Bead;
  call << ".in QMMM";
  call << "_" << Bead;
  call << ".data QMMM";
  call << "_" << Bead;
  call << ".log QMMMlog.txt";
  return E;
};

double LAMMPSOpt(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
       int Bead)
{
  //Function for optimizing with LAMMPS
  fstream ofile,ifile;
  string dummy; //Generic string
  stringstream call;
  call.copyfmt(cout);
  double E = 0.0;
  int ct;
  //Construct LAMMPS data file
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".data";
  ifile.open("DATA",ios_base::in);
  while (!ifile.eof())
  {
     //Copy the potential line by line
     getline(ifile,dummy);
     call << dummy << '\n';
  }
  ifile.close();
  call << '\n';
  call << "Atoms";
  call << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    call << (Struct[i].id+1);
    call << " 1 "; //Dummy molecule ID
    call << Struct[i].NumTyp << " ";
    call << Struct[i].MP[Bead].q << " ";
    call << Struct[i].P[Bead].x;
    call << " ";
    call << Struct[i].P[Bead].y;
    call << " ";
    call << Struct[i].P[Bead].z;
    call << '\n';
  }
  ofile << call.str();
  ofile.flush();
  ofile.close();
  //Construct input file
  call.str("");
  call << "QMMM";
  call << "_" << Bead;
  call << ".in";
  ofile.open(call.str().c_str(),ios_base::out);
  call.str("");
  call << "atom_style full" << '\n';
  call << "units metal"; //eV,Ang,ps,bar,K
  call << '\n';
  call << "read_data QMMM";
  call << "_" << Bead;
  call << ".data";
  call << '\n' << '\n';
  ifile.open("POTENTIAL",ios_base::in);
  while (!ifile.eof())
  {
     //Copy the potential line by line
     getline(ifile,dummy);
     call << dummy << '\n';
  }
  ifile.close();
  call << '\n';
  if (Nqm > 0)
  {
    //Partition atoms into groups
    call << "group qm id "; //QM and PB
    ct = 0;
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].QMregion or Struct[i].PBregion)
      {
        call << (Struct[i].id+1); //LAMMPS id
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
      if (Struct[i].MMregion or Struct[i].BAregion)
      {
        call << (Struct[i].id+1); //LAMMPS id
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
  call << QMMMOpts.MaxOptSteps << " ";
  //Max force iterations
  call << (3*Natoms*QMMMOpts.MaxOptSteps);
  call << '\n';
  ofile << call.str();
  ofile.flush();
  ofile.close();
  //Run calculation
  call.str("");
  call << "lammps -suffix omp -log QMMM";
  call << "_" << Bead;
  call << "< QMMM";
  call << "_" << Bead;
  call << ".in > QMMMlog";
  call << "_" << Bead;
  call << ".txt";
  GlobalSys = system(call.str().c_str());
  //Extract new geometry
  
  //Clean up files
  call.str("");
  call << "rm -f QMMM";
  call << "_" << Bead;
  call << ".in QMMM";
  call << "_" << Bead;
  call << ".data QMMM";
  call << "_" << Bead;
  call << ".log QMMMlog.txt";
  return E;
};

