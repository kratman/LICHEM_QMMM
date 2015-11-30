/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Path integral functions using QM, MM, and QMMM energies. Calls to wrappers
 are parallel over the number of beads. Other functions are mostly parallel
 over the number of atoms.

*/

//Path integral Monte Carlo functions
double Get_PI_Espring(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Calculate total harmonic PI ring energy
  double E = 0.0;
  double w0 = 1/(QMMMOpts.Beta*hbar);
  w0 *= w0*ToeV*QMMMOpts.Nbeads;
  #pragma omp parallel for reduction(+:E)
  for (int i=0;i<Natoms;i++)
  {
    Struct[i].Ep = 0.0;
    double w = w0*Struct[i].m;
    for (int j=0;j<QMMMOpts.Nbeads;j++)
    {
      //Bead energy, one bond to avoid double counting
      int j2 = j-1;
      if (j2 == -1)
      {
        j2 = QMMMOpts.Nbeads-1; //Ring PBC
      }
      //Calculate displacement with PBC
      double dr2 = CoordDist2(Struct[i].P[j],Struct[i].P[j2]).VecMag();
      Struct[i].Ep += 0.5*w*dr2; //Harmonic energy
    }
    E += Struct[i].Ep; //Save energy
  }
  #pragma omp barrier
  return E;
};

double Get_PI_Epot(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Potential for all beads
  double E = 0.0;
  //Fix parallel for classical MC
  int MCThreads = Nthreads;
  if (QMMMOpts.Nbeads == 1)
  {
    MCThreads = 1;
  }
  //Calculate energy
  #pragma omp parallel for num_threads(MCThreads) reduction(+:E,QMTime,MMTime)
  for (int i=0;i<QMMMOpts.Nbeads;i++)
  {
    //Run the wrappers for all beads
    double Es = 0.0;
    //Timer variables
    int t_qm_start = 0;
    int t_mm_start = 0;
    int Times_qm = 0;
    int Times_mm = 0;
    //Calculate QM energy
    if (Gaussian)
    {
      t_qm_start = (unsigned)time(0);
      Es += GaussianEnergy(Struct,QMMMOpts,i);
      Times_qm += (unsigned)time(0)-t_qm_start;
    }
    if (PSI4)
    {
      t_qm_start = (unsigned)time(0);
      Es += PSI4Energy(Struct,QMMMOpts,i);
      Times_qm += (unsigned)time(0)-t_qm_start;
      //Delete annoying useless files
      GlobalSys = system("rm -f psi.* timer.*");
    }
    if (NWChem)
    {
      t_qm_start = (unsigned)time(0);
      Es += NWChemEnergy(Struct,QMMMOpts,i);
      Times_qm += (unsigned)time(0)-t_qm_start;
    }
    //Calculate MM energy
    if (TINKER)
    {
      t_mm_start = (unsigned)time(0);
      Es += TINKEREnergy(Struct,QMMMOpts,i);
      Times_mm += (unsigned)time(0)-t_mm_start;
    }
    if (AMBER)
    {
      t_mm_start = (unsigned)time(0);
      Es += AMBEREnergy(Struct,QMMMOpts,i);
      Times_mm += (unsigned)time(0)-t_mm_start;
    }
    if (LAMMPS)
    {
      t_mm_start = (unsigned)time(0);
      Es += LAMMPSEnergy(Struct,QMMMOpts,i);
      Times_mm += (unsigned)time(0)-t_mm_start;
    }
    //Add temp variables to the totals
    E += Es;
    QMTime += Times_qm;
    MMTime += Times_mm;
  }
  #pragma omp barrier
  E /= QMMMOpts.Nbeads;
  return E;
};

bool MCMove(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts, double& Emc)
{
  //Function to perform Monte Carlo moves and accept/reject the moves
  bool acc = 0; //Accept or reject
  //Copy Struct
  vector<QMMMAtom> Struct2;
  Struct2 = Struct;
  //Pick random move and apply PBC
  double randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum > (1-CentProb))
  {
    //Move a centroid
    int p;
    bool FrozenAt = 1;
    while (FrozenAt)
    {
      //Make sure the atom is not frozen
      p = (rand()%Natoms);
      if (Struct[p].Frozen == 0)
      {
        FrozenAt = 0;
      }
    }
    double randx = (((double)rand())/((double)RAND_MAX));
    double randy = (((double)rand())/((double)RAND_MAX));
    double randz = (((double)rand())/((double)RAND_MAX));
    double dx = 2*(randx-0.5)*step*CentRatio;
    double dy = 2*(randy-0.5)*step*CentRatio;
    double dz = 2*(randz-0.5)*step*CentRatio;
    //Update positions
    #pragma omp parallel for
    for (int i=0;i<QMMMOpts.Nbeads;i++)
    {
      Struct2[p].P[i].x += dx;
      Struct2[p].P[i].y += dy;
      Struct2[p].P[i].z += dz;
    }
    #pragma omp barrier
  }
  if (randnum < BeadProb)
  {
    //Move all beads in a centroid
    int p;
    bool FrozenAt = 1;
    while (FrozenAt)
    {
      //Make sure the atom is not frozen
      p = (rand()%Natoms);
      if (Struct[p].Frozen == 0)
      {
        FrozenAt = 0;
      }
    }
    for (int i=0;i<QMMMOpts.Nbeads;i++)
    {
      //Randomly displace each bead
      double randx = (((double)rand())/((double)RAND_MAX));
      double randy = (((double)rand())/((double)RAND_MAX));
      double randz = (((double)rand())/((double)RAND_MAX));
      double dx = 2*(randx-0.5)*step;
      double dy = 2*(randy-0.5)*step;
      double dz = 2*(randz-0.5)*step;
      Struct2[p].P[i].x += dx;
      Struct2[p].P[i].y += dy;
      Struct2[p].P[i].z += dz;
    }
  }
  //Initialize energies
  double Eold = QMMMOpts.Eold;
  double Enew = 0;
  //Save box lengths
  double Lxtmp = Lx;
  double Lytmp = Ly;
  double Lztmp = Lz;
  //Attempt a volume move
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < VolProb)
  {
    //Anisotropic volume change
    if (Isotrop == 0)
    {
      randnum = (((double)rand())/((double)RAND_MAX));
      double Lmin = 0.97*Lx;
      double Lmax = 1.03*Lx;
      Lx = Lmin+randnum*(Lmax-Lmin);
      randnum = (((double)rand())/((double)RAND_MAX));
      Lmin = 0.97*Ly;
      Lmax = 1.03*Ly;
      Ly = Lmin+randnum*(Lmax-Lmin);
      randnum = (((double)rand())/((double)RAND_MAX));
      Lmin = 0.97*Lz;
      Lmax = 1.03*Lz;
      Lz = Lmin+randnum*(Lmax-Lmin);
    }
    //Isotropic volume change
    if (Isotrop == 1)
    {
      //Currently assumes that MM cutoffs are safe
      randnum = (((double)rand())/((double)RAND_MAX));
      double expan = 0.97+randnum*0.06;
      Lx *= expan;
      Ly *= expan;
      Lz *= expan;
    }
    //Scale positions
    #pragma omp parallel for
    for (int i=0;i<Natoms;i++)
    {
      for (int j=0;j<QMMMOpts.Nbeads;j++)
      {
        Struct2[i].P[j].x *= Lx/Lxtmp;
        Struct2[i].P[j].y *= Ly/Lytmp;
        Struct2[i].P[j].z *= Lz/Lztmp;
      }
    }
    #pragma omp barrier
    //Add PV energy term
    Enew += QMMMOpts.Press*Lx*Ly*Lz*atm2eV;
  }
  Enew += Get_PI_Epot(Struct2,QMMMOpts);
  Enew += Get_PI_Espring(Struct2,QMMMOpts);
  //Accept or reject
  double dE = Enew-Eold;
  if (QMMMOpts.Ensemble == "NPT")
  {
    //Flip sign on PV terms
    //dE -= 2*QMMMOpts.Press*Lx*Ly*Lz*atm2eV; //Flip Enew PV
    //dE += 2*QMMMOpts.Press*Lxtmp*Lytmp*Lztmp*atm2eV; //Flip Eold PV
    //Add log term
    dE -= Natoms*log((Lx*Ly*Lz)/(Lxtmp*Lytmp*Lztmp))/QMMMOpts.Beta;
  }
  double Prob = exp(-1*dE*QMMMOpts.Beta);
  randnum = (((double)rand())/((double)RAND_MAX));
  if ((dE <= 0) or (randnum < Prob))
  {
    //Accept
    Struct = Struct2;
    Emc = Enew;
    QMMMOpts.Eold = Enew;
    acc = 1;
  }
  else
  {
    //Reject
    Emc = Eold;
    //Revert to old box sizes
    Lx = Lxtmp;
    Ly = Lytmp;
    Lz = Lztmp;
  }
  //Return decision
  return acc;
};

