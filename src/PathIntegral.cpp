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
double SpringEnergy(double k, double r2)
{
  //General harmonic bond for PI rings
  double E = 0.5*k*r2;
  return E;
};

double Get_PI_Espring(vector<QMMMAtom>& parts, QMMMSettings& QMMMOpts)
{
  //Calculate total harmonic PI ring energy
  double E = 0.0;
  double w0 = 1/(QMMMOpts.Beta*hbar);
  w0 *= w0*ToeV;
  #pragma omp parallel for reduction(+:E)
  for (int i=0;i<Natoms;i++)
  {
    parts[i].Ep = 0.0;
    double w = w0*parts[i].m*QMMMOpts.Nbeads;
    for (int j=0;j<QMMMOpts.Nbeads;j++)
    {
      //Bead energy, one bond to avoid double counting
      int j2 = j-1;
      if (j2 == -1)
      {
        j2 = QMMMOpts.Nbeads-1; //Ring PBC
      }
      double dr2 = CoordDist2(parts[i].P[j],parts[i].P[j2]);
      parts[i].Ep += SpringEnergy(w,dr2);
    }
    E += parts[i].Ep;
  }
  #pragma omp barrier
  return E;
};

double Get_PI_Epot(vector<QMMMAtom>& parts, QMMMSettings& QMMMOpts)
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
    if (Gaussian == 1)
    {
      t_qm_start = (unsigned)time(0);
      Es += GaussianEnergy(parts,QMMMOpts,i);
      Times_qm += (unsigned)time(0)-t_qm_start;
    }
    if (PSI4 == 1)
    {
      t_qm_start = (unsigned)time(0);
      Es += PSIEnergy(parts,QMMMOpts,i);
      Times_qm += (unsigned)time(0)-t_qm_start;
      //Clean up annoying useless files
      GlobalSys = system("rm -f psi.*");
    }
    if (NWChem == 1)
    {
      t_qm_start = (unsigned)time(0);
      Es += NWChemEnergy(parts,QMMMOpts,i);
      Times_qm += (unsigned)time(0)-t_qm_start;
    }
    //Calculate MM energy
    if (TINKER == 1)
    {
      t_mm_start = (unsigned)time(0);
      Es += TINKEREnergy(parts,QMMMOpts,i);
      Times_mm += (unsigned)time(0)-t_mm_start;
    }
    if (AMBER == 1)
    {
      t_mm_start = (unsigned)time(0);
      Es += AMBEREnergy(parts,QMMMOpts,i);
      Times_mm += (unsigned)time(0)-t_mm_start;
    }
    if (LAMMPS == 1)
    {
      t_mm_start = (unsigned)time(0);
      Es += LAMMPSEnergy(parts,QMMMOpts,i);
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

bool MCMove(vector<QMMMAtom>& parts, QMMMSettings& QMMMOpts, double& Emc)
{
  bool acc = 0;
  //Copy parts
  vector<QMMMAtom> parts2;
  parts2 = parts;
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
      if (parts[p].Frozen == 0)
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
    #pragma omp parallel for
    for (int i=0;i<QMMMOpts.Nbeads;i++)
    {
      double xp = parts2[p].P[i].x+dx;
      double yp = parts2[p].P[i].y+dy;
      double zp = parts2[p].P[i].z+dz;
      if (PBCon)
      {
        bool check = 1;
        while (check)
        {
          check = 0;
          if (xp > Lx)
          {
            xp -= Lx;
            check = 1;
          }
          if (xp < 0.0)
          {
            xp += Lx;
            check = 1;
          }
          if (yp > Ly)
          {
            yp -= Ly;
            check = 1;
          }
          if (yp < 0.0)
          {
            yp += Ly;
            check = 1;
          }
          if (zp > Lz)
          {
            zp -= Lz;
            check = 1;
          }
          if (zp < 0.0)
          {
            zp += Lz;
            check = 1;
          }
        }
      }
      parts2[p].P[i].x = xp;
      parts2[p].P[i].y = yp;
      parts2[p].P[i].z = zp;
    }
    #pragma omp barrier
  }
  if (randnum < BeadProb)
  {
    //Move a single bead
    int p;
    bool FrozenAt = 1;
    while (FrozenAt)
    {
      //Make sure the atom is not frozen
      p = (rand()%Natoms);
      if (parts[p].Frozen == 0)
      {
        FrozenAt = 0;
      }
    }
    for (int i=0;i<QMMMOpts.Nbeads;i++)
    {
      double randx = (((double)rand())/((double)RAND_MAX));
      double randy = (((double)rand())/((double)RAND_MAX));
      double randz = (((double)rand())/((double)RAND_MAX));
      double dx = 2*(randx-0.5)*step;
      double dy = 2*(randy-0.5)*step;
      double dz = 2*(randz-0.5)*step;
      bool check = 1;
      parts2[p].P[i].x += dx;
      if (PBCon)
      {
        while (check)
        {
          check = 0;
          if (parts2[p].P[i].x > Lx)
          {
            parts2[p].P[i].x -= Lx;
            check = 1;
          }
          if (parts2[p].P[i].x < 0.0)
          {
            parts2[p].P[i].x += Lx;
            check = 1;
          }
        }
      }
      parts2[p].P[i].y += dy;
      if (PBCon)
      {
        check = 1;
        while (check)
        {
          check = 0;
          if (parts2[p].P[i].y > Ly)
          {
            parts2[p].P[i].y -= Ly;
            check = 1;
          }
          if (parts2[p].P[i].y < 0.0)
          {
            parts2[p].P[i].y += Ly;
            check = 1;
          }
        }
      }
      parts2[p].P[i].z += dz;
      if (PBCon)
      {
        check = 1;
        while (check)
        {
          check = 0;
          if (parts2[p].P[i].z > Lz)
          {
            parts2[p].P[i].z -= Lz;
            check = 1;
          }
          if (parts2[p].P[i].z < 0.0)
          {
            parts2[p].P[i].z += Lz;
            check = 1;
          }
        }
      }
    }
  }
  //Calculate energies
  double Eold = QMMMOpts.Eold;
  double Enew = 0;
  double Lxtmp = Lx;
  double Lytmp = Ly;
  double Lztmp = Lz;
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < VolProb)
  {
    if (Isotrop == 0)
    {
      randnum = (((double)rand())/((double)RAND_MAX));
      double Lmin = 0.90*Lx;
      double Lmax = 1.10*Lx;
      Lx = Lmin+randnum*(Lmax-Lmin);
      randnum = (((double)rand())/((double)RAND_MAX));
      Lmin = 0.90*Ly;
      Lmax = 1.10*Ly;
      Ly = Lmin+randnum*(Lmax-Lmin);
      randnum = (((double)rand())/((double)RAND_MAX));
      Lmin = 0.90*Lz;
      Lmax = 1.10*Lz;
      Lz = Lmin+randnum*(Lmax-Lmin);
    }
    if (Isotrop == 1)
    {
      //Currently assumes that MM cutoffs are safe
      randnum = (((double)rand())/((double)RAND_MAX));
      double expan = 0.9+randnum*0.20;
      Lx *= expan;
      Ly *= expan;
      Lz *= expan;
    }
    #pragma omp parallel for
    for (int i=0;i<Natoms;i++)
    {
      for (int j=0;j<QMMMOpts.Nbeads;j++)
      {
        parts2[i].P[j].x *= Lx/Lxtmp;
        parts2[i].P[j].y *= Ly/Lytmp;
        parts2[i].P[j].z *= Lz/Lztmp;
      }
    }
    #pragma omp barrier
    Enew += QMMMOpts.Press*Lx*Ly*Lz*atm2eV;
  }
  Enew += Get_PI_Epot(parts2,QMMMOpts);
  Enew += Get_PI_Espring(parts2,QMMMOpts);
  //Accept or reject
  double dE = Enew-Eold;
  if (QMMMOpts.Ensemble == "NPT")
  {
    dE -= Natoms*log((Lx*Ly*Lz)/(Lxtmp*Lytmp*Lztmp))/QMMMOpts.Beta;
  }
  double Prob = exp(-1*dE*QMMMOpts.Beta);
  randnum = (((double)rand())/((double)RAND_MAX));
  if ((dE <= 0) or (randnum < Prob))
  {
    parts = parts2;
    Emc = Enew;
    QMMMOpts.Eold = Enew;
    acc = 1;
  }
  else
  {
    Emc = Eold;
    Lx = Lxtmp;
    Ly = Lytmp;
    Lz = Lztmp;
  }
  return acc;
};

