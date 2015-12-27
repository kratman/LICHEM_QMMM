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
  #pragma omp parallel for schedule(dynamic) reduction(+:E)
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
  #pragma omp parallel for schedule(dynamic) num_threads(MCThreads) \
          reduction(+:E,QMTime,MMTime)
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
      if (Struct2[p].Frozen == 0)
      {
        FrozenAt = 0;
      }
    }
    double randx = (((double)rand())/((double)RAND_MAX));
    double randy = (((double)rand())/((double)RAND_MAX));
    double randz = (((double)rand())/((double)RAND_MAX));
    double dx = 2*(randx-0.5)*pimcstep*CentRatio;
    double dy = 2*(randy-0.5)*pimcstep*CentRatio;
    double dz = 2*(randz-0.5)*pimcstep*CentRatio;
    //Update positions
    #pragma omp parallel
    {
      #pragma omp for nowait schedule(dynamic)
      for (int i=0;i<QMMMOpts.Nbeads;i++)
      {
        Struct2[p].P[i].x += dx;
      }
      #pragma omp for nowait schedule(dynamic)
      for (int i=0;i<QMMMOpts.Nbeads;i++)
      {
        Struct2[p].P[i].y += dy;
      }
      #pragma omp for nowait schedule(dynamic)
      for (int i=0;i<QMMMOpts.Nbeads;i++)
      {
        Struct2[p].P[i].z += dz;
      }
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
      if (Struct2[p].Frozen == 0)
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
      double dx = 2*(randx-0.5)*pimcstep;
      double dy = 2*(randy-0.5)*pimcstep;
      double dz = 2*(randz-0.5)*pimcstep;
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
      //Assumes that MM cutoffs are safe
      randnum = (((double)rand())/((double)RAND_MAX));
      Lx += 2*(randnum-0.5)*pimcstep;
      randnum = (((double)rand())/((double)RAND_MAX));
      Ly += 2*(randnum-0.5)*pimcstep;
      randnum = (((double)rand())/((double)RAND_MAX));
      Lz += 2*(randnum-0.5)*pimcstep;
    }
    //Isotropic volume change
    if (Isotrop == 1)
    {
      //Assumes that MM cutoffs are safe
      randnum = (((double)rand())/((double)RAND_MAX));
      Lx += 2*(randnum-0.5)*pimcstep;
      Ly += 2*(randnum-0.5)*pimcstep;
      Lz += 2*(randnum-0.5)*pimcstep;
    }
    //Decide how to scale the centroids
    bool ScaleRing = 0; //Shift the ring
    randnum = (((double)rand())/((double)RAND_MAX));
    if (randnum >= 0.5)
    {
      //Evenly scale the size of the ring
      ScaleRing = 1;
    }
    //Scale positions
    if (ScaleRing)
    {
      #pragma omp parallel
      {
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            double shift;
            shift = Struct2[i].P[j].x;
            //Check PBC without wrapping the molecules
            bool check = 1; //Continue the PBC checks
            while (check)
            {
              //Check the value
              check = 0;
              if (shift > Lx)
              {
                shift -= Lx;
                check = 1;
              }
              if (shift < 0)
              {
                shift += Lx;
                check = 1;
              }
            }
            shift = ((Lx/Lxtmp)-1)*shift;
            Struct2[i].P[j].x += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            double shift;
            shift = Struct2[i].P[j].y;
            //Check PBC without wrapping the molecules
            bool check = 1; //Continue the PBC checks
            while (check)
            {
              //Check the value
              check = 0;
              if (shift > Ly)
              {
                shift -= Ly;
                check = 1;
              }
              if (shift < 0)
              {
                shift += Ly;
                check = 1;
              }
            }
            shift = ((Ly/Lytmp)-1)*shift;
            Struct2[i].P[j].y += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            double shift;
            shift = Struct2[i].P[j].z;
            //Check PBC without wrapping the molecules
            bool check = 1; //Continue the PBC checks
            while (check)
            {
              //Check the value
              check = 0;
              if (shift > Lz)
              {
                shift -= Lz;
                check = 1;
              }
              if (shift < 0)
              {
                shift += Lz;
                check = 1;
              }
            }
            shift = ((Lz/Lztmp)-1)*shift;
            Struct2[i].P[j].z += shift;
          }
        }
      }
      #pragma omp barrier
    }
    else
    {
      #pragma omp parallel
      {
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          //Find centroids
          double shift = 0; //Change of position for the centroid
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            shift += Struct2[i].P[j].x; //Add to the position sum
          }
          shift /= QMMMOpts.Nbeads; //Average position
          //Check PBC without wrapping the molecules
          bool check = 1; //Continue the PBC checks
          while (check)
          {
            //Check the value
            check = 0;
            if (shift > Lx)
            {
              shift -= Lx;
              check = 1;
            }
            if (shift < 0)
            {
              shift += Lx;
              check = 1;
            }
          }
          //Calculate the change in position
          shift = ((Lx/Lxtmp)-1)*shift;
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            //Update the position
            Struct2[i].P[j].x += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          //Find centroids
          double shift = 0; //Change of position for the centroid
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            shift += Struct2[i].P[j].y; //Add to the position sum
          }
          shift /= QMMMOpts.Nbeads; //Average position
          //Check PBC without wrapping the molecules
          bool check = 1; //Continue the PBC checks
          while (check)
          {
            //Check the value
            check = 0;
            if (shift > Ly)
            {
              shift -= Ly;
              check = 1;
            }
            if (shift < 0)
            {
              shift += Ly;
              check = 1;
            }
          }
          //Calculate the change in position
          shift = ((Ly/Lytmp)-1)*shift;
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            //Update the position
            Struct2[i].P[j].y += shift;
          }
        }
        #pragma omp for nowait schedule(dynamic)
        for (int i=0;i<Natoms;i++)
        {
          //Find centroids
          double shift = 0; //Change of position for the centroid
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            shift += Struct2[i].P[j].z; //Add to the position sum
          }
          shift /= QMMMOpts.Nbeads; //Average position
          //Check PBC without wrapping the molecules
          bool check = 1; //Continue the PBC checks
          while (check)
          {
            //Check the value
            check = 0;
            if (shift > Lz)
            {
              shift -= Lz;
              check = 1;
            }
            if (shift < 0)
            {
              shift += Lz;
              check = 1;
            }
          }
          //Calculate the change in position
          shift = ((Lz/Lztmp)-1)*shift;
          for (int j=0;j<QMMMOpts.Nbeads;j++)
          {
            //Update the position
            Struct2[i].P[j].z += shift;
          }
        }
      }
      #pragma omp barrier
    }
    //Add PV energy term
    Enew += QMMMOpts.Press*Lx*Ly*Lz*atm2eV;
  }
  Enew += Get_PI_Epot(Struct2,QMMMOpts);
  Enew += Get_PI_Espring(Struct2,QMMMOpts);
  //Accept or reject
  double dE = QMMMOpts.Beta*(Enew-Eold);
  if (QMMMOpts.Ensemble == "NPT")
  {
    //Add log term
    dE -= Natoms*log((Lx*Ly*Lz)/(Lxtmp*Lytmp*Lztmp));
  }
  double Prob = exp(-1*dE);
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum <= Prob)
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

