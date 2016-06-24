/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Functions definitions for the eFF semi-classical electron model.

 Reference for eFF:
 Su et al., Phys. Rev. Lett., 18, 99, 185002, (2007)

*/

//Kintetic and potential energy functions for eFF
double EFFEnergy(QMMMAtom& atom, QMMMElec& elec, int bead)
{
  //Atom-lepton interactions
  double E = 0.0;
  double r = CoordDist2(atom.P[bead],elec.P[bead]).vecMag();
  if (r <= elecCutoff*elecCutoff)
  {
    if (r == 0.0)
    {
      E = coul2eV*atom.MP[bead].q*elec.q*sqrt(8/pi)/elec.rad[bead];
    }
    else
    {
      r = sqrt(r);
      E = (coul2eV*atom.MP[bead].q*elec.q/r)*erf(sqrt2*r/elec.rad[bead]);
    }
  }
  return E;
};

double EFFCorr(QMMMElec& elec1, QMMMElec& elec2, int bead)
{
  //Lepton-lepton interactions
  double E = 0.0;
  double r = CoordDist2(elec1.P[bead],elec2.P[bead]).vecMag();
  //Electrostatic energy
  if (r <= (elecCutoff*elecCutoff))
  {
    r = sqrt(r);
    if (r == 0.0)
    {
      E = coul2eV*elec1.q*elec2.q*sqrt(8/pi);
      E /= sqrt(elec1.rad[bead]*elec1.rad[bead]
           +elec2.rad[bead]*elec2.rad[bead]);
      return hugeNum; //Escape to avoid singularities later
    }
    else
    {
      double radij = elec1.rad[bead]*elec1.rad[bead];
      radij += elec2.rad[bead]*elec2.rad[bead];
      radij = sqrt(radij);
      E = (coul2eV*elec1.q*elec2.q/r);
      E *= erf(sqrt2*r/radij);
    }
    //Pauli repulsion
    if (elec1.typ == elec2.typ)
    {
      //Overlap
      double Sij = 2/((elec1.rad[bead]/elec2.rad[bead])
             +(elec2.rad[bead]/elec1.rad[bead]));
      Sij *= Sij*Sij;
      Sij = sqrt(Sij);
      double tmp = -1*eFFrbar*eFFrbar*r*r;
      tmp /= (elec1.rad[bead]*elec1.rad[bead]+elec2.rad[bead]*elec2.rad[bead]);
      tmp /= eFFsbar*eFFsbar;
      Sij *= exp(tmp);
      //Kinetic energy difference
      double Tij = 1/(elec1.rad[bead]*elec1.rad[bead]);
      Tij += 1/(elec2.rad[bead]*elec2.rad[bead]);
      Tij *= 3/(2*eFFsbar*eFFsbar);
      tmp = 6*eFFsbar*eFFsbar*(elec1.rad[bead]*elec1.rad[bead]
            +elec2.rad[bead]*elec2.rad[bead]);
      tmp -= 4*eFFrbar*eFFrbar*r*r;
      tmp /= eFFsbar*eFFsbar*(elec1.rad[bead]*elec1.rad[bead]
             +elec2.rad[bead]*elec2.rad[bead]);
      tmp /= eFFsbar*eFFsbar*(elec1.rad[bead]*elec1.rad[bead]
             +elec2.rad[bead]*elec2.rad[bead]);
      Tij -= tmp;
      Tij *= har2eV*bohrRad*bohrRad;
      if (elec1.spin[bead] == elec2.spin[bead])
      {
        //Symmetric VB spin-orbital
        double Etmp = Sij*Sij/(1-(Sij*Sij));
        Etmp += (1-eFFRho)*Sij*Sij/(1+(Sij*Sij));
        E += Etmp*Tij;
      }
      else
      {
        //Antisymmetric VB spin orbital
        E += -1*eFFRho*Sij*Sij*Tij/(1+(Sij*Sij));
      }
    }
  }
  return E;
};

double KineticE_eFF(vector<QMMMElec>& elecs, QMMMSettings& QMMMOpts)
{
  //Total electron kinetic energy
  double E = 0;
  #pragma omp parallel for
  for (unsigned int i=0;i<elecs.size();i++)
  {
    elecs[i].Ep = 0.0;
    for (int k=0;k<QMMMOpts.NBeads;k++)
    {
      //Lepton kinetic energy
      double Etmp = 3/(2*elecs[i].rad[k]*elecs[i].rad[k]);
      Etmp *= elecMass/elecs[i].m; //Scale by mass
      Etmp *= har2eV*bohrRad*bohrRad;
      if (scale_eFF)
      {
        //Reduce kinetic energy as the beads increase
        Etmp /= elRtNBeads;
      }
      elecs[i].Ep += Etmp;
    }
  }
  for (unsigned int i=0;i<elecs.size();i++)
  {
    E += elecs[i].Ep;
  }
  return E;
};

double Get_EeFF(vector<QMMMAtom>& QMMMData, vector<QMMMElec>& elecs,
       QMMMSettings& QMMMOpts)
{
  //Total eFF interaction energy
  double E = 0;
  #pragma omp parallel
  {
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<Natoms;i++)
    {
      QMMMData[i].Ep = 0;
      for (unsigned int j=0;j<elecs.size();j++)
      {
        for (int k=0;k<QMMMOpts.NBeads;k++)
        {
          QMMMData[i].Ep += EFFEnergy(QMMMData[i],elecs[j],k);
        }
      }
    }
    #pragma omp for nowait schedule(dynamic)
    for (unsigned int i=0;i<elecs.size();i++)
    {
      elecs[i].Ep = 0.0;
      for (int k=0;k<QMMMOpts.NBeads;k++)
      {
        //Lepton interaction energy
        for (unsigned int j=0;j<i;j++)
        {
          //Lepton-lepton correlation
          elecs[i].Ep += EFFCorr(elecs[i],elecs[j],k);
        }
      }
    }
  }
  #pragma omp barrier
  for (int i=0;i<Natoms;i++)
  {
    E += QMMMData[i].Ep;
  }
  for (unsigned int i=0;i<elecs.size();i++)
  {
    E += elecs[i].Ep;
  }
  E /= QMMMOpts.NBeads; //Removes double counting
  return E;
};

