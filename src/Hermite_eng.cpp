/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM functions for calculating Hermite Gaussian integrals.

 References for integrals:
 Szabo and Ostlund, Modern Quantum Chemistry, (1989)
 Helgaker et al., Molecular Electronic-Structure Theory, (2000)

*/

//Definitions for the GauDen1s class
GauDen1s::GauDen1s(double magi,double widi,double qi,
                   double xi,double yi,double zi)
{
  //Save data given to the constructor
  mag = magi;
  wid = widi;
  q = qi;
  x = xi;
  y = yi;
  z = zi;
  //Convert to a.u.
  wid *= (BohrRad*BohrRad);
  return;
};

GauDen1s::~GauDen1s()
{
  //Generic destructor
  return;
};

double GauDen1s::ChrgNuc(double qpc, Coord pos, double Rcut)
{
  //Function to calculate interactions between nuclei and charges
  double E = 0;
  //Calculate distance
  double rij = 0.0;
  Coord tmppos; //Set positions in Angstroms for PBC check
  tmppos.x = x;
  tmppos.y = y;
  tmppos.z = z;
  rij = CoordDist2(tmppos,pos).vecMag(); //Squared distance in Angstroms
  //Check cutoff
  if (rij <= (Rcut*Rcut))
  {
    //Calculate Coulomb interaction energy
    rij = sqrt(rij)/BohrRad; //Switch to a.u.
    E = q*qpc/rij; //Energy in a.u.
  }
  //Change units
  E *= Har2eV;
  return E;
};

double GauDen1s::NucNuc(GauDen1s gau2, double Rcut)
{
  //Function to calculate interactions between nuclei
  double E = 0;
  //Calculate distance
  double rij = 0.0;
  Coord tmppos1; //Set positions in Angstroms for PBC check
  tmppos1.x = x;
  tmppos1.y = y;
  tmppos1.z = z;
  Coord tmppos2; //Set positions in Angstroms for PBC check
  tmppos2.x = gau2.x;
  tmppos2.y = gau2.y;
  tmppos2.z = gau2.z;
  rij = CoordDist2(tmppos1,tmppos2).vecMag(); //Squared distance in Angstroms
  //Check cutoff
  if (rij <= (Rcut*Rcut))
  {
    //Calculate Coulomb interaction energy
    rij = sqrt(rij)/BohrRad; //Switch to a.u.
    E = q*gau2.q/rij; //Energy in a.u.
  }
  //Change units
  E *= Har2eV;
  return E;
};

double GauDen1s::TwoOver(GauDen1s gau2)
{
  //Calculates the overlap of two Gaussian functions
  double Sij = 0.0;
  //Calculate new prefactor
  Sij = pi;
  Sij /= (wid+gau2.wid);
  Sij = pow(Sij,1.5);
  //Calculate distance
  double rij = 0.0;
  Coord tmppos1,tmppos2; //Set positions in Angstroms for PBC check
  tmppos1.x = x;
  tmppos1.y = y;
  tmppos1.z = z;
  tmppos2.x = gau2.x;
  tmppos2.y = gau2.y;
  tmppos2.z = gau2.z;
  rij = CoordDist2(tmppos1,tmppos2).vecMag(); //Squared distance in Angstroms
  rij = sqrt(rij)/BohrRad; //Change to a.u.
  //Calculate overlap
  rij *= -1*wid*gau2.wid*rij;
  rij /= (wid+gau2.wid);
  Sij *= exp(rij);
  //Add original prefactors
  Sij *= mag*gau2.mag;
  //Return overlap
  return Sij;
};

double GauDen1s::OneCoulPC(double qpc, Coord pos, double Rcut)
{
  //Calculates the electron-charge interactions for a point-charge
  double E = 0.0;
  //Calculate distance
  double rij = 0.0;
  Coord tmppos; //Set positions in Angstroms for PBC check
  tmppos.x = x;
  tmppos.y = y;
  tmppos.z = z;
  rij = CoordDist2(tmppos,pos).vecMag(); //Squared distance in Angstroms
  //Check cutoff
  if (rij <= (Rcut*Rcut))
  {
    //Calculate Coulomb interaction energy
    rij = sqrt(rij)/BohrRad; //Switch to a.u.
    E = erf(sqrt(wid)*rij);
    E *= -1*mag*qpc/rij; //Negative due to electron charge
  }
  //Change units
  E *= Har2eV;
  return E;
};

double GauDen1s::OneCoulNuc(GauDen1s gau2, double Rcut)
{
  //Calculates the electron-charge interactions for a point-charge
  double E = 0.0;
  //Calculate distance
  double rij = 0.0;
  Coord tmppos1,tmppos2; //Set positions in Angstroms for PBC check
  tmppos1.x = x;
  tmppos1.y = y;
  tmppos1.z = z;
  tmppos2.x = gau2.x;
  tmppos2.y = gau2.y;
  tmppos2.z = gau2.z;
  rij = CoordDist2(tmppos1,tmppos2).vecMag(); //Squared distance in Angstroms
  //Check cutoff
  if (rij <= (Rcut*Rcut))
  {
    //Calculate Coulomb interaction energy
    rij = sqrt(rij)/BohrRad; //Switch to a.u.
    E = erf(sqrt(wid)*rij);
    E *= -1*mag*gau2.q/rij; //Negative due to electron charge
  }
  //Change units
  E *= Har2eV;
  return E;
};

double GauDen1s::TwoCoul(GauDen1s gau2, double Rcut)
{
  //Calculates the electron-electron interactions
  double E = 0.0;
  //Calculate distance
  double rij = 0.0;
  Coord tmppos1,tmppos2; //Set positions in Angstroms for PBC check
  tmppos1.x = x;
  tmppos1.y = y;
  tmppos1.z = z;
  tmppos2.x = gau2.x;
  tmppos2.y = gau2.y;
  tmppos2.z = gau2.z;
  rij = CoordDist2(tmppos1,tmppos2).vecMag(); //Squared distance in Angstroms
  //Check cutoff
  if (rij <= (Rcut*Rcut))
  {
    //Calculate Coulomb interaction energy
    rij = sqrt(rij)/BohrRad; //Switch to a.u.
    //Calculate new gaussian exponent
    double a = 0.0;
    a = (wid*gau2.wid);
    a /= (wid+gau2.wid);
    //Calculate energy
    E = erf(sqrt(a)*rij);
    E *= mag*gau2.mag/rij; //Double negative
  }
  //Change units
  E *= Har2eV;
  return E;
};

//Definitions for the HermGau class
HermGau::HermGau(double ci, double a, int ix, int iy, int iz,
                 double xi, double yi, double zi)
{
  //Save data given to the constructor
  mag = ci;
  alpha = a;
  powx = ix;
  powy = iy;
  powz = iz;
  x = xi;
  y = yi;
  z = zi;
  return;
};

HermGau::~HermGau()
{
  //Generic destructor
  return;
};

double HermGau::Coeff()
{
  //Return the coefficient
  return mag;
};

double HermGau::XPos()
{
  //Return the x position
  return x;
};

double HermGau::YPos()
{
  //Return the y position
  return y;
};

double HermGau::ZPos()
{
  //Return the z position
  return z;
};

double HermGau::Alpha()
{
  //Return the Gaussian coefficient (width)
  return alpha;
};

int HermGau::XPow()
{
  //Return the Hermite power in the x direction
  return powx;
};

int HermGau::YPow()
{
  //Return the Hermite power in the y direction
  return powy;
};

int HermGau::ZPow()
{
  //Return the Hermite power in the z direction
  return powz;
};

double HermGau::Value(double xi, double yi, double zi)
{
  //Return the value at point (xi,yi,zi)
  double Val = 0; //Final value
  //Calculate distance
  double Xij = (xi-x)/BohrRad; //X distance (a.u.)
  double Yij = (yi-y)/BohrRad; //Y distance (a.u.)
  double Zij = (zi-z)/BohrRad; //Z distance (a.u.)
  //Scale Xij,Yij,Zij by alpha
  Xij *= alpha*Xij; //Alpha*Xij^2
  Yij *= alpha*Yij; //Alpha*Yij^2
  Zij *= alpha*Zij; //Alpha*Zij^2
  //Calculate the value of the basis function
  int Ni,signct;
  double Xval = 0; //Value of the X component
  Ni = ((int)floor(((double)powx)/2)); //Loop length
  signct = 1; //Flips the sign during the loop
  for (int i=0;i<Ni;i++)
  {
    //Calculate value for iteration i
    double valtmp; //Temporary storage
    valtmp = signct*pow(2*Xij,powx-(2*i));
    valtmp /= LICHEMFactorial(i)*LICHEMFactorial(powx-(2*i));
    //Update sum and sign
    Xval += valtmp;
    signct *= -1; //Change sign
  }
  Xval *= exp(-1*Xij); //Add Gaussian
  Xval *= LICHEMFactorial(powx);
  double Yval = 0; //Value of the Y component
  Ni = ((int)floor(((double)powy)/2)); //Loop length
  signct = 1; //Flips the sign during the loop
  for (int i=0;i<Ni;i++)
  {
    //Calculate value for iteration i
    double valtmp; //Temporary storage
    valtmp = signct*pow(2*Yij,powy-(2*i));
    valtmp /= LICHEMFactorial(i)*LICHEMFactorial(powy-(2*i));
    //Update sum and sign
    Yval += valtmp;
    signct *= -1; //Change sign
  }
  Yval *= exp(-1*Yij); //Add Gaussian
  Yval *= LICHEMFactorial(powy);
  double Zval = 0; //Value of the Z component
  Ni = ((int)floor(((double)powz)/2)); //Loop length
  signct = 1; //Flips the sign during the loop
  for (int i=0;i<Ni;i++)
  {
    //Calculate value for iteration i
    double valtmp; //Temporary storage
    valtmp = signct*pow(2*Zij,powz-(2*i));
    valtmp /= LICHEMFactorial(i)*LICHEMFactorial(powz-(2*i));
    //Update sum and sign
    Zval += valtmp;
    signct *= -1; //Change sign
  }
  Zval *= exp(-1*Zij); //Add Gaussian
  Zval *= LICHEMFactorial(powz);
  //Combine values from x,y,z
  Val = mag*Xval*Yval*Zval;
  return Val;
};

//Functions for calculating Gaussian integrals
double BoysFunc(int n, double x)
{
  //Recursive Boys function
  double Val = 0.0;
  if (n == 0)
  {
    //Zero order Boys function
    Val = sqrt(pi/(4*x))*erf(sqrt(x));
    return Val;
  }
  //Recursively calculate value
  Val = (((2*n-1)*BoysFunc(n-1,x))-exp(-1*x))/(2*x);
  return Val;
};

double HermCoul2e(HermGau& Gi, HermGau& Gj)
{
  //Recursive two electron Coulomb integral
  double Eij = 0; //Energy
  //Combine Gaussians with the Gaussian product rule
  double anew = Gi.Alpha()+Gj.Alpha(); //New Gaussian coefficient
  int powx = Gi.XPow()+Gj.XPow(); //New X power
  int powy = Gi.YPow()+Gj.YPow(); //New Y power
  int powz = Gi.ZPow()+Gj.ZPow(); //New Z power
  //Update magnitude based on the separation
  double mu = Gi.Alpha()*Gj.Alpha()/anew; //Smearing parameter
  Coord posi,posj; //Temporary storage for positions
  posi.x = Gi.XPos();
  posi.y = Gi.YPos();
  posi.z = Gi.ZPos();
  posj.x = Gj.XPos();
  posj.y = Gj.YPos();
  posj.z = Gj.ZPos();
  Coord Disp = CoordDist2(posi,posj); //Calculate distances
  double Xij = Disp.x/BohrRad; //X distance (a.u.)
  double Yij = Disp.y/BohrRad; //Y distance (a.u.)
  double Zij = Disp.z/BohrRad; //Z distance (a.u.)
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians (a.u.)
  double newmag = Gi.Coeff()*Gj.Coeff(); //Product of old coefficients
  newmag *= exp(-mu*Rij2); //Scale based on distance
  //Create product Gaussian
  HermGau Gij(newmag,anew,powx,powy,powz,Xij,Yij,Zij);
  //Calculate integrals
  double Ix = 0; //Integral in the x direction
  if (Gij.XPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> HermMags; //Magnitude of the Gaussians
    vector<int> HermOrders; //Order of the Hermite function
    vector<int> BoysOrders; //Order of the Boys function
    HermMags.push_back(1.0);
    HermOrders.push_back(Gij.XPow());
    BoysOrders.push_back(0);
    //Recursion
    bool ContRecurs = 1; //Keeps the while loop going
    while (ContRecurs)
    {
      //Stop recursion unless orders are greater than zero
      ContRecurs = 0;
      //Create temporary arrays
      vector<double> NewMags; //New Gaussian magnitudes
      vector<int> NewOrders; //New Hermite orders
      vector<int> NewBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<HermOrders.size();i++)
      {
        //Create new Hermites
        if (HermOrders[i] > 0)
        {
          ContRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = HermOrders[i]-2;
          coeffi *= HermMags[i];
          if ((HermOrders[i]-2) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-2);
            NewBoys.push_back(BoysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.XPos();
          coeffi *= HermMags[i];
          if ((HermOrders[i]-1) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-1);
            NewBoys.push_back(BoysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      HermMags = NewMags;
      HermOrders = NewOrders;
      BoysOrders = NewBoys;
    }
    //Calculate X integral
    for (unsigned int i=0;i<BoysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.Alpha()),BoysOrders[i]);
      Itmp *= BoysFunc(BoysOrders[i],(Gij.Alpha()*Rij2));
      Ix += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.Alpha()),0);
    Itmp *= BoysFunc(0,(Gij.Alpha()*Rij2));
    Ix += Itmp; //Update integral
  }
  double Iy = 0; //Integral in the y direction
  if (Gij.YPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> HermMags; //Magnitude of the Gaussians
    vector<int> HermOrders; //Order of the Hermite function
    vector<int> BoysOrders; //Order of the Boys function
    HermMags.push_back(1.0);
    HermOrders.push_back(Gij.YPow());
    BoysOrders.push_back(0);
    //Recursion
    bool ContRecurs = 1; //Keeps the while loop going
    while (ContRecurs)
    {
      //Stop recursion unless orders are greater than zero
      ContRecurs = 0;
      //Create temporary arrays
      vector<double> NewMags; //New Gaussian magnitudes
      vector<int> NewOrders; //New Hermite orders
      vector<int> NewBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<HermOrders.size();i++)
      {
        //Create new Hermites
        if (HermOrders[i] > 0)
        {
          ContRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = HermOrders[i]-2;
          coeffi *= HermMags[i];
          if ((HermOrders[i]-2) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-2);
            NewBoys.push_back(BoysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.YPos();
          coeffi *= HermMags[i];
          if ((HermOrders[i]-1) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-1);
            NewBoys.push_back(BoysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      HermMags = NewMags;
      HermOrders = NewOrders;
      BoysOrders = NewBoys;
    }
    //Calculate Y integral
    for (unsigned int i=0;i<BoysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.Alpha()),BoysOrders[i]);
      Itmp *= BoysFunc(BoysOrders[i],(Gij.Alpha()*Rij2));
      Iy += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.Alpha()),0);
    Itmp *= BoysFunc(0,(Gij.Alpha()*Rij2));
    Iy += Itmp; //Update integral
  }
  double Iz = 0; //Integral in the z direction
  if (Gij.ZPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> HermMags; //Magnitude of the Gaussians
    vector<int> HermOrders; //Order of the Hermite function
    vector<int> BoysOrders; //Order of the Boys function
    HermMags.push_back(1.0);
    HermOrders.push_back(Gij.ZPow());
    BoysOrders.push_back(0);
    //Recursion
    bool ContRecurs = 1; //Keeps the while loop going
    while (ContRecurs)
    {
      //Stop recursion unless orders are greater than zero
      ContRecurs = 0;
      //Create temporary arrays
      vector<double> NewMags; //New Gaussian magnitudes
      vector<int> NewOrders; //New Hermite orders
      vector<int> NewBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<HermOrders.size();i++)
      {
        //Create new Hermites
        if (HermOrders[i] > 0)
        {
          ContRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = HermOrders[i]-2;
          coeffi *= HermMags[i];
          if ((HermOrders[i]-2) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-2);
            NewBoys.push_back(BoysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.ZPos();
          coeffi *= HermMags[i];
          if ((HermOrders[i]-1) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-1);
            NewBoys.push_back(BoysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      HermMags = NewMags;
      HermOrders = NewOrders;
      BoysOrders = NewBoys;
    }
    //Calculate Z integral
    for (unsigned int i=0;i<BoysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.Alpha()),BoysOrders[i]);
      Itmp *= BoysFunc(BoysOrders[i],(Gij.Alpha()*Rij2));
      Iz += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.Alpha()),0);
    Itmp *= BoysFunc(0,(Gij.Alpha()*Rij2));
    Iz += Itmp; //Update integral
  }
  //Combine integrals
  Eij = Ix*Iy*Iz; //Combine the integrals
  Eij *= Gij.Coeff(); //Scale by magnitude
  //Change units and return
  Eij *= Har2eV;
  return Eij;
};

double HermCoul1e(HermGau& Gi, double qj, Coord& Posj)
{
  //Recursive one electron Coulomb integral
  double Eij = 0; //Energy
  //Create a temporary Gaussian
  double newmag = Gi.Coeff(); //Magnitude
  double anew = Gi.Alpha(); //Gaussian coefficient
  int powx = Gi.XPow(); //X power
  int powy = Gi.YPow(); //Y power
  int powz = Gi.ZPow(); //Z power
  Coord posi; //Temporary storage for positions
  posi.x = Gi.XPos();
  posi.y = Gi.YPos();
  posi.z = Gi.ZPos();
  Coord Disp = CoordDist2(posi,Posj); //Calculate distances
  double Xij = Disp.x/BohrRad; //X distance (a.u.)
  double Yij = Disp.y/BohrRad; //Y distance (a.u.)
  double Zij = Disp.z/BohrRad; //Z distance (a.u.)
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians (a.u.)
  HermGau Gij(newmag,anew,powx,powy,powz,Xij,Yij,Zij);
  //Calculate integrals
  double Ix = 0; //Integral in the x direction
  if (Gij.XPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> HermMags; //Magnitude of the Gaussians
    vector<int> HermOrders; //Order of the Hermite function
    vector<int> BoysOrders; //Order of the Boys function
    HermMags.push_back(1.0);
    HermOrders.push_back(Gij.XPow());
    BoysOrders.push_back(0);
    //Recursion
    bool ContRecurs = 1; //Keeps the while loop going
    while (ContRecurs)
    {
      //Stop recursion unless orders are greater than zero
      ContRecurs = 0;
      //Create temporary arrays
      vector<double> NewMags; //New Gaussian magnitudes
      vector<int> NewOrders; //New Hermite orders
      vector<int> NewBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<HermOrders.size();i++)
      {
        //Create new Hermites
        if (HermOrders[i] > 0)
        {
          ContRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = HermOrders[i]-2;
          coeffi *= HermMags[i];
          if ((HermOrders[i]-2) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-2);
            NewBoys.push_back(BoysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.XPos();
          coeffi *= HermMags[i];
          if ((HermOrders[i]-1) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-1);
            NewBoys.push_back(BoysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      HermMags = NewMags;
      HermOrders = NewOrders;
      BoysOrders = NewBoys;
    }
    //Calculate X integral
    for (unsigned int i=0;i<BoysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.Alpha()),BoysOrders[i]);
      Itmp *= BoysFunc(BoysOrders[i],(Gij.Alpha()*Rij2));
      Ix += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.Alpha()),0);
    Itmp *= BoysFunc(0,(Gij.Alpha()*Rij2));
    Ix += Itmp; //Update integral
  }
  double Iy = 0; //Integral in the y direction
  if (Gij.YPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> HermMags; //Magnitude of the Gaussians
    vector<int> HermOrders; //Order of the Hermite function
    vector<int> BoysOrders; //Order of the Boys function
    HermMags.push_back(1.0);
    HermOrders.push_back(Gij.YPow());
    BoysOrders.push_back(0);
    //Recursion
    bool ContRecurs = 1; //Keeps the while loop going
    while (ContRecurs)
    {
      //Stop recursion unless orders are greater than zero
      ContRecurs = 0;
      //Create temporary arrays
      vector<double> NewMags; //New Gaussian magnitudes
      vector<int> NewOrders; //New Hermite orders
      vector<int> NewBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<HermOrders.size();i++)
      {
        //Create new Hermites
        if (HermOrders[i] > 0)
        {
          ContRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = HermOrders[i]-2;
          coeffi *= HermMags[i];
          if ((HermOrders[i]-2) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-2);
            NewBoys.push_back(BoysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.YPos();
          coeffi *= HermMags[i];
          if ((HermOrders[i]-1) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-1);
            NewBoys.push_back(BoysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      HermMags = NewMags;
      HermOrders = NewOrders;
      BoysOrders = NewBoys;
    }
    //Calculate Y integral
    for (unsigned int i=0;i<BoysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.Alpha()),BoysOrders[i]);
      Itmp *= BoysFunc(BoysOrders[i],(Gij.Alpha()*Rij2));
      Iy += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.Alpha()),0);
    Itmp *= BoysFunc(0,(Gij.Alpha()*Rij2));
    Iy += Itmp; //Update integral
  }
  double Iz = 0; //Integral in the z direction
  if (Gij.ZPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> HermMags; //Magnitude of the Gaussians
    vector<int> HermOrders; //Order of the Hermite function
    vector<int> BoysOrders; //Order of the Boys function
    HermMags.push_back(1.0);
    HermOrders.push_back(Gij.ZPow());
    BoysOrders.push_back(0);
    //Recursion
    bool ContRecurs = 1; //Keeps the while loop going
    while (ContRecurs)
    {
      //Stop recursion unless orders are greater than zero
      ContRecurs = 0;
      //Create temporary arrays
      vector<double> NewMags; //New Gaussian magnitudes
      vector<int> NewOrders; //New Hermite orders
      vector<int> NewBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<HermOrders.size();i++)
      {
        //Create new Hermites
        if (HermOrders[i] > 0)
        {
          ContRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = HermOrders[i]-2;
          coeffi *= HermMags[i];
          if ((HermOrders[i]-2) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-2);
            NewBoys.push_back(BoysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.ZPos();
          coeffi *= HermMags[i];
          if ((HermOrders[i]-1) >= 0)
          {
            NewMags.push_back(coeffi);
            NewOrders.push_back(HermOrders[i]-1);
            NewBoys.push_back(BoysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      HermMags = NewMags;
      HermOrders = NewOrders;
      BoysOrders = NewBoys;
    }
    //Calculate Z integral
    for (unsigned int i=0;i<BoysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.Alpha()),BoysOrders[i]);
      Itmp *= BoysFunc(BoysOrders[i],(Gij.Alpha()*Rij2));
      Iz += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.Alpha()),0);
    Itmp *= BoysFunc(0,(Gij.Alpha()*Rij2));
    Iz += Itmp; //Update integral
  }
  //Combine integrals
  Eij = Ix*Iy*Iz; //Combine the integrals
  Eij *= Gij.Coeff()*qj; //Scale by magnitude
  //Change units and return
  Eij *= Har2eV;
  return Eij;
};

double HermOverlap(HermGau& Gi, HermGau& Gj)
{
  //Recursive two electron overlap integral
  double Sij = 0; //Overlap
  //Combine Gaussians with the Gaussian product rule
  double anew = Gi.Alpha()+Gj.Alpha(); //New Gaussian coefficient
  int powx = Gi.XPow()+Gj.XPow(); //New X power
  int powy = Gi.YPow()+Gj.YPow(); //New Y power
  int powz = Gi.ZPow()+Gj.ZPow(); //New Z power
  //Update magnitude based on the separation
  double mu = Gi.Alpha()*Gj.Alpha()/anew; //Smearing parameter
  Coord posi,posj; //Temporary storage for positions
  posi.x = Gi.XPos();
  posi.y = Gi.YPos();
  posi.z = Gi.ZPos();
  posj.x = Gj.XPos();
  posj.y = Gj.YPos();
  posj.z = Gj.ZPos();
  Coord Disp = CoordDist2(posi,posj); //Calculate distances
  double Xij = Disp.x/BohrRad; //X distance (a.u.)
  double Yij = Disp.y/BohrRad; //Y distance (a.u.)
  double Zij = Disp.z/BohrRad; //Z distance (a.u.)
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians (a.u.)
  double newmag = Gi.Coeff()*Gj.Coeff(); //Product of old coefficients
  newmag *= exp(-mu*Rij2); //Scale based on distance
  //Create product Gaussian
  HermGau Gij(newmag,anew,powx,powy,powz,Xij,Yij,Zij);
  //Calculate integral
  
  Sij *= Gij.Coeff(); //Scale by magnitude
  //Return overlap
  return Sij;
};

