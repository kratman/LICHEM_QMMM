/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM functions for Calculating Hermite Gaussian integrals.

 References for integrals:
 Szabo and Ostlund, Modern Quantum Chemistry, (1989)
 Helgaker et al., Molecular Electronic-Structure Theory, (2000)

*/

//Definitions for the HermGau class
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
  double Xij = Gi.XPos()-Gj.XPos(); //X distance
  double Yij = Gi.YPos()-Gj.YPos(); //Y distance
  double Zij = Gi.ZPos()-Gj.ZPos(); //Z distance
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians
  double newmag = Gi.Coeff()*Gj.Coeff(); //Product of old coefficients
  newmag *= exp(-mu*Rij2); //Scale based on distance
  //Create product Gaussian
  HermGau Gij(newmag,anew,powx,powy,powz,Xij,Yij,Zij);
  //Calculate integral
  
  Eij *= Gij.Coeff(); //Scale by magnitude
  //Change units and return
  Eij *= Har2eV;
  return Eij;
};

double HermCoul1e(HermGau& Gi, double qj, Coord& Pos)
{
  //Recursive one electron Coulomb integral
  double Eij = 0; //Energy
  //Create a temporary Gaussian
  double newmag = Gi.Coeff(); //Magnitude
  double anew = Gi.Alpha(); //Gaussian coefficient
  int powx = Gi.XPow(); //X power
  int powy = Gi.YPow(); //Y power
  int powz = Gi.ZPow(); //Z power
  double Xij = Gi.XPos()-Pos.x; //X distance
  double Yij = Gi.YPos()-Pos.y; //Y distance
  double Zij = Gi.ZPos()-Pos.z; //Z distance
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians
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
  double xcent = 0.5*(Gi.XPos()+Gj.XPos()); //Center of the Gaussians
  double ycent = 0.5*(Gi.YPos()+Gj.YPos()); //Center of the Gaussians
  double zcent = 0.5*(Gi.ZPos()+Gj.ZPos()); //Center of the Gaussians
  double anew = Gi.Alpha()+Gj.Alpha(); //New Gaussian coefficient
  int powx = Gi.XPow()+Gj.XPow(); //New X power
  int powy = Gi.YPow()+Gj.YPow(); //New Y power
  int powz = Gi.ZPow()+Gj.ZPow(); //New Z power
  //Update magnitude based on the separation
  double mu = Gi.Alpha()*Gj.Alpha()/anew; //Smearing parameter
  double Xij = Gi.XPos()-Gj.XPos(); //X distance
  double Yij = Gi.YPos()-Gj.YPos(); //Y distance
  double Zij = Gi.ZPos()-Gj.ZPos(); //Z distance
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians
  double newmag = Gi.Coeff()*Gj.Coeff(); //Product of old coefficients
  newmag *= exp(-mu*Rij2); //Scale based on distance
  //Create product Gaussian
  HermGau Gij(newmag,anew,powx,powy,powz,xcent,ycent,zcent);
  //Calculate integral
  
  Sij *= Gij.Coeff(); //Scale by magnitude
  //Return overlap
  return Sij;
};

