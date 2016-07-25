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

//Definitions for the HermGau class
HermGau::HermGau(double ci, double a, int ix, int iy, int iz,
                 double xi, double yi, double zi)
{
  //Save data given to the constructor
  mag_ = ci;
  alpha_ = a;
  powX_ = ix;
  powY_ = iy;
  powZ_ = iz;
  x_ = xi;
  y_ = yi;
  z_ = zi;
  return;
};

HermGau::~HermGau()
{
  //Generic destructor
  return;
};

double HermGau::coeff()
{
  //Return the coefficient
  return mag_;
};

double HermGau::xPos()
{
  //Return the x position
  return x_;
};

double HermGau::yPos()
{
  //Return the y position
  return y_;
};

double HermGau::zPos()
{
  //Return the z position
  return z_;
};

double HermGau::getAlpha()
{
  //Return the Gaussian coefficient (width)
  return alpha_;
};

int HermGau::xPow()
{
  //Return the Hermite power in the x direction
  return powX_;
};

int HermGau::yPow()
{
  //Return the Hermite power in the y direction
  return powY_;
};

int HermGau::zPow()
{
  //Return the Hermite power in the z direction
  return powZ_;
};

double HermGau::value(double xi, double yi, double zi)
{
  //Return the value at point (xi,yi,zi)
  double val = 0; //Final value
  //Calculate distance
  double Xij = (xi-x_)/bohrRad; //X distance (a.u.)
  double Yij = (yi-y_)/bohrRad; //Y distance (a.u.)
  double Zij = (zi-z_)/bohrRad; //Z distance (a.u.)
  //Scale Xij,Yij,Zij by alpha
  Xij *= alpha_*Xij; //Alpha*Xij^2
  Yij *= alpha_*Yij; //Alpha*Yij^2
  Zij *= alpha_*Zij; //Alpha*Zij^2
  //Calculate the value of the basis function
  int Ni,signCt;
  double xVal = 0; //Value of the X component
  Ni = ((int)floor(((double)powX_)/2)); //Loop length
  signCt = 1; //Flips the sign during the loop
  for (int i=0;i<Ni;i++)
  {
    //Calculate value for iteration i
    double valTemp; //Temporary storage
    valTemp = signCt*pow(2*Xij,powX_-(2*i));
    valTemp /= LICHEMFactorial(i)*LICHEMFactorial(powX_-(2*i));
    //Update sum and sign
    xVal += valTemp;
    signCt *= -1; //Change sign
  }
  xVal *= exp(-1*Xij); //Add Gaussian
  xVal *= LICHEMFactorial(powX_);
  double yVal = 0; //Value of the Y component
  Ni = ((int)floor(((double)powY_)/2)); //Loop length
  signCt = 1; //Flips the sign during the loop
  for (int i=0;i<Ni;i++)
  {
    //Calculate value for iteration i
    double valTemp; //Temporary storage
    valTemp = signCt*pow(2*Yij,powY_-(2*i));
    valTemp /= LICHEMFactorial(i)*LICHEMFactorial(powY_-(2*i));
    //Update sum and sign
    yVal += valTemp;
    signCt *= -1; //Change sign
  }
  yVal *= exp(-1*Yij); //Add Gaussian
  yVal *= LICHEMFactorial(powY_);
  double zVal = 0; //Value of the Z component
  Ni = ((int)floor(((double)powZ_)/2)); //Loop length
  signCt = 1; //Flips the sign during the loop
  for (int i=0;i<Ni;i++)
  {
    //Calculate value for iteration i
    double valTemp; //Temporary storage
    valTemp = signCt*pow(2*Zij,powZ_-(2*i));
    valTemp /= LICHEMFactorial(i)*LICHEMFactorial(powZ_-(2*i));
    //Update sum and sign
    zVal += valTemp;
    signCt *= -1; //Change sign
  }
  zVal *= exp(-1*Zij); //Add Gaussian
  zVal *= LICHEMFactorial(powZ_);
  //Combine values from x,y,z
  val = mag_*xVal*yVal*zVal;
  return val;
};

//Functions for calculating Gaussian integrals
double BoysFunc(int n, double x)
{
  //Recursive Boys function
  double val = 0.0;
  if (n == 0)
  {
    //Zero order Boys function
    val = sqrt(pi/(4*x))*erf(sqrt(x));
    return val;
  }
  //Recursively calculate value
  val = (((2*n-1)*BoysFunc(n-1,x))-exp(-1*x))/(2*x);
  return val;
};

double HermCoul2e(HermGau& Gi, HermGau& Gj)
{
  //Recursive two electron Coulomb integral
  double Eij = 0; //Energy
  //Combine Gaussians with the Gaussian product rule
  double aNew = Gi.getAlpha()+Gj.getAlpha(); //New Gaussian coefficient
  int powX = Gi.xPow()+Gj.xPow(); //New X power
  int powY = Gi.yPow()+Gj.yPow(); //New Y power
  int powZ = Gi.zPow()+Gj.zPow(); //New Z power
  //Update magnitude based on the separation
  double mu = Gi.getAlpha()*Gj.getAlpha()/aNew; //Smearing parameter
  Coord posi,posj; //Temporary storage for positions
  posi.x = Gi.xPos();
  posi.y = Gi.yPos();
  posi.z = Gi.zPos();
  posj.x = Gj.xPos();
  posj.y = Gj.yPos();
  posj.z = Gj.zPos();
  Coord disp = CoordDist2(posi,posj); //Calculate distances
  double Xij = disp.x/bohrRad; //X distance (a.u.)
  double Yij = disp.y/bohrRad; //Y distance (a.u.)
  double Zij = disp.z/bohrRad; //Z distance (a.u.)
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians (a.u.)
  double newMag = Gi.coeff()*Gj.coeff(); //Product of old coefficients
  newMag *= exp(-mu*Rij2); //Scale based on distance
  //Create product Gaussian
  HermGau Gij(newMag,aNew,powX,powY,powZ,Xij,Yij,Zij);
  //Calculate integrals
  double Ix = 0; //Integral in the x direction
  if (Gij.xPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> hermMags; //Magnitude of the Gaussians
    vector<int> hermOrders; //Order of the Hermite function
    vector<int> boysOrders; //Order of the Boys function
    hermMags.push_back(1.0);
    hermOrders.push_back(Gij.xPow());
    boysOrders.push_back(0);
    //Recursion
    bool contRecurs = 1; //Keeps the while loop going
    while (contRecurs)
    {
      //Stop recursion unless orders are greater than zero
      contRecurs = 0;
      //Create temporary arrays
      vector<double> newMags; //New Gaussian magnitudes
      vector<int> newOrders; //New Hermite orders
      vector<int> newBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<hermOrders.size();i++)
      {
        //Create new Hermites
        if (hermOrders[i] > 0)
        {
          contRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = hermOrders[i]-2;
          coeffi *= hermMags[i];
          if ((hermOrders[i]-2) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-2);
            newBoys.push_back(boysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.xPos();
          coeffi *= hermMags[i];
          if ((hermOrders[i]-1) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-1);
            newBoys.push_back(boysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      hermMags = newMags;
      hermOrders = newOrders;
      boysOrders = newBoys;
    }
    //Calculate X integral
    for (unsigned int i=0;i<boysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.getAlpha()),boysOrders[i]);
      Itmp *= BoysFunc(boysOrders[i],(Gij.getAlpha()*Rij2));
      Ix += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.getAlpha()),0);
    Itmp *= BoysFunc(0,(Gij.getAlpha()*Rij2));
    Ix += Itmp; //Update integral
  }
  double Iy = 0; //Integral in the y direction
  if (Gij.yPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> hermMags; //Magnitude of the Gaussians
    vector<int> hermOrders; //Order of the Hermite function
    vector<int> boysOrders; //Order of the Boys function
    hermMags.push_back(1.0);
    hermOrders.push_back(Gij.yPow());
    boysOrders.push_back(0);
    //Recursion
    bool contRecurs = 1; //Keeps the while loop going
    while (contRecurs)
    {
      //Stop recursion unless orders are greater than zero
      contRecurs = 0;
      //Create temporary arrays
      vector<double> newMags; //New Gaussian magnitudes
      vector<int> newOrders; //New Hermite orders
      vector<int> newBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<hermOrders.size();i++)
      {
        //Create new Hermites
        if (hermOrders[i] > 0)
        {
          contRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = hermOrders[i]-2;
          coeffi *= hermMags[i];
          if ((hermOrders[i]-2) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-2);
            newBoys.push_back(boysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.yPos();
          coeffi *= hermMags[i];
          if ((hermOrders[i]-1) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-1);
            newBoys.push_back(boysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      hermMags = newMags;
      hermOrders = newOrders;
      boysOrders = newBoys;
    }
    //Calculate Y integral
    for (unsigned int i=0;i<boysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.getAlpha()),boysOrders[i]);
      Itmp *= BoysFunc(boysOrders[i],(Gij.getAlpha()*Rij2));
      Iy += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.getAlpha()),0);
    Itmp *= BoysFunc(0,(Gij.getAlpha()*Rij2));
    Iy += Itmp; //Update integral
  }
  double Iz = 0; //Integral in the z direction
  if (Gij.zPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> hermMags; //Magnitude of the Gaussians
    vector<int> hermOrders; //Order of the Hermite function
    vector<int> boysOrders; //Order of the Boys function
    hermMags.push_back(1.0);
    hermOrders.push_back(Gij.zPow());
    boysOrders.push_back(0);
    //Recursion
    bool contRecurs = 1; //Keeps the while loop going
    while (contRecurs)
    {
      //Stop recursion unless orders are greater than zero
      contRecurs = 0;
      //Create temporary arrays
      vector<double> newMags; //New Gaussian magnitudes
      vector<int> newOrders; //New Hermite orders
      vector<int> newBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<hermOrders.size();i++)
      {
        //Create new Hermites
        if (hermOrders[i] > 0)
        {
          contRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = hermOrders[i]-2;
          coeffi *= hermMags[i];
          if ((hermOrders[i]-2) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-2);
            newBoys.push_back(boysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.zPos();
          coeffi *= hermMags[i];
          if ((hermOrders[i]-1) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-1);
            newBoys.push_back(boysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      hermMags = newMags;
      hermOrders = newOrders;
      boysOrders = newBoys;
    }
    //Calculate Z integral
    for (unsigned int i=0;i<boysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.getAlpha()),boysOrders[i]);
      Itmp *= BoysFunc(boysOrders[i],(Gij.getAlpha()*Rij2));
      Iz += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.getAlpha()),0);
    Itmp *= BoysFunc(0,(Gij.getAlpha()*Rij2));
    Iz += Itmp; //Update integral
  }
  //Combine integrals
  Eij = Ix*Iy*Iz; //Combine the integrals
  Eij *= Gij.coeff(); //Scale by magnitude
  //Change units and return
  Eij *= har2eV;
  return Eij;
};

double HermCoul1e(HermGau& Gi, double qj, Coord& Posj)
{
  //Recursive one electron Coulomb integral
  double Eij = 0; //Energy
  //Create a temporary Gaussian
  double newMag = Gi.coeff(); //Magnitude
  double aNew = Gi.getAlpha(); //Gaussian coefficient
  int powX = Gi.xPow(); //X power
  int powY = Gi.yPow(); //Y power
  int powZ = Gi.zPow(); //Z power
  Coord posi; //Temporary storage for positions
  posi.x = Gi.xPos();
  posi.y = Gi.yPos();
  posi.z = Gi.zPos();
  Coord disp = CoordDist2(posi,Posj); //Calculate distances
  double Xij = disp.x/bohrRad; //X distance (a.u.)
  double Yij = disp.y/bohrRad; //Y distance (a.u.)
  double Zij = disp.z/bohrRad; //Z distance (a.u.)
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians (a.u.)
  HermGau Gij(newMag,aNew,powX,powY,powZ,Xij,Yij,Zij);
  //Calculate integrals
  double Ix = 0; //Integral in the x direction
  if (Gij.xPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> hermMags; //Magnitude of the Gaussians
    vector<int> hermOrders; //Order of the Hermite function
    vector<int> boysOrders; //Order of the Boys function
    hermMags.push_back(1.0);
    hermOrders.push_back(Gij.xPow());
    boysOrders.push_back(0);
    //Recursion
    bool contRecurs = 1; //Keeps the while loop going
    while (contRecurs)
    {
      //Stop recursion unless orders are greater than zero
      contRecurs = 0;
      //Create temporary arrays
      vector<double> newMags; //New Gaussian magnitudes
      vector<int> newOrders; //New Hermite orders
      vector<int> newBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<hermOrders.size();i++)
      {
        //Create new Hermites
        if (hermOrders[i] > 0)
        {
          contRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = hermOrders[i]-2;
          coeffi *= hermMags[i];
          if ((hermOrders[i]-2) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-2);
            newBoys.push_back(boysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.xPos();
          coeffi *= hermMags[i];
          if ((hermOrders[i]-1) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-1);
            newBoys.push_back(boysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      hermMags = newMags;
      hermOrders = newOrders;
      boysOrders = newBoys;
    }
    //Calculate X integral
    for (unsigned int i=0;i<boysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.getAlpha()),boysOrders[i]);
      Itmp *= BoysFunc(boysOrders[i],(Gij.getAlpha()*Rij2));
      Ix += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.getAlpha()),0);
    Itmp *= BoysFunc(0,(Gij.getAlpha()*Rij2));
    Ix += Itmp; //Update integral
  }
  double Iy = 0; //Integral in the y direction
  if (Gij.yPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> hermMags; //Magnitude of the Gaussians
    vector<int> hermOrders; //Order of the Hermite function
    vector<int> boysOrders; //Order of the Boys function
    hermMags.push_back(1.0);
    hermOrders.push_back(Gij.yPow());
    boysOrders.push_back(0);
    //Recursion
    bool contRecurs = 1; //Keeps the while loop going
    while (contRecurs)
    {
      //Stop recursion unless orders are greater than zero
      contRecurs = 0;
      //Create temporary arrays
      vector<double> newMags; //New Gaussian magnitudes
      vector<int> newOrders; //New Hermite orders
      vector<int> newBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<hermOrders.size();i++)
      {
        //Create new Hermites
        if (hermOrders[i] > 0)
        {
          contRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = hermOrders[i]-2;
          coeffi *= hermMags[i];
          if ((hermOrders[i]-2) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-2);
            newBoys.push_back(boysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.yPos();
          coeffi *= hermMags[i];
          if ((hermOrders[i]-1) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-1);
            newBoys.push_back(boysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      hermMags = newMags;
      hermOrders = newOrders;
      boysOrders = newBoys;
    }
    //Calculate Y integral
    for (unsigned int i=0;i<boysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.getAlpha()),boysOrders[i]);
      Itmp *= BoysFunc(boysOrders[i],(Gij.getAlpha()*Rij2));
      Iy += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.getAlpha()),0);
    Itmp *= BoysFunc(0,(Gij.getAlpha()*Rij2));
    Iy += Itmp; //Update integral
  }
  double Iz = 0; //Integral in the z direction
  if (Gij.zPow() > 0)
  {
    //Aspherical Hermite Gaussians
    vector<double> hermMags; //Magnitude of the Gaussians
    vector<int> hermOrders; //Order of the Hermite function
    vector<int> boysOrders; //Order of the Boys function
    hermMags.push_back(1.0);
    hermOrders.push_back(Gij.zPow());
    boysOrders.push_back(0);
    //Recursion
    bool contRecurs = 1; //Keeps the while loop going
    while (contRecurs)
    {
      //Stop recursion unless orders are greater than zero
      contRecurs = 0;
      //Create temporary arrays
      vector<double> newMags; //New Gaussian magnitudes
      vector<int> newOrders; //New Hermite orders
      vector<int> newBoys; //New Boys function orders
      //Loop over Hermite functions
      for (unsigned int i=0;i<hermOrders.size();i++)
      {
        //Create new Hermites
        if (hermOrders[i] > 0)
        {
          contRecurs = 1; //Continue recursion
          double coeffi; //Temp. coefficient storage
          //First Hermite
          coeffi = hermOrders[i]-2;
          coeffi *= hermMags[i];
          if ((hermOrders[i]-2) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-2);
            newBoys.push_back(boysOrders[i]+1);
          }
          //Second Hermite
          coeffi = Gij.zPos();
          coeffi *= hermMags[i];
          if ((hermOrders[i]-1) >= 0)
          {
            newMags.push_back(coeffi);
            newOrders.push_back(hermOrders[i]-1);
            newBoys.push_back(boysOrders[i]+1);
          }
        }
      }
      //Save new magnitudes and orders
      hermMags = newMags;
      hermOrders = newOrders;
      boysOrders = newBoys;
    }
    //Calculate Z integral
    for (unsigned int i=0;i<boysOrders.size();i++)
    {
      double Itmp; //Temp. storage for integrals
      Itmp = pow((-1*Gij.getAlpha()),boysOrders[i]);
      Itmp *= BoysFunc(boysOrders[i],(Gij.getAlpha()*Rij2));
      Iz += Itmp; //Update integral
    }
  }
  else
  {
    //Spherical Hermite Gaussian
    double Itmp; //Temp. storage for integrals
    Itmp = pow((-1*Gij.getAlpha()),0);
    Itmp *= BoysFunc(0,(Gij.getAlpha()*Rij2));
    Iz += Itmp; //Update integral
  }
  //Combine integrals
  Eij = Ix*Iy*Iz; //Combine the integrals
  Eij *= Gij.coeff()*qj; //Scale by magnitude
  //Change units and return
  Eij *= har2eV;
  return Eij;
};

double HermOverlap(HermGau& Gi, HermGau& Gj)
{
  //Recursive two electron overlap integral
  double Sij = 0; //Overlap
  //Combine Gaussians with the Gaussian product rule
  double aNew = Gi.getAlpha()+Gj.getAlpha(); //New Gaussian coefficient
  int powX = Gi.xPow()+Gj.xPow(); //New X power
  int powY = Gi.yPow()+Gj.yPow(); //New Y power
  int powZ = Gi.zPow()+Gj.zPow(); //New Z power
  //Update magnitude based on the separation
  double mu = Gi.getAlpha()*Gj.getAlpha()/aNew; //Smearing parameter
  Coord posi,posj; //Temporary storage for positions
  posi.x = Gi.xPos();
  posi.y = Gi.yPos();
  posi.z = Gi.zPos();
  posj.x = Gj.xPos();
  posj.y = Gj.yPos();
  posj.z = Gj.zPos();
  Coord disp = CoordDist2(posi,posj); //Calculate distances
  double Xij = disp.x/bohrRad; //X distance (a.u.)
  double Yij = disp.y/bohrRad; //Y distance (a.u.)
  double Zij = disp.z/bohrRad; //Z distance (a.u.)
  double Rij2 = Xij*Xij+Yij*Yij+Zij*Zij; //Distance between Gaussians (a.u.)
  double newMag = Gi.coeff()*Gj.coeff(); //Product of old coefficients
  newMag *= exp(-mu*Rij2); //Scale based on distance
  //Create product Gaussian
  HermGau Gij(newMag,aNew,powX,powY,powZ,Xij,Yij,Zij);
  //Calculate integral
  
  Sij *= Gij.coeff(); //Scale by magnitude
  //Return overlap
  return Sij;
};

