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
  //Combine Gaussians
  
  //Calculate integral
  
  Eij *= Gi.Coeff()*Gj.Coeff(); //Scale by magnitude
  //Change units and return
  Eij *= Har2eV;
  return Eij;
};

double HermCoul1e(HermGau& Gi, double qj, Coord& Pos)
{
  //Recursive one electron Coulomb integral
  double Eij = 0; //Energy
  //Calculate integral
  
  Eij *= Gi.Coeff()*qj; //Scale by magnitude
  //Change units and return
  Eij *= Har2eV;
  return Eij;
};

double HermOverlap(HermGau& Gi, HermGau& Gj)
{
  //Recursive two electron overlap integral
  double Sij = 0; //Overlap
  //Combine Gaussians
  
  //Calculate integral
  
  Sij *= Gi.Coeff()*Gj.Coeff(); //Scale by magnitude
  //Return overlap
  return Sij;
};

