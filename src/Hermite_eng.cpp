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

