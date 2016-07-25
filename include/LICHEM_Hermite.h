/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM classes for calculating Hermite Gaussian integrals.

*/

//Make including safe
#ifndef LICHEM_HERMITE_HEADERS
#define LICHEM_HERMITE_HEADERS

//Gaussian basis set classes
class HermGau
{
  //Class for Hermite Gaussians
  private:
    //Gaussian properties
    double mag_; //Coefficient in front of the Gaussian
    double alpha_; //Exponent (width)
    int powX_; //Power in the x direction
    int powY_; //Power in the y direction
    int powZ_; //Power in the z direction
    //Position
    double x_; //X position in Angstroms
    double y_; //Y position in Angstroms
    double z_; //Z position in Angstroms
  public:
    //Constructor
    HermGau(double,double,int,int,int,double,double,double);
    //Destructor
    ~HermGau();
    //Functions to return private data
    double coeff(); //Return the coefficient (magnitude)
    double xPos(); //Return the x position
    double yPos(); //Return the y position
    double zPos(); //Return the z position
    double getAlpha(); //Return the Gaussian width
    int xPow(); //Return the Hermite power in the x direction
    int yPow(); //Return the Hermite power in the y direction
    int zPow(); //Return the Hermite power in the z direction
    //Functions for calculations
    double value(double,double,double); //Return the magnitude at (x,y,z)
};

#endif

