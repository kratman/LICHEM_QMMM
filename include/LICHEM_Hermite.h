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
class GauDen1s
{
  //Simple 1s Gaussian class
  private:
    //Properties
    double mag; //Magnitude/population (prefactor)
    double wid; //Width in a.u.
    double q; //Nuclear charge
    double x; //X position
    double y; //Y position
    double z; //Z position
  public:
    //Constructor
    GauDen1s(double,double,double,double,double,double);
    //Destructor
    ~GauDen1s();
    //Point-charge interactions
    double ChrgNuc(double,Coord,double); //Nuclei-charge electrostatic term
    double NucNuc(GauDen1s,double); //Nuclei-nuclei electrostatic term
    //Electron density integrals
    double TwoOver(GauDen1s); //Density-density overlap
    double OneCoulPC(double,Coord,double); //Density-charge (MM)
    double OneCoulNuc(GauDen1s,double); //Density-nucleus
    double TwoCoul(GauDen1s,double); //Density-density Coulomb repulsion
};

class HermGau
{
  //Class for Hermite Gaussians
  private:
    //Gaussian properties
    double mag; //Coefficient in front of the Gaussian
    double alpha; //Exponent (width)
    int powx; //Power in the x direction
    int powy; //Power in the y direction
    int powz; //Power in the z direction
    //Position
    double x; //X position in Angstroms
    double y; //Y position in Angstroms
    double z; //Z position in Angstroms
  public:
    //Constructor
    HermGau(double,double,int,int,int,double,double,double);
    //Destructor
    ~HermGau();
    //Functions to return private data
    double Coeff(); //Return the coefficient (magnitude)
    double XPos(); //Return the x position
    double YPos(); //Return the y position
    double ZPos(); //Return the z position
    double Alpha(); //Return the Gaussian width
    int XPow(); //Return the Hermite power in the x direction
    int YPow(); //Return the Hermite power in the y direction
    int ZPow(); //Return the Hermite power in the z direction
    //Functions for calculations
    double Value(double,double,double); //Return the magnitude at (x,y,z)
};

#endif

