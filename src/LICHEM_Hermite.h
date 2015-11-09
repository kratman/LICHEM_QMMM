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
#ifndef HERMITE_HEADERS
#define HERMITE_HEADERS

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
    GauDen1s(double magi,double widi,double qi,double xi,double yi,double zi)
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
    }
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
    HermGau(double ci, double a, int ix, int iy, int iz,
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
      //Convert to a.u.
      alpha *= (BohrRad*BohrRad);
      return;
    }
    double Coeff(); //Return the coefficient (magnitude)
    double XPos(); //Return the x position
    double YPos(); //Return the y position
    double ZPos(); //Return the z position
    double Alpha(); //Return the Gaussian width
    int XPow(); //Return the Hermite power in the x direction
    int YPow(); //Return the Hermite power in the y direction
    int ZPow(); //Return the Hermite power in the z direction
};

#endif

