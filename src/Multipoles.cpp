/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 LICHEM functions for manipulating multipoles.

 NB: Multipole routines are being held back until the LICHEM paper
 is published.

 Reference for TINKER and frame of reference rotations:
 

 References for conversion to point-charges:
 

*/

//TINKER routines
void ExtractTINKpoles(vector<QMMMAtom>& Struct, int Bead)
{
  //Parses TINKER parameter files to find multipoles and local frames
  
  return;
};

void RotateTINKCharges(vector<QMMMAtom>& Struct, int Bead)
{
  //Switches from the local frame of reference to the global frame
  //of reference
  
  return;
};

void WriteTINKMpole(vector<QMMMAtom>& Struct, fstream& ofile, int i, int Bead)
{
  //Write a new multipole definition for pseudo-atoms

  return;
};

//General routines
void WriteChargeFile(vector<QMMMAtom>& Struct, int Bead)
{
  //Function to write a file for the MM charges
  
  return;
};

void ExtractGlobalPoles(int& argc, char**& argv)
{
  //Function to print the multipoles in the global frame
  
  return;
};

RedMpole Cart2SphHarm(Mpole& poles)
{
  //Converts Cartesian multipoles to spherical harmonic multipoles
  RedMpole SHpoles; //Spherical harmonic multipoles

  return SHpoles;
};

OctCharges SphHarm2Charges(RedMpole poles)
{
  //Converts spherical harmonic multipoles to point-charges
  OctCharges PCgrid; //New point-charge multipoles

  return PCgrid;
};

