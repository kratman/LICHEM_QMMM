/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE functions for manipulating multipoles.

 Citation for TINKER and frame of reference rotations:
 

 Citations for conversion to point charges:
 Stone, The Theory of Intermolecular Forces, 2013
 Devereux et al., J. Chem. Theory Comp., 10, 10, 4229, 2014

*/

//TINKER routines
void ExtractTINKpoles(vector<QMMMAtom>& Struct)
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

//General routines
RedMpole Cart2SphHarm(Mpole& poles)
{
  //Converts cartesian multipoles to spherical harmonic multipoles
  RedMpole SHpoles; //Spherical harmonic multipoles

  return SHpoles;
};

OctCharges SphHarm2Charges(RedMpole poles)
{
  //Converts spherical harmonic multipoles to point charges
  OctCharges PCgrid; //New point charge multipoles

  return PCgrid;
};
