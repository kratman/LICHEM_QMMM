/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 FLUKE functions for manipulating multipoles.

*/

RedMpole CartMP2SphHarm(Mpole& poles)
{
  //Converts cartesian multipoles to spherical harmonic multipoles
  RedMpole SHpoles; //Spherical harmonic multipoles

  return SHpoles;
};

OctCharges SphHarmMP2Charges(RedMpole& poles)
{
  //Converts spherical harmonic multipoles to point charges
  OctCharges PCgrid; //New point charge multipoles

  return PCgrid;
};
