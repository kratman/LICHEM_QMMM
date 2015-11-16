/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Compiled databases for looking up basis sets

 Reference for basis sets:
 

*/

//Basis set definitions
vector<HermGau> HermBasis(string Typ, string BasName)
{
  //Function to set specific Hermite basis sets
  bool BadBasis = 0; //Flag to exit if no basis set is found
  vector<HermGau> NewBasis;
  //NB: Organized by basis name, then element
  
  //Check for errors
  if (BadBasis)
  {
    cerr << "Error: Basis set " << BasName;
    cerr << " is not defined for atom " << Typ;
    cerr << "!!!" << '\n' << '\n';
    cerr.flush();
    exit(0);
  }
  //Return basis set if it was found in the database
  return NewBasis;
};

