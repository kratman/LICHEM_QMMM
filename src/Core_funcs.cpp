/*

##############################################################################
#                                                                            #
#              FLUKE: Fields Layered Under Kohn-sham Electrons               #
#                             By: Eric G. Kratz                              #
#                                                                            #
##############################################################################

 Primary functions and classes for FLUKE. This must be the first cpp
 file imported into main().

*/

//Core utility functions
void PrintFancyTitle()
{
  cout << '\n';
  cout << "#######################################";
  cout << "########################################";
  cout << '\n';
  cout << "#                                      ";
  cout << "                                       #";
  cout << '\n';
  cout << "#               ";
  cout << "FLUKE: Fields Layered Under Kohn-Sham Electrons";
  cout << "               #";
  cout << '\n';
  cout << "#                                      ";
  cout << "                                       #";
  cout << '\n';
  cout << "#######################################";
  cout << "########################################";
  cout << '\n' << '\n';
  cout.flush();
  return;
};

bool CheckFile(const string& file)
{
  //Checks if a file exists
  struct stat buffer;
  if (stat(file.c_str(),&buffer) != -1)
  {
    return 1;
  }
  return 0;
};

int Get_threads()
{
  //Function to count the number of allowed threads
  int ct = 0;
  #pragma omp parallel reduction(+:ct)
  ct += 1;
  #pragma omp barrier
  return ct;
};

double Bohring(double ri)
{
  //Convert ri (Bohr) to Angstroms
  double r = 0;
  r = ri/BohrRad;
  return r;
};

double CoordDist2(Coord& a, Coord& b)
{
  //Displacements
  double dx = a.x-b.x;
  double dy = a.y-b.y;
  double dz = a.z-b.z;
  //PBC
  if (PBCon == 1)
  {
    bool check = 1;
    while (check)
    {
      check = 0;
      if (abs(dx) > (0.5*Lx))
      {
        dx = Lx-abs(dx);
        check = 1;
      }
      if (abs(dy) > (0.5*Ly))
      {
        dy = Ly-abs(dy);
        check = 1;
      }
      if (abs(dz) > (0.5*Lz))
      {
        dz = Lz-abs(dz);
        check = 1;
      }
    }
  }
  //Squared radius
  double r2 = dx*dx+dy*dy+dz*dz;
  return r2;
};

string Typing(int Z)
{
  //Function to convert nuclear charges to atom types
  //Messy, but it is rarely used
  string AtTyp = "";
  vector<string> Typs; //List of atomic symbols
  Typs.push_back("H");
  Typs.push_back("He");
  Typs.push_back("Li");
  Typs.push_back("Be");
  Typs.push_back("B");
  Typs.push_back("C");
  Typs.push_back("N");
  Typs.push_back("O");
  Typs.push_back("F");
  Typs.push_back("Ne");
  Typs.push_back("Na");
  Typs.push_back("Mg");
  Typs.push_back("Al");
  Typs.push_back("Si");
  Typs.push_back("P");
  Typs.push_back("S");
  Typs.push_back("Cl");
  Typs.push_back("Ar");
  Typs.push_back("K");
  Typs.push_back("Ca");
  Typs.push_back("Sc");
  Typs.push_back("Ti");
  Typs.push_back("V");
  Typs.push_back("Cr");
  Typs.push_back("Mn");
  Typs.push_back("Fe");
  Typs.push_back("Co");
  Typs.push_back("Ni");
  Typs.push_back("Cu");
  Typs.push_back("Zn");
  Typs.push_back("Ga");
  Typs.push_back("Ge");
  Typs.push_back("As");
  Typs.push_back("Se");
  Typs.push_back("Br");
  Typs.push_back("Kr");
  Typs.push_back("Rb");
  Typs.push_back("Sr");
  Typs.push_back("Y");
  Typs.push_back("Zr");
  Typs.push_back("Nb");
  Typs.push_back("Mo");
  Typs.push_back("Tc");
  Typs.push_back("Ru");
  Typs.push_back("Rh");
  Typs.push_back("Pd");
  Typs.push_back("Ag");
  Typs.push_back("Cd");
  Typs.push_back("In");
  Typs.push_back("Sn");
  Typs.push_back("Sb");
  Typs.push_back("Te");
  Typs.push_back("I");
  Typs.push_back("Xe");
  //Find atom type
  AtTyp = Typs[Z-1];
  return AtTyp;
};

int RevTyping(string AtName)
{
  //Function to convert atom types to nuclear charges
  //Messy, but it is rarely used
  int Z = 0;
  vector<string> Typs; //List of atomic symbols
  Typs.push_back("H");
  Typs.push_back("He");
  Typs.push_back("Li");
  Typs.push_back("Be");
  Typs.push_back("B");
  Typs.push_back("C");
  Typs.push_back("N");
  Typs.push_back("O");
  Typs.push_back("F");
  Typs.push_back("Ne");
  Typs.push_back("Na");
  Typs.push_back("Mg");
  Typs.push_back("Al");
  Typs.push_back("Si");
  Typs.push_back("P");
  Typs.push_back("S");
  Typs.push_back("Cl");
  Typs.push_back("Ar");
  Typs.push_back("K");
  Typs.push_back("Ca");
  Typs.push_back("Sc");
  Typs.push_back("Ti");
  Typs.push_back("V");
  Typs.push_back("Cr");
  Typs.push_back("Mn");
  Typs.push_back("Fe");
  Typs.push_back("Co");
  Typs.push_back("Ni");
  Typs.push_back("Cu");
  Typs.push_back("Zn");
  Typs.push_back("Ga");
  Typs.push_back("Ge");
  Typs.push_back("As");
  Typs.push_back("Se");
  Typs.push_back("Br");
  Typs.push_back("Kr");
  Typs.push_back("Rb");
  Typs.push_back("Sr");
  Typs.push_back("Y");
  Typs.push_back("Zr");
  Typs.push_back("Nb");
  Typs.push_back("Mo");
  Typs.push_back("Tc");
  Typs.push_back("Ru");
  Typs.push_back("Rh");
  Typs.push_back("Pd");
  Typs.push_back("Ag");
  Typs.push_back("Cd");
  Typs.push_back("In");
  Typs.push_back("Sn");
  Typs.push_back("Sb");
  Typs.push_back("Te");
  Typs.push_back("I");
  Typs.push_back("Xe");
  #pragma omp parallel for num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      Z = i+1;
    }
  }
  #pragma omp barrier
  return Z;
};

vector<int> TraceBoundary(vector<QMMMAtom>& Struct, int AtID)
{
  bool BondError = 0;
  vector<int> BoundAtoms;
  //Add atoms bonded to atom "AtID"
  for (unsigned int i=0;i<Struct[AtID].Bonds.size();i++)
  {
    int BondID = Struct[AtID].Bonds[i];
    if (Struct[BondID].BAregion)
    {
      BoundAtoms.push_back(BondID);
    }
    if (Struct[BondID].PAregion and (BondID != AtID))
    {
      //Two PAs are connected and this system will fail
      BondError = 1;
    }
  }
  //Check find other boundary atoms bonded to the initial set
  bool MoreFound = 1; //More boundary atoms were found
  while (MoreFound and (!BondError))
  {
    MoreFound = 0; //Break the loop
    vector<int> tmp;
    for (unsigned int i=0;i<BoundAtoms.size();i++)
    {
      int BAID = BoundAtoms[i];
      for (unsigned int j=0;j<Struct[BAID].Bonds.size();j++)
      {
        int BondID = Struct[BAID].Bonds[j];
        //Check if it is on the list
        bool IsThere = 0;
        for (unsigned int k=0;k<BoundAtoms.size();k++)
        {
          if (BondID == BoundAtoms[k])
          {
            IsThere = 1;
          }
        }
        if (!IsThere)
        {
          if (Struct[BondID].BAregion)
          {
            MoreFound = 1; //Keep going
            tmp.push_back(BondID);
          }
          if (Struct[BondID].PAregion and (BondID != AtID))
          {
            //Two PAs are connected and this system will fail
            BondError = 1;
          }
        }
      }
    }
    //Add them to the list
    for (unsigned int i=0;i<tmp.size();i++)
    {
      bool IsThere = 0;
      for (unsigned int j=0;j<BoundAtoms.size();j++)
      {
        //Avoid adding the atom twice
        if (tmp[i] == BoundAtoms[j])
        {
          IsThere = 1;
        }
      }
      if (!IsThere)
      {
        BoundAtoms.push_back(tmp[i]);
      }
    }
  }
  if (BondError)
  {
    cerr << "Error: Two pseudoatoms are connected through boudary atoms!!!";
    cerr << '\n';
    cerr << " The connections prevent FLUKE from correctly updating";
    cerr << " the charges" << '\n' << '\n';
    cerr.flush();
    exit(0);
  }
  return BoundAtoms;
};

