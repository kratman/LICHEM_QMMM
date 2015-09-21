/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Primary functions for LICHEM.

*/

//Periodic Table functions
PeriodicTable::PeriodicTable()
{
  //Set atomic properties
  //Note: The widths were calculated for the all-electron density
  Typs.push_back("H");
  GauWids.push_back(0.33046422838582712); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("He");
  GauWids.push_back(0.44937774826454929); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Li");
  GauWids.push_back(0.24891493951151511); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Be");
  GauWids.push_back(0.28393617417828898); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("B");
  GauWids.push_back(0.31688766703304322); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("C");
  GauWids.push_back(0.34748552120700499); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("N");
  GauWids.push_back(0.38230385698895192); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("O");
  GauWids.push_back(0.40455289847810627); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("F");
  GauWids.push_back(0.43227871563688025); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ne");
  GauWids.push_back(0.45940792416326015); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Na");
  GauWids.push_back(0.28603607723272795); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mg");
  GauWids.push_back(0.29359655766949572); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Al");
  GauWids.push_back(0.30283227512360922); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Si");
  GauWids.push_back(0.31462038964458622); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("P");
  GauWids.push_back(0.33301416291877084); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("S");
  GauWids.push_back(0.35056486830985728); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cl");
  GauWids.push_back(0.37229368179032268); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ar");
  GauWids.push_back(0.39498740756144213); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("K");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ca");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sc");
  GauWids.push_back(0.29098058180517811); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ti");
  GauWids.push_back(0.29703979129006947); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("V");
  GauWids.push_back(0.30365511638085596); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cr");
  GauWids.push_back(0.30996489546441913); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mn");
  GauWids.push_back(0.32461806746013033); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Fe");
  GauWids.push_back(0.33074715374383107); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Co");
  GauWids.push_back(0.33482470199473979); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ni");
  GauWids.push_back(0.33983030212559556); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cu");
  GauWids.push_back(0.34432982070592233); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Zn");
  GauWids.push_back(0.34411081636918234); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ga");
  GauWids.push_back(0.33702898576385676); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ge");
  GauWids.push_back(0.33935669473845503); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("As");
  GauWids.push_back(0.34818738929821641); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Se");
  GauWids.push_back(0.35888462095130624); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Br");
  GauWids.push_back(0.37243540751493021); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Kr");
  GauWids.push_back(0.3880309254592651); //Method: PBE0/aug-cc-pVQZ
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Y");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Zr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Nb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mo");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tc");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ru");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rh");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ag");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("In");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sn");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Te");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("I");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Xe");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cs");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ba");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("La");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ce");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Nd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Eu");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Gd");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Dy");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ho");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Er");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Yb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Lu");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Hf");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ta");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("W");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Re");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Os");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ir");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pt");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Au");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Hg");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Tl");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pb");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Bi");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Po");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("At");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rn");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Fr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ra");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Ac");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Th");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pa");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("U");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Np");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Pu");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Am");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Bk");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Cf");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Es");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Fm");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Md");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("No");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Lr");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Rf");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Db");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Sg");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Bh");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Hs");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
  Typs.push_back("Mt");
  GauWids.push_back(0.3); //Default value
  CovRadii.push_back(1.0); //Default value
};

string PeriodicTable::Typing(int Z)
{
  //Function to convert nuclear charges to atom types
  return Typs[Z-1];
};

int PeriodicTable::RevTyping(string AtName)
{
  //Function to convert atom types to nuclear charges
  int Z = 0;
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

double PeriodicTable::GetGauWid(string AtName)
{
  //Function to find the 1s Gaussian width of an atom
  double gauwid = 0;
  #pragma omp parallel for num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      gauwid = GauWids[i];
    }
  }
  #pragma omp barrier
  return gauwid;
};

double PeriodicTable::GetCovRadius(string AtName)
{
  //Function to find the 1s Gaussian width of an atom
  double radius = 0;
  #pragma omp parallel for num_threads(Ncpus)
  for (unsigned int i=0;i<Typs.size();i++)
  {
    if (AtName == Typs[i])
    {
      //Save value
      radius = CovRadii[i];
    }
  }
  #pragma omp barrier
  return radius;
};

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
  cout << "#                 ";
  cout << "LICHEM: Layered Interacting CHEmical Models";
  cout << "                 #";
  cout << '\n';
  cout << "#                                      ";
  cout << "                                       #";
  cout << '\n';
  cout << "#                      ";
  cout << "Symbiotic Computational Chemistry";
  cout << "                      #";
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

int FindMaxThreads()
{
  //Function to count the number of allowed threads
  int ct = 0; //Generic counter
  #pragma omp parallel reduction(+:ct)
  ct += 1; //Add one for each thread
  #pragma omp barrier
  //Return total count
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
  if (PBCon)
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
    if (Struct[BondID].PBregion and (BondID != AtID))
    {
      //Two PBs are connected and this system will fail
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
          if (Struct[BondID].PBregion and (BondID != AtID))
          {
            //Two PBs are connected and this system will fail
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
    cerr << "Error: Two pseudo-bonds are connected through boudary atoms!!!";
    cerr << '\n';
    cerr << " The connections prevent LICHEM from correctly updating";
    cerr << " the charges" << '\n' << '\n';
    cerr.flush();
    exit(0);
  }
  return BoundAtoms;
};

bool Bonded(vector<QMMMAtom>& Struct, int Atom1, int Atom2)
{
  //Function to check if two atoms are 1-2 connected
  bool IsBound = 0;
  for (unsigned int i=0;i<Struct[Atom2].Bonds.size();i++)
  {
    if (Atom1 == Struct[Atom2].Bonds[i])
    {
      IsBound = 1;
    }
  }
  return IsBound;
};

bool Angled(vector<QMMMAtom>& Struct, int Atom1, int Atom3)
{
  //Function to check if two atoms are 1-3 connected
  bool IsBound = 0;
  for (unsigned int i=0;i<Struct[Atom1].Bonds.size();i++)
  {
    int Atom2 = Struct[Atom1].Bonds[i];
    for (unsigned int j=0;j<Struct[Atom2].Bonds.size();j++)
    {
      if (Atom3 == Struct[Atom2].Bonds[j])
      {
        IsBound = 1;
      }
    }
  }
  return IsBound;
};

bool Dihedraled(vector<QMMMAtom>& Struct, int Atom1, int Atom4)
{
  //Function to check if two atoms are 1-4 connected
  bool IsBound = 0;
  for (unsigned int i=0;i<Struct[Atom1].Bonds.size();i++)
  {
    int Atom2 = Struct[Atom1].Bonds[i];
    for (unsigned int j=0;j<Struct[Atom2].Bonds.size();j++)
    {
      int Atom3 = Struct[Atom2].Bonds[j];
      for (unsigned int k=0;k<Struct[Atom3].Bonds.size();k++)
      {
        if (Atom4 == Struct[Atom3].Bonds[k])
        {
          IsBound = 1;
        }
      }
    }
  }
  return IsBound;
};

