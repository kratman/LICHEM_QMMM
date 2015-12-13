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

double LICHEMFactorial(int n)
{
  //Calculate a factorial
  double val = 1;
  while (n > 0)
  {
    val *= n;
    n -= 1;
  }
  return val;
};

string LICHEMFormDouble(double InpVal, int wid)
{
  //Resizes a double to a set number of characters
  //NB: This was a product of my frustration with stream settings
  stringstream oldvalue;
  string newvalue;
  //Initialize settings
  oldvalue.str("");
  oldvalue << fixed;
  oldvalue.precision(wid);
  //Save input value to the string
  oldvalue << InpVal;
  newvalue = oldvalue.str();
  int Nchars = newvalue.length();
  //Resize string
  if (Nchars > wid)
  {
    //Delete characters
    newvalue.erase(newvalue.begin()+wid,newvalue.end());
  }
  else
  {
    //Pad with zeros
    int diff = wid-newvalue.length();
    for (int i=0;i<diff;i++)
    {
      //Add a zero
      newvalue += "0";
    }
  }
  return newvalue;
};

string LICHEMFormFloat(double InpVal, int wid)
{
  //Resizes a double to a set number of characters
  //NB: This was a product of my frustration with stream settings
  stringstream oldvalue;
  string newvalue;
  //Initialize settings
  oldvalue.str("");
  oldvalue << fixed;
  oldvalue.precision(wid);
  //Save input value to the string
  oldvalue << InpVal;
  newvalue = oldvalue.str();
  int Nchars = newvalue.length();
  //Resize string
  if (Nchars > wid)
  {
    //Delete characters
    newvalue.erase(newvalue.begin()+wid,newvalue.end());
  }
  else
  {
    //Pad with zeros
    int diff = wid-newvalue.length();
    for (int i=0;i<diff;i++)
    {
      //Add a zero
      newvalue += "0";
    }
  }
  return newvalue;
};

bool CheckFile(const string& file)
{
  //Checks if a file exists
  struct stat buffer;
  if (stat(file.c_str(),&buffer) != -1)
  {
    //Yep...
    return 1;
  }
  //Nope...
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

Coord CoordDist2(Coord& a, Coord& b)
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
  //Save displacements
  Coord DispAB; //Distance between A and B
  DispAB.x = dx;
  DispAB.y = dy;
  DispAB.z = dz;
  return DispAB;
};

//Functions to check connectivity
vector<int> TraceBoundary(vector<QMMMAtom>& Struct, int AtID)
{
  //Function to find all boundary atoms connected to a pseudobond atom
  bool BondError = 0; //Checks if the molecular structure "breaks" the math
  vector<int> BoundAtoms; //Final list of atoms
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
  //Check for errors
  if (BondError)
  {
    cerr << "Error: Two pseudo-bonds are connected through boudary atoms!!!";
    cerr << '\n';
    cerr << " The connections prevent LICHEM from correctly updating";
    cerr << " the charges" << '\n' << '\n';
    cerr.flush();
    //Quit to avoid unphysical results
    exit(0);
  }
  //Return list if there are no errors
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

//Structure correction functions
void PBCCenter(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Move the system to the center of the simulation box
  double avgx = 0;
  double avgy = 0;
  double avgz = 0;
  #pragma omp parallel for schedule(dynamic) reduction(+:avgx,avgy,avgz)
  for (int i=0;i<Natoms;i++)
  {
    //Loop over all beads
    double centx = 0;
    double centy = 0;
    double centz = 0;
    for (int j=0;j<QMMMOpts.Nbeads;j++)
    {
      //Update local average postion
      centx += Struct[i].P[j].x;
      centy += Struct[i].P[j].y;
      centz += Struct[i].P[j].z;
    }
    //Upate full averages
    avgx += centx;
    avgy += centy;
    avgz += centz;
  }
  //Convert sums to averages
  avgx /= Natoms*QMMMOpts.Nbeads;
  avgy /= Natoms*QMMMOpts.Nbeads;
  avgz /= Natoms*QMMMOpts.Nbeads;
  //Move atoms to the center of the box
  #pragma omp parallel for schedule(dynamic)
  for (int i=0;i<Natoms;i++)
  {
    //Loop over all beads
    for (int j=0;j<QMMMOpts.Nbeads;j++)
    {
      //Move bead to the center
      Struct[i].P[j].x -= avgx;
      Struct[i].P[j].y -= avgy;
      Struct[i].P[j].z -= avgz;
    }
  }
  //Return with updated structure
  return;
};

void PrintLapin()
{
  //Print a nice picture
  stringstream lapin;
  lapin.str("");
  //Add obfuscated text
  lapin << '\n';
  lapin << "                   /^\\";
  lapin << "     /^\\";
  lapin << '\n';
  lapin << "                   |  \\";
  lapin << "   /  |";
  lapin << '\n';
  lapin << "                   |  |";
  lapin << "   |  |";
  lapin << '\n';
  lapin << "                   |";
  lapin << "  |___|  |";
  lapin << '\n';
  lapin << "                  /";
  lapin << "           \\";
  lapin << '\n';
  lapin << "                 /";
  lapin << "   &";
  lapin << "     &";
  lapin << "   \\";
  lapin << '\n';
  lapin << "                 |";
  lapin << "             |";
  lapin << '\n';
  lapin << "                 |";
  lapin << "     .";
  lapin << " .";
  lapin << "     |";
  lapin << '\n';
  lapin << "                  \\";
  lapin << "   \\___/";
  lapin << "   /";
  lapin << '\n';
  lapin << "                   \\_";
  lapin << "       _/";
  lapin << '\n';
  lapin << "                     |";
  lapin << "     |";
  lapin << '\n' << '\n';
  lapin << "         ______";
  lapin << "       ______";
  lapin << "       ______";
  lapin << '\n';
  lapin << "        /";
  lapin << "      \\";
  lapin << "     /";
  lapin << "      \\";
  lapin << "     /";
  lapin << "      \\";
  lapin << '\n';
  lapin << "       /";
  lapin << "        \\";
  lapin << "   /";
  lapin << "        \\";
  lapin << "   /";
  lapin << "        \\";
  lapin << '\n';
  lapin << "      |";
  lapin << "----------";
  lapin << "|";
  lapin << " |";
  lapin << "          |";
  lapin << " |";
  lapin << "**********";
  lapin << "|";
  lapin << '\n';
  lapin << "      [";
  lapin << "          ]";
  lapin << " [";
  lapin << "%%%%%%%%%%";
  lapin << "]";
  lapin << " [";
  lapin << "          ]";
  lapin << '\n';
  lapin << "      |";
  lapin << "----------";
  lapin << "|";
  lapin << " |";
  lapin << "          |";
  lapin << " |";
  lapin << "**********";
  lapin << "|";
  lapin << '\n';
  lapin << "       \\";
  lapin << "        /";
  lapin << "   \\";
  lapin << "        /";
  lapin << "   \\";
  lapin << "        /";
  lapin << '\n';
  lapin << "        \\";
  lapin << "______";
  lapin << "/";
  lapin << "     \\";
  lapin << "______";
  lapin << "/";
  lapin << "     \\";
  lapin << "______";
  lapin << "/";
  lapin << '\n' << '\n' << '\n';
  //Print text and return
  cout << lapin.str();
  cout.flush();
  return;
};

