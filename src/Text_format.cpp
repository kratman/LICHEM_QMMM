/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Text formatting functions for LICHEM.

 NB: These functions were written to make LICHEM text formatting consistent
 while not requiring additional libraries.

*/

//Number formatting functions
template<typename T>
string LICHEMFormFloat(T InpVal, int wid)
{
  //Resizes a floating-point number to a set number of characters
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
  //Check for a decimal place
  bool HasDot = 0;
  for (unsigned int i=0;i<newvalue.length();i++)
  {
    if (newvalue[i] == '.')
    {
      HasDot = 1;
    }
  }
  if (!HasDot)
  {
    //Fix integers
    if (LICHEMCount(newvalue) < (wid-2))
    {
      //Add a decimal point
      newvalue += ".";
    }
    else
    {
      //Escape if the integer is too long
      return newvalue;
    }
  }
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
    int diff = wid-Nchars;
    for (int i=0;i<diff;i++)
    {
      //Add a zero
      newvalue += "0";
    }
  }
  return newvalue;
};

//String formatting functions
template<typename T> int LICHEMCount(T origval)
{
  //Function to count the number of characters in a string or int
  stringstream origtext; //Stream to convert origval to a string
  origtext << origval;
  string inttext = origtext.str();
  //Count characters and return
  int Nchar = inttext.length();
  return Nchar;
};

void LICHEMLowerText(string& origtext)
{
  //Function to switch a string to lowercase letters
  stringstream newtext;
  newtext.str("");
  for (unsigned int i=0;i<origtext.length();i++)
  {
    //Check and change case
    if (origtext[i] == 'A')
    {
      newtext << "a";
    }
    else if (origtext[i] == 'B')
    {
      newtext << "b";
    }
    else if (origtext[i] == 'C')
    {
      newtext << "c";
    }
    else if (origtext[i] == 'D')
    {
      newtext << "d";
    }
    else if (origtext[i] == 'E')
    {
      newtext << "e";
    }
    else if (origtext[i] == 'F')
    {
      newtext << "f";
    }
    else if (origtext[i] == 'G')
    {
      newtext << "g";
    }
    else if (origtext[i] == 'H')
    {
      newtext << "h";
    }
    else if (origtext[i] == 'I')
    {
      newtext << "i";
    }
    else if (origtext[i] == 'J')
    {
      newtext << "j";
    }
    else if (origtext[i] == 'K')
    {
      newtext << "k";
    }
    else if (origtext[i] == 'L')
    {
      newtext << "l";
    }
    else if (origtext[i] == 'M')
    {
      newtext << "m";
    }
    else if (origtext[i] == 'N')
    {
      newtext << "n";
    }
    else if (origtext[i] == 'O')
    {
      newtext << "o";
    }
    else if (origtext[i] == 'P')
    {
      newtext << "p";
    }
    else if (origtext[i] == 'Q')
    {
      newtext << "q";
    }
    else if (origtext[i] == 'R')
    {
      newtext << "r";
    }
    else if (origtext[i] == 'S')
    {
      newtext << "s";
    }
    else if (origtext[i] == 'T')
    {
      newtext << "t";
    }
    else if (origtext[i] == 'U')
    {
      newtext << "u";
    }
    else if (origtext[i] == 'V')
    {
      newtext << "v";
    }
    else if (origtext[i] == 'W')
    {
      newtext << "w";
    }
    else if (origtext[i] == 'X')
    {
      newtext << "x";
    }
    else if (origtext[i] == 'Y')
    {
      newtext << "y";
    }
    else if (origtext[i] == 'Z')
    {
      newtext << "z";
    }
    else
    {
      //Numbers, spaces, and other characters
      newtext << origtext[i];
    }
  }
  //Overwrite original string
  origtext = newtext.str();
  return;
};

void LICHEMUpperText(string& origtext)
{
  //Function to switch a string to uppercase letters
  stringstream newtext;
  newtext.str("");
  for (unsigned int i=0;i<origtext.length();i++)
  {
    //Check and change case
    if (origtext[i] == 'a')
    {
      newtext << "A";
    }
    else if (origtext[i] == 'b')
    {
      newtext << "B";
    }
    else if (origtext[i] == 'c')
    {
      newtext << "C";
    }
    else if (origtext[i] == 'd')
    {
      newtext << "D";
    }
    else if (origtext[i] == 'e')
    {
      newtext << "E";
    }
    else if (origtext[i] == 'f')
    {
      newtext << "F";
    }
    else if (origtext[i] == 'g')
    {
      newtext << "G";
    }
    else if (origtext[i] == 'h')
    {
      newtext << "H";
    }
    else if (origtext[i] == 'i')
    {
      newtext << "I";
    }
    else if (origtext[i] == 'j')
    {
      newtext << "J";
    }
    else if (origtext[i] == 'k')
    {
      newtext << "K";
    }
    else if (origtext[i] == 'l')
    {
      newtext << "L";
    }
    else if (origtext[i] == 'm')
    {
      newtext << "M";
    }
    else if (origtext[i] == 'n')
    {
      newtext << "N";
    }
    else if (origtext[i] == 'o')
    {
      newtext << "O";
    }
    else if (origtext[i] == 'p')
    {
      newtext << "P";
    }
    else if (origtext[i] == 'q')
    {
      newtext << "Q";
    }
    else if (origtext[i] == 'r')
    {
      newtext << "R";
    }
    else if (origtext[i] == 's')
    {
      newtext << "S";
    }
    else if (origtext[i] == 't')
    {
      newtext << "T";
    }
    else if (origtext[i] == 'u')
    {
      newtext << "U";
    }
    else if (origtext[i] == 'v')
    {
      newtext << "V";
    }
    else if (origtext[i] == 'w')
    {
      newtext << "W";
    }
    else if (origtext[i] == 'x')
    {
      newtext << "X";
    }
    else if (origtext[i] == 'y')
    {
      newtext << "Y";
    }
    else if (origtext[i] == 'z')
    {
      newtext << "Z";
    }
    else
    {
      //Numbers, spaces, and other characters
      newtext << origtext[i];
    }
  }
  //Overwrite original string
  origtext = newtext.str();
  return;
};

void LICHEMFixSciNot(string& origtext)
{
  //Function change D scientific notation to E notation
  for (unsigned int i=0;i<origtext.length();i++)
  {
    //Check and change case
    if ((origtext[i] == 'D') or (origtext[i] == 'd'))
    {
      origtext[i] = 'e';
    }
    if (origtext[i] == 'E')
    {
      origtext[i] = 'e';
    }
  }
  return;
};

