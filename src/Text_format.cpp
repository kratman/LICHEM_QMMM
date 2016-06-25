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
string LICHEMFormFloat(T inpVal, int wid)
{
  //Resizes a floating-point number to a set number of characters
  //NB: This was a product of my frustration with stream settings
  stringstream oldValue;
  string newValue;
  //Initialize settings
  oldValue.str("");
  oldValue << fixed;
  oldValue.precision(wid);
  //Save input value to the string
  oldValue << inpVal;
  newValue = oldValue.str();
  //Check for a decimal place
  bool hasDot = 0;
  for (unsigned int i=0;i<newValue.length();i++)
  {
    if (newValue[i] == '.')
    {
      hasDot = 1;
    }
  }
  if (!hasDot)
  {
    //Fix integers
    if (LICHEMCount(newValue) < (wid-2))
    {
      //Add a decimal point
      newValue += ".";
    }
    else
    {
      //Escape if the integer is too long
      return newValue;
    }
  }
  int NChars = newValue.length();
  //Resize string
  if (NChars > wid)
  {
    //Delete characters
    newValue.erase(newValue.begin()+wid,newValue.end());
  }
  else
  {
    //Pad with zeros
    int diff = wid-NChars;
    for (int i=0;i<diff;i++)
    {
      //Add a zero
      newValue += "0";
    }
  }
  return newValue;
};

//String formatting functions
template<typename T> int LICHEMCount(T origVal)
{
  //Function to count the number of characters in a string or int
  stringstream origText; //Stream to convert origval to a string
  origText << origVal;
  string intText = origText.str();
  //Count characters and return
  int NChar = intText.length();
  return NChar;
};

void LICHEMLowerText(string& origText)
{
  //Function to switch a string to lowercase letters
  stringstream newText;
  newText.str("");
  for (unsigned int i=0;i<origText.length();i++)
  {
    //Check and change case
    if (origText[i] == 'A')
    {
      newText << "a";
    }
    else if (origText[i] == 'B')
    {
      newText << "b";
    }
    else if (origText[i] == 'C')
    {
      newText << "c";
    }
    else if (origText[i] == 'D')
    {
      newText << "d";
    }
    else if (origText[i] == 'E')
    {
      newText << "e";
    }
    else if (origText[i] == 'F')
    {
      newText << "f";
    }
    else if (origText[i] == 'G')
    {
      newText << "g";
    }
    else if (origText[i] == 'H')
    {
      newText << "h";
    }
    else if (origText[i] == 'I')
    {
      newText << "i";
    }
    else if (origText[i] == 'J')
    {
      newText << "j";
    }
    else if (origText[i] == 'K')
    {
      newText << "k";
    }
    else if (origText[i] == 'L')
    {
      newText << "l";
    }
    else if (origText[i] == 'M')
    {
      newText << "m";
    }
    else if (origText[i] == 'N')
    {
      newText << "n";
    }
    else if (origText[i] == 'O')
    {
      newText << "o";
    }
    else if (origText[i] == 'P')
    {
      newText << "p";
    }
    else if (origText[i] == 'Q')
    {
      newText << "q";
    }
    else if (origText[i] == 'R')
    {
      newText << "r";
    }
    else if (origText[i] == 'S')
    {
      newText << "s";
    }
    else if (origText[i] == 'T')
    {
      newText << "t";
    }
    else if (origText[i] == 'U')
    {
      newText << "u";
    }
    else if (origText[i] == 'V')
    {
      newText << "v";
    }
    else if (origText[i] == 'W')
    {
      newText << "w";
    }
    else if (origText[i] == 'X')
    {
      newText << "x";
    }
    else if (origText[i] == 'Y')
    {
      newText << "y";
    }
    else if (origText[i] == 'Z')
    {
      newText << "z";
    }
    else
    {
      //Numbers, spaces, and other characters
      newText << origText[i];
    }
  }
  //Overwrite original string
  origText = newText.str();
  return;
};

void LICHEMUpperText(string& origText)
{
  //Function to switch a string to uppercase letters
  stringstream newText;
  newText.str("");
  for (unsigned int i=0;i<origText.length();i++)
  {
    //Check and change case
    if (origText[i] == 'a')
    {
      newText << "A";
    }
    else if (origText[i] == 'b')
    {
      newText << "B";
    }
    else if (origText[i] == 'c')
    {
      newText << "C";
    }
    else if (origText[i] == 'd')
    {
      newText << "D";
    }
    else if (origText[i] == 'e')
    {
      newText << "E";
    }
    else if (origText[i] == 'f')
    {
      newText << "F";
    }
    else if (origText[i] == 'g')
    {
      newText << "G";
    }
    else if (origText[i] == 'h')
    {
      newText << "H";
    }
    else if (origText[i] == 'i')
    {
      newText << "I";
    }
    else if (origText[i] == 'j')
    {
      newText << "J";
    }
    else if (origText[i] == 'k')
    {
      newText << "K";
    }
    else if (origText[i] == 'l')
    {
      newText << "L";
    }
    else if (origText[i] == 'm')
    {
      newText << "M";
    }
    else if (origText[i] == 'n')
    {
      newText << "N";
    }
    else if (origText[i] == 'o')
    {
      newText << "O";
    }
    else if (origText[i] == 'p')
    {
      newText << "P";
    }
    else if (origText[i] == 'q')
    {
      newText << "Q";
    }
    else if (origText[i] == 'r')
    {
      newText << "R";
    }
    else if (origText[i] == 's')
    {
      newText << "S";
    }
    else if (origText[i] == 't')
    {
      newText << "T";
    }
    else if (origText[i] == 'u')
    {
      newText << "U";
    }
    else if (origText[i] == 'v')
    {
      newText << "V";
    }
    else if (origText[i] == 'w')
    {
      newText << "W";
    }
    else if (origText[i] == 'x')
    {
      newText << "X";
    }
    else if (origText[i] == 'y')
    {
      newText << "Y";
    }
    else if (origText[i] == 'z')
    {
      newText << "Z";
    }
    else
    {
      //Numbers, spaces, and other characters
      newText << origText[i];
    }
  }
  //Overwrite original string
  origText = newText.str();
  return;
};

void LICHEMFixSciNot(string& origText)
{
  //Function change D scientific notation to E notation
  for (unsigned int i=0;i<origText.length();i++)
  {
    //Check and change case
    if ((origText[i] == 'D') or (origText[i] == 'd'))
    {
      origText[i] = 'e';
    }
    if (origText[i] == 'E')
    {
      origText[i] = 'e';
    }
  }
  return;
};

