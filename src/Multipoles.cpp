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

 Reference for TINKER and frame of reference rotations:
 Ponder, TINKER - Software Tools for Molecular Design

 References for conversion to point-charges:
 Stone, The Theory of Intermolecular Forces, (2013)
 Devereux et al., J. Chem. Theory Comp., 10, 10, 4229, (2014)

*/

//TINKER routines
void ExtractTINKpoles(vector<QMMMAtom>& Struct, int Bead)
{
  //Parses TINKER parameter files to find multipoles and local frames
  string dummy; //Generic string
  fstream inFile,outFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Write XYZ data
    outFile << setw(6) << (Struct[i].id+1);
    outFile << " ";
    outFile << setw(3) << Struct[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(Struct[i].P[Bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(Struct[i].P[Bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(Struct[i].P[Bead].z,16);
    outFile << " ";
    outFile << setw(4) << Struct[i].numTyp;
    for (unsigned int j=0;j<Struct[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (Struct[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Write poledit input
  call.str("");
  call << "LICHM_" << Bead << ".txt";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << "2" << '\n';
  outFile << "LICHM_" << Bead << ".xyz" << '\n';
  outFile << '\n';
  outFile.flush();
  outFile.close();
  //Run poledit
  call.str("");
  call << "poledit < LICHM_" << Bead << ".txt > LICHM_" << Bead << ".out";
  globalSys = system(call.str().c_str());
  //Extract multipole frames
  call.str("");
  call << "LICHM_" << Bead << ".out";
  inFile.open(call.str().c_str(),ios_base::in);
  while (!inFile.eof())
  {
    //Parse file line by line
    getline(inFile,dummy);
    stringstream line(dummy);
    line >> dummy;
    if (dummy == "Multipoles")
    {
      line >> dummy;
      //Second check
      if (dummy == "With")
      {
        line >> dummy;
        //Third check
        if (dummy == "Altered")
        {
          //Read multipoles and frames
          for (int i=0;i<Natoms;i++)
          {
            //Clear junk
            getline(inFile,dummy);
            getline(inFile,dummy);
            getline(inFile,dummy);
            //Check if poles exist
            getline(inFile,dummy);
            stringstream line(dummy);
            line >> dummy;
            if (dummy == "Local")
            {
              //Collect definition
              line >> dummy >> Struct[i].MP[Bead].type;
              line >> Struct[i].MP[Bead].atom1;
              line >> Struct[i].MP[Bead].atom2;
              line >> Struct[i].MP[Bead].atom3;
              //Correct numbering
              Struct[i].MP[Bead].atom1 -= 1;
              Struct[i].MP[Bead].atom2 -= 1;
              Struct[i].MP[Bead].atom3 -= 1;
              //Collect charge
              inFile >> dummy >> Struct[i].MP[Bead].q;
              //Collect static dipole
              inFile >> dummy;
              inFile >> Struct[i].MP[Bead].Dx;
              inFile >> Struct[i].MP[Bead].Dy;
              inFile >> Struct[i].MP[Bead].Dz;
              //Initialize induced dipole
              Struct[i].MP[Bead].IDx = 0;
              Struct[i].MP[Bead].IDy = 0;
              Struct[i].MP[Bead].IDz = 0;
              //Collect quadrupole
              inFile >> dummy;
              inFile >> Struct[i].MP[Bead].Qxx;
              inFile >> Struct[i].MP[Bead].Qxy;
              inFile >> Struct[i].MP[Bead].Qyy;
              inFile >> Struct[i].MP[Bead].Qxz;
              inFile >> Struct[i].MP[Bead].Qyz;
              inFile >> Struct[i].MP[Bead].Qzz;
            }
            else
            {
              //Initialize the "blank" multipole
              Struct[i].MP[Bead].type = "None";
              Struct[i].MP[Bead].Dx = 0;
              Struct[i].MP[Bead].Dy = 0;
              Struct[i].MP[Bead].Dz = 0;
              Struct[i].MP[Bead].IDx = 0;
              Struct[i].MP[Bead].IDy = 0;
              Struct[i].MP[Bead].IDz = 0;
              Struct[i].MP[Bead].Qxx = 0;
              Struct[i].MP[Bead].Qyy = 0;
              Struct[i].MP[Bead].Qzz = 0;
              Struct[i].MP[Bead].Qxy = 0;
              Struct[i].MP[Bead].Qxz = 0;
              Struct[i].MP[Bead].Qyz = 0;
            }
            //Clear more junk
            getline(inFile,dummy);
          }
        }
      }
    }
  }
  inFile.close();
  //Clean up files
  call.str("");
  call << "rm -f LICHM_" << Bead << ".txt LICHM_";
  call << Bead << ".key LICHM_" << Bead << ".xyz LICHM_";
  call << Bead << ".out";
  globalSys = system(call.str().c_str());
  return;
};

void RotateTINKCharges(vector<QMMMAtom>& Struct, int Bead)
{
  //Switches from the local frame of reference to the global frame
  //of reference
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (int i=0;i<Natoms;i++)
  {
    Vector3d vecX,vecY,vecZ; //Local frame vectors
    //Initialize vectors in the global frame
    vecX(0) = 1;
    vecX(1) = 0;
    vecX(2) = 0;
    vecY(0) = 0;
    vecY(1) = 1;
    vecY(2) = 0;
    vecZ(0) = 0;
    vecZ(1) = 0;
    vecZ(2) = 1;
    //Rotate the charges
    if (Struct[i].MMregion)
    {
      //Find current orientation
      double x,y,z;
      if (Struct[i].MP[Bead].type == "Bisector")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Fill in z vector
        vecZ += vecX;
        vecZ.normalize();
        //Find x vector by subtracting overlap
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "Z-then-X")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find x vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "Z-Bisect")
      {
        //Find first vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Find third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom3].P[Bead].z;
        vecY(0) = -1*x; //Correct the direction
        vecY(1) = -1*y; //Correct the direction
        vecY(2) = -1*z; //Correct the direction
        vecY.normalize();
        //Combine vectors
        vecX += vecY;
        vecX.normalize();
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "3-Fold")
      {
        //First vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom3].P[Bead].z;
        vecY(0) = -1*x; //Correct the direction
        vecY(1) = -1*y; //Correct the direction
        vecY(2) = -1*z; //Correct the direction
        vecY.normalize();
        //Combine vectors and normalize
        vecZ += vecX+vecY;
        vecZ.normalize();
        //Find second axis by subtracting overlap and normalizing
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "Z-Only")
      {
        //Primary vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Use a global axis for the second vector
        vecX(0) = 1.0;
        vecX(1) = 0.0;
        vecX(2) = 0.0;
        if (vecZ.dot(vecX) > 0.85)
        {
          //Switch to y axis if overlap is large
          vecX(0) = 0.0;
          vecX(1) = 1.0;
        }
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "None")
      {
        Struct[i].PC[Bead].q1 = 0;
        Struct[i].PC[Bead].q2 = 0;
        Struct[i].PC[Bead].q3 = 0;
        Struct[i].PC[Bead].q4 = 0;
        Struct[i].PC[Bead].q5 = 0;
        Struct[i].PC[Bead].q6 = 0;
      }
      //Fill in y vector
      vecY = vecX.cross(vecZ);
      vecY.normalize();
      //Rotate to the global frame
      Mpole newPoles;
      //Add monopoles
      newPoles.q = Struct[i].MP[Bead].q;
      //Rotate dipoles
      newPoles.Dx = 0; //X component
      newPoles.Dx += Struct[i].MP[Bead].Dx*vecX(0);
      newPoles.Dx += Struct[i].MP[Bead].Dy*vecY(0);
      newPoles.Dx += Struct[i].MP[Bead].Dz*vecZ(0);
      newPoles.Dy = 0; //Y component
      newPoles.Dy += Struct[i].MP[Bead].Dx*vecX(1);
      newPoles.Dy += Struct[i].MP[Bead].Dy*vecY(1);
      newPoles.Dy += Struct[i].MP[Bead].Dz*vecZ(1);
      newPoles.Dz = 0; //Z component
      newPoles.Dz += Struct[i].MP[Bead].Dx*vecX(2);
      newPoles.Dz += Struct[i].MP[Bead].Dy*vecY(2);
      newPoles.Dz += Struct[i].MP[Bead].Dz*vecZ(2);
      //Add induced dipoles (Already in global frame)
      newPoles.Dx += Struct[i].MP[Bead].IDx;
      newPoles.Dy += Struct[i].MP[Bead].IDy;
      newPoles.Dz += Struct[i].MP[Bead].IDz;
      newPoles.IDx = 0;
      newPoles.IDy = 0;
      newPoles.IDz = 0;
      //Rotate quadrupoles (This looks awful, but it works)
      //NB: This is a hard coded matrix rotation
      newPoles.Qxx = 0; //XX component
      newPoles.Qxx += vecX(0)*vecX(0)*Struct[i].MP[Bead].Qxx;
      newPoles.Qxx += vecX(0)*vecY(0)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxx += vecX(0)*vecZ(0)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxx += vecY(0)*vecX(0)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxx += vecY(0)*vecY(0)*Struct[i].MP[Bead].Qyy;
      newPoles.Qxx += vecY(0)*vecZ(0)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecX(0)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxx += vecZ(0)*vecY(0)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecZ(0)*Struct[i].MP[Bead].Qzz;
      newPoles.Qxy = 0; //XY component
      newPoles.Qxy += vecX(0)*vecX(1)*Struct[i].MP[Bead].Qxx;
      newPoles.Qxy += vecX(0)*vecY(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxy += vecX(0)*vecZ(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxy += vecY(0)*vecX(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxy += vecY(0)*vecY(1)*Struct[i].MP[Bead].Qyy;
      newPoles.Qxy += vecY(0)*vecZ(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecX(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxy += vecZ(0)*vecY(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecZ(1)*Struct[i].MP[Bead].Qzz;
      newPoles.Qxz = 0; //XZ component
      newPoles.Qxz += vecX(0)*vecX(2)*Struct[i].MP[Bead].Qxx;
      newPoles.Qxz += vecX(0)*vecY(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxz += vecX(0)*vecZ(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxz += vecY(0)*vecX(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxz += vecY(0)*vecY(2)*Struct[i].MP[Bead].Qyy;
      newPoles.Qxz += vecY(0)*vecZ(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecX(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxz += vecZ(0)*vecY(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecZ(2)*Struct[i].MP[Bead].Qzz;
      newPoles.Qyy = 0; //YY component
      newPoles.Qyy += vecX(1)*vecX(1)*Struct[i].MP[Bead].Qxx;
      newPoles.Qyy += vecX(1)*vecY(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyy += vecX(1)*vecZ(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyy += vecY(1)*vecX(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyy += vecY(1)*vecY(1)*Struct[i].MP[Bead].Qyy;
      newPoles.Qyy += vecY(1)*vecZ(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecX(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyy += vecZ(1)*vecY(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecZ(1)*Struct[i].MP[Bead].Qzz;
      newPoles.Qyz = 0; //YZ component
      newPoles.Qyz += vecX(1)*vecX(2)*Struct[i].MP[Bead].Qxx;
      newPoles.Qyz += vecX(1)*vecY(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyz += vecX(1)*vecZ(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyz += vecY(1)*vecX(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyz += vecY(1)*vecY(2)*Struct[i].MP[Bead].Qyy;
      newPoles.Qyz += vecY(1)*vecZ(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecX(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyz += vecZ(1)*vecY(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecZ(2)*Struct[i].MP[Bead].Qzz;
      newPoles.Qzz = 0; //ZZ component
      newPoles.Qzz += vecX(2)*vecX(2)*Struct[i].MP[Bead].Qxx;
      newPoles.Qzz += vecX(2)*vecY(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qzz += vecX(2)*vecZ(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qzz += vecY(2)*vecX(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qzz += vecY(2)*vecY(2)*Struct[i].MP[Bead].Qyy;
      newPoles.Qzz += vecY(2)*vecZ(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecX(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qzz += vecZ(2)*vecY(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecZ(2)*Struct[i].MP[Bead].Qzz;
      //switch to point-charges
      Struct[i].PC[Bead] = SphHarm2Charges(Cart2SphHarm(newPoles));
      //Translate charges to the atom's location in the global frame
      Struct[i].PC[Bead].x1 += Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y1 += Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z1 += Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x2 += Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y2 += Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z2 += Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x3 += Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y3 += Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z3 += Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x4 += Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y4 += Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z4 += Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x5 += Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y5 += Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z5 += Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x6 += Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y6 += Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z6 += Struct[i].P[Bead].z;
    }
    else
    {
      //Set other charges to zero
      //NB: This makes sure that variables are not undefined
      Struct[i].PC[Bead].q1 = 0;
      Struct[i].PC[Bead].q2 = 0;
      Struct[i].PC[Bead].q3 = 0;
      Struct[i].PC[Bead].q4 = 0;
      Struct[i].PC[Bead].q5 = 0;
      Struct[i].PC[Bead].q6 = 0;
      Struct[i].PC[Bead].x1 = Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y1 = Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z1 = Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x2 = Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y2 = Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z2 = Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x3 = Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y3 = Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z3 = Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x4 = Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y4 = Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z4 = Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x5 = Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y5 = Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z5 = Struct[i].P[Bead].z;
      Struct[i].PC[Bead].x6 = Struct[i].P[Bead].x;
      Struct[i].PC[Bead].y6 = Struct[i].P[Bead].y;
      Struct[i].PC[Bead].z6 = Struct[i].P[Bead].z;
    }
  }
  return;
};

void WriteTINKMpole(vector<QMMMAtom>& Struct, fstream& outFile, int i, int Bead)
{
  //Write a new multipole definition for pseudo-bonds and QM atoms
  outFile << "multipole -"; //Negative sign defines the frame with atom IDs
  outFile << (Struct[i].id+1) << " ";
  //Print frame
  if (Struct[i].MP[Bead].type == "Z-then-X")
  {
    outFile << (Struct[i].MP[Bead].atom1+1) << " ";
    outFile << (Struct[i].MP[Bead].atom2+1) << " ";
    outFile << (Struct[i].MP[Bead].atom3+1) << " ";
  }
  if (Struct[i].MP[Bead].type == "Bisector")
  {
    outFile << "-"; //Defines the bisectors
    outFile << (Struct[i].MP[Bead].atom1+1) << " ";
    outFile << (Struct[i].MP[Bead].atom2+1) << " ";
    outFile << (Struct[i].MP[Bead].atom3+1) << " ";
  }
  if (Struct[i].MP[Bead].type == "Z-Bisector")
  {
    outFile << (Struct[i].MP[Bead].atom1+1) << " ";
    outFile << "-"; //Defines the bisectors
    outFile << (Struct[i].MP[Bead].atom2+1) << " ";
    outFile << (Struct[i].MP[Bead].atom3+1) << " ";
  }
  if (Struct[i].MP[Bead].type == "Z-Only")
  {
    outFile << (Struct[i].MP[Bead].atom1+1) << " ";
    outFile << "0" << " ";
    outFile << "0" << " ";
  }
  if (Struct[i].MP[Bead].type == "None")
  {
    outFile << "0" << " ";
    outFile << "0" << " ";
    outFile << "0" << " ";
  }
  //Print only the point-charge
  outFile << Struct[i].MP[Bead].q << '\n';
  //Print dummy dipoles
  outFile << "        0.0 0.0 0.0" << '\n';
  //Print dummy quadrupoles
  outFile << "        0.0" << '\n';
  outFile << "        0.0 0.0" << '\n';
  outFile << "        0.0 0.0 0.0" << '\n';
  outFile << '\n';
  return;
};

//General routines
void WriteChargeFile(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                     int Bead)
{
  //Function to write a file for the MM charges
  stringstream call; //Generic stream
  fstream outFile; //Stream for the charge file
  bool firstCharge = 1; //Always write the first charge
  //Find the center of mass
  Coord QMCOM; //QM region center of mass
  if (PBCon or QMMMOpts.useLREC)
  {
    QMCOM = FindQMCOM(Struct,QMMMOpts,Bead);
  }
  //Initialize charges
  if (AMOEBA)
  {
    if (TINKER)
    {
      //Set up current multipoles
      RotateTINKCharges(Struct,Bead);
    }
  }
  //Write charge file
  if (Gaussian or NWChem)
  {
    //Save file
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    outFile.open(call.str().c_str(),ios_base::out);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        //Check PBC (minimum image convention)
        Coord distCent; //Distance from QM COM
        double xShft = 0;
        double yShft = 0;
        double zShft = 0;
        if (PBCon or QMMMOpts.useLREC)
        {
          //Initialize displacements
          double dx,dy,dz; //Starting displacements
          dx = Struct[i].P[Bead].x-QMCOM.x;
          dy = Struct[i].P[Bead].y-QMCOM.y;
          dz = Struct[i].P[Bead].z-QMCOM.z;
          distCent = CoordDist2(Struct[i].P[Bead],QMCOM);
          //Calculate the shift in positions
          //NB: Generally this work out to be +/- {Lx,Ly,Lz}
          if (PBCon)
          {
            xShft = distCent.x-dx;
            yShft = distCent.y-dy;
            zShft = distCent.z-dz;
          }
        }
        //Check for long-range corrections
        double scrq = 1;
        if (QMMMOpts.useLREC)
        {
          //Use the long-range correction
          double rcom = 0; //Distance from center of mass
          //Calculate the distance from the center of mass
          rcom = distCent.vecMag();
          if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
          {
            //Scale the charge
            rcom = sqrt(rcom);
            double scrqA,scrqB; //Temporary variables
            //Calculate temp. variables
            scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
            scrqB = -3*scrqA*scrqA;
            scrqA *= 2*scrqA*scrqA;
            //Combine temp. variables
            scrqA += scrqB+1;
            //Set the scale factor
            scrq -= pow(scrqA,QMMMOpts.LRECPow);
          }
          else
          {
            //Delete the charge
            scrq = 0;
          }
        }
        if ((scrq > 0) or firstCharge)
        {
          if (CHRG)
          {
            //Add charges
            firstCharge = 0; //Skips writing the remaining zeros
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].x+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].y+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].z+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16);
            outFile << '\n';
          }
          if (AMOEBA)
          {
            //Add multipoles
            firstCharge = 0; //Skips writing the remaining zeros
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x1+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y1+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z1+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x2+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y2+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z2+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x3+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y3+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z3+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x4+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y4+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z4+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x5+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y5+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z5+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x6+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y6+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z6+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16);
            outFile << '\n';
          }
        }
      }
    }
  }
  if (PSI4)
  {
    //Save file
    call.str("");
    call << "MMCharges_" << Bead << ".txt";
    outFile.open(call.str().c_str(),ios_base::out);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        //Check PBC (minimum image convention)
        Coord distCent; //Distance from QM COM
        double xShft = 0;
        double yShft = 0;
        double zShft = 0;
        if (PBCon or QMMMOpts.useLREC)
        {
          //Initialize displacements
          double dx,dy,dz; //Starting displacements
          dx = Struct[i].P[Bead].x-QMCOM.x;
          dy = Struct[i].P[Bead].y-QMCOM.y;
          dz = Struct[i].P[Bead].z-QMCOM.z;
          distCent = CoordDist2(Struct[i].P[Bead],QMCOM);
          //Calculate the shift in positions
          //NB: Generally this work out to be +/- {Lx,Ly,Lz}
          if (PBCon)
          {
            xShft = distCent.x-dx;
            yShft = distCent.y-dy;
            zShft = distCent.z-dz;
          }
        }
        //Check for long-range corrections
        double scrq = 1;
        if (QMMMOpts.useLREC)
        {
          //Use the long-range correction
          double rcom = 0; //Distance from center of mass
          //Calculate the distance from the center of mass
          rcom = distCent.vecMag();
          if (rcom <= (QMMMOpts.LRECCut*QMMMOpts.LRECCut))
          {
            //Scale the charge
            rcom = sqrt(rcom);
            double scrqA,scrqB; //Temporary variables
            //Calculate temp. variables
            scrqA = (QMMMOpts.LRECCut-rcom)/QMMMOpts.LRECCut;
            scrqB = -3*scrqA*scrqA;
            scrqA *= 2*scrqA*scrqA;
            //Combine temp. variables
            scrqA += scrqB+1;
            //Set the scale factor
            scrq -= pow(scrqA,QMMMOpts.LRECPow);
          }
          else
          {
            //Delete the charge
            scrq = 0;
          }
        }
        if (scrq > 0)
        {
          if (CHRG)
          {
            //Add charges
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].x+xShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].y+yShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].P[Bead].z+zShft,16);
            outFile << ")" << '\n';
          }
          if (AMOEBA)
          {
            //Add multipoles
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x1+xShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y1+yShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z1+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x2+xShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y2+yShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z2+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x3+xShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y3+yShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z3+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x4+xShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y4+yShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z4+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x5+xShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y5+yShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z5+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].x6+xShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].y6+yShft,16) << ",";
            outFile << LICHEMFormFloat(Struct[i].PC[Bead].z6+zShft,16);
            outFile << ")" << '\n';
          }
        }
      }
    }
  }
  //Write to files
  outFile.flush();
  outFile.close();
  //Return to the QM calculations
  return;
};

void ExtractGlobalPoles(int& argc, char**& argv)
{
  //Function to print the multipoles in the global frame
  cout << fixed;
  cout.precision(12);
  fstream xyzFile,connectFile,regionFile;
  vector<QMMMAtom> Struct;
  QMMMSettings QMMMOpts;
  stringstream call; //Stream for system calls and reading/writing files
  string dummy;
  for (int i=0;i<argc;i++)
  {
    //Read file names
    dummy = string(argv[i]);
    if (dummy == "-n")
    {
      Ncpus = atoi(argv[i+1]);
    }
    if (dummy == "-x")
    {
      xyzFilename = string(argv[i+1]);
      xyzFile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-c")
    {
      conFilename = string(argv[i+1]);
      connectFile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-r")
    {
      regFilename = string(argv[i+1]);
      regionFile.open(argv[i+1],ios_base::in);
    }
  }
  //Make sure input files can be read
  bool doQuit = 0;
  if (!xyzFile.good())
  {
    cout << "Error: Could not open xyz file.";
    cout << '\n';
    cout.flush();
    doQuit = 1;
  }
  if (!connectFile.good())
  {
    cout << "Error: Could not open connectivity file.";
    cout << '\n';
    cout.flush();
    doQuit = 1;
  }
  if (!regionFile.good())
  {
    cout << "Error: Could not open region file.";
    cout << '\n';
    cout.flush();
    doQuit = 1;
  }
  if (doQuit)
  {
    //Quit with an error
    exit(0);
  }
  ReadLICHEMInput(xyzFile,connectFile,regionFile,Struct,QMMMOpts);
  //Assume there is a single bead
  int Bead = 0;
  //Rotate multipoles
  if (TINKER)
  {
    #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
    for (int i=0;i<Natoms;i++)
    {
      Vector3d vecX,vecY,vecZ; //Local frame vectors
      //Initialize vectors in the global frame
      vecX(0) = 1;
      vecX(1) = 0;
      vecX(2) = 0;
      vecY(0) = 0;
      vecY(1) = 1;
      vecY(2) = 0;
      vecZ(0) = 0;
      vecZ(1) = 0;
      vecZ(2) = 1;
      //Find current orientation and rotate multipoles
      double x,y,z;
      if (Struct[i].MP[Bead].type == "Bisector")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Fill in z vector
        vecZ += vecX;
        vecZ.normalize();
        //Find x vector by subtracting overlap
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "Z-then-X")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find x vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "Z-Bisect")
      {
        //Find first vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Find third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom3].P[Bead].z;
        vecY(0) = -1*x; //Correct the direction
        vecY(1) = -1*y; //Correct the direction
        vecY(2) = -1*z; //Correct the direction
        vecY.normalize();
        //Combine vectors
        vecX += vecY;
        vecX.normalize();
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "3-Fold")
      {
        //First vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom2].P[Bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom3].P[Bead].z;
        vecY(0) = -1*x; //Correct the direction
        vecY(1) = -1*y; //Correct the direction
        vecY(2) = -1*z; //Correct the direction
        vecY.normalize();
        //Combine vectors and normalize
        vecZ += vecX+vecY;
        vecZ.normalize();
        //Find second axis by subtracting overlap and normalizing
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (Struct[i].MP[Bead].type == "Z-Only")
      {
        //Primary vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].atom1].P[Bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Use a global axis for the second vector
        vecX(0) = 1.0;
        vecX(1) = 0.0;
        vecX(2) = 0.0;
        if (vecZ.dot(vecX) > 0.85)
        {
          //Switch to y axis if overlap is large
          vecX(0) = 0.0;
          vecX(1) = 1.0;
        }
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      //Fill in y vector
      vecY = vecX.cross(vecZ);
      vecY.normalize();
      //Rotate to the global frame
      Mpole newPoles;
      //Add monopoles
      newPoles.q = Struct[i].MP[Bead].q;
      //Rotate dipoles
      newPoles.Dx = 0; //X component
      newPoles.Dx += Struct[i].MP[Bead].Dx*vecX(0);
      newPoles.Dx += Struct[i].MP[Bead].Dy*vecY(0);
      newPoles.Dx += Struct[i].MP[Bead].Dz*vecZ(0);
      newPoles.Dy = 0; //Y component
      newPoles.Dy += Struct[i].MP[Bead].Dx*vecX(1);
      newPoles.Dy += Struct[i].MP[Bead].Dy*vecY(1);
      newPoles.Dy += Struct[i].MP[Bead].Dz*vecZ(1);
      newPoles.Dz = 0; //Z component
      newPoles.Dz += Struct[i].MP[Bead].Dx*vecX(2);
      newPoles.Dz += Struct[i].MP[Bead].Dy*vecY(2);
      newPoles.Dz += Struct[i].MP[Bead].Dz*vecZ(2);
      //Add induced dipoles (Already in global frame)
      newPoles.Dx += Struct[i].MP[Bead].IDx;
      newPoles.Dy += Struct[i].MP[Bead].IDy;
      newPoles.Dz += Struct[i].MP[Bead].IDz;
      newPoles.IDx = 0;
      newPoles.IDy = 0;
      newPoles.IDz = 0;
      //Rotate quadrupoles (This looks awful, but it works)
      newPoles.Qxx = 0; //XX component
      newPoles.Qxx += vecX(0)*vecX(0)*Struct[i].MP[Bead].Qxx;
      newPoles.Qxx += vecX(0)*vecY(0)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxx += vecX(0)*vecZ(0)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxx += vecY(0)*vecX(0)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxx += vecY(0)*vecY(0)*Struct[i].MP[Bead].Qyy;
      newPoles.Qxx += vecY(0)*vecZ(0)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecX(0)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxx += vecZ(0)*vecY(0)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecZ(0)*Struct[i].MP[Bead].Qzz;
      newPoles.Qxy = 0; //XY component
      newPoles.Qxy += vecX(0)*vecX(1)*Struct[i].MP[Bead].Qxx;
      newPoles.Qxy += vecX(0)*vecY(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxy += vecX(0)*vecZ(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxy += vecY(0)*vecX(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxy += vecY(0)*vecY(1)*Struct[i].MP[Bead].Qyy;
      newPoles.Qxy += vecY(0)*vecZ(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecX(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxy += vecZ(0)*vecY(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecZ(1)*Struct[i].MP[Bead].Qzz;
      newPoles.Qxz = 0; //XZ component
      newPoles.Qxz += vecX(0)*vecX(2)*Struct[i].MP[Bead].Qxx;
      newPoles.Qxz += vecX(0)*vecY(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxz += vecX(0)*vecZ(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxz += vecY(0)*vecX(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qxz += vecY(0)*vecY(2)*Struct[i].MP[Bead].Qyy;
      newPoles.Qxz += vecY(0)*vecZ(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecX(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qxz += vecZ(0)*vecY(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecZ(2)*Struct[i].MP[Bead].Qzz;
      newPoles.Qyy = 0; //YY component
      newPoles.Qyy += vecX(1)*vecX(1)*Struct[i].MP[Bead].Qxx;
      newPoles.Qyy += vecX(1)*vecY(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyy += vecX(1)*vecZ(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyy += vecY(1)*vecX(1)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyy += vecY(1)*vecY(1)*Struct[i].MP[Bead].Qyy;
      newPoles.Qyy += vecY(1)*vecZ(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecX(1)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyy += vecZ(1)*vecY(1)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecZ(1)*Struct[i].MP[Bead].Qzz;
      newPoles.Qyz = 0; //YZ component
      newPoles.Qyz += vecX(1)*vecX(2)*Struct[i].MP[Bead].Qxx;
      newPoles.Qyz += vecX(1)*vecY(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyz += vecX(1)*vecZ(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyz += vecY(1)*vecX(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qyz += vecY(1)*vecY(2)*Struct[i].MP[Bead].Qyy;
      newPoles.Qyz += vecY(1)*vecZ(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecX(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qyz += vecZ(1)*vecY(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecZ(2)*Struct[i].MP[Bead].Qzz;
      newPoles.Qzz = 0; //ZZ component
      newPoles.Qzz += vecX(2)*vecX(2)*Struct[i].MP[Bead].Qxx;
      newPoles.Qzz += vecX(2)*vecY(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qzz += vecX(2)*vecZ(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qzz += vecY(2)*vecX(2)*Struct[i].MP[Bead].Qxy;
      newPoles.Qzz += vecY(2)*vecY(2)*Struct[i].MP[Bead].Qyy;
      newPoles.Qzz += vecY(2)*vecZ(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecX(2)*Struct[i].MP[Bead].Qxz;
      newPoles.Qzz += vecZ(2)*vecY(2)*Struct[i].MP[Bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecZ(2)*Struct[i].MP[Bead].Qzz;
      //Save the global frame multipoles
      Struct[i].MP[Bead] = newPoles;
    }
  }
  //Calculate the total charge
  double totChg = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus) \
          reduction(+:totChg)
  for (int i=0;i<Natoms;i++)
  {
    totChg += Struct[i].MP[Bead].q;
  }
  //Write multipoles to the screen
  cout << '\n';
  cout << "Global frame multipole moments (a.u.):";
  cout << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    cout << "Atom " << i << ": ";
    cout << Struct[i].MMTyp << '\n';
    cout << "Position (x,y,z):" << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << Struct[i].P[Bead].x << " ";
    cout << Struct[i].P[Bead].y << " ";
    cout << Struct[i].P[Bead].z << '\n';
    cout << "Charge and dipole (x,y,z):" << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << Struct[i].MP[Bead].q << " ";
    cout << Struct[i].MP[Bead].Dx << " ";
    cout << Struct[i].MP[Bead].Dy << " ";
    cout << Struct[i].MP[Bead].Dz << '\n';
    cout << "Quadrupole tensor (x,y,z):" << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << Struct[i].MP[Bead].Qxx << " ";
    cout << Struct[i].MP[Bead].Qxy << " ";
    cout << Struct[i].MP[Bead].Qxz << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << Struct[i].MP[Bead].Qxy << " ";
    cout << Struct[i].MP[Bead].Qyy << " ";
    cout << Struct[i].MP[Bead].Qyz << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << Struct[i].MP[Bead].Qxz << " ";
    cout << Struct[i].MP[Bead].Qyz << " ";
    cout << Struct[i].MP[Bead].Qzz << '\n';
    cout << '\n';
  }
  cout << "Total charge: " << totChg;
  cout << '\n' << '\n';
  //Quit
  exit(0);
  return;
};

RedMpole Cart2SphHarm(Mpole& pole)
{
  //Converts Cartesian multipoles to spherical harmonic multipoles
  RedMpole SHPole; //Spherical harmonic multipoles
  //Diagonalize quadrupole moment tensor
  Matrix3d QPole;
  QPole(0,0) = pole.Qxx;
  QPole(0,1) = pole.Qxy;
  QPole(0,2) = pole.Qxz;
  QPole(1,0) = pole.Qxy;
  QPole(1,1) = pole.Qyy;
  QPole(1,2) = pole.Qyz;
  QPole(2,0) = pole.Qxz;
  QPole(2,1) = pole.Qyz;
  QPole(2,2) = pole.Qzz;
  //Change out of a.u.
  QPole *= bohrRad*bohrRad; //NB: TINKER also divides by 3
  EigenSolver<Matrix3d> QTensor; //There might be a better method
  QTensor.compute(QPole);
  Vector3d SHTensor;
  SHTensor = QTensor.eigenvalues().real();
  //Save vector
  Matrix3d vec = QTensor.eigenvectors().real();
  SHPole.vecX(0) = vec(0,0);
  SHPole.vecX(1) = vec(1,0);
  SHPole.vecX(2) = vec(2,0);
  SHPole.vecY(0) = vec(0,1);
  SHPole.vecY(1) = vec(1,1);
  SHPole.vecY(2) = vec(2,1);
  SHPole.vecZ(0) = vec(0,2);
  SHPole.vecZ(1) = vec(1,2);
  SHPole.vecZ(2) = vec(2,2);
  //Normalize vectors (Probably not needed)
  SHPole.vecX.normalize();
  SHPole.vecY.normalize();
  SHPole.vecZ.normalize();
  //Convert to spherical harmonics and rotate dipoles
  SHPole.Q00 = pole.q;
  SHPole.Q11c = 0; //X component
  SHPole.Q11c += pole.Dx*SHPole.vecX(0);
  SHPole.Q11c += pole.Dy*SHPole.vecX(1);
  SHPole.Q11c += pole.Dz*SHPole.vecX(2);
  SHPole.Q11c *= bohrRad; //Change out of a.u.
  SHPole.Q11s = 0; //Y component
  SHPole.Q11s += pole.Dx*SHPole.vecY(0);
  SHPole.Q11s += pole.Dy*SHPole.vecY(1);
  SHPole.Q11s += pole.Dz*SHPole.vecY(2);
  SHPole.Q11s *= bohrRad; //Change out of a.u.
  SHPole.Q10 = 0; //Z component
  SHPole.Q10 += pole.Dx*SHPole.vecZ(0);
  SHPole.Q10 += pole.Dy*SHPole.vecZ(1);
  SHPole.Q10 += pole.Dz*SHPole.vecZ(2);
  SHPole.Q10 *= bohrRad; //Change out of a.u.
  SHPole.Q22c = (SHTensor(0)-SHTensor(1))/sqrt(3); //Diagonal Qxx-Qyy
  SHPole.Q20 = SHTensor(2); //Diagonal Qzz
  return SHPole;
};

OctCharges SphHarm2Charges(RedMpole pole)
{
  //Converts spherical harmonic multipoles to point-charges
  OctCharges PCGrid; //New point-charge multipoles
  double pd = 0.25*bohrRad; //Positive displacement of the charges
  double nd = -1*pd; //Negative of the displacement
  //Charge in the +x direction
  PCGrid.q1 = pole.Q00/6;
  PCGrid.q1 += pole.Q11c/(2*pd);
  PCGrid.q1 -= pole.Q20/(6*pd*pd);
  PCGrid.q1 += pole.Q22c/(2*sqrt(3)*pd*pd);
  PCGrid.x1 = pd*pole.vecX(0);
  PCGrid.y1 = pd*pole.vecX(1);
  PCGrid.z1 = pd*pole.vecX(2);
  //Charge in the +y direction
  PCGrid.q2 = pole.Q00/6;
  PCGrid.q2 += pole.Q11s/(2*pd);
  PCGrid.q2 -= pole.Q20/(6*pd*pd);
  PCGrid.q2 -= pole.Q22c/(2*sqrt(3)*pd*pd);
  PCGrid.x2 = pd*pole.vecY(0);
  PCGrid.y2 = pd*pole.vecY(1);
  PCGrid.z2 = pd*pole.vecY(2);
  //Charge in the +z direction
  PCGrid.q3 = pole.Q00/6;
  PCGrid.q3 += pole.Q10/(2*pd);
  PCGrid.q3 += pole.Q20/(3*pd*pd);
  PCGrid.x3 = pd*pole.vecZ(0);
  PCGrid.y3 = pd*pole.vecZ(1);
  PCGrid.z3 = pd*pole.vecZ(2);
  //Charge in the -x direction
  PCGrid.q4 = pole.Q00/6;
  PCGrid.q4 -= pole.Q11c/(2*pd);
  PCGrid.q4 -= pole.Q20/(6*pd*pd);
  PCGrid.q4 += pole.Q22c/(2*sqrt(3)*pd*pd);
  PCGrid.x4 = nd*pole.vecX(0);
  PCGrid.y4 = nd*pole.vecX(1);
  PCGrid.z4 = nd*pole.vecX(2);
  //Charge in the -y direction
  PCGrid.q5 = pole.Q00/6;
  PCGrid.q5 -= pole.Q11s/(2*pd);
  PCGrid.q5 -= pole.Q20/(6*pd*pd);
  PCGrid.q5 -= pole.Q22c/(2*sqrt(3)*pd*pd);
  PCGrid.x5 = nd*pole.vecY(0);
  PCGrid.y5 = nd*pole.vecY(1);
  PCGrid.z5 = nd*pole.vecY(2);
  //Charge in the -z direction
  PCGrid.q6 = pole.Q00/6;
  PCGrid.q6 -= pole.Q10/(2*pd);
  PCGrid.q6 += pole.Q20/(3*pd*pd);
  PCGrid.x6 = nd*pole.vecZ(0);
  PCGrid.y6 = nd*pole.vecZ(1);
  PCGrid.z6 = nd*pole.vecZ(2);
  //Return
  return PCGrid;
};

