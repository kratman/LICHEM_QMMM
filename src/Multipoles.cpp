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
void ExtractTINKpoles(vector<QMMMAtom>& QMMMData, int bead)
{
  //Parses TINKER parameter files to find multipoles and local frames
  string dummy; //Generic string
  fstream inFile,outFile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << bead << ".xyz";
  outFile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  outFile << Natoms << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Write XYZ data
    outFile << setw(6) << (QMMMData[i].id+1);
    outFile << " ";
    outFile << setw(3) << QMMMData[i].MMTyp;
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].x,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].y,16);
    outFile << " ";
    outFile << LICHEMFormFloat(QMMMData[i].P[bead].z,16);
    outFile << " ";
    outFile << setw(4) << QMMMData[i].numTyp;
    for (unsigned int j=0;j<QMMMData[i].bonds.size();j++)
    {
      outFile << " "; //Avoids trailing spaces
      outFile << setw(6) << (QMMMData[i].bonds[j]+1);
    }
    outFile << '\n';
  }
  outFile.flush();
  outFile.close();
  //Write poledit input
  call.str("");
  call << "LICHM_" << bead << ".txt";
  outFile.open(call.str().c_str(),ios_base::out);
  outFile << "2" << '\n';
  outFile << "LICHM_" << bead << ".xyz" << '\n';
  outFile << '\n';
  outFile.flush();
  outFile.close();
  //Run poledit
  call.str("");
  call << "poledit < LICHM_" << bead << ".txt > LICHM_" << bead << ".out";
  globalSys = system(call.str().c_str());
  //Extract multipole frames
  call.str("");
  call << "LICHM_" << bead << ".out";
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
              line >> dummy >> QMMMData[i].MP[bead].type;
              line >> QMMMData[i].MP[bead].atom1;
              line >> QMMMData[i].MP[bead].atom2;
              line >> QMMMData[i].MP[bead].atom3;
              //Correct numbering
              QMMMData[i].MP[bead].atom1 -= 1;
              QMMMData[i].MP[bead].atom2 -= 1;
              QMMMData[i].MP[bead].atom3 -= 1;
              //Collect charge
              inFile >> dummy >> QMMMData[i].MP[bead].q;
              //Collect static dipole
              inFile >> dummy;
              inFile >> QMMMData[i].MP[bead].Dx;
              inFile >> QMMMData[i].MP[bead].Dy;
              inFile >> QMMMData[i].MP[bead].Dz;
              //Initialize induced dipole
              QMMMData[i].MP[bead].IDx = 0;
              QMMMData[i].MP[bead].IDy = 0;
              QMMMData[i].MP[bead].IDz = 0;
              //Collect quadrupole
              inFile >> dummy;
              inFile >> QMMMData[i].MP[bead].Qxx;
              inFile >> QMMMData[i].MP[bead].Qxy;
              inFile >> QMMMData[i].MP[bead].Qyy;
              inFile >> QMMMData[i].MP[bead].Qxz;
              inFile >> QMMMData[i].MP[bead].Qyz;
              inFile >> QMMMData[i].MP[bead].Qzz;
            }
            else
            {
              //Initialize the "blank" multipole
              QMMMData[i].MP[bead].type = "None";
              QMMMData[i].MP[bead].Dx = 0;
              QMMMData[i].MP[bead].Dy = 0;
              QMMMData[i].MP[bead].Dz = 0;
              QMMMData[i].MP[bead].IDx = 0;
              QMMMData[i].MP[bead].IDy = 0;
              QMMMData[i].MP[bead].IDz = 0;
              QMMMData[i].MP[bead].Qxx = 0;
              QMMMData[i].MP[bead].Qyy = 0;
              QMMMData[i].MP[bead].Qzz = 0;
              QMMMData[i].MP[bead].Qxy = 0;
              QMMMData[i].MP[bead].Qxz = 0;
              QMMMData[i].MP[bead].Qyz = 0;
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
  call << "rm -f LICHM_" << bead << ".txt LICHM_";
  call << bead << ".key LICHM_" << bead << ".xyz LICHM_";
  call << bead << ".out";
  globalSys = system(call.str().c_str());
  return;
};

void RotateTINKCharges(vector<QMMMAtom>& QMMMData, int bead)
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
    if (QMMMData[i].MMregion)
    {
      //Find current orientation
      double x,y,z;
      if (QMMMData[i].MP[bead].type == "Bisector")
      {
        //Find z vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
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
      if (QMMMData[i].MP[bead].type == "Z-then-X")
      {
        //Find z vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find x vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (QMMMData[i].MP[bead].type == "Z-Bisect")
      {
        //Find first vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find second vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Find third vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].z;
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
      if (QMMMData[i].MP[bead].type == "3-Fold")
      {
        //First vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Second vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Third vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].z;
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
      if (QMMMData[i].MP[bead].type == "Z-Only")
      {
        //Primary vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
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
      if (QMMMData[i].MP[bead].type == "None")
      {
        QMMMData[i].PC[bead].q1 = 0;
        QMMMData[i].PC[bead].q2 = 0;
        QMMMData[i].PC[bead].q3 = 0;
        QMMMData[i].PC[bead].q4 = 0;
        QMMMData[i].PC[bead].q5 = 0;
        QMMMData[i].PC[bead].q6 = 0;
      }
      //Fill in y vector
      vecY = vecX.cross(vecZ);
      vecY.normalize();
      //Rotate to the global frame
      MPole newPoles;
      //Add monopoles
      newPoles.q = QMMMData[i].MP[bead].q;
      //Rotate dipoles
      newPoles.Dx = 0; //X component
      newPoles.Dx += QMMMData[i].MP[bead].Dx*vecX(0);
      newPoles.Dx += QMMMData[i].MP[bead].Dy*vecY(0);
      newPoles.Dx += QMMMData[i].MP[bead].Dz*vecZ(0);
      newPoles.Dy = 0; //Y component
      newPoles.Dy += QMMMData[i].MP[bead].Dx*vecX(1);
      newPoles.Dy += QMMMData[i].MP[bead].Dy*vecY(1);
      newPoles.Dy += QMMMData[i].MP[bead].Dz*vecZ(1);
      newPoles.Dz = 0; //Z component
      newPoles.Dz += QMMMData[i].MP[bead].Dx*vecX(2);
      newPoles.Dz += QMMMData[i].MP[bead].Dy*vecY(2);
      newPoles.Dz += QMMMData[i].MP[bead].Dz*vecZ(2);
      //Add induced dipoles (Already in global frame)
      newPoles.Dx += QMMMData[i].MP[bead].IDx;
      newPoles.Dy += QMMMData[i].MP[bead].IDy;
      newPoles.Dz += QMMMData[i].MP[bead].IDz;
      newPoles.IDx = 0;
      newPoles.IDy = 0;
      newPoles.IDz = 0;
      //Rotate quadrupoles (This looks awful, but it works)
      //NB: This is a hard coded matrix rotation
      newPoles.Qxx = 0; //XX component
      newPoles.Qxx += vecX(0)*vecX(0)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qxx += vecX(0)*vecY(0)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxx += vecX(0)*vecZ(0)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxx += vecY(0)*vecX(0)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxx += vecY(0)*vecY(0)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qxx += vecY(0)*vecZ(0)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecX(0)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxx += vecZ(0)*vecY(0)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecZ(0)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qxy = 0; //XY component
      newPoles.Qxy += vecX(0)*vecX(1)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qxy += vecX(0)*vecY(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxy += vecX(0)*vecZ(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxy += vecY(0)*vecX(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxy += vecY(0)*vecY(1)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qxy += vecY(0)*vecZ(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecX(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxy += vecZ(0)*vecY(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecZ(1)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qxz = 0; //XZ component
      newPoles.Qxz += vecX(0)*vecX(2)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qxz += vecX(0)*vecY(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxz += vecX(0)*vecZ(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxz += vecY(0)*vecX(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxz += vecY(0)*vecY(2)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qxz += vecY(0)*vecZ(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecX(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxz += vecZ(0)*vecY(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecZ(2)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qyy = 0; //YY component
      newPoles.Qyy += vecX(1)*vecX(1)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qyy += vecX(1)*vecY(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyy += vecX(1)*vecZ(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyy += vecY(1)*vecX(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyy += vecY(1)*vecY(1)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qyy += vecY(1)*vecZ(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecX(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyy += vecZ(1)*vecY(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecZ(1)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qyz = 0; //YZ component
      newPoles.Qyz += vecX(1)*vecX(2)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qyz += vecX(1)*vecY(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyz += vecX(1)*vecZ(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyz += vecY(1)*vecX(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyz += vecY(1)*vecY(2)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qyz += vecY(1)*vecZ(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecX(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyz += vecZ(1)*vecY(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecZ(2)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qzz = 0; //ZZ component
      newPoles.Qzz += vecX(2)*vecX(2)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qzz += vecX(2)*vecY(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qzz += vecX(2)*vecZ(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qzz += vecY(2)*vecX(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qzz += vecY(2)*vecY(2)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qzz += vecY(2)*vecZ(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecX(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qzz += vecZ(2)*vecY(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecZ(2)*QMMMData[i].MP[bead].Qzz;
      //switch to point-charges
      QMMMData[i].PC[bead] = SphHarm2Charges(Cart2SphHarm(newPoles));
      //Translate charges to the atom's location in the global frame
      QMMMData[i].PC[bead].x1 += QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y1 += QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z1 += QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x2 += QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y2 += QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z2 += QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x3 += QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y3 += QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z3 += QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x4 += QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y4 += QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z4 += QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x5 += QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y5 += QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z5 += QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x6 += QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y6 += QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z6 += QMMMData[i].P[bead].z;
    }
    else
    {
      //Set other charges to zero
      //NB: This makes sure that variables are not undefined
      QMMMData[i].PC[bead].q1 = 0;
      QMMMData[i].PC[bead].q2 = 0;
      QMMMData[i].PC[bead].q3 = 0;
      QMMMData[i].PC[bead].q4 = 0;
      QMMMData[i].PC[bead].q5 = 0;
      QMMMData[i].PC[bead].q6 = 0;
      QMMMData[i].PC[bead].x1 = QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y1 = QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z1 = QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x2 = QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y2 = QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z2 = QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x3 = QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y3 = QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z3 = QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x4 = QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y4 = QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z4 = QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x5 = QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y5 = QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z5 = QMMMData[i].P[bead].z;
      QMMMData[i].PC[bead].x6 = QMMMData[i].P[bead].x;
      QMMMData[i].PC[bead].y6 = QMMMData[i].P[bead].y;
      QMMMData[i].PC[bead].z6 = QMMMData[i].P[bead].z;
    }
  }
  return;
};

void WriteTINKMPole(vector<QMMMAtom>& QMMMData, fstream& outFile, int i,
                    int bead)
{
  //Write a new multipole definition for pseudo-bonds and QM atoms
  outFile << "multipole -"; //Negative sign defines the frame with atom IDs
  outFile << (QMMMData[i].id+1) << " ";
  //Print frame
  if (QMMMData[i].MP[bead].type == "Z-then-X")
  {
    outFile << (QMMMData[i].MP[bead].atom1+1) << " ";
    outFile << (QMMMData[i].MP[bead].atom2+1) << " ";
    outFile << (QMMMData[i].MP[bead].atom3+1) << " ";
  }
  if (QMMMData[i].MP[bead].type == "Bisector")
  {
    outFile << "-"; //Defines the bisectors
    outFile << (QMMMData[i].MP[bead].atom1+1) << " ";
    outFile << (QMMMData[i].MP[bead].atom2+1) << " ";
    outFile << (QMMMData[i].MP[bead].atom3+1) << " ";
  }
  if (QMMMData[i].MP[bead].type == "Z-Bisector")
  {
    outFile << (QMMMData[i].MP[bead].atom1+1) << " ";
    outFile << "-"; //Defines the bisectors
    outFile << (QMMMData[i].MP[bead].atom2+1) << " ";
    outFile << (QMMMData[i].MP[bead].atom3+1) << " ";
  }
  if (QMMMData[i].MP[bead].type == "Z-Only")
  {
    outFile << (QMMMData[i].MP[bead].atom1+1) << " ";
    outFile << "0" << " ";
    outFile << "0" << " ";
  }
  if (QMMMData[i].MP[bead].type == "None")
  {
    outFile << "0" << " ";
    outFile << "0" << " ";
    outFile << "0" << " ";
  }
  //Print only the point-charge
  outFile << QMMMData[i].MP[bead].q << '\n';
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
void WriteChargeFile(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts,
                     int bead)
{
  //Function to write a file for the MM charges
  stringstream call; //Generic stream
  fstream outFile; //Stream for the charge file
  bool firstCharge = 1; //Always write the first charge
  //Find the center of mass
  Coord QMCOM; //QM region center of mass
  if (PBCon or QMMMOpts.useLREC)
  {
    QMCOM = FindQMCOM(QMMMData,QMMMOpts,bead);
  }
  //Initialize charges
  if (AMOEBA)
  {
    if (TINKER)
    {
      //Set up current multipoles
      RotateTINKCharges(QMMMData,bead);
    }
  }
  //Write charge file
  if (Gaussian or NWChem)
  {
    //Save file
    call.str("");
    call << "MMCharges_" << bead << ".txt";
    outFile.open(call.str().c_str(),ios_base::out);
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].MMregion)
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
          dx = QMMMData[i].P[bead].x-QMCOM.x;
          dy = QMMMData[i].P[bead].y-QMCOM.y;
          dz = QMMMData[i].P[bead].z-QMCOM.z;
          distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
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
            outFile << LICHEMFormFloat(QMMMData[i].P[bead].x+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].P[bead].y+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].P[bead].z+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].MP[bead].q*scrq,16);
            outFile << '\n';
          }
          if (AMOEBA)
          {
            //Add multipoles
            firstCharge = 0; //Skips writing the remaining zeros
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x1+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y1+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z1+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q1*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x2+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y2+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z2+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q2*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x3+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y3+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z3+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q3*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x4+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y4+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z4+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q4*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x5+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y5+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z5+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q5*scrq,16);
            outFile << '\n';
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x6+xShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y6+yShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z6+zShft,16);
            outFile << " ";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q6*scrq,16);
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
    call << "MMCharges_" << bead << ".txt";
    outFile.open(call.str().c_str(),ios_base::out);
    for (int i=0;i<Natoms;i++)
    {
      if (QMMMData[i].MMregion)
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
          dx = QMMMData[i].P[bead].x-QMCOM.x;
          dy = QMMMData[i].P[bead].y-QMCOM.y;
          dz = QMMMData[i].P[bead].z-QMCOM.z;
          distCent = CoordDist2(QMMMData[i].P[bead],QMCOM);
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
            outFile << LICHEMFormFloat(QMMMData[i].MP[bead].q*scrq,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].P[bead].x+xShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].P[bead].y+yShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].P[bead].z+zShft,16);
            outFile << ")" << '\n';
          }
          if (AMOEBA)
          {
            //Add multipoles
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q1*scrq,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x1+xShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y1+yShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z1+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q2*scrq,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x2+xShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y2+yShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z2+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q3*scrq,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x3+xShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y3+yShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z3+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q4*scrq,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x4+xShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y4+yShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z4+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q5*scrq,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x5+xShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y5+yShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z5+zShft,16);
            outFile << ")" << '\n';
            outFile << "Chrgfield.extern.addCharge(";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].q6*scrq,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].x6+xShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].y6+yShft,16);
            outFile << ",";
            outFile << LICHEMFormFloat(QMMMData[i].PC[bead].z6+zShft,16);
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
  vector<QMMMAtom> QMMMData;
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
  ReadLICHEMInput(xyzFile,connectFile,regionFile,QMMMData,QMMMOpts);
  //Assume there is a single bead
  int bead = 0;
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
      if (QMMMData[i].MP[bead].type == "Bisector")
      {
        //Find z vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
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
      if (QMMMData[i].MP[bead].type == "Z-then-X")
      {
        //Find z vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find x vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Subtract overlap and normalize
        vecX -= vecZ*(vecX.dot(vecZ));
        vecX.normalize();
      }
      if (QMMMData[i].MP[bead].type == "Z-Bisect")
      {
        //Find first vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Find second vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Find third vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].z;
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
      if (QMMMData[i].MP[bead].type == "3-Fold")
      {
        //First vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
        vecZ(0) = -1*x; //Correct the direction
        vecZ(1) = -1*y; //Correct the direction
        vecZ(2) = -1*z; //Correct the direction
        vecZ.normalize();
        //Second vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom2].P[bead].z;
        vecX(0) = -1*x; //Correct the direction
        vecX(1) = -1*y; //Correct the direction
        vecX(2) = -1*z; //Correct the direction
        vecX.normalize();
        //Third vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom3].P[bead].z;
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
      if (QMMMData[i].MP[bead].type == "Z-Only")
      {
        //Primary vector
        x = QMMMData[i].P[bead].x;
        x -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].x;
        y = QMMMData[i].P[bead].y;
        y -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].y;
        z = QMMMData[i].P[bead].z;
        z -= QMMMData[QMMMData[i].MP[bead].atom1].P[bead].z;
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
      MPole newPoles;
      //Add monopoles
      newPoles.q = QMMMData[i].MP[bead].q;
      //Rotate dipoles
      newPoles.Dx = 0; //X component
      newPoles.Dx += QMMMData[i].MP[bead].Dx*vecX(0);
      newPoles.Dx += QMMMData[i].MP[bead].Dy*vecY(0);
      newPoles.Dx += QMMMData[i].MP[bead].Dz*vecZ(0);
      newPoles.Dy = 0; //Y component
      newPoles.Dy += QMMMData[i].MP[bead].Dx*vecX(1);
      newPoles.Dy += QMMMData[i].MP[bead].Dy*vecY(1);
      newPoles.Dy += QMMMData[i].MP[bead].Dz*vecZ(1);
      newPoles.Dz = 0; //Z component
      newPoles.Dz += QMMMData[i].MP[bead].Dx*vecX(2);
      newPoles.Dz += QMMMData[i].MP[bead].Dy*vecY(2);
      newPoles.Dz += QMMMData[i].MP[bead].Dz*vecZ(2);
      //Add induced dipoles (Already in global frame)
      newPoles.Dx += QMMMData[i].MP[bead].IDx;
      newPoles.Dy += QMMMData[i].MP[bead].IDy;
      newPoles.Dz += QMMMData[i].MP[bead].IDz;
      newPoles.IDx = 0;
      newPoles.IDy = 0;
      newPoles.IDz = 0;
      //Rotate quadrupoles (This looks awful, but it works)
      newPoles.Qxx = 0; //XX component
      newPoles.Qxx += vecX(0)*vecX(0)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qxx += vecX(0)*vecY(0)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxx += vecX(0)*vecZ(0)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxx += vecY(0)*vecX(0)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxx += vecY(0)*vecY(0)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qxx += vecY(0)*vecZ(0)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecX(0)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxx += vecZ(0)*vecY(0)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxx += vecZ(0)*vecZ(0)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qxy = 0; //XY component
      newPoles.Qxy += vecX(0)*vecX(1)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qxy += vecX(0)*vecY(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxy += vecX(0)*vecZ(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxy += vecY(0)*vecX(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxy += vecY(0)*vecY(1)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qxy += vecY(0)*vecZ(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecX(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxy += vecZ(0)*vecY(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxy += vecZ(0)*vecZ(1)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qxz = 0; //XZ component
      newPoles.Qxz += vecX(0)*vecX(2)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qxz += vecX(0)*vecY(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxz += vecX(0)*vecZ(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxz += vecY(0)*vecX(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qxz += vecY(0)*vecY(2)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qxz += vecY(0)*vecZ(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecX(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qxz += vecZ(0)*vecY(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qxz += vecZ(0)*vecZ(2)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qyy = 0; //YY component
      newPoles.Qyy += vecX(1)*vecX(1)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qyy += vecX(1)*vecY(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyy += vecX(1)*vecZ(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyy += vecY(1)*vecX(1)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyy += vecY(1)*vecY(1)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qyy += vecY(1)*vecZ(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecX(1)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyy += vecZ(1)*vecY(1)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyy += vecZ(1)*vecZ(1)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qyz = 0; //YZ component
      newPoles.Qyz += vecX(1)*vecX(2)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qyz += vecX(1)*vecY(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyz += vecX(1)*vecZ(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyz += vecY(1)*vecX(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qyz += vecY(1)*vecY(2)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qyz += vecY(1)*vecZ(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecX(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qyz += vecZ(1)*vecY(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qyz += vecZ(1)*vecZ(2)*QMMMData[i].MP[bead].Qzz;
      newPoles.Qzz = 0; //ZZ component
      newPoles.Qzz += vecX(2)*vecX(2)*QMMMData[i].MP[bead].Qxx;
      newPoles.Qzz += vecX(2)*vecY(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qzz += vecX(2)*vecZ(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qzz += vecY(2)*vecX(2)*QMMMData[i].MP[bead].Qxy;
      newPoles.Qzz += vecY(2)*vecY(2)*QMMMData[i].MP[bead].Qyy;
      newPoles.Qzz += vecY(2)*vecZ(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecX(2)*QMMMData[i].MP[bead].Qxz;
      newPoles.Qzz += vecZ(2)*vecY(2)*QMMMData[i].MP[bead].Qyz;
      newPoles.Qzz += vecZ(2)*vecZ(2)*QMMMData[i].MP[bead].Qzz;
      //Save the global frame multipoles
      QMMMData[i].MP[bead] = newPoles;
    }
  }
  //Calculate the total charge
  double totChg = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus) \
          reduction(+:totChg)
  for (int i=0;i<Natoms;i++)
  {
    totChg += QMMMData[i].MP[bead].q;
  }
  //Write multipoles to the screen
  cout << '\n';
  cout << "Global frame multipole moments (a.u.):";
  cout << '\n' << '\n';
  for (int i=0;i<Natoms;i++)
  {
    cout << "Atom " << i << ": ";
    cout << QMMMData[i].MMTyp << '\n';
    cout << "Position (x,y,z):" << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << QMMMData[i].P[bead].x << " ";
    cout << QMMMData[i].P[bead].y << " ";
    cout << QMMMData[i].P[bead].z << '\n';
    cout << "Charge and dipole (x,y,z):" << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << QMMMData[i].MP[bead].q << " ";
    cout << QMMMData[i].MP[bead].Dx << " ";
    cout << QMMMData[i].MP[bead].Dy << " ";
    cout << QMMMData[i].MP[bead].Dz << '\n';
    cout << "Quadrupole tensor (x,y,z):" << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << QMMMData[i].MP[bead].Qxx << " ";
    cout << QMMMData[i].MP[bead].Qxy << " ";
    cout << QMMMData[i].MP[bead].Qxz << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << QMMMData[i].MP[bead].Qxy << " ";
    cout << QMMMData[i].MP[bead].Qyy << " ";
    cout << QMMMData[i].MP[bead].Qyz << '\n';
    cout << " "; //Sometimes you just need a little space
    cout << QMMMData[i].MP[bead].Qxz << " ";
    cout << QMMMData[i].MP[bead].Qyz << " ";
    cout << QMMMData[i].MP[bead].Qzz << '\n';
    cout << '\n';
  }
  cout << "Total charge: " << totChg;
  cout << '\n' << '\n';
  //Quit
  exit(0);
  return;
};

RedMPole Cart2SphHarm(MPole& pole)
{
  //Converts Cartesian multipoles to spherical harmonic multipoles
  RedMPole SHPole; //Spherical harmonic multipoles
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

OctCharges SphHarm2Charges(RedMPole pole)
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

