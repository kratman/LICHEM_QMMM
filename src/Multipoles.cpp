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
  fstream ifile,ofile; //Generic file streams
  stringstream call; //Stream for system calls and reading/writing files
  //Create TINKER xyz file from the structure
  call.str("");
  call << "LICHM_" << Bead << ".xyz";
  ofile.open(call.str().c_str(),ios_base::out);
  //Write atoms to the xyz file
  ofile << Natoms << '\n';
  for (int i=0;i<Natoms;i++)
  {
    //Write XYZ data
    ofile << setw(6) << (Struct[i].id+1);
    ofile << " ";
    ofile << setw(3) << Struct[i].MMTyp;
    ofile << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].x,16);
    ofile << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].y,16);
    ofile << " ";
    ofile << LICHEMFormFloat(Struct[i].P[Bead].z,16);
    ofile << " ";
    ofile << setw(4) << Struct[i].NumTyp;
    for (unsigned int j=0;j<Struct[i].Bonds.size();j++)
    {
      ofile << " "; //Avoids trailing spaces
      ofile << setw(6) << (Struct[i].Bonds[j]+1);
    }
    ofile << '\n';
  }
  ofile.flush();
  ofile.close();
  //Write poledit input
  call.str("");
  call << "LICHM_" << Bead << ".txt";
  ofile.open(call.str().c_str(),ios_base::out);
  ofile << "2" << '\n';
  ofile << "LICHM_" << Bead << ".xyz" << '\n';
  ofile << '\n';
  ofile.flush();
  ofile.close();
  //Run poledit
  call.str("");
  call << "poledit < LICHM_" << Bead << ".txt > LICHM_" << Bead << ".out";
  GlobalSys = system(call.str().c_str());
  //Extract multipole frames
  call.str("");
  call << "LICHM_" << Bead << ".out";
  ifile.open(call.str().c_str(),ios_base::in);
  while (!ifile.eof())
  {
    //Parse file line by line
    getline(ifile,dummy);
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
            getline(ifile,dummy);
            getline(ifile,dummy);
            getline(ifile,dummy);
            //Check if poles exist
            getline(ifile,dummy);
            stringstream line(dummy);
            line >> dummy;
            if (dummy == "Local")
            {
              //Collect definition
              line >> dummy >> Struct[i].MP[Bead].Type;
              line >> Struct[i].MP[Bead].Atom1;
              line >> Struct[i].MP[Bead].Atom2;
              line >> Struct[i].MP[Bead].Atom3;
              //Correct numbering
              Struct[i].MP[Bead].Atom1 -= 1;
              Struct[i].MP[Bead].Atom2 -= 1;
              Struct[i].MP[Bead].Atom3 -= 1;
              //Collect charge
              ifile >> dummy >> Struct[i].MP[Bead].q;
              //Collect static dipole
              ifile >> dummy;
              ifile >> Struct[i].MP[Bead].Dx;
              ifile >> Struct[i].MP[Bead].Dy;
              ifile >> Struct[i].MP[Bead].Dz;
              //Initialize induced dipole
              Struct[i].MP[Bead].IDx = 0;
              Struct[i].MP[Bead].IDy = 0;
              Struct[i].MP[Bead].IDz = 0;
              //Collect quadrupole
              ifile >> dummy;
              ifile >> Struct[i].MP[Bead].Qxx;
              ifile >> Struct[i].MP[Bead].Qxy;
              ifile >> Struct[i].MP[Bead].Qyy;
              ifile >> Struct[i].MP[Bead].Qxz;
              ifile >> Struct[i].MP[Bead].Qyz;
              ifile >> Struct[i].MP[Bead].Qzz;
            }
            else
            {
              //Initialize the "blank" multipole
              Struct[i].MP[Bead].Type = "None";
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
            getline(ifile,dummy);
          }
        }
      }
    }
  }
  ifile.close();
  //Clean up files
  call.str("");
  call << "rm -f LICHM_" << Bead << ".txt LICHM_";
  call << Bead << ".key LICHM_" << Bead << ".xyz LICHM_";
  call << Bead << ".out";
  GlobalSys = system(call.str().c_str());
  return;
};

void RotateTINKCharges(vector<QMMMAtom>& Struct, int Bead)
{
  //Switches from the local frame of reference to the global frame
  //of reference
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
  for (int i=0;i<Natoms;i++)
  {
    Vector3d Vecx,Vecy,Vecz; //Local frame vectors
    //Initialize vectors in the global frame
    Vecx(0) = 1;
    Vecx(1) = 0;
    Vecx(2) = 0;
    Vecy(0) = 0;
    Vecy(1) = 1;
    Vecy(2) = 0;
    Vecz(0) = 0;
    Vecz(1) = 0;
    Vecz(2) = 1;
    //Rotate the charges
    if (Struct[i].MMregion)
    {
      //Find current orientation
      double x,y,z;
      if (Struct[i].MP[Bead].Type == "Bisector")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Fill in z vector
        Vecz += Vecx;
        Vecz.normalize();
        //Find x vector by subtracting overlap
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "Z-then-X")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Find x vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Subtract overlap and normalize
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "Z-Bisect")
      {
        //Find first vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Find second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Find third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom3].P[Bead].z;
        Vecy(0) = -1*x; //Correct the direction
        Vecy(1) = -1*y; //Correct the direction
        Vecy(2) = -1*z; //Correct the direction
        Vecy.normalize();
        //Combine vectors
        Vecx += Vecy;
        Vecx.normalize();
        //Subtract overlap and normalize
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "3-Fold")
      {
        //First vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom3].P[Bead].z;
        Vecy(0) = -1*x; //Correct the direction
        Vecy(1) = -1*y; //Correct the direction
        Vecy(2) = -1*z; //Correct the direction
        Vecy.normalize();
        //Combine vectors and normalize
        Vecz += Vecx+Vecy;
        Vecz.normalize();
        //Find second axis by subtracting overlap and normalizing
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "Z-Only")
      {
        //Primary vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Random vector
        Vecx(0) = 2*((((double)rand())/((double)RAND_MAX))-0.5);
        Vecx(1) = 2*((((double)rand())/((double)RAND_MAX))-0.5);
        Vecx(2) = 2*((((double)rand())/((double)RAND_MAX))-0.5);
        Vecx.normalize();
        //Subtract overlap and normalize
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "None")
      {
        Struct[i].PC[Bead].q1 = 0;
        Struct[i].PC[Bead].q2 = 0;
        Struct[i].PC[Bead].q3 = 0;
        Struct[i].PC[Bead].q4 = 0;
        Struct[i].PC[Bead].q5 = 0;
        Struct[i].PC[Bead].q6 = 0;
      }
      //Fill in y vector
      Vecy = Vecx.cross(Vecz);
      Vecy.normalize();
      //Rotate to the global frame
      Mpole NewPoles;
      //Add monopoles
      NewPoles.q = Struct[i].MP[Bead].q;
      //Rotate dipoles
      NewPoles.Dx = 0; //X component
      NewPoles.Dx += Struct[i].MP[Bead].Dx*Vecx(0);
      NewPoles.Dx += Struct[i].MP[Bead].Dy*Vecy(0);
      NewPoles.Dx += Struct[i].MP[Bead].Dz*Vecz(0);
      NewPoles.Dy = 0; //Y component
      NewPoles.Dy += Struct[i].MP[Bead].Dx*Vecx(1);
      NewPoles.Dy += Struct[i].MP[Bead].Dy*Vecy(1);
      NewPoles.Dy += Struct[i].MP[Bead].Dz*Vecz(1);
      NewPoles.Dz = 0; //Z component
      NewPoles.Dz += Struct[i].MP[Bead].Dx*Vecx(2);
      NewPoles.Dz += Struct[i].MP[Bead].Dy*Vecy(2);
      NewPoles.Dz += Struct[i].MP[Bead].Dz*Vecz(2);
      //Add induced dipoles (Already in global frame)
      NewPoles.Dx += Struct[i].MP[Bead].IDx;
      NewPoles.Dy += Struct[i].MP[Bead].IDy;
      NewPoles.Dz += Struct[i].MP[Bead].IDz;
      NewPoles.IDx = 0;
      NewPoles.IDy = 0;
      NewPoles.IDz = 0;
      //Rotate quadrupoles (This looks awful, but it works)
      //NB: This is a hard coded matrix rotation
      NewPoles.Qxx = 0; //XX component
      NewPoles.Qxx += Vecx(0)*Vecx(0)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qxx += Vecx(0)*Vecy(0)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxx += Vecx(0)*Vecz(0)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxx += Vecy(0)*Vecx(0)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxx += Vecy(0)*Vecy(0)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qxx += Vecy(0)*Vecz(0)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxx += Vecz(0)*Vecx(0)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxx += Vecz(0)*Vecy(0)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxx += Vecz(0)*Vecz(0)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qxy = 0; //XY component
      NewPoles.Qxy += Vecx(0)*Vecx(1)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qxy += Vecx(0)*Vecy(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxy += Vecx(0)*Vecz(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxy += Vecy(0)*Vecx(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxy += Vecy(0)*Vecy(1)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qxy += Vecy(0)*Vecz(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxy += Vecz(0)*Vecx(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxy += Vecz(0)*Vecy(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxy += Vecz(0)*Vecz(1)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qxz = 0; //XZ component
      NewPoles.Qxz += Vecx(0)*Vecx(2)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qxz += Vecx(0)*Vecy(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxz += Vecx(0)*Vecz(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxz += Vecy(0)*Vecx(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxz += Vecy(0)*Vecy(2)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qxz += Vecy(0)*Vecz(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxz += Vecz(0)*Vecx(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxz += Vecz(0)*Vecy(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxz += Vecz(0)*Vecz(2)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qyy = 0; //YY component
      NewPoles.Qyy += Vecx(1)*Vecx(1)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qyy += Vecx(1)*Vecy(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyy += Vecx(1)*Vecz(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyy += Vecy(1)*Vecx(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyy += Vecy(1)*Vecy(1)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qyy += Vecy(1)*Vecz(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyy += Vecz(1)*Vecx(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyy += Vecz(1)*Vecy(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyy += Vecz(1)*Vecz(1)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qyz = 0; //YZ component
      NewPoles.Qyz += Vecx(1)*Vecx(2)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qyz += Vecx(1)*Vecy(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyz += Vecx(1)*Vecz(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyz += Vecy(1)*Vecx(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyz += Vecy(1)*Vecy(2)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qyz += Vecy(1)*Vecz(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyz += Vecz(1)*Vecx(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyz += Vecz(1)*Vecy(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyz += Vecz(1)*Vecz(2)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qzz = 0; //ZZ component
      NewPoles.Qzz += Vecx(2)*Vecx(2)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qzz += Vecx(2)*Vecy(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qzz += Vecx(2)*Vecz(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qzz += Vecy(2)*Vecx(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qzz += Vecy(2)*Vecy(2)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qzz += Vecy(2)*Vecz(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qzz += Vecz(2)*Vecx(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qzz += Vecz(2)*Vecy(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qzz += Vecz(2)*Vecz(2)*Struct[i].MP[Bead].Qzz;
      //switch to point-charges
      Struct[i].PC[Bead] = SphHarm2Charges(Cart2SphHarm(NewPoles));
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

void WriteTINKMpole(vector<QMMMAtom>& Struct, fstream& ofile, int i, int Bead)
{
  //Write a new multipole definition for pseudo-bonds and QM atoms
  ofile << "multipole -"; //Negative sign defines the frame with atom IDs
  ofile << (Struct[i].id+1) << " ";
  //Print frame
  if (Struct[i].MP[Bead].Type == "Z-then-X")
  {
    ofile << (Struct[i].MP[Bead].Atom1+1) << " ";
    ofile << (Struct[i].MP[Bead].Atom2+1) << " ";
    ofile << (Struct[i].MP[Bead].Atom3+1) << " ";
  }
  if (Struct[i].MP[Bead].Type == "Bisector")
  {
    ofile << "-"; //Defines the bisectors
    ofile << (Struct[i].MP[Bead].Atom1+1) << " ";
    ofile << (Struct[i].MP[Bead].Atom2+1) << " ";
    ofile << (Struct[i].MP[Bead].Atom3+1) << " ";
  }
  if (Struct[i].MP[Bead].Type == "Z-Bisector")
  {
    ofile << (Struct[i].MP[Bead].Atom1+1) << " ";
    ofile << "-"; //Defines the bisectors
    ofile << (Struct[i].MP[Bead].Atom2+1) << " ";
    ofile << (Struct[i].MP[Bead].Atom3+1) << " ";
  }
  if (Struct[i].MP[Bead].Type == "Z-Only")
  {
    ofile << (Struct[i].MP[Bead].Atom1+1) << " ";
    ofile << "0" << " ";
    ofile << "0" << " ";
  }
  if (Struct[i].MP[Bead].Type == "None")
  {
    ofile << "0" << " ";
    ofile << "0" << " ";
    ofile << "0" << " ";
  }
  //Print only the point-charge
  ofile << Struct[i].MP[Bead].q << '\n';
  //Print dummy dipoles
  ofile << "        0.0 0.0 0.0" << '\n';
  //Print dummy quadrupoles
  ofile << "        0.0" << '\n';
  ofile << "        0.0 0.0" << '\n';
  ofile << "        0.0 0.0 0.0" << '\n';
  ofile << '\n';
  return;
};

//General routines
void WriteChargeFile(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts,
                     int Bead)
{
  //Function to write a file for the MM charges
  stringstream call; //Generic stream
  fstream ofile; //Stream for the charge file
  bool FirstCharge = 1; //Always write the first charge
  //Find the center of mass
  Coord QMCOM; //QM region center of mass
  if (PBCon or QMMMOpts.UseLREC)
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
    ofile.open(call.str().c_str(),ios_base::out);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        //Check PBC (minimum image convention)
        Coord DistCent; //Distance from QM COM
        double xshft = 0;
        double yshft = 0;
        double zshft = 0;
        if (PBCon or QMMMOpts.UseLREC)
        {
          //Initialize displacements
          double dx,dy,dz; //Starting displacements
          dx = Struct[i].P[Bead].x-QMCOM.x;
          dy = Struct[i].P[Bead].y-QMCOM.y;
          dz = Struct[i].P[Bead].z-QMCOM.z;
          DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
          //Calculate the shift in positions
          //NB: Generally this work out to be +/- {Lx,Ly,Lz}
          if (PBCon)
          {
            xshft = DistCent.x-dx;
            yshft = DistCent.y-dy;
            zshft = DistCent.z-dz;
          }
        }
        //Check for long-range corrections
        double scrq = 1;
        if (QMMMOpts.UseLREC)
        {
          //Use the long-range correction
          double rcom = 0; //Distance from center of mass
          //Calculate the distance from the center of mass
          rcom = DistCent.VecMag();
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
            scrq -= scrqA*scrqA*scrqA;
          }
          else
          {
            //Delete the charge
            scrq = 0;
          }
        }
        if ((scrq > 0) or FirstCharge)
        {
          if (CHRG)
          {
            //Add charges
            FirstCharge = 0; //Skips writing the remaining zeros
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].x+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].y+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].z+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16);
            ofile << '\n';
          }
          if (AMOEBA)
          {
            //Add multipoles
            FirstCharge = 0; //Skips writing the remaining zeros
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x1+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y1+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z1+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x2+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y2+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z2+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x3+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y3+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z3+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x4+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y4+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z4+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x5+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y5+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z5+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16);
            ofile << '\n';
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x6+xshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y6+yshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z6+zshft,16);
            ofile << " ";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16);
            ofile << '\n';
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
    ofile.open(call.str().c_str(),ios_base::out);
    for (int i=0;i<Natoms;i++)
    {
      if (Struct[i].MMregion)
      {
        //Check PBC (minimum image convention)
        Coord DistCent; //Distance from QM COM
        double xshft = 0;
        double yshft = 0;
        double zshft = 0;
        if (PBCon or QMMMOpts.UseLREC)
        {
          //Initialize displacements
          double dx,dy,dz; //Starting displacements
          dx = Struct[i].P[Bead].x-QMCOM.x;
          dy = Struct[i].P[Bead].y-QMCOM.y;
          dz = Struct[i].P[Bead].z-QMCOM.z;
          DistCent = CoordDist2(Struct[i].P[Bead],QMCOM);
          //Calculate the shift in positions
          //NB: Generally this work out to be +/- {Lx,Ly,Lz}
          if (PBCon)
          {
            xshft = DistCent.x-dx;
            yshft = DistCent.y-dy;
            zshft = DistCent.z-dz;
          }
        }
        //Check for long-range corrections
        double scrq = 1;
        if (QMMMOpts.UseLREC)
        {
          //Use the long-range correction
          double rcom = 0; //Distance from center of mass
          //Calculate the distance from the center of mass
          rcom = DistCent.VecMag();
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
            scrq -= scrqA*scrqA*scrqA;
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
            ofile << "Chrgfield.extern.addCharge(";
            ofile << LICHEMFormFloat(Struct[i].MP[Bead].q*scrq,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].x+xshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].y+yshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].P[Bead].z+zshft,16);
            ofile << ")" << '\n';
          }
          if (AMOEBA)
          {
            //Add multipoles
            ofile << "Chrgfield.extern.addCharge(";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q1*scrq,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x1+xshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y1+yshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z1+zshft,16);
            ofile << ")" << '\n';
            ofile << "Chrgfield.extern.addCharge(";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q2*scrq,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x2+xshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y2+yshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z2+zshft,16);
            ofile << ")" << '\n';
            ofile << "Chrgfield.extern.addCharge(";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q3*scrq,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x3+xshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y3+yshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z3+zshft,16);
            ofile << ")" << '\n';
            ofile << "Chrgfield.extern.addCharge(";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q4*scrq,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x4+xshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y4+yshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z4+zshft,16);
            ofile << ")" << '\n';
            ofile << "Chrgfield.extern.addCharge(";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q5*scrq,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x5+xshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y5+yshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z5+zshft,16);
            ofile << ")" << '\n';
            ofile << "Chrgfield.extern.addCharge(";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].q6*scrq,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].x6+xshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].y6+yshft,16) << ",";
            ofile << LICHEMFormFloat(Struct[i].PC[Bead].z6+zshft,16);
            ofile << ")" << '\n';
          }
        }
      }
    }
  }
  //Write to files
  ofile.flush();
  ofile.close();
  //Return to the QM calculations
  return;
};

void ExtractGlobalPoles(int& argc, char**& argv)
{
  //Function to print the multipoles in the global frame
  cout << fixed;
  cout.precision(12);
  fstream xyzfile,connectfile,regionfile;
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
      xyzfilename = string(argv[i+1]);
      xyzfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-c")
    {
      confilename = string(argv[i+1]);
      connectfile.open(argv[i+1],ios_base::in);
    }
    if (dummy == "-r")
    {
      regfilename = string(argv[i+1]);
      regionfile.open(argv[i+1],ios_base::in);
    }
  }
  //Make sure input files can be read
  bool DoQuit = 0;
  if (!xyzfile.good())
  {
    cout << "Error: Could not open xyz file.";
    cout << '\n';
    cout.flush();
    DoQuit = 1;
  }
  if (!connectfile.good())
  {
    cout << "Error: Could not open connectivity file.";
    cout << '\n';
    cout.flush();
    DoQuit = 1;
  }
  if (!regionfile.good())
  {
    cout << "Error: Could not open region file.";
    cout << '\n';
    cout.flush();
    DoQuit = 1;
  }
  if (DoQuit)
  {
    //Quit with an error
    exit(0);
  }
  ReadLICHEMInput(xyzfile,connectfile,regionfile,Struct,QMMMOpts);
  //Assume there is a single bead
  int Bead = 0;
  //Rotate multipoles
  if (TINKER)
  {
    #pragma omp parallel for schedule(dynamic) num_threads(Ncpus)
    for (int i=0;i<Natoms;i++)
    {
      Vector3d Vecx,Vecy,Vecz; //Local frame vectors
      //Initialize vectors in the global frame
      Vecx(0) = 1;
      Vecx(1) = 0;
      Vecx(2) = 0;
      Vecy(0) = 0;
      Vecy(1) = 1;
      Vecy(2) = 0;
      Vecz(0) = 0;
      Vecz(1) = 0;
      Vecz(2) = 1;
      //Find current orientation and rotate multipoles
      double x,y,z;
      if (Struct[i].MP[Bead].Type == "Bisector")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Fill in z vector
        Vecz += Vecx;
        Vecz.normalize();
        //Find x vector by subtracting overlap
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "Z-then-X")
      {
        //Find z vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Find x vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Subtract overlap and normalize
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "Z-Bisect")
      {
        //Find first vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Find second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Find third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom3].P[Bead].z;
        Vecy(0) = -1*x; //Correct the direction
        Vecy(1) = -1*y; //Correct the direction
        Vecy(2) = -1*z; //Correct the direction
        Vecy.normalize();
        //Combine vectors
        Vecx += Vecy;
        Vecx.normalize();
        //Subtract overlap and normalize
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "3-Fold")
      {
        //First vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Second vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom2].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom2].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom2].P[Bead].z;
        Vecx(0) = -1*x; //Correct the direction
        Vecx(1) = -1*y; //Correct the direction
        Vecx(2) = -1*z; //Correct the direction
        Vecx.normalize();
        //Third vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom3].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom3].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom3].P[Bead].z;
        Vecy(0) = -1*x; //Correct the direction
        Vecy(1) = -1*y; //Correct the direction
        Vecy(2) = -1*z; //Correct the direction
        Vecy.normalize();
        //Combine vectors and normalize
        Vecz += Vecx+Vecy;
        Vecz.normalize();
        //Find second axis by subtracting overlap and normalizing
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      if (Struct[i].MP[Bead].Type == "Z-Only")
      {
        //Primary vector
        x = Struct[i].P[Bead].x-Struct[Struct[i].MP[Bead].Atom1].P[Bead].x;
        y = Struct[i].P[Bead].y-Struct[Struct[i].MP[Bead].Atom1].P[Bead].y;
        z = Struct[i].P[Bead].z-Struct[Struct[i].MP[Bead].Atom1].P[Bead].z;
        Vecz(0) = -1*x; //Correct the direction
        Vecz(1) = -1*y; //Correct the direction
        Vecz(2) = -1*z; //Correct the direction
        Vecz.normalize();
        //Random vector
        Vecx(0) = 2*((((double)rand())/((double)RAND_MAX))-0.5);
        Vecx(1) = 2*((((double)rand())/((double)RAND_MAX))-0.5);
        Vecx(2) = 2*((((double)rand())/((double)RAND_MAX))-0.5);
        Vecx.normalize();
        //Subtract overlap and normalize
        Vecx -= Vecz*(Vecx.dot(Vecz));
        Vecx.normalize();
      }
      //Fill in y vector
      Vecy = Vecx.cross(Vecz);
      Vecy.normalize();
      //Rotate to the global frame
      Mpole NewPoles;
      //Add monopoles
      NewPoles.q = Struct[i].MP[Bead].q;
      //Rotate dipoles
      NewPoles.Dx = 0; //X component
      NewPoles.Dx += Struct[i].MP[Bead].Dx*Vecx(0);
      NewPoles.Dx += Struct[i].MP[Bead].Dy*Vecy(0);
      NewPoles.Dx += Struct[i].MP[Bead].Dz*Vecz(0);
      NewPoles.Dy = 0; //Y component
      NewPoles.Dy += Struct[i].MP[Bead].Dx*Vecx(1);
      NewPoles.Dy += Struct[i].MP[Bead].Dy*Vecy(1);
      NewPoles.Dy += Struct[i].MP[Bead].Dz*Vecz(1);
      NewPoles.Dz = 0; //Z component
      NewPoles.Dz += Struct[i].MP[Bead].Dx*Vecx(2);
      NewPoles.Dz += Struct[i].MP[Bead].Dy*Vecy(2);
      NewPoles.Dz += Struct[i].MP[Bead].Dz*Vecz(2);
      //Add induced dipoles (Already in global frame)
      NewPoles.Dx += Struct[i].MP[Bead].IDx;
      NewPoles.Dy += Struct[i].MP[Bead].IDy;
      NewPoles.Dz += Struct[i].MP[Bead].IDz;
      NewPoles.IDx = 0;
      NewPoles.IDy = 0;
      NewPoles.IDz = 0;
      //Rotate quadrupoles (This looks awful, but it works)
      NewPoles.Qxx = 0; //XX component
      NewPoles.Qxx += Vecx(0)*Vecx(0)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qxx += Vecx(0)*Vecy(0)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxx += Vecx(0)*Vecz(0)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxx += Vecy(0)*Vecx(0)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxx += Vecy(0)*Vecy(0)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qxx += Vecy(0)*Vecz(0)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxx += Vecz(0)*Vecx(0)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxx += Vecz(0)*Vecy(0)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxx += Vecz(0)*Vecz(0)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qxy = 0; //XY component
      NewPoles.Qxy += Vecx(0)*Vecx(1)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qxy += Vecx(0)*Vecy(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxy += Vecx(0)*Vecz(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxy += Vecy(0)*Vecx(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxy += Vecy(0)*Vecy(1)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qxy += Vecy(0)*Vecz(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxy += Vecz(0)*Vecx(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxy += Vecz(0)*Vecy(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxy += Vecz(0)*Vecz(1)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qxz = 0; //XZ component
      NewPoles.Qxz += Vecx(0)*Vecx(2)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qxz += Vecx(0)*Vecy(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxz += Vecx(0)*Vecz(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxz += Vecy(0)*Vecx(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qxz += Vecy(0)*Vecy(2)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qxz += Vecy(0)*Vecz(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxz += Vecz(0)*Vecx(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qxz += Vecz(0)*Vecy(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qxz += Vecz(0)*Vecz(2)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qyy = 0; //YY component
      NewPoles.Qyy += Vecx(1)*Vecx(1)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qyy += Vecx(1)*Vecy(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyy += Vecx(1)*Vecz(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyy += Vecy(1)*Vecx(1)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyy += Vecy(1)*Vecy(1)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qyy += Vecy(1)*Vecz(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyy += Vecz(1)*Vecx(1)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyy += Vecz(1)*Vecy(1)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyy += Vecz(1)*Vecz(1)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qyz = 0; //YZ component
      NewPoles.Qyz += Vecx(1)*Vecx(2)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qyz += Vecx(1)*Vecy(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyz += Vecx(1)*Vecz(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyz += Vecy(1)*Vecx(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qyz += Vecy(1)*Vecy(2)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qyz += Vecy(1)*Vecz(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyz += Vecz(1)*Vecx(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qyz += Vecz(1)*Vecy(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qyz += Vecz(1)*Vecz(2)*Struct[i].MP[Bead].Qzz;
      NewPoles.Qzz = 0; //ZZ component
      NewPoles.Qzz += Vecx(2)*Vecx(2)*Struct[i].MP[Bead].Qxx;
      NewPoles.Qzz += Vecx(2)*Vecy(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qzz += Vecx(2)*Vecz(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qzz += Vecy(2)*Vecx(2)*Struct[i].MP[Bead].Qxy;
      NewPoles.Qzz += Vecy(2)*Vecy(2)*Struct[i].MP[Bead].Qyy;
      NewPoles.Qzz += Vecy(2)*Vecz(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qzz += Vecz(2)*Vecx(2)*Struct[i].MP[Bead].Qxz;
      NewPoles.Qzz += Vecz(2)*Vecy(2)*Struct[i].MP[Bead].Qyz;
      NewPoles.Qzz += Vecz(2)*Vecz(2)*Struct[i].MP[Bead].Qzz;
      //Save the global frame multipoles
      Struct[i].MP[Bead] = NewPoles;
    }
  }
  //Calculate the total charge
  double TotChg = 0;
  #pragma omp parallel for schedule(dynamic) num_threads(Ncpus) \
          reduction(+:TotChg)
  for (int i=0;i<Natoms;i++)
  {
    TotChg += Struct[i].MP[Bead].q;
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
  cout << "Total charge: " << TotChg;
  cout << '\n' << '\n';
  //Quit
  exit(0);
  return;
};

RedMpole Cart2SphHarm(Mpole& pole)
{
  //Converts Cartesian multipoles to spherical harmonic multipoles
  RedMpole SHpole; //Spherical harmonic multipoles
  //Diagonalize quadrupole moment tensor
  Matrix3d Qpole;
  Qpole(0,0) = pole.Qxx;
  Qpole(0,1) = pole.Qxy;
  Qpole(0,2) = pole.Qxz;
  Qpole(1,0) = pole.Qxy;
  Qpole(1,1) = pole.Qyy;
  Qpole(1,2) = pole.Qyz;
  Qpole(2,0) = pole.Qxz;
  Qpole(2,1) = pole.Qyz;
  Qpole(2,2) = pole.Qzz;
  //Change out of a.u.
  Qpole *= BohrRad*BohrRad; //NB: TINKER also divides by 3
  EigenSolver<Matrix3d> Qtensor; //There might be a better method
  Qtensor.compute(Qpole);
  Vector3d SHtensor;
  SHtensor = Qtensor.eigenvalues().real();
  //Save vector
  Matrix3d Vec = Qtensor.eigenvectors().real();
  SHpole.Vecx(0) = Vec(0,0);
  SHpole.Vecx(1) = Vec(1,0);
  SHpole.Vecx(2) = Vec(2,0);
  SHpole.Vecy(0) = Vec(0,1);
  SHpole.Vecy(1) = Vec(1,1);
  SHpole.Vecy(2) = Vec(2,1);
  SHpole.Vecz(0) = Vec(0,2);
  SHpole.Vecz(1) = Vec(1,2);
  SHpole.Vecz(2) = Vec(2,2);
  //Normalize vectors (Probably not needed)
  SHpole.Vecx.normalize();
  SHpole.Vecy.normalize();
  SHpole.Vecz.normalize();
  //Convert to spherical harmonics and rotate dipoles
  SHpole.Q00 = pole.q;
  SHpole.Q11c = 0; //X component
  SHpole.Q11c += pole.Dx*SHpole.Vecx(0);
  SHpole.Q11c += pole.Dy*SHpole.Vecx(1);
  SHpole.Q11c += pole.Dz*SHpole.Vecx(2);
  SHpole.Q11c *= BohrRad; //Change out of a.u.
  SHpole.Q11s = 0; //Y component
  SHpole.Q11s += pole.Dx*SHpole.Vecy(0);
  SHpole.Q11s += pole.Dy*SHpole.Vecy(1);
  SHpole.Q11s += pole.Dz*SHpole.Vecy(2);
  SHpole.Q11s *= BohrRad; //Change out of a.u.
  SHpole.Q10 = 0; //Z component
  SHpole.Q10 += pole.Dx*SHpole.Vecz(0);
  SHpole.Q10 += pole.Dy*SHpole.Vecz(1);
  SHpole.Q10 += pole.Dz*SHpole.Vecz(2);
  SHpole.Q10 *= BohrRad; //Change out of a.u.
  SHpole.Q22c = (SHtensor(0)-SHtensor(1))/sqrt(3); //Diagonal Qxx-Qyy
  SHpole.Q20 = SHtensor(2); //Diagonal Qzz
  return SHpole;
};

OctCharges SphHarm2Charges(RedMpole pole)
{
  //Converts spherical harmonic multipoles to point-charges
  OctCharges PCgrid; //New point-charge multipoles
  double pd = 0.25*BohrRad; //Positive displacement of the charges
  double nd = -1*pd; //Negative of the displacement
  //Charge in the +x direction
  PCgrid.q1 = pole.Q00/6;
  PCgrid.q1 += pole.Q11c/(2*pd);
  PCgrid.q1 -= pole.Q20/(6*pd*pd);
  PCgrid.q1 += pole.Q22c/(2*sqrt(3)*pd*pd);
  PCgrid.x1 = pd*pole.Vecx(0);
  PCgrid.y1 = pd*pole.Vecx(1);
  PCgrid.z1 = pd*pole.Vecx(2);
  //Charge in the +y direction
  PCgrid.q2 = pole.Q00/6;
  PCgrid.q2 += pole.Q11s/(2*pd);
  PCgrid.q2 -= pole.Q20/(6*pd*pd);
  PCgrid.q2 -= pole.Q22c/(2*sqrt(3)*pd*pd);
  PCgrid.x2 = pd*pole.Vecy(0);
  PCgrid.y2 = pd*pole.Vecy(1);
  PCgrid.z2 = pd*pole.Vecy(2);
  //Charge in the +z direction
  PCgrid.q3 = pole.Q00/6;
  PCgrid.q3 += pole.Q10/(2*pd);
  PCgrid.q3 += pole.Q20/(3*pd*pd);
  PCgrid.x3 = pd*pole.Vecz(0);
  PCgrid.y3 = pd*pole.Vecz(1);
  PCgrid.z3 = pd*pole.Vecz(2);
  //Charge in the -x direction
  PCgrid.q4 = pole.Q00/6;
  PCgrid.q4 -= pole.Q11c/(2*pd);
  PCgrid.q4 -= pole.Q20/(6*pd*pd);
  PCgrid.q4 += pole.Q22c/(2*sqrt(3)*pd*pd);
  PCgrid.x4 = nd*pole.Vecx(0);
  PCgrid.y4 = nd*pole.Vecx(1);
  PCgrid.z4 = nd*pole.Vecx(2);
  //Charge in the -y direction
  PCgrid.q5 = pole.Q00/6;
  PCgrid.q5 -= pole.Q11s/(2*pd);
  PCgrid.q5 -= pole.Q20/(6*pd*pd);
  PCgrid.q5 -= pole.Q22c/(2*sqrt(3)*pd*pd);
  PCgrid.x5 = nd*pole.Vecy(0);
  PCgrid.y5 = nd*pole.Vecy(1);
  PCgrid.z5 = nd*pole.Vecy(2);
  //Charge in the -z direction
  PCgrid.q6 = pole.Q00/6;
  PCgrid.q6 -= pole.Q10/(2*pd);
  PCgrid.q6 += pole.Q20/(3*pd*pd);
  PCgrid.x6 = nd*pole.Vecz(0);
  PCgrid.y6 = nd*pole.Vecz(1);
  PCgrid.z6 = nd*pole.Vecz(2);
  //Return
  return PCgrid;
};

