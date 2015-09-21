/*

###############################################################################
#                                                                             #
#                 LICHEM: Layered Interacting CHEmical Models                 #
#                              By: Eric G. Kratz                              #
#                                                                             #
#                      Symbiotic Computational Chemistry                      #
#                                                                             #
###############################################################################

 Functions for data analysis and trajectory processing.

*/

//Trajectory analysis functions
void Print_traj(vector<QMMMAtom>& parts, fstream& traj, QMMMSettings& QMMMOpts)
{
  //Function to print the trajectory or restart files for all beads
  stringstream call; //Only used to save traj stream settings
  call.copyfmt(traj); //Save settings
  traj.precision(8); //Adjust printing
  int Ntot = QMMMOpts.Nbeads*Natoms; //Total number of particles
  traj << Ntot << '\n' << '\n'; //Print number of particles and a blank line
  for (int i=0;i<Natoms;i++)
  {
    //Print all replicas of atom i
    for (int j=0;j<QMMMOpts.Nbeads;j++)
    {
      traj << setw(3) << parts[i].QMTyp << " ";
      traj << setw(10) << parts[i].P[j].x << " ";
      traj << setw(10) << parts[i].P[j].y << " ";
      traj << setw(10) << parts[i].P[j].z << '\n';
    }
  }
  traj.flush(); //Force printing
  traj.copyfmt(call); //Replace settings
  return;
};

void BurstTraj(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Function to split reaction path and path-integral trajectory frames
  int ct; //Generic counter
  stringstream call; //Stream for system calls and reading/writing files
  fstream burstfile;
  string dummy; //Generic string
  //Open new split trajectory file
  call.str("");
  call << "BurstStruct.xyz";
  ct = 0; //Start counting at the second file
  while (CheckFile(call.str()))
  {
    //Avoids overwriting files
    ct += 1;
    call.str("");
    call << "BurstStruct_";
    call << ct << ".xyz";
  }
  burstfile.open(call.str().c_str(),ios_base::out);
  //Print trajectory
  burstfile.precision(8);
  for (int j=0;j<QMMMOpts.Nbeads;j++)
  {
    //Print all atoms in replica j
    burstfile << Natoms; //Number of atoms
    burstfile << '\n' << '\n'; //Print blank comment line
    for (int i=0;i<Natoms;i++)
    {
      //Print data for atom i
      burstfile << setw(3) << Struct[i].QMTyp << " ";
      burstfile << setw(10) << Struct[i].P[j].x << " ";
      burstfile << setw(10) << Struct[i].P[j].y << " ";
      burstfile << setw(10) << Struct[i].P[j].z << '\n';
    }
  }
  burstfile.flush(); //Print trajectory
  burstfile.close(); //File is nolonger needed
  return;
};

//Trajectory manipulation functions
void KabschRotation(MatrixXd& A, MatrixXd& B, int MatSize)
{
  //Function to translate/rotate two structures for maximum overlap
  MatrixXd CoVar(MatSize,MatSize); //Covariance matrix
  MatrixXd RotMat(3,3); //Rotation matrix
  //Translate to the centroid
  double Ax = 0; //Average x position of matrix A
  double Ay = 0; //Average y position of matrix A
  double Az = 0; //Average z position of matrix A
  double Bx = 0; //Average x position of matrix B
  double By = 0; //Average y position of matrix B
  double Bz = 0; //Average z position of matrix B
  #pragma omp parallel for reduction(+:Ax,Ay,Az,Bx,By,Bz)
  for (int i=0;i<MatSize;i++)
  {
    //Update sum of the atomic positions
    Ax += A(i,0);
    Ay += A(i,1);
    Az += A(i,2);
    Bx += B(i,0);
    By += B(i,1);
    Bz += B(i,2);
  }
  #pragma omp barrier
  //Take average
  Ax /= MatSize;
  Ay /= MatSize;
  Az /= MatSize;
  Bx /= MatSize;
  By /= MatSize;
  Bz /= MatSize;
  #pragma omp parallel for
  for (int i=0;i<MatSize;i++)
  {
    //Move A and B to (0,0,0)
    A(i,0) -= Ax;
    A(i,1) -= Ay;
    A(i,2) -= Az;
    B(i,0) -= Bx;
    B(i,1) -= By;
    B(i,2) -= Bz;
  }
  #pragma omp barrier
  //Calculate covariance matrix
  CoVar = A.transpose()*B;
  //Compute rotation matrix and identity matrix
  JacobiSVD<MatrixXd> SVDMat(CoVar,ComputeFullU|ComputeFullV);
  MatrixXd DetMat = SVDMat.matrixV()*SVDMat.matrixU().transpose();
  int SignVal = DetMat.determinant();
  if (SignVal < 0)
  {
    SignVal = -1;
  }
  else
  {
    SignVal = 1;
  }
  MatrixXd Ident(3,3);
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<3;j++)
    {
      Ident(i,j) = 0;
    }
    Ident(i,i) = 1;
  }
  Ident(2,2) *= SignVal;
  //Find optimal rotation matrix
  RotMat = SVDMat.matrixV()*Ident*SVDMat.matrixU().transpose();
  //Rotate matrix A
  A *= RotMat;
  //Return the modified positions
  return;
};

VectorXd KabschDisplacement(MatrixXd& A, MatrixXd& B, int MatSize)
{
  //Returns the distance between two superimposed structures
  VectorXd Dist(MatSize);
  //Rotate structures
  KabschRotation(A,B,MatSize);
  //Calculate displacement
  #pragma omp parallel for
  for (int i=0;i<(3*MatSize);i++)
  {
    //Find the correct location in the arrays
    int Direc = i%3; //Find the remainder: x=0,y=1,z=2
    int AtID = (i-Direc)/3; //Find array index
    //Displacement
    Dist(i) = A(AtID,Direc)-B(AtID,Direc);
  }
  #pragma omp barrier
  //Return array
  return Dist;
};

