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

 Reference for Kabsch algorithm:
 Kabsch, Acta Crystallogr. Sect. A, 32, 922, (1976)
 Kabsch, Acta Crystallogr. Sect. A, 34, 827, (1978)

*/

/*!
  \ingroup Analysis
*/
//! \{

//Trajectory analysis functions

//! \brief Prints the current structure to a trajectory file.
//! \param QMMMData - All QMMM atomic data for the simulation
//! \param traj - File stream for all trajectory output
//! \param QMMMOpts - Simulation settings
void Print_traj(vector<QMMMAtom>& QMMMData, fstream& traj,
                QMMMSettings& QMMMOpts)
{
  //Function to print the trajectory or restart files for all beads
  stringstream call; //Only used to save traj stream settings
  //Print XYZ file
  int Ntot = QMMMOpts.NBeads*Natoms; //Total number of particles
  traj << Ntot << '\n' << '\n'; //Print number of particles and a blank line
  //Loop over the atoms in the structure
  for (int i=0;i<Natoms;i++)
  {
    //Print all replicas of atom i
    for (int j=0;j<QMMMOpts.NBeads;j++)
    {
      traj << setw(3) << left << QMMMData[i].QMTyp << " ";
      traj << LICHEMFormFloat(QMMMData[i].P[j].x,16) << " ";
      traj << LICHEMFormFloat(QMMMData[i].P[j].y,16) << " ";
      traj << LICHEMFormFloat(QMMMData[i].P[j].z,16) << '\n';
    }
  }
  //Write data and return
  traj.flush(); //Force printing
  return;
};

//! \brief Splits a multireplica trajectory into individual frames.
//! \param QMMMData - Simulation trajectory data
//! \param QMMMOpts - Simulation settings
void BurstTraj(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts)
{
  //Function to split reaction path and path-integral trajectory frames
  int ct; //Generic counter
  stringstream call; //Stream for system calls and reading/writing files
  fstream burstFile;
  string dummy; //Generic string
  //Open new trajectory file
  call.str("");
  call << "BurstStruct.xyz";
  ct = 0; //Start counting at the second file
  while (CheckFile(call.str()))
  {
    //Avoids overwriting files
    ct += 1; //Increase file counter
    call.str(""); //Change file name
    call << "BurstStruct_";
    call << ct << ".xyz";
  }
  burstFile.open(call.str().c_str(),ios_base::out);
  //Print trajectory
  for (int j=0;j<QMMMOpts.NBeads;j++)
  {
    //Print all atoms in replica j
    burstFile << Natoms; //Number of atoms
    burstFile << '\n' << '\n'; //Print blank comment line
    for (int i=0;i<Natoms;i++)
    {
      //Print data for atom i
      burstFile << setw(3) << left << QMMMData[i].QMTyp << " ";
      burstFile << LICHEMFormFloat(QMMMData[i].P[j].x,16) << " ";
      burstFile << LICHEMFormFloat(QMMMData[i].P[j].y,16) << " ";
      burstFile << LICHEMFormFloat(QMMMData[i].P[j].z,16) << '\n';
    }
  }
  //Write data and return
  burstFile.flush(); //Print trajectory
  burstFile.close(); //File is nolonger needed
  return;
};

//Trajectory manipulation functions

//! \brief Reads reactant, TS, and product data to create a path.
//! \param argc - Initial number of arguments passed to LICHEM
//! \param argv - Initial argument values passed to LICHEM
void PathLinInterpolate(int& argc, char**& argv)
{
  //Linearly interpolate a path from reactant and product geometries
  //NB: Transition state structures are optional
  stringstream call; //Stream for system calls and reading/writing files
  string dummy; //Generic string
  string reactFilename, tsFilename, prodFilename; //File names
  fstream reactFile, tsFile, prodFile, pathFile; //File streams
  int Nbeads = 3; //Default to react, ts, and prod
  bool includeTS = false; //Flag to use the TS structure
  bool doQuit = false; //Quit with an error
  reactFilename = "N/A";
  tsFilename = "N/A";
  prodFilename = "N/A";
  //Read settings
  cout << "Reading LICHEM input:";
  for (int i=0;i<argc;i++)
  {
    dummy = string(argv[i]);
    //Check number of beads
    if (dummy == "-b")
    {
      stringstream file;
      file << argv[i+1]; //Save to the stream
      file >> Nbeads; //Change to an int
    }
    //Check reactant file
    if (dummy == "-r")
    {
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open reactant file!!!";
        cout << '\n';
        doQuit = true;
      }
      reactFilename = file.str();
      reactFile.open(argv[i+1],ios_base::in);
      cout << " " << argv[i+1];
    }
    //Check transition state file
    if (dummy == "-t")
    {
      includeTS = true; //Do a 3 point interpolation
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open transition state file!!!";
        cout << '\n';
        doQuit = true;
      }
      tsFilename = file.str();
      tsFile.open(argv[i+1],ios_base::in);
      cout << " " << argv[i+1];
    }
    //Check product file
    if (dummy == "-p")
    {
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open product file!!!";
        cout << '\n';
        doQuit = true;
      }
      prodFilename = file.str();
      prodFile.open(argv[i+1],ios_base::in);
      cout << " " << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if (!CheckFile(reactFilename))
  {
    //Missing flag
    cout << "Error: Missing reactant file!!!";
    cout << '\n' << '\n';
    doQuit = true;
  }
  if ((!CheckFile(tsFilename)) && (includeTS))
  {
    //Missing flag
    cout << "Error: Missing transition state file!!!";
    cout << '\n' << '\n';
    doQuit = true;
  }
  if (!CheckFile(prodFilename))
  {
    //Missing flag
    cout << "Error: Missing product file!!!";
    cout << '\n' << '\n';
    doQuit = true;
  }
  if (doQuit)
  {
    //Exit with errors
    exit(0);
  }
  //Read geometries
  vector<string> atTyps; //Element names
  vector<Coord> reactPOS; //Reactant coordinates
  vector<Coord> transPOS; //Transition state coordinates
  vector<Coord> prodPOS; //Product coordinates
  reactFile >> Natoms; //Read number of atoms
  for (int i=0;i<Natoms;i++)
  {
    string temptyp;
    Coord temppos;
    reactFile >> temptyp;
    reactFile >> temppos.x;
    reactFile >> temppos.y;
    reactFile >> temppos.z;
    atTyps.push_back(temptyp);
    reactPOS.push_back(temppos);
  }
  reactFile.close();
  if (includeTS)
  {
    tsFile >> Natoms; //Read number of atoms
    for (int i=0;i<Natoms;i++)
    {
      string temptyp;
      Coord temppos;
      tsFile >> temptyp;
      tsFile >> temppos.x;
      tsFile >> temppos.y;
      tsFile >> temppos.z;
      atTyps.push_back(temptyp);
      transPOS.push_back(temppos);
    }
    tsFile.close();
  }
  prodFile >> Natoms; //Read number of atoms
  for (int i=0;i<Natoms;i++)
  {
    string temptyp;
    Coord temppos;
    prodFile >> temptyp;
    prodFile >> temppos.x;
    prodFile >> temppos.y;
    prodFile >> temppos.z;
    atTyps.push_back(temptyp);
    prodPOS.push_back(temppos);
  }
  prodFile.close();
  //Check for more errors
  if (reactPOS.size() != prodPOS.size())
  {
    //Incorrect number of atoms in the reactant or product
    cout << "Error: Different number of atoms for the reactant";
    cout << " and product!!!";
    cout << '\n' << '\n';
    exit(0);
  }
  else if (includeTS)
  {
    if (reactPOS.size() != transPOS.size())
    {
      //Incorrect number of atoms in the reactant or product
      cout << "Error: Different number of atoms for the reactant";
      cout << " and transition state!!!";
      cout << '\n' << '\n';
      exit(0);
    }
    if (transPOS.size() != prodPOS.size())
    {
      //Incorrect number of atoms in the reactant or product
      cout << "Error: Different number of atoms for the transition state";
      cout << " and product!!!";
      cout << '\n' << '\n';
      exit(0);
    }
  }
  //Interpolate between points
  pathFile.open("BeadStartStruct.xyz",ios_base::out);
  pathFile << (Natoms*Nbeads) << '\n' << '\n';
  if (includeTS)
  {
    //Linear interpolation between the react, ts, and prod structures
    if ((Nbeads%2) != 1)
    {
      //Adjust number of beads
      Nbeads += 1;
      cout << "Warning: A three structure interpolation requires an odd";
      cout << " number of points." << '\n';
      cout << " Nbeads increased to " << Nbeads << '\n' << '\n';
    }
    for (int i=0;i<Natoms;i++)
    {
      //Loop over first half of the beads
      for (int j=0;j<((Nbeads-1)/2);j++)
      {
        //Print element
        pathFile << atTyps[i] << " ";
        //Print interpolated coordinates
        double x,y,z;
        x = reactPOS[i].x;
        x += (2*j*transPOS[i].x)/(Nbeads-1);
        x -= (2*j*reactPOS[i].x)/(Nbeads-1);
        pathFile << LICHEMFormFloat(x,16) << " ";
        y = reactPOS[i].y;
        y += (2*j*transPOS[i].y)/(Nbeads-1);
        y -= (2*j*reactPOS[i].y)/(Nbeads-1);
        pathFile << LICHEMFormFloat(y,16) << " ";
        z = reactPOS[i].z;
        z += (2*j*transPOS[i].z)/(Nbeads-1);
        z -= (2*j*reactPOS[i].z)/(Nbeads-1);
        pathFile << LICHEMFormFloat(z,16) << '\n';
      }
      //Loop over second half of the beads
      for (int j=0;j<(((Nbeads-1)/2)+1);j++)
      {
        //Print element
        pathFile << atTyps[i] << " ";
        //Print interpolated coordinates
        double x,y,z;
        x = transPOS[i].x;
        x += (2*j*prodPOS[i].x)/(Nbeads-1);
        x -= (2*j*transPOS[i].x)/(Nbeads-1);
        pathFile << LICHEMFormFloat(x,16) << " ";
        y = transPOS[i].y;
        y += (2*j*prodPOS[i].y)/(Nbeads-1);
        y -= (2*j*transPOS[i].y)/(Nbeads-1);
        pathFile << LICHEMFormFloat(y,16) << " ";
        z = transPOS[i].z;
        z += (2*j*prodPOS[i].z)/(Nbeads-1);
        z -= (2*j*transPOS[i].z)/(Nbeads-1);
        pathFile << LICHEMFormFloat(z,16) << '\n';
      }
    }
  }
  else
  {
    //Linear interpolation between the react and prod structures
    for (int i=0;i<Natoms;i++)
    {
      //Loop over beads
      for (int j=0;j<Nbeads;j++)
      {
        //Print element
        pathFile << atTyps[i] << " ";
        //Print interpolated coordinates
        double x,y,z;
        x = reactPOS[i].x;
        x += (j*prodPOS[i].x)/(Nbeads-1);
        x -= (j*reactPOS[i].x)/(Nbeads-1);
        pathFile << LICHEMFormFloat(x,16) << " ";
        y = reactPOS[i].y;
        y += (j*prodPOS[i].y)/(Nbeads-1);
        y -= (j*reactPOS[i].y)/(Nbeads-1);
        pathFile << LICHEMFormFloat(y,16) << " ";
        z = reactPOS[i].z;
        z += (j*prodPOS[i].z)/(Nbeads-1);
        z -= (j*reactPOS[i].z)/(Nbeads-1);
        pathFile << LICHEMFormFloat(z,16) << '\n';
      }
    }
  }
  pathFile.flush();
  pathFile.close();
  //Exit LICHEM
  exit(0);
  return;
};

//! \brief Reads a multireplica trajectory file and separates the replicas.
//! \param argc - Initial number of arguments passed to LICHEM
//! \param argv - Initial argument values passed to LICHEM
void SplitPathTraj(int& argc, char**& argv)
{
  //Function to separate a reaction path frame into a trajectory
  stringstream call; //Stream for system calls and reading/writing files
  string dummy; //Generic string
  string pathFilename; //Name of the merged trajectory file
  fstream pathFile, burstFile; //File streams
  int ct; //Generic counter
  int Nbeads = 1; //Default to a single bead
  int frameID = 0; //The trajectory frame that will be separated
  bool doQuit = false; //Quit with an error
  pathFilename = "N/A";
  //Read settings
  cout << "Reading LICHEM input:";
  for (int i=0;i<argc;i++)
  {
    dummy = string(argv[i]);
    //Check number of beads
    if (dummy == "-b")
    {
      stringstream file;
      file << argv[i+1]; //Save to the stream
      file >> Nbeads; //Change to an int
    }
    //Check frame
    if (dummy == "-f")
    {
      stringstream file;
      file << argv[i+1]; //Save to the stream
      file >> frameID; //Change to an int
    }
    //Read reaction path trajectory file name
    if (dummy == "-p")
    {
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open trajectory file!!!";
        cout << '\n';
        doQuit = true;
      }
      pathFilename = file.str();
      pathFile.open(argv[i+1],ios_base::in);
      cout << " " << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Print other settings
  cout << "Number of beads: " << Nbeads << '\n';
  cout << "Frame ID: " << frameID << '\n';
  cout << '\n';
  //Check for errors
  if (!CheckFile(pathFilename))
  {
    //Missing flag
    cout << "Error: Missing trajectory file!!!";
    cout << '\n' << '\n';
    doQuit = true;
  }
  if (doQuit)
  {
    //Exit with errors
    exit(0);
  }
  //Open new trajectory file
  call.str("");
  call << "BurstStruct.xyz";
  ct = 0; //Start counting at the second file
  while (CheckFile(call.str()))
  {
    //Avoids overwriting files
    ct += 1; //Increase file counter
    call.str(""); //Change file name
    call << "BurstStruct_";
    call << ct << ".xyz";
  }
  burstFile.open(call.str().c_str(),ios_base::out);
  cout << "Trajectory output: " << call.str();
  cout << '\n' << '\n';
  //Read the number of atoms
  getline(pathFile,dummy); //Read the first line of the file
  call.str(dummy); //Save to a stream
  call >> Natoms; //Read number of atoms
  if ((Natoms%Nbeads) != 0)
  {
    //Check the number of particles
    cout << "Error: Number of beads does not match the number of";
    cout << " particles!!!";
    cout << '\n' << '\n';
    exit(0);
  }
  else
  {
    Natoms /= Nbeads;
  }
  //Clear first comment line
  getline(pathFile,dummy); //Read and discard junk
  //Move to the correct frame
  for (int i=0;i<frameID;i++)
  {
    for (int j=0;j<(Natoms*Nbeads+2);j++)
    {
      string line;
      getline(pathFile,line); //Read and discard junk
    }
  }
  //Create string array for the frames
  vector<string> allFrames;
  for (int i=0;i<Nbeads;i++)
  {
    //Each element is a long string holding a single frame
    stringstream line;
    line.str("");
    line << Natoms << '\n' << '\n';
    allFrames.push_back(line.str());
  }
  //Separate coordinates
  for (int i=0;i<Natoms;i++)
  {
    for (int j=0;j<Nbeads;j++)
    {
      string line;
      getline(pathFile,line); //Read a line
      allFrames[j] += line; //Append to the frame
      allFrames[j] += '\n'; //Terminate the line
    }
  }
  //Write the output file
  for (int i=0;i<Nbeads;i++)
  {
    burstFile << allFrames[i];
  }
  //Close files
  pathFile.close();
  burstFile.flush();
  burstFile.close();
  //Exit
  exit(0);
  return;
};

//! \brief Translates and rotates two structures for maximum overlap.
//! \param A - Matrix of coordinates for structure A
//! \param B - Matrix of coordinates for structure B
//! \param matSize - Number of atoms in structure A/B
void KabschRotation(MatrixXd& A, MatrixXd& B, int matSize)
{
  //Function to translate/rotate two structures for maximum overlap
  MatrixXd coVar; //Covariance matrix
  MatrixXd rotMat; //Rotation matrix
  //Calculate the center of the structures
  double Ax = 0; //Average x position of matrix A
  double Ay = 0; //Average y position of matrix A
  double Az = 0; //Average z position of matrix A
  double Bx = 0; //Average x position of matrix B
  double By = 0; //Average y position of matrix B
  double Bz = 0; //Average z position of matrix B
  #pragma omp parallel
  {
    //Update sum of the atomic positions
    #pragma omp for nowait schedule(dynamic) reduction(+:Ax)
    for (int i=0;i<matSize;i++)
    {
      Ax += A(i,0);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Ay)
    for (int i=0;i<matSize;i++)
    {
      Ay += A(i,1);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Az)
    for (int i=0;i<matSize;i++)
    {
      Az += A(i,2);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Bx)
    for (int i=0;i<matSize;i++)
    {
      Bx += B(i,0);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:By)
    for (int i=0;i<matSize;i++)
    {
      By += B(i,1);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Bz)
    for (int i=0;i<matSize;i++)
    {
      Bz += B(i,2);
    }
  }
  #pragma omp barrier
  //Take average
  Ax /= matSize;
  Ay /= matSize;
  Az /= matSize;
  Bx /= matSize;
  By /= matSize;
  Bz /= matSize;
  //Translate centroids
  #pragma omp parallel
  {
    //Move A and B to (0,0,0)
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<matSize;i++)
    {
      A(i,0) -= Ax;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<matSize;i++)
    {
      A(i,1) -= Ay;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<matSize;i++)
    {
      A(i,2) -= Az;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<matSize;i++)
    {
      B(i,0) -= Bx;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<matSize;i++)
    {
      B(i,1) -= By;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<matSize;i++)
    {
      B(i,2) -= Bz;
    }
  }
  #pragma omp barrier
  //Calculate covariance matrix
  coVar = (B.transpose())*A;
  //Compute SVD and identity matrix
  JacobiSVD<MatrixXd> SVDMat(coVar,ComputeFullU|ComputeFullV);
  MatrixXd detMat = SVDMat.matrixV()*(SVDMat.matrixU().transpose());
  double signVal = detMat.determinant();
  if (signVal < 0)
  {
    signVal = -1;
  }
  else
  {
    signVal = 1;
  }
  MatrixXd Ident(3,3);
  Ident.setIdentity();
  Ident(2,2) *= signVal; //Change sign for rotation
  //Find optimal rotation matrix
  rotMat = SVDMat.matrixV()*Ident*(SVDMat.matrixU().transpose());
  //Rotate matrix A
  B *= rotMat;
  //Return the modified positions
  return;
};

//! \brief Calculates the displacement between two structures.
//! \param A - Matrix of coordinates for structure A
//! \param B - Matrix of coordinates for structure B
//! \param matSize - Number of atoms in structure A/B
//! \return dist - Total distance between the two structures
VectorXd KabschDisplacement(MatrixXd& A, MatrixXd& B, int matSize)
{
  //Returns the distance between two superimposed structures
  VectorXd dist(3*matSize);
  //Rotate structures
  KabschRotation(A,B,matSize);
  //Calculate displacement
  #pragma omp parallel for schedule(dynamic)
  for (int i=0;i<(3*matSize);i++)
  {
    //Find the correct location in the arrays
    int direc = i%3; //Find the remainder: x=0,y=1,z=2
    int atID = (i-direc)/3; //Find array index
    //Calculate displacement
    dist(i) = A(atID,direc)-B(atID,direc);
  }
  //Return array
  return dist;
};

//Physical property analysis functions

//! \brief Calculates the density of a periodic structure.
//! \param QMMMData - Simulation trajectory data
//! \param QMMMOpts - Simulation settings
//! \return rho - Density of the simulation box
double LICHEMDensity(vector<QMMMAtom>& QMMMData, QMMMSettings& QMMMOpts)
{
  //Function to calculate the density for periodic calculations
  double rho = 0;
  //Sum the masses
  #pragma omp parallel for schedule(dynamic) reduction(+:rho)
  for (int i=0;i<Natoms;i++)
  {
    rho += QMMMData[i].m;
  }
  //Divide by the volume
  rho /= (Lx*Ly*Lz);
  //Change units to SI
  rho *= (amu2kg*m2Ang*m2Ang*m2Ang);
  //Change to g/cm^3
  rho /= 1000; //m^3 to cm^3
  //Return density
  return rho;
};

//! \brief Calculates the frequencies of the QM atoms in the system.
//! \param QMMMData - Simulation trajectory data
//! \param QMMMHess - Input Hessian matrix
//! \param QMMMOpts - Simulation settings
//! \param bead - Replica used to calculate the frequencies
//! \param remCt - The number of low frequency modes removed in the analysis
//! \return QMMMFreqs - Array containing the harmonic frequencies
VectorXd LICHEMFreq(vector<QMMMAtom>& QMMMData, MatrixXd& QMMMHess,
                    QMMMSettings& QMMMOpts, int bead, int& remCt)
{
  //Function to perform a QMMM frequency analysis
  double projTol = 0.65; //Amount of overlap to remove a mode
  double zeroTol = 1.00; //Smallest possible frequency (cm^-1)
  int transRotCt = 0; //Number of deleted translation and rotational modes
  //Define variables
  int Ndof = 3*(Nqm+Npseudo); //Degrees of freedom
  //Define arrays
  EigenSolver<MatrixXd> freqAnalysis;
  VectorXd QMMMFreqs(Ndof); //Vibrational frequencies (cm^-1)
  MatrixXd QMMMNormModes(Ndof,Ndof); //Normal modes
  MatrixXd freqMatrix(Ndof,Ndof); //Diagonal frequency matrix
  VectorXd transX(Ndof); //X translation
  VectorXd transY(Ndof); //Y translation
  VectorXd transZ(Ndof); //Z translation
  VectorXd rotX(Ndof); //X rotation
  VectorXd rotY(Ndof); //Y rotation
  VectorXd rotZ(Ndof); //Z rotation
  //Initialize arrays
  QMMMFreqs.setZero();
  QMMMNormModes.setZero();
  freqMatrix.setZero();
  transX.setZero();
  transY.setZero();
  transZ.setZero();
  rotX.setZero();
  rotY.setZero();
  rotZ.setZero();
  //Collect QM and PB masses
  vector<double> masses;
  for (int i=0;i<Natoms;i++)
  {
    //Locate QM and PB atoms
    if (QMMMData[i].QMRegion || QMMMData[i].PBRegion)
    {
      //Switch to a.u. and save mass
      double massVal = QMMMData[i].m/elecMass;
      masses.push_back(massVal); //X component
      masses.push_back(massVal); //Y component
      masses.push_back(massVal); //Z component
    }
  }
  //Mass scale the Hessian matrix
  #pragma omp parallel for
  for (int i=0;i<Ndof;i++)
  {
    //Update diagonal matrix elements
    QMMMHess(i,i) /= masses[i];
    //Update off-diagonal matrix elements
    for (int j=0;j<i;j++)
    {
      //Mass scale
      QMMMHess(i,j) /= sqrt(masses[i]*masses[j]);
      //Apply symmetry
      QMMMHess(j,i) = QMMMHess(i,j);
    }
  }
  //Diagonalize Hessian matrix
  freqAnalysis.compute(QMMMHess);
  QMMMFreqs = freqAnalysis.eigenvalues().real();
  QMMMNormModes = freqAnalysis.eigenvectors().real();
  //Create model translation and rotation modes
  #pragma omp parallel for
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Translational modes
    transX(3*i) = sqrt(masses[3*i]);
    transY(3*i+1) = sqrt(masses[3*i+1]);
    transZ(3*i+2) = sqrt(masses[3*i+2]);
  }
  transX.normalize();
  transY.normalize();
  transZ.normalize();
  //Remove translation and rotation
  transRotCt = 0; //Use as a counter
  #pragma omp parallel for reduction(+:transRotCt)
  for (int i=0;i<Ndof;i++)
  {
    bool isTransRot = false;
    double dotTest; //Saves overlap
    //Locate translations
    dotTest = abs(transX.dot(QMMMNormModes.col(i)));
    if (dotTest > projTol)
    {
      isTransRot = true;
    }
    dotTest = abs(transY.dot(QMMMNormModes.col(i)));
    if (dotTest > projTol)
    {
      isTransRot = true;
    }
    dotTest = abs(transZ.dot(QMMMNormModes.col(i)));
    if (dotTest > projTol)
    {
      isTransRot = true;
    }
    //Locate rotations
    dotTest = abs(rotX.dot(QMMMNormModes.col(i)));
    if (dotTest > projTol)
    {
      isTransRot = true;
    }
    dotTest = abs(rotY.dot(QMMMNormModes.col(i)));
    if (dotTest > projTol)
    {
      isTransRot = true;
    }
    dotTest = abs(rotZ.dot(QMMMNormModes.col(i)));
    if (dotTest > projTol)
    {
      isTransRot = true;
    }
    //Locate small frequencies
    double smallFreq; //Temporary storage
    smallFreq = zeroTol/har2Wavenum; //Convert tolerance to a.u.
    smallFreq *= smallFreq; //Square tolerance
    if (abs(QMMMFreqs(i)) < smallFreq)
    {
      isTransRot = true;
    }
    //Adjust frequencies
    if (isTransRot && (!QMMM))
    {
      //Remove frequency
      QMMMFreqs(i);
      freqMatrix(i,i) = 0;
      transRotCt += 1;
    }
    else
    {
      //Save frequency
      freqMatrix(i,i) = QMMMFreqs(i);
    }
  }
  if (transRotCt > 0)
  {
    //Remove unwanted frequencies
    QMMMHess = QMMMNormModes*freqMatrix*QMMMNormModes.inverse();
    freqAnalysis.compute(QMMMHess);
    QMMMFreqs = freqAnalysis.eigenvalues().real();
    QMMMNormModes = freqAnalysis.eigenvectors().real();
  }
  //Take the square root and keep the sign
  #pragma omp parallel for
  for (int i=0;i<Ndof;i++)
  {
    //Save sign
    int freqSign = 1;
    if (QMMMFreqs(i) < 0)
    {
      freqSign = -1;
    }
    //Remove sign
    QMMMFreqs(i) = abs(QMMMFreqs(i));
    //Take square root
    QMMMFreqs(i) = sqrt(QMMMFreqs(i));
    //Replace sign
    QMMMFreqs(i) *= freqSign;
  }
  //Change units
  QMMMFreqs *= har2Wavenum;
  //Remove negligible frequencies
  transRotCt = 0; //Reset counter
  #pragma omp parallel for reduction(+:transRotCt)
  for (int i=0;i<Ndof;i++)
  {
    //Delete frequencies below the tolerance
    if (abs(QMMMFreqs(i)) < zeroTol)
    {
      transRotCt += 1;
      QMMMFreqs(i) = 0;
    }
  }
  //Write normal modes
  if (QMMMOpts.printNormModes)
  {
    if (!QMMM)
    {
      //Write all normal modes
      WriteModes(QMMMData,0,QMMMFreqs,QMMMNormModes,QMMMOpts,bead);
    }
    else
    {
      //Write modes for imaginary frequencies
      WriteModes(QMMMData,1,QMMMFreqs,QMMMNormModes,QMMMOpts,bead);
    }
  }
  //Return frequencies
  remCt = transRotCt;
  return QMMMFreqs;
};

//! \brief Prints trajectory files for animating the normal modes.
//! \param QMMMData - Simulation trajectory data
//! \param imagOnly - Flag to print only the imaginary frequencies
//! \param freqs - Input harmonic frequencies
//! \param normModes - Input normal modes
//! \param QMMMOpts - Simulation settings
//! \param bead - Replica used in the animation
void WriteModes(vector<QMMMAtom>& QMMMData, bool imagOnly, VectorXd& freqs,
                MatrixXd& normModes, QMMMSettings& QMMMOpts, int bead)
{
  //Function to write normal modes
  int Nframes = 5; //Trajectory files have (2N+1) frames
  double ampFrac = 0.50; //Fractional amplitude of the mode
  stringstream call; //Generic stream
  fstream modeFile; //Normal mode output file
  int Ndof = 3*(Nqm+Npseudo); //Number of vibrational modes
  int ct; //Generic counter
  double curAmp; //Stores the current amplitude for each frame
  for (int i=0;i<Ndof;i++)
  {
    //Check print options
    if (((!imagOnly) || (freqs(i) < 0)) && (freqs(i) != 0))
    {
      //Print normal mode
      call.str("");
      //Create file
      call << "NormModes_" << i << ".xyz";
      modeFile.open(call.str().c_str(),ios_base::out);
      //Write file
      for (int j=0;j<(2*Nframes+1);j++)
      {
        //Loop over frames
        curAmp = -1*ampFrac; //Current amplitude
        curAmp += (j*ampFrac)/Nframes; //Move along the trajectory
        ct = 0; //Keeps track of the atom IDs
        modeFile << (Nqm+Npseudo) << '\n' << '\n';
        for (int k=0;k<Natoms;k++)
        {
          if (QMMMData[k].QMRegion || QMMMData[k].PBRegion)
          {
            //Write element
            modeFile << QMMMData[k].QMTyp << " ";
            //Write X component
            modeFile << QMMMData[k].P[bead].x+(curAmp*normModes(ct,i));
            modeFile << " ";
            ct += 1;
            //Write Y component
            modeFile << QMMMData[k].P[bead].y+(curAmp*normModes(ct,i));
            modeFile << " ";
            ct += 1;
            //Write Z component
            modeFile << QMMMData[k].P[bead].z+(curAmp*normModes(ct,i));
            modeFile << '\n';
            ct += 1;
          }
        }
      }
      modeFile.flush();
      modeFile.close();
    }
  }
  return;
};

//End of file group
//! \}

