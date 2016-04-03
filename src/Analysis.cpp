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

//Trajectory analysis functions
void Print_traj(vector<QMMMAtom>& Struct, fstream& traj,
                QMMMSettings& QMMMOpts)
{
  //Function to print the trajectory or restart files for all beads
  stringstream call; //Only used to save traj stream settings
  //Print XYZ file
  int Ntot = QMMMOpts.Nbeads*Natoms; //Total number of particles
  traj << Ntot << '\n' << '\n'; //Print number of particles and a blank line
  //Loop over the atoms in the structure
  for (int i=0;i<Natoms;i++)
  {
    //Print all replicas of atom i
    for (int j=0;j<QMMMOpts.Nbeads;j++)
    {
      traj << setw(3) << left << Struct[i].QMTyp << " ";
      traj << LICHEMFormFloat(Struct[i].P[j].x,16) << " ";
      traj << LICHEMFormFloat(Struct[i].P[j].y,16) << " ";
      traj << LICHEMFormFloat(Struct[i].P[j].z,16) << '\n';
    }
  }
  //Write data and return
  traj.flush(); //Force printing
  return;
};

void BurstTraj(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Function to split reaction path and path-integral trajectory frames
  int ct; //Generic counter
  stringstream call; //Stream for system calls and reading/writing files
  fstream burstfile;
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
  burstfile.open(call.str().c_str(),ios_base::out);
  //Print trajectory
  for (int j=0;j<QMMMOpts.Nbeads;j++)
  {
    //Print all atoms in replica j
    burstfile << Natoms; //Number of atoms
    burstfile << '\n' << '\n'; //Print blank comment line
    for (int i=0;i<Natoms;i++)
    {
      //Print data for atom i
      burstfile << setw(3) << left << Struct[i].QMTyp << " ";
      burstfile << LICHEMFormFloat(Struct[i].P[j].x,16) << " ";
      burstfile << LICHEMFormFloat(Struct[i].P[j].y,16) << " ";
      burstfile << LICHEMFormFloat(Struct[i].P[j].z,16) << '\n';
    }
  }
  //Write data and return
  burstfile.flush(); //Print trajectory
  burstfile.close(); //File is nolonger needed
  return;
};

//Trajectory manipulation functions
void PathLinInterpolate(int& argc, char**& argv)
{
  //Linearly interpolate a path from reactant and product geometries
  //NB: Transition state structures are optional
  stringstream call; //Stream for system calls and reading/writing files
  string dummy; //Generic string
  string reactfilename, tsfilename, prodfilename; //File names
  fstream reactfile, tsfile, prodfile, pathfile; //File streams
  int Nbeads = 3; //Default to react, ts, and prod
  bool IncludeTS = 0; //Flag to use the TS structure
  bool DoQuit = 0; //Quit with an error
  //Read settings
  cout << "Reading LICHEM input: ";
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
        DoQuit = 1;
      }
      reactfilename = file.str();
      reactfile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
    }
    //Check transition state file
    if (dummy == "-t")
    {
      IncludeTS = 1; //Do a 3 point interpolation
      stringstream file;
      file << argv[i+1];
      if (!CheckFile(file.str()))
      {
        cout << "Error: Could not open transition state file!!!";
        cout << '\n';
        DoQuit = 1;
      }
      tsfilename = file.str();
      tsfile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
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
        DoQuit = 1;
      }
      prodfilename = file.str();
      prodfile.open(argv[i+1],ios_base::in);
      cout << argv[i+1];
    }
  }
  cout << '\n' << '\n'; //Terminate output
  //Error check
  if (!CheckFile(reactfilename))
  {
    //Missing flag
    cout << "Error: Missing reactant file!!!";
    cout << '\n' << '\n';
    DoQuit = 1;
  }
  if ((!CheckFile(tsfilename)) and (IncludeTS))
  {
    //Missing flag
    cout << "Error: Missing transition state file!!!";
    cout << '\n' << '\n';
    DoQuit = 1;
  }
  if (!CheckFile(prodfilename))
  {
    //Missing flag
    cout << "Error: Missing product file!!!";
    cout << '\n' << '\n';
    DoQuit = 1;
  }
  if (DoQuit)
  {
    //Exit with errors
    exit(0);
  }
  //Read geometries
  vector<string> AtTyps; //Element names
  vector<Coord> ReactPOS; //Reactant coordinates
  vector<Coord> TransPOS; //Transition state coordinates
  vector<Coord> ProdPOS; //Product coordinates
  reactfile >> Natoms; //Read number of atoms
  for (int i=0;i<Natoms;i++)
  {
    string temptyp;
    Coord temppos;
    reactfile >> temptyp;
    reactfile >> temppos.x;
    reactfile >> temppos.y;
    reactfile >> temppos.z;
    AtTyps.push_back(temptyp);
    ReactPOS.push_back(temppos);
  }
  reactfile.close();
  if (IncludeTS)
  {
    tsfile >> Natoms; //Read number of atoms
    for (int i=0;i<Natoms;i++)
    {
      string temptyp;
      Coord temppos;
      tsfile >> temptyp;
      tsfile >> temppos.x;
      tsfile >> temppos.y;
      tsfile >> temppos.z;
      AtTyps.push_back(temptyp);
      TransPOS.push_back(temppos);
    }
    tsfile.close();
  }
  prodfile >> Natoms; //Read number of atoms
  for (int i=0;i<Natoms;i++)
  {
    string temptyp;
    Coord temppos;
    prodfile >> temptyp;
    prodfile >> temppos.x;
    prodfile >> temppos.y;
    prodfile >> temppos.z;
    AtTyps.push_back(temptyp);
    ProdPOS.push_back(temppos);
  }
  prodfile.close();
  //Check for more errors
  if (ReactPOS.size() != ProdPOS.size())
  {
    //Incorrect number of atoms in the reactant or product
    cout << "Error: Different number of atoms for the reactant";
    cout << " and product!!!";
    cout << '\n' << '\n';
    exit(0);
  }
  else if (IncludeTS)
  {
    if (ReactPOS.size() != TransPOS.size())
    {
      //Incorrect number of atoms in the reactant or product
      cout << "Error: Different number of atoms for the reactant";
      cout << " and transition state!!!";
      cout << '\n' << '\n';
      exit(0);
    }
    if (TransPOS.size() != ProdPOS.size())
    {
      //Incorrect number of atoms in the reactant or product
      cout << "Error: Different number of atoms for the transition state";
      cout << " and product!!!";
      cout << '\n' << '\n';
      exit(0);
    }
  }
  //Interpolate between points
  pathfile.open("BeadStartStruct.xyz",ios_base::out);
  pathfile << (Natoms*Nbeads) << '\n' << '\n';
  if (IncludeTS)
  {
    //Linear interpolation between the react, ts, and prod structures
    
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
        pathfile << AtTyps[i];
        //Print interpolated coordinates
        double x,y,z;
        x = ReactPOS[i].x;
        x += (j*ProdPOS[i].x)/(Nbeads-1);
        x -= (j*ReactPOS[i].x)/(Nbeads-1);
        pathfile << LICHEMFormFloat(x,16) << " ";
        y = ReactPOS[i].y;
        y += (j*ProdPOS[i].y)/(Nbeads-1);
        y -= (j*ReactPOS[i].y)/(Nbeads-1);
        pathfile << LICHEMFormFloat(y,16) << " ";
        z = ReactPOS[i].z;
        z += (j*ProdPOS[i].z)/(Nbeads-1);
        z -= (j*ReactPOS[i].z)/(Nbeads-1);
        pathfile << LICHEMFormFloat(z,16) << '\n';

      }
    }
  }
  pathfile.flush();
  pathfile.close();
  //Exit LICHEM
  exit(0);
  return;
};

void ReorderQMPBBA(int& argc, char**& argv)
{
  //Rewrite the xyzfile to put QM, PB, and BA atoms first
  stringstream call; //Stream for system calls and reading/writing files
  string dummy; //Generic string
  //Analyze structure
  
  //Print structure
  
  //Exit LICHEM
  exit(0);
  return;
};

void KabschRotation(MatrixXd& A, MatrixXd& B, int MatSize)
{
  //Function to translate/rotate two structures for maximum overlap
  MatrixXd CoVar; //Covariance matrix
  MatrixXd RotMat; //Rotation matrix
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
    for (int i=0;i<MatSize;i++)
    {
      Ax += A(i,0);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Ay)
    for (int i=0;i<MatSize;i++)
    {
      Ay += A(i,1);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Az)
    for (int i=0;i<MatSize;i++)
    {
      Az += A(i,2);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Bx)
    for (int i=0;i<MatSize;i++)
    {
      Bx += B(i,0);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:By)
    for (int i=0;i<MatSize;i++)
    {
      By += B(i,1);
    }
    #pragma omp for nowait schedule(dynamic) reduction(+:Bz)
    for (int i=0;i<MatSize;i++)
    {
      Bz += B(i,2);
    }
  }
  #pragma omp barrier
  //Take average
  Ax /= MatSize;
  Ay /= MatSize;
  Az /= MatSize;
  Bx /= MatSize;
  By /= MatSize;
  Bz /= MatSize;
  //Translate centroids
  #pragma omp parallel
  {
    //Move A and B to (0,0,0)
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<MatSize;i++)
    {
      A(i,0) -= Ax;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<MatSize;i++)
    {
      A(i,1) -= Ay;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<MatSize;i++)
    {
      A(i,2) -= Az;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<MatSize;i++)
    {
      B(i,0) -= Bx;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<MatSize;i++)
    {
      B(i,1) -= By;
    }
    #pragma omp for nowait schedule(dynamic)
    for (int i=0;i<MatSize;i++)
    {
      B(i,2) -= Bz;
    }
  }
  #pragma omp barrier
  //Calculate covariance matrix
  CoVar = (B.transpose())*A;
  //Compute SVD and identity matrix
  JacobiSVD<MatrixXd> SVDMat(CoVar,ComputeFullU|ComputeFullV);
  MatrixXd DetMat = SVDMat.matrixV()*(SVDMat.matrixU().transpose());
  double SignVal = DetMat.determinant();
  if (SignVal < 0)
  {
    SignVal = -1;
  }
  else
  {
    SignVal = 1;
  }
  MatrixXd Ident(3,3);
  Ident.setIdentity();
  Ident(2,2) *= SignVal; //Change sign for rotation
  //Find optimal rotation matrix
  RotMat = SVDMat.matrixV()*Ident*(SVDMat.matrixU().transpose());
  //Rotate matrix A
  B *= RotMat;
  //Return the modified positions
  return;
};

VectorXd KabschDisplacement(MatrixXd& A, MatrixXd& B, int MatSize)
{
  //Returns the distance between two superimposed structures
  VectorXd Dist(3*MatSize);
  //Rotate structures
  KabschRotation(A,B,MatSize);
  //Calculate displacement
  #pragma omp parallel for schedule(dynamic)
  for (int i=0;i<(3*MatSize);i++)
  {
    //Find the correct location in the arrays
    int Direc = i%3; //Find the remainder: x=0,y=1,z=2
    int AtID = (i-Direc)/3; //Find array index
    //Calculate displacement
    Dist(i) = A(AtID,Direc)-B(AtID,Direc);
  }
  //Return array
  return Dist;
};

//Physical property analysis functions
double LICHEMDensity(vector<QMMMAtom>& Struct, QMMMSettings& QMMMOpts)
{
  //Function to calculate the density for periodic calculations
  double rho = 0;
  //Sum the masses
  #pragma omp parallel for schedule(dynamic) reduction(+:rho)
  for (int i=0;i<Natoms;i++)
  {
    rho += Struct[i].m;
  }
  //Divide by the volume (SI)
  rho /= (Lx*Ly*Lz);
  //Change to g/cm^3 and return
  rho *= amu2kg; //amu to kg
  rho *= m2Ang*m2Ang*m2Ang/1000; //A^3 to cm^3
  return rho;
};

VectorXd LICHEMFreq(vector<QMMMAtom>& Struct, MatrixXd& QMMMHess,
                    QMMMSettings& QMMMOpts, int Bead, int& remct)
{
  //Function to perform a QMMM frequency analysis
  double ProjTol = 0.65; //Amount of overlap to remove a mode
  double ZeroTol = 1.00; //Smallest possible frequency (cm^-1)
  int transrotct = 0; //Number of deleted translation and rotational modes
  //Define variables
  int Ndof = 3*(Nqm+Npseudo); //Degrees of freedom
  //Define arrays
  EigenSolver<MatrixXd> FreqAnalysis;
  VectorXd QMMMFreqs(Ndof); //Vibrational frequencies (cm^-1)
  MatrixXd QMMMNormModes(Ndof,Ndof); //Normal modes
  MatrixXd FreqMatrix(Ndof,Ndof); //Diagonal frequency matrix
  VectorXd TransX(Ndof); //X translation
  VectorXd TransY(Ndof); //Y translation
  VectorXd TransZ(Ndof); //Z translation
  VectorXd RotX(Ndof); //X rotation
  VectorXd RotY(Ndof); //Y rotation
  VectorXd RotZ(Ndof); //Z rotation
  //Initialize arrays
  QMMMFreqs.setZero();
  QMMMNormModes.setZero();
  FreqMatrix.setZero();
  TransX.setZero();
  TransY.setZero();
  TransZ.setZero();
  RotX.setZero();
  RotY.setZero();
  RotZ.setZero();
  //Collect QM and PB masses
  vector<double> Masses;
  for (int i=0;i<Natoms;i++)
  {
    //Locate QM and PB atoms
    if (Struct[i].QMregion or Struct[i].PBregion)
    {
      //Switch to a.u. and save mass
      double massval = Struct[i].m/ElecMass;
      Masses.push_back(massval); //X component
      Masses.push_back(massval); //Y component
      Masses.push_back(massval); //Z component
    }
  }
  //Mass scale the Hessian matrix
  #pragma omp parallel for
  for (int i=0;i<Ndof;i++)
  {
    //Update diagonal matrix elements
    QMMMHess(i,i) /= Masses[i];
    //Update off-diagonal matrix elements
    for (int j=0;j<i;j++)
    {
      //Mass scale
      QMMMHess(i,j) /= sqrt(Masses[i]*Masses[j]);
      //Apply symmetry
      QMMMHess(j,i) = QMMMHess(i,j);
    }
  }
  //Diagonalize Hessian matrix
  FreqAnalysis.compute(QMMMHess);
  QMMMFreqs = FreqAnalysis.eigenvalues().real();
  QMMMNormModes = FreqAnalysis.eigenvectors().real();
  //Create model translation and rotation modes
  #pragma omp parallel for
  for (int i=0;i<(Nqm+Npseudo);i++)
  {
    //Translational modes
    TransX(3*i) = sqrt(Masses[3*i]);
    TransY(3*i+1) = sqrt(Masses[3*i+1]);
    TransZ(3*i+2) = sqrt(Masses[3*i+2]);
  }
  TransX.normalize();
  TransY.normalize();
  TransZ.normalize();
  //Remove translation and rotation
  transrotct = 0; //Use as a counter
  #pragma omp parallel for reduction(+:transrotct)
  for (int i=0;i<Ndof;i++)
  {
    bool IsTransRot = 0;
    double DotTest; //Saves overlap
    //Locate translations
    DotTest = abs(TransX.dot(QMMMNormModes.col(i)));
    if (DotTest > ProjTol)
    {
      IsTransRot = 1;
    }
    DotTest = abs(TransY.dot(QMMMNormModes.col(i)));
    if (DotTest > ProjTol)
    {
      IsTransRot = 1;
    }
    DotTest = abs(TransZ.dot(QMMMNormModes.col(i)));
    if (DotTest > ProjTol)
    {
      IsTransRot = 1;
    }
    //Locate rotations
    DotTest = abs(RotX.dot(QMMMNormModes.col(i)));
    if (DotTest > ProjTol)
    {
      IsTransRot = 1;
    }
    DotTest = abs(RotY.dot(QMMMNormModes.col(i)));
    if (DotTest > ProjTol)
    {
      IsTransRot = 1;
    }
    DotTest = abs(RotZ.dot(QMMMNormModes.col(i)));
    if (DotTest > ProjTol)
    {
      IsTransRot = 1;
    }
    //Locate small frequencies
    double SmallFreq; //Temporary storage
    SmallFreq = ZeroTol/Har2wavenum; //Convert tolerance to a.u.
    SmallFreq *= SmallFreq; //Square tolerance
    if (abs(QMMMFreqs(i)) < SmallFreq)
    {
      IsTransRot = 1;
    }
    //Adjust frequencies
    if (IsTransRot and (!QMMM))
    {
      //Remove frequency
      QMMMFreqs(i);
      FreqMatrix(i,i) = 0;
      transrotct += 1;
    }
    else
    {
      //Save frequency
      FreqMatrix(i,i) = QMMMFreqs(i);
    }
  }
  if (transrotct > 0)
  {
    //Remove unwanted frequencies
    QMMMHess = QMMMNormModes*FreqMatrix*QMMMNormModes.inverse();
    FreqAnalysis.compute(QMMMHess);
    QMMMFreqs = FreqAnalysis.eigenvalues().real();
    QMMMNormModes = FreqAnalysis.eigenvectors().real();
  }
  //Take the square root and keep the sign
  #pragma omp parallel for
  for (int i=0;i<Ndof;i++)
  {
    //Save sign
    int freqsign = 1;
    if (QMMMFreqs(i) < 0)
    {
      freqsign = -1;
    }
    //Remove sign
    QMMMFreqs(i) = abs(QMMMFreqs(i));
    //Take square root
    QMMMFreqs(i) = sqrt(QMMMFreqs(i));
    //Replace sign
    QMMMFreqs(i) *= freqsign;
  }
  //Change units
  QMMMFreqs *= Har2wavenum;
  //Remove negligible frequencies
  transrotct = 0; //Reset counter
  #pragma omp parallel for reduction(+:transrotct)
  for (int i=0;i<Ndof;i++)
  {
    //Delete frequencies below the tolerance
    if (abs(QMMMFreqs(i)) < ZeroTol)
    {
      transrotct += 1;
      QMMMFreqs(i) = 0;
    }
  }
  //Write all normal modes
  if ((!QMMM) and QMMMOpts.PrintNormModes)
  {
    WriteModes(Struct,0,QMMMFreqs,QMMMNormModes,QMMMOpts,Bead);
  }
  //Write modes for imaginary frequencies
  else if (QMMM)
  {
    WriteModes(Struct,1,QMMMFreqs,QMMMNormModes,QMMMOpts,Bead);
  }
  //Return frequencies
  remct = transrotct;
  return QMMMFreqs;
};

void WriteModes(vector<QMMMAtom>& Struct, bool ImagOnly, VectorXd& Freqs,
                MatrixXd& NormModes, QMMMSettings& QMMMOpts, int Bead)
{
  //Function to write normal modes
  int Nframes = 5; //Trajectory files have (2N+1) frames
  double AmpFrac = 0.50; //Fractional amplitude of the mode
  stringstream call; //Generic stream
  fstream modefile; //Normal mode output file
  int Ndof = 3*(Nqm+Npseudo); //Number of vibrational modes
  int ct; //Generic counter
  double CurAmp; //Stores the amplitude for each frame
  for (int i=0;i<Ndof;i++)
  {
    //Check print options
    if (((!ImagOnly) or (Freqs(i) < 0)) and (Freqs(i) != 0))
    {
      //Print normal mode
      call.str("");
      //Create file
      call << "NormModes_" << i << ".xyz";
      modefile.open(call.str().c_str(),ios_base::out);
      //Write file
      for (int j=0;j<(2*Nframes+1);j++)
      {
        //Loop over frames
        CurAmp = -1*AmpFrac; //Current amplitude
        CurAmp += (j*AmpFrac)/Nframes; //Move along the trajectory
        ct = 0; //Keeps track of the atom IDs
        modefile << (Nqm+Npseudo) << '\n' << '\n';
        for (int k=0;k<Natoms;k++)
        {
          if (Struct[k].QMregion or Struct[k].PBregion)
          {
            //Write element
            modefile << Struct[k].QMTyp << " ";
            //Write X component
            modefile << Struct[k].P[Bead].x+(CurAmp*NormModes(ct,i));
            modefile << " ";
            ct += 1;
            //Write Y component
            modefile << Struct[k].P[Bead].y+(CurAmp*NormModes(ct,i));
            modefile << " ";
            ct += 1;
            //Write Z component
            modefile << Struct[k].P[Bead].z+(CurAmp*NormModes(ct,i));
            modefile << '\n';
            ct += 1;
          }
        }
      }
      modefile.flush();
      modefile.close();
    }
  }
  return;
};

