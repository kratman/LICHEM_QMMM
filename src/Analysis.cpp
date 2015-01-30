

//Trajectory analysis functions
void BurstTraj(vector<QMMMAtom>& Struct, string& filename,
     QMMMSettings& QMMMOpts)
{
  //Function to split reaction path and path-integral trajectory frames
  int ct; //Generic counter
  stringstream call;
  fstream burstfile;
  string dummy; //Generic string
  //Open new split trajectory file
  call.str("");
  call << "BurstStruct.xyz";
  ct = 1; //Start counting at the second file
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
  for (int j=0;j<QMMMOpts.Nbeads;j++)
  {
    burstfile << Natoms; //Number of atoms
    burstfile << '\n' << '\n';
    for (int i=0;i<Natoms;i++)
    {
      burstfile << Struct[i].QMTyp << " ";
      burstfile << Struct[i].P[j].x << " ";
      burstfile << Struct[i].P[j].y << " ";
      burstfile << Struct[i].P[j].z << '\n';
    }
  }
  return;
};
