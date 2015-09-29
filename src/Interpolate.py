#!/usr/bin/python

#Temporary script for interpolating reaction pathways

#Modules
import sys

#Check input and print help
if (len(sys.argv) != 4):
  line = '\n'
  line += "Usage: ./Inperpolate.py Npoints Reactant.xyz Product.xyz"
  line += '\n'
  print(line)
  exit(0)

#Set number of points
Npt = int(sys.argv[1])

#Open files
rfile = open(sys.argv[2],"r")
pfile = open(sys.argv[3],"r")
ofile = open("BeadStartStruct.xyz","w")

#Create structures
React = []
Prod = []
Typs = []

#Check number of atoms
Nr = rfile.readline()
Nr = Nr.strip()
Nr = int(Nr)
Np = pfile.readline() 
Np = Np.strip()
Np = int(Np)
if (Np != Nr):
  line = '\n'
  line += "Error: Number of atoms does not match!"
  line += '\n'
  print(line)
  exit(0)
else:
  Natoms = Nr

#Clear junk
rfile.readline()
pfile.readline()

#Read structures
for i in range(Natoms):
  tmp1 = []
  tmp2 = []
  dummy = rfile.readline()
  dummy = dummy.split()
  tmp1.append(dummy[0])
  tmp2.append(float(dummy[1]))
  tmp2.append(float(dummy[2]))
  tmp2.append(float(dummy[3]))
  React.append(tmp2)
  tmp2 = []
  dummy = pfile.readline()
  dummy = dummy.split()
  tmp1.append(dummy[0])
  tmp2.append(float(dummy[1]))
  tmp2.append(float(dummy[2]))
  tmp2.append(float(dummy[3]))
  Prod.append(tmp2)  
  Typs.append(tmp1)

#Check atom order
for i in range(Natoms):
  if (Typs[i][0] != Typs[i][1]):
    line = '\n'
    line += "Error: Reactant and product atoms are not in the same order!"
    line += '\n'
    print(line)
    exit(0)

#Print some output
line = '\n'
line += "Interpolating structure..."
print(line)

#Interpolate structure
Npts = float(Npt)
line = ""
line += `Natoms*Npt`
line += '\n'
line += '\n'
for i in range(Natoms):
  for j in range(Npt):
    attyp = Typs[i][0]
    x = React[i][0]
    x += (j/(Npts-1))*Prod[i][0]
    x -= (j/(Npts-1))*React[i][0]
    y = React[i][1]
    y += (j/(Npts-1))*Prod[i][1]
    y -= (j/(Npts-1))*React[i][1]
    z = React[i][2]
    z += (j/(Npts-1))*Prod[i][2]
    z -= (j/(Npts-1))*React[i][2]
    line += attyp+" "
    line += `x`+" "
    line += `y`+" "
    line += `z`+'\n'
ofile.write(line)

#Print some output
line = '\n'
line += "Done."
line += '\n'
print(line)

#Quit
ofile.close()
exit(0)
