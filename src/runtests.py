###################################################
#                                                 #
#   LICHEM: Layered Interacting CHEmical Models   #
#                                                 #
#        Symbiotic Computational Chemistry        #
#                                                 #
###################################################

#LICHEM semi-automated test suite 

#Modules
import subprocess
import sys
import os

#Classes
class ClrSet:
  #Unicode colors
  Norm = '\033[0m'
  Red = '\033[91m'
  Bold = '\033[1m'
  Green = '\033[92m'
  Blue = '\033[94m'
  Cyan = '\033[36m'
  #Set colors
  TFail = Bold+Red #Highlight failed tests
  TPass = Bold+Green #Highlight passed tests
  Reset = Norm #Reset to defaults

#Print title
line = '\n'
line += "###################################################"
line += '\n'
line += "#                                                 #"
line += '\n'
line += "#   LICHEM: Layered Interacting CHEmical Models   #"
line += '\n'
line += "#                                                 #"
line += '\n'
line += "#        Symbiotic Computational Chemistry        #"
line += '\n'
line += "#                                                 #"
line += '\n'
line += "###################################################"
line += '\n'
line += '\n'
line += "Running LICHEM tests..."
line += '\n'
print(line)

#Read arguments
DryRun = 0 #Only check packages
if (len(sys.argv) < 4):
  #Print help if arguments are missing
  line = ""
  line += "Usage: python runtests.py Ncpus QMPackage MMPackage"
  line += '\n'
  line += '\n'
  #Identify QM wrappers
  line += "Available QM wrappers:"
  line += '\n'
  line += " PSI4: "
  cmd = "which psi4"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  line += QMbin
  line += '\n'
  line += " Gaussian: "
  cmd = "which g09"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  line += QMbin
  line += '\n'
  line += " NWChem: "
  cmd = "which nwchem"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  line += QMbin
  line += '\n'
  line += '\n'
  #Identify MM wrappers
  line += "Available MM wrappers:"
  line += '\n'
  line += " TINKER: "
  cmd = "which analyze"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  line += MMbin
  line += '\n'
  line += " LAMMPS: "
  cmd = "which lammps"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  line += MMbin
  line += '\n'
  line += " AMBER: "
  cmd = "which pmemd"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  line += MMbin
  line += '\n'
  print(line)
  exit(0)
Ncpus = int(sys.argv[1]) #Threads
QMPack = sys.argv[2] #QM wrapper for calculations
MMPack = sys.argv[3] #MM wrapper for calculations
if (len(sys.argv) > 4):
  if ((sys.argv[4] == "Dry") or (sys.argv[4] == "dry")):
    #Quit early
    DryRun = 1
QMbin = ""
MMbin = ""

#Check packages
BadQM = 1
if ((QMPack == "PSI4") or (QMPack == "psi4") or (QMPack == "Psi4")):
  QMPack = "PSI4"
  cmd = "which psi4"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  BadQM = 0
if ((QMPack == "Gaussian") or (QMPack == "gaussian") or (QMPack == "g09")):
  QMPack = "Gaussian"
  cmd = "which g09"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  BadQM = 0
if ((QMPack == "NWChem") or (QMPack == "nwchem") or (QMPack == "NWCHEM")):
  QMPack = "NWChem"
  cmd = "which nwchem"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  BadQM = 0
if (BadQM == 1):
  #Quit with an error
  line = '\n'
  line += "Error: QM package name '"
  line += QMPack
  line += "' not recognized." 
  line += '\n'
  print(line)
  exit(0)
BadMM = 1
if ((MMPack == "TINKER") or (MMPack == "tinker") or (MMPack == "Tinker")):
  MMPack = "TINKER"
  cmd = "which analyze"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  BadMM = 0
if ((MMPack == "LAMMPS") or (MMPack == "lammps") or (MMPack == "Lammps")):
  MMPack = "LAMMPS"
  cmd = "which lammps"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  BadMM = 0
if ((MMPack == "AMBER") or (MMPack == "amber") or (MMPack == "Amber")):
  MMPack = "AMBER"
  cmd = "which pmemd" #Check
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  BadMM = 0
if (BadMM == 1):
  #Quit with error
  line = '\n'
  line += "Error: MM package name '"
  line += MMPack
  line += "' not recognized."
  line += '\n'
  print(line)
  exit(0)

#Print settings
line = "Settings:"
line += '\n'
line += " Threads: "
line += `Ncpus`
line += '\n'
line += " QM package: "
line += QMPack
line += '\n'
line += " Binary: "
line += QMbin
line += '\n'
line += " MM package: "
line += MMPack
line += '\n'
line += " Binary: "
line += MMbin
line += '\n'
print(line)

#Escape for dry runs
if (DryRun == 1):
  #Quit without an error
  line = "Dry run completed."
  line += '\n'
  print(line)
  exit(0)

#Escape if binaries not found
if ((QMbin == "N/A") or (MMbin == "N/A")):
  #Quit with an error
  line = "Error: Missing binaries."
  line += '\n'
  print(line)
  exit(0)

#Run tests
DirPath = "./"
cmd = "lichem -n "
cmd += `Ncpus`
cmd += " "

#Set path based on packages
if (QMPack == "PSI4"):
  DirPath += "PSI4_"
if (QMPack == "Gaussian"):
  DirPath += "Gau_"
if (QMPack == "NWChem"):
  DirPath += "NWChem_"
DirPath += MMPack
DirPath += "/"

#Change directory
os.chdir(DirPath)

#Complete commands
cmd += "-x xyzfile.xyz "
cmd += "-r regions.inp "
cmd += "-c connect.inp "
cmd += "-o trash.xyz "
cmd += "> tests.out"

#Run calculations
subprocess.call(cmd,shell=True)

#Start printing results
line = "Results:"
line += '\n'

#Check QMMM energy results
PassEnergy = 0
cmd = ""
cmd += "grep -e"
cmd += ' "QMMM energy: " '
cmd += "tests.out"
QMMMEnergy = subprocess.check_output(cmd,shell=True)
QMMMEnergy = QMMMEnergy.split()
QMMMEnergy = float(QMMMEnergy[2])
QMMMEnergy = round(QMMMEnergy,6)
if (QMPack == "PSI4"):
  #Check again saved energy
  if (QMMMEnergy == round(-2077.771998247412,6)):
    PassEnergy = 1
if (QMPack == "Gaussian"):
  #Check again saved energy
  if (QMMMEnergy == round(-2077.775625181200,6)):
    PassEnergy = 1
if (QMPack == "NWChem"):
  #Check again saved energy
  if (QMMMEnergy == round(-2077.775529980780,6)):
    PassEnergy = 1
line += " QMMM energy (water dimer): "
if (PassEnergy == 1):
  line += ClrSet.TPass+"Pass"+ClrSet.Reset
  line += '\n'
else:
  line += ClrSet.TFail+"Fail"+ClrSet.Reset
  line += '\n'

#Print all results
print(line)

#Clean up files
cmd = ""
cmd += "rm -f trash.xyz tests.out"
if (QMPack == "Gaussian"):
  #Remove checkpoint files
  cmd += " *.chk"
subprocess.call(cmd,shell=True)

#Quit
os.chdir("../")
line = '\n'
line += "Done."
line += '\n'
print(line)
exit(0)
