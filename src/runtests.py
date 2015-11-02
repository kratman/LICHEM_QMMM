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
line += "***************************************************"
line += '\n'
line += "*                                                 *"
line += '\n'
line += "*   LICHEM: Layered Interacting CHEmical Models   *"
line += '\n'
line += "*                                                 *"
line += '\n'
line += "*        Symbiotic Computational Chemistry        *"
line += '\n'
line += "*                                                 *"
line += '\n'
line += "***************************************************"
line += '\n'
print(line)

#Read arguments
DryRun = 0 #Only check packages
if (len(sys.argv) < 4):
  #Print help if arguments are missing
  line = ""
  line += "Usage:"
  line += '\n'
  line += " user:$ ./runtests Ncpus All"
  line += '\n'
  line += "  or "
  line += '\n'
  line += " user:$ ./runtests Ncpus QMPackage MMPackage"
  line += '\n'
  line += "  or "
  line += '\n'
  line += " user:$ ./runtests Ncpus QMPackage MMPackage dry"
  line += '\n'
  line += '\n'
  #Find LICHEM
  line += "LICHEM binary: "
  cmd = "which lichem"
  try:
    #Find path
    LICHEMbin = subprocess.check_output(cmd,shell=True)
    LICHEMbin = ClrSet.TPass+LICHEMbin.strip()+ClrSet.Reset
  except:
    LICHEMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += LICHEMbin
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
    QMbin = ClrSet.TPass+QMbin.strip()+ClrSet.Reset
  except:
    QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += QMbin
  line += '\n'
  line += " Gaussian: "
  cmd = "which g09"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = ClrSet.TPass+QMbin.strip()+ClrSet.Reset
  except:
    QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += QMbin
  line += '\n'
  line += " NWChem: "
  cmd = "which nwchem"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = ClrSet.TPass+QMbin.strip()+ClrSet.Reset
  except:
    QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
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
    MMbin = ClrSet.TPass+MMbin.strip()+ClrSet.Reset
  except:
    MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += MMbin
  line += '\n'
  line += " LAMMPS: "
  cmd = "which lammps"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = ClrSet.TPass+MMbin.strip()+ClrSet.Reset
  except:
    MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += MMbin
  line += '\n'
  line += " AMBER: "
  cmd = "which pmemd"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = ClrSet.TPass+MMbin.strip()+ClrSet.Reset
  except:
    MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += MMbin
  line += '\n'
  print(line)
  exit(0)

line = "Running LICHEM tests..."
line += '\n'
print(line)

Ncpus = int(sys.argv[1]) #Threads
QMPack = sys.argv[2] #QM wrapper for calculations
QMPack = QMPack.lower()
MMPack = sys.argv[3] #MM wrapper for calculations
MMPack = MMPack.lower()
if (len(sys.argv) > 4):
  if ((sys.argv[4] == "Dry") or (sys.argv[4] == "dry")):
    #Quit early
    DryRun = 1

LICHEMbin = ""
QMbin = ""
MMbin = ""
#Check packages and identify missing binaries
BadLICHEM = 0
cmd = "which lichem"
try:
  #Find path
  LICHEMbin = subprocess.check_output(cmd,shell=True)
  LICHEMbin = LICHEMbin.strip()
except:
  LICHEMbin = "N/A"
  BadLICHEM = 1
if (BadLICHEM == 1):
  #Quit with an error
  line = ""
  line += "Error: LICHEM binary not found!"
  line += '\n'
  print(line)
  exit(0)
BadQM = 1
if ((QMPack == "psi4") or (QMPack == "psi")):
  QMPack = "PSI4"
  cmd = "which psi4"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  BadQM = 0
if ((QMPack == "gaussian") or (QMPack == "g09")):
  QMPack = "Gaussian"
  cmd = "which g09"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = QMbin.strip()
  except:
    QMbin = "N/A"
  BadQM = 0
if (QMPack == "nwchem"):
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
if (MMPack == "tinker"):
  MMPack = "TINKER"
  cmd = "which analyze"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  BadMM = 0
if (MMPack == "lammps"):
  MMPack = "LAMMPS"
  cmd = "which lammps"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = MMbin.strip()
  except:
    MMbin = "N/A"
  BadMM = 0
if (MMPack == "amber"):
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
line += " LICHEM binary: "
line += LICHEMbin
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

#Set path based on packages
DirPath = ""
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

#Start printing results
print("Results:")

#Check QM energy
line = ""
PassEnergy = 0
cmd = "lichem -n "
cmd += `Ncpus`
cmd += " "
cmd += "-x xyzfile.xyz "
cmd += "-r qmreg.inp "
cmd += "-c connect.inp "
cmd += "-o trash.xyz "
cmd += "> tests.out " #Capture stdout
cmd += "2>&1" #Capture stderr
subprocess.call(cmd,shell=True) #Run calculations
cmd = ""
cmd += "grep -e"
cmd += ' "QM energy: " '
cmd += "tests.out"
try:
  QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
  QMMMEnergy = QMMMEnergy.split()
  QMMMEnergy = float(QMMMEnergy[2])
  QMMMEnergy = round(QMMMEnergy,5)
except:
  QMMMEnergy = 0.0
if (QMPack == "PSI4"):
  #Check against saved energy
  if (QMMMEnergy == round(-4155.414685776608,5)):
    PassEnergy = 1
if (QMPack == "Gaussian"):
  #Check against saved energy
  if (QMMMEnergy == round(-4155.413972324945,5)):
    PassEnergy = 1
if (QMPack == "NWChem"):
  #Check against saved energy
  if (QMMMEnergy == round(-4155.414734834268,5)):
    PassEnergy = 1
line += " QM energy: "
if (PassEnergy == 1):
  line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
else:
  line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
cmd = ""
cmd += "grep -e"
cmd += ' "Total wall time: " '
cmd += "tests.out"
try:
  RunTime = subprocess.check_output(cmd,shell=True) #Get run time
  RunTime = RunTime.split()
  RunTime = " "+RunTime[3]+" "+RunTime[4]
except:
  RunTime = " N/A"
line += RunTime
print(line)

#Clean up files
cmd = ""
cmd += "rm -f tinker.key tests.out trash.xyz"
if (QMPack == "Gaussian"):
  #Remove checkpoint files
  cmd += " *.chk"
if (QMPack == "PSI4"):
  #Remove checkpoint files
  cmd += " timer.* psi.*"
subprocess.call(cmd,shell=True)

#TINKER tests

#Check MM energy
line = ""
PassEnergy = 0
cmd = "cp pchrg.key tinker.key"
subprocess.call(cmd,shell=True) #Copy key file
cmd = "lichem -n "
cmd += `Ncpus`
cmd += " "
cmd += "-x xyzfile.xyz "
cmd += "-r mmreg.inp "
cmd += "-c connect.inp "
cmd += "-o trash.xyz "
cmd += "> tests.out " #Capture stdout
cmd += "2>&1" #Capture stderr
subprocess.call(cmd,shell=True) #Run calculations
cmd = ""
cmd += "grep -e"
cmd += ' "MM energy: " '
cmd += "tests.out"
try:
  QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
  QMMMEnergy = QMMMEnergy.split()
  QMMMEnergy = float(QMMMEnergy[2])
  QMMMEnergy = round(QMMMEnergy,5)
except:
  QMMMEnergy = 0.0
#Check against saved energy
if (QMMMEnergy == round(-0.259690353622,5)):
  PassEnergy = 1
line += " MM energy: "
if (PassEnergy == 1):
  line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
else:
  line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
cmd = ""
cmd += "grep -e"
cmd += ' "Total wall time: " '
cmd += "tests.out"
try:
  RunTime = subprocess.check_output(cmd,shell=True) #Get run time
  RunTime = RunTime.split()
  RunTime = " "+RunTime[3]+" "+RunTime[4]
except:
  RunTime = " N/A"
line += RunTime
print(line)

#Clean up files
cmd = ""
cmd += "rm -f tinker.key tests.out trash.xyz"
if (QMPack == "Gaussian"):
  #Remove checkpoint files
  cmd += " *.chk"
if (QMPack == "PSI4"):
  #Remove checkpoint files
  cmd += " timer.* psi.*"
subprocess.call(cmd,shell=True)

#Check QMMM point-charge energy results
line = ""
PassEnergy = 0
cmd = "cp pchrg.key tinker.key"
subprocess.call(cmd,shell=True) #Copy key file
cmd = "lichem -n "
cmd += `Ncpus`
cmd += " "
cmd += "-x xyzfile.xyz "
cmd += "-r pchrgreg.inp "
cmd += "-c connect.inp "
cmd += "-o trash.xyz "
cmd += "> tests.out " #Capture stdout
cmd += "2>&1" #Capture stderr
subprocess.call(cmd,shell=True) #Run calculations
cmd = ""
cmd += "grep -e"
cmd += ' "QMMM energy: " '
cmd += "tests.out"
try:
  QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
  QMMMEnergy = QMMMEnergy.split()
  QMMMEnergy = float(QMMMEnergy[2])
  QMMMEnergy = round(QMMMEnergy,5)
except:
  QMMMEnergy = 0.0
if (QMPack == "PSI4"):
  #Check against saved energy
  if (QMMMEnergy == round(-2077.868475298003,5)):
    PassEnergy = 1
if (QMPack == "Gaussian"):
  #Check against saved energy
  if (QMMMEnergy == round(-2077.868129577162,5)):
    PassEnergy = 1
if (QMPack == "NWChem"):
  #Check against saved energy
  if (QMMMEnergy == round(-2077.868519717375,5)):
    PassEnergy = 1
line += " QMMM energy: "
if (PassEnergy == 1):
  line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
else:
  line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
cmd = ""
cmd += "grep -e"
cmd += ' "Total wall time: " '
cmd += "tests.out"
try:
  RunTime = subprocess.check_output(cmd,shell=True) #Get run time
  RunTime = RunTime.split()
  RunTime = " "+RunTime[3]+" "+RunTime[4]
except:
  RunTime = " N/A"
line += RunTime
print(line)

#Clean up files
cmd = ""
cmd += "rm -f tinker.key tests.out trash.xyz"
if (QMPack == "Gaussian"):
  #Remove checkpoint files
  cmd += " *.chk"
if (QMPack == "PSI4"):
  #Remove checkpoint files
  cmd += " timer.* psi.*"
subprocess.call(cmd,shell=True)

#Check QMMM polarizable energy results
line = ""
PassEnergy = 0
cmd = "cp pol.key tinker.key"
subprocess.call(cmd,shell=True) #Copy key file
cmd = "lichem -n "
cmd += `Ncpus`
cmd += " "
cmd += "-x xyzfile.xyz "
cmd += "-r polreg.inp "
cmd += "-c connect.inp "
cmd += "-o trash.xyz "
cmd += "> tests.out " #Capture stdout
cmd += "2>&1" #Capture stderr
subprocess.call(cmd,shell=True) #Run calculations
cmd = ""
cmd += "grep -e"
cmd += ' "QMMM energy: " '
cmd += "tests.out"
try:
  QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
  QMMMEnergy = QMMMEnergy.split()
  QMMMEnergy = float(QMMMEnergy[2])
  QMMMEnergy = round(QMMMEnergy,5)
except:
  QMMMEnergy = 0.0
if (QMPack == "PSI4"):
  #Check against saved energy
  if (QMMMEnergy == round(-2077.771998494848,5)):
    PassEnergy = 1
if (QMPack == "Gaussian"):
  #Check against saved energy
  if (QMMMEnergy == round(-2077.775628449552,5)):
    PassEnergy = 1
if (QMPack == "NWChem"):
  #Check against saved energy
  if (QMMMEnergy == round(-2077.775532434397,5)):
    PassEnergy = 1
line += " Pol. QMMM energy: "
if (PassEnergy == 1):
  line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
else:
  line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
cmd = ""
cmd += "grep -e"
cmd += ' "Total wall time: " '
cmd += "tests.out"
try:
  RunTime = subprocess.check_output(cmd,shell=True) #Get run time
  RunTime = RunTime.split()
  RunTime = " "+RunTime[3]+" "+RunTime[4]
except:
  RunTime = " N/A"
line += RunTime
line += '\n'

#Clean up files
cmd = ""
cmd += "rm -f tinker.key tests.out trash.xyz"
if (QMPack == "Gaussian"):
  #Remove checkpoint files
  cmd += " *.chk"
if (QMPack == "PSI4"):
  #Remove checkpoint files
  cmd += " timer.* psi.*"
subprocess.call(cmd,shell=True)

#Print all results
print(line)

#Quit
os.chdir("../")
line = "Done."
line += '\n'
print(line)
exit(0)
