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
import time
import sys
import os

#Development settings
#NB: These settings should not be modified
UpdateResults = 0 #Flag to print energies to update tests

#Start timer and counters immediately
StartTime = time.time()
passct = 0
failct = 0

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
AllTests = 0 #Run all tests at once
if (len(sys.argv) == 3):
  if ((sys.argv[2] == "All") or (sys.argv[2] == "all")):
    #Automatically run all tests
    AllTests = 1
if (len(sys.argv) < 4):
  line = ""
  if (AllTests == 0):
    #Print help if arguments are missing
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
  if (AllTests == 0):
    #Quit
    exit(0)

Ncpus = int(sys.argv[1]) #Threads
if (AllTests == 0):
  QMPack = sys.argv[2] #QM wrapper for calculations
  QMPack = QMPack.lower()
  MMPack = sys.argv[3] #MM wrapper for calculations
  MMPack = MMPack.lower()
  if (len(sys.argv) > 4):
    if ((sys.argv[4] == "Dry") or (sys.argv[4] == "dry")):
      #Quit early
      DryRun = 1

#Initialize variables
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
if (AllTests == 0):
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
if (AllTests == 0):
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
else:
  line += " Mode: All tests"
  line += '\n'
if (DryRun == 1):
  line += " Mode: Dry run"
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
if (((QMbin == "N/A") or (MMbin == "N/A")) and (AllTests == 0)):
  #Quit with an error
  line = "Error: Missing binaries."
  line += '\n'
  print(line)
  exit(0)

line = "***************************************************"
line += '\n'
line += '\n'
line += "Running LICHEM tests..."
line += '\n'
print(line)

#Make a list of tests
QMTests = []
MMTests = []
if (AllTests == 1):
  QMTests.append("PSI4")
  QMTests.append("Gaussian")
  QMTests.append("NWChem")
  MMTests.append("TINKER")
  MMTests.append("LAMMPS")
  MMTests.append("AMBER")
else:
  QMTests.append(QMPack)
  MMTests.append(MMPack)

#Loop over tests
for qmtest in QMTests:
  for mmtest in MMTests:
    #Set packages
    QMPack = qmtest
    MMPack = mmtest

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
    line = QMPack+"/"+MMPack
    line += " results:"
    print(line)

    #Delete line to avoid bugs
    line = ""

    #Check HF energy
    if ((QMPack == "PSI4") or (QMPack == "Gaussian")):
      PassEnergy = 0
      cmd = "lichem -n "
      cmd += `Ncpus`
      cmd += " "
      cmd += "-x waterdimer.xyz "
      cmd += "-r hfreg.inp "
      cmd += "-c watercon.inp "
      cmd += "-o trash.xyz "
      cmd += "> tests.out " #Capture stdout
      cmd += "2>&1" #Capture stderr
      subprocess.call(cmd,shell=True) #Run calculations
      cmd = ""
      cmd += "grep -e"
      cmd += ' "QM energy: " ' #Find final energy
      cmd += "tests.out"
      SavedEnergy = "Crashed..."
      try:
        #Safely check energy
        QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
        QMMMEnergy = QMMMEnergy.split()
        QMMMEnergy = float(QMMMEnergy[2])
        SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
        QMMMEnergy = round(QMMMEnergy,5)
      except:
        #Calculation failed
        QMMMEnergy = 0.0
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-4136.9303981392,5)):
          PassEnergy = 1
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-4136.9317704519,5)):
          PassEnergy = 1
      line += " HF energy:           "
      if (PassEnergy == 1):
        line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
        passct += 1
      else:
        line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
        failct += 1
      cmd = ""
      cmd += "grep -e"
      cmd += ' "Total wall time: " ' #Find run time
      cmd += "tests.out"
      try:
        RunTime = subprocess.check_output(cmd,shell=True) #Get run time
        RunTime = RunTime.split()
        RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
      except:
        RunTime = " N/A"
      line += RunTime
      if (UpdateResults == 1):
        line += ", "
        line += SavedEnergy
      print(line)

      #Clean up files
      cmd = ""
      cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
      cmd += " BeadStartStruct.xyz BurstStruct.xyz"
      if (QMPack == "Gaussian"):
        #Remove checkpoint files
        cmd += " *.chk"
      if (QMPack == "PSI4"):
        #Remove checkpoint files
        cmd += " timer.* psi.* *.32 *.180"
      if (QMPack == "NWChem"):
        #Remove checkpoint files
        cmd += " *.movecs"
      subprocess.call(cmd,shell=True)

      #Delete line to avoid bugs
      line = ""

    #Check DFT energy
    PassEnergy = 0
    cmd = "lichem -n "
    cmd += `Ncpus`
    cmd += " "
    cmd += "-x waterdimer.xyz "
    cmd += "-r pbereg.inp "
    cmd += "-c watercon.inp "
    cmd += "-o trash.xyz "
    cmd += "> tests.out " #Capture stdout
    cmd += "2>&1" #Capture stderr
    subprocess.call(cmd,shell=True) #Run calculations
    cmd = ""
    cmd += "grep -e"
    cmd += ' "QM energy: " ' #Find final energy
    cmd += "tests.out"
    SavedEnergy = "Crashed..."
    try:
      #Safely check energy
      QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
      QMMMEnergy = QMMMEnergy.split()
      QMMMEnergy = float(QMMMEnergy[2])
      SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
      QMMMEnergy = round(QMMMEnergy,5)
    except:
      #Calculation failed
      QMMMEnergy = 0.0
    if (QMPack == "PSI4"):
      #Check against saved energy
      if (QMMMEnergy == round(-4154.1683659877,5)):
        PassEnergy = 1
    if (QMPack == "Gaussian"):
      #Check against saved energy
      if (QMMMEnergy == round(-4154.1676114324,5)):
        PassEnergy = 1
    if (QMPack == "NWChem"):
      #Check against saved energy
      if (QMMMEnergy == round(-4154.1683939169,5)):
        PassEnergy = 1
    line += " PBE0 energy:         "
    if (PassEnergy == 1):
      line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
      passct += 1
    else:
      line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
      failct += 1
    cmd = ""
    cmd += "grep -e"
    cmd += ' "Total wall time: " ' #Find run time
    cmd += "tests.out"
    try:
      RunTime = subprocess.check_output(cmd,shell=True) #Get run time
      RunTime = RunTime.split()
      RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
    except:
      RunTime = " N/A"
    line += RunTime
    if (UpdateResults == 1):
      line += ", "
      line += SavedEnergy
    print(line)

    #Clean up files
    cmd = ""
    cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
    cmd += " BeadStartStruct.xyz BurstStruct.xyz"
    if (QMPack == "Gaussian"):
      #Remove checkpoint files
      cmd += " *.chk"
    if (QMPack == "PSI4"):
      #Remove checkpoint files
      cmd += " timer.* psi.* *.32 *.180"
    if (QMPack == "NWChem"):
      #Remove checkpoint files
      cmd += " *.movecs"
    subprocess.call(cmd,shell=True)

    #Delete line to avoid bugs
    line = ""

    #Check CCSD energy
    if (QMPack == "PSI4"):
      PassEnergy = 0
      cmd = "lichem -n "
      cmd += `Ncpus`
      cmd += " "
      cmd += "-x waterdimer.xyz "
      cmd += "-r ccsdreg.inp "
      cmd += "-c watercon.inp "
      cmd += "-o trash.xyz "
      cmd += "> tests.out " #Capture stdout
      cmd += "2>&1" #Capture stderr
      subprocess.call(cmd,shell=True) #Run calculations
      cmd = ""
      cmd += "grep -e"
      cmd += ' "QM energy: " ' #Find final energy
      cmd += "tests.out"
      SavedEnergy = "Crashed..."
      try:
        #Safely check energy
        QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
        QMMMEnergy = QMMMEnergy.split()
        QMMMEnergy = float(QMMMEnergy[2])
        SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
        QMMMEnergy = round(QMMMEnergy,5)
      except:
        #Calculation failed
        QMMMEnergy = 0.0
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-4147.730483706,5)):
          PassEnergy = 1
      line += " CCSD energy:         "
      if (PassEnergy == 1):
        line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
        passct += 1
      else:
        line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
        failct += 1
      cmd = ""
      cmd += "grep -e"
      cmd += ' "Total wall time: " ' #Find run time
      cmd += "tests.out"
      try:
        RunTime = subprocess.check_output(cmd,shell=True) #Get run time
        RunTime = RunTime.split()
        RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
      except:
        RunTime = " N/A"
      line += RunTime
      if (UpdateResults == 1):
        line += ", "
        line += SavedEnergy
      print(line)

      #Clean up files
      cmd = ""
      cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
      cmd += " BeadStartStruct.xyz BurstStruct.xyz"
      if (QMPack == "Gaussian"):
        #Remove checkpoint files
        cmd += " *.chk"
      if (QMPack == "PSI4"):
        #Remove checkpoint files
        cmd += " timer.* psi.* *.32 *.180"
      if (QMPack == "NWChem"):
        #Remove checkpoint files
        cmd += " *.movecs"
      subprocess.call(cmd,shell=True)

      #Delete line to avoid bugs
      line = ""

    #Check PM6 energy
    if (QMPack == "Gaussian"):
      PassEnergy = 0
      cmd = "lichem -n "
      cmd += `Ncpus`
      cmd += " "
      cmd += "-x waterdimer.xyz "
      cmd += "-r pm6reg.inp "
      cmd += "-c watercon.inp "
      cmd += "-o trash.xyz "
      cmd += "> tests.out " #Capture stdout
      cmd += "2>&1" #Capture stderr
      subprocess.call(cmd,shell=True) #Run calculations
      cmd = ""
      cmd += "grep -e"
      cmd += ' "QM energy: " ' #Find final energy
      cmd += "tests.out"
      SavedEnergy = "Crashed..."
      try:
        #Safely check energy
        QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
        QMMMEnergy = QMMMEnergy.split()
        QMMMEnergy = float(QMMMEnergy[2])
        SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
        QMMMEnergy = round(QMMMEnergy,5)
      except:
        #Calculation failed
        QMMMEnergy = 0.0
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-4.8623027634995,5)):
          PassEnergy = 1
      line += " PM6 energy:          "
      if (PassEnergy == 1):
        line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
        passct += 1
      else:
        line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
        failct += 1
      cmd = ""
      cmd += "grep -e"
      cmd += ' "Total wall time: " ' #Find run time
      cmd += "tests.out"
      try:
        RunTime = subprocess.check_output(cmd,shell=True) #Get run time
        RunTime = RunTime.split()
        RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
      except:
        RunTime = " N/A"
      line += RunTime
      if (UpdateResults == 1):
        line += ", "
        line += SavedEnergy
      print(line)

      #Clean up files
      cmd = ""
      cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
      cmd += " BeadStartStruct.xyz BurstStruct.xyz"
      if (QMPack == "Gaussian"):
        #Remove checkpoint files
        cmd += " *.chk"
      if (QMPack == "PSI4"):
        #Remove checkpoint files
        cmd += " timer.* psi.* *.32 *.180"
      if (QMPack == "NWChem"):
        #Remove checkpoint files
        cmd += " *.movecs"
      subprocess.call(cmd,shell=True)

      #Delete line to avoid bugs
      line = ""

    #Check NEB optimization
    PassEnergy = 0
    cmd = "cp methflbeads.xyz BeadStartStruct.xyz"
    subprocess.call(cmd,shell=True) #Copy restart file
    cmd = "lichem -n "
    cmd += `Ncpus`
    cmd += " "
    cmd += "-x methfluor.xyz "
    cmd += "-r nebreg.inp "
    cmd += "-c methflcon.inp "
    cmd += "-o trash.xyz "
    cmd += "> tests.out " #Capture stdout
    cmd += "2>&1" #Capture stderr
    subprocess.call(cmd,shell=True) #Run calculations
    cmd = ""
    cmd += "grep -e"
    cmd += ' "Opt. step: 2 " ' #Find final energy
    cmd += "tests.out"
    SavedEnergy = "Crashed..."
    try:
      #Safely check energy
      QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
      QMMMEnergy = QMMMEnergy.split()
      QMMMEnergy = float(QMMMEnergy[11])
      SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
      QMMMEnergy = round(QMMMEnergy,5)
    except:
      #Calculation failed
      QMMMEnergy = 0.0
    if (QMPack == "PSI4"):
      #Check against saved energy
      if (QMMMEnergy == round(-6511.0580192214,5)):
        PassEnergy = 1
    if (QMPack == "Gaussian"):
      #Check against saved energy
      if (QMMMEnergy == round(-6511.0567955964,5)):
        PassEnergy = 1
    if (QMPack == "NWChem"):
      #Check against saved energy
      if (QMMMEnergy == round(-6511.0579547077,5)):
        PassEnergy = 1
    line += " NEB TS energy:       "
    if (PassEnergy == 1):
      line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
      passct += 1
    else:
      line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
      failct += 1
    cmd = ""
    cmd += "grep -e"
    cmd += ' "Total wall time: " ' #Find run time
    cmd += "tests.out"
    try:
      RunTime = subprocess.check_output(cmd,shell=True) #Get run time
      RunTime = RunTime.split()
      RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
    except:
      RunTime = " N/A"
    line += RunTime
    if (UpdateResults == 1):
      line += ", "
      line += SavedEnergy
    print(line)

    #Clean up files
    cmd = ""
    cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
    cmd += " BeadStartStruct.xyz BurstStruct.xyz"
    if (QMPack == "Gaussian"):
      #Remove checkpoint files
      cmd += " *.chk"
    if (QMPack == "PSI4"):
      #Remove checkpoint files
      cmd += " timer.* psi.* *.32 *.180"
    if (QMPack == "NWChem"):
      #Remove checkpoint files
      cmd += " *.movecs"
    subprocess.call(cmd,shell=True)

    #Delete line to avoid bugs
    line = ""

    #TINKER tests
    if (MMPack == "TINKER"):
      #Check MM energy
      line = ""
      PassEnergy = 0
      cmd = "cp pchrg.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      cmd = "lichem -n "
      cmd += `Ncpus`
      cmd += " "
      cmd += "-x waterdimer.xyz "
      cmd += "-r mmreg.inp "
      cmd += "-c watercon.inp "
      cmd += "-o trash.xyz "
      cmd += "> tests.out " #Capture stdout
      cmd += "2>&1" #Capture stderr
      subprocess.call(cmd,shell=True) #Run calculations
      cmd = ""
      cmd += "grep -e"
      cmd += ' "MM energy: " ' #Find final energy
      cmd += "tests.out"
      SavedEnergy = "Crashed..."
      try:
        #Safely check energy
        QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
        QMMMEnergy = QMMMEnergy.split()
        QMMMEnergy = float(QMMMEnergy[2])
        SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
        QMMMEnergy = round(QMMMEnergy,5)
      except:
        #Calculation failed
        QMMMEnergy = 0.0
      #Check against saved energy
      if (QMMMEnergy == round(-0.2596903536223,5)):
        PassEnergy = 1
      line += " TIP3P energy:        "
      if (PassEnergy == 1):
        line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
        passct += 1
      else:
        line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
        failct += 1
      cmd = ""
      cmd += "grep -e"
      cmd += ' "Total wall time: " ' #Find run time
      cmd += "tests.out"
      try:
        RunTime = subprocess.check_output(cmd,shell=True) #Get run time
        RunTime = RunTime.split()
        RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
      except:
         RunTime = " N/A"
      line += RunTime
      if (UpdateResults == 1):
        line += ", "
        line += SavedEnergy
      print(line)

      #Clean up files
      cmd = ""
      cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
      cmd += " BeadStartStruct.xyz BurstStruct.xyz"
      if (QMPack == "Gaussian"):
        #Remove checkpoint files
        cmd += " *.chk"
      if (QMPack == "PSI4"):
        #Remove checkpoint files
        cmd += " timer.* psi.* *.32 *.180"
      if (QMPack == "NWChem"):
        #Remove checkpoint files
        cmd += " *.movecs"
      subprocess.call(cmd,shell=True)

      #Check MM energy
      line = ""
      PassEnergy = 0
      cmd = "cp pol.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      cmd = "lichem -n "
      cmd += `Ncpus`
      cmd += " "
      cmd += "-x waterdimer.xyz "
      cmd += "-r solvreg.inp "
      cmd += "-c watercon.inp "
      cmd += "-o trash.xyz "
      cmd += "> tests.out " #Capture stdout
      cmd += "2>&1" #Capture stderr
      subprocess.call(cmd,shell=True) #Run calculations
      cmd = ""
      cmd += "grep -e"
      cmd += ' "MM energy: " ' #Find final energy
      cmd += "tests.out"
      SavedEnergy = "Crashed..."
      try:
        #Safely check energy
        QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
        QMMMEnergy = QMMMEnergy.split()
        QMMMEnergy = float(QMMMEnergy[2])
        SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
        QMMMEnergy = round(QMMMEnergy,5)
      except:
        #Calculation failed
        QMMMEnergy = 0.0
      #Check against saved energy
      if (QMMMEnergy == round(-1.2549403662026,5)):
        PassEnergy = 1
      line += " AMOEBA/GK energy:    "
      if (PassEnergy == 1):
        line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
        passct += 1
      else:
        line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
        failct += 1
      cmd = ""
      cmd += "grep -e"
      cmd += ' "Total wall time: " ' #Find run time
      cmd += "tests.out"
      try:
        RunTime = subprocess.check_output(cmd,shell=True) #Get run time
        RunTime = RunTime.split()
        RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
      except:
         RunTime = " N/A"
      line += RunTime
      if (UpdateResults == 1):
        line += ", "
        line += SavedEnergy
      print(line)

      #Clean up files
      cmd = ""
      cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
      cmd += " BeadStartStruct.xyz BurstStruct.xyz"
      if (QMPack == "Gaussian"):
        #Remove checkpoint files
        cmd += " *.chk"
      if (QMPack == "PSI4"):
        #Remove checkpoint files
        cmd += " timer.* psi.* *.32 *.180"
      if (QMPack == "NWChem"):
        #Remove checkpoint files
        cmd += " *.movecs"
      subprocess.call(cmd,shell=True)

      #Check QMMM point-charge energy results
      line = ""
      PassEnergy = 0
      cmd = "cp pchrg.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      cmd = "lichem -n "
      cmd += `Ncpus`
      cmd += " "
      cmd += "-x waterdimer.xyz "
      cmd += "-r pchrgreg.inp "
      cmd += "-c watercon.inp "
      cmd += "-o trash.xyz "
      cmd += "> tests.out " #Capture stdout
      cmd += "2>&1" #Capture stderr
      subprocess.call(cmd,shell=True) #Run calculations
      cmd = ""
      cmd += "grep -e"
      cmd += ' "QMMM energy: " ' #Find final energy
      cmd += "tests.out"
      SavedEnergy = "Crashed..."
      try:
        #Safely check energy
        QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
        QMMMEnergy = QMMMEnergy.split()
        QMMMEnergy = float(QMMMEnergy[2])
        SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
        QMMMEnergy = round(QMMMEnergy,5)
      except:
        #Calculation failed
        QMMMEnergy = 0.0
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.2021947277,5)):
          PassEnergy = 1
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.2018207808,5)):
          PassEnergy = 1
      if (QMPack == "NWChem"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.2022117306,5)):
          PassEnergy = 1
      line += " PBE0/TIP3P energy:   "
      if (PassEnergy == 1):
        line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
        passct += 1
      else:
        line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
        failct += 1
      cmd = ""
      cmd += "grep -e"
      cmd += ' "Total wall time: " ' #Find run time
      cmd += "tests.out"
      try:
        RunTime = subprocess.check_output(cmd,shell=True) #Get run time
        RunTime = RunTime.split()
        RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
      except:
        RunTime = " N/A"
      line += RunTime
      if (UpdateResults == 1):
        line += ", "
        line += SavedEnergy
      print(line)

      #Clean up files
      cmd = ""
      cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
      cmd += " BeadStartStruct.xyz BurstStruct.xyz"
      if (QMPack == "Gaussian"):
        #Remove checkpoint files
        cmd += " *.chk"
      if (QMPack == "PSI4"):
        #Remove checkpoint files
        cmd += " timer.* psi.* *.32 *.180"
      if (QMPack == "NWChem"):
        #Remove checkpoint files
        cmd += " *.movecs"
      subprocess.call(cmd,shell=True)

      #Check QMMM polarizable energy results
      line = ""
      PassEnergy = 0
      cmd = "cp pol.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      cmd = "lichem -n "
      cmd += `Ncpus`
      cmd += " "
      cmd += "-x waterdimer.xyz "
      cmd += "-r polreg.inp "
      cmd += "-c watercon.inp "
      cmd += "-o trash.xyz "
      cmd += "> tests.out " #Capture stdout
      cmd += "2>&1" #Capture stderr
      subprocess.call(cmd,shell=True) #Run calculations
      cmd = ""
      cmd += "grep -e"
      cmd += ' "QMMM energy: " ' #Find final energy
      cmd += "tests.out"
      SavedEnergy = "Crashed..."
      try:
        #Safely check energy
        QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
        QMMMEnergy = QMMMEnergy.split()
        QMMMEnergy = float(QMMMEnergy[2])
        SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
        QMMMEnergy = round(QMMMEnergy,5)
      except:
        #Calculation failed
        QMMMEnergy = 0.0
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.1114201829,5)):
          PassEnergy = 1
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.1090319595,5)):
          PassEnergy = 1
      if (QMPack == "NWChem"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.1094168459,5)):
          PassEnergy = 1
      line += " PBE0/AMOEBA energy:  "
      if (PassEnergy == 1):
        line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
        passct += 1
      else:
        line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
        failct += 1
      cmd = ""
      cmd += "grep -e"
      cmd += ' "Total wall time: " ' #Find run time
      cmd += "tests.out"
      try:
        RunTime = subprocess.check_output(cmd,shell=True) #Get run time
        RunTime = RunTime.split()
        RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
      except:
        RunTime = " N/A"
      line += RunTime
      if (UpdateResults == 1):
        line += ", "
        line += SavedEnergy
      print(line)

      #Clean up files
      cmd = ""
      cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
      cmd += " BeadStartStruct.xyz BurstStruct.xyz"
      if (QMPack == "Gaussian"):
        #Remove checkpoint files
        cmd += " *.chk"
      if (QMPack == "PSI4"):
        #Remove checkpoint files
        cmd += " timer.* psi.* *.32 *.180"
      if (QMPack == "NWChem"):
        #Remove checkpoint files
        cmd += " *.movecs"
      subprocess.call(cmd,shell=True)

      #Check pseudobond optimizations
      if ((QMPack == "Gaussian") or (QMPack == "NWChem")):
        line = ""
        PassEnergy = 0
        cmd = "cp pbopt.key tinker.key"
        subprocess.call(cmd,shell=True) #Copy key file
        cmd = "cp pbbasis.txt BASIS"
        subprocess.call(cmd,shell=True) #Copy BASIS set file
        cmd = "lichem -n "
        cmd += `Ncpus`
        cmd += " "
        cmd += "-x alkyl.xyz "
        cmd += "-r pboptreg.inp "
        cmd += "-c alkcon.inp "
        cmd += "-o trash.xyz "
        cmd += "> tests.out " #Capture stdout
        cmd += "2>&1" #Capture stderr
        subprocess.call(cmd,shell=True) #Run calculations
        cmd = ""
        cmd += "grep -e"
        cmd += ' "Opt. step: 2 " ' #Find final energy
        cmd += "tests.out"
        SavedEnergy = "Crashed..."
        try:
          #Safely check energy
          QMMMEnergy = subprocess.check_output(cmd,shell=True) #Get results
          QMMMEnergy = QMMMEnergy.split()
          QMMMEnergy = float(QMMMEnergy[6])
          SavedEnergy = "Energy: "+`QMMMEnergy` #Save it for later
          QMMMEnergy = round(QMMMEnergy,5)
        except:
          #calculation failed
          QMMMEnergy = 0.0
        if (QMPack == "Gaussian"):
          #Check against saved energy
          if (QMMMEnergy == round(-3015.0548490566,5)):
            PassEnergy = 1
        if (QMPack == "NWChem"):
          #Check against saved energy
          if (QMMMEnergy == round(-3015.2278310975,5)):
            PassEnergy = 1
        line += " DFP/Pseudobonds:     "
        if (PassEnergy == 1):
          line += ClrSet.TPass+"Pass"+ClrSet.Reset+","
          passct += 1
        else:
          line += ClrSet.TFail+"Fail"+ClrSet.Reset+","
          failct += 1
        cmd = ""
        cmd += "grep -e"
        cmd += ' "Total wall time: " ' #Find run time
        cmd += "tests.out"
        try:
          RunTime = subprocess.check_output(cmd,shell=True) #Get run time
          RunTime = RunTime.split()
          RunTime = " "+('%.4f'%round(float(RunTime[3]),4))+" "+RunTime[4]
        except:
          RunTime = " N/A"
        line += RunTime
        if (UpdateResults == 1):
          line += ", "
          line += SavedEnergy
        print(line)

        #Clean up files
        cmd = ""
        cmd += "rm -f BASIS tinker.key tests.out trash.xyz"
        cmd += " BeadStartStruct.xyz BurstStruct.xyz"
        if (QMPack == "Gaussian"):
          #Remove checkpoint files
          cmd += " *.chk"
        if (QMPack == "PSI4"):
          #Remove checkpoint files
          cmd += " timer.* psi.* *.32 *.180"
        if (QMPack == "NWChem"):
          #Remove checkpoint files
          cmd += " *.movecs"
        subprocess.call(cmd,shell=True)

    #Print blank line and change directory
    line = ""
    print(line)
    os.chdir("../")

#Start printing the statistics
line = ""
line += "***************************************************"
line += '\n'
line += '\n'
line += "Statistics:"
line += '\n'
line += " Tests passed: "+`passct`+'\n'
line += " Tests failed: "+`failct`+'\n'

#Stop timer
EndTime = time.time()
TotalTime = (EndTime-StartTime)

#Find the correct units
TimeUnits = " seconds"
if (TotalTime > 60):
  TotalTime /= 60.0
  TimeUnits = " minutes"
  if (TotalTime > 60):
    TotalTime /= 60.0
    TimeUnits = " hours"
    if (TotalTime > 24):
      TotalTime /= 24
      TimeUnits = " days"

#Finish printing the statistics
line += " Total run time: "+('%.2f'%round(TotalTime,2))+TimeUnits+'\n'
line += '\n'
line += "***************************************************"
line += '\n'
line += '\n'
line += "Done."
line += '\n'
print(line)

#Quit
exit(0)

