###################################################
#                                                 #
#   LICHEM: Layered Interacting CHEmical Models   #
#                                                 #
#        Symbiotic Computational Chemistry        #
#                                                 #
###################################################

#LICHEM semi-automated test suite

###
#  Usage:
#
#    user:$ ./runtests Ncpus All
#      or
#    user:$ ./runtests Ncpus QMPackage MMPackage
#     or
#    user:$ ./runtests Ncpus QMPackage MMPackage dry
####

#Modules
import subprocess
import time
import sys
import os

#Start timer immediately
startTime = time.time()

#Initialize globals
TTxtLen = 30 #Number of characters for the test name
passCt = 0 #Number of tests passed
failCt = 0 #Number of tests failed

#Development settings
#NB: Modified by the Makefile
updateResults = 0 #Flag to print energies to update tests
forceAll = 0 #Flag to force it to do tests even if they will fail

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

#Functions
def RunLICHEM(xName,rName,cName):
  #Submit LICHEM jobs
  cmd = "lichem -n "
  cmd += str(Ncpus)
  cmd += " "
  cmd += "-x "
  cmd += xName
  cmd += " "
  cmd += "-r "
  cmd += rName
  cmd += " "
  cmd += "-c "
  cmd += cName
  cmd += " "
  cmd += "-o trash.xyz "
  cmd += "> tests.out " #Capture stdout
  cmd += "2>&1" #Capture stderr
  subprocess.call(cmd,shell=True) #Run calculations
  return

def CleanFiles():
  #Delete junk files
  cleanCmd = "rm -f"
  #Remove LICHEM files
  cleanCmd += " BASIS tests.out trash.xyz"
  cleanCmd += " BeadStartStruct.xyz BurstStruct.xyz"
  #Remove TINKER files
  cleanCmd += " tinker.key"
  #Remove LAMMPS files
  cleanCmd += ""
  #Remove AMBER files
  cleanCmd += ""
  #Remove Gaussian files
  cleanCmd += " *.chk"
  #Remove PSI4 files
  cleanCmd += " timer.* psi.* *.32 *.180"
  #Remove NWChem files
  cleanCmd += " *.movecs"
  #Delete the files
  subprocess.call(cleanCmd,shell=True)
  return

def RecoverEnergy(txtLabel,itemNum):
  #Recover the energy from LICHEM output
  cmd = ""
  cmd += "grep -e "
  cmd += '"'
  cmd += txtLabel
  cmd += ' "'
  cmd += " tests.out"
  savedResult = "Crashed..."
  try:
    #Safely check energy
    finalEnergy = subprocess.check_output(cmd,shell=True) #Get results
    finalEnergy = finalEnergy.decode('utf-8').split()
    finalEnergy = float(finalEnergy[itemNum])
    savedResult = "Energy: "+str(finalEnergy) #Save it for later
    finalEnergy = round(finalEnergy,5)
  except:
    #Calculation failed
    finalEnergy = 0.0
  return finalEnergy,savedResult

def RecoverFreqs():
  #Recover a list of frequencies
  cmd = ""
  cmd += "sed '/Usage Statistics/,$d' tests.out | "
  cmd += "sed -n '/Frequencies:/,$p' | "
  cmd += "sed '/Frequencies:/d'"
  try:
    #Safely check energy
    freqList = []
    tmpFreqs = subprocess.check_output(cmd,shell=True) #Get results
    tmpFreqs = tmpFreqs.decode('utf-8').strip()
    tmpFreqs = tmpFreqs.split()
    for freqVal in tmpFreqs:
      freqList.append(float(freqVal))
  except:
    #Calculation failed
    freqList = []
  return freqList

def AddPass(tName,testPass,txtLn):
  #Add a colored pass or fail message
  global passCt
  global failCt
  global TTxtLen
  #Add the name of the test
  tName = " "+tName
  deltaTxt = TTxtLen-len(tName)
  if (deltaTxt > 0):
    #Make the test name consistent with TTxtLen
    for i in range(deltaTxt):
      tName += " "
  else:
    #Update bad length
    TTxtLen -= deltaTxt
    TTxtLen += 1
    deltaTxt = TTxtLen-len(tname)
    for i in range(deltaTxt):
      tName += " "
  #Label as pass or fail
  txtLn += tName
  if (testPass == 1):
    txtLn += ClrSet.TPass+"Pass"+ClrSet.Reset+","
    passCt += 1
  else:
    txtLn += ClrSet.TFail+"Fail"+ClrSet.Reset+","
    failCt += 1
  return txtLn

def AddRunTime(txtLn):
  #Collect the LICHEM run time and add it to a string
  cmd = ""
  cmd += "grep -e"
  cmd += ' "Total wall time: " ' #Find run time
  cmd += "tests.out"
  try:
    runTime = subprocess.check_output(cmd,shell=True) #Get run time
    runTime = runTime.decode('utf-8').split()
    runTime = " "+('%.4f'%round(float(runTime[3]),4))+" "+runTime[4]
  except:
    runTime = " N/A"
  txtLn += runTime
  return txtLn

def AddEnergy(devOpt,txtLn,enVal):
  if (devOpt == 1):
    txtLn += ", "
    txtLn += enVal
  return txtLn

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
dryRun = 0 #Only check packages
allTests = 0 #Run all tests at once
if (len(sys.argv) == 3):
  if ((sys.argv[2]).lower() == "all"):
    #Automatically run all tests
    allTests = 1
if (len(sys.argv) < 4):
  line = ""
  if (allTests == 0):
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
    LICHEMbin = ClrSet.TPass+LICHEMbin.decode('utf-8').strip()+ClrSet.Reset
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
    QMbin = ClrSet.TPass+QMbin.decode('utf-8').strip()+ClrSet.Reset
  except:
    QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += QMbin
  line += '\n'
  line += " Gaussian: "
  cmd = "which g09"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = ClrSet.TPass+QMbin.decode('utf-8').strip()+ClrSet.Reset
  except:
    QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += QMbin
  line += '\n'
  line += " NWChem: "
  cmd = "which nwchem"
  try:
    #Find path
    QMbin = subprocess.check_output(cmd,shell=True)
    QMbin = ClrSet.TPass+QMbin.decode('utf-8').strip()+ClrSet.Reset
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
    MMbin = ClrSet.TPass+MMbin.decode('utf-8').strip()+ClrSet.Reset
  except:
    MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += MMbin
  line += '\n'
  line += " LAMMPS: "
  cmd = "which lammps"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = ClrSet.TPass+MMbin.decode('utf-8').strip()+ClrSet.Reset
  except:
    MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += MMbin
  line += '\n'
  line += " AMBER: "
  cmd = "which pmemd"
  try:
    #Find path
    MMbin = subprocess.check_output(cmd,shell=True)
    MMbin = ClrSet.TPass+MMbin.decode('utf-8').strip()+ClrSet.Reset
  except:
    MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
  line += MMbin
  line += '\n'
  print(line)
  if (allTests == 0):
    #Quit
    exit(0)

Ncpus = int(sys.argv[1]) #Threads
if (allTests == 0):
  QMPack = sys.argv[2] #QM wrapper for calculations
  QMPack = QMPack.lower()
  MMPack = sys.argv[3] #MM wrapper for calculations
  MMPack = MMPack.lower()
  if (len(sys.argv) > 4):
    if ((sys.argv[4]).lower() == "dry"):
      #Quit early
      dryRun = 1

#Initialize variables
LICHEMbin = ""
QMbin = ""
MMbin = ""

#Check packages and identify missing binaries
badLICHEM = 0
cmd = "which lichem"
try:
  #Find path
  LICHEMbin = subprocess.check_output(cmd,shell=True)
  LICHEMbin = LICHEMbin.decode('utf-8').strip()
except:
  LICHEMbin = "N/A"
  badLICHEM = 1
if (badLICHEM == 1):
  #Quit with an error
  line = ""
  line += "Error: LICHEM binary not found!"
  line += '\n'
  print(line)
  exit(0)
if (allTests == 0):
  badQM = 1
  if ((QMPack == "psi4") or (QMPack == "psi")):
    QMPack = "PSI4"
    cmd = "which psi4"
    try:
      #Find path
      QMbin = subprocess.check_output(cmd,shell=True)
      QMbin = QMbin.decode('utf-8').strip()
    except:
      QMbin = "N/A"
    badQM = 0
  if ((QMPack == "gaussian") or (QMPack == "g09")):
    QMPack = "Gaussian"
    cmd = "which g09"
    try:
      #Find path
      QMbin = subprocess.check_output(cmd,shell=True)
      QMbin = QMbin.decode('utf-8').strip()
    except:
      QMbin = "N/A"
    badQM = 0
  if (QMPack == "nwchem"):
    QMPack = "NWChem"
    cmd = "which nwchem"
    try:
      #Find path
      QMbin = subprocess.check_output(cmd,shell=True)
      QMbin = QMbin.decode('utf-8').strip()
    except:
      QMbin = "N/A"
    badQM = 0
  if (badQM == 1):
    #Quit with an error
    line = '\n'
    line += "Error: QM package name '"
    line += QMPack
    line += "' not recognized."
    line += '\n'
    print(line)
    exit(0)
  badMM = 1
  if (MMPack == "tinker"):
    MMPack = "TINKER"
    cmd = "which analyze"
    try:
      #Find path
      MMbin = subprocess.check_output(cmd,shell=True)
      MMbin = MMbin.decode('utf-8').strip()
    except:
      MMbin = "N/A"
    badMM = 0
  if (MMPack == "lammps"):
    MMPack = "LAMMPS"
    cmd = "which lammps"
    try:
      #Find path
      MMbin = subprocess.check_output(cmd,shell=True)
      MMbin = MMbin.decode('utf-8').strip()
    except:
      MMbin = "N/A"
    badMM = 0
  if (MMPack == "amber"):
    MMPack = "AMBER"
    cmd = "which pmemd" #Check
    try:
      #Find path
      MMbin = subprocess.check_output(cmd,shell=True)
      MMbin = MMbin.decode('utf-8').strip()
    except:
      MMbin = "N/A"
    badMM = 0
  if (badMM == 1):
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
line += str(Ncpus)
line += '\n'
if (allTests == 0):
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
  if (forceAll == 1):
    line += " Mode: Development"
  else:
    line += " Mode: All tests"
  line += '\n'
if (dryRun == 1):
  line += " Mode: Dry run"
  line += '\n'
print(line)

#Escape for dry runs
if (dryRun == 1):
  #Quit without an error
  line = "Dry run completed."
  line += '\n'
  print(line)
  exit(0)

#Escape if binaries not found
if (((QMbin == "N/A") or (MMbin == "N/A")) and (allTests == 0)):
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
if (allTests == 1):
  #Safely add PSI4
  cmd = "which psi4"
  try:
    #Run PSI4 tests
    packBin = subprocess.check_output(cmd,shell=True)
    QMTests.append("PSI4")
  except:
    #Skip tests that will fail
    if (forceAll == 1):
      QMTests.append("PSI4")
  #Safely add Gaussian
  cmd = "which g09"
  try:
    #Run Gaussian tests
    packBin = subprocess.check_output(cmd,shell=True)
    QMTests.append("Gaussian")
  except:
    #Skip tests that will fail
    if (forceAll == 1):
      QMTests.append("Gaussian")
  #Safely add NWChem
  cmd = "which nwchem"
  try:
    #Run NWChem tests
    packBin = subprocess.check_output(cmd,shell=True)
    QMTests.append("NWChem")
  except:
    #Skip tests that will fail
    if (forceAll == 1):
      QMTests.append("NWChem")
  #Safely add TINKER
  cmd = "which analyze"
  try:
    #Run TINKER tests
    packBin = subprocess.check_output(cmd,shell=True)
    MMTests.append("TINKER")
  except:
    #Skip tests that will fail
    if (forceAll == 1):
      MMTests.append("TINKER")
  #Safely add lammps
  cmd = "which lammps"
  try:
    #Run LAMMPS tests
    packBin = subprocess.check_output(cmd,shell=True)
    MMTests.append("LAMMPS")
  except:
    #Skip tests that will fail
    if (forceAll == 1):
      MMTests.append("LAMMPS")
  #Safely add AMBER
  cmd = "which pmemd"
  try:
    #Run AMBER tests
    packBin = subprocess.check_output(cmd,shell=True)
    MMTests.append("AMBER")
  except:
    #Skip tests that will fail
    if (forceAll == 1):
      MMTests.append("AMBER")
else:
  #Add only the specified packages
  QMTests.append(QMPack)
  MMTests.append(MMPack)

#NB: Tests are in the following order:
#     1) HF energy
#     2) PBE0 energy
#     3) CCSD energy
#     4) PM6 energy
#     5) Frequencies
#     6) NEB TS energy
#     7) TIP3P energy
#     8) AMOEBA/GK energy
#     9) PBE0/TIP3P energy
#    10) PBE0/AMOEBA energy
#    11) DFP/Pseudobonds

#Loop over tests
for qmTest in QMTests:
  for mmTest in MMTests:
    #Set packages
    QMPack = qmTest
    MMPack = mmTest

    #Set path based on packages
    dirPath = ""
    if (QMPack == "PSI4"):
      dirPath += "PSI4_"
    if (QMPack == "Gaussian"):
      dirPath += "Gau_"
    if (QMPack == "NWChem"):
      dirPath += "NWChem_"
    dirPath += MMPack
    dirPath += "/"

    #Change directory
    os.chdir(dirPath)

    #Start printing results
    line = QMPack+"/"+MMPack
    line += " results:"
    print(line)

    #Check HF energy
    if ((QMPack == "PSI4") or (QMPack == "Gaussian")):
      line = ""
      passEnergy = 0
      RunLICHEM("waterdimer.xyz","hfreg.inp","watercon.inp")
      QMMMEnergy,savedEnergy = RecoverEnergy("QM energy:",2)
      #Check result
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-4136.9303981392,5)):
          passEnergy = 1
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-4136.9317704519,5)):
          passEnergy = 1
      line = AddPass("HF energy:",passEnergy,line)
      line = AddRunTime(line)
      line = AddEnergy(updateResults,line,savedEnergy)
      print(line)
      CleanFiles() #Clean up files

    #Check DFT energy
    line = ""
    passEnergy = 0
    RunLICHEM("waterdimer.xyz","pbereg.inp","watercon.inp")
    QMMMEnergy,savedEnergy = RecoverEnergy("QM energy:",2)
    #Check result
    if (QMPack == "PSI4"):
      #Check against saved energy
      if (QMMMEnergy == round(-4154.1683659877,5)):
        passEnergy = 1
    if (QMPack == "Gaussian"):
      #Check against saved energy
      if (QMMMEnergy == round(-4154.1676114324,5)):
        passEnergy = 1
    if (QMPack == "NWChem"):
      #Check against saved energy
      if (QMMMEnergy == round(-4154.1683939169,5)):
        passEnergy = 1
    line = AddPass("PBE0 energy:",passEnergy,line)
    line = AddRunTime(line)
    line = AddEnergy(updateResults,line,savedEnergy)
    print(line)
    CleanFiles() #Clean up files

    #Check CCSD energy
    if (QMPack == "PSI4"):
      line = ""
      passEnergy = 0
      RunLICHEM("waterdimer.xyz","ccsdreg.inp","watercon.inp")
      QMMMEnergy,savedEnergy = RecoverEnergy("QM energy:",2)
      #Check result
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-4147.730483706,5)):
          passEnergy = 1
      line = AddPass("CCSD energy:",passEnergy,line)
      line = AddRunTime(line)
      line = AddEnergy(updateResults,line,savedEnergy)
      print(line)
      CleanFiles() #Clean up files

    #Check PM6 energy
    if (QMPack == "Gaussian"):
      line = ""
      passEnergy = 0
      RunLICHEM("waterdimer.xyz","pm6reg.inp","watercon.inp")
      QMMMEnergy,savedEnergy = RecoverEnergy("QM energy:",2)
      #Check result
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-4.8623027634995,5)):
          passEnergy = 1
      line = AddPass("PM6 energy:",passEnergy,line)
      line = AddRunTime(line)
      line = AddEnergy(updateResults,line,savedEnergy)
      print(line)
      CleanFiles() #Clean up files

    #Check imaginary frequencies
    line = ""
    passEnergy = 0
    RunLICHEM("methfluor.xyz","freqreg.inp","methflcon.inp")
    QMMMFreqs = RecoverFreqs()
    QMMMEnergy = 5e100 #Huge number
    #Sort frequencies
    for freqVal in QMMMFreqs:
      #Find lowest frequency
      if (freqVal < QMMMEnergy):
        QMMMEnergy = freqVal
        savedEnergy = "Freq:   "+str(freqVal)
    #Check for errors
    if (QMMMEnergy > 1e100):
      savedEnergy = "Crashed..."
    #Check results
    if (QMPack == "PSI4"):
      #Check against saved frequency
      if (round(QMMMEnergy,0) == round(-496.58664,0)):
        passEnergy = 1
    if (QMPack == "Gaussian"):
      #Check against saved frequency
      if (round(QMMMEnergy,0) == round(-496.73073,0)):
        passEnergy = 1
    if (QMPack == "NWChem"):
      #Check against saved frequency
      if (round(QMMMEnergy,0) == round(-496.79703,0)):
        passEnergy = 1
    line = AddPass("Frequencies:",passEnergy,line)
    line = AddRunTime(line)
    line = AddEnergy(updateResults,line,savedEnergy)
    print(line)
    CleanFiles() #Clean up files

    #Check NEB optimization
    line = ""
    PassEnergy = 0
    cmd = "cp methflbeads.xyz BeadStartStruct.xyz"
    subprocess.call(cmd,shell=True) #Copy restart file
    RunLICHEM("methfluor.xyz","nebreg.inp","methflcon.inp")
    QMMMEnergy,savedEnergy = RecoverEnergy("Opt. step: 2",11)
    #Check result
    if (QMPack == "PSI4"):
      #Check against saved energy
      if (QMMMEnergy == round(-6511.0580192214,5)):
        passEnergy = 1
    if (QMPack == "Gaussian"):
      #Check against saved energy
      if (QMMMEnergy == round(-6511.0567955964,5)):
        passEnergy = 1
    if (QMPack == "NWChem"):
      #Check against saved energy
      if (QMMMEnergy == round(-6511.0579547077,5)):
        passEnergy = 1
    line = AddPass("NEB TS energy:",passEnergy,line)
    line = AddRunTime(line)
    line = AddEnergy(updateResults,line,savedEnergy)
    print(line)
    CleanFiles() #Clean up files

    #TINKER tests
    if (MMPack == "TINKER"):
      #Check MM energy
      line = ""
      passEnergy = 0
      cmd = "cp pchrg.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      RunLICHEM("waterdimer.xyz","mmreg.inp","watercon.inp")
      QMMMEnergy,savedEnergy = RecoverEnergy("MM energy:",2)
      #Check result
      if (QMMMEnergy == round(-0.2596903536223,5)):
        #Check against saved energy
        passEnergy = 1
      line = AddPass("TIP3P energy:",passEnergy,line)
      line = AddRunTime(line)
      line = AddEnergy(updateResults,line,savedEnergy)
      print(line)
      CleanFiles() #Clean up files

      #Check MM energy
      line = ""
      passEnergy = 0
      cmd = "cp pol.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      RunLICHEM("waterdimer.xyz","solvreg.inp","watercon.inp")
      QMMMEnergy,savedEnergy = RecoverEnergy("MM energy:",2)
      #Check result
      if (QMMMEnergy == round(-1.2549403662026,5)):
        #Check against saved energy
        passEnergy = 1
      line = AddPass("AMOEBA/GK energy:",passEnergy,line)
      line = AddRunTime(line)
      line = AddEnergy(updateResults,line,savedEnergy)
      print(line)
      CleanFiles() #Clean up files

      #Check QMMM point-charge energy results
      line = ""
      passEnergy = 0
      cmd = "cp pchrg.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      RunLICHEM("waterdimer.xyz","pchrgreg.inp","watercon.inp")
      QMMMEnergy,savedEnergy = RecoverEnergy("QMMM energy:",2)
      #Check result
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.2021947277,5)):
          passEnergy = 1
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.2018207808,5)):
          passEnergy = 1
      if (QMPack == "NWChem"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.2022117306,5)):
          passEnergy = 1
      line = AddPass("PBE0/TIP3P energy:",passEnergy,line)
      line = AddRunTime(line)
      line = AddEnergy(updateResults,line,savedEnergy)
      print(line)
      CleanFiles() #Clean up files

      #Check QMMM polarizable energy results
      line = ""
      passEnergy = 0
      cmd = "cp pol.key tinker.key"
      subprocess.call(cmd,shell=True) #Copy key file
      RunLICHEM("waterdimer.xyz","polreg.inp","watercon.inp")
      QMMMEnergy,savedEnergy = RecoverEnergy("QMMM energy:",2)
      #Check result
      if (QMPack == "PSI4"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.1114201829,5)):
          passEnergy = 1
      if (QMPack == "Gaussian"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.1090319595,5)):
          passEnergy = 1
      if (QMPack == "NWChem"):
        #Check against saved energy
        if (QMMMEnergy == round(-2077.1094168459,5)):
          passEnergy = 1
      line = AddPass("PBE0/AMOEBA energy:",passEnergy,line)
      line = AddRunTime(line)
      line = AddEnergy(updateResults,line,savedEnergy)
      print(line)
      CleanFiles() #Clean up files

      #Check pseudobond optimizations
      if ((QMPack == "Gaussian") or (QMPack == "NWChem")):
        line = ""
        passEnergy = 0
        cmd = "cp pbopt.key tinker.key"
        subprocess.call(cmd,shell=True) #Copy key file
        cmd = "cp pbbasis.txt BASIS"
        subprocess.call(cmd,shell=True) #Copy BASIS set file
        RunLICHEM("alkyl.xyz","pboptreg.inp","alkcon.inp")
        QMMMEnergy,savedEnergy = RecoverEnergy("Opt. step: 2",6)
        #Check result
        if (QMPack == "Gaussian"):
          #Check against saved energy
          if (QMMMEnergy == round(-3015.0548490566,5)):
            passEnergy = 1
        if (QMPack == "NWChem"):
          #Check against saved energy
          if (QMMMEnergy == round(-3015.2278310975,5)):
            passEnergy = 1
        line = AddPass("DFP/Pseudobonds:",passEnergy,line)
        line = AddRunTime(line)
        line = AddEnergy(updateResults,line,savedEnergy)
        print(line)
        CleanFiles() #Clean up files

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
line += " Tests passed: "+str(passCt)+'\n'
line += " Tests failed: "+str(failCt)+'\n'

#Stop timer
endTime = time.time()
totalTime = (endTime-startTime)

#Find the correct units
timeUnits = " seconds"
if (totalTime > 60):
  totalTime /= 60.0
  timeUnits = " minutes"
  if (totalTime > 60):
    totalTime /= 60.0
    timeUnits = " hours"
    if (totalTime > 24):
      totalTime /= 24.0
      timeUnits = " days"

#Finish printing the statistics
line += " Total run time: "+('%.2f'%round(totalTime,2))+timeUnits+'\n'
line += '\n'
line += "***************************************************"
line += '\n'
line += '\n'
line += "Done."
line += '\n'
print(line)

#Quit
exit(0)

