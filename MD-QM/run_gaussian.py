#!/hpc2n/eb/software/Python/2.7.18-GCCcore-9.3.0/bin/python

# by Maximilian Scheurer (mscheurer@ks.uiuc.edu), Dec 2016

##################  Imports

from sys import argv as sargv
from sys import exit
from os.path import dirname
import subprocess as sp
import re
import numpy as np


################## Command file parameters

gaussianConfigLines1 = """\
%CPU=0-3
%chk=geom_optim.chk
%mem=16GB
# B3LYP/6-31G* Charge Force Integral=(FineGrid,Acc2E=10) Constants=2006

test structure 

1 1
"""

# if you specified an excited state ( TD(NStates=2,Root=1) ) in the com file and want to give its gradient to NAMD,
# specify the number of the excited state here,
# so that the python script can read the corresponding excited state energy
# if no TD calculation is requested, set excitedState to 0
excitedState=0


################## Processing of file names
gaussianWhitespace = "\n\n"

inputFilename = sargv[1]

#print(inputFilename)

directory = dirname(inputFilename)

# Prepares the name of the configuration file based on full path to data file provided by NAMD
gaussianInFileName = directory + "/"
gaussianInFileName += "qmmm.com"

# Name of the Gaussian log file
gaussianOutFileName = gaussianInFileName +".log"

# Prepares the file name for the file which will be read by NAMD
finalResFileName = inputFilename
finalResFileName += ".result"

################## Reading and parsing NAMD's data ; Writing gaussian's Input

# Reads NAMD data
infile = open(inputFilename,"r")

line = infile.readline()

# Gets number of atoms in the Quantum Chemistry region (= QM atoms + Link atoms)
numQMatms = int(line.split()[0])
# Gets number of point charges
numPntChr = int(line.split()[1].replace("\n",""))

# print("numQMatms:",numQMatms,"; numPntChr",numPntChr)

# stores all lines written to gaussian's input file
outLinesQM = []

# Identation
ident = "  "

lineIndx = 1
pointChargeList = []
pointChargeDict = {}
charges = []
for line in infile:

    posx = line.split()[0]
    posy = line.split()[1]
    posz = line.split()[2]

    if lineIndx <= numQMatms:

        # gaussian's format requires the fileds to be ordered begining with the
        # atom's element symbol, and followed by the XYZ coordinates.

        element = line.split()[3].replace("\n","")

        outLinesQM.append(" ".join([element,posx,posy,posz]) + "\n")

    else:
        #output linebreak to separate atoms from charges
        if lineIndx == numQMatms+1:
            outLinesQM.append("\n")

        # gaussian's format requires the fields to be ordered begining with the
        # XYZ, and followed by the charge .
        pos = " ".join(line.split()[0:3])
        charge = line.split()[3]
        charges.append(charge)
        if pos in pointChargeDict:
        	print "double occurence: ", pos, pointChargeDict[pos] , " with new charge ", charge
        	pointChargeDict[pos] += float(charge)
        	print "new merged charge: ", pointChargeDict[pos]
        else:
        	pointChargeDict[pos] = float(charge)


    lineIndx += 1

cnp = np.array(charges,dtype=float)
print "Sum of the charges: ", np.sum(cnp)

for k in pointChargeDict:
	c = pointChargeDict[k]
	p = k.split()
	pointChargeList.append([p[0],p[1],p[2],str(c)])
	outLinesQM.append(" ".join([p[0],p[1],p[2],'{0:.16f}'.format(c)]) + "\n")

# print len(pointChargeList)
infile.close()

###

with open(gaussianInFileName,"w") as outQMFile:

    outQMFile.write(gaussianConfigLines1)

    for line in outLinesQM:
        outQMFile.write(line)

    outQMFile.write(gaussianWhitespace)

################## Run gaussian

# We first move the shell to the target directory where calculations are to be
# performed
cmdline = "cd " + directory + "; "
# Then we run gaussian with our output file receiving all standard output.

# we probably need to set some environment variables:
import subprocess, os
current_env = os.environ.copy()
current_env["GAUSS_EXEDIR"] = "/hpc2n/eb/software/gaussian/16.C.01-AVX2/g16/"
# subprocess.Popen(my_command, env=my_env)

cmdline += "/hpc2n/eb/software/gaussian/16.C.01-AVX2/g16/g16 "
cmdline += gaussianInFileName + " " + gaussianOutFileName

# print "command:", cmdline
proc = sp.Popen(args=cmdline, shell=True, env=current_env)
proc.wait()

########### READING Gaussian output
excitedStateEnergy=0

tmpOutFile = open(gaussianOutFileName,"r")

qmCharges = []
gradient = []

# Bohr radius for conversion
a0 = 0.52917721067
conversionHB = 627.509469/a0
selfEnergyGaussian = 0.0

# Iterates ultil we find the section of output that contains atomic partial
# charges for QM atoms
chargeSection = False

gradientSection = False

scfenergyFound = False

excitedFound = False

iterate = True

selfEnergyFound = False

while iterate:

    line = tmpOutFile.readline()

    if line.find("Mulliken charges:") != -1:
        chargeSection = True
        # Skips a dividing line
        line = tmpOutFile.readline()
        continue

    if not scfenergyFound:
        if line.find("SCF Done") != -1:
            reg = re.search("SCF Done:(.+)=(.*)([-+]?\d*\.\d+|\d+)(.+)A.U.",line)
            scfenergy = 627.509469 * float(reg.group(2).strip())
            print "SCF energy: ", scfenergy
            scfenergyFound = True

    if not excitedFound and excitedState != 0:
        if line.find("Excited State") != -1:
            line = line.strip().replace(":","").split()
            if int(line[2]) != excitedState:
                continue
            line =tmpOutFile.readline()
            # check if we really requested the excited state that the gradient will be calculated for!
            while line.find("This state for optimization and/or second-order correction.") == -1:
                line = tmpOutFile.readline()
            line = tmpOutFile.readline()
            line = line.strip()
            reg = re.search("([-+]?\d*\.\d+|\d+)(.+)",line)
            excitedStateEnergy = 627.509469 * float(reg.group(1))
            print "Excited State energy: " , excitedStateEnergy
            excitedFound = True

    if line.find("Center     Atomic                   Forces (Hartrees/Bohr)") != -1:
        gradientSection = True
        #skip 2 lines
        line = tmpOutFile.readline()
        line = tmpOutFile.readline()
        line = tmpOutFile.readline()

    if gradientSection:
        line = line.strip().split()
        # don't switch the sign of the number, as Gaussian prints F = -dE/dx, so we can directly pass the numbers to NAMD
        gradient.append( [ str( float(line[2])*conversionHB ), str( float(line[3])*conversionHB ), str( float(line[4])*conversionHB ) ])
        if len(gradient) == numQMatms:
            gradientSection = False
            iterate = False

    if chargeSection:
        qmCharges.append(line.split()[2].replace("\n","").strip())

    if len(qmCharges) == numQMatms:
        chargeSection = False

    if selfEnergyFound != True and line.find("Self energy of the charges") != -1:
    	line = line.strip().split()
    	selfEnergyGaussian = float(line[-2])
        selfEnergyFound = True

tmpOutFile.close()

finFile = open(finalResFileName,"w")

if excitedState == 0:
	finFile.write(str( scfenergy - selfEnergyGaussian*627.509469 ) + "\n")
else:
	finFile.write(str( excitedStateEnergy - selfEnergyGaussian*627.509469 ) + "\n")

for i in range(numQMatms):

    finFile.write(" ".join(gradient[i]) + " " + qmCharges[i] + "\n")

finFile.close()


##########

exit(0)
