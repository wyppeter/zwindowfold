# ZWindowFold
# Sliding-window RNA-folding algorithm to determine the z-score of ordered structure in the RNA sequence
# Identifies proposed structural domains
# v2019.12.13
# Peter Wang
# Bartel Lab / Simon Lab
# Requires: RNAfold command line installed

# python ZWindowFold.py [input.fa] [output.csv] [WINDOWSIZE] [WINDOWSHIFT]

# Import modules
import sys
from os import remove as osremove
import subprocess
from random import sample as randomsample
from re import findall as refindall
import numpy

# Get args
# Input FASTA file
try:
    inputFASTA = sys.argv[1]
except IndexError:
    inputFASTA = "input.fa"
# Output CSV file
try:
    outputCSV = sys.argv[2]
except IndexError:
    outputCSV = "output.csv"
# Configure window size (do 50, 150, 500)
try:
    WINDOWSIZE = int(sys.argv[3])
except IndexError:
    WINDOWSIZE = 150
# Configure window interval (10, 20, 100)
try:
    WINDOWSHIFT = int(sys.argv[4])
except IndexError:
    WINDOWSHIFT = 10

# Define config and constants
NUMRAN = 20  # Number of shuffling controls to test
AVGOF = 5  # Number of repeats for z-score mean

def RNAfolddG(rnaseq):
    """ Wrapper to call RNAfold """
    # returns dG
    foldOutput = subprocess.run(
            ["RNAfold", "--noPS"],
            input = rnaseq,
            capture_output = True, text = True
            ).stdout
    return float(foldOutput.split("\n")[1].split()[-1].strip("()"))

# Define functions
def zscore(vals):
    # Calculate z-score from a list of values
    # First item is positive, the rest are shuffled controls
    zval = (vals[0] - numpy.mean(vals[1:])) / numpy.std(vals)
    if numpy.isnan(zval):
        return 0
    else:
        return zval

# Get inputs and seq storage structures
with open(inputFASTA, "r") as inFile:
    seq = "".join([l.strip() for l in inFile.readlines()[1:]]).upper()
SEQLENGTH = len(seq)   # total length of seq
WINDOWPOS = [WINDOWSHIFT * ind for ind in range((SEQLENGTH - WINDOWSIZE) // WINDOWSHIFT + 1)]  # all the starting positions of the windows

print("Analyzing %s..." % inputFASTA)

# Run this n times to get average z-score for each window
outDictSet = []  # saving the outDicts for each trial as each listing in this list
outDict = {}  # holder for each trial
for trial in range(AVGOF):
    print("Run#%s in progress" % str(trial+1))

    # Set up region-sequence_list dictionary
    seqWindows = {pos:[seq[pos:pos+WINDOWSIZE]] for pos in WINDOWPOS}

    # Generate random shuffled strings and attach to the 1-member lists of seq values in the dictionary above
    for segmentList in seqWindows.values():
        for _ in range(NUMRAN):
            segmentList += ["".join(randomsample(segmentList[0], WINDOWSIZE))]

    # Here comes the folding; extract the free energies, and calculate z-scores
    # outDict is a dict, window pos (start) as keys, values being a list of dG floats (first is WT, rest are shuffled)
    outDict = {}
    for pos in WINDOWPOS:
        for windowSeq in seqWindows[pos]:
            dGthis = RNAfolddG(windowSeq)
            outDict[pos] = outDict.get(pos,[]) + [dGthis]

    # Load into big data structure
    outDictSet += [outDict]

# Get minimum free energies of WT seq windows (first item in each list of dGs)
energies = {indx:outDictSet[0][indx][0] for indx in sorted(outDictSet[0].keys())}

print("Printing output...")

# Output
outputCSVFile = open(outputCSV, "w")
# Header
outputCSVFile.write("pos,range,%s,mean,stdev,minE\n" % (",".join([("trial" + str(n+1)) for n in range(AVGOF)])))
# Each row is: position (mid-point), range, trial#_z, trial#_z, ..., mean, stdev, minE
for ind in sorted(energies.keys()):

    windowran = str(ind + 1) + "-" + str(ind + WINDOWSIZE)  # min - max (1-based, inclusive)

    outputCSVFile.write(",".join([
        str(int(ind + WINDOWSIZE/2)),  # midpoint
        windowran,
        (",".join(["%.10f"%(zscore(out[ind])) for out in outDictSet])),
        ("%.10f" % numpy.mean(
            [zscore(out[ind]) for out in outDictSet]
            )),
        ("%.10f" % numpy.std(
            [zscore(out[ind]) for out in outDictSet]
            )),
        ("%.10f" % energies[ind])
        ]) + "\n")
outputCSVFile.close()

print("Done!")
