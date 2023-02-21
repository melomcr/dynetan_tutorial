######################################################################
######### Dynamical Network Analysis - Tutorial for command line interface #####
######################################################################

# Load the python package
import dynetan as dna
from dynetan.toolkit import getSelFromNode

import os
import networkx as nx
from itertools import islice
import warnings

# Here we turn off `notebookMode` to activate progress bars in the
# command line.
dnap = dna.proctraj.DNAproc(notebookMode=False)

######################################################################
######### User defined values #######
######################################################################

# Path where input files will searched and results be written.
workDir = "./TutorialData/"

# PSF file name
psfFile = os.path.join(workDir, "decarboxylase.0.psf")

# DCD file name
dcdFiles = [os.path.join(workDir, "decarboxylase.1.dcd")]
# dcdFiles = [os.path.join(workDir, "decarboxylase.1.short.dcd")]

ligandSegID = "OMP"

# Segment IDs for regions that will be studied.
segIDs = ["OMP","ENZY"]

# Residue name for solvent molecule(s)
h2oName = ["TIP3"]

# Number of windows created from full simulation.
numWinds = 4

# Sampled frames per window
numSampledFrames = 10

usrNodeGroups = {}

usrNodeGroups["TIP3"] = {}
usrNodeGroups["TIP3"]["OH2"] = set("OH2 H1 H2".split())

usrNodeGroups["OMP"] = {}
usrNodeGroups["OMP"]["N1"] = set("N1 C2 O2 N3 C4 O4 C5 C6 C7 OA OB".split())
usrNodeGroups["OMP"]["P"] = set("P OP1 OP2 OP3 O5' C5' C4' O4' C1' C3' C2' O2' O3'".split())

#################################
### Extra configuration

# Cutoff for contact map (In Angstroms)
cutoffDist = 4.5

# Minimum contact persistance (In ratio of total trajectory frames)
contactPersistence = 0.75

#################################
### Load info to object

dnap.setNumWinds(numWinds)
dnap.setNumSampledFrames(numSampledFrames)
dnap.setCutoffDist(cutoffDist)
dnap.setContactPersistence(contactPersistence)
dnap.setSolvNames(h2oName)
dnap.setSegIDs(segIDs)

dnap.setNodeGroups(usrNodeGroups)

######################################################################
######### Load topology and trajectory files #######
######################################################################

print("Loading topology file {} and trajectory file(s) {}.".format(psfFile, dcdFiles))

dnap.loadSystem(psfFile, dcdFiles)

print("System loaded.")

# We can access the trajectory data directly.
print("MDAnalysis universe:", dnap.getU().trajectory)

######################################################################
######### Prepare system for network calculations  #######
######################################################################

dnap.checkSystem()

dnap.selectSystem(withSolvent=True)

dnap.prepareNetwork()

print("Aligning trajectory...")
dnap.alignTraj()

######################################################################
######### Determine contacts and calculate correlations  #######
######################################################################

print("Finding contacts...")
dnap.findContacts(stride=1, verbose=True)

# This may be necessary for systems with low default recursion limits.
import sys
print("Recursion limit:", sys.getrecursionlimit())
sys.setrecursionlimit(3000)
print("New recursion limit:", sys.getrecursionlimit())

print("Filtering contacts...")
dnap.filterContacts(notSameRes=True, notConsecutiveRes=False, removeIsolatedNodes=True)

dnap.calcCor(ncores=4)

######################################################################
######### Determine cartesian distances  #######
######################################################################

dnap.calcCartesian(backend="serial")

######################################################################
######### Determine Network Properties #######
######################################################################

dnap.calcGraphInfo()

# Basic information of the network as interpreted as a graph.
print("Graph with {} nodes and {} edges".format(len(dnap.nxGraphs[0].nodes),
                                                len(dnap.nxGraphs[0].edges)))

# Both density and transitivity are scaled from 0 to 1
for win in range(dnap.numWinds):
    print("----- Window {} -----".format(win))
    print("Density:", round( nx.density(dnap.nxGraphs[win]), 4) )
    print("Transitivity:", round( nx.transitivity(dnap.nxGraphs[win]), 4) )
    print()

from operator import itemgetter

# We can check the nodes that have the most connections in each window.
for win in range(dnap.numWinds):
    print("----- Window {} -----".format(win))

    sorted_degree = sorted(dnap.getDegreeDict(win).items(), key=itemgetter(1), reverse=True)

    print("Top 5 nodes by degree: [node --> degree : selection]")
    for n,d in sorted_degree[:5]:
        print("{0:>4} --> {1:>2} : {2}".format(n, d, getSelFromNode(n, dnap.nodesAtmSel)))

    print()

# calculate optimal paths
print("Calculating optimal paths...")
dnap.calcOptPaths(ncores=4)

print("Calculating edge betweenness...")
# calculate betweenness values
dnap.calcBetween(ncores=4)

print("Here are the top 5 pairs of nodes based on Betweenness values, "
      "compared to their correlation values (in Window 0):")
for k, v in islice(dnap.btws[0].items(),5):
    node_pair = k
    btw = round(v, 3)
    corr = round(dnap.corrMatAll[0, k[0], k[1]], 3)
    print(f"\tNodes {node_pair} have betweenness {btw} and correlation {corr}.")

dnap.calcEigenCentral()

dnap.calcCommunities()

print("Here are the top 5 communities based on number of nodes:")
# Sort communities based on number of nodes
for comIndx in islice(dnap.nodesComm[0]["commOrderSize"], 5):
    print("Modularity Class {0:>2}: {1:>3} nodes.".format(comIndx,
                                                          len(dnap.nodesComm[0]["commNodes"][comIndx])))

print("Here are the top 5 communities based on Eigenvector Centrality:")
# Sort communities based on the node with the highest eigenvector centrality
for comIndx in islice(dnap.nodesComm[0]["commOrderEigenCentr"], 5):
    print("Modularity Class {0} ({1} nodes) Sorted by Eigenvector Centrality:".format(
                                                                    comIndx,
                                                                len(dnap.nodesComm[0]["commNodes"][comIndx])))
    for node in dnap.nodesComm[0]["commNodes"][comIndx][:5]:
        print("Name: {0:>4} | Degree: {1:>2} | Eigenvector Centrality: {2}".format(
            node, dnap.nxGraphs[win].nodes[node]['degree'], dnap.nxGraphs[win].nodes[node]['eigenvector']))
    print()

######################################################################
######### Determine Interface Properties #######
######################################################################

# This method returns the number of unique nodes in interface node pairs.
num_nodes = dnap.interfaceAnalysis(selAstr="segid ENZY", selBstr="segid OMP")
print(f"There are {num_nodes} nodes making contacts along the complex interface.")

######################################################################
######### Save data for analysis and plots #######
######################################################################

pathToData = "./TutorialResults/"

fileNameRoot = "dnaData"

fullPathRoot = os.path.join(pathToData, fileNameRoot)

dnap.saveData(fullPathRoot)

# This function will save a reduced DCD trajectory with the heavy atoms used for network analysis
# A smaller trajectory can be created by choosing a "stride" that sub-samples the original trajectory.
# This function will also produce a PDB file so that information on atoms and residues can be loaded to
#    visualization software such as VMD.

dcd_stride = 1
num_atoms = dnap.workU.atoms.n_atoms
num_frames = len(dnap.workU.trajectory[::dcd_stride])

print(f"We will save {num_atoms} heavy atoms and {num_frames} frames.")

# MDAnalysis may print warnings regarding missing data fields, such as altLocs,
# icodes, occupancies, or tempfactor, which provide information commonly found
# in PDB files. The warnings are for your information, and in the context of
# this tutorial, they are expected and do not indicate a problem. We will silence
# such "UserWarning"s for clarity.
warnings.filterwarnings("ignore", category=UserWarning)

dnap.saveReducedTraj(fullPathRoot, stride=dcd_stride)

######################################################################
######### END OF TUTORIAL #######
######################################################################

print("\n\n END OF TUTORIAL \n\n")
