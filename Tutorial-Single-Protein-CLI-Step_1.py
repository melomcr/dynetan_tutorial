######################################################################
######### Dynamical Network Analysis - Tutorial for command line interface #####
######################################################################

# Load the python package
import dynetan as dna
from dynetan.toolkit import getSelFromNode, getNodeFromSel

import os, sys
import networkx as nx

from operator import itemgetter
from itertools import islice

# Here we turn off `notebookMode` to activate progress bars in the
# command line.
dnap = dna.proctraj.DNAproc(notebookMode=False)

######################################################################
######### User defined values #######
######################################################################

#################################
### File names and paths

# Path where input files will searched and results be written.
workDir = "./TutorialData/"

# PSF file name
psfFile = os.path.join(workDir, "decarboxylase.0.psf")

# DCD file name
dcdFiles = [os.path.join(workDir, "decarboxylase.1.dcd")]
# dcdFiles = [os.path.join(workDir, "decarboxylase.1.short.dcd")]

# Path where output files will saved.
pathToData = "./TutorialResults/"

fileNameRoot = "dnaData_SinglePtn"

#################################
### System-specific values

# Segment IDs for regions that will be studied.
segIDs = ["ENZY"]

# Residue name for solvent molecule(s)
h2oName = ["TIP3"]

# Network Analysis will auomatically create one node per protein residue (in the alpha carbon).
# This selection can be modified for one or all amino acids by manually defining node groups.

# For non-proteic and non-canonical residues, the user must specify atom(s) that will represent a node.
# We also need to know the heavy atoms that compose each node group.

usrNodeGroups = {}

usrNodeGroups["TIP3"] = {}
usrNodeGroups["TIP3"]["OH2"] = set("OH2 H1 H2".split())

#usrNodeGroups["OMP"] = {}
#usrNodeGroups["OMP"]["N1"] = set("N1 C2 O2 N3 C4 O4 C5 C6 C7 OA OB".split())
#usrNodeGroups["OMP"]["P"] = set("P OP1 OP2 OP3 O5' C5' C4' O4' C1' C3' C2' O2' O3'".split())

#################################
### Analysis values

# Number of windows created from full simulation.
numWinds = 4

# Sampled frames per window
numSampledFrames = 10

# Number of sampled frames for automatic selection of solvent and ions.
# numAutoFrames = numSampledFrames*numWinds

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
dnap.seth2oName(h2oName)
dnap.setSegIDs(segIDs)

dnap.setUsrNodeGroups(usrNodeGroups)

fullPathRoot = os.path.join(pathToData, fileNameRoot)

######################################################################
######### Load topology and trajectory files #######
######################################################################

print("Loading topology file {} and trajectory file(s) {}.".format(psfFile,dcdFiles))

dnap.loadSystem(psfFile,dcdFiles)

print("System loaded.")

# We can access the trajectory data directly.
print("MDAnalysis universe:",dnap.getU().trajectory)

######################################################################
######### Prepare system for network calculations  #######
######################################################################

dnap.checkSystem()

# In this case, we are only selecting the protein for analysis. This will exclude
# all water and ions, and any ligands from the system.
selectionStr = "segid ENZY and (not (name H* or name [123]H*))"

dnap.selectSystem(userSelStr=selectionStr)

dnap.prepareNetwork()

print("Aligning trajectory...")
dnap.alignTraj()

######################################################################
######### Determine contacts and calculate correlations  #######
######################################################################

print("Finding contacts...")
dnap.findContacts(stride=1, verbose=True)

# This may be necessary for systems with low default recursion limits.
print("Recursion limit:",sys.getrecursionlimit())
sys.setrecursionlimit(3000)
print("New recursion limit:",sys.getrecursionlimit())

print("Filtering contacts...")
dnap.filterContacts(notSameRes=True, notConsecutiveRes=False, removeIsolatedNodes=True)

dnap.calcCor(ncores=1)

######################################################################
######### Determine cartesian distances  #######
######################################################################

dnap.calcCartesian(backend="serial")

######################################################################
######### Determine Network Properties #######
######################################################################

dnap.calcGraphInfo()

# Basic information of the network as interpreted as a graph.
print("Graph with {} nodes and {} edges".format(len(dnap.nxGraphs[0].nodes), len(dnap.nxGraphs[0].edges)) )

# Both density and transitivity are scaled from 0 to 1
for win in range(dnap.numWinds):
    print("----- Window {} -----".format(win))
    print("Density:", round( nx.density(dnap.nxGraphs[win]), 4) )
    print("Transitivity:", round( nx.transitivity(dnap.nxGraphs[win]), 4) )
    print()

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
dnap.calcOptPaths(ncores=1)

print("Calculating edge betweeness...")
# calculate betweeness values
dnap.calcBetween(ncores=1)

print("Here are the top 5 pairs of nodes based on Betweeness values, compared to their correlation values (in Window 0):")
for k,v in islice(dnap.btws[0].items(),5):
    print("\tNodes {} have betweenes {} and correlation {}.".format(k, 
                                                                  round(v,3), 
                                                                  round(dnap.corrMatAll[0, k[0], k[1]], 3) ) )

dnap.calcEigenCentral()

dnap.calcCommunities()

print("Here are the top 5 communities based on number of nodes:")
# Sort communities based on number of nodes
for comIndx in islice(dnap.nodesComm[0]["commOrderSize"], 5):
    print("Modularity Class {0:>2}: {1:>3} nodes.".format(comIndx, len(dnap.nodesComm[0]["commNodes"][comIndx])))

print("Here are the top 5 communities based on Eigenvector Centrality:")
# Sort communities based on the node with highest eigenvector centrality
for comIndx in islice(dnap.nodesComm[0]["commOrderEigenCentr"], 5):
    print("Modularity Class {0} ({1} nodes) Sorted by Eigenvector Centrality:".format(
                                                                    comIndx, 
                                                                len(dnap.nodesComm[0]["commNodes"][comIndx])))
    for node in dnap.nodesComm[0]["commNodes"][comIndx][:5]:
        print("Name: {0:>4} | Degree: {1:>2} | Eigenvector Centrality: {2}".format(
            node, dnap.nxGraphs[win].nodes[node]['degree'], dnap.nxGraphs[win].nodes[node]['eigenvector']))
    print()

######################################################################
######### Save data for analysis and plots #######
######################################################################

dnap.saveData(fullPathRoot)

# This function will save a reduced DCD trajectory with the heavy atoms used for network analysis
# A smaller trajectory can be created by choosing a "stride" that sub-samples the original trajectory.
# This function will also produce a PDB file so that information on atoms and residues can be loaded to
#    visualization software such as VMD.

dcdstride = 1

print("We will save {} heavy atoms and {} frames.".format(dnap.workU.atoms.n_atoms, 
                                                          len(dnap.workU.trajectory[::dcdstride]) ))

dnap.saveReducedTraj(fullPathRoot, stride = dcdstride)

print("It is normal to find warnings from MDanalysis due to missing PDB attributes, such as occupancies. This is not important for DyNetAn since we are mainly looking for atom IDs and positions. The rest of the relevant information has been stored in a topology file.")

######################################################################
######### END OF TUTORIAL #######
######################################################################

print("\n\n END OF TUTORIAL \n\n")
