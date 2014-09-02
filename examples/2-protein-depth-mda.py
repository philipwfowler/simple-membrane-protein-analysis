#! /usr/bin/env python

# load the MDAnalysis and numpy modules
import MDAnalysis
import numpy

# load the trajectory
u = MDAnalysis.Universe("trajectory-files/peptso-1a.gro","trajectory-files/peptso-1a-100ns-dt1ns.xtc")

# open a file for writing the data to
OUTPUT = open("dat/2-protein-depth-mda.dat",'w')

# iterate through the trajectory, frame by frame
for timestep in u.trajectory:

    # identify the protein
    protein = u.selectAtoms("protein")

    # identify the bilayer
    bilayer = u.selectAtoms("resname POPC")

    proteinZ = protein.centerOfMass()[2]
    bilayerZ = bilayer.centerOfMass()[2]
    
    proteinCoord = proteinZ - bilayerZ
    
    # what frame are we on?
    frame = u.trajectory.frame

    # so we can see what is happening, write the output to the screen for the moment
    print frame,proteinCoord

    # write the data for this timestep to the file (the bit in the brackets formats the data nicely)
    print >> OUTPUT, "%7i %7.3f" % (frame,proteinCoord)
    
# close the file    
OUTPUT.close