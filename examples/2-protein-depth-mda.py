#! /usr/bin/env python

# load the MDAnalysis and numpy modules
import MDAnalysis
import numpy

# load trajectory
u = MDAnalysis.Universe("trajectory-files/peptso-1a.gro","trajectory-files/peptso-1a-100ns-dt1ns.xtc")

# open a file for writing the data to
OUTPUT = open("dat/2-protein-depth-mda.dat",'w')

# loop over all frames
for frame in u.trajectory:

    # identify the protein
    protein = u.selectAtoms("protein")

    # identify the bilayer
    bilayer = u.selectAtoms("resname POPC")

    # find out the z coordinate of the centre of mass of the protein
    proteinZ = protein.centerOfMass()[2]
    
    # find out the z coordinate of the centre of mass of the bilayer
    bilayerZ = bilayer.centerOfMass()[2]
    
    # calculate the relative depth of the protein
    proteinCoord = proteinZ - bilayerZ
    
    # find out what frame we are on
    frame = u.trajectory.frame

    # so we can see what is happening, write the output to the screen for the moment
    print frame,proteinCoord

    # write out the frame number and the relative depth of the protein to the file (the bit in the brackets formats the data nicely)
    print >> OUTPUT, "%7i %7.3f" % (frame,proteinCoord)
    
# close the file    
OUTPUT.close