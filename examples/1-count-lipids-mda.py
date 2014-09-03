#! /usr/bin/env python

# load the MDAnalysis module
import MDAnalysis

# load trajectory
u = MDAnalysis.Universe("trajectory-files/peptso-1a.gro","trajectory-files/peptso-1a-100ns-dt1ns.xtc")

# open a file for writing the data to
OUTPUT = open("dat/1-count-lipids-mda.dat",'w')

# loop over all frames
for timestep in u.trajectory:

    # identify all lipid atoms within 0.36nm of any protein atom
    lipids = u.selectAtoms("resname POPC and around 3.6 protein")
    
    # count the lipid atoms
    lipidNumber = lipids.numberOfAtoms()

    # find out what frame we are on
    frame = u.trajectory.frame

    # so we can see what is happening, write the output to the screen for the moment
    print frame,lipidNumber

    # write out the frame number and number of lipid atoms to file (the bit in the brackets formats the data nicely)
    print >> OUTPUT, "%7i %7i" % (frame,lipidNumber)
    
# close the file    
OUTPUT.close