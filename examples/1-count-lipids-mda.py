#! /usr/bin/env python

# load the MDAnalysis module
import MDAnalysis

# load the trajectory
u = MDAnalysis.Universe("trajectory-files/peptso-1a.gro","trajectory-files/peptso-1a-100ns-dt1ns.xtc")

# open a file for writing the data to
OUTPUT = open("dat/1-mda-lipid-contacts.dat",'w')

# iterate through the trajectory, frame by frame
for timestep in u.trajectory:

    # identify the lipid atoms within 0.36 nm of the protein
    lipids = u.selectAtoms("resname POPC and around 3.6 protein")
    
    # how many atoms is that?
    lipidNumber = lipids.numberOfAtoms()

    # what frame are we on?
    frame = u.trajectory.frame

    # so we can see what is happening, write the output to the screen for the moment
    print frame,lipidNumber

    # write the data for this timestep to the file (the bit in the brackets formats the data nicely)
    print >> OUTPUT, "%6i %7i" % (frame,lipidNumber)
    
# close the file    
OUTPUT.close