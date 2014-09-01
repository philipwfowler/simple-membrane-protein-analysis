#! /usr/bin/env python

# load the MDAnalysis and numpy modules
import MDAnalysis
import numpy

# load the trajectory
u = MDAnalysis.Universe("trajectory-files/peptso-1a.gro","trajectory-files/peptso-1a-100ns-dt1ns.xtc")

# define the species using a dictionary (will make it easier to generalise)
speciesList = {'phosphate':'name P1','lipids':'resname POPC','headgroups':'resname POPC and name N','protein':'protein'}

# make a list of the bin edges
bins = range(-80,80,2)

# what frames do we want to analyse?
start = 50
end = 100

#loop over the molecular species considered  
for species in speciesList:

    # create an empty list
    coordinates = []

    # open file to store the results
    OUTPUT = open("dat/2-density-profile-" + species + "-mda.dat",'w')

    # loop over all frames
    for timestep in u.trajectory:

        # find out frame we are on
        frame = u.trajectory.frame 

        # check that the frame lies within the region we wish to analyse
        if frame >= start and frame < end: 

            # choose the atoms 
            atoms = u.selectAtoms(speciesList[species])

            # select the bilayer
            bilayer = u.selectAtoms("resname POPC")

            # find out the centre of mass of the bilayer along z (needed for normalisation)
            bilayerZ = bilayer.centerOfMass()[2]

            # loop over the number of species
            for i in atoms:

                # calculate the normalised value of z (i.e. relative to the centre of the bilayer)
                z = i.position[2] - bilayerZ

                # append it to the list
                coordinates.append(z)

    # now, using numpy, calculate the histogram
    values,edges = numpy.histogram(coordinates,bins)    

    # and write it to a file (zip is a bit of python magic that iterates through two lists at the same time)
    for (i,j) in zip(edges,values):
        print >> OUTPUT, "%5i %7.3e" % (i,j)
            
    # close the file    
    OUTPUT.close