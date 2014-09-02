# 3. Calculate the density profile of the different species along the membrane normal

Hopefully we've convinced ourselves that our transporter protein is embedded into the lipid bilayer after ~50 ns and isn't moving about a lot. Now we'd like to look at where the protein sits along the membrane normal in relation to the different components of the phospholipids. 

Let's think about what we want to do. In a nutshell, we want to build up histograms along the membrane normal (i.e. $z$) of say, the protein, the phosphate groups, the lipid tails, the headgroups and water. We might also want to split the protein down into individual transmembrane helices, but let's see what the result looks like first. So each frame, we want to measure the position along *z* of each of our chosen components, subtract the centre of mass of the membrane to normalise it and add it to a list. When we've parsed all our chosen frames, we then histogram each list and write out the results to file. So the pseudo-code looks something like

	load trajectory
	loop over the molecular species considered  
		open file to store the results
		create an empty list
		loop over all frames
			measure the coordinates of the species
			loop over the number of species
				append the coordinates to the list
		histogram the list of coordinates
		write to the file
		close the file	

This is a *lot* easier in python than in VMD because python has, in addition to MDAnalysis, some other very useful modules, specifically numpy and scipy, which provide all sorts of useful numerical and scientific functionality. So, I'm only going to show you this in python. As before, let's prototype but only look at one species (the phosphate atoms) from one frame. First let's load the trajectory.

	In [2]: import MDAnalysis
	
	In [3]: u = MDAnalysis.Universe("example-trajectory/peptso-1a.gro","example-trajectory/peptso-1a-100ns-dt1ns.xtc"

The phosphate atoms can be identified by

	In [5]: species = u.selectAtoms("name P1")

	In [6]: species
	Out[6]: <AtomGroup with 205 atoms>

Which is right as there are 205 lipids (check by running `tail trajectory-files/peptso-1a.gro` in a Terminal). As you have to append to an existing list, we shall need to create an empty list first

	In [7]: coordinates = [] 

Now we can loop over the phosphates.

	In [12]: for i in species:
	     z=i.position[2]
	     coordinates.append(z)

You need to take the 3rd (0-based!) element of position as this is *z*; the 1st and 2nd are *x* and *y*. If you now type `print coordinates` you should get a list of z values of the phosphate atoms 205 elements long (`print len(coordinates)`). Now numpy has a method that, if you give it a list of values and a list of bin edges, will calculate the histogram for you - [numpy.histogram](http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html). First we need to make a list of bin edges which should be symmetric and say twice as big as the bilayer which will be 35-40nm wide. Hence:

	In [17]: bins = range(-80,80,2)

Now we can do the numpy magic in one line (after importing numpy, of course)

	In [18]: import numpy

	In [19]: values,edges = numpy.histogram(coordinates,bins)

That is the heart of what we need; everything else is adding loops over frames and different species and writing the results to disc etc.		