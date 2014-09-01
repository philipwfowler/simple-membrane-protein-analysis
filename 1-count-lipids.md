1. Calculate how many lipid atoms are within a set distance of the protein over time
================================

Why do we want to do this? Well, often one of the trickiest things in any membrane protein simulation is getting the protein in the membrane in the first place. There are lots of different ways of doing this and they all perturb the membrane in some way. So we need to check that the protein has become adequately embedded in the membrane during the simulation. Note that this is *separate* to looking at whether the structure of the protein has equilibrated - we shall not cover that here. This is isn't especially sophisticated, we just want to count how many lipid atoms are within a set distance of any protein atom. You can choose the distance to be what you want, some people use 0.3 nm, or you could use 0.36 nm which is about 150% the van der Waals radius of a hydrogen atom. We'll do this every timestep and then plot a graph and see if it reaches a plateau.

So, in pseudo-code we want to do something like 

    load trajectory
    open a file to store the results
    loop over all frames
        select all lipid atoms within 0.36nm of a protein atom
        count the lipid atoms
        find out what frame we are on
        write out the frame number and number of lipid atoms to disc
    close the file
    
I'm going to show you how to do this in both VMD and MDAnalysis to show you that simple tasks can often be done more than one way and there are advantages and disadvantages to each. 

**1.1 MDAnalysis**

Let's start with python and MDAnalysis. One advantage of MDAnalysis is we can prototype in python (if you have it installed). So

	$ ipython
	Python 2.7.6 (default, Nov 18 2013, 15:12:51) 
	Type "copyright", "credits" or "license" for more information.
	
	IPython 2.0.0 -- An enhanced Interactive Python.
	?         -> Introduction and overview of IPython's features.
	%quickref -> Quick reference.
	help      -> Python's own help system.
	object?   -> Details about 'object', use 'object??' for extra details.
	
	In [1]:

First we need to load the MDAnalysis module

	In [2]: import MDAnalysis
	
Now we follow the pseudocode and load the trajectory

	In [3]: u = MDAnalysis.Universe("example-trajectory/peptso-1a.gro","example-trajectory/peptso-1a-100ns-dt1ns.xtc"

All the coordinates across all the frames are now stored in the MDAnalysis Universe object which I've called u. 

	In [4]: u
	Out[4]: <Universe with 76245 atoms>

To see what methods and properties u has type u. and then hit TAB in ipython 

	In [5]: u.
	u.SYSTEM       u.coord        u.load_new     u.selectAtoms  
	u.atoms        u.dimensions   u.residues     u.trajectory   
	u.bonds        u.filename     u.segments     u.universe	

Since we are doing this interactively to begin with, it is easiest to just work on the current frame and then when we've got our logic sorted wrap a loop around it. Let's check we are on the first frame.

	In [5]: u.trajectory.frame
	Out[5]: 1

Ok, good. So to do this we need to select the atoms which in MDAnalysis is done using the `selectAtoms` universe method.

	In [8]: lipids = u.selectAtoms("resname POPC and around 3.6 protein")

	In [9]: lipids
	Out[9]: <AtomGroup with 1619 atoms>

The syntax that `selectAtoms` understands is based on CHARMM and, as we shall see, is subtly different to that understood by VMD. Although we can understand from the output that we have found 1619 atoms, that is not in a format that we could write to disc. So, let's look what methods and properties lipids has

	In [12]: lipids.

I haven't put the output as it would take up a lot of space. You should see a method called `lipids.numberOfAtoms`. Let's see what that does

	In [15]: lipids.numberOfAtoms()
	Out[15]: 1619

Hurray, an integer. That is basically it; we know what frame we are on and how to count the number of lipids within a specified distance of a protein. In the supplied python script (1-count-lipids-mdanalysis.py) I've just added some comments, the loop over time steps and the logic to write the numbers to a file and the screen (i.e. STDOUT).  you can run it either by typing

	python 1-count-lipids-mdanalysis.py

or, since I have made the file executable using `chmod`  you can run it directly 

	./1-count-lipids-mdanalysis.py

It should write the data to the screen (STDOUT)	

**1.2 VMD**

Now let's try and do exactly the same thing in VMD. 

It is a bit harder to prototype here, but let's give it a shot. 

Although some analysis is done in VMD using little GUIs we are not going to do that. We are going to write some code in Tcl, which is the script language VMD exposes that lets you do this sort of thing. Tcl can feel a bit strange if you are used to programming C, python, perl etc as it has no equals sign. You'll see what I mean. First, load VMD, either by clicking the icon or typing `vmd` into your terminal.

**1.3 Extension exercises**

- Which is faster? Why do you think that is? Do you think it will always be faster?
- Write down the disadvantages and advantages of each approach
- Add command line arguments (using an appropriate python module like getopt) to generalise the MDAnalysis version. Turn the following into arguments:
	- the name/path of the input PDB/GRO file
	- the name/path of the input trajectory file
	- the distance 
	- the name/path of the output file
	- and perhaps the string used to identify the lipids (here it is "resname POPC"")