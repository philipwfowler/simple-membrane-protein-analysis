# 1. Calculate how many lipid atoms are within a set distance of the protein over time

Why do we want to do this? Well, often one of the trickiest things in any membrane protein simulation is getting the protein in the membrane in the first place. There are lots of different ways of doing this and they all perturb the membrane in some way. So we need to check that the protein has become adequately embedded in the membrane during the simulation. Note that this is *separate* to looking at whether the structure of the protein has equilibrated - we shall not cover that here. This analysis is isn't especially sophisticated, we just want to count how many lipid atoms are within a set distance of any protein atom. You can choose the distance to be what you want, some people use 0.3 nm, or you could try and base it on a physical argument and for example use 0.36 nm which is about 150% the van der Waals radius of a hydrogen atom. We'll do this for every frame and then plot a graph and see if it reaches a plateau.

So, in [pseudo-code](http://xkcd.com/1185/) we want to do something like 

    load trajectory
    open a file for writing the data to
    loop over all frames
        identify all lipid atoms within 0.36nm of any protein atom
        count the lipid atoms
        find out what frame we are on
        write out the frame number and number of lipid atoms to file
    close the file
    
I'm going to show you how to do this in both [VMD](http://www.ks.uiuc.edu/Research/vmd/) and [MDAnalysis](https://code.google.com/p/mdanalysis/) to demonstrate that simple tasks can often be done using both packages and there are advantages and disadvantages to each. 

## 1.1 MDAnalysis

Let's start with python and [MDAnalysis](https://code.google.com/p/mdanalysis/). First move to the `examples/` subdirectory.

	cd examples/

One advantage of [MDAnalysis](https://code.google.com/p/mdanalysis/) is we can prototype our method in [ipython](http://ipython.org) (if you have it installed). So

	$ ipython
	Python 2.7.6 (default, Nov 18 2013, 15:12:51) 
	Type "copyright", "credits" or "license" for more information.
	
	IPython 2.0.0 -- An enhanced Interactive Python.
	?         -> Introduction and overview of IPython's features.
	%quickref -> Quick reference.
	help      -> Python's own help system.
	object?   -> Details about 'object', use 'object??' for extra details.
	
	In [1]:

First we need to load the [MDAnalysis](https://code.google.com/p/mdanalysis/) module

	In [2]: import MDAnalysis
	
Now we follow the pseudocode and load the trajectory

	In [3]: u = MDAnalysis.Universe("trajectory-files/peptso-1a.gro","trajectory-files/peptso-1a-100ns-dt1ns.xtc")

All the coordinates across all the frames are now stored in the MDAnalysis [Universe](https://code.google.com/p/mdanalysis/wiki/Universe) object which I've called u. 

	In [4]: u
	Out[4]: <Universe with 76245 atoms>

To see what methods and properties u has type u. and then hit `TAB` in ipython 

	In [5]: u.
	u.SYSTEM       u.coord        u.load_new     u.selectAtoms  
	u.atoms        u.dimensions   u.residues     u.trajectory   
	u.bonds        u.filename     u.segments     u.universe	

Since we are doing this interactively to begin with, it is easiest to just work on the current frame and then when we've got our logic sorted wrap a loop around it. Let's check we are on the first frame.

	In [5]: u.trajectory.frame
	Out[5]: 1

Ok, good. So to do this we need to select the atoms which in [MDAnalysis](https://code.google.com/p/mdanalysis/) is done using the `selectAtoms` universe method.

	In [6]: lipids = u.selectAtoms("resname POPC and around 3.6 protein")

	In [7]: lipids
	Out[7]: <AtomGroup with 1619 atoms>

The syntax that `selectAtoms` understands is based on CHARMM and, as we shall see, is subtly different to that understood by [VMD](http://www.ks.uiuc.edu/Research/vmd/). Although we can understand from the output that we have found 1619 atoms, that is not in a format that we could write to disc. So, let's look what methods and properties lipids has

	In [8]: lipids.

I haven't put the output as it would take up a lot of space. You should see a method called `lipids.numberOfAtoms`. Let's see what that does

	In [8]: lipids.numberOfAtoms()
	Out[8]: 1619

Hurray, an integer - we know how to handle them. That is basically it; we know what frame we are on and how to count the number of lipids within a specified distance of a protein. In the supplied python script ([1-count-lipids-mdanalysis.py](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/examples/1-count-lipids-mda.py)) I've added some comments that match the pseudo-code, the loop over all frames and the logic to write the numbers to both a file and the screen (i.e. STDOUT).  You can run it either by typing

	python 1-count-lipids-mdanalysis.py

or, since I have made the file executable using `chmod`  you can run it directly 

	./1-count-lipids-mdanalysis.py

It should write the data to the screen (STDOUT) and also to a file called `dat/1-count-lipids-mda.dat`. Now you can plot these data using your preferred graphing package (I'd recommend gnuplot). Does the protein appear to be embedded to you?

![Graph of something](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/images/graph-1-count-lipids-mda.png)	

## 1.2 VMD

Now let's try and do exactly the same thing in [VMD](http://www.ks.uiuc.edu/Research/vmd/). Although you may have done some analysis in [VMD](http://www.ks.uiuc.edu/Research/vmd/) using a bespoke GUI we are not going to do that. We are going to write some code in `Tcl`, which is the scripting language [VMD](http://www.ks.uiuc.edu/Research/vmd/) exposes that lets you do this sort of thing. `Tcl` can feel a bit strange if you are used to programming `C`, `python`, `perl` etc as, for example, it has no equals sign. You'll see what I mean. First, load [VMD](http://www.ks.uiuc.edu/Research/vmd/), either by clicking the icon or typing `vmd` into your terminal.

It is a bit harder to prototype in [VMD](http://www.ks.uiuc.edu/Research/vmd/), but let's try it. You'll need to have the `Tk Console` window open: go to `Extensions` then `Tk Console` and something that looks a bit like a Terminal will appear, but with coloured text on a white background. Again, we shall just follow the pseudo-code, so first, let's load in the trajectory.

	(examples) 1 % mol load gro trajectory-files/peptso-1a.gro	0	>Main< (examples) 2 % mol addfile trajectory-files/peptso-1a-100ns-dt1ns.xtc type xtc first 0 last -1 waitfor all	0

You should now see see the waters, lipids and protein in the VMD Display. Seeing what is going on is a key advantage of using [VMD](http://www.ks.uiuc.edu/Research/vmd/). Again, we shall just work on one frame and wrap everything in a loop later. Here we bump into an annoying inconsistency between the two packages: [VMD](http://www.ks.uiuc.edu/Research/vmd/) loads the GRO/PDB file and makes that frame 0, then it loads the trajectory (frames 1-101) and it leaves you at the last frame, 101. [MDAnalysis](https://code.google.com/p/mdanalysis/) just uses the GRO/PDB file to parse the trajectory which is stored in  frames 1-101 and you start at frame 1. For consistency I'm going to return to the start using the arrows in the `Main` window (or you could type `animate goto 1` in the `Tk Console`). Now to see how many lipid atoms are within 0.36 nm of the protein:

	>Main< (examples) 7 % set lipids [atomselect top "resname POPC and within 3.6 of protein"]	atomselect1	>Main< (examples) 8 % $lipids num	1619

The same answer as before which is encouraging! Now we can put the logic into a `Tcl` script, add the loop and file output and some comments, as before. To run the script ([1-count-lipids-vmd.tcl](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/examples/1-count-lipids-vmd.tcl)), either type the following the `Tk Console`.

	>Main< (examples) 13 % source 1-count-lipids-vmd.tcl

Or, close [VMD](http://www.ks.uiuc.edu/Research/vmd/) (type `quit` in the `Tk Console`) and from the terminal you can launch [VMD](http://www.ks.uiuc.edu/Research/vmd/) without a Display and load the Tcl file

	$ vmd -dispdev text -e 1-count-lipids-vmd.tcl

Again, you can plot the data to see if the protein is embedded: 

![Graph of something](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/images/graph-1-count-lipids-vmd.png)


It looks the same as we got using [MDAnalysis](https://code.google.com/p/mdanalysis/), but how similar are the results? Try using `diff` which is a very useful Linux command that compares two text files.

	$ diff dat/1-count-lipids-mda.dat dat/1-count-lipids-vmd.dat 
	0a1
	>       0     730

Yes both methods produce identical results - the only difference is that, as we expected, the [VMD](http://www.ks.uiuc.edu/Research/vmd/) file contains a frame 0 which is the results for the GRO file.

## 1.3 Extensions

- Which is faster? Why do you think that is? Do you think it will always be faster?
- Write down the disadvantages and advantages of each approach
- Convert the script to work on a [coarse-grained MARTINI simulation](http://md.chem.rug.nl/cgmartini/). What distance should you use then?
- Add command line arguments (using an appropriate python module like [`getopt`](https://docs.python.org/2/library/getopt.html)) to generalise the [MDAnalysis](https://code.google.com/p/mdanalysis/) version. Turn the following into arguments:
	- the name/path of the input PDB/GRO file
	- the name/path of the input trajectory file
	- the distance 
	- the name/path of the output file
	- and perhaps the string used to identify the lipids (here it is "resname POPC")

Now we are ready to move onto the [next exercise](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/2-protein-depth.md).	