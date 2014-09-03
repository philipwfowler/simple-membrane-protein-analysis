# 2. Measure the depth of the centre of mass of the protein relative to the centre of mass of the bilayer over time

Let's now look at how the protein is moving in the bilayer, specifically whether it is sinking or rising relative to the plane of the membrane. Caveat: we are going to *assume* that the membrane normal is aligned with the *z* axis (this is common practice). Our patch of bilayer is fairly small, so this is reasonable but it isn't always true! Let's start by thinking how the pseudocode might look

    load trajectory
    open a file for writing the data to
    loop over all frames
        find out the z coordinate of the centre of mass of the protein
        find out the z coordinate of the centre of mass of the bilayer
        calculate the relative depth of the protein
        find out what frame we are on
        write out the frame number and the relative depth of the protein to the file
    close the file

Like before, we shall do this analysis using both [VMD](http://www.ks.uiuc.edu/Research/vmd/) and python/[MDAnalysis](https://code.google.com/p/mdanalysis/) and compare the results. To avoid too many gear-changes in our brains, let's do [VMD](http://www.ks.uiuc.edu/Research/vmd/) first and then try [MDAnalysis](https://code.google.com/p/mdanalysis/). 

## 2.1 VMD 

Fire up [VMD](http://www.ks.uiuc.edu/Research/vmd/) as before and open the `Tk Console` window. First we need to load the trajectory, just as before.

	(examples) 1 % mol load gro trajectory-files/peptso-1a.gro	0	>Main< (examples) 2 % mol addfile trajectory-files/peptso-1a-100ns-dt1ns.xtc type xtc first 0 last -1 waitfor all	0

[VMD](http://www.ks.uiuc.edu/Research/vmd/) has a series of analysis commands that start with the keyword [`measure`](http://www.ks.uiuc.edu/Research/vmd/current/ug/node136.html). They tend to be quite fast as they have (mostly) been optimised in C and some of them will even run on your GPU, if you have the right drivers etc installed. The one we are going to use is `measure centre` which returns the geometric centre of an atom selection.

	>Main< (examples) 4 % set protein [atomselect top "protein"]	atomselect1 	>Main< (examples) 5 % measure center $protein weight mass	42.8310661315918 40.41354751586914 49.69547653198242

Let's store the centre of mass in a variable.

	>Main< (examples) 8 % set proteinCoord [measure center $protein weight mass]	42.8310661315918 40.41354751586914 49.69547653198242

Now Tcl being Tcl, lists are a bit different. To pick out the third value we shall use [lindex](http://www.tcl.tk/man/tcl8.4/TclCmd/lindex.htm) (remember, like python, lists in Tcl are [0-based](http://www.xkcd.com/163/)).

	>Main< (examples) 10 % lindex $proteinCoord 2	49.69547653198242

That's the kernel of what we need to this in [VMD](http://www.ks.uiuc.edu/Research/vmd/) (we can simply repeat the commands for the bilayer and subtract one from the other to get the relative depth of the protein). Now let's write a program ([2-protein-depth-vmd.tcl](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/examples/2-protein-depth-vmd.tcl)). 	
You can run the program as before by either loading it directly into the `Tk Console` window by typing

	>Main< (examples) 11 % source 2-protein-depth-vmd.tcl
	
Or in the Terminal loading [VMD](http://www.ks.uiuc.edu/Research/vmd/)

	$  vmd -dispdev text -e 2-protein-depth-vmd.tcl 
	
Plotting the resulting data shows the protein is moving by a few Angstroms relative to the bilayer. We wouldn't expect a transmembrane protein to move that much (a peripheral protein might though). At this stage we would probably examine different parts of the protein, for example, we might look to see if there were any differences between the N- and C-terminal halves of the protein, or between individual transmembrane helices.

![Graph of something](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/images/graph-2-protein-depth-vmd.png)

## 2.2 MDAnalysis.

Back to python (well, ipython to start with). Let's do the same calculation using [MDAnalysis](https://code.google.com/p/mdanalysis/). Again we need to import the module and load the coordinates.

	In [1]: import MDAnalysis

	In [2]: u = MDAnalysis.Universe("example-trajectory/peptso-1a.gro","example-trajectory/peptso-1a-100ns-dt1ns.xtc"

Again, we will just look at the first frame of the trajectory. First we need to identify the protein.

	In [3]: protein = u.selectAtoms("protein")

Let's see if protein has a helpful method. Type `protein.` in the ipython terminal and hit TAB. You should get a long list. One in particular looks interesting! Try out

	In [4]: protein.centerOfMass()
	Out[4]: array([ 42.83182501,  40.41418755,  49.6961789 ])

Remember we just want the *z* value which is the third element which, since this is python, we access in a more familiar way.

	In [5]: protein.centerOfMass()[2]
	Out[5]: 49.696178896249158

That is basically all we need. Now we can write a python program that puts all this together ([2-protein-depth-mda.py](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/examples/2-protein-depth-mda.py)).

If you run this program and plot the results you should get something that looks similar to the results we got using VMD.

![Graph of something](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/images/graph-2-protein-depth-mda.png)

But are they the same? Analysing the files using `diff`

	diff dat/2-protein-depth-mda.dat dat/2-protein-depth-vmd.dat
	
we find that they are different, but not very different (typically at the second decimal place). One way to see how different is to plot the two datasets, one on the x-axis and the other on the y-axis. 

![Graph of something](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/images/graph-2-protein-depth-comparison.png)

Do you think this discrepancy is important? What do you think is responsible? What might you do to check your hypothesis?

## 2.3 Extension exercises

- Demonstrate the cause of the discrepancy above.
- Generalise the python approach by adding appropriate command line flags (perhaps using the [`getopt`](https://docs.python.org/2/library/getopt.html) module`)
- Examine the relative motion of different components of the protein e.g. 1st six transmembrane helices versus the last six. How will you identify where the helices start and finish?
- Is selecting all the lipid atoms the best approach? What if we only considered the tails? How might we select these in both [VMD](http://www.ks.uiuc.edu/Research/vmd/) and [MDAnalysis](https://code.google.com/p/mdanalysis/)?