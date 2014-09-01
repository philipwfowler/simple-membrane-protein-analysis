# 2. Measure the depth of the centre of mass of the protein relative to the centre of mass of the bilayer over time

Let's now look at how the protein is moving in the bilayer, specifically whether it is sinking or rising relative to the membrane normal. Caveat: we are going to *assume* that the membrane normal is aligned with the *z* axis (this is common practice). Our patch of bilayer is fairly small, so this is reasonable but it isn't always true! Let's think how the pseudocode might look

    load trajectory
    open a file to store the results
    loop over all frames
        find out the z coordinate of the centre of mass of the protein
        find out the z coordinate of the centre of mass of the bilayer
        calculate the difference; this is the relative depth of the protein
        write out the frame number and the relative depth of the protein to disc
    close the file

Like before, we shall do this analysis using both VMD and python/MDAnalysis and compare the results. To avoid too many gear-changes in our brains, let's do VMD first and then try MDAnalysis. Fire up VMD as before and open the `Tk Console` window. First we need to load the trajectory, just as before.

	(examples) 1 % mol load gro trajectory-files/peptso-1a.gro	0	>Main< (examples) 2 % mol addfile trajectory-files/peptso-1a-100ns-dt1ns.xtc type xtc first 0 last -1 waitfor all	0

VMD has a series of analysis commands that start with the keyword [`measure`](http://www.ks.uiuc.edu/Research/vmd/current/ug/node136.html). They tend to be quite fast as they have (mostly) been optimised in C and some of them will even run on your GPU, if you have the right drivers etc installed. The one we are going to use is `measure centre` which returns the geometric centre of an atom selection.

	>Main< (examples) 4 % set protein [atomselect top "protein"]	atomselect1	>Main< (examples) 5 % measure center $protein weight mass	39.187652587890625 37.207035064697266 52.02757263183594

Let's store the centre of mass in a variable.

	>Main< (examples) 8 % set proteinCoord [measure center $protein weight mass]	39.187652587890625 37.207035064697266 52.02757263183594

Now Tcl being Tcl, lists are a bit different. To pick out the third value we shall use [lindex](http://www.tcl.tk/man/tcl8.4/TclCmd/lindex.htm) (remember, like python, lists in Tcl are [0-based](http://www.xkcd.com/163/)).

	>Main< (examples) 10 % lindex $proteinCoord 2	52.02757263183594

That's the kernel. 	
