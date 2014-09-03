# simple-membrane-protein-analysis

This is a short tutorial that will show you how to use both [VMD](http://www.ks.uiuc.edu/Research/vmd/) and [MDAnalysis](https://code.google.com/p/mdanalysis/) to analyse a molecular dynamics simulation of a membrane protein. A sample trajectory is provided (of [a bacterial peptide transporter](http://emboj.embopress.org/content/30/2/417) which, structurally, belongs to the [Major Facilitator Superfamily](http://en.wikipedia.org/wiki/Major_facilitator_superfamily) if you are interested). Although the trajectory we shall use was generated using [GROMACS](http://www.gromacs.org), all the examples can, with a few tweaks, could be applied to trajectories produced by other packages, such as NAMD, CHARMM or AMBER.

The aims are:
- to show you how to run some simple analysis on a membrane protein molecular dynamics simulation
- to encourage you to write your own analysis code, rather than use a 'script' that someone gives you 
- to show you how simple analysis tasks can often be done more than one way. Here I shall use both [VMD](http://www.ks.uiuc.edu/Research/vmd/) and [MDAnalysis](https://code.google.com/p/mdanalysis/) to do the same analysis.
- to introduce you briefly to the importance of writing [good code](http://www.xkcd.com/844/) and [Software Carpentry](http://software-carpentry.org/index.html), specifically to
    - test your code
    - include helpful comments
    - generalise your code so it can be more easily used again
    - fit into a (semi)-automated workflow

We won't be
- exhaustively looking at all the ways you can analysis membrane proteins
- categorically saying which approach is 'best'

You need to understand something about
- membrane proteins
- molecular dynamics simulations
- writing computer code (any language)
- ideally 
	- be familiar with [VMD](http://www.ks.uiuc.edu/Research/vmd/)
	- be familiar with simple python

There are currently three sections to this tutorial. In each of the first two we shall run a simple analysis task and in the third we shall compare VMD and MDAnalysis and list some of the advantages and disadvantages of each to analyse protein simulations.

You will need to have installed
- [git](http://git-scm.com)
- python 2.6.X or 2.7.X, including
    - [numpy](http://www.numpy.org)
    - [scipy](http://www.scipy.org)
    - [MDAnalysis](https://code.google.com/p/mdanalysis/)
- [VMD](http://www.ks.uiuc.edu/Research/vmd/)
- a graphing package. I like [gnuplot](http://gnuplot.sourceforge.net) but XMGrace or Excel are fine too.

To get started you will need to clone this [git](http://git-scm.com) repository (collection of files in a folder) onto your computer. If you are using a Mac or a Linux machine, open a Terminal and check you have [git](http://git-scm.com) by typing

    git

Assuming [git](http://git-scm.com) is installed, you should then copy the address of this repository to your clipboard - click the button with the arrow pointing to the clipboard just below the textbox which has *HTTPS* above it on the right hand side of this page. Then on your machine in the terminal type (or just copy the below)

    git clone https://github.com/philipwfowler/simple-membrane-protein-analysis.git
    
and, volia, you should have a new folder called simple-membrane-protein-analysis/ which contains this file (README.md) along with all the others. 

If you don't have [git](http://git-scm.com) installed, then you will either need to install it or  try using the appropriate [github GUI](http://git-scm.com/downloads/guis) for your platform (includes Windows). If all else fails, you can always download a zip of the repo using the "Download ZIP" button which you can find in the bottom right of this page but this will just give you a copy of the files.

In the first example we shall count the number of lipids within a set distance of the protein during a simulation to assess when the system has equilibrated.  Click [here](https://github.com/philipwfowler/simple-membrane-protein-analysis/blob/master/1-count-lipids.md) to go this example