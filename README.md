# simple-membrane-protein-analysis

This is a short tutorial showing how to use both [VMD](http://www.ks.uiuc.edu/Research/vmd/) and [MDAnalysis](https://code.google.com/p/mdanalysis/) to analyse a molecular dynamics simulation of a membrane protein. A sample trajectory is provided (of a bacterial peptide transporter which, structurally, belongs to the Major Facilitator Superfamily if you are interested). Although the trajectory we shall use was generated using GROMACS, all the examples can, with a few tweaks, could be applied to trajectories produced by other packages, such as NAMD, CHARMM or AMBER.

The aims are:
- to show you how to run some simple analysis on a membrane protein molecular dynamics simulation
- to encourage you to write your own analysis code, rather than use a 'script' that someone gives you 
- to show you how simple analysis tasks can often be done more than one way. Here I shall use either [VMD](http://www.ks.uiuc.edu/Research/vmd/) or [MDAnalysis](https://code.google.com/p/mdanalysis/)
- to introduce you briefly to the importance of writing good code and [Software Carpentry](http://software-carpentry.org/index.html). Specifically to
    - test your code
    - include helpful comments
    - generalise your code so it can be used at different occasions
    - fit into a (semi)-automated workflow

We won't be
- exhaustively looking at all the ways you can analysis membrane proteins
- categorically saying which approach is 'best'

You need to understand something about
- membrane proteins
- molecular dynamics simulations
- writing computer code (any language)
- ideally 
	- be familiar with VMD
	- be familiar with simple python


There are currently four sections to this tutorial
1. calculate how many lipid atoms are within a set distance of the protein over time
2. measure the depth of the centre of mass of the protein relative to the centre of mass of the bilayer over time
3. calculate the density profile of the different species along the membrane normal
4. determine which lipids belong to which leaflet

You will need to have installed
- git
- python 2.6.X or 2.7.X, including
    - numpy
    - scipy
    - MDAnalysis
- [VMD](http://www.ks.uiuc.edu/Research/vmd/)
- a graphing package. I like [gnuplot](http://gnuplot.sourceforge.net) but XMGrace or Excel are fine too.

To get started you will need to clone this git repository (collection of files in a folder) onto your computer. If you are using a Mac or a Linux machine, open a Terminal and check you have git by typing

    git

Assuming git is installed, you should then copy the address of this repository to your clipboard - click the button with the arrow pointing to the clipboard just below the textbox which has *HTTPS* above it on the right hand side of this page. Then on your machine in the terminal type

    git clone https://github.com/philipwfowler/simple-membrane-protein-analysis.git
    
and, volia, you should have a new folder called simple-membrane-protein-analysis/ which contains this file (README.md) along with all the others. 

If you don't have git installed, then you will either need to install it or  try using the appropriate [github GUI](http://git-scm.com/downloads/guis) for your platform (includes Windows). If all else fails, you can always download a zip of the repo using the "Download ZIP" button which you can find in the bottom right of this page.

