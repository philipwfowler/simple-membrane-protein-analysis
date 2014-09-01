simple-membrane-protein-analysis
================================

This is a short tutorial showing how to use both [VMD][vmd] and [MDAnalysis][mda] to analyse a molecular dynamics simulation of a membrane protein. A sample trajectory is provided (of a bacterial peptide transporter which, structurally, belongs to the Major Facilitator Superfamily if you are interested). Although the trajectory we shall use was generated using GROMACS, all the examples can, with a few tweaks, could be applied to trajectories produced by other packages, such as NAMD, CHARMM or AMBER.

The aims are:
- to show you how to run some simple analysis on a membrane protein molecular dynamics simulation
- to encourage you to write your own analysis code, rather than use a 'script' that someone gives you 
- to show you how simple analysis tasks can often be done more than one way. Here I shall use either [VMD][vmd] or [MDAnalysis][mda]
- to introduce you briefly to the importance of writing good code and [Software Carpentry][swc]. Specifically to
    - test your code
    - include helpful comments
    - generalise your code so it can be used at different occasions
    - fit into a (semi)-automated workflow

We won't be
- exhaustively looking at all the ways you can analysis membrane proteins
- categorically saying which approach is 'best'

You need to understand something about
- membrane proteins
- molecular dynamics simulations (

You will need to have installed
- git
- python 2.6.X or 2.7.X, including
    - numpy
    - scipy
    - MDAnalysis
- [VMD][vmd]
- a graphing package. I like [gnuplot](http://gnuplot.sourceforge.net) but XMGrace or Excel are fine too.

[vmd]: http://www.ks.uiuc.edu/Research/vmd/ VMD
[mda]: https://code.google.com/p/mdanalysis/ MDAnalysis
[swc]: http://software-carpentry.org/index.html Software Carpentry
