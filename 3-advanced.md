# 3. More Advanced Techniques

We won't go through any problems in this section but will instead just look at some of the differences between [VMD](http://www.ks.uiuc.edu/Research/vmd/) and [MDAnalysis](https://code.google.com/p/mdanalysis/).

## 3.1 VMD

###Advantages

- You can see what is happening! For example you can check that the text for an `atomselect` command gives the result you expect by creating a Graphical Representation using the same logic.
- Many of the `measure` commands have been written in `C` so are very fast. Some have even been written to use a GPU if one is available and you have the right CUDA drivers installed.
- There are a large number of [examples](http://www.ks.uiuc.edu/Research/vmd/script_library/) 
- [VMD](http://www.ks.uiuc.edu/Research/vmd/) is available for Windows, Linux and Mac.


###Disadvantages

- If you are familiar with `C`, `perl`, `python` or other common languages then learning `Tcl` can be a bit hard. There aren't the books either or number of online resources that `python` has, for example.
- The `Tcl` implementation is part of the [VMD](http://www.ks.uiuc.edu/Research/vmd/) application installed on your machine so you are restricted to the functionality that is included. There is some math functionality but there is no linear algebra or anything more sophisticated. This also means you can't load other modules/packages.
- Because you launch [VMD](http://www.ks.uiuc.edu/Research/vmd/) to run any `Tcl` you write, it is difficult to automated any approach that uses [VMD](http://www.ks.uiuc.edu/Research/vmd/). For example, you can't add command line flags but typically have to generate the `Tcl` code each time. This can make automating an analysis workflow that uses `Tcl` difficult.
- It is more difficult to ""write good code" since there is no testing framework.

### Things it can do that others can't

- Render images
- Some analysis, for example measure the solvent accessible surface area of a protein

## 3.2 python/MDAnalysis

###Advantages

- Because [MDAnalysis](https://code.google.com/p/mdanalysis/) is a python module, you just `import` it into your preferred python implementation. So you can use regular `python` or perhaps `ipython` or even an [`ipython notebook`](http://ipython.org/notebook.html). 
- This also means that once you have the coordinates loaded in and available through the Universe object, you are free to use any other python module to help in your analysis. The most obvious ones are [numpy](http://www.numpy.org) and [scipy](http://www.scipy.org), but there are many that could be useful, e.g. [networkx](https://networkx.github.io) or even [image processing techniques](http://pubs.rsc.org/en/Content/ArticleLanding/2014/FD/c3fd00131h).
- There are [xkcd comics](http://xkcd.com/353/) about python.
- Because you can add command line flag functionality to your program it is easy to integrate it into an automated workflow on a Linux (or Mac) machine.
- python makes writing [good code](http://www.xkcd.com/844/) easy since you can 
	- add unit tests (e.g. [nose](https://nose.readthedocs.org/en/latest/))
	- automatically create documentation (e.g. [python docstrings](http://www.pythonforbeginners.com/basics/python-docstrings))

	
###Disadvantages

- Although many of the [numpy](http://www.numpy.org) methods are optimised in `C`, [MDAnalysis](https://code.google.com/p/mdanalysis/) is often slower than [VMD](http://www.ks.uiuc.edu/Research/vmd/), although this is not always true.
- It is newer and therefore less developed (at present, at least)
- There is no way of 'seeing' what is happening so you can end up using it in conjunction with [VMD](http://www.ks.uiuc.edu/Research/vmd/) anyway.

### Things it can do that others can't

- [MDAnalysis](https://code.google.com/p/mdanalysis/) has a neat function called [LeafletFinder](https://mdanalysis.googlecode.com/svn/trunk/doc/html/documentation_pages/analysis/leaflet.html) that will work out which lipids belong to the top leaflet and which belong to the lower leaflet. This works on highly curved bilayers and is a surprisingly difficult problem.
- Automate the [generation of HOLE profiles](https://pythonhosted.org/MDAnalysis/documentation_pages/analysis/hole.html). HOLE is the standard code for measuring a pore profile through e.g. an ion channel.

That is it. Hope you found this short tutorial helpful. Please send me any comments, or even clone the repo, make some changes and send me a pull request.