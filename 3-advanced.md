# 3. More Advanced Techniques

We won't go through any problems in this section but will instead just look at some of the differences between [VMD](http://www.ks.uiuc.edu/Research/vmd/) and [MDAnalysis](https://code.google.com/p/mdanalysis/).

## 3.1 VMD

###Advantages

- you can see what is happening!
- many of the `measure` commands have been written in `C` so are very fast. Some have even been written to use a GPU if one is available and you have the right CUDA drivers installed.
- there are a large number of [examples](http://www.ks.uiuc.edu/Research/vmd/script_library/) 
- [VMD](http://www.ks.uiuc.edu/Research/vmd/) is available for Windows, Linux and Mac.


###Disadvantages

- If you are familiar with `C`, `perl`, `python` or other common languages then `Tcl` can be a bit hard to learn. There aren't the books either or number of online resources.
- The `Tcl` implementation is part of the [VMD](http://www.ks.uiuc.edu/Research/vmd/) application so you are restricted to the functionality that is included. There is some math functionality but there is no linear algebra or anything more sophisticated. Also you can't use other modules.
- Because you launch [VMD](http://www.ks.uiuc.edu/Research/vmd/) to run any `Tcl` you write, it is difficult to automated any approach that uses [VMD](http://www.ks.uiuc.edu/Research/vmd/). For example, you can't add command line flags but typically have to generate the `Tcl` code each time.
- Difficult to ""write good code" since there is no testing framework.

### Things it can do that others can't

- Render images
- Measure the solvent accessible surface area

## 3.2 python/MDAnalysis

###Advantages

- Because [MDAnalysis](https://code.google.com/p/mdanalysis/) is a python module, it is independent of whichever python implementation you are using. So you can use regular `python` or perhaps `ipython` or even an [`ipython notebook`](http://ipython.org/notebook.html). 
- This also means that once you have the coordinates loaded in and available through the Universe object, you are free to use any other python module to help in your analysis. The most obvious ones are numpy and scipy, but there are many that could be useful, e.g. [networkx](https://networkx.github.io) or even [image processing techniques](http://pubs.rsc.org/en/Content/ArticleLanding/2014/FD/c3fd00131h).
- Because you can add command line flag functionality to your program it is easy to integrate it into an automated workflow.
- python makes writing good code since you can 
	- add unit tests (e.g. [nose](https://nose.readthedocs.org/en/latest/))
	- automatically create documentation (e.g. [python docstrings](http://www.pythonforbeginners.com/basics/python-docstrings))

	
###Disadvantages

- Although many of the numpy methods are optimised in `C`, it can be slower than [VMD](http://www.ks.uiuc.edu/Research/vmd/).
- It is newer and therefore less developed (at present, at least)
- There is no way of 'seeing' what is happening so you can end up using it ini conjunction with [VMD](http://www.ks.uiuc.edu/Research/vmd/).

### Things it can do that others can't

- [MDAnalysis](https://code.google.com/p/mdanalysis/) has a neat function called [LeafletFinder](https://mdanalysis.googlecode.com/svn/trunk/doc/html/documentation_pages/analysis/leaflet.html) that will work out which lipids belong to the top leaflet and which belong to the lower leaflet. This works on highly curved bilayers and is a surprisingly difficult problem.
- Automate the [generation of HOLE profiles](https://pythonhosted.org/MDAnalysis/documentation_pages/analysis/hole.html). HOLE is the standard code for measuring a pore profile through e.g. an ion channel.