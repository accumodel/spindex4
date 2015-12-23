# spindex4
A quad tree spatial index written in python and cython
Some of the differences between other good spatial indices
available in Python (such as Rtree https://pypi.python.org/pypi/Rtree/)
are that:
 - Written in Python with some speedups by adding some Cython and compiling
 - There is a cool visualization method that uses qt to plot the boxes and points
 - There are some other troubleshooting methods
 - Your object key can be a normal hashable object and is not constrained to be an integer
 
Anyone who wants to help expand the spatial index is welcome.  Some areas of work to be done:
 - Complete "remove" method.
 - Set up a better test suite
 - Add a "nearest" method
 
 Dependencies:
  - Cython if you want to modify the .pyx file and recompile
