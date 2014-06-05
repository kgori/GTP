GTP
===

A fork of Megan Owen and Scott Provan's [__java implementation__](http://comet.lehman.cuny.edu/owen/code.html) of their **GTP algorithm** for calculating the geodesic distance between phylogenetic trees, from [__this paper__](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5396323). The original code only allowed for all-all comparisons, this version is modified to allow row-wise comparisons - i.e. the first tree in a file against all other trees. This is activated with the '-r' switch.

Installation
------------

First compile the sources:

    javac polyAlg/*.java distanceAlg1/*.java
    
Then pack the jar file:

    jar cfm gtp.jar manifest.txt polyAlg/*.class distanceAlg1/*.class

Or run the Makefile by typing:

    make

You can clean up the compiled classes with:

    make clean

Usage
-----

Run with ```java -jar gtp.jar [options] trees_file```.

--help outputs this:

    Command line syntax:
    gtp [options] treefile
    Optional arguments:
    	 -d 	 double check results, by computing each distance with the target tree as the starting tree and vice versa; default is false
    	 -h || --help 	 displays this message
    	 -n 	 normalize (vector of the lengths of all edges has length 1)
    	 -o <outfile> 	 store the output in the file <outfile>
    	 -u 	 unrooted trees (default is rooted trees)
    	 -v || --verbose 	 verbose output
    	 -r || --row 	 Compute geodesics between the first tree in the input file and all the others (not an all-all comparison)
