Stereo matching with Loopy Belief Propagation and Graph Cuts
===========================================================
This is a git repository for my qualifying exam to solve
stereo matchin proglem using two MRF based algorithms
(Loopy Belief Propagation and Graph Cuts).

## How to run

You can run this code using the following commands

``` bash

./stereo_reconstruction algorithm[block,lbp,gt_swap,gt_exp] left.png right.png
``` 

Please see the following example.

``` bash
cd MRF_JW/
make
cd ..
make
./stereo_reconstruction block images/Tsukuba/tsukuba-imL.png images/Tsukuba/tsukuba-imR.png
```

## Experimental Results
You can find experimenal results (disparity maps and log file that has informatin like
running time of algorithms, mean energy, and mean squared errors)
in the "stereo_matching_results directory."

December, 2016
Lee, Jangwon
leejang@indiana.edu

