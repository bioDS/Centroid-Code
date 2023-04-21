# Centroid-Code

This repository contains the Python implementation of the Centroid algorithm and a script to generate simulated datasets and recreate plots from the paper.
The simulation part is separated, whereas the Centroid algorithm is part of the python package called `tetres` .

The package has a documentation included as `html{pdf}`, see `tetres/docs/build/html{latex}/index.html{tetres.pdf}`. 
The simulation script contains comments and the functions should be self-explanatory.

# Installation of tetres

Before using the package the c code contained in `/tetres/tetres/trees/` has to be manually compiled by running `make`.
After that install the python package to your Python environment using the `pip install -e '/path/to/repo/tetres/'`.

## Additional dependencies

Some funcitons use R and the following packages need to be installed:
`
beautier,
ape,
phangorn,
phytools
`

BEAST needs to be installed to run the Simulations from `simulationGenerator.py`.

For the error measure plot recreation the following python packages need to be installed
`
Biopython,
scipy,
seaborn,
pandas,
multipledispatch
`

# Running centroid on data

Check out the example folder `/example/` if you want to run the centroid algorithm on your own data

