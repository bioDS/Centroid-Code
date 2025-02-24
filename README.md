# Centroid-Code

This repository contains the Python implementation of the [Centroid algorithm](https://doi.org/10.1101/2023.05.08.539790) and a script to generate simulated datasets and recreate plots from the paper.
The simulation part is separated, whereas the Centroid algorithm is part of the python package called `tetres` .

| :information_source: New version available                             |
|:---------------------------------------------------------|
|Check out the more recent version of this package at [https://github.com/bioDS/tetres](https://github.com/bioDS/tetres)|


The package has a documentation included as `html{pdf}`, see `tetres/docs/build/html{latex}/index.html{tetres.pdf}`.<br>
The simulation script contains comments and the functions should be self-explanatory.


# Installation of tetres

Before using the package the c code contained in `/tetres/tetres/trees/` has to be manually compiled by running `make`.
After that install the python package to your Python environment using the `pip install -e '/path/to/repo/tetres/'`.

## Additional dependencies

Some functions use R and the following packages should be installed:
`
beautier,
ape,
phangorn,
phytools
`

BEAST2.6 (this is because of the R package beautier) needs to be installed to run the Simulations from `simulationGenerator.py`.
You will also need to set the path to BEAST2 and to treeannotator in the script `simulationGenerator.py` line 25 and 27.

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


# Reference

If you publish a paper using this software, please cite<br>
Lars Berling, Lena Collienne, and Alex Gavryushkin<br>
**Estimating the mean in the space of ranked phylogenetic trees**<br>
*Bioinformatics* 2024<br>
[https://doi.org/10.1093/bioinformatics/btae514](https://doi.org/10.1093/bioinformatics/btae514) 
