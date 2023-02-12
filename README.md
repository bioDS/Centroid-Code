# Centroid-Code

This repository contains the Python implementation of the Centroid algorithm and a script to generate simulated datasets.
The simulation part is a separate script, whereas the Centroid algorithm is part of a package called `tetres` .

The package has a documentation included as `html` and `pdf`, see `treeeoclock/docs/build/html{latex}/index.html{tetres.pdf}`. 
The simulation script contains comments and the functions should be self-explanatory.

# Usage

Before using the package the c code contained in `/tetres/tetres/trees/` has to be manually compiled by running `make`.
After that install the python package to your Python environment using the `pip install '/path/to/repo/tetres/'`.
Then the `simulationGenerator.py` script can be executed, it is currently setup to run a single simulation on ten taxa.
To change the setup, the funciton call within the main function of the script has to be altered.

