This is the repository for the manuscript: "Critical dynamics in the spread of focal epileptic seizures: network connectivity, neural excitability and phase transitions
Authors: S Amin Moosavi, Viktor K Jirsa, Wilson Truccolo


The repository contains:

(1) The compressed (zip) datasets from each of the 5 patient-specific connectivity networks required to reproduced the results in this manuscript; Each dataset contains ascii files specifying 

(a) the name of the brain areas (areas.txt)
(b) non-normalized connectivity weights (weights.txt)
(c) tract lengths of white-matter fibers (tract_lengths.txt)

The additional files centres.txt, average_orientations.txt, cortical.txt, are provided, but were not used in this study.

(2) The Matlab codes for the linear stability analysis, and the Python codes for the simulation of patient-specific Epileptor network models. 

(a) to run the Matlab code for stability analysis please unzip the connectivity files for each patient located in /data/ directory.

(b) run /code/main.m for stability analysis of the epileptor model.

(c) run /code/main.ipynb for simulation of the epileptor model.



