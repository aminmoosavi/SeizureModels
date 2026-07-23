Code repository 
Paper: Criticality in probabilistic models of spreading dynamics in brain networks: epileptic seizures
Authors: S Amin Moosavi, Wilson Truccolo
Journal: PLoS Computational Biology (2023)

This part of the repository contains:

(1) Patient specific connectivity data available under the './files/' directory. These datasets have been used and publicly shared before in Moosavi, Jirsa, Truccolo (PLoS ONE, 2022). For each patient-specific network the following is available:

(a) the name of the brain areas (areas.txt)
(b) non-normalized connectivity weights (weights.txt)
(c) tract lengths of white-matter fibers (tract_lengths.txt)

The additional files centres.txt, average_orientations.txt, cortical.txt, are provided, but were not used in this study.

(2) C++ Code for simulation of the probabilistic model via the temporal Gillespie algorithm. 

To run a simulation of the probabilistic model with Gillespie algorithm please follow the steps below.

(a) To clean old compiled files, run the following line in terminal 

make clean

(b) To compile the code, run the following line in terminal

make

(c) To run an example of a simulation with network of subject P1 run the following 

./main1 P1 1 0.45 -0.112 0.0026 20

The 6 inputs to main1 file are clarified below

./main1 <network> <EZ> <w> <E> <E_ez> <R>

The parameters being passed are:

<network>: the patient-specific network. The patient-specific networks are provided in this repository under "./files/" folder. To run simulations of a patient-specific network not provided here, generate three txt files: (i) the interaction weights with N*N entries in N lines, (ii) the tract-lengths with N*N entries in N lines, and (iii) the EZ nodes. Use the same labeling and format as the files given in this repository.  

<EZ>: the index for the EZ node considered to be active in the simulation. Note that some patient-specific networks have multiple EZ nodes. For example in case of P1, there are two EZ nodes: EZ_1 is node 61 and EZ_2 is node 64.

<w>: global coupling connectivity strength.

<E>: excitability of susceptible nodes in the surround; <E> < 0

<Eez>: excitability of a chosen EZ node; <Eez> > 0

<R>: number of stochastic realizations.

We also note that one can adapt the code in main.cpp to sweep a range of w and E parameters. For that set dw and dE as the intervals between two points on the w and E axes, respectively. In addition set Nw and NE for the number of points on the corresponding axes. The simulations will sweep over the Nw by NE grid in the (w,E) space. w = [ <w> : dw : <w>+Nw*dw ] and E = [ <E> : dE : <E>+NE*dE ]

Results of simulations are be saved in .txt files in the directory ./data/<network>/. For each (w,E), the results for all R realizations are saved in three files: (i) the onset times (in seconds); (ii) offset times (in seconds), and (iii) spread size. The file format and labeling:

onsets_w-<w>_E<E>.txt
offsets_w-<w>_E<E>.txt
spread_w-<w>_E<E>.txt

For example for the simulation "./main1 P1 1 0.45 -0.112 0.0026 20", the following files are generated:

onsets_w-0.450000_E-0.112000.txt
offsets_w-0.450000_E-0.112000.txt
spread_w-0.450000_E-0.112000.txt
