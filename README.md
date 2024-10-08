# EIC | Simulation of Semi-Inclusive DIS Physics


This program aims to analyze Semi-Inclusive Deep Inelastic Scattering (SIDIS) processes, involving polarized beams of electrons and protons at high energy to simulate the future innovative measure of the Electron-Ion Collider (EIC). The EIC, set to be located at the Brookhaven National Laboratory (BNL) in New York, promises unprecedented accuracy in studying nuclear components, providing crucial insights into the fundamental quark and gluon constituents of the nuclear structure.

The data are taken from:  root://dtn-ic.jlab.org//work/eic2/EPIC/RECO/24.06.0/epic_craterlake/DIS/NC/18x275/minQ2=1/pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_4.1388.eicrecon.tree.edm4eic.root ./ 
(here the 24.06 is the campaign and the file name is just an example).

____________________________________
FILE NAME:
- "pythia_ep_noradcor_18x275_eic.txt":   Pythia6.4 event generation of electron-proton collisions, the energy beam for the electron is set at 18 GeV while the proton one at 275 GeV. This program was created by the Software & Computing group in the ePIC Campaign1.0.0 some months ago, so it is not mine, there are just a few corrections but nothing valuable. More about the Event generation and Pythia6 on the main portal of the ePIC Collaboration https://eic.github.io/software/pythia6.html.
  
- "pythia_ep_18x275_q2_0.1_1e7_rid.txt":  This file is a reduced output file created to show the program's behaviour quickly. Contains 3000 events with high energy change ($0.1 < Q^2 < 10^7$ $GeV ^2$), instead of 100K events of the original output file.

- "ep_count.py":  This is the Python program used to analyze the output file and reconstruct some important behaviour of the processes.

- "ePIC_task.py":  Just the same program  structure as before, but, with other calculations (a new program was used to have a more performant setup), here there were performed the first task for the ePIC collaboration in the particles identification through the observation of charged particles (still an ongoing program, it is not finished yet).

- "taskA.C" and "taskB.C": This is the principal program of the analysis. It perform the Recontruction and PID for the particle productions. A is for positive particles and B is for the charged ones.

- "run_all.C": Is necessary to run all the main inside the task program.

- "combined.cpp": Combine all the canvas produced by the different main.

- "efficiencyA.cpp" and "efficiencyB.cpp": Calculate the Reconstruction efficiency via MC ID and PID efficiency for positive particles (A) and charged particles (B), and use the combined canvas from 'combined.cpp'.
  
____________________________________
OUTPUT

"ep_count.py" reconstructs different kinematic variables of the system and its output consists of a set of four Figures which contain different graphs:

- Figure 1: Show the Pion, Kaons and Protons production with the variation of four main variables. $Q^2$ encodes information about the energy exchanged in the process, $x_B$ known as 'Bjorken variables', is the fraction of the nucleon's momentum carried by the struct quark (in the parton model), $z$ is a SIDIS observable and show the longitudinal momentum fraction carried by the identified hadron and $P_{hT}$ is the transversal momentum of the identified hadron.
  
- Figure 2: The first line shows the total particle count and the normalized count per event. The second line focuses on the reconstruction  of the $x_B$ vs $Q^2$ graph for Pions and Kaons (due to the low number of events in the reduced file, they will have a very low sensibility and resolution, especially for the Kaons one). The latter line displays the density of the pion momentum in the transverse plane (x,y) and in the longitudinal plane (z, $p_T).

- Figure 3: Contains the normalized particle identification for the initial and final state, with a focus on the unstable particles generated during the process.

- Figure 4: There is a bar graph used to show how the free quarks and diquarks hadronize and three polar diagrams to show the angular propagation of the Pions, Kaons and scattered Electrons, which will be crucial to achieving the best detector design.


Instead, in "ePIC_task.py", there are three sets of graphs:

- Figure 1: Number of particles / total # of positive particles detected, with particles = Pions+, Kaons+ and Protons, in function of $Q^2$, $x$, $z$, $P_{hT}$, Momentum and rapidity ($\eta$).

- Figure 2: The same as before but with negative particles and only with the reconstruction of Pions- and Kaons-

- Figure 3: Show the angular production of the particles, the first line is for the positive case while the second line is for the negative one. Each graph is normalized to 1 to achieve a better resolution and comprehension.


  "task#.C":

  - Production as a function of $Q^2, x_B, z, P_{hT}, \eta, \phi, P_h$ for pion, kaon and proton.
    
    
________________________________________
