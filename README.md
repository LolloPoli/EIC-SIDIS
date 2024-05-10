# EIC | Simulation of Semi-Inclusive DIS Physics


This program aims to analyze Semi-Inclusive Deep Inelastic Scattering (SIDIS) processes, involving polarized beams of electrons and protons at high energy to simulate the future innovative measure of the Electron-Ion Collider (EIC). The EIC, set to be located at the Brookhaven National Laboratory (BNL) in New York, promises unprecedented accuracy in studying nuclear components, providing crucial insights into the fundamental quark and gluon constituents of the nuclear structure.

____________________________________
FILE NAME:
- "pythia_ep_noradcor_18x275_eic.txt":   Pythia6.4 event generation of electron-proton collisions, the energy beam for the electron is set at 18 GeV while the proton one at 275 GeV. This program was created by the Software & Computing group in the ePIC Campaign1.0.0 some months ago, so it IS NOT MINE, there are just a few corrections but nothing valuable. More about the Event generation and Pythia6 on the main portal of the ePIC Collaboration https://eic.github.io/software/pythia6.html.
  
- "pythia_ep_18x275_q2_0.1_1e7_rid.txt":  This file is a reduced output file created to show the program's behaviour quickly. Contains 3000 events with high energy change ($0.1 < Q^2 < 10^7$ $GeV ^2$), instead of 100K events of the original output file.

- "ep_count.py":  This is the Python program used to analyze the output file and reconstruct some important behaviour of the processes.

- "ePIC_task.py:  Just the same program as before, with other calculations (a new program was used to have a more performant setup), here there were performed the first task for the ePIC collaboration in the particles identification through the observation of charged particles.
____________________________________
OUTPUT

"ep_count.py" reconstructs different kinematic variables of the system and its output consists of a set of 4 Figures which contain different graphs:
- Figure 1: Show the Pion, Kaons and Protons production with the variation of four main variables. $Q^2$ encodes information about the energy exchanged in the process, $x_B$ known as 'Bjorken variables', is the fraction of the nucleon's momentum carried by the struct quark (in the parton model), $z$ is a SIDIS observable and show the longitudinal momentum fraction carried by the identified hadron and $P_{hT}$ is the transversal momentum of the identified hadron.
  
- Figure 2: The first line shows the total particle count and the normalized count per event. The second line focuses on the reconstruction  of the $x_B$ vs $Q^2$ graph for Pions and Kaons (due to the low number of events in the reduced file, they will have a very low sensibility and resolution, especially for the Kaons one). The latter line displays the density of the pion momentum in the transverse plane (x,y) and in the longitudinal plane (z, $p_T).

- Figure 3: Contains the normalized particle identification for the initial and final state, with a focus on the unstable particles generated during the process.

- Figure 4: There is a bar graph used to show how the free quarks and diquarks hadronize and three polar diagrams to show the angular propagation of the Pions, Kaons and scattered Electrons, which will be crucial to achieving the best detector design.


Instead, in "ePIC_task.py" three graphs will be:

- Figure 1: Number of particles / total # of positive particles detected, with particles = Pions+, Kaons+ and Protons, in function of $Q^2$, $x$, $z$, $P_{hT}$, Momentum and rapidity ($\eta$).

- Figure 2: The same as before but with negative particles and only with the reconstruction of Pions- and Kaons-

- Figure 3: Show the angular production of the particles, the first line is for the positive case while the second line is for the negative one. Each graph is normalized to 1 to achieve a better resolution and comprehension.
________________________________________
