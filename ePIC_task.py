import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
from scipy.interpolate import make_interp_spline 
from scipy.interpolate import griddata

#_____________________________________________________________________________________________________________________________________________________________________________________

# COPY PASTE OF THE PREVIOUS PROGRAM, TAKING ALL THE DATA SINCE COULD BE USEFUL

# DEFINITION OF A FUNCTION ABLE TO READ MY PYTHIA6.4 OUTPUT AND COLLECT IMPORTANT DATA INTO AN ARRAY (OR MORE)
# THE COUNT OF THE TOTAL NUMBER OF EVENTS IS NECESSARY TO NORMALIZE THE DATA AND HAVE AN ESTIMATE BY EVENT OF THE MEASUREMENTS
def read_pythia_output(filename):
    particles = []
    particles_scatter = []
    event_count = 0  
    last_s_el_px, last_s_el_py, last_s_el_pz, last_s_el_E = 0, 0, 0, 0
    last_s_ph_pz, last_s_ph_E = 0, 0
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() != '' and line[0] != '!':  #  TO SKIP USELESS LINE LIKE COMMENT OR EMPTY LINES
                parts = line.split()
                if parts[0] == '0':                    # IN THE OUTPUT FILE THE LINE OF THE EVENT DESCRIPTION START WITH A '0', SO WILL BE EASY TO COUNT THE NUMBER OF EVENTS
                    event_count += 1
                if len(parts) >= 11 and parts[0] not in ('1', '2'):  # ENSURE THAT THE LINES ARE LONG ENOUGH, THE FIRST TWO LINES ARE REMOVED SINCE REFERS TO THE BEAM PARTICLES
                    try:
                        pdg_code = int(parts[2])      # THE PDG CODE IS RELATED TO THE IDENTITY OF THE PARTICLE 
                        status = int(parts[1])        # 'K(I,1)' OR 'KS' IS THE STATUS AND TELL ME INFORMATION ABOUT THE PARTICLES CONDITIONS
                        # FOR EXAMPLE, IF KS < 10 THE PARTICLES DO NOT DECAY/FRAGMENT INTO OTHER PARTICLES FOR ALL THE EVENT TIME MEASUREMENT
                        # KS > 10 DESCRIBED PARTICLES THAT COULD UNDERGO INTO DIFFERENT TYPE OF DECAYS OR FRAGMENTATION, HENCE, WILL NOT BE PRESENT IN THE FINAL STATE
                        origin = int (parts[3])       # INFORMATION ABOUT THE IDENTITY OF THE PARENT PARTICLE
                        daughter1 = int(parts[4])     # FIRST DAUGHTER
                        daughter2 = int(parts[5])     # LAST DAUGHTER
                        px = float(parts[6])          # MOMENTUM CALCULATION ALONG THE THREE AXIS (x,y,z) [GeV/c]
                        py = float(parts[7])
                        pz = float(parts[8])
                        E = float(parts[9])           # ENERGY OF THE PARTICLE [GeV/c]
                        M = float(parts[10])          # MASS IN [GeV/c^2]
                        phi = float(parts[15]) if len(parts) >= 16 else last_phi  # THE PHI IS COLLECTED FROM THE EVENT LINE, SO IS IMPORTANT TO MAINTAIN THE SAME VALUE FOR ALL THE EVENT FOR ALL PARTICLES CALCULATIONS 
                        last_phi = phi                # UPDATE OF THE LAST CHECKED VALUE
                        real_q2 = float(parts[11]) if len(parts) >= 16 else last_real_q2 
                        last_real_q2 = real_q2
                        event_number = int(parts[0])  # THAT WILL HELP TO RECALL THE EVENT 3,4,5 WHICH CARRY INFORMATION ABOUT THE e-,gamma,p+ SCATTERED PARTICLES
                        # SAME PHI AS BEFORE BUT WITH LOWER MULTIPLICITY TO PERFORM A GRAPH
                        s_phi = float(parts[15]) if event_number == 0 else None
                        # ALSO THE ENERGY OF THE SCATTERED ELECTRON IS CALCULATED TWO TIMES, HERE IS NECESSARY FOR THE GRAPH
                        s_E_electron = float(parts[9]) if event_number == 3 else None
                        # COLLECTION OF THE SCATTERED DATA, AND MULTIPLICATION TO USE THE SAME VALUE FOR ALL THE EVENT
                        # SCATTERED ELECTRON
                        if event_number == 3:         # THE THIRD LINE OF EACH EVENT REPRESENT THE SCATTERED ELECTRON, SO I COLLECTED THE DATA
                            s_el_px = px
                            s_el_py = py
                            s_el_pz = pz
                            s_el_E = E
                        else:                         # I REPRODUCE THE SAME VALUE FOR ALL THE EVENT UNTIL THE NEW EVENT START
                            s_el_px = last_s_el_px
                            s_el_py = last_s_el_py
                            s_el_pz = last_s_el_pz
                            s_el_E = last_s_el_E
                        # SCATTERED PION              # SAME STRUCTURE AS BEFORE, THE FOURTH LINE IS ALWAYS THE SCATTERED PHOTON
                        if event_number == 4:  
                            s_ph_pz = pz
                            s_ph_E = E
                        else:
                            s_ph_pz = last_s_ph_pz
                            s_ph_E = last_s_ph_E
                        # UPDATE OF THE VALUES   
                        last_s_el_px = s_el_px
                        last_s_el_py = s_el_py
                        last_s_el_pz = s_el_pz
                        last_s_el_E = s_el_E
                        last_s_ph_pz = s_ph_pz
                        last_s_ph_E = s_ph_E                      
                        # COLLECTION OF THE DATA
                        particles.append((pdg_code, status, px, py, pz, E, M, phi, real_q2, event_number, origin, daughter1, daughter2))
                        particles_scatter.append((s_el_px, s_el_py, s_el_pz, s_el_E, s_ph_pz, s_ph_E, s_phi, s_E_electron))
                                           
                    except (ValueError, IndexError):
                        pass                          # LINE SKIP
    return particles, particles_scatter, event_count

# DEFINITION OF A FUNCTION TO IDENTIFY THE PARTICLES TYPE WITH THE PDG CODE
def identify_particle(pdg_code):
    particle_dict = {
        11: "Electron",
        -11: "Positron",
        211: "Pion+",
        -211: "Pion-",
        111: "Pion0",
        321: "Kaon+",
        -321: "Kaon-",
        311: "Kaon0",
        -311: "Kaon0bar",
        2212: "Proton+",
        -2212: "Proton-",
        2112: "Neutron",
        -2112: "Neutronbar",
        21: "Gluon",
        22: "Photon",
        2: "Up",
        -2: "Upbar",
        1: "Down",
        -1: "Downbar",
        3: "Strange",
        -3: "Strangebar",
        4: "Charm",
        -4: "Charmbar",
        213: "Rho+",
        -213: "Rho-",
        113: "Rho0",
        2101: "ud0",                 # DIQUARK (ud)_0
        2103: "ud1",                 # DIQUARK (ud)_1
        2203: "uu",                  # DIQUARK (uu)_1
        2224: "Delta++",
        -2224: "Delta--bar",
        2214: "Delta+",
        -2214: "Delta+bar",
        2114: "Delta0",
        -2114: "Delta0bar",
        1114: "Delta-",
        -1114: "Delta-bar",
        3122: "Lambda",
        -3122: "Lambdabar",
        4122: "CLambda+",
        -4122: "CLambda+bar",
        3222: "Epsilon+",
        -3222: "Epsilon+bar",
        3212: "Epsilon0",
        -3212: "Epsilon0bar",
        3112: "Epsilon-",
        -3112: "Epsilon-bar",
        3224: "Epsilon*+",
        -3224: "Epsilon*+bar",
        3214: "Epsilon*0",
        -3214: "Epsilon*0bar",
        3114: "Epsilon*-",
        -3114: "Epsilon*-bar",
        221: "Eta",
        223: "Omega",
        333: "Phi",
        91: "Cluster",             # PARTON SYSTEM IN CLUSTER FRAGMENTATION
        92: "String",              # PARTON SYSTEM IN STRING FRAGMENTATION
        9900110: "Dif",            # DIFRACTIVE pi0/rho0/photon STATE
        9900210: "pi_dif",         #     //     PION     //
        9900220: "w_dif",          #     //     OMEGA    //
        9900330: "phi_dif",        #     //     PHI      //
        9900440: "jpsi_dif",       #     //     J/PSI    //
        9902110: "n_dif",          #     //     NEUTRON  //
        9902210: "p_dif",          #     //     PROTON   //
        130: "K0L",                # SEQUENCE OF STRANGE KAONS
        -130: "K0Lbar",
        310: "K0S",
        -310: "K0Sbar",
        313: "K*0",
        -313: "K*0bar",
        323: "K*+",
        -323: "K*+bar",
        3322: "Xi0",
        -3322: "Xi0bar",
        3312: "Xi-",
        -3312: "Xi-bar",
        3324: "Xi*0",
        -3324: "Xi*0bar",
        3314: "Xi*-",
        -3314: "Xi*-bar",
    }
    return particle_dict.get(pdg_code, "Unknown")

# READ OF THE PYTHIA OUTPUT FILE
filename = "pythia_ep_18x275_q2_0.1_1e7_rid.txt"                           # SMALL PART OF THE DATA, USED TO OBTAIN A FASTER SIMULATION 
#filename = "pythia_ep_noradcor_18x275_q2_0.1_10000000_run1_100K.txt"      # ORIGINAL FILE 
particles, particles_scatter, event_count = read_pythia_output(filename)

# COUNTERS FOR ALL THE CONSIDERED PARTICLE TYPES
electron_count_pos, electron_count_neg = 0, 0
pion_count_pos, pion_count_neg = 0, 0
kaon_count_pos, kaon_count_neg = 0, 0
proton_count_pos, proton_count_neg = 0, 0
neutron_count_pos, neutron_count_neg = 0, 0
rho_count_pos, rho_count_neg= 0, 0
gluon_count_pos, gluon_count_neg = 0, 0
photon_count_pos, photon_count_neg = 0, 0
other_count_pos, other_count_neg = 0, 0
up_count_pos, up_count_neg = 0, 0
down_count_pos, down_count_neg = 0, 0
charm_count_pos, charm_count_neg = 0, 0
strange_count_pos, strange_count_neg = 0, 0
diq_count_pos, diq_count_neg = 0, 0
diq2_count_pos, diq2_count_neg = 0, 0
delta_count_pos, delta_count_neg = 0, 0
lambda_count_pos, lambda_count_neg = 0, 0
epsilon_count_pos, epsilon_count_neg = 0, 0
eta_count_pos, eta_count_neg = 0, 0
omegaphi_count_pos, omegaphi_count_neg = 0, 0
frag_count_pos, frag_count_neg = 0, 0
difractive_count_pos, difractive_count_neg = 0, 0
strangeK_count_pos, strangeK_count_neg = 0, 0
xi_count_pos, xi_count_neg = 0, 0
fragm_up, fragm_down, fragm_uu, fragm_ud = 0, 0, 0, 0


# WE NEED A COUNT FOR ALL THE POSITIVE PARTICLES IN THE FINAL STATE
# ONLY THE POSITIVE PARTICLES WITHOUT DAUGHTERS ARE CONSIDERED
# POSITIVE CASE
total_count_pos = 0
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)       
    event_number = particle[9]  
    dau1 = particle[11]
    dau2 = particle[12]
    if event_number not in (0, 1, 2, 3, 4, 5) and (dau1 == 0 and dau2 == 0):  # THE 0-5 EVENT ARE REMOVED SINCE ARE RELATED TO THE BEAM'S PARTICLES
        if particle_name == "Positron":  
            electron_count_pos+= 1     
            total_count_pos += 1                             
        elif particle_name == "Pion+":
            pion_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Kaon+" or particle_name == "K*+":
            kaon_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Proton+":
            proton_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Up":
            up_count_pos+= 1
            total_count_pos += 1
        elif particle_name =="Downbar":
            down_count_pos+= 1
            total_count_pos += 1
        if particle_name == "Rho+":
            rho_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Strangebar":
            strange_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Charm":
            charm_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "uu":
            diq_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Delta++" or particle_name == "Delta+" or particle_name == "Delta-bar":
            delta_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "CLambda+":
            lambda_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Epsilon+" or particle_name == "Epsilon*+" or particle_name == "Epsilon-bar" or particle_name == "Epsilon*-bar":
            epsilon_count_pos+= 1
            total_count_pos += 1
        elif particle_name == "Xi-" or particle_name == "Xi*-":
            xi_count_pos+= 1
            total_count_pos += 1

# NEGATIVE CASE
# SAME AS BEFORE BUT FOR THE NEGATIVE ONE
total_count_neg = 0
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)       
    event_number = particle[9]  
    dau1 = particle[11]
    dau2 = particle[12]
    if event_number not in (0, 1, 2, 3, 4, 5) and (dau1 == 0 and dau2 == 0): 
        if particle_name == "Electron":  
            electron_count_neg+= 1     
            total_count_neg += 1                             
        elif particle_name == "Pion-":
            pion_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Kaon-" or particle_name == "K*+bar":
            kaon_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Proton-":
            proton_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Upbar":
            up_count_neg+= 1
            total_count_neg += 1
        elif particle_name =="Down":
            down_count_neg+= 1
            total_count_neg += 1
        if particle_name == "Rho-":
            rho_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Strange":
            strange_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Charmbar":
            charm_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Delta++bar" or particle_name == "Delta-" or particle_name == "Delta+bar":
            delta_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "CLambda+bar":
            lambda_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Epsilon-" or particle_name == "Epsilon*+bar" or particle_name == "Epsilon-" or particle_name == "Epsilon*-":
            epsilon_count_neg+= 1
            total_count_neg += 1
        elif particle_name == "Xi-bar" or particle_name == "Xi*-bar":
            xi_count_neg+= 1
            total_count_neg += 1


# NORMALIZATION OF THE COUNTS
def norm_particle_count(particle_c, event_c):
    norm_count = particle_c/event_c
    return norm_count
# POSITIVE COUNT 
norm_pion_pos = norm_particle_count(pion_count_pos, total_count_pos)
norm_kaon_pos = norm_particle_count(kaon_count_pos, total_count_pos)
norm_proton_pos = norm_particle_count(proton_count_pos, total_count_pos)
# NEGATIVE
norm_pion_neg = norm_particle_count(pion_count_neg, total_count_neg)
norm_kaon_neg = norm_particle_count(kaon_count_neg, total_count_neg)

# FUNCTION TO HAVE LINE HISTOGRAMS
def plot_histogram_lines(data, bins, color, label, total_events, linewidth=2):
    hist, bins_edges = np.histogram(data, bins=bins)
    bins_centers = (bins_edges[:-1] + bins_edges[1:]) / 2
    norm_hist = hist / total_events 
    plt.plot(bins_centers, norm_hist, color=color, label=label, linewidth=linewidth)


# DEFINITION OF MY VARIABLES | WITHOUT DISTINCTION OF THE CHARGE
# ENERGY
pion_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"] and particle[11] == 0 and particle[12] == 0]
proton_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
kaon_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"] and particle[11] == 0 and particle[12] == 0]
pion_mass = [particle[6] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"] and particle[11] == 0 and particle[12] == 0]
kaon_mass = [particle[6] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"] and particle[11] == 0 and particle[12] == 0]
proton_mass = [particle[6] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
# PION
pion_x = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"] and particle[11] == 0 and particle[12] == 0]
pion_y = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"] and particle[11] == 0 and particle[12] == 0]
pion_z = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"] and particle[11] == 0 and particle[12] == 0]
# PROTON
proton_x = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
proton_y = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
proton_z = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
# KAON
kaon_x = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"] and particle[11] == 0 and particle[12] == 0]
kaon_y = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"] and particle[11] == 0 and particle[12] == 0]
kaon_z = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"] and particle[11] == 0 and particle[12] == 0]

#_________________________________________________ POSITIVE PARTICLES ________________________________________________________________________________________________________________________________

# PION
pion_x_pos = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Pion+"] and particle[11] == 0 and particle[12] == 0]
pion_y_pos = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Pion+"] and particle[11] == 0 and particle[12] == 0]
pion_z_pos = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Pion+"] and particle[11] == 0 and particle[12] == 0]
pion_E_pos = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Pion+"] and particle[11] == 0 and particle[12] == 0]
# PROTON
proton_x_pos = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Proton+"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
proton_y_pos = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Proton+"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
proton_z_pos = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Proton+"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
proton_E_pos = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Proton+"] and particle[9] != 5 and particle[11] == 0 and particle[12] == 0]
# KAON
kaon_x_pos = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "K*+"] and particle[11] == 0 and particle[12] == 0]
kaon_y_pos = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "K*+"] and particle[11] == 0 and particle[12] == 0]
kaon_z_pos = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "K*+"] and particle[11] == 0 and particle[12] == 0]
kaon_E_pos = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "K*+"] and particle[11] == 0 and particle[12] == 0]

#_________________________________________________ NEGATIVE PARTICLES _______________________________________________________________________________________________________________________________

# PION
pion_x_neg = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Pion-"] and particle[11] == 0 and particle[12] == 0]
pion_y_neg = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Pion-"] and particle[11] == 0 and particle[12] == 0]
pion_z_neg = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Pion-"] and particle[11] == 0 and particle[12] == 0]
pion_E_neg = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Pion-"] and particle[11] == 0 and particle[12] == 0]
# KAON
kaon_x_neg = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Kaon-", "K*+bar"] and particle[11] == 0 and particle[12] == 0]
kaon_y_neg = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Kaon-", "K*+bar"] and particle[11] == 0 and particle[12] == 0]
kaon_z_neg = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Kaon-", "K*+bar"] and particle[11] == 0 and particle[12] == 0]
kaon_E_neg = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Kaon-", "K*+bar"] and particle[11] == 0 and particle[12] == 0]

#____________________________________________________________________________________________________________________________________________________________________________________________________

# ELECTRON
electron_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Electron", "Positron"]] # ALL THE ELECTRON OF THE PROCESS
electron_px =  [particle[2] for particle in particles if identify_particle(particle[0]) in ["Electron", "Positron"]]
electron_py =  [particle[3] for particle in particles if identify_particle(particle[0]) in ["Electron", "Positron"]]
electron_pz =  [particle[4] for particle in particles if identify_particle(particle[0]) in ["Electron", "Positron"]]
lepton_angles = [particle[7] for particle in particles]                  # SCATTERED ANGLE WITH HIGH MULTIPLICITY
Q2 = [particle[8] for particle in particles]
# SCATTERED PARTICLES
sc_electron_px = [particle[2] for particle in particles if particle[9] == 3]
sc_electron_py = [particle[3] for particle in particles if particle[9] == 3]
sc_electron_pz = [particle[4] for particle in particles if particle[9] == 3]
sc_electron_E = [particle[5] for particle in particles if particle[9] == 3]
sc_photon_pz = [particle[4] for particle in particles if particle[9] == 4]
sc_photon_E = [particle[5] for particle in particles if particle[9] == 4]
# ENHANCED COUNTS
scattered_el_E = [particle[3] for particle in particles_scatter[1:]]        # THE FIRST TERM IS REMOVED SINCE IS A '0' GENERATED FROM THE 'if' CYCLE 
scattered_el_px = [particle[0] for particle in particles_scatter[1:]]
scattered_el_py = [particle[1] for particle in particles_scatter[1:]]
scattered_el_pz = [particle[2] for particle in particles_scatter[1:]]
scattered_gamma_E = [particle[5] for particle in particles_scatter[2:]]     # THE FIRST TWO TERMS ARE REMOVED DUE TO A DOUBLE NULL GENERATION OF DATA FROM THE CYCLE
scattered_gamma_pz = [particle[4] for particle in particles_scatter[2:]]
sc_el_ang = [particle[6] for particle in particles_scatter]            
scattered_el_ang = [value for value in sc_el_ang if value is not None]      # FILTERED VALUE OF THE SCATTERED ANGLE, THE 'None' ARE REMOVED FROM THE ARRAY
scattered_el_ang_deg = [np.degrees(a) for a in scattered_el_ang]
sc_el_E = [particle[7] for particle in particles_scatter]
scat_el_E = [value for value in sc_el_E if value is not None]               # ENERGY NEEDED FOR THE GRAPH 
mod_q = [np.sqrt(E**2 + pz**2) for E, pz in zip(scattered_gamma_E, scattered_gamma_pz)]   

# KINEMATIC VARIABLES CONSTRUCTION
#_____________________________________________________________________________ POSITIVE CASE _______________________________________________________________________________________________________

# PION
pion_transv_mom_pos = [np.sqrt(px**2 + py**2) for px, py in zip (pion_x_pos, pion_y_pos)]  
pion_mom_pos = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(pion_x_pos, pion_y_pos, pion_z_pos)]
pion_long_mom_pos = [np.sum(E - pz) for E, pz in zip(pion_E_pos, pion_z_pos)]
pion_angles_pos = [np.arccos(pz/p) for pz, p in zip(pion_z_pos, pion_mom_pos)]            # Z-AXIS
pion_angles_pos_deg = [np.degrees(a) for a in pion_angles_pos]
pion_rapidity_pos = [-np.log(np.tan(t/2)) for t in pion_angles_pos]
z_pion_pos = [(Eph - np.abs(Phz))/(Ee + np.abs(Pze)) for Eph, Phz, Ee, Pze in zip(pion_E_pos, pion_z_pos, scattered_gamma_E, scattered_gamma_pz)]
PhT_pion_pos = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(pion_z_pos, scattered_el_pz, mod_q)]
# PROTON
proton_transv_mom_pos = [np.sqrt(px**2 + py**2) for px, py in zip (proton_x_pos, proton_y_pos)]  
proton_mom_pos = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(proton_x_pos, proton_y_pos, proton_z_pos)]
proton_long_mom_pos = [np.sum(E - pz) for E, pz in zip(proton_E_pos, proton_z_pos)]
proton_angles_pos = [np.arccos(pz/p) for pz, p in zip(proton_z_pos, proton_mom_pos)]      # Z-AXIS
proton_angles_pos_deg = [np.degrees(a) for a in proton_angles_pos]
proton_rapidity_pos = [-np.log(np.tan(t/2)) for t in proton_angles_pos]
z_proton_pos = [(Eph - np.abs(Phz))/(Ee + np.abs(Pze)) for Eph, Phz, Ee, Pze in zip(proton_E_pos, proton_z_pos, scattered_gamma_E, scattered_gamma_pz)]
PhT_proton_pos = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(proton_z_pos, scattered_el_pz, mod_q)]
# KAON
kaon_transv_mom_pos = [np.sqrt(px**2 + py**2) for px, py in zip (kaon_x_pos, kaon_y_pos)]  
kaon_mom_pos = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(kaon_x_pos, kaon_y_pos, kaon_z_pos)]
kaon_long_mom_pos = [np.sum(E - pz) for E, pz in zip(kaon_E_pos, kaon_z_pos)]
kaon_angles_pos = [np.arccos(pz/p) for pz, p in zip(kaon_z_pos, kaon_mom_pos)]            # Z-AXIS
kaon_angles_pos_deg = [np.degrees(a) for a in kaon_angles_pos]
kaon_rapidity_pos = [-np.log(np.tan(t/2)) for t in kaon_angles_pos]
z_kaon_pos = [(Eph - np.abs(Phz))/(Ee + np.abs(Pze)) for Eph, Phz, Ee, Pze in zip(kaon_E_pos, kaon_z_pos, scattered_gamma_E, scattered_gamma_pz)]
PhT_kaon_pos = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(kaon_z_pos, scattered_el_pz, mod_q)]

# DOUBLE ANGLE METHOD | USED FOR THE CONSTRUCTION OF THE BJORKEN VARIABLE 'x' AND THE INELASTICITY 'y'
# PION
y_DA_pi_pos = [np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) for p, t in zip(pion_angles_pos, lepton_angles) if 0.001 <= np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) <= 0.95]
Q2_DA_pi_pos = [np.abs(4*18*18*(1-abs(y))/(np.tan(t/2)**2)) for y, t in zip(y_DA_pi_pos, lepton_angles)]
x_DA_pi_pos = [Q/(4*18*275*y*275) for Q, y in zip(Q2_DA_pi_pos, y_DA_pi_pos) if Q/(4*18*275*y*275)]
# PROTON
y_DA_proton_pos = [np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) for p, t in zip(proton_angles_pos, lepton_angles) if 0.001 <= np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) <= 0.95]
Q2_DA_proton_pos = [np.abs(4*18*18*(1-abs(y))/(np.tan(t/2)**2)) for y, t in zip(y_DA_proton_pos, lepton_angles)]
x_DA_proton_pos = [Q/(4*18*275*y*275) for Q, y in zip(Q2_DA_proton_pos, y_DA_proton_pos)]
# KAON
y_DA_kaon_pos = [np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) for p, t in zip(kaon_angles_pos, lepton_angles) if 0.001 <= np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) <= 0.95]
Q2_DA_kaon_pos = [np.abs(4*18*18*(1-abs(y))/(np.tan(t/2)**2)) for y, t in zip(y_DA_kaon_pos, lepton_angles)]
x_DA_kaon_pos = [Q/(4*18*275*y*275) for Q, y in zip(Q2_DA_kaon_pos, y_DA_kaon_pos)]

#_______________________________________________________________________ NEGATIVE CASE ______________________________________________________________________________________________________________

# PION
pion_transv_mom_neg = [np.sqrt(px**2 + py**2) for px, py in zip (pion_x_neg, pion_y_neg)]  
pion_mom_neg = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(pion_x_neg, pion_y_neg, pion_z_neg)]
pion_long_mom_neg = [np.sum(E - pz) for E, pz in zip(pion_E_neg, pion_z_neg)]
pion_angles_neg = [np.arccos(pz/p) for pz, p in zip(pion_z_neg, pion_mom_neg)]            # Z-AXIS
pion_angles_neg_deg = [np.degrees(a) for a in pion_angles_neg]
pion_rapidity_neg = [-np.log(np.tan(t/2)) for t in pion_angles_neg]
z_pion_neg = [(Eph - np.abs(Phz))/(Ee + np.abs(Pze)) for Eph, Phz, Ee, Pze in zip(pion_E_neg, pion_z_neg, scattered_gamma_E, scattered_gamma_pz)]
PhT_pion_neg = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(pion_z_neg, scattered_el_pz, mod_q)]
# KAON
kaon_transv_mom_neg = [np.sqrt(px**2 + py**2) for px, py in zip (kaon_x_neg, kaon_y_neg)]  
kaon_mom_neg = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(kaon_x_neg, kaon_y_neg, kaon_z_neg)]
kaon_long_mom_neg = [np.sum(E - pz) for E, pz in zip(kaon_E_neg, kaon_z_neg)]
kaon_angles_neg = [np.arccos(pz/p) for pz, p in zip(kaon_z_neg, kaon_mom_neg)]            # Z-AXIS
kaon_angles_neg_deg = [np.degrees(a) for a in kaon_angles_neg]
kaon_rapidity_neg = [-np.log(np.tan(t/2)) for t in kaon_angles_neg]
z_kaon_neg = [(Eph - np.abs(Phz))/(Ee + np.abs(Pze)) for Eph, Phz, Ee, Pze in zip(kaon_E_neg, kaon_z_neg, scattered_gamma_E, scattered_gamma_pz)]
PhT_kaon_neg = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(kaon_z_neg, scattered_el_pz, mod_q)]

# DOUBLE ANGLE METHOD | USED FOR THE CONSTRUCTION OF THE BJORKEN VARIABLE 'x' AND THE INELASTICITY 'y'
# PION
y_DA_pi_neg = [np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) for p, t in zip(pion_angles_neg, lepton_angles) if 0.001 <= np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) <= 0.95]
Q2_DA_pi_neg = [np.abs(4*18*18*(1-abs(y))/(np.tan(t/2)**2)) for y, t in zip(y_DA_pi_neg, lepton_angles)]
x_DA_pi_neg = [Q/(4*18*275*y*275) for Q, y in zip(Q2_DA_pi_neg, y_DA_pi_neg) if Q/(4*18*275*y*275)]
# KAON
y_DA_kaon_neg = [np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) for p, t in zip(kaon_angles_neg, lepton_angles) if 0.001 <= np.abs(np.tan(p/2)/(np.tan(p/2) + np.tan(t/2))) <= 0.95]
Q2_DA_kaon_neg = [np.abs(4*18*18*(1-abs(y))/(np.tan(t/2)**2)) for y, t in zip(y_DA_kaon_neg, lepton_angles)]
x_DA_kaon_neg = [Q/(4*18*275*y*275) for Q, y in zip(Q2_DA_kaon_neg, y_DA_kaon_neg)]


#________________________________________________________________________ POSITIVE GRAPHS ___________________________________________________________________________________________________________

plt.figure(figsize=(12, 8))
# FILTER POSSIBLE INVALID VALUE OF Q2 
# TAKING THE PRODUCTION OF THE POSITIVE PARTICLES vs Q^2
valid_Q2_pi_pos = [q2 for q2 in Q2_DA_pi_pos if q2 > 0]
valid_Q2_kaon_pos = [q2 for q2 in Q2_DA_kaon_pos if q2 > 0]
valid_Q2_proton_pos = [q2 for q2 in Q2_DA_proton_pos if q2 >0]

# IN CASE THERE ARE NO PROBLEMS WITH THE VALUE 
if valid_Q2_pi_pos and valid_Q2_kaon_pos and valid_Q2_proton_pos:
    # GENERATION OF A GREAT FIGURE TO FOUR GRAPHS
    plt.subplot(3,2,1)
    # DEFFINITION OF THE BINS OF Q2
    bins_pi = np.logspace(np.log10(min(valid_Q2_pi_pos)), np.log10(max(valid_Q2_pi_pos)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_Q2_kaon_pos)), np.log10(max(valid_Q2_kaon_pos)), 60)
    bins_proton = np.logspace(np.log10(min(valid_Q2_proton_pos)), np.log10(max(valid_Q2_proton_pos)), 60)
    # HISTOGRAM PLOT | HERE THE FUNCTION TO OBTAIN A LINE HISTOGRAM IS USED
    plot_histogram_lines(Q2_DA_pi_pos, bins=bins_pi, color='skyblue', label=f'Pions ({round(norm_pion_pos,3)})', total_events=total_count_pos)
    plot_histogram_lines(Q2_DA_kaon_pos, bins=bins_kaon, color='lightgreen', label=f'Kaons ({round(norm_kaon_pos,3)})', total_events=total_count_pos)
    plot_histogram_lines(Q2_DA_proton_pos, bins=bins_proton, color='mediumorchid', label=f'Protons ({round(norm_proton_pos,3)})', total_events=total_count_pos)
    plt.xscale('log')
    plt.xlim(1,1e7)
    plt.xlabel('Q^2 [GeV^2]')
    plt.ylabel('Number of Particles')
    plt.title('Number of pi+, P and K+ vs Q^2')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
else:
    print("No valid Q^2 values found")

# x (BJORKEN) VARIABLES OBSERVATION 
valid_X_pi_pos = [x for x in x_DA_pi_pos if x > 0]
valid_X_kaon_pos = [x for x in x_DA_kaon_pos if x > 0]
valid_X_proton_pos = [x for x in x_DA_proton_pos if x >0]

if valid_X_pi_pos and valid_X_kaon_pos and valid_X_proton_pos:
    plt.subplot(3,2,2)
    # BIN DEFINITION
    bins_pi = np.logspace(np.log10(min(valid_X_pi_pos)), np.log10(max(valid_X_pi_pos)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_X_kaon_pos)), np.log10(max(valid_X_kaon_pos)), 60)
    bins_proton = np.logspace(np.log10(min(valid_X_proton_pos)), np.log10(max(valid_X_proton_pos)), 60)
    plot_histogram_lines(x_DA_pi_pos, bins_pi, color='skyblue', label='Pions', total_events=total_count_pos)
    plot_histogram_lines(x_DA_kaon_pos, bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_pos)
    plot_histogram_lines(x_DA_proton_pos, bins_proton, color='mediumorchid', label='Protons', total_events=total_count_pos)
    plt.xscale('log')
    plt.xlabel('x_B')
    plt.xlim(1e-6, 1)
    plt.ylabel('Number of Particles')
    plt.title('Number of pi+, P and K+ vs x')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

else:
    print("No valid Q^2 values found ")

# z OBSERVATION
valid_z_pion_pos = [z for z in z_pion_pos if z > 0]
valid_z_kaon_pos = [z for z in z_kaon_pos if z > 0]
valid_z_proton_pos = [z for z in z_proton_pos if z >0]

if valid_z_pion_pos and valid_z_kaon_pos and valid_z_proton_pos:
    plt.subplot(3,2,3)
    bins_pi = np.logspace(np.log10(min(valid_z_pion_pos)), np.log10(max(valid_z_pion_pos)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_z_kaon_pos)), np.log10(max(valid_z_kaon_pos)), 60)
    bins_proton = np.logspace(np.log10(min(valid_z_proton_pos)), np.log10(max(valid_z_proton_pos)), 60)
    plot_histogram_lines(z_pion_pos, bins=bins_pi, color='skyblue', label='Pions', total_events=total_count_pos)
    plot_histogram_lines(z_kaon_pos, bins=bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_pos)
    plot_histogram_lines(z_proton_pos, bins=bins_proton, color='mediumorchid', label='Protons', total_events=total_count_pos)
    plt.xscale('log')
    plt.xlabel('z')
    plt.xlim(5e-5,1)
    plt.ylabel('Number of Particles')
    plt.title('Number of pi+, P and K+ vs z')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

else:
    print("No valid z values found")

#p_hT OBSERVATION
valid_pT_pi_pos = [pT for pT in PhT_pion_pos]
valid_pT_kaon_pos = [pT for pT in PhT_kaon_pos]
valid_pT_proton_pos = [pT for pT in PhT_proton_pos]

if valid_pT_pi_pos and valid_pT_kaon_pos and valid_pT_proton_pos:
    plt.subplot(3,2,4)
    bins_pi = np.logspace(np.log10(min(valid_pT_pi_pos)), np.log10(max(valid_pT_pi_pos)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_pT_kaon_pos)), np.log10(max(valid_pT_kaon_pos)), 60)
    bins_proton = np.logspace(np.log10(min(valid_pT_proton_pos)), np.log10(max(valid_pT_proton_pos)), 60)
    plot_histogram_lines(PhT_pion_pos, bins=bins_pi, color='skyblue', label='Pions', total_events=total_count_pos)
    plot_histogram_lines(PhT_kaon_pos, bins=bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_pos)
    plot_histogram_lines(PhT_proton_pos, bins=bins_proton, color='mediumorchid', label='Protons', total_events=total_count_pos)
    plt.xscale('log')
    plt.xlim(1e-2,1e5)
    plt.xlabel('PhT [GeV^2]')
    plt.ylabel('Normalized Number of Particles')
    plt.title('Number of pi+, P and K+ vs p_T')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
else:
    print("No valid pT values found")


# RAPIDITY | SET POSITIVE ALONG THE Z-AXIS
plt.subplot(3,2,5)
plot_histogram_lines(pion_rapidity_pos, bins=50, color='skyblue', label='Pions', total_events=total_count_pos)
plot_histogram_lines(kaon_rapidity_pos, bins=50, color='lightgreen', label='Kaons', total_events=total_count_pos)
plot_histogram_lines(proton_rapidity_pos, bins=50, color='mediumorchid', label='Protons', total_events=total_count_pos)
plt.xlim(-4, 4)
plt.xlabel('Rapidity | Angle')
plt.ylabel('Normalized Number of Particles')
plt.title('Number of pi+, P and K+ vs rapidity')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
rapidity_to_degrees = {-4: '178°', -2: '', -0.88: '135°', 0: '90°', 0.88: '45°', 2: '', 4: '2°'}      # TO HAVE A SMALL IDEA ABOUT THE CONNECTION OF THE RAPIDITY AND THE EMISSION ANGLE
tick_positions = [-4, -2, -0.88, 0, 0.88, 2, 4]
tick_labels = ['{}  \n  {}'.format(rapidity, rapidity_to_degrees[rapidity]) for rapidity in tick_positions]
plt.xticks(tick_positions, tick_labels)

valid_p_pi_pos = [p for p in pion_mom_pos]
valid_p_kaon_pos = [p for p in kaon_mom_pos]
valid_p_proton_pos = [p for p in proton_mom_pos]

if valid_p_pi_pos and valid_p_kaon_pos and valid_p_proton_pos:
    plt.subplot(3,2,6)
    bins_pi = np.logspace(np.log10(min(valid_p_pi_pos)), np.log10(max(valid_p_pi_pos)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_p_kaon_pos)), np.log10(max(valid_p_kaon_pos)), 60)
    bins_proton = np.logspace(np.log10(min(valid_p_proton_pos)), np.log10(max(valid_p_proton_pos)), 60)
    plot_histogram_lines(pion_mom_pos, bins=bins_pi, color='skyblue', label='Pions', total_events=total_count_pos)
    plot_histogram_lines(kaon_mom_pos, bins=bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_pos)
    plot_histogram_lines(proton_mom_pos, bins=bins_proton, color='mediumorchid', label='Protons', total_events=total_count_pos)
    plt.xscale('log')
    plt.xlim(1e-1,275)
    plt.xlabel('Momentum [GeV^2]')
    plt.ylabel('Normalized Number of Particles')
    plt.title('Number of pi- and K- vs p_T')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
else:
    print("No valid momentum values found")

plt.tight_layout()

#_______________________________________________________________________ NEGATIVE GRAPHS ___________________________________________________________________________________________________________

plt.figure(figsize=(12, 8))
# FILTER POSSIBLE INVALID VALUE OF Q2 
valid_Q2_pi_neg = [q2 for q2 in Q2_DA_pi_neg if q2 > 0]
valid_Q2_kaon_neg = [q2 for q2 in Q2_DA_kaon_neg if q2 > 0]

# IN CASE THERE ARE NO PROBLEMS WITH THE VALUE 
if valid_Q2_pi_neg and valid_Q2_kaon_neg:
    # GENERATION OF A GREAT FIGURE TO FOUR GRAPHS
    plt.subplot(3,2,1)
    # DEFFINITION OF THE BINS OF Q2
    bins_pi = np.logspace(np.log10(min(valid_Q2_pi_neg)), np.log10(max(valid_Q2_pi_neg)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_Q2_kaon_neg)), np.log10(max(valid_Q2_kaon_neg)), 60)
    # HISTOGRAM PLOT | HERE THE FUNCTION TO OBTAIN A LINE HISTOGRAM IS USED
    plot_histogram_lines(Q2_DA_pi_neg, bins=bins_pi, color='skyblue', label=f'Pions ({round(norm_pion_neg,3)})', total_events=total_count_neg)
    plot_histogram_lines(Q2_DA_kaon_neg, bins=bins_kaon, color='lightgreen', label=f'Kaons ({round(norm_kaon_neg,3)})', total_events=total_count_neg)
    plt.xscale('log')
    plt.xlim(1,1e7)
    plt.xlabel('Q^2 [GeV^2]')
    plt.ylabel('Number of Particles')
    plt.title('Number of pi- and K- vs Q^2')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
else:
    print("No valid Q^2 values found")

# x (BJORKEN) VARIABLES OBSERVATION
valid_X_pi_neg = [x for x in x_DA_pi_neg if x > 0]
valid_X_kaon_neg = [x for x in x_DA_kaon_neg if x > 0]

if valid_X_pi_neg and valid_X_kaon_neg:
    plt.subplot(3,2,2)
    # BIN DEFINITION
    bins_pi = np.logspace(np.log10(min(valid_X_pi_neg)), np.log10(max(valid_X_pi_neg)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_X_kaon_neg)), np.log10(max(valid_X_kaon_neg)), 60)
    plot_histogram_lines(x_DA_pi_neg, bins_pi, color='skyblue', label='Pions', total_events=total_count_neg)
    plot_histogram_lines(x_DA_kaon_neg, bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_neg)
    plt.xscale('log')
    plt.xlabel('x_B')
    plt.xlim(1e-6, 1)
    plt.ylabel('Number of Particles')
    plt.title('Number of pi- and K- vs x')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

else:
    print("No valid Q^2 values found ")

# z OBSERVATION
valid_z_pion_neg = [z for z in z_pion_neg if z > 0]
valid_z_kaon_neg = [z for z in z_kaon_neg if z > 0]

if valid_z_pion_neg and valid_z_kaon_neg:
    plt.subplot(3,2,3)
    bins_pi = np.logspace(np.log10(min(valid_z_pion_neg)), np.log10(max(valid_z_pion_neg)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_z_kaon_neg)), np.log10(max(valid_z_kaon_neg)), 60)
    plot_histogram_lines(z_pion_neg, bins=bins_pi, color='skyblue', label='Pions', total_events=total_count_neg)
    plot_histogram_lines(z_kaon_neg, bins=bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_neg)
    plt.xscale('log')
    plt.xlabel('z')
    plt.xlim(5e-5,1)
    plt.ylabel('Number of Particles')
    plt.title('Number of pi- and K- vs z')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

else:
    print("No valid z values found")

#p_hT OBSERVATION
valid_pT_pi_neg = [pT for pT in PhT_pion_neg]
valid_pT_kaon_neg = [pT for pT in PhT_kaon_neg]

if valid_pT_pi_neg and valid_pT_kaon_neg:
    plt.subplot(3,2,4)
    bins_pi = np.logspace(np.log10(min(valid_pT_pi_neg)), np.log10(max(valid_pT_pi_neg)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_pT_kaon_neg)), np.log10(max(valid_pT_kaon_neg)), 60)
    plot_histogram_lines(PhT_pion_neg, bins=bins_pi, color='skyblue', label='Pions', total_events=total_count_neg)
    plot_histogram_lines(PhT_kaon_neg, bins=bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_neg)
    plt.xscale('log')
    plt.xlim(1e-2,1e5)
    plt.xlabel('PhT [GeV^2]')
    plt.ylabel('Number of Particles')
    plt.title('Number of pi- and K- vs p_T')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
else:
    print("No valid pT values found")

# RAPIDITY
plt.subplot(3,2,5)
plot_histogram_lines(pion_rapidity_neg, bins=50, color='skyblue', label='Pions', total_events=total_count_pos)
plot_histogram_lines(kaon_rapidity_neg, bins=50, color='lightgreen', label='Kaons', total_events=total_count_pos)
plt.xlim(-4, 4)
plt.xlabel('Rapidity | Angle')
plt.ylabel('Number of Particles')
plt.title('Number of pi- and K- vs rapidity')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
rapidity_to_degrees = {-4: '178°', -2: '', -0.88: '135°', 0: '90°', 0.88: '45°', 2: '', 4: '2°'}
tick_positions = [-4, -2, -0.88, 0, 0.88, 2, 4]
tick_labels = ['{}  \n  {}'.format(rapidity, rapidity_to_degrees[rapidity]) for rapidity in tick_positions]
plt.xticks(tick_positions, tick_labels)

plt.tight_layout()

valid_p_pi_neg = [p for p in pion_mom_neg]
valid_p_kaon_neg = [p for p in kaon_mom_neg]

if valid_p_pi_neg and valid_p_kaon_neg:
    plt.subplot(3,2,6)
    bins_pi = np.logspace(np.log10(min(valid_p_pi_neg)), np.log10(max(valid_p_pi_neg)), 60)
    bins_kaon = np.logspace(np.log10(min(valid_p_kaon_neg)), np.log10(max(valid_p_kaon_neg)), 60)
    plot_histogram_lines(pion_mom_neg, bins=bins_pi, color='skyblue', label='Pions', total_events=total_count_neg)
    plot_histogram_lines(kaon_mom_neg, bins=bins_kaon, color='lightgreen', label='Kaons', total_events=total_count_neg)
    plt.xscale('log')
    plt.xlim(1e-1,275)
    plt.xlabel('Momentum [GeV^2]')
    plt.ylabel('Number of Particles')
    plt.title('Number of pi- and K- vs p_T')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
else:
    print("No valid momentum values found")


#____________________________________________________________________ ANGLE PRODUCTION ______________________________________________________________________________________________________________

plt.figure(figsize=(13, 7))
plt.title('Production of Pions, Kaons and Protons in function of the Angle')   # GENERIC TITLE OF THE FIGURE
plt.axis('off')                                                                # TO REMOVE THE BORDER OF THE GRAPH

# POSITIVE PION
plt.subplot(2,3,1, polar=True)
plt.hist(pion_angles_pos_deg, bins=1500, color='skyblue', edgecolor='grey', label=f'Pions+ ({round(norm_pion_pos,3)})', density='True')
plt.xlim(0, np.pi)
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
plt.gca().set_rgrids([0, 0.05, 0.1, 0.15, 0.2, 0.25], labels=[0, "", 0.1, "", 0.2, ""])
plt.legend()

# NEGATIVE PION
plt.subplot(2,3,4, polar=True)
plt.hist(pion_angles_neg_deg, bins=1500, color='skyblue', label=f'Pions- ({round(norm_pion_neg,3)})', edgecolor='grey', density='True')
plt.xlim(0, np.pi)
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
#plt.gca().set_rgrids([0, 0.05, 0.1, 0.15, 0.2])
plt.gca().set_rgrids([0, 0.05, 0.1, 0.15, 0.2, 0.25], labels=[0, "", 0.1, "", 0.2, ""])
plt.legend()

# POSITIVE KAON
plt.subplot(2,3,2, polar=True)
plt.hist(kaon_angles_pos_deg, bins=1500, color='lightgreen', label=f'Kaons+ ({round(norm_kaon_pos,3)})', edgecolor='grey', density='True')
plt.xlim(0, np.pi)
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
plt.gca().set_rgrids([0, 0.05, 0.1, 0.15, 0.2, 0.25], labels=[0, "", 0.1, "", 0.2, ""])
plt.legend()

# NEGATIVE KAON
plt.subplot(2,3,5, polar=True)
plt.hist(kaon_angles_neg_deg, bins=1500, color='lightgreen', label=f'Kaons- ({round(norm_kaon_neg,3)})', edgecolor='grey', density='True')
plt.xlim(0, np.pi)
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
plt.gca().set_rgrids([0, 0.05, 0.1, 0.15, 0.2, 0.25], labels=[0, "", 0.1, "", 0.2, ""])
plt.legend()

# PROTON
plt.subplot(2,3,3, polar=True)
plt.hist(proton_angles_pos_deg, bins=1500, color='mediumorchid', label=f'Protons ({round(norm_proton_pos,3)})', edgecolor='grey', density='True')
plt.xlim(0, np.pi)
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
plt.gca().set_rgrids([0, 1, 2, 3, 4])
plt.legend()

plt.tight_layout()

plt.show()

