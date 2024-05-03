import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
from scipy.interpolate import make_interp_spline 
from scipy.interpolate import griddata


#__________________________________________________________________________________________________________________________________
# DEFINITION OF A FUNCTIO ABLE TO READ MY PYTHIA6.4 OUTPUT AND COLLECT INPORTANTDATA INTO AN ARRAY (OR MORE)
# THE COUNT OF THE TOTAL NUMBER OF EVENT IS NECESSARY TO NORMALIZE THE DATA AND HAVE AN ESTIMATE BY EVENT OF THE MEASUREMENTS
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
#filename = "pythia_ep_18x275_q2_1e-9_1_rid.txt"                           # SMALL DATA WITH LOWER ENEGY RANGE
#filename = "pythia_ep_noradcor_18x275_q2_0.1_10000000_run1_100K.txt"      # ORIGINAL FILE BUT REQUIRE TOO MUCH TIME
#filename = "pythia_ep_noradcor_18x275_q2_0.1_10000000_run1_10K.txt"       # 10K EVENT SET AT HIGH Q2
#filename = "pythia_ep_noradcor_18x275_q2_1e-9_1.0_run1_10K.txt"           # 10K EVENT SET AT LOW Q2
particles, particles_scatter, event_count = read_pythia_output(filename)

# COUNTERS FOR ALL THE CONSIDERED PARTICLE TYPES
electron_count, electron_count_i, electron_count_un, electron_count_f = 0, 0, 0, 0
pion_count, pion_count_i, pion_count_un, pion_count_f = 0, 0, 0, 0
kaon_count, kaon_count_i, kaon_count_un, kaon_count_f = 0, 0, 0, 0
proton_count, proton_count_i, proton_count_un, proton_count_f = 0, 0, 0, 0
neutron_count, neutron_count_i, neutron_count_un, neutron_count_f = 0, 0, 0, 0
rho_count, rho_count_i, rho_count_un, rho_count_f = 0, 0, 0, 0
gluon_count, gluon_count_i, gluon_count_un, gluon_count_f = 0, 0, 0, 0
photon_count, photon_count_i, photon_count_un, photon_count_f = 0, 0, 0, 0
other_count, other_count_i, other_count_un, other_count_f = 0, 0, 0, 0
up_count, up_count_i, up_count_un, up_count_f = 0, 0, 0, 0
down_count, down_count_i, down_count_un, down_count_f = 0, 0, 0, 0
charm_count, charm_count_i, charm_count_un, charm_count_f = 0, 0, 0, 0
strange_count, strange_count_i, strange_count_un, strange_count_f = 0, 0, 0, 0
diq_count, diq_count_i, diq_count_un, diq_count_f = 0, 0, 0, 0
diq2_count, diq2_count_i, diq2_count_un, diq2_count_f = 0, 0, 0, 0
delta_count, delta_count_i, delta_count_un, delta_count_f = 0, 0, 0, 0
lambda_count, lambda_count_i, lambda_count_un, lambda_count_f = 0, 0, 0, 0
epsilon_count, epsilon_count_i, epsilon_count_un, epsilon_count_f = 0, 0, 0, 0
eta_count, eta_count_i, eta_count_un, eta_count_f = 0, 0, 0, 0
omegaphi_count, omegaphi_count_i, omegaphi_count_un, omegaphi_count_f = 0, 0, 0, 0
frag_count, frag_count_i, frag_count_un, frag_count_f = 0, 0, 0, 0
difractive_count, difractive_count_i, difractive_count_un, difractive_count_f = 0, 0, 0, 0
strangeK_count, strangeK_count_i, strangeK_count_un, strangeK_count_f = 0, 0, 0, 0
xi_count, xi_count_i, xi_count_un, xi_count_f = 0, 0, 0, 0
fragm_up, fragm_down, fragm_uu, fragm_ud = 0, 0, 0, 0

# IDENTIFICATION OF THE TOTAL PARTICLES OF THE EVENT
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)
    event_number = particle[9]              # RECALL THE EVENT NUMBER TO REMOVE THE SCATTERED PARTICLES SINCE ARE NOT PRODUCTS OF THE INTERACTION
    if event_number not in (0, 3, 4, 5):    # THE LINE 1 AND 2 ARE REMOVED FROM THE READOUT PROCESS, ALSO THE LINES OF THE EVENT AND SCATTERED PARTICLES ARE REMOVED FROM THE IDENTIFICATION SINCE IS NOT NEEDED
        if particle_name == "Electron" or particle_name == "Positron":  
            electron_count += 1
        elif particle_name == "Pion+" or particle_name == "Pion-":
            pion_count += 1
        elif particle_name == "Pion0":
            pion_count += 1
        elif particle_name == "Kaon+" or particle_name == "Kaon-" or particle_name == "Kaon0" or particle_name == "Kaon0bar":
            kaon_count += 1
        elif particle_name == "Proton+" or particle_name == "Proton-":
            proton_count += 1
        elif particle_name == "Neutron" or particle_name == "Neutronbar":
            neutron_count += 1
        elif particle_name == "Up" or particle_name == "Upbar":
            up_count += 1
        elif particle_name == "Down" or particle_name == "Downbar":
            down_count += 1
        elif particle_name == "Photon":
            photon_count += 1
        elif particle_name == "Gluon":
            gluon_count += 1
        elif particle_name == "Cluster" or particle_name == "String":
            frag_count += 1
        else:
            other_count += 1

# IDENTIFICATION OF THE INITIAL STATE PARTICLES
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)
    origin = particle[10]                                            # THE ORIGIN IS NEEDED TO REGOGNIZE THE PARTICLES WHICH ARE PRODUCED DIRECTLY FROM THE SCATTERED PARTICLES
    event_number = particle[9] 
    if event_number not in (0, 3, 4, 5) and origin in (3, 4, 5):     # FOR THAT REASON THE ORIGIN IS SETTED TO BE 3,4,5
        if particle_name == "Electron" or particle_name == "Positron":
            electron_count_i += 1
        elif particle_name == "Pion+" or particle_name == "Pion-":
            pion_count_i += 1
        elif particle_name == "Pion0":
            pion_count_i += 1
        elif particle_name == "Kaon+" or particle_name == "Kaon-" or particle_name == "Kaon0" or particle_name == "Kaon0bar":
            kaon_count_i += 1
        elif particle_name == "Proton+" or particle_name == "Proton-":
            proton_count_i += 1
        elif particle_name == "Neutron" or particle_name == "Neutronbar":
            neutron_count_i += 1
        elif particle_name == "Up" or particle_name == "Upbar":
            up_count_i += 1
        elif particle_name == "Down" or particle_name == "Downbar":
            down_count_i += 1
        elif particle_name == "Photon":
            photon_count_i += 1
        elif particle_name == "Gluon":
            gluon_count_i += 1
        elif particle_name == "Cluster" or particle_name == "String":
            frag_count_i += 1
        else:
            other_count_i += 1

# ZOOM IN THE OTHER PARTICLES | INITIAL STATE
for particle in particles:    # HERE A ZOOM IN THE 'OTHER' COLUMN IS PERFORMED SINCE BRING AN INTERESTING CONTRIBUTION IN THE INITIAL STATE
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)
    event_number = particle[9]  
    origin = particle[10]
    if event_number not in (0, 3, 4, 5) and origin in (3, 4, 5): 
        if particle_name == "Rho+" or particle_name == "Rho-" or particle_name == "Rho0":
            rho_count_i += 1
        elif particle_name == "Strange" or particle_name == "Strangebar":
            strange_count_i += 1
        elif particle_name == "Charm" or particle_name == "Charmbar":
            charm_count_i += 1
        elif particle_name == "ud0" or particle_name == "ud1":
            diq_count_i += 1
        elif particle_name == "uu":
            diq2_count_i += 1
        elif particle_name == "Delta++" or particle_name == "Delta+" or particle_name == "Delta0" or particle_name == "Delta-":
            delta_count_i += 1
        elif particle_name == "Delta++bar" or particle_name == "Delta+bar" or particle_name == "Delta0bar" or particle_name == "Delta-bar":
            delta_count_i += 1
        elif particle_name == "Lambda" or particle_name == "CLambda+" or particle_name == "Lambdabar" or particle_name == "CLambda+bar":
            lambda_count_i += 1
        elif particle_name == "Epsilon+" or particle_name == "Epsilon0" or particle_name == "Epsilon-" or particle_name == "Epsilon*+" or particle_name == "Epsilon*0" or particle_name == "Epsilon*-":
            epsilon_count_i += 1
        elif particle_name == "Epsilon+bar" or particle_name == "Epsilon0bar" or particle_name == "Epsilon-bar" or particle_name == "Epsilon*+bar" or particle_name == "Epsilon*0bar" or particle_name == "Epsilon*-bar":
            epsilon_count_i += 1
        elif particle_name == "Eta":
            eta_count_i += 1
        elif particle_name == "Omega" or particle_name == "Phi":
            omegaphi_count_i += 1
        elif particle_name == "Cluster" or particle_name == "String":
            frag_count_i += 1
        elif particle_name == "Dif" or particle_name == "pi_dif" or particle_name == "w_diff" or particle_name == "phi_dif" or particle_name == "jpsi_dif" or particle_name == "n_dif" or particle_name == "p_dif":
            difractive_count_i += 1
        elif particle_name == "K0L" or particle_name == "K0S" or particle_name == "K*+" or particle_name == "K*0":
            strangeK_count_i += 1
        elif particle_name == "K0Lbar" or particle_name == "K0Sbar" or particle_name == "K*+bar" or particle_name == "K*0bar":
            strangeK_count_i += 1    
        elif particle_name == "Xi0" or particle_name == "Xi-" or particle_name == "Xi*0" or particle_name == "Xi*-":
            xi_count_i += 1
        elif particle_name == "Xi0bar" or particle_name == "Xi-bar" or particle_name == "Xi*0bar" or particle_name == "Xi*-bar":
            xi_count_i += 1


# IDENTIFICATION OF THE UNSTABLE PARTICLES
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)
    event_number = particle[9]  
    dau1 = particle[11]                                                 # HERE ARE INTRODUCED THE DAUGHTER 1,2 NEEDED TO OBSERVE IF A PARTICLES GENERATE OTHER PARTICLES,
    dau2 = particle[12]                                                 # HENCE, DECAY AND WILL NOT BE PRESENT IN THE FINAL STATE
    if event_number not in (0, 3, 4, 5) and (dau1 != 0 and dau2 != 0):  # SO ALL THE PARTICLES WHICH HAS DAUGHTER ARE CONSIDERED
        if particle_name == "Electron" or particle_name == "Positron":
            electron_count_un += 1
        elif particle_name == "Pion+" or particle_name == "Pion-":
            pion_count_un += 1
        elif particle_name == "Pion0":
            pion_count_un += 1
        elif particle_name == "Kaon+" or particle_name == "Kaon-" or particle_name == "Kaon0" or particle_name == "Kaon0bar":
            kaon_count_un += 1
        elif particle_name == "Proton+" or particle_name == "Proton-":
            proton_count_un += 1
        elif particle_name == "Neutron" or particle_name == "Neutronbar":
            neutron_count_un += 1
        elif particle_name == "Up" or particle_name == "Upbar":
            up_count_un += 1
        elif particle_name == "Down" or particle_name == "Downbar":
            down_count_un += 1
        elif particle_name == "Photon":
            photon_count_un += 1
        elif particle_name == "Gluon":
            gluon_count_un += 1                
        else:
            other_count_un += 1

# ZOOM IN THE OTHER PARTICLES | UNSTABLE PARTICLES
for particle in particles:    # HERE THE 'OTHER' COLUMN IS DOMINAND SO A ZOOM IS MANDATORY TO UNNDERSTAND THE PROCESS BEHAVIOUR
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)
    event_number = particle[9]  
    dau1 = particle[11]
    dau2 = particle[12]
    if event_number not in (0, 3, 4, 5) and (dau1 != 0 and dau2 != 0): 
        if particle_name == "Rho+" or particle_name == "Rho-" or particle_name == "Rho0":
            rho_count_un += 1
        elif particle_name == "Strange" or particle_name == "Strangebar":
            strange_count_un += 1
        elif particle_name == "Charm" or particle_name == "Charmbar":
            charm_count_un += 1
        elif particle_name == "ud0" or particle_name == "ud1" or particle_name == "uu":
            diq_count_un += 1
        elif particle_name == "Delta++" or particle_name == "Delta+" or particle_name == "Delta0" or particle_name == "Delta-":
            delta_count_un += 1
        elif particle_name == "Delta++bar" or particle_name == "Delta+bar" or particle_name == "Delta0bar" or particle_name == "Delta-bar":
            delta_count_un += 1
        elif particle_name == "Lambda" or particle_name == "CLambda+" or particle_name == "Lambdabar" or particle_name == "CLambda+bar":
            lambda_count_un += 1
        elif particle_name == "Epsilon+" or particle_name == "Epsilon0" or particle_name == "Epsilon-" or particle_name == "Epsilon*+" or particle_name == "Epsilon*0" or particle_name == "Epsilon*-":
            epsilon_count_un += 1
        elif particle_name == "Epsilon+bar" or particle_name == "Epsilon0bar" or particle_name == "Epsilon-bar" or particle_name == "Epsilon*+bar" or particle_name == "Epsilon*0bar" or particle_name == "Epsilon*-bar":
            epsilon_count_un += 1
        elif particle_name == "Eta":
            eta_count_un += 1
        elif particle_name == "Omega" or particle_name == "Phi":
            omegaphi_count_un += 1
        elif particle_name == "Cluster" or particle_name == "String":
            frag_count_un += 1
        elif particle_name == "Dif" or particle_name == "pi_dif" or particle_name == "w_diff" or particle_name == "phi_dif" or particle_name == "jpsi_dif" or particle_name == "n_dif" or particle_name == "p_dif":
            difractive_count_un += 1
        elif particle_name == "K0L" or particle_name == "K0S" or particle_name == "K*+" or particle_name == "K*0":
            strangeK_count_un += 1
        elif particle_name == "K0Lbar" or particle_name == "K0Sbar" or particle_name == "K*+bar" or particle_name == "K*0bar":
            strangeK_count_un += 1    
        elif particle_name == "Xi0" or particle_name == "Xi-" or particle_name == "Xi*0" or particle_name == "Xi*-":
            xi_count_un += 1
        elif particle_name == "Xi0bar" or particle_name == "Xi-bar" or particle_name == "Xi*0bar" or particle_name == "Xi*-bar":
            xi_count_un += 1

# IDENTIFICATION OF THE PARTICLES IN THE FINAL STATE
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)                  # IF A PARTICLE DOES NOT DECAY INTO OTHER WILL AND DOES NOT GENERATE OTHER PARTICLES, IT MEANS THAT IT WILL BE 'STABLE'
    event_number = particle[9]                                   # ('STABLE' -> WITH A LIFETIME LONGER THAN THE OBSERVATION TIME OF THE EVENT)
    dau1 = particle[11]                                          # AND THEREFORE WILL BE AN ELEMENT OF THE FINAL STATE
    dau2 = particle[12]
    if event_number not in (0, 3, 4, 5) and (dau1 == 0 or dau2 == 0):   
        if particle_name == "Electron" or particle_name == "Positron":  
            electron_count_f += 1                                      
        elif particle_name == "Pion+" or particle_name == "Pion-":
            pion_count_f += 1
        elif particle_name == "Pion0":
            pion_count_f += 1
        elif particle_name == "Kaon+" or particle_name == "Kaon-" or particle_name == "Kaon0" or particle_name == "Kaon0bar":
            kaon_count_f += 1
        elif particle_name == "Proton+" or particle_name == "Proton-":
            proton_count_f += 1
        elif particle_name == "Neutron" or particle_name == "Neutronbar":
            neutron_count_f += 1
        elif particle_name == "Up" or particle_name == "Upbar":
            up_count_f += 1
        elif particle_name == "Down" or particle_name == "Downbar":
            down_count_f += 1
        elif particle_name == "Photon":
            photon_count_f += 1
        elif particle_name == "Gluon":
            gluon_count_f += 1
        elif particle_name == "Cluster" or particle_name == "String":
            frag_count_f += 1
        else:
            other_count_f += 1

# ZOOM IN THE OTHER PARTICLES | FINAL STATE
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)       # IN THE FINAL STATE THE 'OTHER' COLUMN IS NOT DOMINAT BUT COULD BE INTERESTING TO OBSERVE THE EVOLUTION OF THE OTHER PARTICLES
    event_number = particle[9]  
    dau1 = particle[11]
    dau2 = particle[12]
    if event_number not in (0, 3, 4, 5) and (dau1 == 0 and dau2 == 0): 
        if particle_name == "Rho+" or particle_name == "Rho-" or particle_name == "Rho0":
            rho_count_f += 1
        elif particle_name == "Strange" or particle_name == "Strangebar":
            strange_count_f += 1
        elif particle_name == "Charm" or particle_name == "Charmbar":
            charm_count_f += 1
        elif particle_name == "ud0" or particle_name == "ud1" or particle_name == "uu":
            diq_count_f += 1
        elif particle_name == "Delta++" or particle_name == "Delta+" or particle_name == "Delta0" or particle_name == "Delta-":
            delta_count_f += 1
        elif particle_name == "Delta++bar" or particle_name == "Delta+bar" or particle_name == "Delta0bar" or particle_name == "Delta-bar":
            delta_count_f += 1
        elif particle_name == "Lambda" or particle_name == "CLambda+" or particle_name == "Lambdabar" or particle_name == "CLambda+bar":
            lambda_count_f += 1
        elif particle_name == "Epsilon+" or particle_name == "Epsilon0" or particle_name == "Epsilon-" or particle_name == "Epsilon*+" or particle_name == "Epsilon*0" or particle_name == "Epsilon*-":
            epsilon_count_f += 1
        elif particle_name == "Epsilon+bar" or particle_name == "Epsilon0bar" or particle_name == "Epsilon-bar" or particle_name == "Epsilon*+bar" or particle_name == "Epsilon*0bar" or particle_name == "Epsilon*-bar":
            epsilon_count_f += 1
        elif particle_name == "Eta":
            eta_count_f += 1
        elif particle_name == "Omega" or particle_name == "Phi":
            omegaphi_count_f += 1
        elif particle_name == "Cluster" or particle_name == "String":
            frag_count_f += 1
        elif particle_name == "Dif" or particle_name == "pi_dif" or particle_name == "w_diff" or particle_name == "phi_dif" or particle_name == "jpsi_dif" or particle_name == "n_dif" or particle_name == "p_dif":
            difractive_count_f += 1
        elif particle_name == "K0L" or particle_name == "K0S" or particle_name == "K*+" or particle_name == "K*0":
            strangeK_count_f += 1
        elif particle_name == "K0Lbar" or particle_name == "K0Sbar" or particle_name == "K*+bar" or particle_name == "K*0bar":
            strangeK_count_f += 1    
        elif particle_name == "Xi0" or particle_name == "Xi-" or particle_name == "Xi*0" or particle_name == "Xi*-":
            xi_count_f += 1
        elif particle_name == "Xi0bar" or particle_name == "Xi-bar" or particle_name == "Xi*0bar" or particle_name == "Xi*-bar":
            xi_count_f += 1

# OBSERVATION OF THE QUARKS AND DIQUARKS HADRONIZATION
daughter1_up, daughter2_up = [], []                                              
count_daughters_up = {"pion": 0, "kaon": 0, "proton": 0, "neutron":0, "other":0}       # UP QUARK PRODUCTS
daughter1_down, daughter2_down = [], []
count_daughters_down = {"pion": 0, "kaon": 0, "proton": 0, "neutron":0, "other":0}     # DOWN QUARK PRODUCTS
daughter1_ud, daughter2_ud = [], []
count_daughters_ud = {"pion": 0, "kaon": 0, "proton": 0, "neutron":0, "other":0}       # DIQUARK (ud) PRODUCTS
daughter1_uu, daughter2_uu = [], []
count_daughters_uu = {"pion": 0, "kaon": 0, "proton": 0, "neutron":0, "other":0}       # DIQUARK (uu) PRODUCTS
# COUNTS OF THE HADRONIZATIONS
for particle in particles:
    pdg_code = particle[0]
    particle_name = identify_particle(pdg_code)
    event_number = particle[9]  
    origin = particle[10]
    dau1 = particle[11]
    dau2 = particle[12]
    if event_number not in (0,1,2,3,4,5) and origin in (3,4,5) and (dau1 != 0 and dau2 != 0):
        # UP CYCLE
        if particle_name == "Up" or particle_name == "Upbar": 
            fragm_up += 1
            daughter1_up = [d for d in particles if d[9] == dau1]
            daughter2_up = [d for d in particles if d[9] == dau2]
            # DAUGHTER 1 IDENTIFICATION
            for d1 in daughter1_up:
                pdg_code_d1 = d1[0]
                particle_name_d1 = identify_particle(pdg_code_d1)
                if particle_name_d1 in ["Pion+", "Pioni-", "Pion0"]:
                    count_daughters_up["pion"] += 1
                    break
                elif particle_name_d1 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_up["kaon"] += 1
                    break
                elif particle_name_d1 in ["Proton+", "Proton-"]:
                    count_daughters_up["proton"] += 1
                    break
                elif particle_name_d1 in ["Neutron", "Neutronbar"]:
                    count_daughters_up["neutron"] += 1
                    break
                else:
                    count_daughters_up["other"] += 1
                    break
            # DAUGHTER 2 IDENTIFICATION
            for d2 in daughter2_up:
                pdg_code_d2 = d2[0]
                particle_name_d2 = identify_particle(pdg_code_d2)
                if particle_name_d2 in ["Pion+", "Pion-", "Pion0"]:
                    count_daughters_up["pion"] += 1
                    break
                elif particle_name_d2 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_up["kaon"] += 1
                    break
                elif particle_name_d2 in ["Proton+", "Proton-"]:
                    count_daughters_up["proton"] += 1
                    break
                elif particle_name_d2 in ["Neutron", "Neutronbar"]:
                    count_daughters_up["neutron"] += 1
                    break
                else:
                    count_daughters_up["other"] += 1
                    break

        # DOWN CYCLE
        elif particle_name == "Down" or particle_name == "Downbar":
            fragm_down += 1
            daughter1_down = [d for d in particles if d[9] == dau1]
            daughter2_down = [d for d in particles if d[9] == dau2]
            # DAUGHTER 1 IDENTIFICATION
            for d1 in daughter1_down:
                pdg_code_d1 = d1[0]
                particle_name_d1 = identify_particle(pdg_code_d1)
                if particle_name_d1 in ["Pion+", "Pioni-", "Pion0"]:
                    count_daughters_down["pion"] += 1
                    break
                elif particle_name_d1 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_down["kaon"] += 1
                    break
                elif particle_name_d1 in ["Proton+", "Proton-"]:
                    count_daughters_down["proton"] += 1
                    break
                elif particle_name_d1 in ["Neutron", "Neutronbar"]:
                    count_daughters_down["neutron"] += 1
                    break
                else:
                    count_daughters_down["other"] += 1
                    break
            # DAUGHTER 2 IDENTIFICATION
            for d2 in daughter2_down:
                pdg_code_d2 = d2[0]
                particle_name_d2 = identify_particle(pdg_code_d2)
                if particle_name_d2 in ["Pion+", "Pion-", "Pion0"]:
                    count_daughters_down["pion"] += 1
                    break
                elif particle_name_d2 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_down["kaon"] += 1
                    break
                elif particle_name_d2 in ["Proton+", "Proton-"]:
                    count_daughters_down["proton"] += 1
                    break
                elif particle_name_d2 in ["Neutron", "Neutronbar"]:
                    count_daughters_down["neutron"] += 1
                    break
                else:
                    count_daughters_down["other"] += 1
                    break

        # DIQUARK (ud) CYCLE
        elif particle_name == "ud0" or particle_name == "ud1":
            fragm_ud += 1
            daughter1_ud = [d for d in particles if d[9] == dau1]
            daughter2_ud = [d for d in particles if d[9] == dau2]
            # DAUGHTER 1 IDENTIFICATION
            for d1 in daughter1_ud:
                pdg_code_d1 = d1[0]
                particle_name_d1 = identify_particle(pdg_code_d1)
                if particle_name_d1 in ["Pion+", "Pioni-", "Pion0"]:
                    count_daughters_ud["pion"] += 1
                    break
                elif particle_name_d1 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_ud["kaon"] += 1
                    break
                elif particle_name_d1 in ["Proton+", "Proton-"]:
                    count_daughters_ud["proton"] += 1
                    break
                elif particle_name_d1 in ["Neutron", "Neutronbar"]:
                    count_daughters_ud["neutron"] += 1
                    break
                else:
                    count_daughters_ud["other"] += 1
                    break
            # DAUGHTER 2 IDENTIFICATION
            for d2 in daughter2_ud:
                pdg_code_d2 = d2[0]
                particle_name_d2 = identify_particle(pdg_code_d2)
                if particle_name_d2 in ["Pion+", "Pion-", "Pion0"]:
                    count_daughters_ud["pion"] += 1
                    break
                elif particle_name_d2 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_ud["kaon"] += 1
                    break
                elif particle_name_d2 in ["Proton+", "Proton-"]:
                    count_daughters_ud["proton"] += 1
                    break
                elif particle_name_d2 in ["Neutron", "Neutronbar"]:
                    count_daughters_ud["neutron"] += 1
                    break
                else:
                    count_daughters_ud["other"] += 1
                    break

        # DIQUARK (uu) CYCLE
        elif particle_name == "uu":
            fragm_uu += 1
            daughter1_uu = [d for d in particles if d[9] == dau1]
            daughter2_uu = [d for d in particles if d[9] == dau2]
            # DAUGHTER 1 IDENTIFICATION
            for d1 in daughter1_uu:
                pdg_code_d1 = d1[0]
                particle_name_d1 = identify_particle(pdg_code_d1)
                if particle_name_d1 in ["Pion+", "Pioni-", "Pion0"]:
                    count_daughters_uu["pion"] += 1
                    break
                elif particle_name_d1 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_uu["kaon"] += 1
                    break
                elif particle_name_d1 in ["Proton+", "Proton-"]:
                    count_daughters_uu["proton"] += 1
                    break
                elif particle_name_d1 in ["Neutron", "Neutronbar"]:
                    count_daughters_uu["neutron"] += 1
                    break
                else:
                    count_daughters_uu["other"] += 1
                    break
            # DAUGHTER 2 IDENTIFICATION
            for d2 in daughter2_uu:
                pdg_code_d2 = d2[0]
                particle_name_d2 = identify_particle(pdg_code_d2)
                if particle_name_d2 in ["Pion+", "Pion-", "Pion0"]:
                    count_daughters_uu["pion"] += 1
                    break
                elif particle_name_d2 in ["Kaon+", "Kaon-", "Kaon0", "Kaon0bar"]:
                    count_daughters_uu["kaon"] += 1
                    break
                elif particle_name_d2 in ["Proton+", "Proton-"]:
                    count_daughters_uu["proton"] += 1
                    break
                elif particle_name_d2 in ["Neutron", "Neutronbar"]:
                    count_daughters_uu["neutron"] += 1
                    break
                else:
                    count_daughters_uu["other"] += 1
                    break

# NORMALIZATION OF THE COUNTS
def norm_particle_count(particle_c, event_c):
    norm_count = particle_c/event_c
    return norm_count
# NORMALIZED TOTAL COUNTS
norm_pion = norm_particle_count(pion_count, event_count)
norm_proton = norm_particle_count(proton_count, event_count)
norm_neutron = norm_particle_count(neutron_count, event_count)
norm_kaon = norm_particle_count(kaon_count, event_count)
norm_up = norm_particle_count(up_count, event_count)
norm_down = norm_particle_count(down_count, event_count)
norm_other = norm_particle_count(other_count, event_count)
norm_electron =  norm_particle_count(electron_count, event_count)
norm_photon =  norm_particle_count(photon_count, event_count)
norm_gluon = norm_particle_count(gluon_count, event_count) 
# INITIAL
norm_pion_i = norm_particle_count(pion_count_i, event_count)
norm_proton_i = norm_particle_count(proton_count_i, event_count)
norm_neutron_i = norm_particle_count(neutron_count_i, event_count)
norm_kaon_i = norm_particle_count(kaon_count_i, event_count)
norm_up_i = norm_particle_count(up_count_i, event_count)
norm_down_i = norm_particle_count(down_count_i, event_count)
norm_other_i = norm_particle_count(other_count_i, event_count)
norm_electron_i =  norm_particle_count(electron_count_i, event_count)
norm_photon_i =  norm_particle_count(photon_count_i, event_count)
norm_gluon_i = norm_particle_count(gluon_count_i, event_count) 
# INITIAL | OTHER
norm_rho_i = norm_particle_count(rho_count_i, event_count)
norm_strange_i = norm_particle_count(strange_count_i, event_count)
norm_charm_i = norm_particle_count(charm_count_i, event_count)
norm_diq_i = norm_particle_count(diq_count_i, event_count)
norm_diq2_i = norm_particle_count(diq2_count_i, event_count)
norm_delta_i = norm_particle_count(delta_count_i, event_count)
norm_lambda_i = norm_particle_count(lambda_count_i, event_count)
norm_epsilon_i = norm_particle_count(epsilon_count_i, event_count)
norm_eta_i = norm_particle_count(eta_count_i, event_count)
norm_omegaphi_i = norm_particle_count(omegaphi_count_i, event_count)
norm_frag_i = norm_particle_count(frag_count_i, event_count)
norm_difractive_i = norm_particle_count(difractive_count_i, event_count)
norm_strangeK_i = norm_particle_count(strangeK_count_i, event_count)
norm_xi_i = norm_particle_count(xi_count_i, event_count)
# UNSTABLE
norm_pion_un = norm_particle_count(pion_count_un, event_count)
norm_proton_un = norm_particle_count(proton_count_un, event_count)
norm_neutron_un = norm_particle_count(neutron_count_un, event_count)
norm_kaon_un = norm_particle_count(kaon_count_un, event_count)
norm_up_un = norm_particle_count(up_count_un, event_count)
norm_down_un = norm_particle_count(down_count_un, event_count)
norm_other_un = norm_particle_count(other_count_un, event_count)
norm_electron_un =  norm_particle_count(electron_count_un, event_count)
norm_photon_un =  norm_particle_count(photon_count_un, event_count)
norm_gluon_un = norm_particle_count(gluon_count_un, event_count) 
# UNSTABLE | OTHER
norm_rho_un = norm_particle_count(rho_count_un, event_count)
norm_strange_un = norm_particle_count(strange_count_un, event_count)
norm_charm_un = norm_particle_count(charm_count_un, event_count)
norm_diq_un = norm_particle_count(diq_count_un, event_count)
norm_delta_un = norm_particle_count(delta_count_un, event_count)
norm_lambda_un = norm_particle_count(lambda_count_un, event_count)
norm_epsilon_un = norm_particle_count(epsilon_count_un, event_count)
norm_eta_un = norm_particle_count(eta_count_un, event_count)
norm_omegaphi_un = norm_particle_count(omegaphi_count_un, event_count)
norm_frag_un = norm_particle_count(frag_count_un, event_count)
norm_difractive_un = norm_particle_count(difractive_count_un, event_count)
norm_strangeK_un = norm_particle_count(strangeK_count_un, event_count)
norm_xi_un = norm_particle_count(xi_count_un, event_count)
# FINAL STATE
norm_pion_f = norm_particle_count(pion_count_f, event_count)
norm_proton_f = norm_particle_count(proton_count_f, event_count)
norm_neutron_f = norm_particle_count(neutron_count_f, event_count)
norm_kaon_f = norm_particle_count(kaon_count_f, event_count)
norm_up_f = norm_particle_count(up_count_f, event_count)
norm_down_f = norm_particle_count(down_count_f, event_count)
norm_other_f = norm_particle_count(other_count_f, event_count)
norm_electron_f =  norm_particle_count(electron_count_f, event_count)
norm_photon_f =  norm_particle_count(photon_count_f, event_count)
norm_gluon_f = norm_particle_count(gluon_count_f, event_count) 
# FINAL STATE | OTHER
norm_rho_f = norm_particle_count(rho_count_f, event_count)
norm_strange_f = norm_particle_count(strange_count_f, event_count)
norm_charm_f = norm_particle_count(charm_count_f, event_count)
norm_diq_f = norm_particle_count(diq_count_f, event_count)
norm_delta_f = norm_particle_count(delta_count_f, event_count)
norm_lambda_f = norm_particle_count(lambda_count_f, event_count)
norm_epsilon_f = norm_particle_count(epsilon_count_f, event_count)
norm_eta_f = norm_particle_count(eta_count_f, event_count)
norm_omegaphi_f = norm_particle_count(omegaphi_count_f, event_count)
norm_frag_f = norm_particle_count(frag_count_f, event_count)
norm_difractive_f = norm_particle_count(difractive_count_f, event_count)
norm_strangeK_f = norm_particle_count(strangeK_count_f, event_count)
norm_xi_f = norm_particle_count(xi_count_f, event_count)

# FUNCTION TO HAVE LINE HISTOGRAMS
def plot_histogram_lines(data, bins, color, label, linewidth=2):
    hist, bins_edges = np.histogram(data, bins=bins)
    bins_centers = (bins_edges[:-1] + bins_edges[1:]) / 2
    plt.plot(bins_centers, hist, color=color, label=label, linewidth=linewidth)

# DEFINITION OF MY VARIABLES
# ENERGY
pion_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"]]
pion_charged_en = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-"]]
proton_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5]
kaon_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"]]
pion_mass = [particle[6] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"]]
kaon_mass = [particle[6] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"]]
proton_mass = [particle[6] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5]
# PION
pion_x = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"]]
pion_y = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"]]
pion_z = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Pion+", "Pion-", "Pion0"]]
# PROTON
proton_x = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5]
proton_y = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5]
proton_z = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Proton+", "Proton-"] and particle[9] != 5]
# KAON
kaon_x = [particle[2] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"]]
kaon_y = [particle[3] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"]]
kaon_z = [particle[4] for particle in particles if identify_particle(particle[0]) in ["Kaon+", "Kaon-", "Kaon0"]]
# LEPTON
electron_energies = [particle[5] for particle in particles if identify_particle(particle[0]) in ["Electron", "Positron"]] # LOW MULTIPLICITY
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
sc_phoyon_pz = [particle[4] for particle in particles if particle[9] == 4]
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
# OTHER QUANTITIES DEFINITION 
# PION
pion_transv_mom = [np.sqrt(px**2 + py**2) for px, py in zip (pion_x, pion_y)]  
pion_mom = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(pion_x, pion_y, pion_z)]
pion_long_mom = [np.sum(E - pz) for E, pz in zip(pion_energies, pion_z)]
pion_angles = [np.arccos(pz/p) for pz, p in zip(pion_z, pion_mom)]            # Z-AXIS
pion_angles_deg = [np.degrees(a) for a in pion_angles]
z_pion = [np.abs((Eph - np.abs(Phz))/(Ee - Pze)) for Eph, Phz, Ee, Pze in zip(pion_energies, pion_z, scattered_gamma_E, scattered_gamma_pz)]
PhT_pion = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(pion_z, scattered_el_pz, mod_q)]
# PROTON
proton_transv_mom = [np.sqrt(px**2 + py**2) + 0.0001 for px, py in zip (proton_x, proton_y)]  
proton_mom = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(proton_x, proton_y, proton_z)]
proton_long_mom = [np.sum(E - pz) for E, pz in zip(proton_energies, proton_z)]
proton_angles = [np.arccos(pz/p) for pz, p in zip(proton_z, proton_mom)]      # Z-AXIS
proton_angles_deg = [np.degrees(a) for a in proton_angles]
z_proton = [np.abs((Eph - np.abs(Phz))/(Ee - Pze)) for Eph, Phz, Ee, Pze in zip(proton_energies, proton_z, scattered_gamma_E, scattered_gamma_pz)]
PhT_proton = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(proton_z, scattered_el_pz, mod_q)]
# KAON
kaon_transv_mom = [np.sqrt(px**2 + py**2) for px, py in zip (kaon_x, kaon_y)]  
kaon_mom = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(kaon_x, kaon_y, kaon_z)]
kaon_long_mom = [np.sum(E - pz) for E, pz in zip(kaon_energies, kaon_z)]
kaon_angles = [np.arccos(pz/p) for pz, p in zip(kaon_z, kaon_mom)]            # Z-AXIS
kaon_angles_deg = [np.degrees(a) for a in kaon_angles]
z_kaon = [np.abs((Eph - np.abs(Phz))/(Ee - Pze)) for Eph, Phz, Ee, Pze in zip(kaon_energies, kaon_z, scattered_gamma_E, scattered_gamma_pz)]
PhT_kaon = [np.abs(Phz - ((Phz*z)/mz)) for Phz, z, mz in zip(kaon_z, scattered_el_pz, mod_q)]
# SCATTERED ELECTRON
sc_electron_mom = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(sc_electron_px, sc_electron_py, sc_electron_pz)]
sc_electron_ang = [np.arccos(pz/p) for pz, p in zip(sc_electron_pz, sc_electron_mom)]
sc_electron_ang_deg = [np.degrees(a) for a in sc_electron_ang]

# DOUBLE ANGLE METHOD | USED FOR THE CONSTRUCTION OF THE BJORKEN VARIABLE 'x' AND THE INELASTICITY 'y'
# PION
y_DA_pi = [np.tan(p/2)/(np.tan(p/2) + np.tan(t/2)) for p, t in zip(pion_angles, lepton_angles)]
Q2_DA_pi = [4*18*18*(1-abs(y))/(np.tan(t/2)**2) for y, t in zip(y_DA_pi, lepton_angles)]
x_DA_pi = [Q/(40*275*y*275) for Q, y in zip(Q2_DA_pi, y_DA_pi)]
# PROTON
y_DA_proton = [np.tan(p/2)/(np.tan(p/2) + np.tan(t/2)) for p, t in zip(proton_angles, lepton_angles)]
Q2_DA_proton = [4*18*18*(1-abs(y))/(np.tan(t/2)**2) for y, t in zip(y_DA_proton, lepton_angles)]
x_DA_proton = [Q/(40*275*y*275) for Q, y in zip(Q2_DA_proton, y_DA_proton)]
# KAON
y_DA_kaon = [np.tan(p/2)/(np.tan(p/2) + np.tan(t/2)) for p, t in zip(kaon_angles, lepton_angles)]
Q2_DA_kaon = [4*18*18*(1-abs(y))/(np.tan(t/2)**2) for y, t in zip(y_DA_kaon, lepton_angles)]
x_DA_kaon = [Q/(40*275*y*275) for Q, y in zip(Q2_DA_kaon, y_DA_kaon)]

# THEORETICAL CALCULATION  
# PION
xB_pion = [Q/(2*(Ep*E - px*(-pxp) - py*(-pyp) - pz*pzp)) for Q, Ep, E, px, pxp, py, pyp, pz, pzp in zip (Q2, pion_energies, scattered_gamma_E, pion_x, scattered_el_px, pion_y, scattered_el_py, pion_z, scattered_gamma_pz)]
y_pion = [(Ep*E - px*(-pxp) - py*(-pyp) - pz*pzp)/(18*Ep - (-18)*pz) for Q, Ep, E, px, pxp, py, pyp, pz, pzp in zip (Q2, pion_energies, scattered_gamma_E, pion_x, scattered_el_px, pion_y, scattered_el_py, pion_z, scattered_gamma_pz)]
Q2_pion = [Q*Ep/(Ep) for Q, Ep in zip (Q2, pion_energies)]

#________________________________________________________________________________HISTOGRAM FOR THE PARTICLES BEHAVIOURS_______________________________________________________________________________________________________

# FILTER POSSIBLE INVALID VALUE OF Q2 
valid_Q2_pi = [q2 for q2 in Q2_DA_pi if q2 > 0]
valid_Q2_kaon = [q2 for q2 in Q2_DA_kaon if q2 > 0]
valid_Q2_proton = [q2 for q2 in Q2_DA_proton if q2 >0]

# IN CASE THERE ARE NO PROBLES WITHT THE VALUE 
if valid_Q2_pi and valid_Q2_kaon and valid_Q2_proton:
    # GENERATION OF A GREAT FIGURE TO FOUR GRAPHS
    plt.figure(figsize=(12, 7))
    plt.subplot(2,2,1)
    # DEFFINITION OF THE BINS OF Q2
    bins_pi = np.logspace(np.log10(min(valid_Q2_pi)), np.log10(max(valid_Q2_pi)), 50)
    bins_kaon = np.logspace(np.log10(min(valid_Q2_kaon)), np.log10(max(valid_Q2_kaon)), 50)
    bins_proton = np.logspace(np.log10(min(valid_Q2_proton)), np.log10(max(valid_Q2_proton)), 50)
    # HISTOGRAM PLOT | HERE THE FUNCTION TO OBTAIN A LINE HISTOGRAM IS USED
    plot_histogram_lines(Q2_DA_pi, bins=bins_pi, color='skyblue', label='Pions')
    plot_histogram_lines(Q2_DA_kaon, bins=bins_kaon, color='lightgreen', label='Kaons')
    plot_histogram_lines(Q2_DA_proton, bins=bins_proton, color='mediumorchid', label='Protons')
    plt.xscale('log')
    plt.xlim(0.2,1e6)
    plt.xlabel('Q^2 [GeV^2]')
    plt.ylabel('Number of Particles')
    plt.title('Number of pi, P and K vs Q^2')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
    # IN CASE THE VALUES ARE NOT CORRECTED COULD BE INTERESTING TO UNDERSTAND ITS LOCATION
elif not valid_Q2_proton:
    print("No valid Q^2 values found for protons.")
elif not valid_Q2_kaon:
    print("No valid Q^2 values found foR kaons.")
elif not valid_Q2_pi:
    print("No valid Q^2 values found for pions")
else:
    print("No valid Q^2 values found")

# x (BJORKEN) VARIABLES OBSERVATION
valid_X_pi = [q2 for q2 in x_DA_pi if q2 > 0]
valid_X_kaon = [q2 for q2 in x_DA_kaon if q2 > 0]
valid_X_proton = [q2 for q2 in x_DA_proton if q2 >0]

if valid_X_pi and valid_X_kaon and valid_X_proton:
    plt.subplot(2,2,2)
    # BIN DEFINITION
    bins_pi = np.logspace(np.log10(min(valid_X_pi)), np.log10(max(valid_X_pi)), 50)
    bins_kaon = np.logspace(np.log10(min(valid_X_kaon)), np.log10(max(valid_X_kaon)), 50)
    bins_proton = np.logspace(np.log10(min(valid_X_proton)), np.log10(max(valid_X_proton)), 50)
    plot_histogram_lines(x_DA_pi, bins_pi, color='skyblue', label='Pions')
    plot_histogram_lines(x_DA_kaon, bins_kaon, color='lightgreen', label='Kaons')
    plot_histogram_lines(x_DA_proton, bins_proton, color='mediumorchid', label='Protons')
    plt.xscale('log')
    plt.xlabel('x_B')
    plt.xlim(1e-5, 1)
    plt.ylabel('Number of Particles')
    plt.title('Number of pi, P and K vs x')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

elif not valid_X_proton:
    print("No valid Q^2 values found for protons.")
elif not valid_X_kaon:
    print("No valid Q^2 values found for kaons.")
elif not valid_X_pi:
    print("No valid Q^2 values found for pions")
else:
    print("No valid Q^2 values found for kaons or proton.")

# z OBSERVATION
valid_z_pion = [z for z in z_pion if z > 0]
valid_z_kaon = [z for z in z_kaon if z > 0]
valid_z_proton = [z for z in z_proton if z >0]

if valid_z_pion and valid_z_kaon and valid_z_proton:
    plt.subplot(2,2,3)
    bins_pi = np.logspace(np.log10(min(valid_z_pion)), np.log10(max(valid_z_pion)), 50)
    bins_kaon = np.logspace(np.log10(min(valid_z_kaon)), np.log10(max(valid_z_kaon)), 50)
    bins_proton = np.logspace(np.log10(min(valid_z_proton)), np.log10(max(valid_z_proton)), 50)
    plot_histogram_lines(z_pion, bins=bins_pi, color='skyblue', label='Pions')
    plot_histogram_lines(z_kaon, bins=bins_kaon, color='lightgreen', label='Kaons')
    plot_histogram_lines(z_proton, bins=bins_proton, color='mediumorchid', label='Protons')
    plt.xscale('log')
    plt.xlabel('z')
    plt.xlim(1e-5,1e2)
    plt.ylabel('Number of Particles')
    plt.title('Number of pi, P and K vs z')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()

elif not valid_z_proton:
    print("No valid z values found for protons.")
elif not valid_z_kaon:
    print("No valid z values found for kaons.")
elif not valid_z_pion:
    print("No valid z values found for pions")
else:
    print("No valid z values found for kaons or proton.")

plt.tight_layout()

#p_hT OBSERVATION
valid_pT_pi = [pT for pT in PhT_pion]
valid_pT_kaon = [pT for pT in PhT_kaon]
valid_pT_proton = [pT for pT in PhT_proton]

if valid_pT_pi and valid_pT_kaon and valid_pT_proton:
    plt.subplot(2,2,4)
    bins_pi = np.logspace(np.log10(min(valid_pT_pi)), np.log10(max(valid_pT_pi)), 50)
    bins_kaon = np.logspace(np.log10(min(valid_pT_kaon)), np.log10(max(valid_pT_kaon)), 50)
    bins_proton = np.logspace(np.log10(min(valid_pT_proton)), np.log10(max(valid_pT_proton)), 50)
    plot_histogram_lines(PhT_pion, bins=bins_pi, color='skyblue', label='Pions')
    plot_histogram_lines(PhT_kaon, bins=bins_kaon, color='lightgreen', label='Kaons')
    plot_histogram_lines(PhT_proton, bins=bins_proton, color='mediumorchid', label='Protons')
    plt.xscale('log')
    plt.xlim(1e-2,5e5)
    plt.xlabel('PhT [GeV^2]')
    plt.ylabel('Number of Particles')
    plt.title('Number of pi, P and K vs p_T')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
   
elif not valid_pT_proton:
    print("No valid pT values found for protons.")
elif not valid_pT_kaon:
    print("No valid pT values found for kaons.")
elif not valid_pT_pi:
    print("No valid pT values found for pions")
else:
    print("No valid pT values found for kaons or proton.")

#____________________________________________________________________________GENERATION OF OBSERVATION GRAPH_____________________________________________________________________________________________

# BAR GRAPH FOR THE TOTAL MEASURE OF THE EVENT'S PARTICLES
particle_names = ['Pions', 'Kaons', 'Protons', 'Neutrons', 'Up', 'Down', 'Electrons', 'Photons', 'Others']
particle_counts = [pion_count, kaon_count, proton_count, neutron_count, up_count, down_count, electron_count, photon_count, other_count]

# BIG FIGURE TO CONTAIN 6 GRAPHS
plt.figure(figsize=(12, 8), facecolor = 'white')

# BAR PLOT
plt.subplot(3, 2, 1)
palette = sns.color_palette("Blues", len(particle_names))
plt.bar(particle_names, particle_counts, color=palette, alpha=0.7, width=0.5)
plt.ylabel('Number of Particles')
plt.title('Particle Identification Counts')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(particle_counts):
    plt.text(i, count + 0.05 * max(particle_counts), str(count), ha='center', va='bottom', fontsize=8)
total_particles = sum(particle_counts)
plt.text(0.5, 0.9, f'Total Particles: {total_particles}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)


# PLOT FOR NORMALIZED DATA
norm_particle_names = ['Pions', 'Kaons', 'Protons', 'Neutrons', 'Up', 'Down', 'Electrons', 'Photons', 'Others']
norm_particle_counts = [norm_pion, norm_kaon, norm_proton, norm_neutron, norm_up, norm_down, norm_electron, norm_photon, norm_other]

# PLOT
plt.subplot(3, 2, 2)
norm_palette = sns.color_palette("Reds", len(norm_particle_names))
plt.bar(norm_particle_names, norm_particle_counts, color=norm_palette, alpha=0.7, width=0.5)
plt.ylabel('Normalized number of Particles')
plt.title('Normalized Particle counts')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(norm_particle_counts):
    plt.text(i, count + 0.05 * max(norm_particle_counts), f"{round(count, 2)}", ha='center', va='bottom', fontsize=8)
norm_total_particles = sum(norm_particle_counts)
plt.text(0.5, 0.9, f'Particles per event: {round(norm_total_particles,2)}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)

# ANGULAR DISTRIBUTION OF THE SCATTERED ELECTRONS
plt.subplot(3,2,3)
plt.hist2d(x_DA_kaon, Q2_DA_kaon, bins=[np.logspace(np.log10(1e-5), np.log10(1), 30), np.logspace(np.log10(1e1), np.log10(1e6), 30)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('x_B')
plt.ylabel('Q^2')
plt.title('x_B vs Q2   |   TH   |   Pion')
plt.colorbar(label='Counts')

# PLOT OF x vs Q2 WITH THE DOUBLE ANGLE METHOD
plt.subplot(3,2,4)
plt.hist2d(x_DA_pi, Q2_DA_pi, bins=[np.logspace(np.log10(1e-5), np.log10(1), 30), np.logspace(np.log10(1e1), np.log10(1e6), 30)])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('x_B')
plt.ylabel('Q^2')
plt.title('x_B vs Q2   |   DA   |   Pions')
plt.colorbar(label='Counts')
plt.gca().set_facecolor('white')

# DENSITY OF THE PION'S MOMENTUM IN THE XY PLANE
plt.subplot(3, 2, 5)
sns.kdeplot(x=pion_x, y=pion_y, cmap="Reds", fill=True)
plt.xlabel('p_x [GeV]')
plt.ylabel('p_y [GeV]')
plt.title('Density of the Pions momentum')
plt.xlim(-1.2, 1.2)
plt.ylim(-0.8, 0.8)

# DENSITY OF THE PION'S MOMENTUM IN THE TRANSVERSE AND LONGITUDINAL ORIENTATION
plt.subplot(3,2,6)
sns.kdeplot(x=pion_long_mom, y=pion_transv_mom, cmap="Blues", fill=True)
plt.xlabel('p_L [GeV]')
plt.ylabel('Norm p_T [GeV]')
plt.title('Density of the Transv. vs Norm. Longit. momenta of Pions')
plt.xlim(-0.6,2)
plt.ylim(-0.4,1.35)

plt.tight_layout()  # TO GENERATE ORGANIZED GRAPH

#_____________________________________________________________IDENTIFICATION OF THE PARTICLES IN DIFFERENT MOMENTS OF THE EVENT________________________________________________________________________
'''
# GENERATION OF A LARGE CANVAS TO CONTAIN 6 GRAPHS
plt.figure(figsize=(13, 8))
# NORMALIZATION OF ALL THE DATA COLLECTED
norm_particle_names_i = ['Pions', 'Kaons', 'Protons', 'Neutrons', 'Up', 'Down', 'Gluon', 'Electron', 'Photon', 'Others']
norm_particle_counts_i = [norm_pion_i, norm_kaon_i, norm_proton_i, norm_neutron_i, norm_up_i, norm_down_i, norm_gluon_i, norm_electron_i, norm_photon_i, norm_other_i]
norm_particle_names_oi = ['Rho', 'St. K', 'Strange', 'Charm', 'Delta', 'Lambda', 'Epsilon', 'Eta', 'Omg/Phi', 'Difr', 'Diq (ud)', 'Diq(uu)']
norm_particle_counts_oi = [norm_rho_i, norm_strangeK_i, norm_strange_i, norm_charm_i, norm_delta_i, norm_lambda_i, norm_epsilon_i, norm_eta_i, norm_omegaphi_i, norm_difractive_i, norm_diq_i, norm_diq2_i]
norm_particle_names_un = ['Pions', 'Kaons', 'Protons', 'Neutrons', 'Up', 'Down', 'Gluon', 'Electron', 'Photon', 'Others']
norm_particle_counts_un = [norm_pion_un, norm_kaon_un, norm_proton_un, norm_neutron_un, norm_up_un, norm_down_un, norm_gluon_un, norm_electron_un, norm_photon_un, norm_other_un]
norm_particle_names_f = ['Pions', 'Kaons', 'Protons', 'Neutrons', 'Up', 'Down', 'Gluon', 'Electron', 'Photon', 'Others']
norm_particle_counts_f = [norm_pion_f, norm_kaon_f, norm_proton_f, norm_neutron_f, norm_up_f, norm_down_f, norm_gluon_f, norm_electron_f, norm_photon_f, norm_other_f]
norm_particle_names_o = ['Rho', 'St. K', 'Strange', 'Charm', 'Diquark', 'Delta', 'Lambda', 'Epsilon', 'Eta', 'Omg/Phi', 'Fragm', 'Xi']
norm_particle_counts_o = [norm_rho_un, norm_strangeK_un, norm_strange_un, norm_charm_un, norm_diq_un, norm_delta_un, norm_lambda_un, norm_epsilon_un, norm_eta_un, norm_omegaphi_un, norm_frag_un, norm_xi_un]
norm_particle_names_of = ['Rho', 'St. K', 'Strange', 'Charm', 'Diquark', 'Delta', 'Lambda', 'Epsilon', 'Eta', 'Omg/Phi', 'Fragm', 'Difr']
norm_particle_counts_of = [norm_rho_f, norm_strangeK_f, norm_strange_f, norm_charm_f, norm_diq_f, norm_delta_f, norm_lambda_f, norm_epsilon_f, norm_eta_f, norm_omegaphi_f, norm_frag_f, norm_difractive_f]

# INITIAL STATE
plt.subplot(3,2,1)
norm_palette_i = sns.color_palette("Blues", len(norm_particle_names_i))
plt.bar(norm_particle_names_i, norm_particle_counts_i, color=norm_palette_i, alpha=0.7, width=0.5)
plt.ylabel('Normalized Particle counts')
plt.title('Norm. Particles in the Initial State')
plt.ylabel('Normalized Particle counts')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(norm_particle_counts_i):
    plt.text(i, count + 0.05 * max(norm_particle_counts_i), f"{round(count, 2)}", ha='center', va='bottom', fontsize=8)
norm_total_particles_i = sum(norm_particle_counts_i)
plt.text(0.15, 0.9, f'Particles per event: {round(norm_total_particles_i,2)}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)

# OTHER INITIAL STATE
plt.subplot(3,2,2)
norm_palette_oi = sns.color_palette("Blues", len(norm_particle_names_oi))
plt.bar(norm_particle_names_oi, norm_particle_counts_oi, color=norm_palette_oi, alpha=0.7, width=0.5)
plt.ylabel('Normalized number of Particles')
plt.title('Norm. Other Initial Particles ')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(norm_particle_counts_oi):
    plt.text(i, count + 0.05 * max(norm_particle_counts_oi), f"{round(count, 2)}", ha='center', va='bottom', fontsize=8)
norm_total_particles_oi = sum(norm_particle_counts_oi)
plt.text(0.15, 0.9, f'Particles per event: {round(norm_total_particles_oi,2)}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)

# UNSTABLE STATE
plt.subplot(3,2,3)
norm_palette_un = sns.color_palette("Reds", len(norm_particle_names_un))
plt.bar(norm_particle_names_un, norm_particle_counts_un, color=norm_palette_un, alpha=0.7, width=0.5)
plt.ylabel('Normalized number of Particles')
plt.title('Norm. Unstable/Hadronized Particles')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(norm_particle_counts_un):
    plt.text(i, count + 0.05 * max(norm_particle_counts_un), f"{round(count, 2)}", ha='center', va='bottom', fontsize=8)
norm_total_particles_un = sum(norm_particle_counts_un)
plt.text(0.5, 0.9, f'Particles per event: {round(norm_total_particles_un,2)}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)

# OTHER UNSTABLE STATE
plt.subplot(3,2,4)
norm_palette_o = sns.color_palette("Reds", len(norm_particle_names_o))
plt.bar(norm_particle_names_o, norm_particle_counts_o, color=norm_palette_o, alpha=0.7, width=0.5)
plt.ylabel('Normalized number of Particles')
plt.title('Norm. Other Unstable Particles ')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(norm_particle_counts_o):
    plt.text(i, count + 0.05 * max(norm_particle_counts_o), f"{round(count, 2)}", ha='center', va='bottom', fontsize=8)
norm_total_particles_o = sum(norm_particle_counts_o)
plt.text(0.5, 0.9, f'Particles per event: {round(norm_total_particles_o,2)}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)

# FINAL STATE
plt.subplot(3,2,5)
norm_palette_f = sns.color_palette("Greens", len(norm_particle_names_f))
plt.bar(norm_particle_names_f, norm_particle_counts_f, color=norm_palette_f, alpha=0.7, width=0.5)
plt.ylabel('Normalized number of Particles')
plt.title('Norm. Particles in the Final State')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(norm_particle_counts_f):
    plt.text(i, count + 0.05 * max(norm_particle_counts_f), f"{round(count, 2)}", ha='center', va='bottom', fontsize=8)
norm_total_particles_f = sum(norm_particle_counts_f)
plt.text(0.5, 0.9, f'Particles per event: {round(norm_total_particles_f,2)}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)

#OTHER FINAL STATE
plt.subplot(3,2,6)
norm_palette_of = sns.color_palette("Greens", len(norm_particle_names_of))
plt.bar(norm_particle_names_of, norm_particle_counts_of, color=norm_palette_of, alpha=0.7, width=0.5)
plt.ylabel('Normalized number of Particles')
plt.title('Norm. Other Final Particles')
plt.tick_params(axis='both', which='major', labelsize=8)
for i, count in enumerate(norm_particle_counts_of):
    plt.text(i, count + 0.05 * max(norm_particle_counts_of), f"{round(count, 2)}", ha='center', va='bottom', fontsize=8)
norm_total_particles_of = sum(norm_particle_counts_of)
plt.text(0.5, 0.9, f'Particles per event: {round(norm_total_particles_of,2)}', ha='center', va='center', transform=plt.gca().transAxes, fontsize=8)

plt.tight_layout() 
'''
#________________________________________________________________ANGULAR DISTRIBUTION____________________________________________________________________________________________________

plt.figure(figsize=(14,8))
'''
# SCATTERED ELECTRON
plt.subplot(2,3,1, polar=True)
plt.scatter(scattered_el_ang_deg, scat_el_E, c=scat_el_E, cmap='viridis', alpha=0.7)
plt.colorbar(label='Energy', pad = 0.1)
plt.xlabel('Angles')
plt.title('Angular distribution of the electrons')
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
'''
# HADRONIZATION
plt.subplot(2,2,1)
categories = ["pion", "kaon", "proton", "neutron"]                  # NOT INTERESTED IN THE 'OTHERS' PLOT
counts_up = [count_daughters_up[name] for name in categories]
counts_down = [count_daughters_down[name] for name in categories]
counts_ud = [count_daughters_ud[name] for name in categories]
counts_uu = [count_daughters_uu[name] for name in categories]
# BAR POSITION
x = range(len(categories))
width = 0.2
# COLORS
palette_up = sns.color_palette("Blues")[1]
palette_down = sns.color_palette("Blues")[2]
palette_ud = sns.color_palette("Blues")[3]
palette_uu = sns.color_palette("Blues")[4]
# MULTIPLE BAR GRAPH
plt.bar(x, counts_up, width=width, label='Up Quark', color = palette_up)
plt.bar([i + width for i in x], counts_down, width=width, label='Down Quark', color = palette_down)
plt.bar([i + 2 * width for i in x], counts_ud, width=width, label='Diquark (ud)', color = palette_ud)
plt.bar([i + 3 * width for i in x], counts_uu, width=width, label='Diquark (uu)', color = palette_uu)
plt.xlabel('Categories')
plt.ylabel('Counts')
plt.title('Daughters from parton fragments')
plt.xticks([i + 1.5 * width for i in x], categories)
plt.legend()

# PION ANGOLAR DISTRIBUTIONS
plt.subplot(2,2,2, polar=True)
plt.scatter(pion_angles_deg, pion_energies, c=pion_energies, cmap='viridis', alpha=0.7, s=2)
plt.colorbar(label='Energy [GeV]', pad = 0.1)
plt.xlim(0, np.pi)
plt.title('Angular distribution of the pions')
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
plt.gca().set_rgrids([0, 50, 100, 150, 200, 250])

# KAON ANGOLAR DISTRIBUTIONS
plt.subplot(2,2,3, polar=True)
plt.scatter(kaon_angles_deg, kaon_energies, c=kaon_energies, cmap='viridis', alpha=0.7, s=2)
plt.colorbar(label='Energy [GeV]', pad = 0.1)
plt.xlim(0, np.pi)
plt.title('Angular distribution of the kaons')
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
plt.gca().set_rgrids([0, 50, 100, 150, 200, 250])

# PROTON ANGOLAR DISTRIBUTIONS
plt.subplot(2,2,4, polar=True)
plt.scatter(sc_electron_ang_deg, sc_electron_E, c=sc_electron_E, cmap='viridis', alpha=0.7, s=2)
plt.colorbar(label='Energy [GeV]', pad = 0.1)
plt.xlim(0, np.pi)
plt.title('Angular distribution of the protons')
plt.gca().set_theta_zero_location('E')
plt.gca().set_theta_direction(1)
plt.gca().set_rgrids([0, 3, 6, 9, 12, 15])

plt.tight_layout()

plt.show()


