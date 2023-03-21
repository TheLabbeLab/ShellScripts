import os, math

# ----------------------------------------------- MESS PARSING FUNCTIONS ----------------------------------------------------- #

# -------------------------------------------- Writes MESS Input Header -------------------------------------------------------#
def write_Header():
    MESSFile = open('MESSFile.inp', 'w')
    MESSFile.write('TemperatureList[K]                500. 625. 750. 875. 1000. 1125. 1250. 1375. 1500. 1625. 1750. 1875. 2000.\n')
    MESSFile.write('PressureList[atm]                 0.00001 0.001  0.01  0.1  1  10  100.\n')
    MESSFile.write('!PressureList[bar]                1.\n')
    MESSFile.write('EnergyStepOverTemperature         .2\n')
    MESSFile.write('ExcessEnergyOverTemperature       30\n')
    MESSFile.write('ModelEnergyLimit[kcal/mol]        400\n')
    MESSFile.write('CalculationMethod                 direct\n')
    MESSFile.write('!CalculationMethod                low-eigenvalue !direct\n')
    MESSFile.write('WellCutoff                        10\n')
    MESSFile.write('ChemicalEigenvalueMax             0.2\n')
    MESSFile.write('Model\n  EnergyRelaxation\n    Exponential\n      Factor[1/cm]'+' '*16+'200\n')
    MESSFile.write(' '*6+'Power                       .85\n')
    MESSFile.write(' ' *6+'ExponentCutoff'+'              15\n'+' '*4+'End\n')
    MESSFile.write('  CollisionFrequency\n    LennardJones\n      Epsilons[1/cm]     --FILL HERE-- --FILL HERE--  !N2 and parent\n')
    MESSFile.write('  Sigmas[angstrom]       --FILL HERE-- --FILL HERE--\n  Masses[amu]            --FILL HERE-- --FILL HERE--\n    End\n')
    MESSFile.close()

# --------------------------------- CREATES AND WRITES THE BARRIER SECTION FOR A WELL  --------------------------------------- #
# Assumes all wells are connected to Well 1!
def write_Well_Barrier(TSW_xyz_geom, TSW_freqs, TSW_charge_multiplicity, TSW_rots, rotorEnergies, path_dir, m062xFileName):
    MESSFile = open('MESSFile.inp', 'a')
    # gets path folder name using os library
    pathway = os.path.basename(os.path.normpath(path_dir))
    # extracts the number from the pathway (ie. returns 8 for the P8 folder)
    path_num = ''.join([item for item in pathway if item.isdigit()])
    MESSFile.write('Barrier'+' '*8+'B'+str(int(path_num)-1)+' '*2+'W1'+' '*2+pathway+' '*4+'# '+m062xFileName+'\n')
    MESSFile.write(' '*4+'Variational\n')
    MESSFile.write(' '*8+'RRHO\n')
    MESSFile.write(' '*4+'Geometry[angstrom]'+' '*8+str(len(TSW_xyz_geom))+'\n')
    # writes geometry section
    for line in TSW_xyz_geom:
        line = line.split(',')
        MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
    MESSFile.write(' '*4+'Core RigidRotor\n')
    MESSFile.write(' '*8+'SymmetryFactor'+' '*4+'--- INCLUDE SYMMETRY FACTOR HERE ----\n')
    MESSFile.write(' '*4+'End\n')
    # writes rotor section
    for i in range(len(TSW_rots)):
        max_atom_number = max(TSW_rots[i])
        rotEnergies = rotorEnergies['D'+str(max_atom_number-3)]
        MESSFile.write(' '*4+'Rotor'+' '*8+'Hindered\n')
        MESSFile.write(' '*8+'Group'+' '*21+'numbers!\n')
        MESSFile.write(' '*8+'Axis'+' '*22+str(TSW_rots[i][1])+' '+str(TSW_rots[i][2])+'\n')
        # calculates rotor symmetry
        truncated_rotEnergies = [math.floor(energy*10**2)/10**2 for energy in rotEnergies]
        rotorSymmetry = '1'
        if len(truncated_rotEnergies) == 4:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*2:
                rotorSymmetry = '2'
        if len(truncated_rotEnergies) == 6:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*3:
                rotorSymmetry = '3'
        if len(truncated_rotEnergies) == 8:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*4:
                rotorSymmetry = '4'
        MESSFile.write(' '*8+'Symmetry'+' '*19+rotorSymmetry+'\n')
        rotEnergies = [round(energy, 2) for energy in rotEnergies]
        MESSFile.write(' '*8+'Potential[kcal/mol]'+' '*8+str(int(len(rotEnergies)/int(rotorSymmetry)))+'\n')
        MESSFile.write(' '*9)
        if not rotorSymmetry == '1' :
            for i in range(len(rotEnergies)):
                MESSFile.write(str(rotEnergies[i])+'  ')
                if i == 1:
                    break
        else:
            for rotEnergy in rotEnergies:
                MESSFile.write(str(rotEnergy)+'  ')
        MESSFile.write('\n')
        MESSFile.write(' '*4+'End\n')
    # writes frequecies section substracts number of rotors + 1 imaginary frequency!
    MESSFile.write(' '*4+'Frequencies[1/cm]'+' '*4+str(len(TSW_freqs)-(len(TSW_rots)+1))+'\n     ')
    counter = 1
    for freq in TSW_freqs[1:]:
        MESSFile.write(' '*5)
        MESSFile.write(' '+str(round(freq, 2)))
        if counter == 10:
            MESSFile.write('\n     ')
            counter = 0
        counter += 1
    MESSFile.write('\n'+' '*5+'!------- INCLUDE ROTOR FREQUENCIES HERE ---------!')
    MESSFile.write('\n'+' '*4+'ZeroEnergy[kcal/mol]'+' '*4+'--- ZERO ENERGY VALUE FROM HIGH ENERGY CALCULATIONS ---')
    # change for higher excited states / pulls charge and multplicity from m062x file
    MESSFile.write('\n'+' '*4+'ElectronicLevels[1/cm]'+' '*4+str(1)+'\n'+' '*8+TSW_charge_multiplicity[0]+' '*2+TSW_charge_multiplicity[1]+'\n'+' '*4+'End\n')
    MESSFile.write(' '*4+'Tunneling'+' '*16+'Eckart\n')
    MESSFile.write(' '*4+'ImaginaryFrequency[1/cm]'+' '*4+str(TSW_freqs[0])+'\n') # first freq is list the imaginary frequency!
    MESSFile.write(' '*4+'WellDepth[kcal/mol]'+' '*4+'--- FORWARD BARRIER OF TS ---\n')
    MESSFile.write(' '*4+'WellDepth[kcal/mol]'+' '*4+'--- BACKWARD BARRIER OF TS ---\n'+' '*4+'End\n'+'End\n')
    MESSFile.close()

# ------------------------------------------ CREATES AND WRITES A WELL SECTION ----------------------------------------------- #
# Place well's m062x file inside the coresponding W# directory!
# this function pulls pathway, m062xFileName from script
def write_Well(W_xyz_geom, W_freqs, W_charge_multiplicity, W_rots, rotorEnergies, path_dir, m062xFileName):
    MESSFile = open('MESSFile.inp', 'a')
    # gets path folder name using os library
    pathway = os.path.basename(os.path.normpath(path_dir))
    MESSFile.write('Well'+' '*8+pathway+' '*8+'# '+m062xFileName+'\n')
    MESSFile.write(' '*4+'Species\n')
    MESSFile.write(' '*8+'RRHO\n')
    MESSFile.write(' '*4+'Geometry[angstrom]'+' '*8+str(len(W_xyz_geom))+'\n')
    # writes geometry section
    for line in W_xyz_geom:
        line = line.split(',')
        MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
    MESSFile.write(' '*4+'Core RigidRotor\n')
    MESSFile.write(' '*8+'SymmetryFactor'+' '*4+'--- INCLUDE SYMMETRY FACTOR HERE ----\n')
    MESSFile.write(' '*4+'End\n')
    # writes rotor section
    for i in range(len(W_rots)):
        max_atom_number = max(W_rots[i])
        rotEnergies = rotorEnergies['D'+str(max_atom_number-3)]
        MESSFile.write(' '*4+'Rotor'+' '*8+'Hindered\n')
        MESSFile.write(' '*8+'Group'+' '*21+'numbers!\n')
        MESSFile.write(' '*8+'Axis'+' '*22+str(W_rots[i][1])+' '+str(W_rots[i][2])+'\n')
        # calculates rotor symmetry
        truncated_rotEnergies = [math.floor(energy*10**2)/10**2 for energy in rotEnergies]
        rotorSymmetry = '1'
        if len(truncated_rotEnergies) == 4:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*2:
                rotorSymmetry = '2'
        if len(truncated_rotEnergies) == 6:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*3:
                rotorSymmetry = '3'
        if len(truncated_rotEnergies) == 8:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*4:
                rotorSymmetry = '4'
        MESSFile.write(' '*8+'Symmetry'+' '*19+rotorSymmetry+'\n')
        rotEnergies = [round(energy, 2) for energy in rotEnergies]
        MESSFile.write(' '*8+'Potential[kcal/mol]'+' '*8+str(int(len(rotEnergies)/int(rotorSymmetry)))+'\n')
        MESSFile.write(' '*9)
        if not rotorSymmetry == '1' :
            for i in range(len(rotEnergies)):
                MESSFile.write(str(rotEnergies[i])+'  ')
                if i == 1:
                    break
        else:
            for rotEnergy in rotEnergies:
                MESSFile.write(str(rotEnergy)+'  ')
        MESSFile.write('\n')
        MESSFile.write(' '*4+'End\n')
    # writes frequecies section substracts number of rotors
    MESSFile.write(' '*4+'Frequencies[1/cm]'+' '*4+str(len(W_freqs)-(len(W_rots)))+'\n     ')
    counter = 1
    for freq in W_freqs:
        MESSFile.write(' '*5)
        MESSFile.write(' '+str(round(freq, 2)))
        if counter == 10:
            MESSFile.write('\n     ')
            counter = 0
        counter += 1
    MESSFile.write('\n'+' '*5+'!------- INCLUDE ROTOR FREQUENCIES HERE ---------!')
    MESSFile.write('\n'+' '*4+'ZeroEnergy[kcal/mol]'+' '*4+'--- ZERO ENERGY VALUE FROM HIGH ENERGY CALCULATIONS ---')
    # change for higher excited states / pulls charge and multplicity from m062x file
    MESSFile.write('\n'+' '*4+'ElectronicLevels[1/cm]'+' '*4+str(1)+'\n'+' '*8+W_charge_multiplicity[0]+' '*2+W_charge_multiplicity[1]+'\n'+' '*4+'End\n'+'End\n')
    MESSFile.close()

# ---------------------------------- CREATES AND WRITES THE BARRIER SECTION FOR A TS  ---------------------------------------- #
def write_TS_Barrier(xyz_geom, rots, freqs, rotorEnergies, charge_multiplicity, path_dir, m062xFileName, number_of_wells, h_abstraction):
    MESSFile = open('MESSFile.inp', 'a')
    # gets path folder name using os library
    pathway = os.path.basename(os.path.normpath(path_dir))
    # extracts the number from the pathway (ie. returns 8 for the P8 folder)
    path_num = ''.join([item for item in pathway if item.isdigit()])
    # ------------------------ Naming assumes Barrier connects W1 to products!(ie. B# W1 P#) ----------------------------------#
    if h_abstraction:
        MESSFile.write('Barrier'+' '*8+'B'+str(int(path_num)+number_of_wells)+' '*2+'W1'+' '*2+pathway+' '*4+'# '+m062xFileName+'  +  [H]'+'\n')
    else:
        MESSFile.write('Barrier'+' '*8+'B'+str(int(path_num)+number_of_wells)+' '*2+'W1'+' '*2+pathway+' '*4+'# '+m062xFileName+'\n')
    MESSFile.write(' '*4+'Variational\n')
    MESSFile.write(' '*8+'RRHO\n')
    MESSFile.write(' '*4+'Geometry[angstrom]'+' '*8+str(len(xyz_geom))+'\n')
    # writes geometry section
    for line in xyz_geom:
        line = line.split(',')
        MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
    MESSFile.write(' '*4+'Core RigidRotor\n')
    MESSFile.write(' '*8+'SymmetryFactor'+' '*4+'--- INCLUDE SYMMETRY FACTOR HERE ----\n')
    MESSFile.write(' '*4+'End\n')
    # writes rotor section
    for i in range(len(rots)):
        max_atom_number = max(rots[i])
        rotEnergies = rotorEnergies['D'+str(max_atom_number-3)]
        MESSFile.write(' '*4+'Rotor'+' '*8+'Hindered\n')
        MESSFile.write(' '*8+'Group'+' '*21+'numbers!\n')
        MESSFile.write(' '*8+'Axis'+' '*22+str(rots[i][1])+' '+str(rots[i][2])+'\n')
        # calculates rotor symmetry
        truncated_rotEnergies = [math.floor(energy*10**2)/10**2 for energy in rotEnergies]
        rotorSymmetry = '1'
        if len(truncated_rotEnergies) == 4:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*2:
                rotorSymmetry = '2'
        if len(truncated_rotEnergies) == 6:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*3:
                rotorSymmetry = '3'
        if len(truncated_rotEnergies) == 8:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*4:
                rotorSymmetry = '4'
        MESSFile.write(' '*8+'Symmetry'+' '*19+rotorSymmetry+'\n')
        rotEnergies = [round(energy, 2) for energy in rotEnergies]
        MESSFile.write(' '*8+'Potential[kcal/mol]'+' '*8+str(int(len(rotEnergies)/int(rotorSymmetry)))+'\n')
        MESSFile.write(' '*9)
        if not rotorSymmetry == '1' :
            for i in range(len(rotEnergies)):
                MESSFile.write(str(rotEnergies[i])+'  ')
                if i == 1:
                    break
        else:
            for rotEnergy in rotEnergies:
                MESSFile.write(str(rotEnergy)+'  ')
        MESSFile.write('\n')
        MESSFile.write(' '*4+'End\n')
    # writes frequecies section substracts number of rotors + 1 imaginary frequency!
    MESSFile.write(' '*4+'Frequencies[1/cm]'+' '*4+str(len(freqs)-(len(rots)+1))+'\n     ')
    counter = 1
    for freq in freqs[1:]:
        MESSFile.write(' '*5)
        MESSFile.write(' '+str(round(freq, 2)))
        if counter == 10:
            MESSFile.write('\n     ')
            counter = 0
        counter += 1
    MESSFile.write('\n'+' '*5+'!------- INCLUDE ROTOR FREQUENCIES HERE ---------!')
    MESSFile.write('\n'+' '*4+'ZeroEnergy[kcal/mol]'+' '*4+'--- ZERO ENERGY VALUE FROM HIGH ENERGY CALCULATIONS ---')
    # change for higher excited states / pulls charge and multplicity from m062x file
    MESSFile.write('\n'+' '*4+'ElectronicLevels[1/cm]'+' '*4+str(1)+'\n'+' '*8+charge_multiplicity[0]+' '*2+charge_multiplicity[1]+'\n'+' '*4+'End\n')
    MESSFile.write(' '*4+'Tunneling'+' '*16+'Eckart\n')
    MESSFile.write(' '*4+'ImaginaryFrequency[1/cm]'+' '*4+str(freqs[0])+'\n') # first freq is list the imaginary frequency!
    MESSFile.write(' '*4+'WellDepth[kcal/mol]'+' '*4+'--- FORWARD BARRIER OF TS ---\n')
    MESSFile.write(' '*4+'WellDepth[kcal/mol]'+' '*4+'--- BACKWARD BARRIER OF TS ---\n'+' '*4+'End\n'+'End\n')
    MESSFile.close()
# ------------------------------ WRITES THE BIMOLECULAR SECTION FOR A TS  ---------------------------------------- #
# THIS CODE ASSUMES THE TS CONNECTS TO TWO FRAGMENTS!!! -> CHANGE IF NOT!
# TS log must end with _m062x.log: (ie. TS_CH2O_CCCOCCC_m062x.log)
# Fragment folder that stores must have the name of the TS: (ie. TS_CH2O_CCCOCCC)
# Fragment logs inside the folder must end with _m062x.log: (ie. CH2O_m062x.log and CCCOCCC_m062x.log)
# fragments folder must be named 'fragments'
# ----------------- Access each fragment's log file and extract geometries, rotors and frequencies --------------------------- #

# This function takes in root and fragment's directory paths from from script
def write_TS_Bimolecular(iteration, fragments_naming, fragmentFile, fragments_directory, path_dir, m062xFragmentFileList, fr_xyz_geom, fr_freqs, fr_rots, rotorEnergies, fr_charge_multiplicity, h_abstraction):
    MESSFile = open('MESSFile.inp', 'a')
    # gets path folder name using os library
    pathway = os.path.basename(os.path.normpath(path_dir))
    # makes it so Bimolecular header is writen only once!
    if iteration == 1:
        if h_abstraction:
            MESSFile.write('Bimolecular'+' '*8+pathway+' '*4+'# '+fragments_naming+'  +  [H]'+'\n')
        else:
            MESSFile.write('Bimolecular'+' '*8+pathway+' '*4+'# '+fragments_naming+'\n')
    MESSFile.write(' '*4+'Fragment'+' '*2+fragmentFile.split('_m062x.log')[0]+'\n')
    MESSFile.write(' '*8+'RRHO\n')
    MESSFile.write(' '*4+'Geometry[angstrom]'+' '*8+str(len(fr_xyz_geom))+'\n')
    # writes geometry section
    for line in fr_xyz_geom:
        line = line.split(',')
        MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
    MESSFile.write(' '*4+'Core RigidRotor\n')
    MESSFile.write(' '*8+'SymmetryFactor'+' '*4+'--- INCLUDE SYMMETRY FACTOR HERE ----\n')
    MESSFile.write(' '*4+'End\n')
    # writes rotor section
    for i in range(len(fr_rots)):
        max_atom_number = max(fr_rots[i])
        rotEnergies = rotorEnergies['D'+str(max_atom_number-3)]
        MESSFile.write(' '*4+'Rotor'+' '*8+'Hindered\n')
        MESSFile.write(' '*8+'Group'+' '*21+'numbers!\n')
        MESSFile.write(' '*8+'Axis'+' '*22+str(fr_rots[i][1])+' '+str(fr_rots[i][2])+'\n')
        # calculates rotor symmetry
        truncated_rotEnergies = [math.floor(energy*10**2)/10**2 for energy in rotEnergies]
        rotorSymmetry = '1'
        if len(truncated_rotEnergies) == 4:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*2:
                rotorSymmetry = '2'
        if len(truncated_rotEnergies) == 6:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*3:
                rotorSymmetry = '3'
        if len(truncated_rotEnergies) == 8:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*4:
                rotorSymmetry = '4'
        MESSFile.write(' '*8+'Symmetry'+' '*19+rotorSymmetry+'\n')
        rotEnergies = [round(energy, 2) for energy in rotEnergies]
        MESSFile.write(' '*8+'Potential[kcal/mol]'+' '*8+str(int(len(rotEnergies)/int(rotorSymmetry)))+'\n')
        MESSFile.write(' '*9)
        if not rotorSymmetry == '1' :
            for i in range(len(rotEnergies)):
                MESSFile.write(str(rotEnergies[i])+'  ')
                if i == 1:
                    break
        else:
            for rotEnergy in rotEnergies:
                MESSFile.write(str(rotEnergy)+'  ')
        MESSFile.write('\n')
        MESSFile.write(' '*4+'End\n')
    # writes frequecies section substracts number of rotors!
    MESSFile.write(' '*4+'Frequencies[1/cm]'+' '*4+str(len(fr_freqs)-(len(fr_rots)))+'\n     ')
    counter = 1
    for freq in fr_freqs:
        MESSFile.write(' '*5)
        MESSFile.write(' '+str(round(freq, 2)))
        if counter == 10:
            MESSFile.write('\n     ')
            counter = 0
        counter += 1
    MESSFile.write('\n'+' '*5+'!------- INCLUDE ROTOR FREQUENCIES HERE ---------!')
    MESSFile.write('\n'+' '*4+'ZeroEnergy[kcal/mol]'+' '*4+'0')
    # change for higher excited states / pulls charge and multplicity from m062x file
    MESSFile.write('\n'+' '*4+'ElectronicLevels[1/cm]'+' '*4+str(1)+'\n'+' '*8+fr_charge_multiplicity[0]+' '*2+fr_charge_multiplicity[1]+'\n'+' '*4+'End\n')
    if h_abstraction:
        MESSFile.write(' '*4+'Fragment    H\n'+' '*6+'Atom\n'+' '*8+'Mass[amu]    1\n'+' '*8+'ElectronicLevels[1/cm]    1\n'+' '*6+'    0  2'+'\n    End\n')
    # writes ground energy at the end of the section
    if iteration == len(m062xFragmentFileList):
        MESSFile.write(' '*4+'GroundEnergy[kcal/mol]'+' '*8+' ----Ground Energy Value from High Energy Calcs -----'+'\n'+'End\n')
        MESSFile.close()
    # goes back to fragment directory
    os.chdir(fragments_directory)


# ------------------------------ WRITES THE BARRIER SECTION FOR A BF  ---------------------------------------- #
# Place two fragment m062x files in the pathway folder (ie. P3)
# Fragment files must begin with Bi_... and end with _m062x.log (ie. Bi_COH_E9-27_m062x.log)
# ----------------- Access each fragment's log file and extract geometries, rotors and frequencies --------------------------- #

# fragment#_info contains [xyz_geom, freqs, charge/multiplicity and rotors]
def write_BF_Barrier(fragments_naming, fragment1_info, fr1_rotorEnergies, fragment2_info, fr2_rotorEnergies,  pathway, number_of_wells, m062xFileName, h_abstraction):
    MESSFile = open('MESSFile.inp', 'a')
    # extracts the number from the pathway (ie. returns 8 for the P8 folder)
    path_num = ''.join([item for item in pathway if item.isdigit()])
    # ------------------------ Naming assumes Barrier connects W1 to products!(ie. B# W1 P#) ----------------------------------#
    if h_abstraction:
        MESSFile.write('Barrier'+' '*8+'B'+str(int(path_num)+number_of_wells)+' '*2+'W1'+' '*2+pathway+' '*4+'# '+m062xFileName+'  +  [H]'+'\n')
    else:
        MESSFile.write('Barrier'+' '*8+'B'+str(int(path_num)+number_of_wells)+' '*2+'W1'+' '*2+pathway+' '*4+'# '+fragments_naming+'\n')
    MESSFile.write(' '*4+'RRHO\n')
    MESSFile.write(' '*8+'Stoichiometry'+' '*8+'------ADD MOLECULAR FORMULA------\n')
    MESSFile.write(' '*8+'Core'+' '*4+'PhaseSpaceTheory\n')
    MESSFile.write(' '*4+'FragmentGeometry[angstrom]'+' '*8+str(len(fragment1_info[0]))+'\n')
    # writes geometry section for fragment 1
    for line in fragment1_info[0]:
        line = line.split(',')
        MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
    # writes geometry section for fragment 2
    if h_abstraction:
        MESSFile.write(' '*4+'FragmentGeometry[angstrom]        1\n'+' '*5+'H                   0.000000    0.000000    0.000000\n')
    else:
        MESSFile.write(' '*4+'FragmentGeometry[angstrom]'+' '*8+str(len(fragment2_info[0]))+'\n')
        for line in fragment2_info[0]:
            line = line.split(',')
            MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
    MESSFile.write(' '*4+'SymmetryFactor'+' '*4+'--- INCLUDE SYMMETRY FACTOR HERE ----\n')
    MESSFile.write(' '*8+'PotentialPrefactor[au]'+' '*4+' ----------- ADD POTENTIAL PREFACTOR HERE! ----------\n')
    MESSFile.write(' '*8+'PotentialPowerExponent'+' '*4+'6\n')
    MESSFile.write(' '*4+'End\n')
    # writes rotor section for fragment 1
    for i in range(len(fragment1_info[3])):
        max_atom_number = max(fragment1_info[3][i])
        fr1_rotEnergies = fr1_rotorEnergies['D'+str(max_atom_number-3)]
        MESSFile.write(' '*4+'Rotor'+' '*8+'Hindered\n')
        MESSFile.write(' '*4+'Geometry[angstrom]'+' '*8+str(len(fragment1_info[0]))+'\n')
        for line in fragment1_info[0]:
            line = line.split(',')
            MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
        MESSFile.write(' '*8+'Group'+' '*21+'numbers!\n')
        MESSFile.write(' '*8+'Axis'+' '*22+str(fragment1_info[3][i][1])+' '+str(fragment1_info[3][i][2])+'\n')
        # calculates rotor symmetry
        fr1_truncated_rotEnergies = [math.floor(energy*10**2)/10**2 for energy in fr1_rotEnergies]
        rotorSymmetry = '1'
        if len(fr1_truncated_rotEnergies) == 4:
            if fr1_truncated_rotEnergies == fr1_truncated_rotEnergies[:2]*2:
                rotorSymmetry = '2'
        if len(fr1_truncated_rotEnergies) == 6:
            if fr1_truncated_rotEnergies == fr1_truncated_rotEnergies[:2]*3:
                rotorSymmetry = '3'
        if len(fr1_truncated_rotEnergies) == 8:
            if fr1_truncated_rotEnergies == fr1_truncated_rotEnergies[:2]*4:
                rotorSymmetry = '4'
        MESSFile.write(' '*8+'Symmetry'+' '*19+rotorSymmetry+'\n')
        fr1_rotEnergies = [round(energy, 2) for energy in fr1_rotEnergies]
        MESSFile.write(' '*8+'Potential[kcal/mol]'+' '*8+str(int(len(fr1_rotEnergies)/int(rotorSymmetry)))+'\n')
        MESSFile.write(' '*9)
        if not rotorSymmetry == '1' :
            for i in range(len(fr1_rotEnergies)):
                MESSFile.write(str(fr1_rotEnergies[i])+'  ')
                if i == 1:
                    break
        else:
            for rotEnergy in fr1_rotEnergies:
                MESSFile.write(str(rotEnergy)+'  ')
        MESSFile.write('\n')
        MESSFile.write(' '*4+'End\n')
    if not h_abstraction:
        # writes rotor section for fragment 2
        for i in range(len(fragment2_info[3])):
            max_atom_number = max(fragment2_info[3][i])
            fr2_rotEnergies = fr2_rotorEnergies['D'+str(max_atom_number-3)]
            MESSFile.write(' '*4+'Rotor'+' '*8+'Hindered\n')
            MESSFile.write(' '*4+'Geometry[angstrom]'+' '*8+str(len(fragment2_info[0]))+'\n')
            for line in fragment2_info[0]:
                line = line.split(',')
                MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
            MESSFile.write(' '*8+'Group'+' '*21+'numbers!\n')
            MESSFile.write(' '*8+'Axis'+' '*22+str(fragment2_info[3][i][1])+' '+str(fragment2_info[3][i][2])+'\n')
            # calculates rotor symmetry
            fr2_truncated_rotEnergies = [math.floor(energy*10**2)/10**2 for energy in fr2_rotEnergies]
            rotorSymmetry = '1'
            if len(fr2_truncated_rotEnergies) == 4:
                if fr2_truncated_rotEnergies == fr2_truncated_rotEnergies[:2]*2:
                    rotorSymmetry = '2'
            if len(fr2_truncated_rotEnergies) == 6:
                if fr2_truncated_rotEnergies == fr2_truncated_rotEnergies[:2]*3:
                    rotorSymmetry = '3'
            if len(fr2_truncated_rotEnergies) == 8:
                if fr2_truncated_rotEnergies == fr2_truncated_rotEnergies[:2]*4:
                    rotorSymmetry = '4'
            MESSFile.write(' '*8+'Symmetry'+' '*19+rotorSymmetry+'\n')
            fr2_rotEnergies = [round(energy, 2) for energy in fr2_rotEnergies]
            MESSFile.write(' '*8+'Potential[kcal/mol]'+' '*8+str(int(len(fr2_rotEnergies)/int(rotorSymmetry)))+'\n')
            MESSFile.write(' '*9)
            if not rotorSymmetry == '1' :
                for i in range(len(fr2_rotEnergies)):
                    MESSFile.write(str(fr2_rotEnergies[i])+'  ')
                    if i == 1:
                        break
            else:
                for rotEnergy in fr2_rotEnergies:
                    MESSFile.write(str(rotEnergy)+'  ')
            MESSFile.write('\n')
            MESSFile.write(' '*4+'End\n')
        # sums total number of frequencies and extracts number of rotors for both fragments
        total_frequencies = len(fragment1_info[1])+len(fragment2_info[1])-len(fragment1_info[3])-len(fragment2_info[3])
    else:
        total_frequencies = len(fragment1_info[1])-len(fragment1_info[3])
    MESSFile.write(' '*4+'Frequencies[1/cm]'+' '*4+str(total_frequencies)+'\n     ')
    counter = 1
    for freq in fragment1_info[1]:
        MESSFile.write(' '*5)
        MESSFile.write(' '+str(round(freq, 2)))
        if counter == 10:
            MESSFile.write('\n     ')
            counter = 0
        counter += 1
    MESSFile.write('\n     ')
    counter = 1
    if not h_abstraction:
        for freq in fragment2_info[1]:
            MESSFile.write(' '*5)
            MESSFile.write(' '+str(round(freq, 2)))
            if counter == 10:
                MESSFile.write('\n     ')
                counter = 0
            counter += 1
    MESSFile.write('\n'+' '*5+'!------- INCLUDE ROTOR FREQUENCIES HERE ---------!')
    MESSFile.write('\n'+' '*4+'ZeroEnergy[kcal/mol]'+' '*4+'--- ZERO ENERGY VALUE FROM HIGH ENERGY CALCULATIONS ---')
    # change for higher excited states / pulls charge and multplicity from m062x file
    # for a BF, both fragments have charge and multplicity 0 2
    MESSFile.write('\n'+' '*4+'ElectronicLevels[1/cm]'+' '*4+str(1)+'\n'+' '*8+fragment1_info[2][0]+' '*2+fragment1_info[2][1]+'\n'+' '*4+'End\n')
    MESSFile.close()

# ------------------------------ WRITES THE BIMOLECULAR SECTION FOR A BF  ---------------------------------------- #
# Place two fragment m062x files in the pathway folder (ie. P3)
# Fragment files must begin with Bi_... and end with _m062x.log (ie. Bi_COH_E9-27_m062x.log)
# ----------------- Access each fragment's log file and extract geometries, rotors and frequencies --------------------------- #

# This function takes in root pathway's directory paths from from script
def write_BF_Bimolecular(iteration, root_directory, fragments_naming, fragmentFile, path_dir, m062xFileList, fr_xyz_geom, fr_freqs, fr_rots, rotorEnergies, fr_charge_multiplicity, h_abstraction):
    MESSFile = open('MESSFile.inp', 'a')
    # gets path folder name using os library
    pathway = os.path.basename(os.path.normpath(path_dir))
    # makes it so Bimolecular header is writen only once!
    if iteration == 1:
        if h_abstraction:
            MESSFile.write('Bimolecular'+' '*8+pathway+' '*4+'# '+fragmentFile.split('_m062x.log')[0]+'  +  [H]'+'\n')
        else:
            MESSFile.write('Bimolecular'+' '*8+pathway+' '*4+'# '+fragments_naming+'\n')
    MESSFile.write(' '*4+'Fragment'+' '*2+fragmentFile.split('_m062x.log')[0]+'\n')
    MESSFile.write(' '*8+'RRHO\n')
    MESSFile.write(' '*4+'Geometry[angstrom]'+' '*8+str(len(fr_xyz_geom))+'\n')
    # writes geometry section
    for line in fr_xyz_geom:
        line = line.split(',')
        MESSFile.write(' '*5+line[0]+' '*8+line[1]+' '*4+line[2]+' '*4+line[3]+'\n')
    MESSFile.write(' '*4+'Core RigidRotor\n')
    MESSFile.write(' '*8+'SymmetryFactor'+' '*4+'--- INCLUDE SYMMETRY FACTOR HERE ----\n')
    MESSFile.write(' '*4+'End\n')
    # writes rotor section
    for i in range(len(fr_rots)):
        max_atom_number = max(fr_rots[i])
        rotEnergies = rotorEnergies['D'+str(max_atom_number-3)]
        MESSFile.write(' '*4+'Rotor'+' '*8+'Hindered\n')
        MESSFile.write(' '*8+'Group'+' '*21+'numbers!\n')
        MESSFile.write(' '*8+'Axis'+' '*22+str(fr_rots[i][1])+' '+str(fr_rots[i][2])+'\n')
        # calculates rotor symmetry
        truncated_rotEnergies = [math.floor(energy*10**2)/10**2 for energy in rotEnergies]
        rotorSymmetry = '1'
        if len(truncated_rotEnergies) == 4:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*2:
                rotorSymmetry = '2'
        if len(truncated_rotEnergies) == 6:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*3:
                rotorSymmetry = '3'
        if len(truncated_rotEnergies) == 8:
            if truncated_rotEnergies == truncated_rotEnergies[:2]*4:
                rotorSymmetry = '4'
        MESSFile.write(' '*8+'Symmetry'+' '*19+rotorSymmetry+'\n')
        rotEnergies = [round(energy, 2) for energy in rotEnergies]
        MESSFile.write(' '*8+'Potential[kcal/mol]'+' '*8+str(int(len(rotEnergies)/int(rotorSymmetry)))+'\n')
        MESSFile.write(' '*9)
        if not rotorSymmetry == '1' :
            for i in range(len(rotEnergies)):
                MESSFile.write(str(rotEnergies[i])+'  ')
                if i == 1:
                    break
        else:
            for rotEnergy in rotEnergies:
                MESSFile.write(str(rotEnergy)+'  ')
        MESSFile.write('\n')
        MESSFile.write(' '*4+'End\n')
    # writes frequecies section substracts number of rotors!
    MESSFile.write(' '*4+'Frequencies[1/cm]'+' '*4+str(len(fr_freqs)-(len(fr_rots)))+'\n     ')
    counter = 1
    for freq in fr_freqs:
        MESSFile.write(' '*5)
        MESSFile.write(' '+str(round(freq, 2)))
        if counter == 10:
            MESSFile.write('\n     ')
            counter = 0
        counter += 1
    MESSFile.write('\n'+' '*5+'!------- INCLUDE ROTOR FREQUENCIES HERE ---------!')
    MESSFile.write('\n'+' '*4+'ZeroEnergy[kcal/mol]'+' '*4+'--- ZERO ENERGY VALUE FROM HIGH ENERGY CALCULATIONS ---')
    # change for higher excited states / pulls charge and multplicity from m062x file
    MESSFile.write('\n'+' '*4+'ElectronicLevels[1/cm]'+' '*4+str(1)+'\n'+' '*8+fr_charge_multiplicity[0]+' '*2+fr_charge_multiplicity[1]+'\n'+' '*4+'End\n')
    # Adds H atom fragment section
    if h_abstraction:
        MESSFile.write(' '*4+'Fragment    H\n'+' '*6+'Atom\n'+' '*8+'Mass[amu]    1\n'+' '*8+'ElectronicLevels[1/cm]    1\n'+' '*6+'    0  2'+'\n    End\n')
    # writes ground energy at the end of the section
    if iteration == len(m062xFileList):
        MESSFile.write(' '*4+'GroundEnergy[kcal/mol]'+' '*8+' ----Ground Energy Value from High Energy Calcs -----'+'\n'+'End\n')
        MESSFile.close()
    # goes back to BF pathway directory
    os.chdir(root_directory+'\\'+pathway)
