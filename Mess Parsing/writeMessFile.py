from messParsing import *
from geometryInfo import *
from get_energies import *
from os.path import basename

# Energies Dictionary, will contain all energies for each molecule to be added into energies spredsheet.
energiesDictionary = {}
all_molecules = []

# Condition to write excel spreadsheet
writeExcel = str(input("Do you want to create an Excel spreadsheet? (y/n): "))
if writeExcel == 'y':
    createExcel = True
else:
    createExcel = False

## -- Adds energies from m062x and high energy calcs into the energiesDictionary -- ##
def add_energies(pathway, FileList, non_ts_prods, non_ts, energiesDictionary, createExcel=False):
    num_fragments = len([file for file in FileList if file.endswith('_m062x.log') and not file.startswith('TS') and pathway.startswith('P')])
    FileList = sorted(FileList)
    for filename in FileList:
        file = open(filename, 'r')
        filestring = file.read().replace('\n', '').replace(' ', '')
        file.close()
        if filename.endswith('E.log'):
            # gets energy from CCSD(T)
            E_ccsdt = filestring.split('\CCSD(T)')[1].split('\\')[0]
            E_ccsdt = float(re.findall(r'-?\d*\.\d*', E_ccsdt)[0])
            # gets T1 diagnostic
            if 'T1Diagnostic' in filestring:
                t1diag = filestring.split('T1Diagnostic')[-1]
                t1diag = float(re.findall(r'-?\d*\.\d*', t1diag)[0])
        if filename.endswith('E1.log'):
            E1_ccsdt = filestring.split('\CCSD(T)')[1].split('\\')[0]
            E1_ccsdt = float(re.findall(r'-?\d*\.\d*', E1_ccsdt)[0])
            if 'T1Diagnostic' in filestring:
                t1diag = filestring.split('T1Diagnostic')[-1]
                t1diag = float(re.findall(r'-?\d*\.\d*', t1diag)[0])
        if filename.endswith('m062x.log'):
            zeroPoint = filestring.split('\ZeroPoint')[-1]
            zeroPoint = float(re.findall(r'-?\d*\.\d*', zeroPoint)[0])
            molecule  = filename.split('_m062x.log')[0]
            if createExcel:
                energiesDictionary[molecule] = [zeroPoint, E_ccsdt, E1_ccsdt, t1diag]
                # this condition catches single fragments (assumed to be accompanied by H)
                if num_fragments == 1:
                    energiesDictionary['H'] = [0, -0.499278, -0.499809811, 0]
                if not molecule.startswith('TS') and basename(os.getcwd()) != 'fragments':
                    non_ts_prods.append(molecule)
                if not molecule.startswith('TS'):
                    non_ts.append(molecule)
                all_molecules.append(molecule)

    return non_ts_prods, non_ts
## -- Writes MESS file -- ##
def writeMessFile():
    global non_ts_prods, non_ts
    non_ts = []
    non_ts_prods = []
    #get current working directory (root)
    root_directory = os.getcwd()
    write_Header()

    # list that includes all reaction pathways (W1, P1, P2 ... etc.). Makes sure that wells (W) go first!
    pathways = []
    for item in os.listdir():
        pathway_path = os.getcwd()+'\\'+item
        if os.path.isdir(pathway_path) and item[0] in ['P', 'W']:
            pathways.append(item)

    pathways.sort(key = lambda pathways: pathways[0], reverse=True)
    pathways.sort(key = lambda pathways: int(pathways[1:]))

    # counts the number of wells in the mechanism, excludes parent well
    number_of_wells = sum('W' in path for path in pathways) - 1
    lastParse = False
    for pathway in pathways:
        print('Working on: ', pathway)
        if pathway == pathways[-1]:
            lastParse = True
        path_dir = root_directory+'\\'+pathway
        #looks into the TS folder and extract xyz geom from m062x file
        os.chdir(path_dir)

        m062xFileList = [filename for filename in os.listdir() if filename.endswith("m062x.log")]
        m062xFileName = m062xFileList[0]

        # writes well section
        if pathway.startswith('W'):
            # adds energies from all fragment files into 'energies' dictionary
            try:
                non_ts_prods, non_ts = add_energies(pathway, os.listdir(), non_ts_prods, non_ts, energiesDictionary, createExcel)
            except PermissionError:
                pass

            if len(m062xFileList)>1:
                # extracts Well's filename for Well Section
                m062xFileName   = [file for file in m062xFileList if file.startswith('W')]
                m062xFileName = m062xFileName[0]

            # setting bond threshold for non-TS molecules
            bond_thresh = 1.2
            ## -- Topology Function Calls For Wells -- ##
            # get xyz geometry and filename
            W_xyz_geom = get_xyzGeom(m062xFileName)
            # get frequencies from m062x file
            W_freqs, W_charge_multiplicity = get_frequencies(m062xFileName)
            # get rotors from geometry and get_geometry (rotors script) functions
            W_rots = get_rotors(W_xyz_geom, m062xFileName, bond_thresh)
            # finds corresponding z-matrix dihedral for each rotor
            W_dihedrals = ['D'+str(max(rot)-3) for rot in W_rots]

            print(W_dihedrals, m062xFileName)
            rotorEnergies = get_peaks_valleys(W_dihedrals, m062xFileName)

            # Change back to root directory and write in Mess input there!
            os.chdir(root_directory)
            write_Well(W_xyz_geom, W_freqs, W_charge_multiplicity, W_rots, rotorEnergies, path_dir, m062xFileName)

            if len(m062xFileList)>1:
                # extracts name of Well's TS for Well's Barrier Section
                TSWellFileName   = [file for file in m062xFileList if file.startswith('TS')]
                TSWellFileName = TSWellFileName[0]

                #looks into the TS folder and extract xyz geom from m062x file
                os.chdir(path_dir)

                # change threshold to 1.4 for all files here bc they are TS!
                bond_thresh = 1.4
                ## -- Topology Function Calls For Well's TS -- ##
                TSW_xyz_geom = get_xyzGeom(TSWellFileName)
                TSW_freqs, TSW_charge_multiplicity = get_frequencies(TSWellFileName)
                TSW_rots = get_rotors(TSW_xyz_geom, TSWellFileName, bond_thresh, TS=True)
                # finds corresponding z-matrix dihedral for each rotor
                TSW_dihedrals = ['D'+str(max(rot)-3) for rot in TSW_rots]

                rotorEnergies = get_peaks_valleys(TSW_dihedrals, TSWellFileName)

                # Change back to root directory and write in Mess input there!
                os.chdir(root_directory)
                write_Well_Barrier(TSW_xyz_geom, TSW_freqs, TSW_charge_multiplicity, TSW_rots, rotorEnergies, path_dir, m062xFileName)

                continue
            continue

        if m062xFileName.startswith('TS'):
            TS = True
        else:
            TS = False
            bond_thresh = 1.2

        # This section writes transition states
        if TS:
            # adds energies from all fragment files into 'energies' dictionary
            try:
                non_ts_prods, non_ts = add_energies(pathway, os.listdir(), non_ts_prods, non_ts, energiesDictionary, createExcel)
            except PermissionError:
                pass

            fragments_directory = root_directory+'\\'+pathway+'\\fragments'
            # changes to fragments' directory to write bimolecular section
            os.chdir(fragments_directory)

            # adds energies from all fragment files into 'energies' dictionary
            non_ts_prods, non_ts = add_energies(pathway, os.listdir(), non_ts_prods, non_ts, energiesDictionary, createExcel)

            m062xFragmentFileList = [filename for filename in os.listdir() if filename.endswith("m062x.log")]
            if len(m062xFragmentFileList)>1:
                h_abstraction = False
                # gets fragments' names from file for parsing
                frag1 = m062xFragmentFileList[0].split('_m062x.log')[0]
                frag2 = m062xFragmentFileList[1].split('_m062x.log')[0]
                # takes into account the cases with 3 products!
                if len(m062xFragmentFileList) > 2:
                    frag3 = m062xFragmentFileList[2].split('_m062x.log')[0]
                    fragments_naming = frag1+' '*2+'+'+' '*2+frag2+' '*2+'+'+frag3
                else:
                    fragments_naming = frag1+' '*2+'+'+' '*2+frag2
            # Accounts for cases where there is only one product m062x file, H- abstractions!
            else:
                h_abstraction = True
                non_ts.append('H')
                fragments_naming = m062xFragmentFileList[0].split('_m062x.log')[0]
            ## -- Topology Function Calls For Fragments -- ##
            # change threshold to 1.2 for all fragment files here bc they are not a TS!
            bond_thresh = 1.2
            iteration = 1
            for fragmentFile in m062xFragmentFileList:
                fr_xyz_geom = get_xyzGeom(fragmentFile)
                fr_freqs, fr_charge_multiplicity = get_frequencies(fragmentFile)
                fr_rots = get_rotors(fr_xyz_geom, fragmentFile, bond_thresh)
                # finds corresponding z-matrix dihedral for each rotor
                fr_dihedrals = ['D'+str(max(rot)-3) for rot in fr_rots]

                print(fr_dihedrals, m062xFragmentFileList)
                rotorEnergies = get_peaks_valleys(fr_dihedrals, fragmentFile)

                # go back to root directory
                os.chdir(root_directory)
                # Writes TS Bimolecular Section from topology ooutputs! (calling this function opens/closes the file!)
                write_TS_Bimolecular(iteration, fragments_naming, fragmentFile, fragments_directory, path_dir, m062xFragmentFileList, fr_xyz_geom, fr_freqs, fr_rots, rotorEnergies, fr_charge_multiplicity, h_abstraction)
                iteration += 1
            # resetting h-abstraction boolean
            h_abstraction = False
            #looks into the TS folder and extract xyz geom from m062x file
            os.chdir(path_dir)

            # change threshold to 1.4 for all files here bc they are TS!
            bond_thresh = 1.4
            ## -- Topology Function Calls For TS -- ##
            TS_xyz_geom = get_xyzGeom(m062xFileName)
            TS_freqs, TS_charge_multiplicity = get_frequencies(m062xFileName)
            TS_rots = get_rotors(TS_xyz_geom, m062xFileName, bond_thresh, TS=True)
            # finds corresponding z-matrix dihedral for each rotor
            TS_dihedrals = ['D'+str(max(rot)-3) for rot in TS_rots]

            print(TS_dihedrals, m062xFileName)
            rotorEnergies = get_peaks_valleys(TS_dihedrals, m062xFileName)

            # Change back to root directory and write in Mess input!
            os.chdir(root_directory)
            # Writes TS Barrier Section from topology outputs! (calling this function closes the file!)
            write_TS_Barrier(TS_xyz_geom, TS_rots, TS_freqs, rotorEnergies, TS_charge_multiplicity, path_dir, m062xFileName, number_of_wells, h_abstraction)

        # This section writes BFs
        else:
            # goes back into the pathway's directory
            os.chdir(path_dir)
            # adds energies from all files into 'energies' dicionary
            non_ts_prods, non_ts = add_energies(pathway, os.listdir(), non_ts_prods, non_ts, energiesDictionary, createExcel)

            if len(m062xFileList)>1:
                h_abstraction = False
                # gets fragments' names from file for parsing
                frag1 = m062xFileList[0].split('_m062x.log')[0]
                frag2 = m062xFileList[1].split('_m062x.log')[0]
                # takes into account the cases with 3 products!
                if len(m062xFileList) > 2:
                    frag3 = m062xFileList[2].split('_m062x.log')[0]
                    fragments_naming = frag1+' '*2+'+'+' '*2+frag2+' '*2+'+'+frag3
                else:
                    fragments_naming = frag1+' '*2+'+'+' '*2+frag2
            else:
                h_abstraction = True
                non_ts_prods.append('H')
                non_ts.append('H')
                fragments_naming = m062xFileList[0].split('_m062x.log')[0]

            # Function calls for Bimolecular Section of Bond Fissions
            iteration = 1
            for fragmentFile in m062xFileList:
                fr_xyz_geom = get_xyzGeom(fragmentFile)
                fr_freqs, fr_charge_multiplicity = get_frequencies(fragmentFile)
                fr_rots = get_rotors(fr_xyz_geom, fragmentFile, bond_thresh)
                # finds corresponding z-matrix dihedral for each rotor
                fr_dihedrals = ['D'+str(max(rot)-3) for rot in fr_rots]

                rotorEnergies = get_peaks_valleys(fr_dihedrals, fragmentFile)

                # go back to root directory
                os.chdir(root_directory)
                # Writes TS Bimolecular Section from topology ooutputs! (calling this function opens/closes the file!)
                write_BF_Bimolecular(iteration, root_directory, fragments_naming, fragmentFile, path_dir, m062xFileList, fr_xyz_geom, fr_freqs, fr_rots, rotorEnergies, fr_charge_multiplicity, h_abstraction)
                iteration += 1

            # Change back to root directory and create Mess input there (TS Bimolecular Sections)!
            os.chdir(root_directory)
            # Change back into pathway's directory
            os.chdir(path_dir)

            # gets geometry, frequency and rotor information for the two fragments inside the BF pathway
            # fragment 1 information
            fr1_xyz_geom = get_xyzGeom(m062xFileList[0])
            fr1_freqs, fr1_charge_multiplicity = get_frequencies(m062xFileList[0])
            fr1_rots = get_rotors(fr1_xyz_geom, m062xFileList[0], bond_thresh)
            # finds corresponding z-matrix dihedral for each rotor
            fr1_dihedrals = ['D'+str(max(rot)-3) for rot in fr1_rots]

            fr1_rotorEnergies = get_peaks_valleys(fr1_dihedrals, m062xFileList[0])

            fragment1_info = [fr1_xyz_geom, fr1_freqs, fr1_charge_multiplicity, fr1_rots]

            if len(m062xFileList)>1:
                # fragment 2 information
                fr2_xyz_geom = get_xyzGeom(m062xFileList[1])
                fr2_freqs, fr2_charge_multiplicity = get_frequencies(m062xFileList[1])
                fr2_rots = get_rotors(fr2_xyz_geom, m062xFileList[1], bond_thresh)
                # finds corresponding z-matrix dihedral for each rotor
                fr2_dihedrals = ['D'+str(max(rot)-3) for rot in fr2_rots]

                fr2_rotorEnergies = get_peaks_valleys(fr2_dihedrals, m062xFileList[1])

                fragment2_info = [fr2_xyz_geom, fr2_freqs, fr2_charge_multiplicity, fr2_rots]
            # Change back to root directory and create Mess input there (TS Bimolecular Sections)!
            os.chdir(root_directory)
            # Writes BF Barrier Section from topology outputs! (calling this function closes the file!)
            if h_abstraction:
                fragment2_info = []
            write_BF_Barrier(fragments_naming, fragment1_info, fr1_rotorEnergies, fragment2_info, fr2_rotorEnergies, pathway, number_of_wells, m062xFileName, h_abstraction)

        # change back to root dir for next pathway
        os.chdir(root_directory)
        if lastParse:
            MESSFile = open('MESSFile.inp', 'a')
            MESSFile.write('End\n')
            MESSFile.close()

writeMessFile()

# Additional data structures needed for creating energies spreadsheet
ts_energies = {}
for key in energiesDictionary:
    if key not in non_ts_prods:
        ts_energies[key] = energiesDictionary[key]

dissociations = [prod for prod in non_ts_prods if not prod.startswith('W')]
# print('dissociations: ', dissociations)
print(f"EnergiesDictionary {energiesDictionary}")


all_ts_molecules = [molecule for molecule in all_molecules if molecule not in dissociations]
ts_prods = {}
for i in range(len(all_molecules)):
    try:
        if all_molecules[i].startswith('TS'):
            # this logic is why it doesn't work for more than two molecule products
            if not all_molecules[i+1].startswith('TS'):
                if all_molecules[i+1].startswith('W'):
                    ts_prods[all_molecules[i]] = [all_molecules[i+1]]
                else:
                    if not all_molecules[i+2].startswith('TS'):
                        ts_prods[all_molecules[i]] = [all_molecules[i+1], all_molecules[i+2]]
                    else:
                        ts_prods[all_molecules[i]] = [all_molecules[i+1]]

    except IndexError:
        ts_prods[all_molecules[i]] = [all_molecules[i+1]]

# print('all ts molecules', all_ts_molecules)
# print('ts_prods', ts_prods)
# print('ts_energies', ts_energies)

non_repeated_non_ts = sorted(set(non_ts), key=non_ts.index)
# print('non repeated non ts', non_repeated_non_ts)

if createExcel:
    write_energies_spreadsheet(energiesDictionary, non_repeated_non_ts, ts_prods, dissociations)
else:
    pass

print('Seems like it worked!')
# end of program
