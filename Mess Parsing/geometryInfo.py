# Function for extracting XYZ Geometries from m062x files
import os
import re
from difflib import SequenceMatcher
from get_geometry import get_geom, get_bond_graph, get_bonds, get_angles, get_torsions, get_all_rotors

## ------------------------------------------------- TOPOLOGY FUNCTIONS ------------------------------------------------------##

# get xyz geometry from the m062x file
def get_xyzGeom(m062xFileName):
    #open and strip \n from file, locate geom using regex and remove duplicates by using non-duplicate keys feature from dicts.
    #Regex: \d (digit), [A-Z] (caps letter), -? (0 or more appearances of '-'), \. (.), * (0 or more appearances of preceding call)
    m062xFile = open(m062xFileName, 'r')
    m062xFileString = m062xFile.read().replace('\n', '').replace(' ', '')
    xyz_geom = list(dict.fromkeys(re.findall(r'[A-Z],-?\d\.\d*,-?\d\.\d*,-?\d\.\d*', m062xFileString)))
    idxs_toDelete = []
    for i in range(len(xyz_geom)):
        for j in range(len(xyz_geom)):
            # library used to compare how 2 strings match
            # used to delete repeated coordinates that differ to themeselves by a few orders of magnitude (redundant!)
            # ratios > 0.8 are accurate by 6 decimal places!
            if SequenceMatcher(a=xyz_geom[i], b=xyz_geom[j]).ratio() > 0.8 and SequenceMatcher(a=xyz_geom[i], b=xyz_geom[j]).ratio() < 1.0:
                if sorted([i, j]) not in idxs_toDelete:
                    idxs_toDelete.append([i, j])
    # deletes repeated coordinates!
    idxs_toDelete = [max(item) for item in idxs_toDelete]
    for index in sorted(idxs_toDelete, reverse=True):
        del xyz_geom[index]

    m062xFile.close()
    return xyz_geom

# slightly different from rotor script, returns lists of rotors via atom numbering
def get_rotor_dihedrals(torsions, bond_graph, bond_graph_woTSthresh=[]):
    axes = []
    rotors = []
    rotor_rows = []
    maxAtomNumber = []

    rotor_axes = get_all_rotors(bond_graph, bond_graph_woTSthresh)

    for torsion in torsions:
        ax = sorted(torsion[1:3])
        if (ax in rotor_axes) and (ax not in axes) and (max(torsion[:4]) not in maxAtomNumber):
            rotors.append(torsion[:4])
            axes.append(ax)
            maxAtomNumber.append(max(torsion[:4]))

    return rotors

# pulls functions from get get_geomtry to get list of molecule's rotors
def get_rotors(xyz_geom, m062xFileName, bond_thresh, TS=False):
    #converts xyz_geom to xyz_array format compatible with get_rotors_gjf.py script!
    xyz_array = []
    xyz_array.append([len(xyz_geom)])
    xyz_array.append([m062xFileName])
    for i in range(len(xyz_geom)):
        line = [float(element) if len(element)>1 else element for element in xyz_geom[i].split(',')]
        xyz_array.append(line)
    ## -- Using functions from get_geometry script -- ##
    #read in geometry, determine bonded topology
    geom = get_geom(xyz_array)
    bond_graph, double_bonds = get_bond_graph(geom, bond_thresh)
    if TS:
        bond_thresh = 1.2
        bond_graph_woTSthresh, double_bonds = get_bond_graph(geom, bond_thresh)
    #calculate bond lengths, angles, and torsions
    bonds = get_bonds(geom, bond_graph)
    angles = get_angles(geom, bond_graph)
    torsions = get_torsions(geom, bond_graph)
    if TS:
        rotors = get_rotor_dihedrals(torsions, bond_graph, bond_graph_woTSthresh)
    else:
        rotors = get_rotor_dihedrals(torsions, bond_graph)
    # filtering out double bonds from possible rotors
    rotors = [rotor for rotor in rotors if rotor[1:3] not in double_bonds]
    #readjusting indices
    for rotor in rotors:
        for i in range(len(rotor)):
            rotor[i] += 1

    return rotors

#extracts frequencies as well as charge and multiplicity! - Taken from Katie's code
def get_frequencies(m062xFilename):
    m062xFile = open(m062xFilename, 'r')
    m062xFileLines = m062xFile.readlines()
    m062xFile.close()

    all_frequencies = ""
    for line in m062xFileLines:
        # finds frequencies
        if "Frequencies --" in line:
            freq_line=line.strip('\n')
            freq_line=freq_line.strip('Frequencies --')+ '   '
            freq_line=freq_line.replace('             ', '   ')
            #stores all frequenices in one string
            all_frequencies += freq_line
        # finds charge and multiplicity
        if 'Multiplicity' in line:
            charge_multiplicity = [item for item in line if item.isdigit()]

    #'all_frequencies' contains repeated elements
    all_frequencies = [float(freq) for freq in all_frequencies.split()]
    frequencies = []
    #adding the imaginary (negative) frequency to the completed lise (Only for transition states!!!)
    frequencies.append(all_frequencies[0])
    #eliminating duplicates for a complete frequency list, starts from second element!
    for i in range(len(all_frequencies[1:])):
        # try block used for some cases were frequencies are not repeated!
        try:
            if all_frequencies[1:][i+1] < all_frequencies[1:][i]:
                frequencies.append(all_frequencies[1:][i])
                break
        except IndexError:
            frequencies.append(all_frequencies[1:][i])
            break
        else:
            frequencies.append(all_frequencies[1:][i])

    return frequencies, charge_multiplicity
