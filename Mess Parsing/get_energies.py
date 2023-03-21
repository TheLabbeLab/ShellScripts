import xlsxwriter
import numpy as np
import re
from scipy.signal import find_peaks

# Creates and populates energies spreadsheet from input energies dictionary
def write_energies_spreadsheet(energiesDictionary, non_repeated_non_ts, ts_prods, dissociations):
    # Create a workbook and add a worksheet.
    workbook = xlsxwriter.Workbook('Energies_SpreadSheet.xlsx')
    worksheet = workbook.add_worksheet('Energetics')
    # Some data we want to write to the worksheet.
    column_titles = ['m062X/cc-pVTZ', 'CCSD(T)/cc-pVDZ', 'CCSD(T)/cc-pVTZ', 'T1 diag', 'E(hartrees)']
    # Start from the first cell. Rows and columns are zero indexed.
    row = 0
    col = 0

    # Iterate over the data and write it out row by row.
    for column_title in (column_titles):
        col += 1
        worksheet.write(row, col, column_title)
    col = 0
    row += 1
    worksheet.write(row, col, 'Reactants/Products')

    row += 1

    # writes all product energies into excel file
    for prod in non_repeated_non_ts:
        worksheet.write(row, col, prod)
        col = 1
        for energy in energiesDictionary[prod]:
            worksheet.write(row, col, energy)
            col += 1
        # extrapolation of ground energies to inf zeta limit
        worksheet.write(row, col, '=B'+str(row+1)+'+D'+str(row+1)+'+((D'+str(row+1)+'-C'+str(row+1)+')*0.46286)')
        row += 1
        col = 0

    row += 1
    worksheet.write(row, col, 'Transition States')
    worksheet.write(row, col+6, 'E_fwd (hartrees)')
    worksheet.write(row, col+7, 'E_rev (hartrees)')
    worksheet.write(row, col+8, 'E_fwd (kcal/mol)')
    worksheet.write(row, col+9, 'E_rev (kcal/mol)')

    # writes all transition state energies into excel file
    row += 1
    for key in energiesDictionary.keys():
        if key.startswith('TS'):
            temp = set(ts_prods[key])
            prod_indices = [i for i, val in enumerate(non_repeated_non_ts) if val in temp]
            h_index = [i for i, val in enumerate(non_repeated_non_ts) if val == 'H']
            worksheet.write(row, col, key)
            col = 1
            for energy in energiesDictionary[key]:
                worksheet.write(row, col, energy)
                col += 1
            worksheet.write(row, col, '=B'+str(row+1)+'+D'+str(row+1)+'+((D'+str(row+1)+'-C'+str(row+1)+')*0.46286)')
            col += 1
            worksheet.write(row, col, '=F'+str(row+1)+'-$F$3')
            col += 1
            if len(prod_indices)==2:
                worksheet.write(row, col, '=F'+str(row+1)+'-(F'+str(prod_indices[0]+3)+'+F'+str(prod_indices[1]+3)+')')
            else:
                if not list(temp)[0].startswith('W'):
                    worksheet.write(row, col, '=F'+str(row+1)+'-(F'+str(prod_indices[0]+3)+'+F'+str(h_index[0]+3)+')')
                else:
                    worksheet.write(row, col, '=F'+str(row+1)+'-(F'+str(prod_indices[0]+3)+')')
            col += 1
            worksheet.write(row, col, '=G'+str(row+1)+'*627.51')
            col += 1
            worksheet.write(row, col, '=H'+str(row+1)+'*627.51')

            row += 1
            col = 0
    # writes all dissociation energies into excel file
    row += 1
    worksheet.write(row, col, 'Dissociations')
    worksheet.write(row, col+6, 'E_fwd (hartrees)')
    worksheet.write(row, col+7, 'E_rev (hartrees)')
    worksheet.write(row, col+8, 'E_fwd (kcal/mol)')
    worksheet.write(row, col+9, 'E_rev (kcal/mol)')

    row += 1
    for i in range(0, len(dissociations), 2):
        dissoc_indices = [j for j, val in enumerate(non_repeated_non_ts) if val in dissociations[i:i+2]]
        worksheet.write(row, col, dissociations[i]+' + '+dissociations[i+1])
        col += 6
        worksheet.write(row, col,  '=(F'+str(dissoc_indices[0]+3)+'+F'+str(dissoc_indices[1]+3)+')-F3')
        col += 1
        worksheet.write(row, col, '=-G'+str(row+1))
        col += 1
        worksheet.write(row, col, '=G'+str(row+1)+'*627.51')
        col += 1
        worksheet.write(row, col, '=H'+str(row+1)+'*627.51')

        col = 0
        row += 1

    workbook.close()

# Takes in rotor dihedrals, navigates to corresponding rotor output and returns peaks and valleys using SciPy' find_peaks
def get_peaks_valleys(dihedrals, m062xFileName):
    moleculeName = m062xFileName.split('_m062x.log')[0]
    # dictionary, conatins dihedrals as keys and extreme energies lists as values
    rotorEnergies = {}
    degrees = np.arange(0, 36, 1)
    for dihedral in dihedrals:
        peaks_and_valleys = []
        rotorFileName = moleculeName + '_' + dihedral + '_rotor.log'
        rotorFile = open(rotorFileName, 'r')
        rotorFileString = rotorFile.read()
        rotorFile.close()
        rotorFileString = rotorFileString.replace('\n', '').replace(' ', '')
        # conditional for a different command in rotor header
        if "#scan" in rotorFileString:
            # use regex to find scanned energies at bottom of the file
            print(rotorFileName)
            energies = re.findall(r'\\HF=[-?\d*\.\d*,]*', rotorFileString)[0].split(',')
            energies[0] = energies[0].split('=')[-1]
        else:
            rotorFile = open(rotorFileName, 'r')
            rotorFileString = rotorFile.readlines()
            rotorFile.close()
            # use regex to find scanned energies throughout the file
            print(rotorFileName)
            search_optimization_line = re.compile(r'Optimization completed')
            search_SCF_Done_line = re.compile(r'SCF Done')
            num_lines = len(rotorFileString)
            curr_line_num = 0
            energies = []
            while curr_line_num < num_lines:
                curr_line = str.strip(rotorFileString[curr_line_num])
                if search_SCF_Done_line.search(curr_line):
                    SCF_Done_Line = curr_line
                if search_optimization_line.search(curr_line):
                    SCF_Done_line_contents = SCF_Done_Line.split()
                    energies.append(str(SCF_Done_line_contents[4]))
                curr_line_num = curr_line_num + 1

        # adjusts energies with respect to initial energy and converts to kcal/mol
        energies = np.array([float(energy)-float(energies[0]) for energy in energies][:-1])*627.51
        peaks_and_valleys.append(energies[0])

        # finds peaks and valleys using scipy.signal.find_peaks
        peak_indices, _ = find_peaks(energies)
        valley_indices, _  = find_peaks(-1*energies)

        # creates list with all extreme energies
        for i in range(len(peak_indices)):
            try:
                peaks_and_valleys.append(energies[peak_indices[i]])
                peaks_and_valleys.append(energies[valley_indices[i]])
            except IndexError:
                pass
        # populates dicionary of dihedral-energies
        rotorEnergies[dihedral] =  peaks_and_valleys

    return rotorEnergies
