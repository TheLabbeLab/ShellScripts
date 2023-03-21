
import warnings
warnings.filterwarnings(action='ignore') 
import os
import shutil
import numpy as np
import yaml
import time
import colorama
import CodeFunctions as cf

colorama.init(autoreset = True)

print(colorama.Fore.BLUE  + colorama.Style.BRIGHT+ """
====================================
|| Title: Arrheinus Fit Generator ||
|| Author: Pray Shah              ||
|| EMAIL: prsh1291@colorado.edu   ||
====================================
""")

base_dir = os.getcwd() 
if os.path.exists("Plots"):
    shutil.rmtree("Plots")
    os.mkdir("Plots")
else:
    os.mkdir("Plots")

FileName = input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter File Name with extension: ")

filePath = os.path.join(base_dir, FileName)

while True:
    
    if input(colorama.Fore.GREEN + colorama.Style.BRIGHT + "\nDo you want to input from 'UserInput.yaml' file? (y/n): ")[0].lower() == 'y':
        UserInput = yaml.load(open('UserInput.yaml', 'r'), Loader=yaml.Loader)
        Temps = np.array(UserInput['Temperature'])
        Pressures = np.array(UserInput['Pressure'])
        ChemkinNameDict = UserInput['Species']
        MappedKeyword = list(UserInput['Species'].keys())
    
    else:
        print("Need to specify temperautre start, stop and step size:")

        Temps = np.arange(int(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Please enter start temperature: ")), 
                        int(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Please enter end temperature: ")) + 1, 
                        int(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Please enter step size: ")))

        print("\nEnter the pressures for which the arrhenius fits will be generated")
        print(colorama.Fore.YELLOW + colorama.Style.BRIGHT +  "Specify all pressures in one line and use ',' in between pressures. Do not use any other delimiter!\n")

        Pressures = np.array(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter Pressures: ").replace(" ","").split(','))
        
        print(colorama.Fore.YELLOW + colorama.Style.BRIGHT + "Enter one at a time!")

        check = 'y' 
        MappedKeyword = []
        SpeciesName = []

        while check == 'y':    
            MappedKeyword.append(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter the keyword: "))
            SpeciesName.append(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter the species name: "))
            check = input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Have more to map? (y/n): " )[0].lower()

        ChemkinNameDict = dict(zip(MappedKeyword, SpeciesName))
    
    print(colorama.Fore.YELLOW + colorama.Style.NORMAL + f"\nTemperatures: {Temps}")
    print(colorama.Fore.YELLOW + colorama.Style.NORMAL + f"\nPressures: {Pressures}")
    print(colorama.Fore.YELLOW + colorama.Style.NORMAL + '\nSpecies Mapped:')
    for key, value in ChemkinNameDict.items():
        print(colorama.Fore.YELLOW + f'{key} : {value}')
    
    check = input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "\nIs the order correct? (y/n): ")[0].lower()
    
    if check == 'y':
        break
    else:
        print(colorama.Fore.RED + "Re-Enter the species map!")

MappedKeyword = np.array(MappedKeyword)

Parent_HighP, Parent_rates = cf.MessExtraction.messExtract(filePath, Temps, Pressures, MappedKeyword) 

print(colorama.Fore.YELLOW + colorama.Style.BRIGHT + "\nData Extracted for requested pressures and temperatures successfully!\n")

# Weighted Rates

WeightedRatesUserInp = input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Do you want to weight the rates for a specific temperature range? (y/n): ").lower()
if WeightedRatesUserInp == 'y': 
    WeightStartEndTemps = np.array([int(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter the inital temperature: ")), int(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter the final temperature: "))])
    Repitition = int(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter the repitions you want of the rates for give range of temperature: ")) - 1
    
    StartIndex  = np.where(WeightStartEndTemps[0] == Parent_HighP[0][0])[0][0]
    EndIndex    = np.where(WeightStartEndTemps[1] == Parent_HighP[0][0])[0][0] + 1

    SubArray = Parent_HighP[:,:,StartIndex:EndIndex]
    SubArray = np.repeat(SubArray, Repitition, axis=2)
    Parent_HighP = np.concatenate((Parent_HighP, SubArray), axis = 2)
    Parent_HighP = Parent_HighP[:,:,Parent_HighP[0][0].argsort()]

    SubArray = Parent_rates[:,:,StartIndex:EndIndex]
    SubArray = np.repeat(SubArray, Repitition, axis=2)
    Parent_rates = np.concatenate((Parent_rates, SubArray), axis = 2)
    Parent_rates = Parent_rates[:,:,Parent_rates[0][0].argsort()]

# Tolerance

Tolerance = int(input("\nPlease enter the error tolerance in percentage: " ))

# Arrhenius Fit

file = [open("Arrhenius_Fit_Chemkin_Format_ModArr.inp","w"), open("Arrhenius_Fit_Chemkin_Format_DoubleArr.inp","w"), open("Arrhenius_Fit_Chemkin_Format.inp", "w")]

while True:
    Index = 0
    ReactantKeyword = ''
    ReactantKeyword = input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "\nEnter the reactant keyword for which you want rates (e.g: W1, P1 or P6): ")
    if ReactantKeyword not in MappedKeyword:
        print(colorama.Fore.RED + colorama.Style.BRIGHT + "\nError while entering Reactant keyword!")
        print(colorama.Fore.RED + colorama.Style.BRIGHT + "Keyword not found!")
        print(colorama.Fore.RED + colorama.Style.BRIGHT + "Re-Enter the keyword")
    else:
        Bimolec = False
        if '+' in ChemkinNameDict[ReactantKeyword]:
            print(colorama.Fore.RED + colorama.Style.BRIGHT + "\nReaction is Bimolecular")
            Bimolec = True
        else:
            print(colorama.Fore.RED + colorama.Style.BRIGHT + "\nReaction is Unimolecular")
            Bimolec = False        
            
        RatePlotDict = {}
        Index = np.where(MappedKeyword == ReactantKeyword)[0][0]
        MappedForPlot = np.delete(MappedKeyword, [Index])
        for i in range(0, len(MappedForPlot)):
            RatePlotDict[i + 1] = f'{ChemkinNameDict[ReactantKeyword]} = {ChemkinNameDict[MappedForPlot[i]]}'
    
        print(colorama.Fore.YELLOW + "\nMapped dictionary generated:")
        for key, value in RatePlotDict.items():
            print(colorama.Fore.BLUE + colorama.Style.BRIGHT + f"{key} : " + colorama.Fore.YELLOW + f"{value}")

        print(colorama.Fore.YELLOW + colorama.Style.BRIGHT + "\nSpecify all keys in one line and use ',' in between keys. Do not use any other delimiter!\n")
        SpeciesIndex = np.array(input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "Enter the key number from the above dictionary for which you want Arrhenius expression: " ).replace(" ","").split(',')).astype(int)

        print(colorama.Fore.BLUE + colorama.Style.BRIGHT + """\n
        ===========
        ||Writing||
        ===========
        """)
        StartTime = time.time()
        if WeightedRatesUserInp == 'y':
            cf.WriteChemkinFile(file, RatePlotDict, Index, Pressures, SpeciesIndex, Parent_HighP, Parent_rates, Tolerance, WeightStartEndTemps, Bimolec).WriteWeightedRates()
        else:
            cf.WriteChemkinFile(file, RatePlotDict, Index, Pressures, SpeciesIndex, Parent_HighP, Parent_rates, Tolerance, Bimolec=Bimolec).WriteRates()
        EndTime = time.time()
        print(colorama.Back.GREEN + f"\nWrite Time: {EndTime - StartTime} seconds")
        if input(colorama.Fore.GREEN  + colorama.Style.BRIGHT + "\nDo you want to print any more rates? (y/n): ")[0].lower() == 'y':
            continue
        else:
            break

file[0].close()
file[1].close()
file[2].close()

print(colorama.Fore.BLUE + colorama.Style.BRIGHT + "\nArrhenius Expressions Printed!!\n")