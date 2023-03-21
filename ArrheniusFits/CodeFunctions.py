import os 
import scipy.optimize
import numpy as np
import re
import matplotlib.pyplot as plt
import datetime

# MESS Data Extraction Function
class MessExtraction:

    def messExtract(filePath, T, P, keywords, pUnit='atm'):   

        keywords = np.append(keywords, "Temperature-Pressure")

        file = open(filePath,'r')
        content = file.readlines()
        file.close()

        #-------------------------------------------------------------------------------------------------------------------------------#

        # Lines to find High Pressure Data

        HighP_start = re.compile(r"Temperature-Species Rate Tables")
        HighP_end = re.compile(r"Capture/Escape Rate Coefficients")
        line_num = 0

        while line_num < len(content):
            line = str.strip(content[line_num])
            if HighP_start.search(line):
                # print(line, line_num)
                HighP_start_num = line_num
            if HighP_end.search(line):
                # print(line, line_num)
                HighP_end_num = line_num
                break
            line_num = line_num + 1


        HighP_content = content[HighP_start_num : HighP_end_num]


        rates = []
        for key_num in range(0, len(keywords) - 1):

            keyword_start = keywords[key_num]

            keyword_end = keywords[key_num + 1]
            # print(keyword_end)

            line_num = 0
            # start_num = 0
            for line_num in range(0, len(HighP_content)):
                line = str.strip(HighP_content[line_num]).split()
                for keyword in line:
                    if keyword == keyword_start + '->' + keywords[0]:
                        start_num  =  line_num + 1
                        break
                if keyword == keyword_start + '->' + keywords[0]:
                    break

            rates.append(HighP_content[start_num: start_num + len(T)])
            
        for i in range(0, len(rates)):
            for j in range(0, len(rates[i])):
                rates[i][j] = str.strip(rates[i][j]).split()
                # for k in range(0, len(rates[i][j])):
                    # rates[i][j][k] = str.strip(rates[i][j][k]).split()

            # print(rates[1])

        # print(rates[0:2])

        rates = np.array(rates)
        rates[rates=='***'] = ['nan']
        rates = rates.astype(float)
        # print(rates[0])   

        # print(rates.shape)

        hpRates = np.transpose(rates, (0,2,1))


        #----------------------------------------------------------------------------------------------------------------------------------------#


        # Lines to find Temperature-Species Rate Tables

        temp_press_start = re.compile(r"Temperature-Species Rate Tables")
        temp_press_end = re.compile(r"Temperature-Pressure Rate Tables")
        line_num = 0

        while line_num < len(content):
            line = str.strip(content[line_num])
            if temp_press_start.search(line):
                # print(line, line_num)
                temp_press_start_num = line_num
            if temp_press_end.search(line):
                # print(line, line_num)
                temp_press_end_num = line_num
                break
            line_num = line_num + 1


        # print(temp_press_start_num, temp_press_end_num)

        temp_press_content = content[temp_press_start_num : temp_press_end_num + 1]


        # print(temp_press_content[-1])

        rates =[]
        start_num = 0
        for key_num in range(0, len(keywords) - 1):

            keyword_start = keywords[key_num]

            keyword_end = keywords[key_num + 1]
            # print(keyword_end)

            line_num = 0

            for line_num in range(0, len(temp_press_content)):
                line = str.strip(temp_press_content[line_num]).split()
                for keyword in line:
                    if keyword == keyword_start + '->':
                        start_num  =  line_num - 2
                        break
                if keyword == keyword_start + '->':
                    break

            for line_num in range(0, len(temp_press_content)):
                line = str.strip(temp_press_content[line_num]).split()
                for keyword in line:
                    if keyword == keyword_end + '->':
                        end_num  =  line_num - 3
                        break
                    elif keyword == "Temperature-Pressure":
                        end_num = line_num - 3
                        break
                if keyword == keyword_end + '->':
                    break
                elif keyword == "Temperature-Pressure":
                    break

            # print(start_num)
            # print(end_num)        
            # print(temp_press_content[start_num])
            # print(temp_press_content[end_num - 1])

            data = temp_press_content[start_num: end_num]
            # print(data)

            i = 0

            for p in P:
                # print(p)
                for i in range(0, len(data)):
                    if re.compile('Pressure = ' + str(p) + ' '+ pUnit).search(str.strip(data[i])):
                        line_num = i
                        rates.append(data[i + 3: i + 3 + len(T)])
                        # print(str.strip(data[i])) 
            # print(rates[0])

            # df = pd.DataFrame(rates[1])

        for i in range(0, len(rates)):
            for j in range(0, len(rates[i])):
                    rates[i][j] = str.strip(rates[i][j]).split()
                # for k in range(0, len(rates[i][j])):
                    # rates[i][j][k] = str.strip(rates[i][j][k]).split()

            # print(rates[1])

        # print(rates[0:2])

        rates = np.array(rates)
        rates[rates=='***'] = ['nan']
        rates = rates.astype(float)
        # print(rates[0])   

        # print(rates.shape)

        tpRates = np.transpose(rates, (0,2,1))
        # print(rates[55])
        # print(rates.shape)

        # rates_M1M_parent = rates

        # print(rates_M1M_parent.shape)
        # rates =[]
        return hpRates, tpRates

# Fit Functions
class FittingFunctions:
    def __init__(self, Bimolec):
        self.Bimolec = Bimolec
        self.Avagadro = 6.0221408e+23 # For converting Bimolecular rates from cm3/s-molecules to cm3 / s-mol
        
    def Linear(self, Ts, a, b):
        K = a*Ts + b
        if self.Bimolec == True:
            K = self.Avagadro * K
            return np.log(K)
        else:
            return np.log(K)

    def Arrhenius(self, Ts, A, Ea): #Ea given in cal/mol
        K = A*np.exp(-Ea/(1.9872*Ts))
        
        if self.Bimolec == True:
            K = self.Avagadro * K
            return np.log(K)
        else:
            return np.log(K)

    def ModArrhenius(self, Ts, A, n, Ea): #Ea given in cal/mol
        K = A*Ts**n*np.exp(-Ea/(1.9872*Ts))
        
        if self.Bimolec == True:
            K  = self.Avagadro * K
            return np.log(K)
        else:
            return np.log(K)

    def DoubArrhenius(self, Ts, A1, n1, Ea1, A2, n2, Ea2): #Ea given in cal/mol
        K1 = A1*Ts**n1*np.exp(-Ea1/(1.9872*Ts))
        K2 = A2*Ts**n2*np.exp(-Ea2/(1.9872*Ts))
        K = K1 + K2
        if self.Bimolec == True:
            K = self.Avagadro * K
            return np.log(K)
        else:
            return np.log(K)

# Percentage Error
class ErrorFunction:
    def ErrorPercentage(actual, calculated, precision = 3):
        return np.round(np.sum(abs(actual - calculated)/actual)*100/len(actual), precision), np.round((actual - calculated)*100/actual, precision)

# Plotting Class
class PlotFits:
    def __init__(self, Temp, Rates, RateName, Pressure, WeightsIndices=[] , Bimolec = False):
        self.Temp = Temp
        self.Rates = Rates
        self.RateName = RateName
        self.Pressure = Pressure
        self.WeightsIndices = WeightsIndices
        self.Bimolec = Bimolec

    def ArrheniusFitPlot(self):

        fig, ax = plt.subplots(2, figsize=(25,25), dpi=300)# for fits at different pressures

        ax[0].plot(1000/self.Temp, self.Rates, 'b-', linewidth=5, label=f"MESS Rate")# Arrh fit

        # ------------------------------------------------------------------------------------------------------- Mod Arrh fit -------------------------------------------------------------------------------------------------------------- #
        try:
            ModArr_popt, ModArr_pcov = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).ModArrhenius, self.Temp, np.log(self.Rates), maxfev=500000)
            ax[0].plot(1000/self.Temp, np.exp(FittingFunctions(self.Bimolec).ModArrhenius(self.Temp, ModArr_popt[0], ModArr_popt[1], ModArr_popt[2])), 'g--', linewidth=5, label=f'Modified Arrhenius Fit')
            ModArrMeanError, ModArrError = ErrorFunction.ErrorPercentage(self.Rates, np.exp(FittingFunctions(self.Bimolec).ModArrhenius(self.Temp, ModArr_popt[0], ModArr_popt[1], ModArr_popt[2])))

        except RuntimeError:
            ModArr_popt, ModArr_pcov, ModArrMeanError, ModArrError = [0,0,0], [0,0,0], np.inf, np.ones(len(self.Rates)) * np.inf
            

        ax[1].plot(1000/self.Temp, ModArrError, 'g--', linewidth=5, label='Modified Arrhenius Fit Error')

        # ------------------------------------------------------------------------------------------------------ Double Arrh fit --------------------------------------------------------------------------------------------------------------- #
        try:
            half = int(len(self.Temp)/2)
            popt1 ,pcov1 = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).ModArrhenius, self.Temp[:half], np.log(self.Rates[:half]), maxfev=500000)
            popt2 ,pcov2 = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).ModArrhenius, self.Temp[half:], np.log(self.Rates[half:]), maxfev=500000)
            init_guess = [popt1[0], popt1[1], popt1[2], popt2[0], popt2[1], popt2[2]]
            DoubleArr_popt, DoubleArr_pcov = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).DoubArrhenius, self.Temp, np.log(self.Rates), p0=init_guess, maxfev=500000)

            ax[0].plot(1000/self.Temp, np.exp(FittingFunctions(self.Bimolec).DoubArrhenius(self.Temp, DoubleArr_popt[0], DoubleArr_popt[1], DoubleArr_popt[2], 
                                                                        DoubleArr_popt[3], DoubleArr_popt[4], DoubleArr_popt[5])), 'r--', linewidth=5, label=f'Double Arrhenius Fit')

            DoubleArrMeanError, DoubleError = ErrorFunction.ErrorPercentage(self.Rates, np.exp(FittingFunctions(self.Bimolec).DoubArrhenius(self.Temp, DoubleArr_popt[0], DoubleArr_popt[1], DoubleArr_popt[2], 
                                                                                                            DoubleArr_popt[3], DoubleArr_popt[4], DoubleArr_popt[5])))
        except RuntimeError:
            DoubleArr_popt, DoubleArr_pcov, DoubleArrMeanError, DoubleError = [0,0,0,0,0,0], [0,0,0,0,0,0], np.inf, np.ones(len(self.Rates)) * np.inf

        ax[1].plot(1000/self.Temp, DoubleError, 'r--', linewidth=5, label='Double Arrhenius Fit Error')

        # ------------------------------------------------------------------------------------------------------ Plotting Parameters --------------------------------------------------------------------------------------------------------------- #

        ax[0].legend(fontsize=25)
        ax[0].set_title("Rates Vs 1000/Temperature", fontsize=25)
        ax[0].set_xlabel("1000/Temperature [K]", fontsize=25)
        ax[0].set_ylabel("k [1/s]", fontsize=25)
        ax[0].set_yscale("log")
        ax[0].tick_params(axis='both', labelsize=25)
        ax[0].grid()

        ax[1].legend(fontsize=25)
        ax[1].set_title("Percentage Error Vs 1000/Temperautre", fontsize=25)
        ax[1].set_xlabel("1000/Temperature [K]", fontsize=25)
        ax[1].set_ylabel("Error %", fontsize=25)
        ax[1].tick_params(axis='both', labelsize=25)
        ax[1].grid()

        fig.suptitle(f"""{self.RateName}
        Pressure = {self.Pressure}""" , fontsize=40)

        plt.savefig(os.path.join("Plots",f'{self.RateName}_P={self.Pressure}.jpg'), facecolor='w', dpi=300)
        plt.close()
        return ((ModArr_popt, ModArr_pcov, ModArrMeanError, ModArrError), (DoubleArr_popt, DoubleArr_pcov, DoubleArrMeanError, DoubleError))

    # ------------------------------------------------------------------------------------------------------ Weighted  --------------------------------------------------------------------------------------------------------------- #

    def WeightedArrheniusFitPlot(self):

        fig, ax = plt.subplots(2, figsize=(25,25), dpi=300)# for fits at different pressures

        ax[0].plot(1000/self.Temp, self.Rates, 'b-', linewidth=5, label=f"MESS Rate")# Arrh fit

        # ------------------------------------------------------------------------------------------------------- Mod Arrh fit -------------------------------------------------------------------------------------------------------------- #
        try:
            ModArr_popt, ModArr_pcov = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).ModArrhenius, self.Temp, np.log(self.Rates), maxfev=500000)
            ax[0].plot(1000/self.Temp, np.exp(FittingFunctions(self.Bimolec).ModArrhenius(self.Temp, ModArr_popt[0], ModArr_popt[1], ModArr_popt[2])), 'g--', linewidth=5, label=f'Modified Arrhenius Fit')
            
            ModArrMeanError, ModArrError = ErrorFunction.ErrorPercentage(self.Rates, np.exp(FittingFunctions(self.Bimolec).ModArrhenius(self.Temp, ModArr_popt[0], ModArr_popt[1], ModArr_popt[2])))

            TailoredModArrMeanError, TailoredModArrError = ErrorFunction.ErrorPercentage(self.Rates[self.WeightsIndices[0] : self.WeightsIndices[1]],
                                                            np.exp(FittingFunctions(self.Bimolec).ModArrhenius(self.Temp[self.WeightsIndices[0] : self.WeightsIndices[1]], 
                                                                                                ModArr_popt[0], ModArr_popt[1], ModArr_popt[2])))
        except RuntimeError:
            print("Got Runtime Error")
            ModArr_popt, ModArr_pcov, ModArrMeanError, ModArrError = [0,0,0], [0,0,0], np.inf, np.ones(len(self.Rates)) * np.inf
        
            TailoredModArrMeanError, TailoredModArrError = np.inf, np.ones(len(self.Rates[self.WeightsIndices[0] : self.WeightsIndices[1]])) * np.inf

        ax[1].plot(1000/self.Temp, ModArrError, 'g--', linewidth=5, label='Modified Arrhenius Fit Error')

        # ------------------------------------------------------------------------------------------------------ Double Arrh fit --------------------------------------------------------------------------------------------------------------- #
        try:
            half = int(len(self.Temp)/2)
            popt1 ,pcov1 = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).ModArrhenius, self.Temp[:half], np.log(self.Rates[:half]), maxfev=500000)
            popt2 ,pcov2 = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).ModArrhenius, self.Temp[half:], np.log(self.Rates[half:]), maxfev=500000)
            init_guess = [popt1[0], popt1[1], popt1[2], popt2[0], popt2[1], popt2[2]]
            DoubleArr_popt, DoubleArr_pcov = scipy.optimize.curve_fit(FittingFunctions(self.Bimolec).DoubArrhenius, self.Temp, np.log(self.Rates), p0=init_guess, maxfev=500000)

            ax[0].plot(1000/self.Temp, np.exp(FittingFunctions(self.Bimolec).DoubArrhenius(self.Temp, DoubleArr_popt[0], DoubleArr_popt[1], DoubleArr_popt[2], 
                                                                        DoubleArr_popt[3], DoubleArr_popt[4], DoubleArr_popt[5])), 'r--', linewidth=5, label=f'Double Arrhenius Fit')

            DoubleArrMeanError, DoubleError = ErrorFunction.ErrorPercentage(self.Rates, np.exp(FittingFunctions(self.Bimolec).DoubArrhenius(self.Temp, DoubleArr_popt[0], DoubleArr_popt[1], DoubleArr_popt[2], 
                                                                                                            DoubleArr_popt[3], DoubleArr_popt[4], DoubleArr_popt[5])))

            TailoredDoubleArrMeanError, TailoredDoubleError = ErrorFunction.ErrorPercentage(self.Rates[self.WeightsIndices[0] : self.WeightsIndices[1]],
                                                    np.exp(FittingFunctions(self.Bimolec).DoubArrhenius(self.Temp[self.WeightsIndices[0] : self.WeightsIndices[1]], DoubleArr_popt[0], DoubleArr_popt[1], 
                                                                                        DoubleArr_popt[2], DoubleArr_popt[3], DoubleArr_popt[4], DoubleArr_popt[5])))
        except RuntimeError:
            DoubleArr_popt, DoubleArr_pcov, DoubleArrMeanError, DoubleError = [0,0,0,0,0,0], [0,0,0,0,0,0], np.inf, np.ones(len(self.Rates)) * np.inf

            TailoredDoubleArrMeanError, TailoredDoubleError = np.inf, np.ones(len(self.Rates[self.WeightsIndices[0] : self.WeightsIndices[1]])) * np.inf

        ax[1].plot(1000/self.Temp, DoubleError, 'r--', linewidth=5, label='Double Arrhenius Fit Error')

        # ------------------------------------------------------------------------------------------------------ Plotting Parameters --------------------------------------------------------------------------------------------------------------- #

        ax[0].legend(fontsize=25)
        ax[0].set_title("Rates Vs 1000/Temperature", fontsize=25)
        ax[0].set_xlabel("1000/Temperature [K]", fontsize=25)
        ax[0].set_ylabel("k [1/s]", fontsize=25)
        ax[0].set_yscale("log")
        ax[0].tick_params(axis='both', labelsize=25)
        ax[0].grid()

        ax[1].legend(fontsize=25)
        ax[1].set_title("Percentage Error Vs 1000/Temperautre", fontsize=25)
        ax[1].set_xlabel("1000/Temperature [K]", fontsize=25)
        ax[1].set_ylabel("Error %", fontsize=25)
        ax[1].tick_params(axis='both', labelsize=25)
        ax[1].grid()

        fig.suptitle(f"""{self.RateName}
        Pressure = {self.Pressure}""" , fontsize=40)

        if len(self.WeightsIndices) != 0:
            ax[0].axvspan(1000/self.Temp[self.WeightsIndices[0]], 1000/self.Temp[self.WeightsIndices[1] - 1], color='yellow', alpha=0.5)
            ax[1].axvspan(1000/self.Temp[self.WeightsIndices[0]], 1000/self.Temp[self.WeightsIndices[1] - 1], color='yellow', alpha=0.5)
        
        plt.savefig(os.path.join("Plots",f'{self.RateName}_P={self.Pressure}.jpg'), facecolor='w', dpi=300)
        plt.close()

        return ((ModArr_popt, ModArr_pcov, ModArrMeanError, ModArrError,TailoredModArrMeanError, TailoredModArrError), 
                (DoubleArr_popt, DoubleArr_pcov, DoubleArrMeanError, DoubleError, TailoredDoubleArrMeanError, TailoredDoubleError))

# Chemkin File Write Function
class WriteChemkinFile:

    def __init__(self, file, RatePlotDict, HighPIndex, PlottingPressure, SpeciesIndex, HighPressureRates, TempSpeciesRates, Tolerance, WeightStartEndTemps=[], Bimolec = False):
        self.file = file
        self.RatePlotDict = RatePlotDict
        self.HighPIndex = HighPIndex
        self.PlottingPressure = PlottingPressure
        self.SpeciesIndex = SpeciesIndex
        self.HighPressureRates = HighPressureRates
        self.TempSpeciesRates = TempSpeciesRates
        self.Tolerance = Tolerance
        self.WeightStartEndTemps = WeightStartEndTemps
        self.Bimolec = Bimolec

    def __CheckTemps(self,Temp):
            
        try:
            self.WeightsIndices = np.array([np.where(Temp == self.WeightStartEndTemps[0])[0][0], np.where(Temp == self.WeightStartEndTemps[1])[0][-1] + 1])
        except IndexError:
            self.WeightsIndices = np.array([np.where(Temp == self.WeightStartEndTemps[0])[0][0], -1])
                    
        return self.WeightsIndices

    def WriteRates(self):
    
        SpeciesNumber = self.HighPressureRates.shape[0]
        PressureNumber = int(self.TempSpeciesRates.shape[0] / SpeciesNumber)

        PressureShift = PressureNumber*self.HighPIndex

        for RateIndex in self.SpeciesIndex:
            
            # High Pressure:
            # if all data set are NaN then prints 0 0 0.
            if len(np.where(np.isnan ( self.HighPressureRates [self.HighPIndex] [RateIndex + 1]) == True)[0]) == len(self.HighPressureRates [self.HighPIndex] [RateIndex + 1]):
                self.file[0].write(f"{self.RatePlotDict[RateIndex]} 0.000000E+00 0 0.000000\n")
                self.file[1].write(f"{self.RatePlotDict[RateIndex]} 0.000000E+00 0 0.000000\n")
                self.file[2].write(f"{self.RatePlotDict[RateIndex]} 0.000000E+00 0 0.000000\n")

            # elif checks where NaN begins in data set and terminates data when the NaN begins. 
            elif len(np.where(np.isnan(self.HighPressureRates[self.HighPIndex][RateIndex + 1]) == True)[0]) != 0:
                self.Rates = self.HighPressureRates [self.HighPIndex] [RateIndex]    [:np.where(np.isnan( self.HighPressureRates [self.HighPIndex] [RateIndex + 1] ) == True)[0][0]]
                Temp  = self.HighPressureRates     [0]          [0]        [:np.where(np.isnan( self.HighPressureRates [self.HighPIndex] [RateIndex + 1] ) == True)[0][0]]
                ModArr, DoubleArr = PlotFits(Temp, self.Rates, self.RatePlotDict[RateIndex], 'High', self.Bimolec).ArrheniusFitPlot()
                
                self.file[0].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n"), 
                self.file[1].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")
                self.file[2].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")

            # else considers whole data.
            else:
                self.Rates = self.HighPressureRates[self.HighPIndex][RateIndex + 1]
                Temp  = self.HighPressureRates    [0]            [0]
                ModArr, DoubleArr = PlotFits(Temp, self.Rates, self.RatePlotDict[RateIndex], 'High', self.Bimolec).ArrheniusFitPlot()
                
                self.file[0].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n"), 
                self.file[1].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")      
                self.file[2].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")

            # All pressrues:
            for PressureIndex in range(PressureShift, self.PlottingPressure.shape[0] + PressureShift):
                
                if len(np.where(np.isnan( self.TempSpeciesRates [PressureIndex] [RateIndex]) == True)[0]) != 0:
                    self.Rates = self.TempSpeciesRates[PressureIndex][RateIndex]    [:np.where(np.isnan(self.TempSpeciesRates [PressureIndex] [RateIndex]) == True)[0][0]]
                    Temp  = self.TempSpeciesRates     [0]           [0]        [:np.where(np.isnan(self.TempSpeciesRates [PressureIndex] [RateIndex]) == True)[0][0]]
                else:
                    self.Rates = self.TempSpeciesRates [PressureIndex] [RateIndex]
                    Temp  = self.TempSpeciesRates     [0]            [0]

                ModArr, DoubleArr = PlotFits(Temp, self.Rates, self.RatePlotDict[RateIndex], self.PlottingPressure[PressureIndex - PressureShift], self.Bimolec).ArrheniusFitPlot()

                self.file[0].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}/") + 
                                            f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." + 
                                            f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={ModArr[2]}% Max={np.max(ModArr[3])}%.\n")

                self.file[1].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][0]} {DoubleArr[0][1]} {DoubleArr[0][2]}/")+
                                            f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." +
                                            f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={DoubleArr[2]}% Max={np.max(DoubleArr[3])}%." +
                                            f"\nPLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][3]} {DoubleArr[0][4]} {DoubleArr[0][5]}/\n")  
                
                if (np.max(ModArr[3]) < self.Tolerance) or (np.max(ModArr[3]) < np.max(DoubleArr[3])):
                    self.file[2].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}/") + 
                                                f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." + 
                                                f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={ModArr[2]}% Max={np.max(ModArr[3])}%.\n") 
                else:
                    self.file[2].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][0]} {DoubleArr[0][1]} {DoubleArr[0][2]}/")+
                                        f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." +
                                        f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={DoubleArr[2]}% Max={np.max(DoubleArr[3])}%." +
                                        f"\nPLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][3]} {DoubleArr[0][4]} {DoubleArr[0][5]}/\n")

            self.file[0].write('\n')
            self.file[1].write('\n')
            self.file[2].write('\n')

    # ------------------------------------------------------------------------------------------------------ Weighted  --------------------------------------------------------------------------------------------------------------- #

    def WriteWeightedRates(self):

        SpeciesNumber = self.HighPressureRates.shape[0]
        PressureNumber = int(self.TempSpeciesRates.shape[0] / SpeciesNumber)

        PressureShift = PressureNumber*self.HighPIndex

        for RateIndex in self.SpeciesIndex:
            
            # High Pressure:
            # if all data set are NaN then prints 0 0 0.
            if len(np.where(np.isnan ( self.HighPressureRates [self.HighPIndex] [RateIndex + 1]) == True)[0]) == len(self.HighPressureRates [self.HighPIndex] [RateIndex + 1]):
                self.file[0].write(f"{self.RatePlotDict[RateIndex]} 0.000000E+00 0 0.000000\n")
                self.file[1].write(f"{self.RatePlotDict[RateIndex]} 0.000000E+00 0 0.000000\n")
                self.file[2].write(f"{self.RatePlotDict[RateIndex]} 0.000000E+00 0 0.000000\n")

            # elif checks where NaN begins in data set and terminates data when the NaN begins. 
            elif len(np.where(np.isnan(self.HighPressureRates[self.HighPIndex][RateIndex + 1]) == True)[0]) != 0:
                self.Rates = self.HighPressureRates [self.HighPIndex] [RateIndex]    [:np.where(np.isnan( self.HighPressureRates [self.HighPIndex] [RateIndex + 1] ) == True)[0][0]]
                Temp  = self.HighPressureRates     [0]          [0]        [:np.where(np.isnan( self.HighPressureRates [self.HighPIndex] [RateIndex + 1] ) == True)[0][0]]
                self.WeightsIndices = self.__CheckTemps(Temp)
                ModArr, DoubleArr = PlotFits(Temp, self.Rates, self.RatePlotDict[RateIndex], 'High', self.WeightsIndices, self.Bimolec).WeightedArrheniusFitPlot()
                self.file[0].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")
                self.file[1].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")
                self.file[2].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")
            
            # else considers whole data.
            else:
                self.Rates = self.HighPressureRates[self.HighPIndex][RateIndex + 1]
                Temp  = self.HighPressureRates    [0]            [0]
                self.WeightsIndices = self.__CheckTemps(Temp)
                ModArr, DoubleArr = PlotFits(Temp, self.Rates, self.RatePlotDict[RateIndex], 'High', self.WeightsIndices, self.Bimolec).WeightedArrheniusFitPlot()
                self.file[0].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")
                self.file[1].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")    
                self.file[2].write(f"{self.RatePlotDict[RateIndex]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}\n")  
                
            # All pressrues:
            for PressureIndex in range(PressureShift, self.PlottingPressure.shape[0] + PressureShift):
                
                if len(np.where(np.isnan( self.TempSpeciesRates [PressureIndex] [RateIndex]) == True)[0]) != 0:
                    self.Rates = self.TempSpeciesRates[PressureIndex][RateIndex]    [:np.where(np.isnan(self.TempSpeciesRates [PressureIndex] [RateIndex]) == True)[0][0]]
                    Temp  = self.TempSpeciesRates     [0]           [0]        [:np.where(np.isnan(self.TempSpeciesRates [PressureIndex] [RateIndex]) == True)[0][0]]
                else:
                    self.Rates = self.TempSpeciesRates [PressureIndex] [RateIndex]
                    Temp  = self.TempSpeciesRates     [0]            [0]
                
                self.WeightsIndices = self.__CheckTemps(Temp)
                ModArr, DoubleArr = PlotFits(Temp, self.Rates, self.RatePlotDict[RateIndex], self.PlottingPressure[PressureIndex - PressureShift], self.WeightsIndices, self.Bimolec).WeightedArrheniusFitPlot()
                
                if self.WeightsIndices[1] != -1:
                    self.WeightsIndices[1] -= 1

                self.file[0].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}/") + 
                                            f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." + 
                                            f" Fit tailored for {Temp[self.WeightsIndices[0]]}K - {Temp[self.WeightsIndices[1]]}K MAPE={ModArr[4]}% Max={np.max(ModArr[5])}%." + 
                                            f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={ModArr[2]}% Max={np.max(ModArr[3])}%.\n")

                self.file[1].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][0]} {DoubleArr[0][1]} {DoubleArr[0][2]}/")+
                                            f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." +
                                            f" Fit tailored for {Temp[self.WeightsIndices[0]]}K - {Temp[self.WeightsIndices[1]]}K MAPE={DoubleArr[4]}% Max={np.max(DoubleArr[5])}%." +
                                            f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={DoubleArr[2]}% Max={np.max(DoubleArr[3])}%." +
                                            f"\nPLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][3]} {DoubleArr[0][4]} {DoubleArr[0][5]}/\n")  
                
                if (np.max(ModArr[5]) < self.Tolerance) or (np.max(ModArr[5]) < np.max(DoubleArr[5])):
                    self.file[2].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {ModArr[0][0]} {ModArr[0][1]} {ModArr[0][2]}/") + 
                                                f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." + 
                                                f" Fit tailored for {Temp[self.WeightsIndices[0]]}K - {Temp[self.WeightsIndices[1]]}K MAPE={ModArr[4]}% Max={np.max(ModArr[5])}%." + 
                                                f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={ModArr[2]}% Max={np.max(ModArr[3])}%.\n")  
                else:
                    self.file[2].write("{:<60}".format(f"PLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][0]} {DoubleArr[0][1]} {DoubleArr[0][2]}/")+
                                        f"!Labbe Lab {datetime.datetime.now().strftime('%B %d, %Y')}." +
                                        f" Fit tailored for {Temp[self.WeightsIndices[0]]}K - {Temp[self.WeightsIndices[1]]}K MAPE={DoubleArr[4]}% Max={np.max(DoubleArr[5])}%." +
                                        f" Total data range {np.min(Temp)}K - {np.max(Temp)}K MAPE={DoubleArr[2]}% Max={np.max(DoubleArr[3])}%." +
                                        f"\nPLOG/ {self.PlottingPressure[PressureIndex - PressureShift]} {DoubleArr[0][3]} {DoubleArr[0][4]} {DoubleArr[0][5]}/\n")

            self.file[0].write('\n')
            self.file[1].write('\n')
            self.file[2].write('\n')

if __name__ == '__main__':
    print(f"{os.path.basename(__file__)} is built as a module and has no functionality to run on its own")