#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 07:44:54 2021

@ Original AMID Author: Marc M. E. Cormier

@ Current AMIDR Author: Mitchell Ball

"""

import pandas as pd
import numpy as np
import sys
from scipy.optimize import curve_fit, fsolve
from scipy import stats
from pathlib import Path
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings(action = 'ignore')

plt.rc('lines', markersize = 4, linewidth = 0.75)
plt.rc('axes', grid = True, labelsize = 14)
plt.rc('axes.grid', which = 'both')
plt.rc('grid', color = 'grey')
plt.rc('xtick.minor', bottom = True, top = False)
plt.rc('ytick.minor', left = True, right = False)
plt.rc('xtick', bottom = True, top = False, direction = 'out', labelsize = 10)
plt.rc('ytick', left = True, right = False, direction = 'out', labelsize = 10)
plt.rc('legend', frameon = False, fontsize = 10, columnspacing = 1.0, handletextpad = 0.5, handlelength = 1.4)
plt.rc('errorbar', capsize = 2)
plt.rc('savefig', dpi = 300)

RATES = np.array([0.01, 0.05, 0.1, 0.2, 1/3, 0.5, 1, 2, 2.5, 5, 10, 20, 40, 80, 160, 320, 640, 1280])

COLUMNS = ['Time', 'Cycle', 'Step', 'Current', 'Potential', 'Capacity', 'Prot_step']
UNITS = ['(h)', None, None, '(mA)', '(V)', '(mAh)', None]

SHAPES = ['sphere', 'plane']

class BIOCONVERT():
    
    def __init__(self, path, form_files, d_files, c_files, cellname, export_data = True, export_fig = True):
                
        # Acquire header info from first file
        all_files = []
        all_files.extend(form_files)
        all_files.extend(d_files)
        all_files.extend(c_files)
        firstFileLoc = Path(path) / all_files[0]
        with open(firstFileLoc, 'r') as f:
            
            # Read beginning of first file to discover lines in header
            f.readline()
            hlinenum = int(f.readline().strip().split()[-1])
            header = f.readlines()[:hlinenum-3]
            
            # Acquire protocol name, capacity, active mass, and time started
            protText = re.search('Loaded Setting File : (.+)', ''.join(header))
            protName = protText.group(1).strip()
            
            massText = re.search('Mass of active material : (\d+.?\d+) (.+)', ''.join(header))
            massVal = float(massText.group(1))
            massUnit = massText.group(2).strip()
            if massUnit != 'mg':
                print('Mass Unit: ' + massUnit)
                print('Please edit first file to express mass in mg so that specific capacity is accurately calculated')
            
            capacityText = re.search('Battery capacity : (\d+.?\d+) (.+)', ''.join(header))
            capacityVal = float(capacityText.group(1))
            capacityUnit = capacityText.group(2).strip()
            if capacityUnit != 'mA.h':
                print('Capacity Unit: ' + capacityUnit + '100')
                print('Please edit first file to express capacity in mA.h so that specific capacity and rates are accurately calculated')
            capacityUnit = capacityUnit.replace('.h', 'Hr')
            
            startText = re.search('Technique started on : (.+)', ''.join(header))
            startTime = startText.group(1).strip()
            
            # Write header text
            csvHeader = '[Summary]\nCell: ' + cellname + '\nFirst Protocol: ' + protName \
            + '\nMass (' + massUnit + '): ' + str(massVal) + '\nCapacity (' + capacityUnit + '): ' + str(capacityVal) \
            + '\nStarted: ' + startTime + '\n[End Summary]\n[Data]\n'
        
        # Generate complete csv
        df = pd.DataFrame({})
        
        # Generate form dataframe if data available
        if form_files:
            dfForm = pd.DataFrame({})
            
            # Read and combine form file data
            for f in form_files:
                formFileLoc = Path(path) / f
                with open(formFileLoc, 'r') as f:
            
                    # Read beginning of form file to discover lines in header
                    f.readline()
                    hlinenum = int(f.readline().strip().split()[-1]) - 1
                    
                # Read file into dataframe and convert to UHPC format
                dfTempForm = pd.read_csv(formFileLoc, skiprows = hlinenum, sep = '\t', encoding_errors = 'replace')
                dfTempForm = dfTempForm[['mode', 'time/s', 'I/mA', 'Ewe-Ece/V', 'Ewe/V', '(Q-Qo)/mA.h', 'Ns']]
                
                # Convert to NVX initial step convention
                dfTempForm['Ns'] = dfTempForm['Ns'] + 1
                
                # Add last previous capacity, time, and step number to current data
                if not(dfForm.empty):
                    dfTempForm['time/s'] = dfTempForm['time/s'] + dfForm['time/s'].iat[-1]
                    dfTempForm['(Q-Qo)/mA.h'] = dfTempForm['(Q-Qo)/mA.h'] + dfForm['(Q-Qo)/mA.h'].iat[-1]
                    dfTempForm['Ns'] = dfTempForm['Ns'] + dfForm['Ns'].iat[-1]
                
                # Concatenate
                dfForm = pd.concat([dfForm, dfTempForm])
            
            # Convert to hours, base units, NVX labels, and NVX rest step convention while retaining order
            dfForm['time/s'] = dfForm['time/s'] / 3600
            dfForm['I/mA'] = dfForm['I/mA'] / 1000
            dfForm['(Q-Qo)/mA.h'] = dfForm['(Q-Qo)/mA.h'] / 1000
            dfForm['mode'] = dfForm['mode'].replace(3, 0)
            
            dfForm.rename(columns = {'mode':'Step Type', 
                                   'time/s':'Run Time (h)', 
                                   'I/mA':'Current (A)', 
                                   'Ewe-Ece/V':'Potential vs. Counter (V)', 
                                   'Ewe/V':'Potential (V)', 
                                   '(Q-Qo)/mA.h':'Capacity (Ah)', 
                                   'Ns':'Step Number'}, 
                          inplace = True)
            
            # Generate form file
            if export_data:
                pathFileForm = Path(path) / (cellname + ' Form.csv')
                with open(pathFileForm, 'w') as f:
                    f.write(csvHeader)
                
                print("Formation data exporting to:\n" + str(pathFileForm))
                dfForm.to_csv(pathFileForm, mode = 'a', index = False)
            
            # Concatenate
            df = pd.concat([df, dfForm])
        
        # Generate D dataframe if data available
        if d_files:
            dfD = pd.DataFrame({})
            
            # Read, V average, and combine d file data
            for file in d_files:
                dFileLoc = Path(path) / file
                with open(dFileLoc, 'r') as f:
            
                    # Read beginning of d file to discover lines in header
                    f.readline()
                    hlinenum = int(f.readline().strip().split()[-1]) - 1
                    
                # Read file into dataframe and convert to UHPC format
                dfTempD = pd.read_csv(dFileLoc, skiprows = hlinenum, sep = '\t', encoding_errors = 'replace')
                dfTempD = dfTempD[['mode', 'time/s', 'I/mA', 'Ewe-Ece/V', 'Ewe/V', '(Q-Qo)/mA.h', 'Ns', 'control/mA']]
    
                # Convert 0A CC to NVX rest steps and trim off control I column
                dfTempD['mode'].mask(dfTempD['control/mA'] == 0, 0, inplace = True)
                dfTempD['mode'].mask(dfTempD['mode'] == 3, 0, inplace = True)
                dfTempD.drop(columns = ['control/mA'], inplace = True)
 
                dfTempDSteps = dfTempD.drop_duplicates(subset = ['Ns'], ignore_index = True)
                
                # Detect if test ended prematurely
                if dfTempDSteps['mode'].iloc[-1] != 0:
                    
                    # Add dummy step at end to prevent overindexing
                    dummyD = dfTempD.loc[dfTempD.index[-1]:dfTempD.index[-1]].copy()
                    dummyD['Ns'] = dummyD['Ns'] + 1
                    dummyD['mode'] = 0
                    dfTempD = pd.concat([dfTempD, dummyD], ignore_index = True)
                    dfTempDSteps = pd.concat([dfTempDSteps, dummyD], ignore_index = True)
                    
                # Iterate over each pulse starting with their preceeding OCV V rest step    
                for i in dfTempDSteps.index:
                    if i != dfTempDSteps.index[-1]:
                        if dfTempDSteps['mode'][i] == 0 and dfTempDSteps['mode'][i+1] != 0:
                            
                            # Average together all points of the rest step before a pulse into 1 point (OCV V)
                            ocvFinSel = dfTempD['Ns'] == dfTempDSteps['Ns'][i]
                            ocvFinVals = dfTempD[ocvFinSel].mean(axis = 0)
                            dfTempD[ocvFinSel] = ocvFinVals
                
                            # Determine nAvg, the number of datapoints to average together so that there are 10 points in the first step
                            nAvg = int(sum(dfTempD['Ns'] == dfTempDSteps['Ns'][i+1])/10)
                            
                            # Iterate over each CC step in pulse
                            j = 1
                            while dfTempDSteps['mode'][i+j] != 0:
                                pulseInd = dfTempD[dfTempD['Ns'] == dfTempDSteps['Ns'][i+j]].index[0]
                                
                                # Give OCV V to first point in pulse else remove first point in CC step
                                if j == 1:
                                    dfTempD['Ewe-Ece/V'][pulseInd] = ocvFinVals['Ewe-Ece/V']
                                    dfTempD['Ewe/V'][pulseInd] = ocvFinVals['Ewe/V']
                                else:
                                    dfTempD.drop([pulseInd], inplace = True)
                                    
                                    # Prints step and indice of datapoint removed [Default Commented Out]
                                    #print(dfTempDSteps['Ns'][i+j], pulseInd)
                                
                                # Iterate over sets of nAvg datapoints within a CC step skipping the first and remainder datapoints
                                if pulseInd + nAvg <= dfTempD.index[-1]:
                                    while dfTempD['Ns'][pulseInd + nAvg] == dfTempDSteps['Ns'][i+j]:
                                        
                                        # Average together nAvg points into 1 point
                                        CCPointVals = dfTempD.loc[pulseInd + 1:pulseInd + nAvg].mean(axis = 0)
                                        dfTempD.loc[pulseInd + 1:pulseInd + nAvg] = CCPointVals.values
                                        pulseInd = pulseInd + nAvg 
                                        
                                        if pulseInd + nAvg > dfTempD.index[-1]:
                                            break
                                    
                                # Drop all remainder datapoints
                                nextStepInd = dfTempD[dfTempD['Ns'] == dfTempDSteps['Ns'][i+j+1]].index[0]
                                dfTempD.drop(dfTempD.loc[pulseInd + 1:nextStepInd - 1].index, inplace = True)
                                
                                # Prints step and range of indices of datapoints removed [Default Commented Out]
                                #if pulseInd + 1 != nextStepInd:
                                    #print(dfTempDSteps['Ns'][i+j], pulseInd + 1, '-', nextStepInd - 1)
                                
                                j = j + 1
                    else:
                        
                        # Calculate OCV V for end of final pulse
                        if dfTempDSteps['mode'][i] == 0 and dfTempDSteps['mode'][i-1] == 0 and dfTempDSteps['mode'][i-2] == 0 and dfTempDSteps['mode'][i-3] == 1:
                            
                            # Average together all points of the rest step before a pulse into 1 point (OCV V)
                            ocvFinSel = dfTempD['Ns'] == dfTempDSteps['Ns'][i]
                            ocvFinVals = dfTempD[ocvFinSel].mean(axis = 0)
                            dfTempD[ocvFinSel] = ocvFinVals
                        
                        # Label steps after last OCV V as CC to prevent analysis (Test stopped prematurely)
                        else:
                            print(file, 'ended prematurely. Labeling last pulse as unfinished to prevent analysis.')
                            for i in range(len(dfTempDSteps.index)):
                                if dfTempDSteps['mode'].iat[-i-1] == 0:
                                    dfTempDSteps['mode'].iat[-i-1] = 1
                                    dfTempD['mode'][dfTempD['Ns'] == dfTempDSteps['Ns'].iat[-i-1]] = 1
                                else:
                                    break
                        
                # Remove initial rest step series 
                for i in range(len(dfTempDSteps.index)):
                    if dfTempDSteps['mode'][i] == 0:
                        dfTempD.drop(dfTempD.loc[dfTempD['Ns'] == dfTempDSteps['Ns'][i]].index, inplace = True)
                        dfTempDSteps.drop(i, inplace = True)
                    else:
                        dfTempDSteps.reset_index(drop = True, inplace = True)
                        break
                
                # Remove duplicates to simplify to one datapoint per averaging
                dfTempD.drop_duplicates(inplace = True)
                
                # Combine all rest steps except for the last step in a series into 1 step
                for i in range(len(dfTempDSteps.index)):
                    if i != 0 and i != dfTempDSteps.index[-1]:
                        if dfTempDSteps['mode'][i] == 0 and dfTempDSteps['mode'][i-1] == 0 and dfTempDSteps['mode'][i+1] == 0:
                            dfTempD['Ns'][dfTempD['Ns'] == dfTempDSteps['Ns'][i]] = dfTempDSteps['Ns'][i-1]
                            dfTempDSteps['Ns'][i] = dfTempDSteps['Ns'][i-1]
                
                # Combine all CC steps in a series
                for i in range(len(dfTempDSteps.index)):
                    if i != 0:
                        if dfTempDSteps['mode'][i] != 0 and dfTempDSteps['mode'][i-1] != 0:
                            dfTempD['Ns'][dfTempD['Ns'] == dfTempDSteps['Ns'][i]] = dfTempDSteps['Ns'][i-1]
                            dfTempDSteps['Ns'][i] = dfTempDSteps['Ns'][i-1]
                            
                # Label last step as rest step if not already (Test stopped prematurely)
                dfTempD['mode'].iloc[-1] = 0
                dfTempDSteps['mode'].iloc[-1] = 0
                
                # Relabel steps with continuous integers starting from 1
                newDSteps = dfTempD.drop_duplicates(subset = ['Ns'], ignore_index = True)
                dfTempD['Ns'].replace(newDSteps['Ns'].values, newDSteps.index+1, inplace = True)
                dfTempDSteps['Ns'].replace(newDSteps['Ns'].values, newDSteps.index+1, inplace = True)
                
                # Prints dataset after all transformations [Default Commented Out]
                #print(dfTempD)
                #print(dfTempD.loc[9500:13570])
                #print(dfTempDSteps[0:50], dfTempDSteps[50:100], dfTempDSteps[100:150])
                
                # Add last previous capacity, time, and step number to current data
                if not(dfD.empty):
                    dfTempD['time/s'] = dfTempD['time/s'] + dfD['time/s'].iat[-1]
                    dfTempD['(Q-Qo)/mA.h'] = dfTempD['(Q-Qo)/mA.h'] + dfD['(Q-Qo)/mA.h'].iat[-1]
                    dfTempD['Ns'] = dfTempD['Ns'] + dfD['Ns'].iat[-1]
                
                # Concatenate 
                dfD = pd.concat([dfD, dfTempD])
                
            # Convert to hours, base units, and NVX labels while retaining order
            dfD['time/s'] = dfD['time/s'] / 3600
            dfD['I/mA'] = dfD['I/mA'] / 1000
            dfD['(Q-Qo)/mA.h'] = dfD['(Q-Qo)/mA.h'] / 1000
            
            dfD.rename(columns = {'mode':'Step Type', 
                                'time/s':'Run Time (h)', 
                                'I/mA':'Current (A)', 
                                'Ewe-Ece/V':'Potential vs. Counter (V)', 
                                'Ewe/V':'Potential (V)', 
                                '(Q-Qo)/mA.h':'Capacity (Ah)', 
                                'Ns':'Step Number'}, 
                       inplace = True)

            # Add last capacity to output file
            if not(df.empty):
                dfD['Capacity (Ah)'] = dfD['Capacity (Ah)'] + df['Capacity (Ah)'].iat[-1]
                
            # Generate D File
            if export_data:
                pathFileD = Path(path) / (cellname + ' Discharge.csv')
                with open(pathFileD, 'w') as f:
                    f.write(csvHeader)      
                    
                print("Discharge data exporting to:\n" + str(pathFileD))
                dfD.to_csv(pathFileD, mode = 'a', index = False)
                
            # Add last time, and step number to graphs
            if not(df.empty):
                dfD['Run Time (h)'] = dfD['Run Time (h)'] + df['Run Time (h)'].iat[-1]
                dfD['Step Number'] = dfD['Step Number'] + df['Step Number'].iat[-1]
                    
            # Concatenate
            df = pd.concat([df, dfD])
            
        # Generate C dataframe if data available
        if c_files:
            dfC = pd.DataFrame({})
            
            # Read, V average, and combine c file data
            for file in c_files:
                cFileLoc = Path(path) / file
                with open(cFileLoc, 'r') as f:
            
                    # Read beginning of c file to discover lines in header
                    f.readline()
                    hlinenum = int(f.readline().strip().split()[-1]) - 1
                    
                # Read file into dataframe and convert to UHPC format
                dfTempC = pd.read_csv(cFileLoc, skiprows = hlinenum, sep = '\t', encoding_errors = 'replace')
                dfTempC = dfTempC[['mode', 'time/s', 'I/mA', 'Ewe-Ece/V', 'Ewe/V', '(Q-Qo)/mA.h', 'Ns', 'control/mA']]
    
                # Convert 0A CC to NVX rest steps and trim off control I column
                dfTempC['mode'].mask(dfTempC['control/mA'] == 0, 0, inplace = True)
                dfTempC.drop(columns = ['control/mA'], inplace = True)
                
                dfTempCSteps = dfTempC.drop_duplicates(subset = ['Ns'], ignore_index = True)     
                
                # Detect if test ended prematurely
                if dfTempCSteps['mode'].iloc[-1] != 0:
                    
                    # Add dummy step at end to prevent overindexing
                    dummyC = dfTempC.loc[dfTempC.index[-1]:dfTempC.index[-1]].copy()
                    dummyC['Ns'] = dummyC['Ns'] + 1
                    dummyC['mode'] = 0
                    dfTempC = pd.concat([dfTempC, dummyC], ignore_index = True)
                    dfTempCSteps = pd.concat([dfTempCSteps, dummyC], ignore_index = True)
                
                # Iterate over each pulse starting with their preceeding OCV V rest step 
                for i in dfTempCSteps.index:
                    if i != dfTempCSteps.index[-1]:
                        if dfTempCSteps['mode'][i] == 0 and dfTempCSteps['mode'][i+1] != 0:
                            
                            # Average together all points of the rest step before a pulse into 1 point (OCV V)
                            ocvFinSel = dfTempC['Ns'] == dfTempCSteps['Ns'][i]
                            ocvFinVals = dfTempC[ocvFinSel].mean(axis = 0)
                            dfTempC[ocvFinSel] = ocvFinVals
                
                            # Determine nAvg, the number of datapoints to average together so that there are 10 points in the first step
                            nAvg = int(sum(dfTempC['Ns'] == dfTempCSteps['Ns'][i+1])/10)
                            
                            # Iterate over each CC step in pulse
                            j = 1
                            while dfTempCSteps['mode'][i+j] != 0:
                                pulseInd = dfTempC[dfTempC['Ns'] == dfTempCSteps['Ns'][i+j]].index[0]
                                
                                # Give OCV V to first point in pulse else remove first point in CC step
                                if j == 1:
                                    dfTempC['Ewe-Ece/V'][pulseInd] = ocvFinVals['Ewe-Ece/V']
                                    dfTempC['Ewe/V'][pulseInd] = ocvFinVals['Ewe/V']
                                else:
                                    dfTempC.drop([pulseInd], inplace = True)
                                    
                                    # Prints step and indice of datapoint removed [Default Commented Out]
                                    #print(dfTempCSteps['Ns'][i+j], pulseInd)
                                
                                # Iterate over sets of nAvg datapoints within a CC step skipping the first and remainder datapoints
                                if pulseInd + nAvg <= dfTempC.index[-1]:
                                    while dfTempC['Ns'][pulseInd + nAvg] == dfTempCSteps['Ns'][i+j]:
                                        
                                        # Average together nAvg points into 1 point
                                        CCPointVals = dfTempC.loc[pulseInd + 1:pulseInd + nAvg].mean(axis = 0)
                                        dfTempC.loc[pulseInd + 1:pulseInd + nAvg] = CCPointVals.values
                                        pulseInd = pulseInd + nAvg 
                                        
                                        if pulseInd + nAvg > dfTempC.index[-1]:
                                            break
                                    
                                # Drop all remainder datapoints
                                nextStepInd = dfTempC[dfTempC['Ns'] == dfTempCSteps['Ns'][i+j+1]].index[0]
                                dfTempC.drop(dfTempC.loc[pulseInd + 1:nextStepInd - 1].index, inplace = True)
                                
                                # Prints step and range of indices of datapoints removed [Default Commented Out]
                                #if pulseInd + 1 != nextStepInd:
                                    #print(dfTempCSteps['Ns'][i+j], pulseInd + 1, '-', nextStepInd - 1)
                                
                                j = j + 1
                    else:
                        
                        # Calculate OCV V for end of final pulse
                        if dfTempCSteps['mode'][i] == 0 and dfTempCSteps['mode'][i-1] == 0 and dfTempCSteps['mode'][i-2] == 0 and dfTempCSteps['mode'][i-3] == 1:
                            
                            # Average together all points of the rest step before a pulse into 1 point (OCV V)
                            ocvFinSel = dfTempC['Ns'] == dfTempCSteps['Ns'][i]
                            ocvFinVals = dfTempC[ocvFinSel].mean(axis = 0)
                            dfTempC[ocvFinSel] = ocvFinVals
                        
                        # Label steps after last OCV V as CC to prevent analysis (Test stopped prematurely)
                        else:
                            print(file, 'ended prematurely. Labeling last pulse as unfinished to prevent analysis.')
                            for i in range(len(dfTempCSteps.index)):
                                if dfTempCSteps['mode'].iat[-i-1] == 0:
                                    dfTempCSteps['mode'].iat[-i-1] = 1
                                    dfTempC['mode'][dfTempC['Ns'] == dfTempCSteps['Ns'].iat[-i-1]] = 1
                                else:
                                    break
                        
                # Remove initial rest step series 
                for i in range(len(dfTempCSteps.index)):
                    if dfTempCSteps['mode'][i] == 0:
                        dfTempC.drop(dfTempC.loc[dfTempC['Ns'] == dfTempCSteps['Ns'][i]].index, inplace = True)
                        dfTempCSteps.drop(i, inplace = True)
                    else:
                        dfTempCSteps.reset_index(drop = True, inplace = True)
                        break
                
                # Remove duplicates to simplify to one datapoint per averaging
                dfTempC.drop_duplicates(inplace = True)
                
                # Combine all rest steps in a series except for the last step into 1 step
                for i in range(len(dfTempCSteps.index)):
                    if i != 0 and i != dfTempCSteps.index[-1]:
                        if dfTempCSteps['mode'][i] == 0 and dfTempCSteps['mode'][i-1] == 0 and dfTempCSteps['mode'][i+1] == 0:
                            dfTempC['Ns'][dfTempC['Ns'] == dfTempCSteps['Ns'][i]] = dfTempCSteps['Ns'][i-1]
                            dfTempCSteps['Ns'][i] = dfTempCSteps['Ns'][i-1]
                
                # Combine all CC steps in a series
                for i in range(len(dfTempCSteps.index)):
                    if i != 0:
                        if dfTempCSteps['mode'][i] != 0 and dfTempCSteps['mode'][i-1] != 0:
                            dfTempC['Ns'][dfTempC['Ns'] == dfTempCSteps['Ns'][i]] = dfTempCSteps['Ns'][i-1]
                            dfTempCSteps['Ns'][i] = dfTempCSteps['Ns'][i-1]
                
                # Label last step as rest step if not already (Test stopped prematurely)
                dfTempC['mode'].iloc[-1] = 0
                dfTempCSteps['mode'].iloc[-1] = 0
                
                # Relabel steps with continuous integers starting from 1
                newCSteps = dfTempC.drop_duplicates(subset = ['Ns'], ignore_index = True)
                dfTempC['Ns'].replace(newCSteps['Ns'].values, newCSteps.index+1, inplace = True)
                dfTempCSteps['Ns'].replace(newCSteps['Ns'].values, newCSteps.index+1, inplace = True)
                
                # Add last previous capacity, time, and step number to current data
                if not(dfC.empty):
                    dfTempC['time/s'] = dfTempC['time/s'] + dfC['time/s'].iat[-1]
                    dfTempC['(Q-Qo)/mA.h'] = dfTempC['(Q-Qo)/mA.h'] + dfC['(Q-Qo)/mA.h'].iat[-1]
                    dfTempC['Ns'] = dfTempC['Ns'] + dfC['Ns'].iat[-1]
                
                # Concatenate 
                dfC = pd.concat([dfC, dfTempC])
                
            # Convert to hours, base units, and NVX labels while retaining order
            dfC['time/s'] = dfC['time/s'] / 3600
            dfC['I/mA'] = dfC['I/mA'] / 1000
            dfC['(Q-Qo)/mA.h'] = dfC['(Q-Qo)/mA.h'] / 1000
            
            dfC.rename(columns = {'mode':'Step Type', 
                                'time/s':'Run Time (h)', 
                                'I/mA':'Current (A)', 
                                'Ewe-Ece/V':'Potential vs. Counter (V)', 
                                'Ewe/V':'Potential (V)', 
                                '(Q-Qo)/mA.h':'Capacity (Ah)', 
                                'Ns':'Step Number'}, 
                       inplace = True)
            
            # Add last capacity to output file
            if not(df.empty):
                dfC['Capacity (Ah)'] = dfC['Capacity (Ah)'] + df['Capacity (Ah)'].iat[-1]
    
            # Generate C file
            if export_data:
                pathFileC = Path(path) / (cellname + ' Charge.csv')
                with open(pathFileC, 'w') as f:
                    f.write(csvHeader)
                            
                print("Charge data exporting to:\n" + str(pathFileC))
                dfC.to_csv(pathFileC, mode = 'a', index = False)
            
            # Add last time, and step number to graphs
            if not(df.empty):
                dfC['Run Time (h)'] = dfC['Run Time (h)'] + df['Run Time (h)'].iat[-1]
                dfC['Step Number'] = dfC['Step Number'] + df['Step Number'].iat[-1]
                    
            # Concatenate
            df = pd.concat([df, dfC])
        
        # Generate full file
        if export_data:
            pathFile = Path(path) / (cellname + ' All.csv')
            with open(pathFile, 'w') as f:
                f.write(csvHeader)
                
            df.to_csv(pathFile, mode = 'a', index = False)
        
        # Generate complete graph
        fig, axs = plt.subplots(nrows = 1, ncols = 2, sharey = True, figsize = (6, 3), gridspec_kw = {'wspace':0.0})
        
        if form_files:
            axs[0].plot(dfForm['Run Time (h)'], dfForm['Potential vs. Counter (V)'], 'k--', label = 'Formation vs. $E_c$')
            axs[0].plot(dfForm['Run Time (h)'], dfForm['Potential (V)'], 'k-', label = 'Formation vs. $E_r$')
        
        if d_files:
            axs[0].plot(dfD['Run Time (h)'], dfD['Potential vs. Counter (V)'], 'r--', label = 'Discharge vs. $E_c$')
            axs[0].plot(dfD['Run Time (h)'], dfD['Potential (V)'], 'r-', label = 'Discharge vs. $E_r$')
        
        if c_files:
            axs[0].plot(dfC['Run Time (h)'], dfC['Potential vs. Counter (V)'], 'b--', label = 'Charge vs. $E_c$')
            axs[0].plot(dfC['Run Time (h)'], dfC['Potential (V)'], 'b-', label = 'Charge vs. $E_r$')

        axs[0].set_xlabel('Time (h)')
        axs[0].set_ylabel('Voltage (V)')
        axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0].grid(which = 'minor', color = 'lightgrey')
        
        if form_files:
            axs[1].plot(dfForm['Capacity (Ah)']*1000000/massVal, dfForm['Potential vs. Counter (V)'], 'k--', label = 'Formation vs. $E_c$')
            axs[1].plot(dfForm['Capacity (Ah)']*1000000/massVal, dfForm['Potential (V)'], 'k-', label = 'Formation vs. $E_r$')
        
        if d_files:
            axs[1].plot(dfD['Capacity (Ah)']*1000000/massVal, dfD['Potential vs. Counter (V)'], 'r--', label = 'Discharge vs. $E_c$')
            axs[1].plot(dfD['Capacity (Ah)']*1000000/massVal, dfD['Potential (V)'], 'r-', label = 'Discharge vs. $E_r$')
        
        if c_files:
            axs[1].plot(dfC['Capacity (Ah)']*1000000/massVal, dfC['Potential vs. Counter (V)'], 'b--', label = 'Charge vs. $E_c$')
            axs[1].plot(dfC['Capacity (Ah)']*1000000/massVal, dfC['Potential (V)'], 'b-', label = 'Charge vs. $E_r$')
        
        axs[1].set_xlabel('Capacity (mAh/g)')
        axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1].grid(which = 'minor', color = 'lightgrey')
        
        plt.legend(bbox_to_anchor = (1.0, 0.5), loc = 'center left')
        
        if export_fig:
            figname = Path(path) / '{} Protocol.jpg'.format(cellname)
            print(figname)
            plt.savefig(figname, bbox_inches = 'tight')
        
        plt.show()
        plt.close()
        
        if d_files and c_files:
            dfDOCV = dfD[dfD['Step Type'] != 0].drop_duplicates(['Step Number'])
            dfCOCV = dfC[dfC['Step Type'] != 0].drop_duplicates(['Step Number'])
            
            fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize = (6, 3), gridspec_kw = {'wspace':0.0})
            axs.plot(dfDOCV['Capacity (Ah)']*1000000/massVal, dfDOCV['Potential (V)'], 'r.-', label = 'Discharge Relaxed vs. $E_r$')
            axs.plot(dfCOCV['Capacity (Ah)']*1000000/massVal, dfCOCV['Potential (V)'], 'b.-', label = 'Charge Relaxed vs. $E_r$')
            axs.set_xlabel('Capacity (mAh/g)')
            axs.set_ylabel('Voltage (V)')
            axs.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs.grid(which = 'minor', color = 'lightgrey')
            plt.legend(frameon = True)
            
            if export_fig:
                figname = Path(path) / '{} Relax Match.jpg'.format(cellname)
                print(figname)
                plt.savefig(figname, bbox_inches = 'tight')
            
            plt.show()
            plt.close()
        
class AMIDR():
    
    def __init__(self, path, uhpc_file, single_pulse, export_data = True, export_fig = None, use_input_cap = True, 
                 capacitance_corr = False, fcap_min = 0.0, spliced = False, force2e = False, parselabel = None):
        
        if not(export_fig is None): ('There is no figure to export. Feel free to neglect this argument.')
        
        self.single_p = single_pulse
        self.capacitance_corr = capacitance_corr
        self.fcap_min = fcap_min
        if parselabel is None:
            self.cell_label = ('.').join(uhpc_file.split('.')[0:-1])
        else:
            self.cell_label = ('.').join(uhpc_file.split('.')[0:-1]) + '-' + parselabel
        print(self.cell_label)
        self.dst = Path(path) / self.cell_label
        # If does not exist, create dir.
        if self.dst.is_dir() is False and export_data is True:
            self.dst.mkdir()
            print('Create directory: {}'.format(self.dst))
        self.src = Path(path)
        
        self.uhpc_file = self.src / uhpc_file
             
        with open(self.uhpc_file, 'r') as f:
            lines = f.readlines()
        nlines = len(lines)
        headlines = []
        for i in range(nlines):
            headlines.append(lines[i])
            l = lines[i].strip().split()
            if l[0][:6] == '[Data]':
                nskip = i+1
                break
        
        header = ''.join(headlines)
        del lines
                
        # find mass and theoretical cap using re on header str
        m = re.search('Mass\s+\(.*\):\s+(\d+)?\.\d+', header)
        m = m.group(0).split()
        mass_units = m[1][1:-2]
        if mass_units == 'mg':
            self.mass = float(m[-1]) / 1000
        else:
            self.mass = float(m[-1])
        
        m = re.search('Capacity\s+(.*):\s+(\d+)?\.\d+', header)
        m = m.group(0).split()
        cap_units = m[1][1:-2]
        if cap_units == 'mAHr':
            self.input_cap = float(m[-1]) / 1000
        else:
            self.input_cap = float(m[-1])
            
        m = re.search('Cell: .+?(?=,|\\n)', header)
        m = m.group(0).split()
        self.cellname = " ".join(m[1:])
        
        #self.cellname = headlines[1][-1]
        #self.mass = float(headlines[4][-1]) / 1000
        #self.input_cap = float(headlines[5][-1]) / 1000  # Convert to Ah
        #if headlines[10][0] == '[Data]':
        #    hlinenum = 11
            #hline = f.readline()
        #    hline = headlines[10]
        #else:
        #    hlinenum = 12
            #f.readline()
            #hline = f.readline()
        #    hline = headlines[11]
               
        print('Working on cell: {}'.format(self.cellname))
        print('Positive electrode active mass: {} g'.format(self.mass))
        print('Input cell capacity: {} Ah'.format(round(self.input_cap, 10)))
        
        self.df = pd.read_csv(self.uhpc_file, header = nskip)
        
        self.df.rename(columns = {'Capacity (Ah)': 'Capacity', 
                                  'Potential (V)': 'Potential', 
                                  'Potential vs. Counter (V)':'Label Potential', 
                                  'Run Time (h)': 'Time', 
                                  'Time (h)': 'Time', 
                                  'Current (A)': 'Current', 
                                  'Cycle Number': 'Cycle', 
                                  'Meas I (A)': 'Current', 
                                  'Step Type': 'Step', 
                                  'Prot.Step': 'Prot_step', 
                                  'Step Number': 'Prot_step'}, 
                                   inplace = True)
        #print(self.df.columns)
        #print(self.df.Step.unique())
        
        if single_pulse == True and spliced == True:
            sys.exit("single_pulse cannot operate on spliced files. Manually clean up your spliced file and select spliced = false")
        
        # Add Prot_step column if column does not yet exist or spliced file is used.
        if 'Prot_step' not in self.df.columns or spliced == True:
            s = self.df.Step
            self.df['Prot_step'] = s.ne(s.shift()).cumsum() - 1

        #if hline[-4:] == 'Flag':
        #    self.df = self.df.rename(columns = {'Flag':'Prot.Step'})
        #    i = self.df.Step
        #    self.df['Prot.Step'] = i.ne(i.shift()).cumsum() - 1
            
        #self.df.columns = COLUMNS
        
        # Adjust data where time is not monotonically increasing.   
        t = self.df['Time'].values
        cap = self.df['Capacity'].values
        dt = t[1:] - t[:-1]
        inds = np.where(dt < 0.0)[0]
        if len(inds) > 0:
            print('Indices being adjusted due to time non-monotonicity: {}'.format(inds))
            self.df['Time'][inds+1] = (t[inds] + t[inds+2])/2
            self.df['Capacity'][inds+1] = (cap[inds] + cap[inds+2])/2
        # Adjust data where potential is negative.
        inds = self.df.index[self.df['Potential'] < 0.0].tolist()
        if len(inds) > 0:
            print('Indices being adjusted due to negative voltage: {}'.format(inds))
            self.df['Potential'][inds] = (t[inds-1] + t[inds+1])/2
        
        if 'Label Potential' not in self.df:
            self.df['Label Potential'] = self.df['Potential'].copy()
        elif force2e:
            self.df['Potential'] = self.df['Label Potential'].copy()
            print('3-electrode data detected. Ignoring working potential and using complete cell potential for everything. [NOT RECCOMMENDED]')
        else:
            print('3-electrode data detected. Using working potential for calculations and complete cell potential for labelling.')
        
        #plt.plot(self.df['Capacity'], self.df['Potential'])
        
        self.sigdf = self._find_sigcurves()
        #plt.plot(self.sigdf['Capacity'], self.sigdf['Potential'])
        self.sc_stepnums = self.sigdf['Prot_step'].unique()
        self.capacity = self.sigdf['Capacity'].max() - self.sigdf['Capacity'].min()
        self.spec_cap = self.capacity / self.mass
        if use_input_cap:
            self.capacity = self.input_cap
        else:
            print('Specific Capacity achieved by signature curves: {0:.2f} mAh/g'.format(self.spec_cap*1000))
            print('Using {:.8f} Ah to compute rates.'.format(self.capacity))

        self.caps, self.cumcaps, self.volts, self.fcaps, self.rates, self.eff_rates, self.currs, \
        self.ir, self.dqdv, self.resistdrop, self.icaps, self.avg_caps, self.ivolts, self.cvolts, \
        self.avg_volts, self.dvolts, self.vlabels = self._parse_sigcurves()
        
        self.nvolts = len(self.caps)
        
        if export_data:
            caprate_fname = self.dst / '{0} Parsed.xlsx'.format(self.cell_label)
            writer = pd.ExcelWriter(caprate_fname)
            for i in range(self.nvolts):
                if self.single_p is False:
                    caprate_df = pd.DataFrame(data = {'Specific Capacity': self.cumcaps[i],
                                                    'Fractional Capacity': self.fcaps[i],
                                                    'Effective C/n Rate': self.eff_rates[i],
                                                    'C/n Rate': self.rates[i]})
                else:
                    caprate_df = pd.DataFrame(data = {'Specific Capacity': self.caps[i],
                                                    'Voltage': self.volts[i],
                                                    'Fractional Capacity': self.fcaps[i],
                                                    'qi/I': self.eff_rates[i]})
                caprate_df.to_excel(writer, sheet_name = self.vlabels[i], index = False)
            print("Parsed data exporting to:\n" + str(caprate_fname))
            writer.save()
            writer.close()

    def _find_sigcurves(self):
        """
        Use control "step" to find sequence of charge/discharge - OCV 
        characteristic of signature curves.
        """
        newdf = self.df.drop_duplicates(subset = ['Step', 'Prot_step'])
        steps = newdf['Step'].values
        prosteps = newdf['Prot_step'].values
        ocv_inds = np.where(steps == 0)[0]
        
        if self.single_p is False:
            #print(ocv_inds)
            # Require a min of 3 OCV steps with the same step before and after
            # to qualify as a signature curve.
            #print(steps[2], steps[6])
            for i in range(len(ocv_inds)):
                #print(ocv_inds[i], steps[ocv_inds[i] - 1], steps[ocv_inds[i+2] + 1])
                if steps[ocv_inds[i] - 1] == steps[ocv_inds[i+2] + 1]:
                    first_sig_step = prosteps[ocv_inds[i] - 1]
                    break
            
            #last_sig_step = None
            for i in range(len(ocv_inds)):
                ind = -i - 1
                if steps[ocv_inds[ind] + 1] != steps[ocv_inds[ind] - 1]:
                    last_sig_step = prosteps[ocv_inds[ind] - 1]
                    break
                    
                elif steps[ocv_inds[ind] + 1] != steps[ocv_inds[ind] + 2]:
                    last_sig_step = prosteps[ocv_inds[ind] + 1]
                    break
                    
                #print(ocv_inds[-i-1], steps[ocv_inds[-i-1] - 1], steps[ocv_inds[-i-1] + 1])
                #if len(steps) > ocv_inds[-i-1] + 3:
                #    if steps[ocv_inds[-i-1]] != steps[ocv_inds[-i-1] + 2]:
                #        last_sig_step = prosteps[ocv_inds[-i-1] + 1]
                #        break
                #if (steps[ocv_inds[-i-1] - 1] != steps[ocv_inds[-i-1] + 1]):
                #    last_sig_step = prosteps[ocv_inds[-i-1] - 1]
                #    break
            #print(i)
            if i == len(ocv_inds) - 1:
                last_sig_step = prosteps[ocv_inds[-1] + 1]
        else:
            #single_pulse sigcurves selection
            for i in range(len(ocv_inds)):
                if i+1 == len(ocv_inds):
                    print("No adjacent OCV steps detected. Protocol is likely not single_pulse.")
                    break
                if ocv_inds[i] == ocv_inds[i+1] - 1:
                    first_sig_step = prosteps[ocv_inds[i] - 1]
                    break
            for i in range(len(ocv_inds)):
                if ocv_inds[-i] == ocv_inds[-i-1] + 1:
                    last_sig_step = prosteps[ocv_inds[-i]]
                    break    
        
        print('First signature curve step: {}'.format(first_sig_step))
        print('Last signature curve step: {}'.format(last_sig_step))
        
        sigdf = self.df.loc[(self.df['Prot_step'] >= first_sig_step) & (self.df['Prot_step'] <= last_sig_step)]
        
        return sigdf
    
    def plot_protocol(self, xlims = None, ylims = None, export_data = None, export_fig = True):
        
        if not(export_data is None): ('There is no data to export. Feel free to neglect this argument.')

        fig, axs = plt.subplots(nrows = 1, ncols = 2, sharey = True,
                                figsize = (6, 3), gridspec_kw = {'wspace':0.0})
        axs[0].plot(self.df['Time'], self.df['Label Potential'], 'k-')
        axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0].grid(which = 'minor', color = 'lightgrey')
        axs[0].set_xlabel('Time (h)')
        axs[0].set_ylabel('Voltage (V)')
        axs[1].set_xlabel('Specific Capacity \n (mAh/g)')
        axs[1].tick_params(which = 'both', axis = 'y', length = 0)
        axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1].grid(which = 'minor', color = 'lightgrey')
        #axs[0].tick_params(direction = 'in', top = True, right = True)
        
        # plot signature curves first if first
        if self.sc_stepnums[0] == 1:
            axs[1].plot(self.sigdf['Capacity']*1000/self.mass, self.sigdf['Label Potential'],
                color = 'red',
                label = 'Signature Curves')
        
        stepnums = self.df['Prot_step'].unique()
        #print(stepnums)
        fullsteps = np.setdiff1d(stepnums, self.sc_stepnums)
        #print(fullsteps)
        #print(self.sc_stepnums)
        # Need to set prop cycle
        colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(fullsteps)+1))
        c = 0
        for i in range(len(fullsteps)):

            stepdf = self.df.loc[self.df['Prot_step'] == fullsteps[i]]
            avgcurr = stepdf['Current'].mean()
            if avgcurr > 0.0:
                cyclabel = 'Charge'
            else:
                cyclabel = 'Discharge'
                
            if stepdf['Step'].values[0] == 0:
                label = 'OCV'
            else:
                avgcurr = np.absolute(avgcurr)
                minarg = np.argmin(np.absolute(RATES - self.capacity/avgcurr))
                rate = RATES[minarg]
                label = 'C/{0} {1}'.format(int(rate), cyclabel)
            
            axs[1].plot(stepdf['Capacity']*1000/self.mass, stepdf['Label Potential'],
                        color = colors[c],
                        label = label)
            
            c = c + 1
            
            # if the next step is the start of sigcurves, plot sigcurves
            if fullsteps[i] == self.sc_stepnums[0] - 1:
                #print('plotting sig curves...')
                axs[1].plot(self.sigdf['Capacity']*1000/self.mass, self.sigdf['Label Potential'],
                        color = 'red',
                        label = 'Signature Curves')
        
        plt.legend(bbox_to_anchor = (1.0, 0.5), loc = 'center left')
        
        if xlims is not None:
            axs[1].set_xlim(xlims[0], xlims[1])
        if ylims is not None:
            axs[0].set_ylim(ylims[0], ylims[1])
        
        if export_fig:
            figname = self.dst / '{} Protocol.jpg'.format(self.cell_label)
            print(figname)
            plt.savefig(self.dst / '{} Protocol.jpg'.format(self.cell_label), bbox_inches = 'tight')
            
        plt.show()
        plt.close()
    
    def plot_caps(self, export_data = None, export_fig = True):
        
        if not(export_data is None): ('There is no data to export. Feel free to neglect this argument.')
        
        fig, axs = plt.subplots(nrows = 2, ncols = 1, sharex = True, figsize = (3, 6), gridspec_kw = {'hspace':0.0})
        colors = plt.get_cmap('viridis')(np.linspace(0, 1, self.nvolts))
        
        for i in range(self.nvolts):
            axs[0].semilogx(self.eff_rates[i], self.cumcaps[i],
                            color = colors[i])
            axs[1].semilogx(self.eff_rates[i], self.fcaps[i],
                            color = colors[i], label = self.vlabels[i])
        
        if self.single_p:
            axs[1].set_xlabel('$q_{i}/I$ (h)')
        else:
            axs[1].set_xlabel('n$\mathregular{_{eff}}$ in C/n$\mathregular{_{eff}}$')
        axs[1].set_ylabel('$τ$')
        axs[0].set_ylabel('Specific Capacity \n (mAh/g)')
        axs[1].tick_params(axis = 'x', length = 0)
        axs[0].xaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[0].xaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0].grid(which = 'minor', color = 'lightgrey')
        axs[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1].grid(which = 'minor', color = 'lightgrey')
        plt.legend(bbox_to_anchor = (1.0, 1.0), loc = 'center left', ncol = 1 + self.nvolts//25)
        
        if export_fig:
            figname = self.dst / '{} Parsed.jpg'.format(self.cell_label)
            print(figname)
            plt.savefig(figname, bbox_inches = 'tight')

        plt.show()
        plt.close()

    def _parse_sigcurves(self):

        sigs = self.sigdf.loc[self.sigdf['Step'] != 0]
        #capacity = sigs['Capacity'].max() - sigs['Capacity'].min()
        #print('Specific Capacity: {} mAh'.format(capacity))
        Vstart = np.around(sigs['Label Potential'].values[0], decimals = 3)
        Vend = np.around(sigs['Label Potential'].values[-1], decimals = 3)
        print('Starting voltage: {:.3f} V'.format(Vstart))
        print('Ending voltage: {:.3f} V'.format(Vend))
        
        sigsteps = sigs['Prot_step'].unique()
        nsig = len(sigsteps)
        print('Found {} charge or discharge steps in signature curve sequences.'.format(nsig))
        caps = []
        cumcaps = []
        volts = []
        fcaps = []
        rates = []
        initcap = []
        cutcap = []
        initvolts = []
        cutvolts = []
        currs = []
        ir = []
        dqdv = []
        resistdrop = []
        eff_rates = []
        
        if self.single_p is False:
            for i in range(nsig):
                step = sigs.loc[sigs['Prot_step'] == sigsteps[i]]
                pulsecaps = step['Capacity'].values
                pulsevolts = step['Potential'].values
                currents = np.absolute(step['Current'].values)
                rate = self.capacity / np.average(currents)
                minarg = np.argmin(np.absolute(RATES - rate))
                
                # slice first and last current values if possible.
                # if less than 4(NVX) or 5(UHPC) data points, immediate voltage cutoff reached, omit step.
                if len(currents) > 3:
                    if pulsevolts[-2] == np.around(pulsevolts[-2], decimals = 2):
                        currents = currents[1:-1]
                        cvoltind = -2
                    elif len(currents) > 4:
                        if pulsevolts[-3] == np.around(pulsevolts[-3], decimals = 2):
                            currents = currents[1:-1]
                            cvoltind = -3
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
                    
                # determine dqdv based on the measurements before the voltage cutoff
                diffq = (pulsecaps[cvoltind-2] - pulsecaps[cvoltind-1]) / (pulsevolts[cvoltind-2] - pulsevolts[cvoltind-1])
                
                #if (np.amax(pulsecaps) - np.amin(pulsecaps))/self.mass < 5e-5:
                #    continue
            
                if caps == []:
                    caps.append([np.amax(pulsecaps) - np.amin(pulsecaps)])
                    volts.append([np.amax(pulsevolts) - np.amin(pulsevolts)])
                    rates.append([RATES[minarg]])
                    #initcutvolt = np.around(pulsevolts[0], decimals = 3)
                    initcap.append([pulsecaps[0]])
                    cutcap.append([pulsecaps[cvoltind]])
                    initvolts.append([pulsevolts[0]])
                    cutvolts.append([pulsevolts[cvoltind]])
                    currs.append([np.average(currents)])
                    ir.append([np.absolute(pulsevolts[0] - pulsevolts[1])])
                    dqdv.append([diffq])
                    resistdrop.append([ir[-1][-1]/currs[-1][-1]])
                else:
                    #if np.amax(currents) < currs[-1][-1]:
                    if pulsevolts[cvoltind] == cutvolts[-1][-1]:
                        caps[-1].append(np.amax(pulsecaps) - np.amin(pulsecaps))
                        volts[-1].append(np.amax(pulsevolts) - np.amin(pulsevolts))
                        rates[-1].append(RATES[minarg])
                        cutcap[-1].append(pulsecaps[cvoltind])
                        cutvolts[-1].append(pulsevolts[cvoltind])
                        currs[-1].append(np.average(currents))
                        ir[-1].append(np.absolute(pulsevolts[0] - pulsevolts[1]))
                        dqdv[-1].append(diffq)
                        resistdrop[-1].append(ir[-1][-1]/currs[-1][-1])
                    else:
                        if np.absolute(pulsevolts[-2] - cutvolts[-1][-1]) < 0.001:
                            continue
                        #print(np.average(currents), pulsevolts[-2])
                        caps.append([np.amax(pulsecaps) - np.amin(pulsecaps)])
                        volts.append([np.amax(pulsevolts) - np.amin(pulsevolts)])
                        rates.append([RATES[minarg]])
                        initcap.append([pulsecaps[0]])
                        cutcap.append([pulsecaps[cvoltind]])
                        initvolts.append([pulsevolts[0]])
                        cutvolts.append([pulsevolts[cvoltind]])
                        currs.append([np.average(currents)])
                        ir.append([np.absolute(pulsevolts[0] - pulsevolts[1])])
                        dqdv.append([diffq])
                        resistdrop.append([ir[-1][-1]/currs[-1][-1]])
            
            nvolts = len(caps)
            for i in range(nvolts):
                fcaps.append(np.cumsum(caps[i]) / np.sum(caps[i]))
                cumcaps.append(np.cumsum(caps[i]))
                
                eff_rates.append(cumcaps[i][-1]/currs[i])
                
                # Remove data where capacity is too small due to IR
                # i.e., voltage cutoff was reached immediately.
                inds = np.where(fcaps[i] < self.fcap_min)[0]
                if len(inds) > 0:
                    caps[i] = np.delete(caps[i], inds)
                    volts[i] = np.delete(volts[i], inds)
                    cumcaps[i] = np.delete(cumcaps[i], inds)
                    fcaps[i] = np.delete(fcaps[i], inds)
                    eff_rates[i] = np.delete(eff_rates[i], inds)
                    rates[i] = np.delete(rates[i], inds)
                    cutcap[i] = np.delete(cutcap[i], inds)
                    cutvolts[i] = np.delete(cutvolts[i], inds)
                    currs[i] = np.delete(currs[i], inds)
                    ir[i] = np.delete(ir[i], inds)
                    dqdv[i] = np.delete(dqdv[i], inds)
                    resistdrop[i] = np.delete(resistdrop[i], inds)
                    print("Current removed due to being below fcap min")
        
            if self.capacitance_corr == True:
                print("Capacitance correction cannot be applied to multi-pulse AMID data. Data is being analyzed without capacitance correction.")
        else:
            # idcaps is the idealized capacity for a given voltage based upon dqdv
            # cumcurrs is the cumulative averge current for determining where C should be calculated
            idcaps = []
            cumcurrs = []
            voltsAct = []
            time = []
            for i in range(nsig):
                step = sigs.loc[sigs['Prot_step'] == sigsteps[i]]
                pulsecaps = step['Capacity'].values
                pulsevolts = step['Potential'].values
                lpulsevolts = step['Label Potential'].values
                currents = np.absolute(step['Current'].values)
                runtime = step['Time'].values
                
                # Collect succeeding OCV steps (1 OCV or 2 OCV) to calculate dqdv
                ocvstep = self.sigdf.loc[self.sigdf['Prot_step'] == sigsteps[i] + 2]
                if ocvstep['Step'].values[0] != 0:
                    ocvstep = self.sigdf.loc[self.sigdf['Prot_step'] == sigsteps[i] + 1]
                ocvpulsecaps = ocvstep['Capacity'].values
                ocvvolts = ocvstep['Potential'].values
                ocvlvolts = ocvstep['Label Potential'].values
                
                ir.append([np.absolute(pulsevolts[1] - pulsevolts[0])])
                                
                initcap.append([pulsecaps[0]])
                cutcap.append([ocvpulsecaps[-1]]) 
                initvolts.append([lpulsevolts[0]])
                cutvolts.append([ocvlvolts[-1]]) 
                
                dqdv.append([(pulsecaps[0] - ocvpulsecaps[-1])/(pulsevolts[0] - ocvvolts[-1])])
                
                time.append(np.absolute(runtime[1:] - runtime[0]))
                caps.append(np.absolute(pulsecaps[1:] - pulsecaps[0]))
                cumcaps.append(np.absolute(pulsecaps[1:] - pulsecaps[0]))
                volts.append(np.absolute(pulsevolts[1:] - pulsevolts[0]))
                idcaps.append(np.absolute(dqdv[-1][0]*(pulsevolts[1:] - pulsevolts[0])))
                currs.append(currents[1:])
                voltsAct.append(pulsevolts[1:])
                
                cumcurrs.append([])
                for j in range(len(currents[1:])):
                    cumcurrs[-1].append(np.average(currents[1:j+2]))
                    minarg = np.argmin(np.absolute(RATES - self.capacity / cumcurrs[-1][j]))
                    if j == 0:
                        rates.append([RATES[minarg]])
                    else:
                        rates[-1].append(RATES[minarg])
                
                resistdrop.append([ir[-1][-1]/currs[-1][0]])
 
            nvolts = len(caps)
            
            if self.capacitance_corr == True:
                #DL capacitance is calculated from the first 5 consistent current datapoints in the lowest V pulse.
                lowVind = cutvolts.index(min(cutvolts))

                A = np.ones((8, 2))
                y = np.zeros(8)
                n = 0
                for i in range(len(currs[lowVind])):
                    if n == 8:
                        break
                    if abs(currs[lowVind][i]/currs[lowVind][i+1] - 1) < 0.01:
                        A[n][0] = voltsAct[lowVind][i+1]
                        y[n] = caps[lowVind][i+1]
                        n = n + 1
                    else:
                        A = np.ones((8, 2))
                        y = np.zeros(8)
                        n = 0
                
                capacitance = abs(np.linalg.lstsq(A, y)[0][0])
                print('Double layer capacitance found at lowest V pulse: {:.2f} nF'.format(1.0e9*capacitance))
                
                rohm = np.power(10, stats.mode(np.round(np.log10(resistdrop), 2))[0])[0][0]
                print('Logarithmic mode of ohmic resistance over all pulses: {:.2f} Ω'.format(rohm))
                for i in range(nvolts):
                    dlcaps = []
                    for j in range(len(voltsAct[i])):
                        if voltsAct[i][0] > voltsAct[i][-1]:
                            if voltsAct[i][0] + ir[i][0] - currs[i][j]*rohm > voltsAct[i][j]:
                                dlcaps.append(capacitance*((voltsAct[i][0] + ir[i][0] - currs[i][j]*rohm) - voltsAct[i][j]))
                            else:
                                dlcaps.append(0)
                        else:
                            if voltsAct[i][0] - ir[i][0] + currs[i][j]*rohm < voltsAct[i][j]:
                                dlcaps.append(capacitance*(voltsAct[i][j] - (voltsAct[i][0] - ir[i][0] + currs[i][j]*rohm)))
                            else:
                                dlcaps.append(0)

                    caps[i] = caps[i] - dlcaps
                    # if caps is calculated as negative, this datapoint is effectively thrown out (caps set to 0)
                    for j in range(len(caps[i])):
                        if caps[i][j] < 0:
                            caps[i][j] = 0
                    
                    idcaps[i] = idcaps[i] - dlcaps #- dqdv[i][0]*currs[i]*rohm 
                    # if idcaps is calculated as negative or zero, this datapoint is effectively thrown out (caps set to nan)
                    for j in range(len(caps[i])):
                        if idcaps[i][j] <= 0:
                            idcaps[i][j] = float('NaN')
                    
                    #cumcurrs[i] = cumcurrs[i] - dlcaps/time[i] # disabled as it may amplify error if near 0
                    # if cumulative current is calculated as negative or zero, this datapoint is effectively thrown out (caps set to nan)
                    for j in range(len(caps[i])):
                        if cumcurrs[i][j] <= 0:
                            cumcurrs[i][j] = float('NaN')
            
            for i in range(nvolts):
                fcaps.append(caps[i]/idcaps[i])
                eff_rates.append(idcaps[i]/cumcurrs[i])
                
                # outlier repair: if fcap or eff_rates is NaN, make it equal to the succeeding point (or previous if last point).
                for j in range(len(caps[i])):
                    if fcaps[i][-j - 1] != fcaps[i][-j - 1] or eff_rates[i][-j - 1] != eff_rates[i][-j - 1]:
                        fcaps[i][-j - 1] = fcaps[i][-j]
                        eff_rates[i][-j - 1] = eff_rates[i][-j]
                        
            if self.fcap_min != 0.0:
                print("Fractional capacity exclusion cannot be applied to single-pulse AMIDR data. Data is being analyzed without fractional capacity exclusion.")
        
        print('Found {} signature curves.'.format(nvolts))
        
        ivolts = np.zeros(nvolts)
        cvolts = np.zeros(nvolts)
        icaps = np.zeros(nvolts)
        ccaps = np.zeros(nvolts)
        for i in range(nvolts):
            ivolts[i] = np.average(initvolts[i])
            cvolts[i] = np.average(cutvolts[i])
            icaps[i] = np.average(initcap[i])
            ccaps[i] = cutcap[i][-1]
        
        with np.printoptions(precision = 3):
            avg_caps = (icaps + ccaps)/2
            # Get midpoint voltage for each range.
            #avg_volts[0] = (initcutvolt + cvolts[0])/2
            #avg_volts[1:] = (cvolts[:-1] + cvolts[1:])/2
            avg_volts = (ivolts + cvolts)/2                
            dvolts = np.zeros(nvolts)
            #dvolts[0] = np.absolute(initcutvolt - cvolts[0])
            #dvolts[1:] = np.absolute(cvolts[:-1] - cvolts[1:])
            dvolts = np.absolute(ivolts - cvolts)
            if self.single_p is False:
                print('Midpoint capacities (mAh/g): {}'.format((avg_caps*1000/self.mass)))
                print('Cutoff voltages: {}'.format(cvolts))
                print('Midpoint voltages: {}'.format(avg_volts))
                with np.printoptions(precision = 4):
                    print('Voltage intervals widths: {}'.format(dvolts))
            # Make voltage interval labels for legend.
            #vlabels = ['{0:.3f} V - {1:.3f} V'.format(initcutvolt, cvolts[0])]
            #vlabels = vlabels + ['{0:.3f} V - {1:.3f} V'.format(cvolts[i], cvolts[i+1]) for i in range(nvolts-1)]
            vlabels = ['{0:.3f} V - {1:.3f} V'.format(ivolts[i], cvolts[i]) for i in range(nvolts)]
            print('Voltage interval labels: {}'.format(vlabels))
            print('Found {} voltage intervals.'.format(nvolts))
        
        iadj = 0
        for i in range(nvolts):    
            if len([j for j in fcaps[i-iadj] if j>0.001]) < 4:
                print("{} removed due to not having 4 or more datapoints with relative change in capacity ($τ$) above 0.001".format(vlabels[i-iadj]))
                caps.pop(i-iadj)
                cumcaps.pop(i-iadj)
                volts.pop(i-iadj)
                fcaps.pop(i-iadj)
                eff_rates.pop(i-iadj)
                rates.pop(i-iadj)
                currs.pop(i-iadj)
                ir.pop(i-iadj)
                dqdv.pop(i-iadj)
                resistdrop.pop(i-iadj)
                icaps = np.delete(icaps, i-iadj)
                avg_caps = np.delete(avg_caps, i-iadj)
                ivolts = np.delete(ivolts, i-iadj)
                avg_volts = np.delete(avg_volts, i-iadj)
                dvolts = np.delete(dvolts, i-iadj)
                vlabels = np.delete(vlabels, i-iadj)
                iadj = iadj + 1
        nvolts = len(caps)
        
        speccaps = []
        speccumcaps = []
        for i in range(nvolts):
            speccaps.append(1000*np.array(caps[i])/self.mass)
            speccumcaps.append(1000*np.array(cumcaps[i])/self.mass)

        return speccaps, speccumcaps, volts, fcaps, rates, eff_rates, currs, ir, dqdv, resistdrop, icaps, avg_caps, ivolts, cvolts, avg_volts, dvolts, vlabels 
       
    def fit_atlung(self, r, R_corr, ionsat_inputs = [], micR_input = 4.9, ftol = 5e-14, D_bounds = [1e-17, 1e-8], D_guess = 1.0e-11, 
                   fcapadj_bounds = [1.0, 1.5], fcapadj_guess = 1.0, P_bounds = [1e-6, 1e1], P_guess = 1.0e-2, remove_out_of_bounds = True,
                   shape = 'sphere', nalpha = 4000, nQ = 4000, export_data = True, export_fig = True, fitlabel = None):

        self.r = r
        self.R_corr = R_corr
        
        if fitlabel is None:
            cell_label = self.cell_label
        else:
            cell_label = self.cell_label + '-' + fitlabel
        
        if shape not in SHAPES:
            print('The specified shape {0} is not supported.'.format(shape))
            print('Supported shapes are: {0}. Defaulting to sphere.'.format(SHAPES))
            
        # Get geometric constants according to particle shape.
        if shape == 'sphere':
            self.alphas = []
            for i in np.arange(4, 4*nalpha):
                g = lambda a: a/np.tan(a) - 1
                sol = fsolve(g, i)
                self.alphas.append(sol)
            self.alphas = np.unique(np.around(self.alphas, 8))**2
            self.alphas = self.alphas[:nalpha]
            A, B = 3, 5

        elif shape == 'plane':
            self.alphas = (np.arange(1, nalpha+1)*np.pi)**2
            A, B = 1, 3
                
        # Solve for tau vs Q
        if self.R_corr is False:
            print("Optimum Parameters: {}".format("Log(Dc) fCapAdj"))
            Q_arr = np.logspace(-3, 2, nQ)
            tau_sol = np.zeros(nQ)
            tau_guess = 0.5
            for i in range(nQ):
                Q = Q_arr[i]
                func = lambda tau: tau - 1 + (1/(A*Q))*(1/B - 2*(np.sum(np.exp(-self.alphas*tau*Q)/self.alphas)))
                tau_sol[i] = fsolve(func, tau_guess, factor = 1.)
        elif self.single_p is False:
            print("Optimum Parameters: {}".format("Log(Dc) fCapAdj Log(P) Log(P/Dc)"))
        else:
            print("Optimum Parameters: {}".format("Log(Dc) Log(P) Log(P/Dc)"))
                
        dconst = np.zeros(self.nvolts, dtype = float)
        dtconst = np.zeros(self.nvolts, dtype = float)
        pconst = np.zeros(self.nvolts, dtype = float)
        resist = np.zeros(self.nvolts, dtype = float)
        dqdv = np.zeros(self.nvolts, dtype = float)
        sigma = np.zeros(self.nvolts, dtype = float)
        fit_err = np.zeros(self.nvolts, dtype = float)
        cap_max = np.zeros(self.nvolts, dtype = float)
        cap_min = np.zeros(self.nvolts, dtype = float)
        cap_span = np.zeros(self.nvolts, dtype = float)

        for j in range(self.nvolts):
            z = np.ones(len(self.fcaps[j]))
            #fcap = np.array(self.fcaps[j])
            fcap = np.array(self.fcaps[j])
            #self._max_cap = self.cumcaps[j][-1]
            #print('Max cap: {} mAh/g'.format(self._max_cap))
            rates = np.array(self.eff_rates[j])
            #I = np.array(self.currs[j])*1000
            #print("Currents: {} mA".format(I))
            #self._dqdv = np.average(self.dqdV[j][-1])*1000/self.mass
            
            if self.single_p is False:
                # selects the dqdv of C/40 discharge/charge or nearest to C/40
                act_rates = self.capacity / np.array(self.currs[j])
                minarg = np.argmin(np.absolute(40 - act_rates))
                dqdv[j] = self.dqdv[j][minarg]
            else:
                dqdv[j] = self.dqdv[j][0]
                # constrains fcapadj to 1
                fcapadj_bounds = [1.0, 1.0000001]

            #print("dq/dV: {} Ah/V".format(dqdv[j]))
            
            if self.R_corr is False:
                C = np.sum(self.ir[j])
                weights = (C - self.ir[j]) / np.sum(C - self.ir[j])
                bounds = ([np.log10(D_bounds[0]), fcapadj_bounds[0]],
                          [np.log10(D_bounds[1]), fcapadj_bounds[1]])
                p0 = [np.log10(D_guess), fcapadj_guess]    
            else:
                bounds = ([np.log10(D_bounds[0]), fcapadj_bounds[0], np.log10(P_bounds[0])],
                          [np.log10(D_bounds[1]), fcapadj_bounds[1], np.log10(P_bounds[1])])
                p0 = [np.log10(D_guess), fcapadj_guess, np.log10(P_guess)]
            
            if shape == 'sphere':
                if self.R_corr is False:
                    popt, pcov = curve_fit(self._spheres, (fcap, rates), z, p0 = p0,
                               bounds = bounds, sigma = weights,
                               method = 'trf', max_nfev = 5000, x_scale = [1.0, 1.0],
                               ftol = ftol, xtol = None, gtol = None, loss = 'soft_l1', f_scale = 1.0)
                    with np.printoptions(precision = 4):
                        print("{}: {}".format(self.vlabels[j], popt))
                else:
                    p0opt = [p0[2] - p0[0], p0[1], p0[2]]
                    boundsopt = [[bounds[0][2] - bounds[1][0], bounds[0][1], bounds[0][2]], 
                                 [bounds[1][2] - bounds[0][0], bounds[1][1], bounds[1][2]]]
                    
                    popt, pcov = curve_fit(self._spheres_R_corr, (fcap, rates), z, p0 = p0opt,
                               bounds = boundsopt,
                               method = 'trf', max_nfev = 5000, x_scale = [1.0, 1.0, 1.0],
                               ftol = ftol, xtol = None, gtol = None, loss = 'soft_l1', f_scale = 1.0)
                    popt = np.array([popt[2] - popt[0], popt[1], popt[2], popt[0]])
                    with np.printoptions(precision = 3):
                        if self.single_p is False:
                            print("{}: {}".format(self.vlabels[j], popt))
                        else:
                            if remove_out_of_bounds and (round(popt[2], 5) == round(bounds[0][2], 5) or round(popt[2], 5) == round(bounds[1][2], 5)):
                                print("{}: {}".format(self.vlabels[j], 'No fit within P bounds'))
                                popt = np.array([float('NaN'), float('NaN'), float('NaN'), float('NaN')])
                            else:
                                print("{}: {}".format(self.vlabels[j], np.array([popt[0], popt[2], popt[3]])))
                    pconst[j] = 10**popt[2]
                    Q_arr = np.logspace(-6, 2, nQ)
                    tau_sol = np.zeros(nQ)
                    tau_guess = 0.5
                    for i in range(nQ):
                        Q = Q_arr[i]
                        func = lambda tau: tau - 1 + (1/(A*Q))*(1/B - 2*(np.sum(np.exp(-self.alphas*tau*Q)/self.alphas))) + 10**popt[2]/Q if 10**popt[2]<Q else tau

                        tau_sol[i] = fsolve(func, tau_guess, factor = 1.)
                        if tau_sol[i] < 0:
                            tau_sol[i] = 0
                    
            if shape == 'plane':
                popt, pcov = curve_fit(self._planes, (fcap, rates), z, p0 = p0,
                           bounds = bounds, sigma = weights,
                           method = 'trf', max_nfev = 5000, x_scale = [1e-11, 1.0],
                           ftol = ftol, xtol = None, gtol = None, loss = 'soft_l1', f_scale = 1.0)
            
            #(Q_arr, tau_sol, '-k', label = 'Atlung - {}'.format(shape))
            
            sigma[j] = np.sqrt(np.diag(pcov))[0]
            dconst[j] = 10**popt[0]
            Qfit = 3600*rates*dconst[j]/r**2
            tau_fit = fcap/popt[1]
            
            cap_max[j] = tau_fit[-1]
            cap_min[j] = tau_fit[0]
            cap_span[j] = tau_fit[-1] - tau_fit[0]
            
            # Get difference between fitted values and theoretical Atlung curve to get fit_err.
            error = np.zeros(len(Qfit), dtype = float)
            for k in range(len(Qfit)):
                dQ = np.absolute(Q_arr - Qfit[k])
                minarg = np.argmin(dQ)
                error[k] = np.absolute(tau_fit[k] - tau_sol[minarg])
            if R_corr is False:
                fit_err[j] = np.sum(weights*error)
            else:
                fit_err[j] = np.sqrt(np.average((error/cap_max[j])**2))
            
            plt.figure(figsize = (3, 3))
            plt.semilogx(Qfit, tau_fit, 'or', markersize = 2, label = 'Experimental')
            if max(tau_fit) < 0.01:
                plt.ylim(0, 0.01)
            elif max(tau_fit) < 0.1:
                plt.ylim(0, 0.1)
            else:
                plt.ylim(0, 1)
            plt.semilogx(Q_arr, tau_sol, '-k', label = 'Model')
            plt.xlabel('$Q$')
            plt.ylabel('$τ$')
            plt.legend(loc = 'upper left', frameon = True)
            if Qfit[0] < 1.0e-4 or Qfit[1] < 1.0e-3:
                plt.xlim(1.0e-6, 1.0e0)
                plt.xticks(10.**np.arange(-6, 1))
            else:
                plt.xlim(1.0e-4, 1.0e2)
                plt.xticks(10.**np.arange(-4, 3))
            plt.gca().xaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
            plt.gca().xaxis.set_major_locator(ticker.LogLocator(numticks = 10))
            plt.gca().yaxis.set_minor_locator(ticker.AutoMinorLocator())
            plt.gca().grid(which = 'minor', color = 'lightgrey')
            if export_fig:
                figname = self.dst / '{0} {1:.3f} V ({2}).jpg'.format(cell_label, self.avg_volts[j], shape)
                if not(np.isnan(popt[0])):
                    plt.savefig(figname, bbox_inches = 'tight')
            
            plt.close()
            
        if export_fig:
            print("Figures not shown saved to:\n" + str(self.dst))
    
        if R_corr is False:
            DV_df = pd.DataFrame(data = {'Voltage': self.avg_volts, 'D': dconst})
        else:
            # get resist from pconst
            resist = pconst*self.r**2/(3600*dconst*dqdv)
            #print("Resist: {} V/A".format(resist))
            
            # get resist from ir drop
            rdrop = []
            for i in range(len(resist)):
                if self.single_p is False:
                    rdrop.append(np.average(self.resistdrop[i]))
            #        resistdrop.append([])
            #        for j in range(len(self.ir[i])):
            #            resistdrop[i][j].append(self.ir[i][j]/self.currs[i][j])
                else:
                    rdrop.append(self.resistdrop[i][0])
            
            DV_df = pd.DataFrame(data = {'Voltage (V)': self.avg_volts, 'Initial Voltage (V)': self.ivolts, 'Dc (cm^2/s)': dconst, 'P' : pconst, 'dq/dV (mAh/gV)': [i*1000/self.mass for i in dqdv],
                                       'Rfit (Ohm)' : resist, 'micR (Ohmcm^2)' : A*resist*self.mass/(r*micR_input), 'Rdrop (Ohm)' : rdrop, 'Cap Span' : cap_span, 'Fit Error' : fit_err})
            
        # Calculates Free-path Tracer D from inputs
        if ionsat_inputs:
            temp = ionsat_inputs[0]
            theorcap = ionsat_inputs[1]/1000*self.mass
            
            if ionsat_inputs[2] > max(self.ivolts) or ionsat_inputs[2] < min(self.ivolts) or ionsat_inputs[4] > max(self.ivolts) or ionsat_inputs[4] < min(self.ivolts):
                print('Voltage inputs for calculating Free-path Tracer D are out of range. Tracer D will not be calculated.')
            else:
                volt1 = ionsat_inputs[2]
                newcap1 = ionsat_inputs[3]/1000*self.mass
                volt2 = ionsat_inputs[4]
                newcap2 = ionsat_inputs[5]/1000*self.mass
                volt1a = float('NaN')
                volt1b = float('NaN')
                volt2a = float('NaN')
                volt2b = float('NaN')
                
                for i in range(len(self.ivolts) - 1):
                    if volt1 > self.ivolts[i] and not(volt1a > self.ivolts[i]):
                        volt1a = self.ivolts[i]
                        cap1a = self.icaps[i]
                    if volt1 < self.ivolts[i] and not(volt1b < self.ivolts[i]):
                        volt1b = self.ivolts[i]
                        cap1b = self.icaps[i]
                    if volt2 > self.ivolts[i] and not(volt2a > self.ivolts[i]):
                        volt2a = self.ivolts[i]
                        cap2a = self.icaps[i]
                    if volt2 < self.ivolts[i] and not(volt2b < self.ivolts[i]):
                        volt2b = self.ivolts[i]
                        cap2b = self.icaps[i]
                
                cap1 = cap1a + (volt1 - volt1a)/(volt1b - volt1a)*(cap1b - cap1a)
                cap2 = cap2a + (volt2 - volt2a)/(volt2b - volt2a)*(cap2b - cap2a)
                m = (newcap2 - newcap1)/(cap2 - cap1)
                b = newcap1 - m*cap1
                
                kB = 1.380649E-23
                e = 1.602176634E-19
                pulsecaps = m*self.avg_caps + b
                ipulsecaps = m*self.icaps + b
                socs = 1. - pulsecaps/theorcap
                isocs = 1. - ipulsecaps/theorcap
                dtconst = dconst*kB*temp*dqdv/(e*pulsecaps*(1. - pulsecaps/theorcap))
                
                if R_corr:
                    DV_df['Dt* (cm^2/s)'] = dtconst
                    DV_df['SOC'] = socs
                    DV_df['Initial SOC'] = isocs
                    DV_df = DV_df[['Voltage (V)', 'Initial Voltage (V)', 'SOC', 'Initial SOC', 'Dc (cm^2/s)', 'Dt* (cm^2/s)', 'P', 'dq/dV (mAh/gV)', 'Rfit (Ohm)', 'micR (Ohmcm^2)', 'Rdrop (Ohm)', 'Cap Span', 'Fit Error']]

        if export_data:
            df_filename = self.dst / '{0} Fitted ({1}).xlsx'.format(cell_label, shape)
            print("Fitted data exporting to:\n" + str(df_filename))
            DV_df.to_excel(df_filename, index = False)
        
        with np.printoptions(precision = 3):
            print('Fitted Dc: {}'.format(dconst))
            if self.single_p is False:
                print('Standard deviations from fit: {}'.format(sigma))
                print('Atlung fit error: {}'.format(fit_err))
        
        return self.avg_volts, self.ivolts, dconst, dtconst, fit_err, cap_span, cap_max, cap_min, self.caps, self.ir, self.dvolts, pconst, dqdv, resist, self.resistdrop, self.single_p, self.R_corr, cell_label, self.mass, self.dst
        
    def make_summary_graph(self, fit_data, export_data = None, export_fig = True):
        
        if not(export_data is None): ('There is no figure to export. Feel free to neglect this argument.')
        
        voltage = fit_data[0]
        nvolts = len(voltage)
        ivoltage = fit_data[1]
        dconst = fit_data[2]
        dtconst = fit_data[3]
        fit_err = fit_data[4]
        cap_span = fit_data[5]
        cap_max = fit_data[6]
        cap_min = fit_data[7]
        caps = fit_data[8]
        dV_ir = fit_data[9]
        dvolts = fit_data[10]
        pconst = fit_data[11]
        dqdv = fit_data[12]
        resist = fit_data[13]
        resistdrop = fit_data[14]
        single_p = fit_data[15]
        R_corr = fit_data[16]
        cell_label = fit_data[17]
        mass = fit_data[18]
        dst = fit_data[19]
        
        if single_p is False:
            if R_corr is False:
                fig, axs = plt.subplots(ncols = 1, nrows = 5, figsize = (6, 12), sharex = True,
                gridspec_kw = {'height_ratios': [1, 0.5, 0.5, 0.5, 0.5], 'hspace': 0.0})
                
                axs[0].semilogy(voltage, dconst, 'kx-', label = '$D_{c}$')
                axs[0].semilogy(voltage, dtconst, 'k.:', label = '$D_{t}^{*}$')
                axs[0].set_xlabel('Voltage (V)')
                axs[0].set_ylabel('$D$ (cm$\mathregular{^{2}}$/s)')
                axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
                axs[0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
                axs[0].grid(which = 'minor', color = 'lightgrey')
                if sum(dtconst) != 0:
                    axs[0].legend(frameon = True)
                
                axs[1].semilogy(voltage, fit_err, 'kx-')
                axs[1].set_xlabel('Voltage (V)')
                axs[1].set_ylabel('Fit Error')
                axs[1].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
                axs[1].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
                axs[1].grid(which = 'minor', color = 'lightgrey')
                
                axs[2].set_ylim(0, 1.0)
                axs[2].get_yaxis().set_ticks([0.25, 0.5, 0.75])
                axs[2].set_xlabel('Voltage (V)')
                axs[2].set_ylabel('$τ$ Span')
                axs[2].yaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[2].grid(which = 'minor', color = 'lightgrey')
                axs[2].set_axisbelow(True)
                for j in range(nvolts):
                    axs[2].fill(np.array([voltage[j]-0.01, voltage[j]-0.01, voltage[j]+0.01, voltage[j]+0.01]),
                                np.array([cap_min[j], cap_max[j], cap_max[j], cap_min[j]]),
                                color = 'grey', edgecolor = 'k', linestyle = '-')
                
                for j in range(nvolts):
                    cap_in_step = caps[j]
                    ir = dV_ir[j]
                    axs[3].bar(voltage[j], np.sum(cap_in_step), width = dvolts[j], color = 'grey', alpha = 0.5, edgecolor = 'k', linestyle = '-')
                    for i in range(len(ir)):
                        width = dvolts[j]/len(cap_in_step)
                        center = voltage[j] - dvolts[j]/2 + (i+1/2)*width
                        axs[3].bar(center, cap_in_step[i], width = width, color = 'grey', alpha = 0.5, edgecolor = 'k', linestyle = '-')
                        axs[4].bar(voltage[j], ir[i]/dvolts[j], width = 0.04, color = 'k', alpha = 0.3, edgecolor = 'k', linestyle = '-')
                axs[3].get_xaxis().set_ticks(voltage)
                axs[3].tick_params(axis = 'x', which = 'minor', top = False, bottom = False)
                axs[3].set_xticklabels(['{:.3f}'.format(v) for v in voltage], rotation = 45)
                axs[3].set_xlabel('Voltage (V)')
                axs[3].set_ylabel('Specific\nCapacity\n(mAh/g)')
                axs[3].yaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[3].grid(which = 'minor', color = 'lightgrey')
                axs[3].set_axisbelow(True)
                
                axs[4].get_xaxis().set_ticks(voltage)
                axs[4].tick_params(axis = 'x', which = 'minor', top = False, bottom = False)
                axs[4].set_xticklabels(['{:.3f}'.format(v) for v in voltage], rotation = 45)
                axs[4].set_xlabel('Voltage (V)')
                axs[4].set_ylabel('Fractional\nIR drop')
                axs[4].set_axisbelow(True)
                
            else:
                fig, axs = plt.subplots(ncols = 1, nrows = 6, figsize = (6, 12), sharex = True,
                gridspec_kw = {'height_ratios': [1, 1, 0.5, 0.5, 0.5, 0.5], 'hspace': 0.0})
                
                axs[0].semilogy(voltage, dconst, 'kx-', label = '$D_{c}$')
                axs[0].semilogy(voltage, dtconst, 'k.:', label = '$D_{t}^{*}$')
                axs[0].set_xlabel('Voltage (V)')
                axs[0].set_ylabel('$D$ (cm$\mathregular{^{2}}$/s)')
                axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
                axs[0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
                axs[0].grid(which = 'minor', color = 'lightgrey')
                if sum(dtconst) != 0:
                    axs[0].legend(frameon = True)
                
                axs[1].semilogy(voltage, resist, 'kx-', label = 'Fit R')
                rdavg = []
                rddev = []
                for i in range(nvolts):
                    rdavg.append(np.average(resistdrop[i]))
                    rddev.append(np.std(resistdrop[i]))
                axs[1].errorbar(voltage, rdavg, rddev, fmt = 'k.:', capsize = 3.0, label = 'V Drop R')
                axs[1].set_xlabel('Voltage (V)')
                axs[1].set_ylabel('$R$ (Ω)')
                axs[1].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
                axs[1].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
                axs[1].grid(which = 'minor', color = 'lightgrey')
                axs[1].legend(frameon = True)
                
                axs[2].semilogy(voltage, pconst, 'kx-')
                axs[2].semilogy(voltage, 1.0/15*np.ones(len(voltage)), 'k--', linewidth = 1, label = 'Equivalent Limitation')
                axs[2].set_ylim(2.0e-3, 0.8)
                axs[2].set_xlabel('Voltage (V)')
                axs[2].set_ylabel('$P$')
                axs[2].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
                axs[2].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
                axs[2].grid(which = 'minor', color = 'lightgrey')
                axs[2].legend(frameon = True)
                
                axs[3].semilogy(voltage, fit_err, 'kx-')
                axs[3].set_xlabel('Voltage (V)')
                axs[3].set_ylabel('Fit Error')
                axs[3].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
                axs[3].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
                axs[3].grid(which = 'minor', color = 'lightgrey')
                
                axs[4].set_ylim(0, 1.0)
                axs[4].get_yaxis().set_ticks([0.25, 0.5, 0.75])
                axs[4].set_xlabel('Voltage (V)')
                axs[4].set_ylabel('$τ$ Span')
                axs[4].yaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[4].grid(which = 'minor', color = 'lightgrey')
                axs[4].set_axisbelow(True)
                for j in range(nvolts):
                    axs[4].fill(np.array([voltage[j]-0.01, voltage[j]-0.01, voltage[j]+0.01, voltage[j]+0.01]),
                                np.array([cap_min[j], cap_max[j], cap_max[j], cap_min[j]]),
                                color = 'grey', edgecolor = 'k', linestyle = '-')
                
                for j in range(nvolts):
                    cap_in_step = caps[j]
                    ir = dV_ir[j]
                    axs[5].bar(voltage[j], np.sum(cap_in_step), width = dvolts[j], color = 'grey', alpha = 0.5, edgecolor = 'k', linestyle = '-')
                    for i in range(len(ir)):
                        width = dvolts[j]/len(cap_in_step)
                        center = voltage[j] - dvolts[j]/2 + (i+1/2)*width
                        axs[5].bar(center, cap_in_step[i], width = width, color = 'grey', alpha = 0.5, edgecolor = 'k', linestyle = '-')
                axs[5].get_xaxis().set_ticks(voltage)
                axs[5].tick_params(axis = 'x', which = 'minor', top = False, bottom = False)
                axs[5].set_xticklabels(['{:.3f}'.format(v) for v in voltage], rotation = 45)
                axs[5].set_xlabel('Voltage (V)')
                axs[5].set_ylabel('Specific\nCapacity\n(mAh/g)')
                axs[5].yaxis.set_minor_locator(ticker.AutoMinorLocator())
                axs[5].grid(which = 'minor', color = 'lightgrey')
                axs[5].set_axisbelow(True)

        else:
            fig, axs = plt.subplots(ncols = 1, nrows = 6, figsize = (6, 12), sharex = True,
            gridspec_kw = {'height_ratios': [1, 1, 0.5, 0.5, 0.5, 0.5], 'hspace': 0.0})
            
            axs[0].semilogy(voltage, dconst, 'kx-', label = '$D_{c}$')
            axs[0].semilogy(voltage, dtconst, 'k.:', label = '$D_{t}^{*}$')
            axs[0].set_xlabel('Voltage (V)')
            axs[0].set_ylabel('$D$ (cm$\mathregular{^{2}}$/s)')
            axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs[0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
            axs[0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
            axs[0].grid(which = 'minor', color = 'lightgrey')
            if sum(dtconst) != 0:
                axs[0].legend(frameon = True)
            
            axs[1].semilogy(ivoltage, resist, 'kx-', label = 'Fit R')
            axs[1].semilogy(ivoltage, resistdrop, 'k.:', label = 'V Drop R')
            axs[1].set_xlabel('Voltage (V)')
            axs[1].set_ylabel('$R$ (Ω)')
            axs[1].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
            axs[1].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
            axs[1].grid(which = 'minor', color = 'lightgrey')
            axs[1].legend(frameon = True)
            
            axs[2].semilogy(voltage, pconst, 'kx-')
            axs[2].semilogy(voltage, 1.0/15*np.ones(len(voltage)), 'k--', linewidth = 1.0, label = 'Equivalent Limitation')
            axs[2].set_ylim(2.0e-3, 0.8)
            axs[2].set_xlabel('Voltage (V)')
            axs[2].set_ylabel('$P$')
            axs[2].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
            axs[2].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
            axs[2].grid(which = 'minor', color = 'lightgrey')
            axs[2].legend(frameon = True)
            
            axs[3].semilogy(voltage, fit_err, 'kx-')
            axs[3].set_xlabel('Voltage (V)')
            axs[3].set_ylabel('Fit Error')
            axs[3].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
            axs[3].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
            axs[3].grid(which = 'minor', color = 'lightgrey')
            
            axs[4].set_ylim(0, 1.0)
            axs[4].get_yaxis().set_ticks([0.25, 0.5, 0.75])
            axs[4].set_xlabel('Voltage (V)')
            axs[4].set_ylabel('$τ$ Span')
            axs[4].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs[4].grid(which = 'minor', color = 'lightgrey')
            axs[4].set_axisbelow(True)
            for j in range(nvolts):
                axs[4].fill(np.array([voltage[j]-dvolts[j]/2, voltage[j]-dvolts[j]/2, voltage[j]+dvolts[j]/2, voltage[j]+dvolts[j]/2]),
                            np.array([cap_min[j], cap_max[j], cap_max[j], cap_min[j]]),
                            color = 'grey', edgecolor = 'k', linestyle = '-')
            
            axs[5].plot(voltage, [i*1000/mass for i in dqdv], 'kx-')
            axs[5].set_xlabel('Voltage (V)')
            axs[5].set_ylabel('$dq/dV$ (mAh/gV)')
            axs[5].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axs[5].grid(which = 'minor', color = 'lightgrey')
    
        if export_fig:
            figname = dst / '{0} Summary.jpg'.format(cell_label)
            print(figname)
            plt.savefig(figname, bbox_inches = 'tight')
            
        plt.show()
        plt.close()

    def _spheres(self, X, logD, c_max):
        
        D = 10**logD
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis = 1)
        narr = np.repeat(n.reshape(len(n), 1), len(self.alphas), axis = 1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis = 0)
        
        return c/c_max + ((self.r**2)/(3*3600*n*D))*(1/5 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis = 1)))
    
    def _spheres_R_corr(self, X, logPDivD, c_max, logP):
        
        D = 10**(logP - logPDivD)
        P = 10**logP
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis = 1)
        narr = np.repeat(n.reshape(len(n), 1), len(self.alphas), axis = 1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis = 0)
        
        # Calculates inacessible capacity as 1 + tau if P/Q > 1 AND fcap is less than 0.05 of the largest fcap 
        # by setting n so that P = Q. Otherwise standard AMIDR equation.
        # This avoids the divergent region where tau = 0 but infinite summation error is amplified. 
        result = []
        for i in range(len(c)):
            if P>(3600*n[i]*D)/self.r**2 and c[i]/max(c) < 0.05:
                result.append(c[i]/c_max + 1)
            else:
                result.append(c[i]/c_max + ((self.r**2)/(3*3600*n[i]*D))*(1/5 - 2*(np.sum(np.exp(-a[i]*(carr[i]/c_max)*3600*narr[i]*D/self.r**2)/a[i]))) + P*self.r**2/(3600*n[i]*D))
        
        #return c/c_max + ((self.r**2)/(3*3600*n*D))*(1/5 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis = 1))) + self._dqdv*I*P/self._max_cap
        return result
    
    def _planes(self, X, logD, c_max):
        
        D = 10**logD
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis = 1)
        narr = np.repeat(n.reshape(len(c), 1), len(self.alphas), axis = 1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis = 0)
        
        return c/c_max + ((self.r**2)/(3600*n*D))*(1/3 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis = 1)))
    
    def insert_rate_cap(self, rate_cap):

        new_rate_cap = pd.read_csv(rate_cap, na_values = ['no info', '.'])
        
        self.nvolts = 1
        self.fcaps = [np.array(new_rate_cap['Capacity'])]
        self.eff_rates = [list(new_rate_cap['n in C/n'].values)]
        self.ir = [list(new_rate_cap['Crate'].values)]
        self.vlabels = ['inserted']
        self.avg_volts = [0]
        self.dqdv = [list(new_rate_cap['dqdv'].values)] 
        
class BINAVERAGE():
    
    def __init__(self, path, cells, matname, binsize = 0.025, mincapspan = 0.5, maxdqdVchange = 2, export_data = True, export_fig = True, parselabel = None, fitlabel = None):
        
        # Create new folder if necessary
        if parselabel == None:
            plabel = ''
        else:
            plabel = '-' + parselabel
        if fitlabel == None:
            flabel = ''
        else:
            flabel = '-' + fitlabel
        if export_data or export_fig:
            folder = Path(path) / (matname + plabel + flabel)
            if folder.is_dir() is False:
                folder.mkdir()
        
        # Generate input dataframes
        df = pd.DataFrame({})
        dfD = pd.DataFrame({})
        dfC = pd.DataFrame({})
        
        # Generate plot for individual cells
        fig, axs = plt.subplots(ncols = 2, nrows = 3, figsize = (6, 7.5), sharex = 'col', sharey = 'row',
        gridspec_kw = {'height_ratios': [2, 2, 1], 'hspace': 0.0, 'width_ratios': [1, 1], 'wspace': 0.0})
                
        # Find and read file data into dataframes
        for cell in cells:
            cellpath = Path(path) / cell
            for halfcyclepath in sorted(cellpath.iterdir(), reverse = True):
                if halfcyclepath.is_dir():
                    for fitfile in halfcyclepath.iterdir():
                        if fitfile.is_file() and "harge" + plabel + flabel + " Fitted" in str(fitfile):
                            if "Discharge" in str(fitfile):
                                print("Found discharge data for cell {}".format(Path(cellpath).name))
                                halfcycle = 'Discharge'
                            elif "Charge" in str(fitfile):
                                print("Found charge data for cell {}".format(Path(cellpath).name))
                                halfcycle = 'Charge'
                            else:
                                print("Found mislabeled fit file. Fit files should designate whether they contain charge or discharge data")
                                halfcycle = ''
                            dfnew = pd.read_excel(fitfile)
                            
                            # Remove datapoints where 
                            # dq/dV change between subsequent datapoint is greater than the max dqdV factor or
                            # the capacity span is less the minimum capacity span.
                            baddqdv = ~(((dfnew/dfnew.set_index(dfnew.index + 1))['dq/dV (mAh/gV)'] < maxdqdVchange) & \
                                        ((dfnew/dfnew.set_index(dfnew.index + 1))['dq/dV (mAh/gV)'] > 0) & \
                                        ((dfnew/dfnew.set_index(dfnew.index - 1))['dq/dV (mAh/gV)'] < maxdqdVchange) & \
                                        ((dfnew/dfnew.set_index(dfnew.index - 1))['dq/dV (mAh/gV)'] > 0) & \
                                        ((dfnew.set_index(dfnew.index + 1)/dfnew)['dq/dV (mAh/gV)'] < maxdqdVchange) & \
                                        ((dfnew.set_index(dfnew.index + 1)/dfnew)['dq/dV (mAh/gV)'] > 0) & \
                                        ((dfnew.set_index(dfnew.index - 1)/dfnew)['dq/dV (mAh/gV)'] < maxdqdVchange) & \
                                        ((dfnew.set_index(dfnew.index - 1)/dfnew)['dq/dV (mAh/gV)'] > 0))[1:-1]
                            badcapspan = dfnew['Cap Span'] < mincapspan
                            keep = ~(baddqdv|badcapspan)
                            
                            # Save good fits
                            df = pd.concat([df, dfnew[keep]], ignore_index = True)
                            if halfcycle == 'Discharge':
                                dfD = pd.concat([dfD, dfnew[keep]], ignore_index = True)
                                markerA = 'rx-'
                                markerB = 'r.:'
                                color = 'lightcoral'
                            elif halfcycle == 'Charge':
                                dfC = pd.concat([dfC, dfnew[keep]], ignore_index = True)
                                markerA = 'bx-'
                                markerB = 'b.:'
                                color = 'cornflowerblue'
                            else:
                                markerA = 'kx-'
                                markerB = 'k.:'
                                color = 'grey'
                            
                            # Plot individual cells with outliers removed
                            axs[0, 0].semilogy(dfnew[keep]['Voltage (V)'], dfnew[keep]['Dc (cm^2/s)'], markerA, markersize = 3)
                            axs[0, 0].semilogy(dfnew[keep]['Voltage (V)'], dfnew[keep]['Dt* (cm^2/s)'], markerB, markersize = 1.5)                            
                            axs[0, 1].semilogy(dfnew[keep]['SOC'], dfnew[keep]['Dc (cm^2/s)'], markerA, markersize = 3)
                            axs[0, 1].semilogy(dfnew[keep]['SOC'], dfnew[keep]['Dt* (cm^2/s)'], markerB, markersize = 1.5)
                            axs[1, 0].semilogy(dfnew[keep]['Initial Voltage (V)'], dfnew[keep]['micR (Ohmcm^2)'], markerA, markersize = 3)
                            axs[1, 1].semilogy(dfnew[keep]['Initial SOC'], dfnew[keep]['micR (Ohmcm^2)'], markerA, markersize = 3)
                            axs[2, 0].semilogy(dfnew[keep]['Initial Voltage (V)'], dfnew[keep]['dq/dV (mAh/gV)'], markerA, markersize = 3)
                            axs[2, 1].semilogy(dfnew[keep]['Initial SOC'], dfnew[keep]['dq/dV (mAh/gV)'], markerA, markersize = 3)
                            
                            # Plot individual cell outliers (dqdv)
                            axs[0, 0].semilogy(dfnew[baddqdv]['Voltage (V)'], dfnew[baddqdv]['Dc (cm^2/s)'], marker = '4', color = color, linestyle = 'None', markersize = 3)
                            #axs[0, 0].semilogy(dfnew[baddqdv]['Voltage (V)'], dfnew[baddqdv]['Dt* (cm^2/s)'], marker = '4', color = color, linestyle = 'None', markersize = 3)                            
                            axs[0, 1].semilogy(dfnew[baddqdv]['SOC'], dfnew[baddqdv]['Dc (cm^2/s)'], marker = '4', color = color, linestyle = 'None', markersize = 3)
                            #axs[0, 1].semilogy(dfnew[baddqdv]['SOC'], dfnew[baddqdv]['Dt* (cm^2/s)'], marker = '4', color = color, linestyle = 'None', markersize = 3)
                            axs[1, 0].semilogy(dfnew[baddqdv]['Initial Voltage (V)'], dfnew[baddqdv]['micR (Ohmcm^2)'], marker = '4', color = color, linestyle = 'None', markersize = 3)
                            axs[1, 1].semilogy(dfnew[baddqdv]['Initial SOC'], dfnew[baddqdv]['micR (Ohmcm^2)'], marker = '4', color = color, linestyle = 'None', markersize = 3)
                            axs[2, 0].semilogy(dfnew[baddqdv]['Initial Voltage (V)'], dfnew[baddqdv]['dq/dV (mAh/gV)'], marker = '4', color = color, linestyle = 'None', markersize = 3)
                            axs[2, 1].semilogy(dfnew[baddqdv]['Initial SOC'], dfnew[baddqdv]['dq/dV (mAh/gV)'], marker = '4', color = color, linestyle = 'None', markersize = 3)
        
                            # Plot individual cell outliers (capspan)
                            axs[0, 0].semilogy(dfnew[badcapspan]['Voltage (V)'], dfnew[badcapspan]['Dc (cm^2/s)'], marker = '3', color = color, linestyle = 'None', markersize = 3)
                            #axs[0, 0].semilogy(dfnew[badcapspan]['Voltage (V)'], dfnew[badcapspan]['Dt* (cm^2/s)'], marker = '3', color = color, linestyle = 'None', markersize = 3)                            
                            axs[0, 1].semilogy(dfnew[badcapspan]['SOC'], dfnew[badcapspan]['Dc (cm^2/s)'], marker = '3', color = color, linestyle = 'None', markersize = 3)
                            #axs[0, 1].semilogy(dfnew[badcapspan]['SOC'], dfnew[badcapspan]['Dt* (cm^2/s)'], marker = '3', color = color, linestyle = 'None', markersize = 3)
                            axs[1, 0].semilogy(dfnew[badcapspan]['Initial Voltage (V)'], dfnew[badcapspan]['micR (Ohmcm^2)'], marker = '3', color = color, linestyle = 'None', markersize = 3)
                            axs[1, 1].semilogy(dfnew[badcapspan]['Initial SOC'], dfnew[badcapspan]['micR (Ohmcm^2)'], marker = '3', color = color, linestyle = 'None', markersize = 3)
                            axs[2, 0].semilogy(dfnew[badcapspan]['Initial Voltage (V)'], dfnew[badcapspan]['dq/dV (mAh/gV)'], marker = '3', color = color, linestyle = 'None', markersize = 3)
                            axs[2, 1].semilogy(dfnew[badcapspan]['Initial SOC'], dfnew[badcapspan]['dq/dV (mAh/gV)'], marker = '3', color = color, linestyle = 'None', markersize = 3)
                            
                            # Create data files
                            if export_data:
                                filepath = folder / '{0} {1} {2}{3}{4} Filtered.xlsx'.format(cell, matname, halfcycle, plabel, flabel)
                                
                                print("Filtered data exporting to:\n" + str(filepath))
                                
                                writer = pd.ExcelWriter(filepath)
                                dfnew[keep].to_excel(writer, index = False)
                                writer.save()
                                writer.close()
        
        axs[0, 0].set_ylabel('$D$ (cm$\mathregular{^{2}}$/s)')
        axs[1, 0].set_ylabel('Max $ρ_{c}$ (Ωcm$\mathregular{^{2}}$)')
        axs[2, 0].set_ylabel('$dq/dV$ (mAh/gV)')
        axs[2, 0].set_xlabel('Voltage (V)')
        axs[2, 1].set_xlabel('Ion Saturation')
        
        axs[0, 0].semilogy([], [], 'rx-', label = 'Dch $D_{c}$')
        axs[0, 0].semilogy([], [], 'bx-', label = 'Ch $D_{c}$')
        axs[0, 0].semilogy([], [], 'r.:', label = 'Dch $D_{t}^{*}$')
        axs[0, 0].semilogy([], [], 'b.:', label = 'Ch $D_{t}^{*}$')
        axs[0, 0].legend(frameon = True, ncol = 2)
        
        axs[0, 1].semilogy([], [], marker = '4', color = 'grey', linestyle = 'None', label = '> Max $Δdq/dV$ ($D_{c}$)')
        axs[0, 1].semilogy([], [], marker = '3', color = 'grey', linestyle = 'None', label = '< Min $τ$ Span ($D_{c}$)')
        axs[0, 1].legend(frameon = True)
        axs[0, 1].invert_xaxis()
        
        axs[1, 0].semilogy([], [], 'rx-', label = 'Dch')
        axs[1, 0].semilogy([], [], 'bx-', label = 'Ch')
        axs[1, 0].legend(frameon = True)
        
        axs[1, 1].semilogy([], [], marker = '4', color = 'grey', linestyle = 'None', label = '> Max $Δdq/dV$')
        axs[1, 1].semilogy([], [], marker = '3', color = 'grey', linestyle = 'None', label = '< Min $τ$ Span')
        axs[1, 1].legend(frameon = True)
        
        axs[0, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[0, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[1, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[2, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        
        axs[0, 0].grid(which = 'minor', color = 'lightgrey')
        axs[0, 1].grid(which = 'minor', color = 'lightgrey')
        axs[1, 0].grid(which = 'minor', color = 'lightgrey')
        axs[1, 1].grid(which = 'minor', color = 'lightgrey')
        axs[2, 0].grid(which = 'minor', color = 'lightgrey')
        axs[2, 1].grid(which = 'minor', color = 'lightgrey')
        
        if export_fig:
            figname = folder / '{0}{1}{2} Individual Cells.jpg'.format(matname, plabel, flabel)
            print(figname)
            plt.savefig(figname, bbox_inches = 'tight')
            
        plt.show()
        plt.close()
        
        # Establish bins and output dataframes
        firstbinnum = int(min(min(df['Voltage (V)']), min(df['Initial Voltage (V)']))//binsize)
        lastbinnum = int(max(max(df['Voltage (V)']), max(df['Initial Voltage (V)']))//binsize)

        dfO = pd.DataFrame({'Voltage (V)': np.linspace(firstbinnum*binsize, lastbinnum*binsize, abs(firstbinnum - lastbinnum) + 1) + binsize/2})
        dfOD = pd.DataFrame({'Voltage (V)': np.linspace(firstbinnum*binsize, lastbinnum*binsize, abs(firstbinnum - lastbinnum) + 1) + binsize/2})
        dfOC = pd.DataFrame({'Voltage (V)': np.linspace(firstbinnum*binsize, lastbinnum*binsize, abs(firstbinnum - lastbinnum) + 1) + binsize/2})

        # Fill output dataframes
        for i in dfO.index:
            
            # Select data points with voltages within the bin
            binSel = (df['Voltage (V)'] >= dfO['Voltage (V)'][i] - binsize/2) & \
                     (df['Voltage (V)'] < dfO['Voltage (V)'][i] + binsize/2)
            binSelD = (dfD['Voltage (V)'] >= dfO['Voltage (V)'][i] - binsize/2) & \
                      (dfD['Voltage (V)'] < dfO['Voltage (V)'][i] + binsize/2)
            binSelC = (dfC['Voltage (V)'] >= dfO['Voltage (V)'][i] - binsize/2) & \
                      (dfC['Voltage (V)'] < dfO['Voltage (V)'][i] + binsize/2)
                     
            # Select data points with initial voltages within the bin (for resistances which are plotted at initial cap)
            ibinSel = (df['Initial Voltage (V)'] >= dfO['Voltage (V)'][i] - binsize/2) & \
                      (df['Initial Voltage (V)'] < dfO['Voltage (V)'][i] + binsize/2)
            ibinSelD = (dfD['Initial Voltage (V)'] >= dfO['Voltage (V)'][i] - binsize/2) & \
                       (dfD['Initial Voltage (V)'] < dfO['Voltage (V)'][i] + binsize/2)
            ibinSelC = (dfC['Initial Voltage (V)'] >= dfO['Voltage (V)'][i] - binsize/2) & \
                       (dfC['Initial Voltage (V)'] < dfO['Voltage (V)'][i] + binsize/2)
            
            # Calculate geometric mean and standard deviation for diffusivities and resistances
            dfO.loc[i, 'Dc (cm^2/s)'] = np.exp(np.log(df['Dc (cm^2/s)'][binSel]).mean(axis = 0))
            dfO.loc[i, 'Dc geoSTD'] = np.exp(np.log(df['Dc (cm^2/s)'][binSel]).std(axis = 0))
            dfO.loc[i, 'Dt* (cm^2/s)'] = np.exp(np.log(df['Dt* (cm^2/s)'][binSel]).mean(axis = 0))
            dfO.loc[i, 'Dt* geoSTD'] = np.exp(np.log(df['Dt* (cm^2/s)'][binSel]).std(axis = 0))
            dfO.loc[i, 'Rfit (Ohm)'] = np.exp(np.log(df['Rfit (Ohm)'][ibinSel]).mean(axis = 0))
            dfO.loc[i, 'Rfit geoSTD'] = np.exp(np.log(df['Rfit (Ohm)'][ibinSel]).std(axis = 0))
            dfO.loc[i, 'micR (Ohmcm^2)'] = np.exp(np.log(df['micR (Ohmcm^2)'][ibinSel]).mean(axis = 0))
            dfO.loc[i, 'micR geoSTD'] = np.exp(np.log(df['micR (Ohmcm^2)'][ibinSel]).std(axis = 0))
            dfO.loc[i, 'Rdrop (Ohm)'] = np.exp(np.log(df['Rdrop (Ohm)'][ibinSel]).mean(axis = 0))
            dfO.loc[i, 'Rdrop geoSTD'] = np.exp(np.log(df['Rdrop (Ohm)'][ibinSel]).std(axis = 0))
            dfOD.loc[i, 'Dc (cm^2/s)'] = np.exp(np.log(dfD['Dc (cm^2/s)'][binSelD]).mean(axis = 0))
            dfOD.loc[i, 'Dc geoSTD'] = np.exp(np.log(dfD['Dc (cm^2/s)'][binSelD]).std(axis = 0))
            dfOD.loc[i, 'Dt* (cm^2/s)'] = np.exp(np.log(dfD['Dt* (cm^2/s)'][binSelD]).mean(axis = 0))
            dfOD.loc[i, 'Dt* geoSTD'] = np.exp(np.log(dfD['Dt* (cm^2/s)'][binSelD]).std(axis = 0))
            dfOD.loc[i, 'Rfit (Ohm)'] = np.exp(np.log(dfD['Rfit (Ohm)'][ibinSelD]).mean(axis = 0))
            dfOD.loc[i, 'Rfit geoSTD'] = np.exp(np.log(dfD['Rfit (Ohm)'][ibinSelD]).std(axis = 0))
            dfOD.loc[i, 'micR (Ohmcm^2)'] = np.exp(np.log(dfD['micR (Ohmcm^2)'][ibinSelD]).mean(axis = 0))
            dfOD.loc[i, 'micR geoSTD'] = np.exp(np.log(dfD['micR (Ohmcm^2)'][ibinSelD]).std(axis = 0))
            dfOD.loc[i, 'Rdrop (Ohm)'] = np.exp(np.log(dfD['Rdrop (Ohm)'][ibinSelD]).mean(axis = 0))
            dfOD.loc[i, 'Rdrop geoSTD'] = np.exp(np.log(dfD['Rdrop (Ohm)'][ibinSelD]).std(axis = 0))
            dfOC.loc[i, 'Dc (cm^2/s)'] = np.exp(np.log(dfC['Dc (cm^2/s)'][binSelC]).mean(axis = 0))
            dfOC.loc[i, 'Dc geoSTD'] = np.exp(np.log(dfC['Dc (cm^2/s)'][binSelC]).std(axis = 0))
            dfOC.loc[i, 'Dt* (cm^2/s)'] = np.exp(np.log(dfC['Dt* (cm^2/s)'][binSelC]).mean(axis = 0))
            dfOC.loc[i, 'Dt* geoSTD'] = np.exp(np.log(dfC['Dt* (cm^2/s)'][binSelC]).std(axis = 0))
            dfOC.loc[i, 'Rfit (Ohm)'] = np.exp(np.log(dfC['Rfit (Ohm)'][ibinSelC]).mean(axis = 0))
            dfOC.loc[i, 'Rfit geoSTD'] = np.exp(np.log(dfC['Rfit (Ohm)'][ibinSelC]).std(axis = 0))
            dfOC.loc[i, 'micR (Ohmcm^2)'] = np.exp(np.log(dfC['micR (Ohmcm^2)'][ibinSelC]).mean(axis = 0))
            dfOC.loc[i, 'micR geoSTD'] = np.exp(np.log(dfC['micR (Ohmcm^2)'][ibinSelC]).std(axis = 0))
            dfOC.loc[i, 'Rdrop (Ohm)'] = np.exp(np.log(dfC['Rdrop (Ohm)'][ibinSelC]).mean(axis = 0))
            dfOC.loc[i, 'Rdrop geoSTD'] = np.exp(np.log(dfC['Rdrop (Ohm)'][ibinSelC]).std(axis = 0))
            
            # Calculate arithmetic mean and standard deviation for all else
            dfO.loc[i, 'Voltage STD'] = df['Voltage (V)'][binSel|ibinSel].std(axis = 0)
            dfO.loc[i, 'SOC'] = df['SOC'][binSel|ibinSel].mean(axis = 0)
            dfO.loc[i, 'SOC STD'] = df['SOC'][binSel|ibinSel].std(axis = 0)
            dfO.loc[i, 'dq/dV (mAh/gV)'] = df['dq/dV (mAh/gV)'][binSel].mean(axis = 0)
            dfO.loc[i, 'dq/dV STD'] = df['dq/dV (mAh/gV)'][binSel].std(axis = 0)
            dfO.loc[i, 'Cap Span'] = df['Cap Span'][binSel].mean(axis = 0)
            dfO.loc[i, 'Cap Span STD'] = df['Cap Span'][binSel].std(axis = 0)
            dfO.loc[i, 'Fit Error'] = df['Fit Error'][binSel].mean(axis = 0)
            dfO.loc[i, 'Fit Error STD'] = df['Fit Error'][binSel].std(axis = 0)
            dfOD.loc[i, 'Voltage STD'] = dfD['Voltage (V)'][binSelD|ibinSelD].std(axis = 0)
            dfOD.loc[i, 'SOC'] = dfD['SOC'][binSelD|ibinSelD].mean(axis = 0)
            dfOD.loc[i, 'SOC STD'] = dfD['SOC'][binSelD|ibinSelD].std(axis = 0)
            dfOD.loc[i, 'dq/dV (mAh/gV)'] = dfD['dq/dV (mAh/gV)'][binSelD].mean(axis = 0)
            dfOD.loc[i, 'dq/dV STD'] = dfD['dq/dV (mAh/gV)'][binSelD].std(axis = 0)
            dfOD.loc[i, 'Cap Span'] = dfD['Cap Span'][binSelD].mean(axis = 0)
            dfOD.loc[i, 'Cap Span STD'] = dfD['Cap Span'][binSelD].std(axis = 0)
            dfOD.loc[i, 'Fit Error'] = dfD['Fit Error'][binSelD].mean(axis = 0)
            dfOD.loc[i, 'Fit Error STD'] = dfD['Fit Error'][binSelD].std(axis = 0)
            dfOC.loc[i, 'Voltage STD'] = dfC['Voltage (V)'][binSelC|ibinSelC].std(axis = 0)
            dfOC.loc[i, 'SOC'] = dfC['SOC'][binSelC|ibinSelC].mean(axis = 0)
            dfOC.loc[i, 'SOC STD'] = dfC['SOC'][binSelC|ibinSelC].std(axis = 0)
            dfOC.loc[i, 'dq/dV (mAh/gV)'] = dfC['dq/dV (mAh/gV)'][binSelC].mean(axis = 0)
            dfOC.loc[i, 'dq/dV STD'] = dfC['dq/dV (mAh/gV)'][binSelC].std(axis = 0)
            dfOC.loc[i, 'Cap Span'] = dfC['Cap Span'][binSelC].mean(axis = 0)
            dfOC.loc[i, 'Cap Span STD'] = dfC['Cap Span'][binSelC].std(axis = 0)
            dfOC.loc[i, 'Fit Error'] = dfC['Fit Error'][binSelC].mean(axis = 0)
            dfOC.loc[i, 'Fit Error STD'] = dfC['Fit Error'][binSelC].std(axis = 0)
            
        # Plot bin averaged charge and discharge
        fig, axs = plt.subplots(ncols = 2, nrows = 3, figsize = (6, 7.5), sharex = 'col', sharey = 'row',
        gridspec_kw = {'height_ratios': [2, 2, 1], 'hspace': 0.0, 'width_ratios': [1, 1], 'wspace': 0.0})
        
        axs[0, 1].semilogy([], [], 'r-', label = 'Dch $D_{c}$')
        axs[0, 1].semilogy([], [], 'b-', label = 'Ch $D_{c}$')
        axs[0, 1].semilogy([], [], 'r:', label = 'Dch $D_{t}^{*}$')
        axs[0, 1].semilogy([], [], 'b:', label = 'Ch $D_{t}^{*}$')
        axs[0, 1].legend(frameon = True, ncol = 2)
        
        axs[1, 1].semilogy([], [], 'r-', label = 'Dch')
        axs[1, 1].semilogy([], [], 'b-', label = 'Ch')
        axs[1, 1].legend(frameon = True)
        
        yerrDDc = [dfOD['Dc (cm^2/s)']*(1 - 1/dfOD['Dc geoSTD']), dfOD['Dc (cm^2/s)']*(dfOD['Dc geoSTD'] - 1)]
        yerrDDt = [dfOD['Dt* (cm^2/s)']*(1 - 1/dfOD['Dt* geoSTD']), dfOD['Dt* (cm^2/s)']*(dfOD['Dt* geoSTD'] - 1)]
        axs[0, 0].errorbar(dfOD['Voltage (V)'], dfOD['Dc (cm^2/s)'], yerr = yerrDDc, fmt = 'r-')
        axs[0, 0].errorbar(dfOD['Voltage (V)'], dfOD['Dt* (cm^2/s)'], yerr = yerrDDt, fmt = 'r:', elinewidth = 0.5, markeredgewidth = 0.5)
        axs[0, 1].errorbar(dfOD['SOC'], dfOD['Dc (cm^2/s)'], yerr = yerrDDc, fmt = 'r-')
        axs[0, 1].errorbar(dfOD['SOC'], dfOD['Dt* (cm^2/s)'], yerr = yerrDDt, fmt = 'r:', elinewidth = 0.5, markeredgewidth = 0.5)
        yerrDmicR = [dfOD['micR (Ohmcm^2)']*(1 - 1/dfOD['micR geoSTD']), dfOD['micR (Ohmcm^2)']*(dfOD['micR geoSTD'] - 1)]
        axs[1, 0].errorbar(dfOD['Voltage (V)'], dfOD['micR (Ohmcm^2)'], yerr = yerrDmicR, fmt = 'r-')
        axs[1, 1].errorbar(dfOD['SOC'], dfOD['micR (Ohmcm^2)'], yerr = yerrDmicR, fmt = 'r-')
        axs[2, 0].errorbar(dfOD['Voltage (V)'], dfOD['dq/dV (mAh/gV)'], yerr = dfOD['dq/dV STD'], fmt = 'r-')
        axs[2, 1].errorbar(dfOD['SOC'], dfOD['dq/dV (mAh/gV)'], yerr = dfOD['dq/dV STD'], fmt = 'r-')
        
        yerrCDc = [dfOC['Dc (cm^2/s)']*(1 - 1/dfOC['Dc geoSTD']), dfOC['Dc (cm^2/s)']*(dfOC['Dc geoSTD'] - 1)]
        yerrCDt = [dfOC['Dt* (cm^2/s)']*(1 - 1/dfOC['Dt* geoSTD']), dfOC['Dt* (cm^2/s)']*(dfOC['Dt* geoSTD'] - 1)]
        axs[0, 0].errorbar(dfOC['Voltage (V)'], dfOC['Dc (cm^2/s)'], yerr = yerrCDc, fmt = 'b-')
        axs[0, 0].errorbar(dfOC['Voltage (V)'], dfOC['Dt* (cm^2/s)'], yerr = yerrCDt, fmt = 'b:', elinewidth = 0.5, markeredgewidth = 0.5)
        axs[0, 1].errorbar(dfOC['SOC'], dfOC['Dc (cm^2/s)'], yerr = yerrCDc, fmt = 'b-')
        axs[0, 1].errorbar(dfOC['SOC'], dfOC['Dt* (cm^2/s)'], yerr = yerrCDt, fmt = 'b:', elinewidth = 0.5, markeredgewidth = 0.5)
        yerrCmicR = [dfOC['micR (Ohmcm^2)']*(1 - 1/dfOC['micR geoSTD']), dfOC['micR (Ohmcm^2)']*(dfOC['micR geoSTD'] - 1)]
        axs[1, 0].errorbar(dfOC['Voltage (V)'], dfOC['micR (Ohmcm^2)'],  yerr = yerrCmicR, fmt = 'b-')
        axs[1, 1].errorbar(dfOC['SOC'], dfOC['micR (Ohmcm^2)'], yerr = yerrCmicR, fmt = 'b-')
        axs[2, 0].errorbar(dfOC['Voltage (V)'], dfOC['dq/dV (mAh/gV)'], yerr = dfOC['dq/dV STD'], fmt = 'b-')
        axs[2, 1].errorbar(dfOC['SOC'], dfOC['dq/dV (mAh/gV)'], yerr = dfOC['dq/dV STD'], fmt = 'b-')
        
        axs[0, 0].set_ylabel('$D$ (cm$\mathregular{^{2}}$/s)')
        axs[1, 0].set_ylabel('Max $ρ_{c}$ (Ωcm$\mathregular{^{2}}$)')
        axs[2, 0].set_ylabel('$dq/dV$ (mAh/gV)')
        axs[2, 0].set_xlabel('Voltage (V)')
        axs[2, 1].set_xlabel('Ion Saturation')
        
        axs[0, 1].invert_xaxis()
        
        axs[0, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[0, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[1, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[2, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        
        axs[0, 0].grid(which = 'minor', color = 'lightgrey')
        axs[0, 1].grid(which = 'minor', color = 'lightgrey')
        axs[1, 0].grid(which = 'minor', color = 'lightgrey')
        axs[1, 1].grid(which = 'minor', color = 'lightgrey')
        axs[2, 0].grid(which = 'minor', color = 'lightgrey')
        axs[2, 1].grid(which = 'minor', color = 'lightgrey')
        
        if export_fig:
            figname = folder / '{0}{1}{2} Ch vs Dch.jpg'.format(matname, plabel, flabel)
            print(figname)
            plt.savefig(figname, bbox_inches = 'tight')
            
        plt.show()
        plt.close()
        
        # Plot bin averaged all
        fig, axs = plt.subplots(ncols = 2, nrows = 3, figsize = (6, 7.5), sharex = 'col', sharey = 'row',
        gridspec_kw = {'height_ratios': [2, 2, 1], 'hspace': 0.0, 'width_ratios': [1, 1], 'wspace': 0.0})
        
        axs[0, 1].semilogy([], [], 'k-', label = '$D_{c}$')
        axs[0, 1].semilogy([], [], 'k:', label = '$D_{t}^{*}$')
        axs[0, 1].legend(frameon = True, ncol = 2)
        
        axs[1, 0].semilogy([], [])
        
        yerrDc = [dfO['Dc (cm^2/s)']*(1 - 1/dfO['Dc geoSTD']), dfO['Dc (cm^2/s)']*(dfO['Dc geoSTD'] - 1)]
        yerrDt = [dfO['Dt* (cm^2/s)']*(1 - 1/dfO['Dt* geoSTD']), dfO['Dt* (cm^2/s)']*(dfO['Dt* geoSTD'] - 1)]
        axs[0, 0].errorbar(dfO['Voltage (V)'], dfO['Dc (cm^2/s)'], yerr = yerrDc, fmt = 'k-')
        axs[0, 0].errorbar(dfO['Voltage (V)'], dfO['Dt* (cm^2/s)'], yerr = yerrDt, fmt = 'k:', elinewidth = 0.5, markeredgewidth = 0.5)
        axs[0, 1].errorbar(dfO['SOC'], dfO['Dc (cm^2/s)'], yerr = yerrDc, fmt = 'k-')
        axs[0, 1].errorbar(dfO['SOC'], dfO['Dt* (cm^2/s)'], yerr = yerrDt, fmt = 'k:', elinewidth = 0.5, markeredgewidth = 0.5)
        yerrmicR = [dfO['micR (Ohmcm^2)']*(1 - 1/dfO['micR geoSTD']), dfO['micR (Ohmcm^2)']*(dfO['micR geoSTD'] - 1)]
        axs[1, 0].errorbar(dfO['Voltage (V)'], dfO['micR (Ohmcm^2)'], yerr = yerrmicR, fmt = 'k-')
        axs[1, 1].errorbar(dfO['SOC'], dfO['micR (Ohmcm^2)'], yerr = yerrmicR, fmt = 'k-')
        axs[2, 0].errorbar(dfO['Voltage (V)'], dfO['dq/dV (mAh/gV)'], yerr = dfO['dq/dV STD'], fmt = 'k-')
        axs[2, 1].errorbar(dfO['SOC'], dfO['dq/dV (mAh/gV)'], yerr = dfO['dq/dV STD'], fmt = 'k-')
        
        axs[0, 0].set_ylabel('$D$ (cm$\mathregular{^{2}}$/s)')
        axs[1, 0].set_ylabel('Max $ρ_{c}$ (Ωcm$\mathregular{^{2}}$)')
        axs[2, 0].set_ylabel('$dq/dV$ (mAh/gV)')
        axs[2, 0].set_xlabel('Voltage (V)')
        axs[2, 1].set_xlabel('Ion Saturation')
        
        axs[0, 1].invert_xaxis()
        
        axs[0, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[0, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[1, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[2, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())

        axs[0, 0].grid(which = 'minor', color = 'lightgrey')
        axs[0, 1].grid(which = 'minor', color = 'lightgrey')
        axs[1, 0].grid(which = 'minor', color = 'lightgrey')
        axs[1, 1].grid(which = 'minor', color = 'lightgrey')  
        axs[2, 0].grid(which = 'minor', color = 'lightgrey')
        axs[2, 1].grid(which = 'minor', color = 'lightgrey')  
        
        if export_fig:
            figname = folder / '{0}{1}{2} All.jpg'.format(matname, plabel, flabel)
            print(figname)
            plt.savefig(figname, bbox_inches = 'tight')
            
        plt.show()
        plt.close()
        
        # Create data files
        if export_data:
            dfO = dfO[['Voltage (V)', 'Voltage STD', 'SOC', 'SOC STD', 'Dc (cm^2/s)', 'Dc geoSTD', 'Dt* (cm^2/s)', 'Dt* geoSTD', 'dq/dV (mAh/gV)', 'dq/dV STD', \
                   'Rfit (Ohm)', 'Rfit geoSTD', 'micR (Ohmcm^2)', 'micR geoSTD', 'Rdrop (Ohm)', 'Rdrop geoSTD', 'Cap Span', 'Cap Span STD', 'Fit Error', 'Fit Error STD']]
            dfOD = dfOD[['Voltage (V)', 'Voltage STD', 'SOC', 'SOC STD', 'Dc (cm^2/s)', 'Dc geoSTD', 'Dt* (cm^2/s)', 'Dt* geoSTD', 'dq/dV (mAh/gV)', 'dq/dV STD', \
                     'Rfit (Ohm)', 'Rfit geoSTD', 'micR (Ohmcm^2)', 'micR geoSTD', 'Rdrop (Ohm)', 'Rdrop geoSTD', 'Cap Span', 'Cap Span STD', 'Fit Error', 'Fit Error STD']]
            dfOC = dfOC[['Voltage (V)', 'Voltage STD', 'SOC', 'SOC STD', 'Dc (cm^2/s)', 'Dc geoSTD', 'Dt* (cm^2/s)', 'Dt* geoSTD', 'dq/dV (mAh/gV)', 'dq/dV STD', \
                     'Rfit (Ohm)', 'Rfit geoSTD', 'micR (Ohmcm^2)', 'micR geoSTD', 'Rdrop (Ohm)', 'Rdrop geoSTD', 'Cap Span', 'Cap Span STD', 'Fit Error', 'Fit Error STD']]
            
            filepath = folder / '{0}{1}{2} ({3}).xlsx'.format(matname, plabel, flabel, ', '.join(cells))
            
            print("Bin averaged data exporting to:\n" + str(filepath))
            
            writer = pd.ExcelWriter(filepath)
            dfO.to_excel(writer, sheet_name = 'All', index = False)
            dfOD.to_excel(writer, sheet_name = 'Discharge', index = False)
            dfOC.to_excel(writer, sheet_name = 'Charge', index = False)
            writer.save()
            writer.close()
            
class MATCOMPARE():
    
    def __init__(self, path, mats, export_data = None, export_fig = True):
        
        if not(export_data is None): ('There is no figure to export. Feel free to neglect this argument.')
        
        if len(mats) > 4:
            print('Can only compare up to 4 materials at once. Additional materials will not be plotted.')  
        
        folder = Path(path)
        
        # Generate plot for individual cells
        fig, axs = plt.subplots(ncols = 2, nrows = 3, figsize = (6, 7.5), sharex = 'col', sharey = 'row',
        gridspec_kw = {'height_ratios': [2, 2, 1], 'hspace': 0.0, 'width_ratios': [1, 1], 'wspace': 0.0})
                
        axs[0, 1].semilogy([], [], 'k-', label = '$D_{c}$')
        axs[0, 1].semilogy([], [], 'k:', label = '$D_{t}^{*}$')
        axs[0, 1].legend(frameon = True, ncol = 2)
        
        axs[0, 0].set_ylabel('$D$ (cm$\mathregular{^{2}}$/s)')
        axs[1, 0].set_ylabel('Max $ρ_{c}$ (Ωcm$\mathregular{^{2}}$)')
        axs[2, 0].set_ylabel('$dq/dV$ (mAh/gV)')
        axs[2, 0].set_xlabel('Voltage (V)')
        axs[2, 1].set_xlabel('Ion Saturation')
        
        axs[0, 1].invert_xaxis()
        
        axs[0, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[0, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[0, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axs[1, 0].yaxis.set_minor_locator(ticker.LogLocator(subs = np.arange(1.0, 10.0) * 0.1, numticks = 10))
        axs[1, 0].yaxis.set_major_locator(ticker.LogLocator(numticks = 10))
        axs[2, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())

        axs[0, 0].grid(which = 'minor', color = 'lightgrey')
        axs[0, 1].grid(which = 'minor', color = 'lightgrey')
        axs[1, 0].grid(which = 'minor', color = 'lightgrey')
        axs[1, 1].grid(which = 'minor', color = 'lightgrey') 
        axs[2, 0].grid(which = 'minor', color = 'lightgrey')
        axs[2, 1].grid(which = 'minor', color = 'lightgrey') 
        
        colors = ['r', 'g', 'b', 'k']
        
        # Find and read file data into dataframes
        i = 0
        for mat in mats:
            matpath = folder / mat
            for filepath in matpath.iterdir():
                if mat + ' (' in str(filepath):
                    if i > 3:
                        break
                    
                    print("Found data for material: {}".format(filepath))
                    df = pd.read_excel(filepath, sheet_name = 'All')
                    
                    # Plot dataframe
                    yerrDc = [df['Dc (cm^2/s)']*(1 - 1/df['Dc geoSTD']), df['Dc (cm^2/s)']*(df['Dc geoSTD'] - 1)]
                    yerrDt = [df['Dt* (cm^2/s)']*(1 - 1/df['Dt* geoSTD']), df['Dt* (cm^2/s)']*(df['Dt* geoSTD'] - 1)]
                    axs[0, 0].errorbar(df['Voltage (V)'], df['Dc (cm^2/s)'], yerr = yerrDc, fmt = '-', color = colors[i])
                    axs[0, 0].errorbar(df['Voltage (V)'], df['Dt* (cm^2/s)'], yerr = yerrDt, fmt = ':', color = colors[i], elinewidth = 0.5, markeredgewidth = 0.5)
                    axs[0, 1].errorbar(df['SOC'], df['Dc (cm^2/s)'], yerr = yerrDc, fmt = '-', color = colors[i])
                    axs[0, 1].errorbar(df['SOC'], df['Dt* (cm^2/s)'], yerr = yerrDt, fmt = ':', color = colors[i], elinewidth = 0.5, markeredgewidth = 0.5)
                    yerrmicR = [df['micR (Ohmcm^2)']*(1 - 1/df['micR geoSTD']), df['micR (Ohmcm^2)']*(df['micR geoSTD'] - 1)]
                    axs[1, 0].errorbar(df['Voltage (V)'], df['micR (Ohmcm^2)'], yerr = yerrmicR, color = colors[i], fmt = '-')
                    axs[1, 1].errorbar(df['SOC'], df['micR (Ohmcm^2)'], yerr = yerrmicR, fmt = '-', color = colors[i])
                    axs[1, 1].semilogy([], [], color = colors[i], label = mat)
                    axs[2, 0].errorbar(df['Voltage (V)'], df['dq/dV (mAh/gV)'], yerr = df['dq/dV STD'], color = colors[i], fmt = '-')
                    axs[2, 1].errorbar(df['SOC'], df['dq/dV (mAh/gV)'], yerr = df['dq/dV STD'], fmt = '-', color = colors[i])
                    
                    i = i + 1
                    
        axs[1, 1].legend(frameon = True)            
        
        if export_fig:
            figname = folder / 'Material Comparison ({0}).jpg'.format(', '.join(mats))
            print(figname)
            plt.savefig(figname, bbox_inches = 'tight')
            
        plt.show()
        plt.close()