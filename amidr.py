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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


RATES = np.array([0.01, 0.05, 0.1, 0.2, 1/3, 0.5, 1, 2, 2.5, 5, 10, 20, 40, 80,
                  160, 320, 640, 1280])

COLUMNS = ['Time', 'Cycle', 'Step', 'Current', 'Potential', 'Capacity', 'Prot_step']
UNITS = ['(h)', None, None, '(mA)', '(V)', '(mAh)', None]

SHAPES = ['sphere', 'plane']

class BIOCONVERT():
    
    def __init__(self, path, form_files, d_files, c_files, name, export_fig=True):
        
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
            csvHeader = '[Summary]\nCell: ' + name + '\nFirst Protocol: ' + protName \
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
                dfTempForm = pd.read_csv(formFileLoc, skiprows=hlinenum, sep = '\t', encoding_errors='replace')
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
            
            dfForm.rename(columns={'mode':'Step Type', 
                                   'time/s':'Run Time (h)', 
                                   'I/mA':'Current (A)', 
                                   'Ewe-Ece/V':'Potential vs. Counter (V)', 
                                   'Ewe/V':'Potential (V)', 
                                   '(Q-Qo)/mA.h':'Capacity (Ah)', 
                                   'Ns':'Step Number'}, 
                          inplace=True)
            
            # Add header info into form file
            pathFileForm = Path(path) / (name + ' Form.csv')
            with open(pathFileForm, 'w') as f:
                f.write(csvHeader)
            
            # Add data into form file
            dfForm.to_csv(pathFileForm, mode='a', index=False)
        
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
                dfTempD = pd.read_csv(dFileLoc, skiprows=hlinenum, sep = '\t', encoding_errors='replace')
                dfTempD = dfTempD[['mode', 'time/s', 'I/mA', 'Ewe-Ece/V', 'Ewe/V', '(Q-Qo)/mA.h', 'Ns', 'control/mA']]
    
                # Convert 0A CC to NVX rest steps and trim off control I column
                dfTempD['mode'].mask(dfTempD['control/mA'] == 0, 0, inplace=True)
                dfTempD['mode'].mask(dfTempD['mode'] == 3, 0, inplace=True)
                dfTempD.drop(columns=['control/mA'], inplace=True)
 
                dfTempDSteps = dfTempD.drop_duplicates(subset=['Ns'], ignore_index=True)
                
                # Detect if test ended prematurely
                if dfTempDSteps['mode'].iloc[-1] != 0:
                    
                    # Add dummy step at end to prevent overindexing
                    dummyD = dfTempD.loc[dfTempD.index[-1]:dfTempD.index[-1]].copy()
                    dummyD['Ns'] = dummyD['Ns'] + 1
                    dummyD['mode'] = 0
                    dfTempD = pd.concat([dfTempD, dummyD], ignore_index=True)
                    dfTempDSteps = pd.concat([dfTempDSteps, dummyD], ignore_index=True)
                    
                # Iterate over each pulse starting with their preceeding OCV V rest step    
                for i in dfTempDSteps.index:
                    if i != dfTempDSteps.index[-1]:
                        if dfTempDSteps['mode'][i] == 0 and dfTempDSteps['mode'][i+1] != 0:
                            
                            # Average together all points of the rest step before a pulse into 1 point (OCV V)
                            ocvFinSel = dfTempD['Ns'] == dfTempDSteps['Ns'][i]
                            ocvFinVals = dfTempD[ocvFinSel].mean(axis=0)
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
                                        CCPointVals = dfTempD.loc[pulseInd + 1:pulseInd + nAvg].mean(axis=0)
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
                            ocvFinVals = dfTempD[ocvFinSel].mean(axis=0)
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
                newDSteps = dfTempD.drop_duplicates(subset=['Ns'], ignore_index=True)
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
            
            dfD.rename(columns={'mode':'Step Type', 
                                'time/s':'Run Time (h)', 
                                'I/mA':'Current (A)', 
                                'Ewe-Ece/V':'Potential vs. Counter (V)', 
                                'Ewe/V':'Potential (V)', 
                                '(Q-Qo)/mA.h':'Capacity (Ah)', 
                                'Ns':'Step Number'}, 
                       inplace=True)
                
            # Add header info into d file
            pathFileD = Path(path) / (name + ' Discharge.csv')
            with open(pathFileD, 'w') as f:
                f.write(csvHeader)
            
            # Add data into d file
            dfD.to_csv(pathFileD, mode='a', index=False)
            
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
                dfTempC = pd.read_csv(cFileLoc, skiprows=hlinenum, sep = '\t', encoding_errors='replace')
                dfTempC = dfTempC[['mode', 'time/s', 'I/mA', 'Ewe-Ece/V', 'Ewe/V', '(Q-Qo)/mA.h', 'Ns', 'control/mA']]
    
                # Convert 0A CC to NVX rest steps and trim off control I column
                dfTempC['mode'].mask(dfTempC['control/mA'] == 0, 0, inplace=True)
                dfTempC.drop(columns=['control/mA'], inplace=True)
                
                dfTempCSteps = dfTempC.drop_duplicates(subset=['Ns'], ignore_index=True)     
                
                # Detect if test ended prematurely
                if dfTempCSteps['mode'].iloc[-1] != 0:
                    
                    # Add dummy step at end to prevent overindexing
                    dummyC = dfTempC.loc[dfTempC.index[-1]:dfTempC.index[-1]].copy()
                    dummyC['Ns'] = dummyC['Ns'] + 1
                    dummyC['mode'] = 0
                    dfTempC = pd.concat([dfTempC, dummyC], ignore_index=True)
                    dfTempCSteps = pd.concat([dfTempCSteps, dummyC], ignore_index=True)
                
                # Iterate over each pulse starting with their preceeding OCV V rest step 
                for i in dfTempCSteps.index:
                    if i != dfTempCSteps.index[-1]:
                        if dfTempCSteps['mode'][i] == 0 and dfTempCSteps['mode'][i+1] != 0:
                            
                            # Average together all points of the rest step before a pulse into 1 point (OCV V)
                            ocvFinSel = dfTempC['Ns'] == dfTempCSteps['Ns'][i]
                            ocvFinVals = dfTempC[ocvFinSel].mean(axis=0)
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
                                        CCPointVals = dfTempC.loc[pulseInd + 1:pulseInd + nAvg].mean(axis=0)
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
                            ocvFinVals = dfTempC[ocvFinSel].mean(axis=0)
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
                newCSteps = dfTempC.drop_duplicates(subset=['Ns'], ignore_index=True)
                dfTempC['Ns'].replace(newCSteps['Ns'].values, newCSteps.index+1, inplace = True)
                dfTempCSteps['Ns'].replace(newCSteps['Ns'].values, newCSteps.index+1, inplace = True)
                
                # Prints dataset after all transformations [Default Commented Out]
                #print(dfTempC)
                #print(dfTempC.loc[9500:13570])
                #print(dfTempCSteps[0:50], dfTempCSteps[50:100], dfTempCSteps[100:150])
                
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
            
            dfC.rename(columns={'mode':'Step Type', 
                                'time/s':'Run Time (h)', 
                                'I/mA':'Current (A)', 
                                'Ewe-Ece/V':'Potential vs. Counter (V)', 
                                'Ewe/V':'Potential (V)', 
                                '(Q-Qo)/mA.h':'Capacity (Ah)', 
                                'Ns':'Step Number'}, 
                       inplace=True)
                
            # Add header info into c file
            pathFileC = Path(path) / (name + ' Charge.csv')
            with open(pathFileC, 'w') as f:
                f.write(csvHeader)
            
            # Add data into c file
            dfC.to_csv(pathFileC, mode='a', index=False)
        
        # Combine file data
        if form_files:
            
            # Concatenate
            df = pd.concat([df, dfForm])
            
        if d_files:
            
            # Add last capacity, time, and step number to suceeding file
            if not(df.empty):
                dfD['Run Time (h)'] = dfD['Run Time (h)'] + df['Run Time (h)'].iat[-1]
                dfD['Capacity (Ah)'] = dfD['Capacity (Ah)'] + df['Capacity (Ah)'].iat[-1]
                dfD['Step Number'] = dfD['Step Number'] + df['Step Number'].iat[-1]
                    
            # Concatenate
            df = pd.concat([df, dfD])
            
        if c_files:
            
            # Add last capacity, time, and step number to suceeding file
            if not(df.empty):
                dfC['Run Time (h)'] = dfC['Run Time (h)'] + df['Run Time (h)'].iat[-1]
                dfC['Capacity (Ah)'] = dfC['Capacity (Ah)'] + df['Capacity (Ah)'].iat[-1]
                dfC['Step Number'] = dfC['Step Number'] + df['Step Number'].iat[-1]
                    
            # Concatenate
            df = pd.concat([df, dfC])
        
        # Add Header info into complete file
        pathFile = Path(path) / (name + '.csv')
        with open(pathFile, 'w') as f:
            f.write(csvHeader)
            
        # Add data into complete file
        df.to_csv(pathFile, mode='a', index=False)
        
        # Generate complete graph
        with plt.style.context('grapher'):
        
            fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(10, 4), gridspec_kw={'wspace':0.0})
            
            if form_files:
                axs[0].plot(dfForm['Run Time (h)'], dfForm['Potential vs. Counter (V)'], 'k--', label = 'Formation vs. Counter')
                axs[0].plot(dfForm['Run Time (h)'], dfForm['Potential (V)'], 'k-', label = 'Formation vs. Ref')
            
            if d_files:
                axs[0].plot(dfD['Run Time (h)'], dfD['Potential vs. Counter (V)'], 'r--', label = 'Discharge vs. Counter')
                axs[0].plot(dfD['Run Time (h)'], dfD['Potential (V)'], 'r-', label = 'Discharge vs. Ref')
            
            if c_files:
                axs[0].plot(dfC['Run Time (h)'], dfC['Potential vs. Counter (V)'], 'b--', label = 'Charge vs. Counter')
                axs[0].plot(dfC['Run Time (h)'], dfC['Potential (V)'], 'b-', label = 'Charge vs. Ref')

            axs[0].set_xlabel('Time (h)')
            axs[0].set_ylabel('Voltage (V)')
            
            if form_files:
                axs[1].plot(dfForm['Capacity (Ah)']*1000000/massVal, dfForm['Potential vs. Counter (V)'], 'k--', label = 'Formation vs. Counter')
                axs[1].plot(dfForm['Capacity (Ah)']*1000000/massVal, dfForm['Potential (V)'], 'k-', label = 'Formation vs. Ref')
            
            if d_files:
                axs[1].plot(dfD['Capacity (Ah)']*1000000/massVal, dfD['Potential vs. Counter (V)'], 'r--', label = 'Discharge vs. Counter')
                axs[1].plot(dfD['Capacity (Ah)']*1000000/massVal, dfD['Potential (V)'], 'r-', label = 'Discharge vs. Ref')
            
            if c_files:
                axs[1].plot(dfC['Capacity (Ah)']*1000000/massVal, dfC['Potential vs. Counter (V)'], 'b--', label = 'Charge vs. Counter')
                axs[1].plot(dfC['Capacity (Ah)']*1000000/massVal, dfC['Potential (V)'], 'b-', label = 'Charge vs. Ref')
            
            axs[1].set_xlabel('Capacity (mAh/g)')
            
            plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left')
            
            if export_fig is True:
                plt.savefig(Path(path) / 'complete_protocol_{}.jpg'.format(name))
            
            plt.show()
            
            if d_files and c_files:
                dfDOCV = dfD[dfD['Step Type'] != 0].drop_duplicates(['Step Number'])
                dfCOCV = dfC[dfC['Step Type'] != 0].drop_duplicates(['Step Number'])
                
                fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(10, 4), gridspec_kw={'wspace':0.0})
                axs.plot(dfDOCV['Capacity (Ah)']*1000000/massVal, dfDOCV['Potential vs. Counter (V)'], 'r.-', label = 'Discharge OCV vs. Ref')
                axs.plot(dfCOCV['Capacity (Ah)']*1000000/massVal, dfCOCV['Potential vs. Counter (V)'], 'b.-', label = 'Charge OCV vs. Ref')
                axs.set_xlabel('Capacity (mAh/g)')
                axs.set_ylabel('Voltage (V)')
                plt.legend(loc='lower right')
                
                if export_fig is True:
                    plt.savefig(Path(path) / 'OCV_match_{}.jpg'.format(name))
                
                plt.show()
        
class AMIDR():
    
    def __init__(self, dstpath, srcpath, uhpc_files, cell_label, single_pulse, bytesIO=None, export_data=True, use_input_cap=True, 
                 fcap_min=0.0, capacitance_corr=False, spliced=False, force2e=False):
        
        self.single_p = single_pulse
        self.capacitance_corr = capacitance_corr
        self.fcap_min = fcap_min
        self.cell_label = cell_label
        self.dst = Path(dstpath) / self.cell_label
        # If does not exist, create dir.
        if self.dst.is_dir() is False:
            self.dst.mkdir()
            print('Create directory: {}'.format(self.dst))
        self.src = Path(srcpath)
        
        ### TODO:
        ### Need to modify concatenation tool to accommodate bytesIO
        ### and new parsing method. 
        if type(uhpc_files) is list:
            self.uhpc_file = self.dst / "{}-concatenated.csv".format(self.cell_label)
            # concatenate uhpc files if more than 1 is passed.
            fnames = [self.src / f for f in uhpc_files]
            with open(fnames[0], 'r') as f:
                f.readline()
                cellname = f.readline().strip().split()[-1]
                f.readline()
                f.readline()
                mass = float(f.readline().strip().split()[-1]) / 1000
                capacity = float(f.readline().strip().split()[-1])
                for i in range(4):
                    f.readline()
                if f.readline().strip().split()[0] == '[Data]':
                    hlinenum = 11
                else:
                    hlinenum = 12
            with open(fnames[0], 'r') as f:
                header = f.readlines()[:hlinenum]
                    
            df1 = pd.read_csv(fnames[0], header=hlinenum)
            df2 = pd.read_csv(fnames[1], header=hlinenum)
            df2['Time (h)'] = df2['Time (h)'].values + df1['Time (h)'].values[-1]
            df2['Capacity (Ah)'] = df1['Capacity (Ah)'].values[-1] + df2['Capacity (Ah)'].values
            df2['Prot.Step'] = df2['Prot.Step'].values + df1['Prot.Step'].values[-1]
            df_concat = pd.concat([df1, df2], axis=0)
            with open(self.uhpc_file, 'w') as g:
                for i in range(len(header)):
                    g.write(header[i])
                df_concat.to_csv(path_or_buf=g, mode='a', index=False)
        else:
            self.uhpc_file = self.src / uhpc_files
            
        # need to update bytesIO read.
        if bytesIO is not None:
            headlines = []
            i = 0
            for line in bytesIO:
                headlines.append(line.decode().strip().split())
                if i == 12:
                    break
                i = i + 1
        else:  
            with open(self.uhpc_file, 'r') as f:
                lines = f.readlines()
            nlines = len(lines)
            headlines = []
            for i in range(nlines):
                headlines.append(lines[i])
                l = lines[i].strip().split()
                if l[0][:6] == '[Data]':
                    hline = lines[i+1]
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
        print('Input cell capacity: {} Ah'.format(self.input_cap))
        
        if bytesIO is None:
            self.df = pd.read_csv(self.uhpc_file, header=nskip)
            #self.df = pd.read_csv(self.uhpc_file, header=hlinenum)
        else:
            self.df = pd.read_csv(bytesIO, header=hlinenum)
            
        self.df.rename(columns={'Capacity (Ah)': 'Capacity',
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
                       inplace=True)
        #print(self.df.columns)
        #print(self.df.Step.unique())
        
        if single_pulse == True and spliced == True:
            sys.exit("single_pulse cannot operate on spliced files. Manually clean up your spliced file and select spliced = false")
        
        # Add Prot_step column if column does not yet exist or spliced file is used.
        if 'Prot_step' not in self.df.columns or spliced == True:
            s = self.df.Step
            self.df['Prot_step'] = s.ne(s.shift()).cumsum() - 1

        #if hline[-4:] == 'Flag':
        #    self.df = self.df.rename(columns={'Flag':'Prot.Step'})
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
            print('3-electrode data detected. Ignoring working potential and using complete potential for everything. [NOT RECCOMMENDED]')
        else:
            print('3-electrode data detected. Using working potential for calculations and complete potential for graphs and labels.')
        
        #plt.plot(self.df['Capacity'], self.df['Potential'])
        
        self.sigdf = self._find_sigcurves()
        #plt.plot(self.sigdf['Capacity'], self.sigdf['Potential'])
        self.sc_stepnums = self.sigdf['Prot_step'].unique()
        self.capacity = self.sigdf['Capacity'].max() - self.sigdf['Capacity'].min()
        self.spec_cap = self.capacity / self.mass
        print('Specific Capacity achieved in advanced protocol (signature curves): {0:.2f} mAh/g'.format(self.spec_cap*1000))
        if use_input_cap is True:
            self.capacity = self.input_cap
        print('Using {:.8f} Ah to compute rates.'.format(self.capacity))

        self.caps, self.scaps, self.fcaps, self.rates, self.eff_rates, self.currs, \
        self.ir, self.dqdv, self.resistdrop, self.soccaps, self.ivolts, self.cvolts, self.avg_volts, \
        self.dvolts, self.vlabels = self._parse_sigcurves()
        
        self.nvolts = len(self.caps)
                    
        if export_data is True:
            caprate_fname = self.dst / '{0}_rate-cap.xlsx'.format(self.cell_label)
            writer = pd.ExcelWriter(caprate_fname)
            for i in range(self.nvolts):
                caprate_df = pd.DataFrame(data={'specific_capacity': self.scaps[i],
                                                'fractional_capacity': self.fcaps[i],
                                                'effective_rate': self.eff_rates[i],
                                                'C-rates': self.rates[i]})
                caprate_df.to_excel(writer, sheet_name=self.vlabels[i], index=False)
            writer.save()
            writer.close()
            
        print('Done parsing signature curves.')
        


    def _find_sigcurves(self):
        """
        Use control "step" to find sequence of charge/discharge - OCV 
        characteristic of signature curves.
        """
        newdf = self.df.drop_duplicates(subset=['Step', 'Prot_step'])
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
    
    def plot_protocol(self, xlims=None, ylims=None, export_fig=True, return_fig=False):
        
        with plt.style.context('grapher'):
        
            fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True,
                                    figsize=(10, 4), gridspec_kw={'wspace':0.0})
            axs[0].plot(self.df['Time'], self.df['Label Potential'], 'k-')
            axs[0].set_xlabel('Time (h)')
            axs[0].set_ylabel('Voltage (V)')
            axs[1].set_xlabel('Specific Capacity (mAh/g)')
            #axs[0].tick_params(direction='in', top=True, right=True)
            
            # plot signature curves first if first
            if self.sc_stepnums[0] == 1:
                axs[1].plot(self.sigdf['Capacity']*1000/self.mass, self.sigdf['Label Potential'],
                    color='red',
                    label='Sig Curves')
            
            stepnums = self.df['Prot_step'].unique()
            #print(stepnums)
            fullsteps = np.setdiff1d(stepnums, self.sc_stepnums)
            #print(fullsteps)
            #print(self.sc_stepnums)
            # Need to set prop cycle
            colors = plt.get_cmap('viridis')(np.linspace(0,1,len(fullsteps)+1))
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
                            color=colors[c],
                            label=label)
                
                c = c + 1
                
                # if the next step is the start of sigcurves, plot sigcurves
                if fullsteps[i] == self.sc_stepnums[0] - 1:
                    #print('plotting sig curves...')
                    axs[1].plot(self.sigdf['Capacity']*1000/self.mass, self.sigdf['Label Potential'],
                            color='red',
                            label='Sig Curves')
            
            plt.legend(bbox_to_anchor=(1.0, 0.5), loc='center left')
            
            if xlims is not None:
                axs[1].set_xlim(xlims[0], xlims[1])
            if ylims is not None:
                axs[0].set_ylim(ylims[0], ylims[1])
            
            if export_fig is True:
                plt.savefig(self.dst / 'protocol_vis_{}.jpg'.format(self.cell_label))
            elif return_fig is True:
                return fig
            else:
                plt.show()
        
    
    def plot_caps(self, export_fig=True):
        
        with plt.style.context('grapher'):
            fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True,
                                    figsize=(10, 10), gridspec_kw={'hspace':0.0})
            colors = plt.get_cmap('viridis')(np.linspace(0,1,self.nvolts))
            
            for i in range(self.nvolts):
                axs[0].semilogx(self.eff_rates[i], self.scaps[i],
                                color=colors[i], label=self.vlabels[i])
                axs[1].semilogx(self.eff_rates[i], self.fcaps[i],
                                color=colors[i])
            
            if self.single_p:
                axs[1].set_xlabel('q$_{tot}$/I', fontsize=16)
            else:
                axs[1].set_xlabel('n$_{eff}$ in C/n$_{eff}$', fontsize=16)
            axs[1].set_ylabel('τ', fontsize=16)
            axs[0].set_ylabel('Specific Capacity (mAh/g)', fontsize=16)
            axs[0].legend(frameon=False, bbox_to_anchor=(1.0, 0.0), loc='center left', ncol = 1 + self.nvolts//30, fontsize=12)
            axs[0].tick_params(direction='in', which='both', top=True, right=True)
            axs[1].tick_params(direction='in', which='both', top=True, right=True)
            
            if export_fig is True:
                plt.savefig(self.dst / 'cap-rate_{}.jpg'.format(self.cell_label))
            else:
                plt.show()

    
    def _parse_sigcurves(self):

        sigs = self.sigdf.loc[self.sigdf['Step'] != 0]
        #capacity = sigs['Capacity'].max() - sigs['Capacity'].min()
        #print('Specific Capacity: {} mAh'.format(capacity))
        Vstart = np.around(sigs['Label Potential'].values[0], decimals=3)
        Vend = np.around(sigs['Label Potential'].values[-1], decimals=3)
        print('Starting voltage: {:.3f} V'.format(Vstart))
        print('Ending voltage: {:.3f} V'.format(Vend))
        
        sigsteps = sigs['Prot_step'].unique()
        nsig = len(sigsteps)
        print('Found {} charge or discharge steps in sig curve sequences.'.format(nsig))
        caps = []
        scaps = []
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
                stepcaps = step['Capacity'].values
                volts = step['Potential'].values
                currents = np.absolute(step['Current'].values)
                rate = self.capacity / np.average(currents)
                minarg = np.argmin(np.absolute(RATES - rate))
                
                # slice first and last current values if possible.
                # if less than 4(NVX) or 5(UHPC) data points, immediate voltage cutoff reached, omit step.
                if len(currents) > 3:
                    if volts[-2] == np.around(volts[-2], decimals=2):
                        currents = currents[1:-1]
                        cvoltind = -2
                    elif len(currents) > 4:
                        if volts[-3] == np.around(volts[-3], decimals=2):
                            currents = currents[1:-1]
                            cvoltind = -3
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
                    
                # determine dqdv based on the measurements before the voltage cutoff
                diffq = (stepcaps[cvoltind-2] - stepcaps[cvoltind-1]) / (volts[cvoltind-2] - volts[cvoltind-1])
                
                #if (np.amax(stepcaps) - np.amin(stepcaps))/self.mass < 5e-5:
                #    continue
            
                if caps == []:
                    caps.append([np.amax(stepcaps) - np.amin(stepcaps)])
                    rates.append([RATES[minarg]])
                    #initcutvolt = np.around(volts[0], decimals=3)
                    initcap.append([stepcaps[0]])
                    cutcap.append([stepcaps[cvoltind]])
                    initvolts.append([volts[0]])
                    cutvolts.append([volts[cvoltind]])
                    currs.append([np.average(currents)])
                    ir.append([np.absolute(volts[0] - volts[1])])
                    dqdv.append([diffq])
                    resistdrop.append([ir[-1][-1]/currs[-1][-1]])
                else:
                    #if np.amax(currents) < currs[-1][-1]:
                    if volts[cvoltind] == cutvolts[-1][-1]:
                        caps[-1].append(np.amax(stepcaps) - np.amin(stepcaps))
                        rates[-1].append(RATES[minarg])
                        cutcap[-1].append(stepcaps[cvoltind])
                        cutvolts[-1].append(volts[cvoltind])
                        currs[-1].append(np.average(currents))
                        ir[-1].append(np.absolute(volts[0] - volts[1]))
                        dqdv[-1].append(diffq)
                        resistdrop[-1].append(ir[-1][-1]/currs[-1][-1])
                    else:
                        if np.absolute(volts[-2] - cutvolts[-1][-1]) < 0.001:
                            continue
                        #print(np.average(currents), volts[-2])
                        caps.append([np.amax(stepcaps) - np.amin(stepcaps)])
                        rates.append([RATES[minarg]])
                        initcap.append([stepcaps[0]])
                        cutcap.append([stepcaps[cvoltind]])
                        initvolts.append([volts[0]])
                        cutvolts.append([volts[cvoltind]])
                        currs.append([np.average(currents)])
                        ir.append([np.absolute(volts[0] - volts[1])])
                        dqdv.append([diffq])
                        resistdrop.append([ir[-1][-1]/currs[-1][-1]])
            
            nvolts = len(caps)
            for i in range(nvolts):
                fcaps.append(np.cumsum(caps[i]) / np.sum(caps[i]))
                scaps.append(np.cumsum(caps[i]))
                
                eff_rates.append(scaps[i][-1]/currs[i])
                
                # Remove data where capacity is too small due to IR
                # i.e., voltage cutoff was reached immediately.
                inds = np.where(fcaps[i] < self.fcap_min)[0]
                if len(inds) > 0:
                    caps[i] = np.delete(caps[i], inds)
                    scaps[i] = np.delete(scaps[i], inds)
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
            # icaps is the idealized capacity for a given voltage based upon dqdv
            # currsCum is the cumulative averge current for determining where C should be calculated
            icaps = []
            currsCum = []
            voltsAct = []
            time = []
            for i in range(nsig):
                step = sigs.loc[sigs['Prot_step'] == sigsteps[i]]
                stepcaps = step['Capacity'].values
                volts = step['Potential'].values
                lvolts = step['Label Potential'].values
                currents = np.absolute(step['Current'].values)
                runtime = step['Time'].values
                
                # Collect succeeding OCV steps (1 OCV or 2 OCV) to calculate dqdv
                ocvstep = self.sigdf.loc[self.sigdf['Prot_step'] == sigsteps[i] + 2]
                if ocvstep['Step'].values[0] != 0:
                    ocvstep = self.sigdf.loc[self.sigdf['Prot_step'] == sigsteps[i] + 1]
                ocvstepcaps = ocvstep['Capacity'].values
                ocvvolts = ocvstep['Potential'].values
                ocvlvolts = ocvstep['Label Potential'].values
                
                ir.append([np.absolute(volts[1] - volts[0])])
                
                dqdv.append([(stepcaps[0] - ocvstepcaps[-1])/(volts[0] - ocvvolts[-1])])
                
                time.append(np.absolute(runtime[1:] - runtime[0]))
                caps.append(np.absolute(stepcaps[1:] - stepcaps[0]))
                scaps.append(np.absolute(stepcaps[1:] - stepcaps[0]))
                icaps.append(np.absolute(dqdv[-1][0]*(volts[1:] - volts[0])))
                currs.append(currents[1:])
                voltsAct.append(volts[1:])
                
                currsCum.append([])
                for j in range(len(currents[1:])):
                    currsCum[-1].append(np.average(currents[1:j+2]))
                    minarg = np.argmin(np.absolute(RATES - self.capacity / currsCum[-1][j]))
                    if j == 0:
                        rates.append([RATES[minarg]])
                    else:
                        rates[-1].append(RATES[minarg])
                
                resistdrop.append([ir[-1][-1]/currs[-1][0]])
                
                #if cutvolts == []:
                    #initcutvolt = np.around(lvolts[0], decimals=3)
                initcap.append([stepcaps[0]])
                cutcap.append([ocvstepcaps[-1]]) 
                initvolts.append([lvolts[0]])
                cutvolts.append([ocvlvolts[-1]]) 
                
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
                    
                    icaps[i] = icaps[i] - dlcaps #- dqdv[i][0]*currs[i]*rohm 
                    # if icaps is calculated as negative or zero, this datapoint is effectively thrown out (caps set to nan)
                    for j in range(len(caps[i])):
                        if icaps[i][j] <= 0:
                            icaps[i][j] = float('NaN')
                    
                    #currsCum[i] = currsCum[i] - dlcaps/time[i] # disabled as it may amplify error if near 0
                    # if cumulative current is calculated as negative or zero, this datapoint is effectively thrown out (caps set to nan)
                    for j in range(len(caps[i])):
                        if currsCum[i][j] <= 0:
                            currsCum[i][j] = float('NaN')
            
            for i in range(nvolts):
                fcaps.append(caps[i]/icaps[i])
                eff_rates.append(icaps[i]/currsCum[i])
                
                # outlier repair: if fcap or eff_rates is NaN, make it equal to the succeeding point (or previous if last point).
                for j in range(len(caps[i])):
                    if fcaps[i][-j - 1] != fcaps[i][-j - 1] or eff_rates[i][-j - 1] != eff_rates[i][-j - 1]:
                        fcaps[i][-j - 1] = fcaps[i][-j]
                        eff_rates[i][-j - 1] = eff_rates[i][-j]
        
        print('Found {} signature curves.'.format(nvolts))
        
        ivolts = np.zeros(nvolts)
        cvolts = np.zeros(nvolts)
        icap = np.zeros(nvolts)
        ccap = np.zeros(nvolts)
        for i in range(nvolts):
            ivolts[i] = np.average(initvolts[i])
            cvolts[i] = np.average(cutvolts[i])
            icap[i] = np.average(initcap[i])
            ccap[i] = cutcap[i][-1]
        
        with np.printoptions(precision=3):
            soccaps = (icap + ccap)/2
            print('Midpoint capacities: {}'.format((soccaps*1000/self.mass)))
            print('Cutoff voltages: {}'.format(cvolts))
            # Get midpoint voltage for each range.
            #avg_volt[0] = (initcutvolt + cvolts[0])/2
            #avg_volt[1:] = (cvolts[:-1] + cvolts[1:])/2
            avg_volt = (ivolts + cvolts)/2
            print('Midpoint voltages: {}'.format(avg_volt))
            dvolts = np.zeros(nvolts)
            #dvolts[0] = np.absolute(initcutvolt - cvolts[0])
            #dvolts[1:] = np.absolute(cvolts[:-1] - cvolts[1:])
            dvolts = np.absolute(ivolts - cvolts)
            with np.printoptions(precision=4):
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
                print("{} removed due to not having 4 or more datapoints with fcap above 0.001".format(vlabels[i-iadj]))
                caps.pop(i-iadj)
                scaps.pop(i-iadj)
                fcaps.pop(i-iadj)
                eff_rates.pop(i-iadj)
                rates.pop(i-iadj)
                currs.pop(i-iadj)
                ir.pop(i-iadj)
                dqdv.pop(i-iadj)
                resistdrop.pop(i-iadj)
                soccaps = np.delete(soccaps, i-iadj)
                ivolts = np.delete(ivolts, i-iadj)
                avg_volt = np.delete(avg_volt, i-iadj)
                dvolts = np.delete(dvolts, i-iadj)
                vlabels = np.delete(vlabels, i-iadj)
                iadj = iadj + 1
        nvolts = len(caps)
        
        new_caps = []
        new_scaps = []
        for i in range(nvolts):
            new_caps.append(1000*np.array(caps[i])/self.mass)
            new_scaps.append(1000*np.array(scaps[i])/self.mass)
        #print(new_caps)

        return new_caps, new_scaps, fcaps, rates, eff_rates, currs, ir, dqdv, resistdrop, soccaps, ivolts, cvolts, avg_volt, dvolts, vlabels 
       

    def fit_atlung(self, r, r_corr, tracer_inputs = [], ftol=5e-14, D_bounds=[1e-17, 1e-8], D_guess=1.0e-11, 
                   fcapadj_bounds=[1.0, 1.5], fcapadj_guess=1.0, R_eff_bounds=[1e-6, 1e1], R_eff_guess=1.0e-2, 
                   shape='sphere', nalpha=4000, nQ=4000, export_data=True, export_fig=True, label=None):

        self.r = r
        self.r_corr = r_corr
        
        if shape not in SHAPES:
            print('The specified shape {0} is not supported.'.format(shape))
            print('Supported shapes are: {1}. Defaulting to sphere.'.format(SHAPES))
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
        if self.r_corr is False:
            print("Optimum Parameters: {}".format("Log(Dc) fCapAdj"))
            Q_arr = np.logspace(-3, 2, nQ)
            tau_sol = np.zeros(nQ)
            tau_guess = 0.5
            for i in range(nQ):
                Q = Q_arr[i]
                func = lambda tau: tau - 1 + (1/(A*Q))*(1/B - 2*(np.sum(np.exp(-self.alphas*tau*Q)/self.alphas)))
                tau_sol[i] = fsolve(func, tau_guess, factor=1.)
        elif self.single_p is False:
            print("Optimum Parameters: {}".format("Log(Dc) fCapAdj Log(Reff) Log(Reff/Dc)"))
        else:
            print("Optimum Parameters: {}".format("Log(Dc) Log(Reff) Log(Reff/Dc)"))
                
        dconst = np.zeros(self.nvolts, dtype=float)
        resist_eff = np.zeros(self.nvolts, dtype=float)
        resist = np.zeros(self.nvolts, dtype=float)
        dqdv = np.zeros(self.nvolts, dtype=float)
        sigma = np.zeros(self.nvolts, dtype=float)
        fit_err = np.zeros(self.nvolts, dtype=float)
        cap_max = np.zeros(self.nvolts, dtype=float)
        cap_min = np.zeros(self.nvolts, dtype=float)
        cap_span = np.zeros(self.nvolts, dtype=float)

        for j in range(self.nvolts):
            z = np.ones(len(self.fcaps[j]))
            #fcap = np.array(self.fcaps[j])
            fcap = np.array(self.fcaps[j])
            #self._max_cap = self.scaps[j][-1]
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
                fcapadj_bounds=[1.0, 1.0000001]

            #print("dq/dV: {} Ah/V".format(dqdv[j]))
            
            if self.r_corr is False:
                C = np.sum(self.ir[j])
                weights = (C - self.ir[j]) / np.sum(C - self.ir[j])
                bounds = ([np.log10(D_bounds[0]), fcapadj_bounds[0]],
                          [np.log10(D_bounds[1]), fcapadj_bounds[1]])
                p0 = [np.log10(D_guess), fcapadj_guess]    
            else:
                bounds = ([np.log10(D_bounds[0]), fcapadj_bounds[0], np.log10(R_eff_bounds[0])],
                          [np.log10(D_bounds[1]), fcapadj_bounds[1], np.log10(R_eff_bounds[1])])
                p0 = [np.log10(D_guess), fcapadj_guess, np.log10(R_eff_guess)]
                
            with plt.style.context('grapher'):
                fig = plt.figure()
                
                if shape == 'sphere':
                    if self.r_corr is False:
                        popt, pcov = curve_fit(self._spheres, (fcap, rates), z, p0=p0,
                                   bounds=bounds, sigma=weights,
                                   method='trf', max_nfev=5000, x_scale=[1.0, 1.0],
                                   ftol=ftol, xtol=None, gtol=None, loss='soft_l1', f_scale=1.0)
                        with np.printoptions(precision=4):
                            print("{}: {}".format(self.vlabels[j], popt))
                    else:
                        p0opt = [p0[2] - p0[0], p0[1], p0[2]]
                        boundsopt = [[bounds[0][2] - bounds[1][0], bounds[0][1], bounds[0][2]], 
                                      [bounds[1][2] - bounds[0][0], bounds[1][1], bounds[1][2]]]
                        
                        popt, pcov = curve_fit(self._spheres_r_corr, (fcap, rates), z, p0=p0opt,
                                   bounds=boundsopt,
                                   method='trf', max_nfev=5000, x_scale=[1.0, 1.0, 1.0],
                                   ftol=ftol, xtol=None, gtol=None, loss='soft_l1', f_scale=1.0)
                        popt = np.array([popt[2] - popt[0], popt[1], popt[2], popt[0]])
                        with np.printoptions(precision=3):
                            if self.single_p is False:
                                print("{}: {}".format(self.vlabels[j], popt))
                            else:
                                print("{}: {}".format(self.vlabels[j], np.array([popt[2] - popt[0], popt[2], popt[0]])))
                        resist_eff[j] = 10**popt[2]
                        Q_arr = np.logspace(-6, 2, nQ)
                        tau_sol = np.zeros(nQ)
                        tau_guess = 0.5
                        for i in range(nQ):
                            Q = Q_arr[i]
                            func = lambda tau: tau - 1 + (1/(A*Q))*(1/B - 2*(np.sum(np.exp(-self.alphas*tau*Q)/self.alphas))) + 10**popt[2]/Q if 10**popt[2]<Q else tau

                            tau_sol[i] = fsolve(func, tau_guess, factor=1.)
                            if tau_sol[i] < 0:
                                tau_sol[i] = 0
                        
                if shape == 'plane':
                    popt, pcov = curve_fit(self._planes, (fcap, rates), z, p0=p0,
                               bounds=bounds, sigma=weights,
                               method='trf', max_nfev=5000, x_scale=[1e-11, 1.0],
                               ftol=ftol, xtol=None, gtol=None, loss='soft_l1', f_scale=1.0)
                
                #plt.semilogx(Q_arr, tau_sol, '-k', label='Atlung - {}'.format(shape))
                
                sigma[j] = np.sqrt(np.diag(pcov))[0]
                dconst[j] = 10**popt[0]
                Qfit = 3600*rates*dconst[j]/r**2
                tau_fit = fcap/popt[1]
                
                cap_max[j] = tau_fit[-1]
                cap_min[j] = tau_fit[0]
                cap_span[j] = tau_fit[-1] - tau_fit[0]
                
                # get difference between fitted values and
                # theoretical Atlung curve to get fit_err.
                error = np.zeros(len(Qfit), dtype=float)
                for k in range(len(Qfit)):
                    dQ = np.absolute(Q_arr - Qfit[k])
                    minarg = np.argmin(dQ)
                    error[k] = np.absolute(tau_fit[k] - tau_sol[minarg])
                if r_corr is False:
                    fit_err[j] = np.sum(weights*error)
                else:
                    fit_err[j] = np.sqrt(np.average((error/cap_max[j])**2))
                
                plt.semilogx(Qfit, tau_fit, 'or', markersize=4, label='{0} - {1}'.format(self.cell_label, self.vlabels[j]))
                if max(tau_fit) < 0.01:
                    plt.ylim(0, 0.01)
                elif max(tau_fit) < 0.1:
                    plt.ylim(0, 0.1)
                else:
                    plt.ylim(0, 1)
                plt.semilogx(Q_arr, tau_sol, '-k', label='Atlung - {}'.format(shape))
                plt.xlabel('$Q$')
                plt.ylabel('$τ$')
                plt.legend(frameon=False, loc='upper left', fontsize=10)
                if Qfit[0] < 1.0e-4 or Qfit[1] < 1.0e-3:
                    plt.xlim(1.0e-6, 1.0e0)
                    plt.xticks(10.**np.arange(-6, 1))
                else:
                    plt.xlim(1.0e-4, 1.0e2)
                    plt.xticks(10.**np.arange(-4, 3))
                if export_fig is True:
                    if label is None:
                        figname = self.dst / '{0}_Atlung-{1}_{2:.3f}.jpg'.format(self.cell_label, shape, self.avg_volts[j])
                    else:
                        figname = self.dst / '{0}-{1}_Atlung-{2}_{3:.3f}.jpg'.format(self.cell_label, label, shape, self.avg_volts[j])
                    plt.savefig(figname)
                else:
                    plt.show()
                plt.close()
        
        if r_corr is False:
            DV_df = pd.DataFrame(data={'Voltage': self.avg_volts, 'D': dconst})
        else:
            # get resist from resist_eff
            resist = resist_eff*self.r**2/(3600*dconst*dqdv)
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
            
            DV_df = pd.DataFrame(data={'Voltage (V)': self.avg_volts, 'Initial Voltage (V)': self.ivolts, 'Dc (cm^2/s)': dconst, 'R_eff' : resist_eff,
                                       'dqdV (mAh/gV)': [i*1000/self.mass for i in dqdv], 'Rfit (Ohm)' : resist, 'Rdrop (Ohm)' : rdrop, 'Cap Span' : cap_span, 'Fit Error' : fit_err})
            
            # Calculates Free-path Tracer D from inputs
            if tracer_inputs:
                theorcap = tracer_inputs[0]/1000*self.mass
                initialcap = tracer_inputs[1]/1000*self.mass
                temp = tracer_inputs[2]
                kB = 1.380649E-23
                e = 1.602176634E-19
                pulsecap = initialcap + self.soccaps
                soc = pulsecap/theorcap
                dtconst = dconst*kB*temp*dqdv/(e*pulsecap*(1. - pulsecap/theorcap))
                
                DV_df['Dt* (cm^2/s)'] = dtconst
                DV_df['SOC'] = soc
                DV_df = DV_df[['Voltage (V)', 'Initial Voltage (V)', 'Dc (cm^2/s)', 'Dt* (cm^2/s)', 'R_eff', 'dqdV (mAh/gV)', 'Rfit (Ohm)', 'Rdrop (Ohm)', 'Cap Span', 'Fit Error', 'SOC']]
            else:
                dtconst = []
        
        if export_data is True:
            if label is None:
                df_filename = self.dst / '{0}_D-V_{1}.xlsx'.format(self.cell_label, shape)
            else:
                df_filename = self.dst / '{0}-{1}_D-V_{2}.xlsx'.format(self.cell_label, label, shape)
                
            DV_df.to_excel(df_filename, index=False)
        
        with np.printoptions(precision=3):
            print('Fitted Dc: {}'.format(dconst))
            if self.single_p is False:
                print('Standard deviations from fit: {}'.format(sigma))
                print('Atlung fit error: {}'.format(fit_err))
        
        return self.avg_volts, self.ivolts, dconst, dtconst, fit_err, cap_span, cap_max, cap_min, self.caps, self.ir, self.dvolts, resist_eff, dqdv, resist, self.resistdrop, self.single_p, self.r_corr, self.cell_label, self.mass, self. dst
        
    def make_summary_graph(self, fit_data, export_fig=True, label=None):
        
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
        resist_eff = fit_data[11]
        dqdv = fit_data[12]
        resist = fit_data[13]
        resistdrop = fit_data[14]
        single_p = fit_data[15]
        r_corr = fit_data[16]
        cell_label = fit_data[17]
        mass = fit_data[18]
        dst = fit_data[19]
        
        with plt.style.context('grapher'):
            if single_p is False:
                if r_corr is False:
                    fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(7, 15), sharex=True,
                    gridspec_kw={'height_ratios': [2,1,1,1,1], 'hspace': 0.0})
                    
                    axs[0].semilogy(voltage, dconst, 'ko--', linewidth=0.75, label='{}'.format(cell_label))
                    #axs[0].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    #axs[0].xaxis.set_minor_locator(MultipleLocator(0.1))
                    axs[0].get_xaxis().set_ticks(voltage)
                    axs[0].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[0].set_xlabel('Voltage (V)', fontsize=12)
                    axs[0].set_ylabel('D (cm$^2$/s)', fontsize=12)
                    axs[0].legend(frameon=False, fontsize=12)
                    
                    axs[1].plot(voltage, fit_err, 'ks--', linewidth=0.75)
                    axs[1].set_ylim(0, np.amin([0.5, np.amax(fit_err)]))
                    #axs[1].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    #axs[1].xaxis.set_minor_locator(MultipleLocator(0.1))
                    axs[1].get_xaxis().set_ticks(voltage)
                    axs[1].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[1].yaxis.set_minor_locator(MultipleLocator(0.05))
                    #axs[1].set_ylabel('Average \n fractional \n fit error', fontsize=12)
                    axs[1].set_ylabel('Weighted average \n absolute \n fit error', fontsize=12)
                    
                    axs[2].plot(voltage, cap_span, 'k^--', linewidth=0.75)
                    axs[2].set_ylim(0, 1.0)
                    #axs[2].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    #axs[2].xaxis.set_minor_locator(MultipleLocator(0.1))
                    axs[2].get_xaxis().set_ticks(voltage)
                    axs[2].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[2].get_yaxis().set_ticks([0, 0.25, 0.5, 0.75])
                    axs[2].set_xlabel('Voltage (V)', fontsize=12)
                    axs[2].set_ylabel('Fitted fractional \n capacity span', fontsize=12)
                    for j in range(nvolts):
                        axs[2].fill(np.array([voltage[j]-0.01, voltage[j]-0.01, voltage[j]+0.01, voltage[j]+0.01]),
                                    np.array([cap_min[j], cap_max[j], cap_max[j], cap_min[j]]),
                                    color='k', alpha=0.3, edgecolor='k', linestyle='-')
                    
                    for j in range(nvolts):
                        cap_in_step = caps[j]
                        ir = dV_ir[j]
                        axs[3].bar(voltage[j], np.sum(cap_in_step), width=dvolts[j], color='k', alpha=0.15,
                                   edgecolor='k')
                        for i in range(len(ir)):
                            #width = 0.1/len(cap_in_step)
                            #center = voltage[j] - 0.05 + (i+1/2)*width
                            width = dvolts[j]/len(cap_in_step)
                            center = voltage[j] - dvolts[j]/2 + (i+1/2)*width
                            axs[3].bar(center, cap_in_step[i], width=width, color='k', alpha=0.3)
                            axs[4].bar(voltage[j], ir[i]/dvolts[j], width=0.04, color='k', alpha=0.3)
                    axs[3].set_ylabel('Specific Capacity \n in step (mAh/g)', fontsize=12)
                    axs[3].tick_params(axis='x', which='minor', top=False, bottom=False)
                    #axs[3].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    #axs[4].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    axs[4].get_xaxis().set_ticks(voltage)
                    axs[4].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[4].set_xticklabels(['{:.3f}'.format(v) for v in voltage], rotation=45)
                    axs[4].set_xlabel('Voltage (V)', fontsize=12)
                    axs[4].set_ylabel('Fractional \n IR drop', fontsize=12)
                    
                else:
                    fig, axs = plt.subplots(ncols=1, nrows=6, figsize=(7, 15), sharex=True,
                    gridspec_kw={'height_ratios': [2,1,1,1,1,1], 'hspace': 0.0})
                    fig.suptitle('{}'.format(cell_label))
                
                    axs[0].semilogy(voltage, dconst, 'k+-', linewidth=0.75, label='Chemical D')
                    axs[0].semilogy(voltage, dtconst, 'kx-', linewidth=0.75, label='Free-path Tracer D')
                    #axs[0].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    #axs[0].xaxis.set_minor_locator(MultipleLocator(0.1))
                    axs[0].get_xaxis().set_ticks(voltage)
                    axs[0].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[0].set_xlabel('Voltage (V)', fontsize=12)
                    axs[0].set_ylabel('D (cm$^2$/s)', fontsize=12)
                    axs[0].legend(frameon=False, fontsize=12)
                    
                    axs[1].semilogy(voltage, resist_eff, 'ko-', linewidth=0.75)
                    axs[1].semilogy(voltage, 1.0/15*np.ones(len(voltage)), 'k:', linewidth=1.5)
                    axs[1].set_ylim(1.0e-5, 1.0)
                    axs[1].set_xlabel('Voltage (V)', fontsize=12)
                    axs[1].set_ylabel('R$_{eff}$', fontsize=12)
                    
                    axs[2].semilogy(voltage, resist, 'ko-', linewidth=0.75, label='fit R')
                    rdavg = []
                    rddev = []
                    for i in range(nvolts):
                        rdavg.append(np.average(resistdrop[i]))
                        rddev.append(np.std(resistdrop[i]))
                    axs[2].errorbar(voltage, rdavg, rddev, fmt='k^-', capsize = 3.0, linewidth=0.75, label='Vdrop R')
                    axs[2].xscale = 'log'
                    axs[2].set_ylim(1.0e1, 0.99e5)
                    axs[2].set_xlabel('Voltage (V)', fontsize=12)
                    axs[2].set_ylabel('R ($\Omega$)', fontsize=12)
                    axs[2].legend(frameon=False, fontsize=12)
                    
                    axs[3].plot(voltage, fit_err, 'ko-', linewidth=0.75)
                    axs[3].set_ylim(0, np.amin([0.5, np.amax(fit_err)]))
                    #axs[1].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    #axs[1].xaxis.set_minor_locator(MultipleLocator(0.1))
                    axs[3].get_xaxis().set_ticks(voltage)
                    axs[3].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[3].yaxis.set_minor_locator(MultipleLocator(0.05))
                    #axs[1].set_ylabel('Average \n fractional \n fit error', fontsize=12)
                    axs[3].set_ylabel('Fit error', fontsize=12)
                    
                    axs[4].plot(voltage, cap_span, 'ko-', linewidth=0.75)
                    axs[4].set_ylim(0, 1.0)
                    #axs[2].tick_params(direction='in', which='both', top=True, right=True, labelsize=12)
                    #axs[2].xaxis.set_minor_locator(MultipleLocator(0.1))
                    axs[4].get_xaxis().set_ticks(voltage)
                    axs[4].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[4].get_yaxis().set_ticks([0, 0.25, 0.5, 0.75])
                    axs[4].set_xlabel('Voltage (V)', fontsize=12)
                    axs[4].set_ylabel('Fitted fractional \n capacity span', fontsize=12)
                    for j in range(nvolts):
                        axs[4].fill(np.array([voltage[j]-0.01, voltage[j]-0.01, voltage[j]+0.01, voltage[j]+0.01]),
                                    np.array([cap_min[j], cap_max[j], cap_max[j], cap_min[j]]),
                                    color='k', alpha=0.3, edgecolor='k', linestyle='-')
                    
                    for j in range(nvolts):
                        cap_in_step = caps[j]
                        ir = dV_ir[j]
                        axs[5].bar(voltage[j], np.sum(cap_in_step), width=dvolts[j], color='k', alpha=0.15,
                                   edgecolor='k')
                        for i in range(len(ir)):
                            #width = 0.1/len(cap_in_step)
                            #center = voltage[j] - 0.05 + (i+1/2)*width
                            width = dvolts[j]/len(cap_in_step)
                            center = voltage[j] - dvolts[j]/2 + (i+1/2)*width
                            axs[5].bar(center, cap_in_step[i], width=width, color='k', alpha=0.3)
                            #axs[4].bar(voltage[j], ir[i]/dvolts[j], width=0.04, color='k', alpha=0.3)
                    axs[5].get_xaxis().set_ticks(voltage)
                    axs[5].tick_params(axis='x', which='minor', top=False, bottom=False)
                    axs[5].set_xticklabels(['{:.3f}'.format(v) for v in voltage], rotation=45)
                    axs[5].set_xlabel('Voltage (V)', fontsize=12)
                    axs[5].set_ylabel('Specific Capacity \n in step (mAh/g)', fontsize=12)
                    
            else:
                fig, axs = plt.subplots(ncols=1, nrows=6, figsize=(6, 12), sharex=True,
                gridspec_kw={'height_ratios': [1,1,1,1,1,1], 'hspace': 0.0})
                fig.suptitle('{}'.format(cell_label))
                
                axs[0].semilogy(voltage, dconst, 'k+-', linewidth=0.75, label='Chemical D')
                axs[0].semilogy(voltage, dtconst, 'kx-', linewidth=0.75, label='Free-path Tracer D')
                axs[0].set_xlabel('Voltage (V)', fontsize=12)
                axs[0].set_ylabel('D (cm$^2$/s)', fontsize=12)
                axs[0].legend(frameon=False, fontsize=12)
                
                axs[1].semilogy(voltage, resist_eff, 'k+-', linewidth=0.75)
                axs[1].semilogy(voltage, 1.0/15*np.ones(len(voltage)), 'k:', linewidth=1.5)
                axs[1].set_ylim(1.0e-6, 1.0)
                axs[1].set_xlabel('Voltage (V)', fontsize=12)
                axs[1].set_ylabel('$R_{eff}$', fontsize=12)
                
                axs[2].semilogy(ivoltage, resist, 'k+-', linewidth=0.75, label='fit R')
                axs[2].semilogy(ivoltage, resistdrop, 'kx-', linewidth=0.75, label='Vdrop R')
                axs[2].set_ylim(bottom=1.0e1)
                axs[2].set_xlabel('Voltage (V)', fontsize=12)
                axs[2].set_ylabel('R ($\Omega$)', fontsize=12)
                axs[2].legend(frameon=False, fontsize=12)
                
                axs[3].semilogy(voltage, fit_err, 'k+-', linewidth=0.75)
                #axs[2].set_ylim()
                axs[3].set_xlabel('Voltage (V)', fontsize=12)
                axs[3].set_ylabel('Fit Error', fontsize=12)
                
                axs[4].plot(voltage, cap_span, 'k+-', linewidth=0.75)
                axs[4].set_ylim(0, 1.0)
                axs[4].get_yaxis().set_ticks([0, 0.25, 0.5, 0.75])
                axs[4].set_xlabel('Voltage (V)', fontsize=12)
                axs[4].set_ylabel('Fitted fractional \n capacity span', fontsize=12)
                for j in range(nvolts):
                    axs[4].fill(np.array([voltage[j]-dvolts[j]/2, voltage[j]-dvolts[j]/2, voltage[j]+dvolts[j]/2, voltage[j]+dvolts[j]/2]),
                                np.array([cap_min[j], cap_max[j], cap_max[j], cap_min[j]]),
                                color='k', alpha=0.3, edgecolor='k', linestyle='-')
                
                axs[5].plot(voltage, [i*1000/mass for i in dqdv], 'k+-', linewidth=0.75)
                axs[5].set_xlabel('Voltage (V)', fontsize=12)
                axs[5].set_ylabel('dq/dV (mAh/gV)', fontsize=12)
        
            if export_fig is True:
                if label is None:
                    figstr = 'D-V_{0}.jpg'.format(cell_label)
                else:
                    figstr = 'D-V_{0}_{1}.jpg'.format(cell_label, label)
                diff_figname = dst / figstr
                plt.savefig(diff_figname)
                
            else:
                plt.show()
                

    def _spheres(self, X, logD, c_max):
        
        D = 10**logD
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis=1)
        narr = np.repeat(n.reshape(len(n), 1), len(self.alphas), axis=1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis=0)
        
        return c/c_max + ((self.r**2)/(3*3600*n*D))*(1/5 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis=1)))
    
    def _spheres_r_corr(self, X, logR_effDivD, c_max, logR_eff):
        
        D = 10**(logR_eff - logR_effDivD)
        R_eff = 10**logR_eff
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis=1)
        narr = np.repeat(n.reshape(len(n), 1), len(self.alphas), axis=1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis=0)
        
        # Calculates inacessible capacity as 1 + tau if R_eff/Q > 1 AND fcap is less than 0.05 of the largest fcap 
        # by setting n so that R_eff=Q. Otherwise standard AMIDR equation.
        # This avoids the divergent region where tau=0 but infinite summation error is amplified. 
        result = []
        for i in range(len(c)):
            if R_eff>(3600*n[i]*D)/self.r**2 and c[i]/max(c) < 0.05:
                result.append(c[i]/c_max + 1)
            else:
                result.append(c[i]/c_max + ((self.r**2)/(3*3600*n[i]*D))*(1/5 - 2*(np.sum(np.exp(-a[i]*(carr[i]/c_max)*3600*narr[i]*D/self.r**2)/a[i]))) + R_eff*self.r**2/(3600*n[i]*D))
        
        #return c/c_max + ((self.r**2)/(3*3600*n*D))*(1/5 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis=1))) + self._dqdv*I*R_eff/self._max_cap
        return result
    
    def _planes(self, X, logD, c_max):
        
        D = 10**logD
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis=1)
        narr = np.repeat(n.reshape(len(c), 1), len(self.alphas), axis=1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis=0)
        
        return c/c_max + ((self.r**2)/(3600*n*D))*(1/3 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis=1)))
    
    def insert_rate_cap(self, rate_cap):

        new_rate_cap = pd.read_csv(rate_cap, na_values = ['no info', '.'])
        
        self.nvolts = 1
        self.fcaps = [np.array(new_rate_cap['Capacity'])]
        self.eff_rates = [list(new_rate_cap['n in C/n'].values)]
        self.ir = [list(new_rate_cap['Crate'].values)]
        self.vlabels = ['inserted']
        self.avg_volts = [0]
        self.dqdv = [list(new_rate_cap['dqdv'].values)] 