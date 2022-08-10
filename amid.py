#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 07:44:54 2021

@author: Marc M. E. Cormier
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

class AMID():
    
    def __init__(self, dstpath, srcpath, uhpc_files, cell_label, bytesIO=None, export_data=True, 
                 use_input_cap=False, fcap_min=0.0, single_current=False, capacitance_corr=False, spliced=False):
        
        self.single_curr = single_current
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
        self.cellname = m[-1]
        
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
        
        if single_current == True and spliced == True:
            sys.exit("single_current cannot operate on spliced files. Manually clean up your spliced file and select spliced = false")
        
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

        self.caps, self.scaps, self.fcaps, self.rates, self.eff_rates, \
        self.currs, self.ir, self.dqdv, self.resistdrop, self.cvolts, self.avg_volts, \
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
        
        if self.single_curr is False:
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
            #single_current sigcurves selection
            for i in range(len(ocv_inds)):
                if i+1 == len(ocv_inds):
                    print("No adjacent OCV steps detected. Protocol is likely not single_current")
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
            axs[0].plot(self.df['Time'], self.df['Potential'], 'k-')
            axs[0].set_xlabel('Time (h)')
            axs[0].set_ylabel('Voltage (V)')
            axs[1].set_xlabel('Specific Capacity (mAh/g)')
            #axs[0].tick_params(direction='in', top=True, right=True)
            
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
                
                axs[1].plot(stepdf['Capacity']*1000/self.mass, stepdf['Potential'],
                            color=colors[c],
                            label=label)
                c = c + 1
                
                # if the next step is the start of sigcurves, plot sigcurves
                if fullsteps[i] == self.sc_stepnums[0] - 1:
                    #print('plottting sig curves...')
                    axs[1].plot(self.sigdf['Capacity']*1000/self.mass, self.sigdf['Potential'],
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
                                    figsize=(6, 7), gridspec_kw={'hspace':0.0})
            colors = plt.get_cmap('viridis')(np.linspace(0,1,self.nvolts))
            
            for i in range(self.nvolts):
                axs[0].semilogx(self.eff_rates[i], self.scaps[i],
                                color=colors[i], label=self.vlabels[i])
                axs[1].semilogx(self.eff_rates[i], self.fcaps[i],
                                color=colors[i])
                
            axs[1].set_xlabel('n$_{eff}$ in C/n$_{eff}$', fontsize=16)
            axs[1].set_ylabel('Fractional Capacity', fontsize=16)
            axs[0].set_ylabel('Specific Capacity (mAh/g)', fontsize=16)
            axs[0].legend(frameon=False, bbox_to_anchor=(1.0, 0.0), loc='center left')
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
        Vstart = np.around(sigs['Potential'].values[0], decimals=3)
        Vend = np.around(sigs['Potential'].values[-1], decimals=3)
        print('Starting voltage: {:.3f} V'.format(Vstart))
        print('Ending voltage: {:.3f} V'.format(Vend))
        
        sigsteps = sigs['Prot_step'].unique()
        nsig = len(sigsteps)
        print('Found {} charge or discharge steps in sig curve sequences.'.format(nsig))
        caps = []
        scaps = []
        fcaps = []
        rates = []
        cutvolts = []
        currs = []
        ir = []
        dqdv = []
        resistdrop = []
        eff_rates = []
        
        if self.single_curr is False:
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
                    initcutvolt = np.around(volts[0], decimals=3)
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
                        cutvolts[-1].append(volts[cvoltind])
                        currs[-1].append(np.average(currents))
                        #currs[-1].append(np.amax(currents[1:]))
                        ir[-1].append(np.absolute(volts[0] - volts[1]))
                        dqdv[-1].append(diffq)
                        resistdrop[-1].append(ir[-1][-1]/currs[-1][-1])
                    else:
                        if np.absolute(volts[-2] - cutvolts[-1][-1]) < 0.001:
                            continue
                        #print(np.average(currents), volts[-2])
                        caps.append([np.amax(stepcaps) - np.amin(stepcaps)])
                        rates.append([RATES[minarg]])
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
                #inds = np.where(self.scaps[i] < 0.075)[0]
                inds = np.where(fcaps[i] < self.fcap_min)[0]
                if len(inds) > 0:
                    caps[i] = np.delete(caps[i], inds)
                    scaps[i] = np.delete(scaps[i], inds)
                    fcaps[i] = np.delete(fcaps[i], inds)
                    eff_rates[i] = np.delete(eff_rates[i], inds)
                    rates[i] = np.delete(rates[i], inds)
                    cutvolts[i] = np.delete(cutvolts[i], inds)
                    currs[i] = np.delete(currs[i], inds)
                    ir[i] = np.delete(ir[i], inds)
                    dqdv[i] = np.delete(dqdv[i], inds)
                    resistdrop[i] = np.delete(resistdrop[i], inds)
                    print("Current removed due to being below fcap min")
        
            if self.capacitance_corr == True:
                print("Capacitance correction cannot be applied to multi-current AMID data. Data is being analyzed without capacitance correction")
        else:
            #icaps is the idealized capacity for a given voltage based upon dqdv
            #currsCum is the cumulative averge current for determining where C should be calculated
            icaps = []
            currsCum = []
            voltsAct = []
            time = []
            for i in range(nsig):
                step = sigs.loc[sigs['Prot_step'] == sigsteps[i]]
                #step = step.loc[step['Potential'] > step['Potential'].values[-1] + 0.01]
                stepcaps = step['Capacity'].values
                volts = step['Potential'].values
                currents = np.absolute(step['Current'].values)
                runtime = step['Time'].values
                
                #Collect suceeding OCV steps (1 OCV or 2 OCV) to calculate dqdv
                ocvstep = self.sigdf.loc[self.sigdf['Prot_step'] == sigsteps[i] + 2]
                if ocvstep['Step'].values[0] != 0:
                    ocvstep = self.sigdf.loc[self.sigdf['Prot_step'] == sigsteps[i] + 1]
                ocvstepcaps = ocvstep['Capacity'].values
                ocvvolts = ocvstep['Potential'].values
                
                # if less than 4(NVX) or 5(UHPC) data points, immediate voltage cutoff reached, omit step.
                if len(currents) > 3:
                    if volts[-2] == np.around(volts[-2], decimals=2):
                        cvoltind = -2
                    elif len(currents) > 4:
                        if volts[-3] == np.around(volts[-3], decimals=2):
                            cvoltind = -3
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
                
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
                
                if cutvolts == []:
                    prevstep = self.df.loc[self.df['Prot_step'] == sigsteps[i] - 1]
                    prevvolts = prevstep['Potential'].values
                    initcutvolt = np.around(prevvolts[-1], decimals=3)
                cutvolts.append([ocvvolts[-1]]) 
                
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
                    #if caps is calculated as negative or is first datapoint, this datapoint is effectively thrown out (caps set to 0)
                    for j in range(len(caps[i])):
                        if caps[i][j] < 0 or j == 0:
                            caps[i][j] = 0
                    
                    icaps[i] = icaps[i] - dlcaps - dqdv[i][0]*currs[i]*rohm 
                    #if icaps is calculated as negative or zero or is first datapoint, this datapoint is effectively thrown out (caps set to nan)
                    for j in range(len(caps[i])):
                        if icaps[i][j] <= 0 or j == 0:
                            icaps[i][j] = float('NaN')
                    
                    currsCum[i] = currsCum[i] - dlcaps/time[i]
                    #if cumulative current is calculated as negative or zero or is first datapoint, this datapoint is effectively thrown out (caps set to nan)
                    for j in range(len(caps[i])):
                        if currsCum[i][j] <= 0 or j == 0:
                            currsCum[i][j] = float('NaN')
            
            for i in range(nvolts):
                fcaps.append(caps[i]/icaps[i])
                eff_rates.append(icaps[i]/currsCum[i])
                #outlier repair: if fcap or eff_rates is greater than succeeding points or NaN, reduce it to be equal to lowest of succeeding points.
                for j in range(len(caps[i]) - 1):
                    for k in range(len(caps[i]) - j - 1):
                        if not(fcaps[i][j] <= fcaps[i][j + k + 1]):
                            fcaps[i][j] = fcaps[i][j + k + 1]
                        if not(eff_rates[i][j] <= eff_rates[i][j + k + 1]):
                            eff_rates[i][j] = eff_rates[i][j + k + 1]

        
        print('Found {} signature curves.'.format(nvolts))
        cvolts = np.zeros(nvolts)
        for i in range(nvolts):
            cvolts[i] = np.average(cutvolts[i])
        
                #v1 = np.around(cutvolts[-i-1][-1], decimals=2)
        #if i == 0:
            #cvolts[i] = v1
        #else:
            #v2 = np.around(cutvolts[-i][-1], decimals=2)
            #if v2 == v1:
                #cvolts[i] = 2*cvolts[i-1] - cvolts[i-2]
            #else:
                #cvolts[i] = v1
        #cvolts = cvolts[::-1]
        
        print('Cutoff voltages: {}'.format(cvolts))
        avg_volt = np.zeros(nvolts)
        # Get midpoint voltage for each range.
        avg_volt[0] = (initcutvolt + cvolts[0])/2
        avg_volt[1:] = (cvolts[:-1] + cvolts[1:])/2
        print('Midpoint voltages: {}'.format(avg_volt))
        dvolts = np.zeros(nvolts)
        dvolts[0] = np.absolute(initcutvolt - cvolts[0])
        dvolts[1:] = np.absolute(cvolts[:-1] - cvolts[1:])
        print('Voltage intervals widths: {}'.format(dvolts))
        # Make voltage interval labels for legend.
        vlabels = ['{0:.3f} V - {1:.3f} V'.format(initcutvolt, cvolts[0])]
        vlabels = vlabels + ['{0:.3f} V - {1:.3f} V'.format(cvolts[i], cvolts[i+1]) for i in range(nvolts-1)]
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
                initcutvolt = cutvolts[i-iadj][0]
                cutvolts.pop(i-iadj)
                currs.pop(i-iadj)
                ir.pop(i-iadj)
                dqdv.pop(i-iadj)
                resistdrop.pop(i-iadj)
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

        return new_caps, new_scaps, fcaps, rates, eff_rates, currs, ir, dqdv, resistdrop, cvolts, avg_volt, dvolts, vlabels 
       

    def fit_atlung(self, r, ftol=5e-14, D_bounds=[1e-17, 1e-8], D_guess=1.0e-13, 
                   maxfcap_bounds=[0.95, 1.5], maxfcap_guess=1.0, R_eff_bounds=[1e-6, 1e1], R_eff_guess=1.0e-2,
                   shape='sphere', r_corr=False, nalpha=4000, nQ=4000, export_data=True, export_fig=True, label=None):

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
        if r_corr is False:
            print("Opt params: [Log(D), fCapAdj]")
            Q_arr = np.logspace(-3, 2, nQ)
            tau_sol = np.zeros(nQ)
            tau_guess = 0.5
            for i in range(nQ):
                Q = Q_arr[i]
                func = lambda tau: tau - 1 + (1/(A*Q))*(1/B - 2*(np.sum(np.exp(-self.alphas*tau*Q)/self.alphas)))
                tau_sol[i] = fsolve(func, tau_guess, factor=1.)
        else:
            print("Opt params: [Log(D), fCapAdj, log(Reff), Log(Reff/D)]")
                
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
            
            if self.single_curr is False:
                # selects the dqdv of C/40 discharge/charge or nearest to C/40
                act_rates = self.capacity / np.array(self.currs[j])
                minarg = np.argmin(np.absolute(40 - act_rates))
                dqdv[j] = self.dqdv[j][minarg]
            else:
                dqdv[j] = self.dqdv[j][0]
                # constrains maxfcap to 1
                maxfcap_bounds=[1.0, 1.0000001]

            #print("dq/dV: {} Ah/V".format(dqdv[j]))
            
            if self.r_corr is False:
                C = np.sum(self.ir[j])
                weights = (C - self.ir[j]) / np.sum(C - self.ir[j])
                bounds = ([np.log10(D_bounds[0]), maxfcap_bounds[0]],
                              [np.log10(D_bounds[1]), maxfcap_bounds[1]])
                p0 = [np.log10(D_guess), maxfcap_guess]    
            else:
                bounds = ([np.log10(D_bounds[0]), maxfcap_bounds[0], np.log10(R_eff_bounds[0])],
                          [np.log10(D_bounds[1]), maxfcap_bounds[1], np.log10(R_eff_bounds[1])])
                p0 = [np.log10(D_guess), maxfcap_guess, np.log10(R_eff_guess)]
                
            with plt.style.context('grapher'):
                fig = plt.figure()
                
                if shape == 'sphere':
                    if self.r_corr is False:
                        popt, pcov = curve_fit(self._spheres, (fcap, rates), z, p0=p0,
                                   bounds=bounds, sigma=weights,
                                   method='trf', max_nfev=5000, x_scale=[1.0, 1.0],
                                   ftol=ftol, xtol=None, gtol=None, loss='soft_l1', f_scale=1.0)
                        print("{} Opt params: {}".format(self.vlabels[j], popt))
                    else:
                        p0opt = [p0[2] - p0[0], p0[1], p0[2]]
                        boundsopt = [[bounds[0][2] - bounds[1][0], bounds[0][1], bounds[0][2]], 
                                      [bounds[1][2] - bounds[0][0], bounds[1][1], bounds[1][2]]]
                        
                        popt, pcov = curve_fit(self._spheres_r_corr, (fcap, rates), z, p0=p0opt,
                                   bounds=boundsopt,
                                   method='trf', max_nfev=5000, x_scale=[1.0, 1.0, 1.0],
                                   ftol=ftol, xtol=None, gtol=None, loss='soft_l1', f_scale=1.0)
                        popt = np.array([popt[2] - popt[0], popt[1], popt[2], popt[0]])
                        print("{} Opt params: {}".format(self.vlabels[j], popt))
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
                    plt.semilogx(Q_arr, 0.001*np.ones(len(Q_arr)), ':k')
                elif max(tau_fit) < 0.1:
                    plt.ylim(0, 0.1)
                    plt.semilogx(Q_arr, 0.001*np.ones(len(Q_arr)), ':k')
                plt.semilogx(Q_arr, tau_sol, '-k', label='Atlung - {}'.format(shape))
                plt.xlabel(r'$Q = 3600 n_{eff} D / r^2$')
                plt.ylabel('Fractional Capacity')
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
            #cols = ['Voltage', 'D']
        else:
            # get resist from resist_eff
            resist = resist_eff*self.r**2/(3600*dconst*dqdv)
            #print("Resist: {} V/A".format(resist))
            
            # get resist from ir drop
            rdrop = []
            for i in range(len(resist)):
                if self.single_curr is False:
                    rdrop.append(np.average(self.resistdrop[i]))
            #        resistdrop.append([])
            #        for j in range(len(self.ir[i])):
            #            resistdrop[i][j].append(self.ir[i][j]/self.currs[i][j])
                else:
                    rdrop.append(self.resistdrop[i][0])
            
            DV_df = pd.DataFrame(data={'Voltage (V)': self.avg_volts, 'D (cm^2/s)': dconst, 'R_eff' : resist_eff,
                                       'dqdV (mAh/gV)': [i*1000/self.mass for i in dqdv], 'Rfit (Ohm)' : resist, 'Rdrop (Ohm)' : rdrop})
        
        if export_data is True:
            if label is None:
                df_filename = self.dst / '{0}_D-V_{1}.xlsx'.format(self.cell_label, shape)
            else:
                df_filename = self.dst / '{0}-{1}_D-V_{2}.xlsx'.format(self.cell_label, label, shape)
                
            #DV_df.to_excel(df_filename, columns=cols, index=False)
            DV_df.to_excel(df_filename, index=False)
        
        print('Fitted Dc: {}'.format(dconst))
        print('Standard deviations from fit: {}'.format(sigma))
        print('Atlung fit error: {}'.format(fit_err))
        
        return self.avg_volts, dconst, fit_err, cap_span, cap_max, cap_min, self.caps, self.ir, self.dvolts, resist_eff, dqdv, resist, self.resistdrop
        
    def make_summary_graph(self, fit_data, export_fig=True, label=None):
        
        voltage = fit_data[0]
        nvolts = len(voltage)
        dconst = fit_data[1]
        fit_err = fit_data[2]
        cap_span = fit_data[3]
        cap_max = fit_data[4]
        cap_min = fit_data[5]
        caps = fit_data[6]
        dV_ir = fit_data[7]
        dvolts = fit_data[8]
        resist_eff = fit_data[9]
        dqdv = fit_data[10]
        resist = fit_data[11]
        resistdrop = fit_data[12]
        
        with plt.style.context('grapher'):
            if self.single_curr is False:
                if self.r_corr is False:
                    fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(7, 15), sharex=True,
                    gridspec_kw={'height_ratios': [2,1,1,1,1], 'hspace': 0.0})
                    
                    axs[0].semilogy(voltage, dconst, 'ko--', linewidth=0.75, label='{}'.format(self.cell_label))
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
                    
                    axs[0].semilogy(voltage, dconst, 'ko-', linewidth=0.75, label='{}'.format(self.cell_label))
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
                
                axs[0].semilogy(voltage, dconst, 'k+-', linewidth=0.75, label='{}'.format(self.cell_label))
                axs[0].set_xlabel('Voltage (V)', fontsize=12)
                axs[0].set_ylabel('D (cm$^2$/s)', fontsize=12)
                axs[0].legend(frameon=False, fontsize=12)
                
                axs[1].semilogy(voltage, resist_eff, 'k+-', linewidth=0.75)
                axs[1].semilogy(voltage, 1.0/15*np.ones(len(voltage)), 'k:', linewidth=1.5)
                axs[1].set_ylim(1.0e-6, 1.0)
                axs[1].set_xlabel('Voltage (V)', fontsize=12)
                axs[1].set_ylabel('R$_{eff}$', fontsize=12)
                
                axs[2].semilogy(voltage, resist, 'k+-', linewidth=0.75, label='fit R')
                axs[2].semilogy(voltage, resistdrop, 'k1-', linewidth=0.75, label='Vdrop R')
                axs[2].set_ylim(1.0e1, 0.99e5)
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
                
                axs[5].plot(voltage, [i*1000/self.mass for i in dqdv], 'k+-', linewidth=0.75)
                axs[5].set_xlabel('Voltage (V)', fontsize=12)
                axs[5].set_ylabel('dq/dV (mAh/gV)', fontsize=12)
        
            if export_fig is True:
                if label is None:
                    figstr = 'D-V_{0}.jpg'.format(self.cell_label)
                else:
                    figstr = 'D-V_{0}_{1}.jpg'.format(self.cell_label, label)
                diff_figname = self.dst / figstr
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
        
        #Calculates without points with fcap less than 0.001 if the largest fcap is less than 0.1
        #This allows better fits by treating the resistance growth at low voltages as an instantaneous ohmic resistance
        #when total fcap is small enough that this makes a difference in the fit.
        #Calculates inacessible capacity as 1 + tau if R_eff/Q > 1 AND fcap is less than 0.05 of largest fcap 
        #by setting n so that R_eff=Q. Otherwise standard AMIDR equation.
        #This avoids the divergent region where tau=0 but infinite summation error is amplified. 
        result = []
        for i in range(len(c)):
            if c[i] < 0.001 and max(c) < 0.1:
                result.append(1)
            elif R_eff>(3600*n[i]*D)/self.r**2 and c[i]/max(c) < 0.05:
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