#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 07:44:54 2021

@author: Marc M. E. Cormier
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit, fsolve
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


RATES = np.array([0.01, 0.05, 0.1, 0.2, 1/3, 0.5, 1, 2, 2.5, 5, 10, 20, 40, 80,
                  160, 320, 640, 1280])

COLUMNS = ['Time', 'Cycle', 'Step', 'Current', 'Potential', 'Capacity', 'Prot_step']
UNITS = ['(h)', None, None, '(mA)', '(V)', '(mAh)', None]

SHAPES = ['sphere', 'plane']

class AMID():
    
    def __init__(self, dstpath, srcpath, uhpc_files, cell_label):
        self.cell_label = cell_label
        self.dst = Path(dstpath) / self.cell_label
        # If does not exist, create dir.
        if self.dst.is_dir() is False:
            self.dst.mkdir()
            print('Create directory: {}'.format(self.dst))
        self.src = Path(srcpath)
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
            
        with open(self.uhpc_file, 'r') as f:
            f.readline()
            self.cellname = f.readline().strip().split()[-1]
            f.readline()
            f.readline()
            self.mass = float(f.readline().strip().split()[-1]) / 1000
            self.input_cap = float(f.readline().strip().split()[-1]) / 1000
            for i in range(4):
                f.readline()
            if f.readline().strip().split()[0] == '[Data]':
                hlinenum = 11
                hline = f.readline()
            else:
                hlinenum = 12
                f.readline()
                hline = f.readline()
               
        print('Working on cell: {}'.format(self.cellname))
        print('Positive electrode active mass: {} g'.format(self.mass))
        print('Input cell capacity: {} Ah'.format(self.input_cap))
        
        self.df = pd.read_csv(self.uhpc_file, header=hlinenum)
        # Add Prot.Step column if missing.
        if hline.strip()[-4:] == 'Flag':
            self.df = self.df.rename(columns={'Flag':'Prot.Step'})
            i = self.df.Step
            self.df['Prot.Step'] = i.ne(i.shift()).cumsum() - 1
            
        self.df.columns = COLUMNS
        
        # Remove data where time is not increasing.   
        t = self.df['Time'].values
        dt = t[1:] - t[:-1]
        inds = np.where(dt < 0.0)[0]
        print('Indices being removed to time non-monotonicity: {}'.format(inds))
        self.df = self.df.drop(inds+1)
        inds = self.df.index[self.df['Potential'] < 0.0].tolist()
        print('Indices being removed due to negative voltage: {}'.format(inds))
        self.df = self.df.drop(inds)
        
        self.df.columns = COLUMNS
        #self.df['Capacity'] = self.df['Capacity']
        #self.df['Current'] = self.df['Current']
        
        self.sigdf = self._find_sigcurves()
        self.sc_stepnums = self.sigdf['Prot_step'].unique()
        self.capacity = self.sigdf['Capacity'].max() - self.sigdf['Capacity'].min()
        self.spec_cap = self.capacity / self.mass
        print('Specific Capacity: {0:.2f} mAh'.format(self.spec_cap*1000))
        
        self.caps, self.rates, self.eff_rates, self.currs, self.ir, self.cvolts, self.avg_volts, self.dvolts, self.vlabels = self._parse_sigcurves()
        self.nvolts = len(self.caps)
        self.caps = 1000*np.array(self.caps)/self.mass
        self.scaps = []
        self.fcaps = []
        for i in range(self.nvolts):
            self.scaps.append(np.cumsum(self.caps[i]))
            self.fcaps.append(np.cumsum(self.caps[i]) / np.sum(self.caps[i]))
            
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
        if steps[ocv_inds[0] - 1] == steps[ocv_inds[0] + 1]:
            first_sig_step = prosteps[ocv_inds[0] - 1]
        else:
            first_sig_step = prosteps[ocv_inds[0] + 1]
        #print(first_sig_step)
        last_sig_step = ocv_inds[-1] + 1
        #print(last_sig_step)
        sigdf = self.df.loc[(self.df['Prot_step'] >= first_sig_step) & (self.df['Prot_step'] <= last_sig_step)]
        
        return sigdf
    
    def plot_protocol(self, save=True):
        
        with plt.style.context('grapher'):
        
            fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True,
                                    figsize=(8, 4), gridspec_kw={'wspace':0.0})
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
                
                avgcurr = np.absolute(avgcurr)
                minarg = np.argmin(np.absolute(RATES - self.capacity/avgcurr))
                rate = RATES[minarg]
                axs[1].plot(stepdf['Capacity']*1000/self.mass, stepdf['Potential'],
                            color=colors[c],
                            label='C/{0} {1}'.format(int(rate), cyclabel))
                c = c + 1
                
                # if the next step is the start of sigcurves, plot sigcurves
                if fullsteps[i] == self.sc_stepnums[0] - 1:
                    #print('plottting sig curves...')
                    axs[1].plot(self.sigdf['Capacity']*1000/self.mass, self.sigdf['Potential'],
                                color='red',
                                label='Sig Curves')
            plt.legend()
            
            if save is True:
                plt.savefig(self.dst / 'protocol_vis_{}.jpg'.format(self.cell_label))
            else:
                plt.show()
        
    
    def plot_caps(self, save=True):
        
        with plt.style.context('grapher'):
            fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True,
                                    figsize=(7, 8), gridspec_kw={'hspace':0.0})
            colors = plt.get_cmap('viridis')(np.linspace(0,1,self.nvolts))
            
            for i in range(self.nvolts):
                axs[0].semilogx(self.eff_rates[i], self.scaps[i],
                                color=colors[i], label=self.vlabels[i])
                axs[1].semilogx(self.eff_rates[i], self.fcaps[i],
                                color=colors[i])
                
            axs[1].set_xlabel('n$_{eff}$ in C/n$_{eff}$')
            axs[1].set_ylabel('Fractional Capacity')
            axs[0].set_ylabel('Specific Capacity (mAh/g)')
            axs[0].legend(frameon=False, bbox_to_anchor=(1.0, 0.0), loc='center left')
            axs[0].tick_params(direction='in', which='both', top=True, right=True)
            axs[1].tick_params(direction='in', which='both', top=True, right=True)
            
            if save is True:
                plt.savefig(self.dst / 'cap-rate_{}.jpg'.format(self.cell_label))
            else:
                plt.show()

    
    def _parse_sigcurves(self):

        sigs = self.sigdf.loc[self.sigdf['Step'] != 0]
        #capacity = sigs['Capacity'].max() - sigs['Capacity'].min()
        #print('Specific Capacity: {} mAh'.format(capacity))
        Vstart = np.around(sigs['Potential'].values[0], decimals=2)
        Vend = np.around(sigs['Potential'].values[-1], decimals=2)
        print('Starting voltage: {:.3f} V'.format(Vstart))
        print('Ending voltage: {:.3f}'.format(Vend))
        
        sigsteps = sigs['Prot_step'].unique()
        nsig = len(sigsteps)
        print('Found {} charge or discharge steps in sig curve sequences.'.format(nsig))
        caps = []
        rates = []
        cutvolts = []
        currs = []
        ir = []
        for i in range(nsig):
            step = sigs.loc[sigs['Prot_step'] == sigsteps[i]]
            stepcaps = step['Capacity'].values
            volts = step['Potential'].values
            currents = np.absolute(step['Current'].values)
            rate = self.capacity / np.average(currents)
            minarg = np.argmin(np.absolute(RATES - rate))
            if i == 0:
                caps.append([np.amax(stepcaps) - np.amin(stepcaps)])
                rates.append([RATES[minarg]])
                cutvolts.append([volts[-1]])
                currs.append([np.average(currents)])
                ir.append([np.absolute(volts[0] - volts[1])])
            else:
                if np.amax(currents) < currs[-1][-1]:
                    caps[-1].append(np.amax(stepcaps) - np.amin(stepcaps))
                    rates[-1].append(RATES[minarg])
                    cutvolts[-1].append(volts[-1])
                    currs[-1].append(np.amax(currents))
                    ir[-1].append(np.absolute(volts[0] - volts[1]))
                else:
                    caps.append([np.amax(stepcaps) - np.amin(stepcaps)])
                    rates.append([RATES[minarg]])
                    cutvolts.append([volts[-1]])
                    currs.append([np.average(currents)])
                    ir.append([np.absolute(volts[0] - volts[1])])
                    
        nvolts = len(caps)
        cvolts = np.zeros(nvolts)
        for i in range(len(caps)):
            v1 = np.around(cutvolts[-i-1][-1], decimals=2)
            if i == 0:
                cvolts[i] = v1
            else:
                v2 = np.around(cutvolts[-i][-1], decimals=2)
                if v2 == v1:
                    cvolts[i] = 2*cvolts[i-1] - cvolts[i-2]
                else:
                    cvolts[i] = v1
        cvolts = cvolts[::-1]
        print('Cutoff voltages: {}'.format(cvolts))
        avg_volt = np.zeros(nvolts)
        # Get midpoint voltage for each range.
        avg_volt[0] = (Vstart + cvolts[0])/2
        avg_volt[1:] = (cvolts[:-1] + cvolts[1:])/2
        print('Midpoint voltages: {}'.format(avg_volt))
        dvolts = np.zeros(nvolts)
        dvolts[0] = np.absolute(Vstart - cvolts[0])
        dvolts[1:] = np.absolute(cvolts[:-1] - cvolts[1:])
        print('Voltage intervals widths: {}'.format(dvolts))
        # Make voltage interval labels for legend.
        vlabels = ['{0:.2f} V - {1:.2f} V'.format(Vstart, cvolts[0])]
        vlabels = vlabels + ['{0:.2f} V - {1:.2f} V'.format(cvolts[i], cvolts[i+1]) for i in range(nvolts-1)]
        print('Voltage interval labels: {}'.format(vlabels))
        print('Found {} voltage intervals.'.format(nvolts))
        
        eff_rates = []
        vcaps = np.zeros(nvolts, dtype=float)
        for m in range(nvolts):
            nrates = len(currs[m])
            #nrates = len(rates[m])
            eff_rates.append([])
            vcaps[m] = np.sum(caps[m])
            for n in range(nrates):
                eff_rates[-1].append(vcaps[m]/currs[m][n])           
        
        
        return caps, rates, eff_rates, currs, ir, cvolts, avg_volt, dvolts, vlabels 
       

    def fit_atlung(self, r, ftol=5e-14, shape='sphere', nalpha=150, nQ=2000, save=True, label=None):
        
        self.r = r
        
        if shape not in SHAPES:
            print('The specified shape {0} is not supported.'.format(shape))
            print('Supported shapes are: {1}. Defaulting to sphere.'.format(SHAPES))
        # Get geometric constants according to particle shape.
        if shape == 'sphere':
            self.alphas = []
            for i in np.arange(4, 600):
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
        Q_arr = np.logspace(-3, 2, nQ)
        tau_sol = np.zeros(nQ)
        tau_guess = 0.5
        for j in range(nQ):
            Q = Q_arr[j]
            func = lambda tau: tau - 1 + (1/(A*Q))*(1/B - 2*(np.sum(np.exp(-self.alphas*tau*Q)/self.alphas)))
            tau_sol[j] = fsolve(func, tau_guess, factor=1.)
        
                
        dconst = np.zeros(self.nvolts, dtype=float)
        sigma = np.zeros(self.nvolts, dtype=float)
        fit_err = np.zeros(self.nvolts, dtype=float)
        cap_max = np.zeros(self.nvolts, dtype=float)
        cap_min = np.zeros(self.nvolts, dtype=float)
        cap_span = np.zeros(self.nvolts, dtype=float)

        for j in range(self.nvolts):
            z = np.ones(len(self.scaps[j]))
            #fcap = np.array(self.fcaps[j])
            scap = np.array(self.scaps[j])
            rates = np.array(self.eff_rates[j])
            C = np.sum(self.ir[j])
            weights = (C - self.ir[j]) / np.sum(C - self.ir[j])
            bounds = ([1e-15, 0.9*np.amax(scap)], [1e-10, 1.5*np.amax(scap)])
            p0 = [1e-13, np.amax(scap)]
            with plt.style.context('grapher'):
                fig = plt.figure()
                plt.semilogx(Q_arr, tau_sol, '-k', label='Atlung - {}'.format(shape))
                if shape == 'sphere':
                    popt, pcov = curve_fit(self._spheres, (scap, rates), z, p0=p0,
                               bounds=bounds, sigma=weights,
                               method='trf', max_nfev=5000, x_scale=[1e-11, np.amax(scap)],
                               ftol=ftol, xtol=None, gtol=None, loss='soft_l1', f_scale=1.0)
                if shape == 'plane':
                    popt, pcov = curve_fit(self._planes, (scap, rates), z, p0=p0,
                               bounds=bounds, sigma=weights,
                               method='trf', max_nfev=5000, x_scale=[1e-11, np.amax(scap)],
                               ftol=ftol, xtol=None, gtol=None, loss='soft_l1', f_scale=1.0)
                
                sigma[j] = np.sqrt(np.diag(pcov))[0]
                dconst[j] = popt[0]
                Qfit = 3600*rates*dconst[j]/r**2
                tau_fit = scap/popt[1]
                
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
                fit_err[j] = np.sum(weights*error)
                
                plt.semilogx(Qfit, tau_fit, 'or', label='{0} - {1}'.format(self.cell_label, self.vlabels[j]))
                plt.xlabel(r'$Q = 3600 n D / r^2$')
                plt.ylabel('Fractional Capacity')
                plt.legend(frameon=False, loc='lower right')
                if save is True:
                    if label is None:
                        figname = self.dst / '{0}_Atlung-{1}_{2:.3f}.jpg'.format(self.cell_label, shape, self.avg_volts[j])
                    else:
                        figname = self.dst / '{0}-{1}_Atlung-{2}_{3:.3f}.jpg'.format(self.cell_label, label, shape, self.avg_volts[j])
                    plt.savefig(figname)
                else:
                    plt.show()
                plt.close()
                
        DVdf = pd.DataFrame(data={'Voltage': self.avg_volts, 'D': dconst})
        df_filename = self.dst / '{0}_D-V_{1}.xlsx'.format(self.cell_label, shape)
        DVdf.to_excel(df_filename, columns=['Voltage', 'D'], index=False)
        
        print('Fitted Dc: {}'.format(dconst))
        print('Standard deviations from fit: {}'.format(sigma))
        print('Atlung fit error: {}'.format(fit_err))
        
        return self.avg_volts, dconst, fit_err, cap_span, cap_max, cap_min, self.caps, self.ir, self.dvolts
        
    def make_summary_graph(self, fit_data, save=True, label=None):
        
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
        
        with plt.style.context('grapher'):
            fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(8, 12), sharex=True,
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
            
            if save is True:
                if label is None:
                    figstr = 'D-V_{0}.jpg'.format(self.cell_label)
                else:
                    figstr = 'D-V_{0}_{1}.jpg'.format(self.cell_label, label)
                diff_figname = self.dst / figstr
                plt.savefig(diff_figname)
                
            else:
                plt.show()
                

    def _spheres(self, X, D, c_max):
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis=1)
        narr = np.repeat(n.reshape(len(n), 1), len(self.alphas), axis=1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis=0)
        
        return c/c_max + ((self.r**2)/(3*3600*n*D))*(1/5 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis=1)))
    
    def _planes(self, X, D, c_max):
        
        c, n = X
        carr = np.repeat(c.reshape(len(c), 1), len(self.alphas), axis=1)
        narr = np.repeat(n.reshape(len(c), 1), len(self.alphas), axis=1)
        a = np.repeat(self.alphas.reshape(1, len(self.alphas)), np.shape(carr)[0], axis=0)
        
        return c/c_max + ((self.r**2)/(3600*n*D))*(1/3 - 2*(np.sum(np.exp(-a*(carr/c_max)*3600*narr*D/self.r**2)/a, axis=1)))
    
    