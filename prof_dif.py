#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 17:45:11 2019
Author: Jun Wang
E-mail: jun.wang.ucas@gmail.com

Main Function: Show the difference of templates(OA,PA,SA) to make
               sure analytic noise free template are well fitted.

"""
import psrchive as pr
import argparse
import numpy as np
import matplotlib.pyplot as plt



def main(args):
    global add_tmp,addsm_tmp, paas_tmp    
    add_tmp = args.ad[0]
    addsm_tmp = args.sm[0]
    paas_tmp = args.pa[0]
    
    showdif()



def get_value(archive):
    arch = pr.Archive_load(archive)
    arch.tscrunch()
    arch.dedisperse()
    arch.fscrunch()
    arch.pscrunch()
    arch.remove_baseline()
#    data = arch[0].get_Profile(0.0).get_amps()
    data = arch.get_data()
    return data[0,0,0,:]

    
def parse_arguments():
    parser = argparse.ArgumentParser(description='Compare the template profiles.')
    parser.add_argument('-ad', type=str, nargs=1, help='The added template') 
    parser.add_argument('-sm', type=str, nargs=1, help='The smoothed template')
    parser.add_argument('-pa', type=str, nargs=1, help='The paas template')
    parser.add_argument('-o', '--outplot', action='store_true', help='Creates a plot that show the template difference')    
    args = parser.parse_args()
    return args 



def showdif():
    add_prof = get_value(add_tmp)
    addsm_prof = get_value(addsm_tmp)
    paas_prof = get_value(paas_tmp)
    
    psr_nbin = paas_prof.size
    peak = max([np.max(add_prof), np.max(addsm_prof), np.max(paas_prof)])
    
    
    add_addsm = [add_prof[i] - addsm_prof[i] for i in range(len(add_prof))]
    add_paas = [add_prof[i] - paas_prof[i] for i in range(len(add_prof))]
    addsm_paas = [addsm_prof[i] - paas_prof[i] for i in range(len(addsm_prof))]

    bins = np.linspace(0, 1, psr_nbin)
    
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(18, 9))
    #Plot the difference between the added and added smoothed templates
    axes[0].plot(bins, np.array(add_addsm)/peak, color='orange', linewidth=1.0, linestyle="-")
    axes[0].set_title('OA-SA',fontsize=18)
    
    #Plot the difference between the added and paas templates    
    axes[1].plot(bins, np.array(add_paas)/peak, color='b', linewidth=1.0, linestyle="-")
    axes[1].set_title('OA-PA',fontsize=18)
    
    #Plot the difference between the addsm and paas templates    
    axes[2].plot(bins, np.array(addsm_paas)/peak, color='#003333', linewidth=1.0, linestyle="-")
    axes[2].set_title('SA-PA',fontsize=18)
    
    #Plot the difference between the addsm and paas templates    
    axes[3].plot(bins, add_prof, color='c', linewidth=1.0, linestyle="-")
    axes[3].set_xlabel('Pulse Phase',fontsize=18)
    axes[3].set_title('Pulse Profile',fontsize=18)
    
    plt.xticks(fontsize=18)
    for axis in axes:
        for tick in axis.yaxis.get_major_ticks():
            tick.label1.set_fontsize(16)
        
    #plt.yticks(fontsize=18)

    
    
    fig.text(0.06, 0.2, 'Intensity', ha='center', va='center', rotation='vertical',fontsize=18)
    fig.text(0.06, 0.62, 'Normalized Profile Residuals', ha='center', va='center', rotation='vertical',fontsize=18)
    
    if args.outplot:
        plt.savefig("./out_plot.eps", format='eps', dpi=500)
    plt.show()    
    



    

    
if __name__=="__main__":
    args = parse_arguments()
    main(args)    
