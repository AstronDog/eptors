#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon 24 12:29:17 2020
Author: Jun Wang & Caterina Tiburzi
E-mail: jun.wang.ucas@gmail.comi

Main Function:
    EPTA Pulsar Timing Outlier Rejection Scheme(eptors.py). Mainly used for eliminating incorrect and biased TOAs.    
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.linear_model import HuberRegressor
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from statsmodels import robust
import pandas as pd
import os, argparse, subprocess


def res_obtain(obs, eph):
    os.system('''cat %s | awk '{if ($1 != "C" && $1 != "MODE") print $0}' > %s_tmp'''%(obs,obs))

    os.system('''tempo2 -output general2 -s "tadan {freq} {post} {err}e-6\n" -f %s %s_tmp -npsr 1 \
              -nobs 50000 | grep 'tadan'| awk '{print $2,$3,$4}' > tmp.dat'''%(eph,obs))
    freq, post, err = np.genfromtxt("./tmp.dat")[:,0], np.genfromtxt("./tmp.dat")[:,1], np.genfromtxt("./tmp.dat")[:,2]
    freq = freq.reshape(-1,1)
    return freq, post, err

def show_residual(obs, eph):
    subprocess.call('tempo2 -gr plk -f %s %s -npsr 1 -nobs 50000 -showchisq'%(eph, obs), shell=True)


def trim(obs,eph):
    freq, post, err = res_obtain(obs, eph)

    model = make_pipeline(PolynomialFeatures(2), HuberRegressor()) 
    model.fit(freq, post, huberregressor__sample_weight=np.ravel(1/err))     

    y_pred = model.predict(freq)
    residuals = post - y_pred

    median = np.median(residuals)
    MAD = robust.mad(residuals)
    in_mask = (residuals > median - 3*MAD) & (residuals < median + 3*MAD)
    out_mask = (residuals <= median - 3*MAD) | (residuals >= median + 3*MAD)   


    trim_del_toa,trim_del_err,trim_del_freq = post[out_mask], err[out_mask], freq[out_mask]    
    trim_keep_freq,trim_keep_toa,trim_keep_err = freq[in_mask],post[in_mask],err[in_mask]
 

    df = pd.read_csv("%s_tmp"%(obs), skiprows=1, dtype=str, header=None, delim_whitespace=True)

    df['inliers'] = in_mask
    df = df[df.inliers != False]
    del df['inliers']

    df.to_csv("trim.tim", sep=' ', header=False, index=False)
    os.system("sed -i '1iFORMAT 1' trim.tim")
    #Remove mediators
    os.system("rm tmp.dat  %s_tmp" %(obs))
    
    return trim_del_toa, trim_del_err, trim_del_freq, trim_keep_toa, trim_keep_err, trim_keep_freq

def trim_plot(obs, eph):    
    freq, post, err = res_obtain(obs, eph)
    trim_del_toa, trim_del_err, trim_del_freq, trim_keep_toa, trim_keep_err, trim_keep_freq = trim(obs,eph)
   
    # PLot the fitted data with Huber regression model
    model = make_pipeline(PolynomialFeatures(2), HuberRegressor())
    model.fit(freq, post, huberregressor__sample_weight=np.ravel(1/err))
    x_plot = np.linspace(freq.min(), freq.max())
    y_plot = model.predict(x_plot[:, np.newaxis])
    plt.plot(x_plot,y_plot,'k-')
    
    #Plot the outlier residuals
    plt.errorbar(trim_del_freq, trim_del_toa, yerr=trim_del_err,fmt='r.', label="Outliers")

    #Plot the Inliers
    plt.errorbar(trim_keep_freq,trim_keep_toa,yerr=trim_keep_err,fmt='b.',label="Inliers")
    
    plt.title("Median absolute deviation method")
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Post-fit Residuals (sec)')
    plt.legend()
    plt.show()

def gauss_fit(obs, eph):       
    freq, post, err = res_obtain(obs, eph)
       
    rvu = post/err

    num_bins = 500
    data = plt.hist(rvu, num_bins, facecolor='blue', alpha=0.5)

    def f(x, b, c, d):
        return  b * np.exp(-(x - c)**2.0/(2 * d**2))

    x = [0.5 * (data[1][i] + data[1][i+1]) for i in np.arange(len(data[1])-1)]
    y = data[0]



    ar,mean,rms = np.trapz(y/np.sum(y)*len(y), x),np.mean(x * y / np.sum(y) * len(y)),np.std(x * y / np.sum(y) * len(y))
    popt, pcov = curve_fit(f, x, y, p0=[ar, mean, rms])
    


    in_mask =  (rvu > popt[1] - 4 * abs(popt[2])) & (rvu < popt[1] + 3 * abs(popt[2]))
    out_mask =  (rvu <= popt[1] - 4 * abs(popt[2])) | (rvu >= popt[1] + 3 * abs(popt[2]))
    
    gaus_keep_toa,gaus_keep_err,gaus_keep_freq = post[in_mask],err[in_mask],freq[in_mask]   
    gaus_del_toa,gaus_del_err,gaus_del_freq = post[out_mask],err[out_mask],freq[out_mask]
   

    df = pd.read_csv("%s_tmp"%(obs), skiprows=1, dtype=str, header=None, delim_whitespace=True)
  
    df['inliers'] = in_mask
    df = df[df.inliers != False]
    del df['inliers']

    df.to_csv("trim_gauss.tim", sep=' ', header=False, index=False)
    os.system("sed -i '1iFORMAT 1' trim_gauss.tim")
       
    os.system("rm %s_tmp"%(obs))
    plt.close()
    return gaus_del_freq, gaus_del_toa, gaus_del_err, gaus_keep_toa, gaus_keep_err, gaus_keep_freq


def gauss_plot(obs, eph):
    freq, post, err = res_obtain("trim.tim", eph)
    gaus_del_freq, gaus_del_toa, gaus_del_err, gaus_keep_toa, gaus_keep_err, gaus_keep_freq  =  gauss_fit("trim.tim", eph)    
    

    plt.errorbar(gaus_keep_freq,gaus_keep_toa,yerr=gaus_keep_err,fmt='b.',label="Inliers")
    plt.errorbar(gaus_del_freq, gaus_del_toa, yerr=gaus_del_err,fmt='c.', label="Outliers")
    
    plt.title("Gaussian distribution fitting method")
    plt.ylabel('Timing residual / ToA uncertainty')
    plt.legend()
    plt.show()

    
def showeffect(obs, eph, algo): 
    """
    Plot the final result of the total outlier rejection scheme
    Reserved TOAs are plot in blue.
    ToAs eliminate with MAD scheme are in red.
    ToAs eliminate with Gaussian distribution fit scheme are in red.
    """
    trim_del_toa, trim_del_err, _, trim_keep_toa, trim_keep_err, _ = trim(obs, eph)

        
    if algo == "t":   
        plt.errorbar(np.arange(0, len(trim_keep_toa), 1), trim_keep_toa, yerr=trim_keep_err, fmt='b.', label='Inliers')
    
        plt.errorbar(np.arange(len(trim_keep_toa), len(trim_keep_toa) + len(trim_del_toa), 1), 
                               trim_del_toa, yerr=trim_del_err, fmt='r.', label='Trim-removed outliers')
        plt.xlim((-4, len(trim_del_toa) + len(trim_keep_toa) + 4)) 
    else:
        _,gaus_del_toa, gaus_del_err, gaus_keep_toa, gaus_keep_err,_ = gauss_fit("trim.tim", eph)
        
        plt.errorbar(np.arange(0, len(gaus_keep_toa), 1), gaus_keep_toa, yerr=gaus_keep_err, fmt='b.', label='Inliers')
        
        plt.errorbar(np.arange(len(gaus_keep_toa), len(gaus_keep_toa) + len(trim_del_toa), 1), 
                               trim_del_toa, yerr=trim_del_err, fmt='r.', label='Trim-removed outliers')
        
        plt.errorbar(np.arange(len(gaus_keep_toa) + len(trim_del_toa), 
                               len(gaus_keep_toa) + len(trim_del_toa) + len(gaus_del_toa), 1), 
                               gaus_del_toa, yerr=gaus_del_err, fmt='c.', label='Gauss-removed outliers')
        
        plt.xlim((-4, len(gaus_keep_toa) + len(trim_del_toa) + len(gaus_del_toa) + 4))
        
    plt.title("Outlier removed by each method")
    plt.legend(loc='upper left')
    plt.ylabel('Post-fit Residual (sec)')
    plt.xlabel('TOA number')


    ax = plt.gca() 
    plt.ticklabel_format(axis='y', style='sci') 
    ax.yaxis.major.formatter.set_powerlimits((0,0))    
    plt.show()


def parse_arguments():    
    parser = argparse.ArgumentParser(description="Outlier Rejection for EPTA pulsar timing")
    parser.add_argument('-t', '--tim', type=str, nargs=1, help="Name of the TOA file")
    parser.add_argument('-e', '--ephemeris', type=str, nargs=1, help="Selects parameter file")
    parser.add_argument('-np', '--noplot', action='store_true', dest="noplot", default=False,
                        help='Do not display plots')
    parser.add_argument('-o', '--output', type=str, nargs=1, help="The output file's name.")
    parser.add_argument('-A', '--algorithms', type=str, nargs=1, default="a", choices=['a', 't'], 
                        help="Outlier rejection algorithms, where\n"
                        "a = trim and gauss methods are both applied,\n"
                        "t = only trim method is used")                  

    args = parser.parse_args()
    return args    
    
    
def main(args):    
    obs = args.tim[0]
    eph = args.ephemeris[0]
    noplot = args.noplot
    output = args.output
    algo = args.algorithms[0]

    if algo == "t":
        if noplot:
            trim(obs, eph)
        else:
            show_residual(obs, eph)
            trim(obs, eph)
            trim_plot(obs, eph)
            show_residual("trim.tim", eph)
            showeffect(obs, eph, "t") 
        if output:
            os.system("mv trim.tim %s" %output[0])  
        else:
            pass             
    else:
        if noplot:
            trim(obs, eph)
            gauss_fit("trim.tim", eph)
        else:
            show_residual(obs, eph)
            trim(obs, eph)
            trim_plot(obs, eph)
            gauss_fit("trim.tim", eph)
            gauss_plot("trim.tim", eph)
            show_residual("trim_gauss.tim", eph)
            showeffect(obs, eph, "a")
        if output:
            os.system("mv trim_gauss.tim %s" %output[0])
        else:
            pass

        
        
if __name__=="__main__":
    args = parse_arguments()
    main(args)
