#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:47:14 2019
Author: Jun Wang & Caterina Tiburzi
E-mail: jun.wang.ucas@gmail.com
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


def ncy_check(obs):
    """
    Throw away the extra IPTA flag('-nch') produced with nancay dataset.

    This will not influence TOAs produced from other telescopes.
    """
    # Loading tim file
    df = pd.read_csv(obs, skiprows=1, header=None,dtype=str, delim_whitespace=True)
    # Check the format of the time file
    if df.shape[1] == 23:
        return
    elif df.shape[1] == 25:
        print("Throw away '-nch' flag in nancay tim file")
        df.drop([21, 22], axis=1, inplace=True)
        df.to_csv(obs, sep=' ', index=False, header=False)
        os.system("sed -i '1iFORMAT 1' %s"%(obs))
    else:
        print("Error: please check the format of the tim file")
        exit()
        
def trim(obs,eph,num):
    """
    Basing on Caterina's code, this function is mainly remove outliers 
    with Median Absolute Deviation(MAD) method.

    It will remove the outliers which are 3*MAD away.
    """

    # Obtain the frequency, post-fit residuals and Uncertainties
    os.system('''cat %s | awk '{if ($1 != "C" && $1 != "MODE") print $0}' > %s_tmp'''%(obs,obs))

    os.system("tempo2 -output general2 -f %s %s_tmp -npsr 1 -nobs 60000 -s \
              \"tadan {freq} {post} {err}e-6\n\" | grep 'tadan'| awk \'{print $2,$3,$4}\' \
              > tmp.dat"%(eph,obs))

    freq, post, err = np.genfromtxt("./tmp.dat")[:,0], np.genfromtxt("./tmp.dat")[:,1], np.genfromtxt("./tmp.dat")[:,2]
    freq = freq.reshape(-1,1)

    #print("Now fitting a parabola to the residuals with Huber regressor")       
    #constructing a linear function fitted basing on the Huber regression and applying it
    model = make_pipeline(PolynomialFeatures(2), HuberRegressor()) 
    model.fit(freq, post, huberregressor__sample_weight=np.ravel(1/err)) 
    
    
    #print("Now identifying and rejecting the outliers")
    y_pred = model.predict(freq)
    residuals = post - y_pred


    median = np.median(residuals)
    MAD = robust.mad(residuals)
    in_mask = (residuals > median - num*MAD) & (residuals < median + num*MAD)
    out_mask = (residuals <= median - num*MAD) | (residuals >= median + num*MAD)

    
    #Get the remaining residuals after triming code
    freq_in = freq[in_mask]
    post_in = post[in_mask]
    err_in = err[in_mask]

    
    #Get the outliers marked by triming code
    global post_to
    global err_to
    post_to = post[out_mask]
    err_to = err[out_mask]


    #Plot the orignal residuals
    plt.errorbar(freq,post,yerr=err,fmt='k.')
    
    # PLot the fitted data with Huber regression model
    x_plot = np.linspace(freq.min(), freq.max())
    y_plot = model.predict(x_plot[:, np.newaxis])
    plt.plot(x_plot,y_plot,'r-')
    
    #Plot the inlier residuals
    plt.errorbar(freq_in,post_in,yerr=err_in,fmt='b.')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Post-fit Residuals (sec)')
    plt.show()

    #print("Now getting the ToAs newly and putting them in a dataframe")
    df = pd.read_csv("%s_tmp"%(obs), skiprows=1, dtype=str, header=None, delim_whitespace=True, 
                     names=["psr_name", "freq", "toa", "err", "site", "key1", "frontend", "key2", 
                            "backend", "key3", "unknown", "key4", "bandwidth", "key5", "length", 
                            "key6", "template", "key7", "gof", "key8", "nbin", "key9", "nch"])


    #print("Now applying the boolean mask")
    df['inliers'] = in_mask
    df1 = df[df.inliers != False]
    del df1['inliers']

  
    #print("Now dumping them to file")
    #Output the new tim file after MAD filtered
    df1.to_csv("%s_trimtim"%(obs), sep=' ', header=False, index=False)
    os.system("sed -i '1iFORMAT 1' %s_trimtim"%(obs))
    #Remove mediators
    os.system("rm tmp.dat")
    os.system("rm %s_tmp"%(obs))




def gauss_fit(obs, eph, num):
    """
    Fitting the Normalized Residuals(Residuals divide by TOA uncertainty) with 
    Gaussian distribution and remove the residuals which are 3*\sigma away 
    from the mean.
    """
    
    obs = obs + '_trimtim'
    
    # Obtain the frequency, post-fit residuals and Uncertainties    
    os.system('''cat %s | awk '{if ($1 != "C" && $1 != "MODE") print $0}' > %s_tmp'''%(obs,obs))
    os.system("tempo2 -output general2 -f %s %s_tmp -npsr 1 -nobs 50000 -s \
              \"tadan {freq} {post} {err}e-6\n\" | grep 'tadan'| awk \'{print $2,$3,$4}\' > tmp.dat"%(eph,obs))


    #print("Now fitting a parabola to the residuals with Huber regressor")
    freq, post, err = np.genfromtxt("./tmp.dat")[:,0], np.genfromtxt("./tmp.dat")[:,1], np.genfromtxt("./tmp.dat")[:,2]
    freq = freq.reshape(-1,1)
    
    #Calculate the normalized residuals, which means residuals divided by ToA uncertainties    
    rvu = post/err

    #Plot the distribution of normalized residuals 
    num_bins = 500
    data = plt.hist(rvu, num_bins, facecolor='blue', alpha=0.5)

    #Definite the Gaussian function and fit the normalized residuals
    def f(x, b, c, d):
        return  b * np.exp(-(x - c)**2.0/(2 * d**2))

    x = [0.5 * (data[1][i] + data[1][i+1]) for i in xrange(len(data[1])-1)]
    y = data[0]


    ar = np.trapz(y/np.sum(y)*len(y), x)
    mean = np.mean(x * y / np.sum(y) * len(y))
    rms = np.std(x * y / np.sum(y) * len(y))
    popt, pcov = curve_fit(f, x, y, p0=[ar, mean, rms]) #fit and return the parameters array popt
    
    x_fit = np.linspace(x[0], x[-1], 5000)
    y_fit = f(x_fit, *popt)
    plt.plot(x_fit, y_fit, lw=3, color="r")
    plt.annotate(r'cut-off point',
             xy=(popt[1] + 3 * abs(popt[2]), 0), xycoords='data',
             xytext=(-20, +30), textcoords='offset points', fontsize=12,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    plt.annotate(r'cut-off point',
             xy=(popt[1] - 3 * abs(popt[2]), 0), xycoords='data',
             xytext=(-40, +30), textcoords='offset points', fontsize=12,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    plt.show()

    in_mask =  (rvu > popt[1] - 3 * abs(popt[2])) & (rvu < popt[1] + 3 * abs(popt[2]))
    out_mask =  (rvu <= popt[1] - 3 * abs(popt[2])) | (rvu >= popt[1] + 3 * abs(popt[2]))
    
    
    #Get the outliers marked by Gaussian function, and use it for show_effect function as well
    global post_le
    global err_le
    freq_le = freq[in_mask]
    post_le = post[in_mask]
    err_le = err[in_mask]
    
    
    #Inlier resiudals, used for show_effect function    
    global post_go
    global err_go
    post_go = post[out_mask]
    err_go = err[out_mask]
    

    #Plot Plot the orignal residuals in black and the inlier residuals in cyan
    plt.errorbar(freq,post,yerr=err,fmt='k.')
    plt.errorbar(freq_le,post_le,yerr=err_le,fmt='c.')
    plt.ylabel('Timing residual / ToA uncertainty')
    plt.xlabel('Percentage')
    plt.show()


    df = pd.read_csv("%s_tmp"%(obs), skiprows=1, dtype=str, header=None, delim_whitespace=True, names=["name",
                     "freq", "toa", "err", "site", "key1", "frontend", "key2", "backend", "key3", "unknown", "key4", 
                     "bandwidth", "key5", "length", "key6", "template", "key7", "gof", "key8", "nbin", "key9", "snr"])

    
    #print("Now applying the boolean mask")    
    df['inliers'] = in_mask
    df = df[df.inliers != False]
    del df['inliers']


    #print("Now dumping them to file")
    #Output the new tim file after Gaussian function filtered
    df.to_csv("%s_gauss"%(obs), sep=' ', header=False, index=False)
    os.system("sed -i '1iFORMAT 1' %s_gauss"%(obs))
    
    #Remove mediators    
    os.system("rm %s_tmp"%(obs))

 
    
def showeffect(): 
    """
    Plot the final result of the total outlier rejection scheme
    Reserved TOAs are plot in blue.
    ToAs eliminate with MAD scheme are in red.
    ToAs eliminate with Gaussian distribution fit scheme are in red.
    """
    x_le = np.arange(0, len(post_le), 1 )
    x_to = np.arange(len(post_le),  len(post_le) + len(post_to), 1)
    x_go = np.arange(len(post_le) + len(post_to), len(post_le) + len(post_go) + len(post_to), 1)

    plt.errorbar(x_le, post_le, yerr=err_le, fmt='b.', label='Inliers')
    plt.errorbar(x_to, post_to, yerr=err_to, fmt='r.', label='Trim-removed outliers')
    plt.errorbar(x_go, post_go, yerr=err_go, fmt='c.', label='Gauss-removed outliers')
    plt.legend(loc='upper left')
    plt.ylabel('Post-fit Residual (sec)')
    plt.xlabel('TOA number')

    ax = plt.gca()
    ax.axhline(linewidth=4, color='k')
    xfmt = plt.ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((0, 0))  # Or whatever your limits are . . .
    plt.gca().yaxis.set_major_formatter(xfmt)
    plt.rcParams.update({'font.size': 30, 'font.family': 'serif'})
    plt.show()
    

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Outlier Rejection for EPTA pulsar timing")
    parser.add_argument('-t', '--tim', type=str, nargs=1, help="Name of the TOA file")
    parser.add_argument('-e', '--ephemeris', type=str, nargs=1, help="Current ephemeris")
    #parser.add_argument('-n', '--num', type=int, nargs=1, help="Value of factor n")

    args = parser.parse_args()
    return args    
    
    
def main(args):    
    obs = args.tim[0]
    eph = args.ephemeris[0]
    #num = args.num[0]

    ncy_check(obs)
    subprocess.call('tempo2 -gr plk -f %s %s -npsr 1 -nobs 60000 -showchisq'%(eph, obs), shell=True)
    trim(obs,eph,3)
    #subprocess.call('tempo2 -gr plk -f %s %s_trimtim -npsr 1 -nobs 60000 -showchisq'%(eph, obs), shell=True)
    gauss_fit(obs, eph, 3)
    subprocess.call('tempo2 -gr plk -f %s %s_trimtim_gauss -npsr 1 -nobs 60000 -showchisq'%(eph, obs), shell=True)    
    showeffect()

if __name__=="__main__":
    args = parse_arguments()
    main(args)