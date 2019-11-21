#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 12:28:10 2019
Author: Jun Wang
E-mail: jun.wang.ucas@gmail.com

Main Function: Show the plot of subints and channels vs bins.
               And fould out the reason that caused outliers.
               Not work if un-scrunched file in not in the same direcoty or have postfix other than .p .med .ar or .Tp
    
"""

import pandas as pd
import subprocess
import argparse, os


def show_arch(orig, prod):
    df_orig = pd.read_csv(orig,skiprows=1, dtype=str, header=None, delim_whitespace=True)
    
    df_prod = pd.read_csv(prod,skiprows=1, dtype=str, header=None, delim_whitespace=True)
    
    df_orig = df_orig.append(df_prod)
    df_orig = df_orig.append(df_prod)
    df_orig = df_orig.drop_duplicates(keep=False)

    for psr in df_orig[0]:
        if prod.rsplit('_', 1)[-1] == "trimtim":
            code_flag = "trimming process"
        elif prod.rsplit('_', 1)[-1] == "gauss":
            code_flag = "gaussian fitting process"
        
        print("The acrhive that was removed by %s code: **%s** \nq"%(code_flag, psr))


        if os.path.exists(psr.rsplit('.', 1)[0] + '.ar'):
            psr_raw = psr.rsplit('.', 1)[0] + '.ar'
        elif os.path.exists(psr.rsplit('.', 1)[0] + '.med'):
            psr_raw = psr.rsplit('.', 1)[0] + '.med'
        elif os.path.exists(psr.rsplit('.', 1)[0] + '.p'):
            psr_raw = psr.rsplit('.', 1)[0] + '.p'
        elif os.path.exists(psr.rsplit('.', 1)[0] + '.Tp'):
            psr_raw = psr.rsplit('.', 1)[0] + '.Tp'
        
        subprocess.call("pazi %s"%(psr_raw), shell=True)
        
        ar_rea = raw_input("Please input the problem of this archive: ")
        
        print(ar_rea, type(ar_rea))
        
        with open("reason.txt","a+") as f:
            f.write('The acrhive that was removed by %s code: %s \n'%(code_flag, psr))
            f.write('%s \n'%(ar_rea))



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Outlier Rejection for EPTA pulsar timing")
    parser.add_argument('-t', '--tim', type=str, nargs=1, help="Name of the TOA file")

    args = parser.parse_args()
    
    orig = args.tim[0]
    if orig.rsplit('.', 1)[-1] == "tim":
        prod = args.tim[0] + '_trimtim'
    elif orig.rsplit('.', 1)[-1] == "tim_trimtim":
        prod = args.tim[0] + '_gauss'
    else:
        print("Error")
        exit()
        
    
    show_arch(orig, prod)