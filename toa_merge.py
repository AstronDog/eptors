#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 10:31:33 2019
Author: Jun Wang
E-mail: jun.wang.ucas@gmail.com

Main Function: This code is used to update the origin TOAs with newly created TOAs after re-analysis of the outlier.
    
"""

import pandas as pd
import os, argparse

def df_merge(ori_tim, upd_tim):

    origin_toa = pd.read_csv(ori_tim, skiprows=1, dtype=str, header=None, delim_whitespace=True)

    update_toa = pd.read_csv(upd_tim, skiprows=1, dtype=str, header=None, delim_whitespace=True) 

   

    origin_toa( origin_toa[[0]].merge(update_toa, 'left'))
    
    if args.output == '':
        orig_name = ori_tim.split(':', 1)[1].strip()
        o_name = orig_name + '_new.tim'
    else:
        o_name = args.output
    
    
    
    
    origin_toa.to_csv(o_name, sep=' ', header=False, index=False)
    os.system("sed -i '1iFORMAT 1' test.new.tim")




def parse_arguments():    
    parser = argparse.ArgumentParser(description="Outlier Rejection for EPTA pulsar timing")
    parser.add_argument('-t1', '--oritim', type=str, nargs=1, help="name of the origin TOA file")
    parser.add_argument('-t2', '--updtim', type=str, nargs=1, help="name of the update TOA file")
    parser.add_argument('-o', '--output', type=str, nargs=1, default='', help="Name of the output file.")
    
    
    args = parser.parse_args()
    return args    
    
    
def main(args):    
    ori_tim = args.oritim[0]
    upd_tim = args.updtim[0]
    
    df_merge(ori_tim, upd_tim)

if __name__=="__main__":
    args = parse_arguments()
    main(args)