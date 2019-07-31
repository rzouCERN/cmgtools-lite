#!/bin/bash
#arg1 = 6 for 2016, 7 for 2017, etc.
#arg2 = ttH, TTZ, TTW, etc.
   python ttH-multilepton/ttH_plots.py no 2lss_diff_201$1 -p $2  
