#!/bin/bash
#arg1 = 6 for 2016, 7 for 2017, etc.
#arg2 = ttH, TTZ, TTW, etc.
   python ttH-multilepton/ttH_plots.py no 2lss_diff_201$1 -p $2 
#./diff_plot.sh 6 ttH,TTWW,TTW,TTZ,EWK,Rares,TT,DY,WJets,SingleTop,WW,Flips,Convs
#./diff_plot.sh 7 ttH,tHW,TTWW,TTW,TTZ,WZ,ZZ,Rares,tttt,tZq,TT,DY,WJets,SingleTop,WW,Flips,Convs
#./diff_plot.sh 8 ttH,TTWW,TTW,TTZ,WZ,ZZ,Rares,tttt,tZq,TT,DY,WJets,SingleTop,WW,Flips,Convs 
