#### All events
python ttH-multilepton/ttH_plots.py no 2lss_SR_extr_synch_2016 syncTree_{cname}.root -p ttH_nomatch -U alwaystrue | bash

echo 2lss_SR_extr_synch_2017; python mcNtuple.py  -P ../../macros/synch_2017 --Fs {P}/2_synch -f -j 88 -l 41.4 --s2v -L ttH-multilepton/functionsTTH.cc --tree nanoAODskim --mcc ttH-multilepton/lepchoice-ttH-FO.txt --split-factor=-1   mca-synch.txt ttH-multilepton/2lss_tight.txt  ttH-multilepton/synchTuple.txt --xP '^lep(3|4)_.*' --xP '^(3|4)lep_.*' --xP 'kinMVA_3l_.*'  --binname 2lss    syncTree_2017_{cname}.root -U alwaystrue


rootcp --recreate syncTree_2016_TTH.root:t synchTuples_2016_v3.root:syncTree
rootcp --recreate syncTree_2017_TTH.root:t synchTuples_2017_v3.root:syncTree
rootcp --recreate syncTree_2018_TTH.root:t synchTuples_2018_v3.root:syncTree
/nfs/user/pvischia/tth/


echo 2lss_SR_extr_synch_2016; python mcNtuple.py  -P /nfs/user/pvischia/tth/synch_2016 --Fs {P}/2_synch -f -j 88 -l 41.4 --s2v -L ttH-multilepton/functionsTTH.cc --tree nanoAODskim --mcc ttH-multilepton/lepchoice-ttH-FO.txt --split-factor=-1   mca-synch.txt ttH-multilepton/2lss_tight.txt  ttH-multilepton/synchTuple.txt --xP '^lep(3|4)_.*' --xP '^(3|4)lep_.*' --xP 'kinMVA_3l_.*'  --binname 2lss    syncTree_2017_{cname}.root -U alwaystrue

#### 2lss

# SR
python ttH-multilepton/ttH_plots.py no 2lss_SR_extr_synch_2016 syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_2lSS_SR
#note that the year is missing in some of the commands, i.e. for the just below command, you need to add 2016 after synch_ so you get it running.
# Fake (removing TT cut)
python ttH-multilepton/ttH_plots.py no 2lss_SR_extr_synch_appl syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_2lSS_Fake

# Fake (removing TT cut)
python ttH-multilepton/ttH_plots.py no 2lss_SR_extr_synch_flip syncTree_{cname}.root -p ttH -E mm -I mm | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_2lSS_Flip


#### 3l

# SR
python ttH-multilepton/ttH_plots.py no 3l_extr_SR_synch syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_3l_SR


# Fake (removing TTT cut)
python ttH-multilepton/ttH_plots.py no 3l_extr_SR_synch_appl syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_3l_Fake

#### WZ
# CR
python ttH-multilepton/ttH_plots.py no cr_wz_extr_synch syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_WZctrl_SR
# Fake (removing TTT cut)
python ttH-multilepton/ttH_plots.py no cr_wz_extr_synch_appl syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_WZctrl_Fake

#### TTZ
# CR
python ttH-multilepton/ttH_plots.py no cr_ttz_extr_synch syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_ttZctrl_SR
# Fake (removing TTT cut)
python ttH-multilepton/ttH_plots.py no cr_ttz_extr_synch_appl syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_ttZctrl_Fake


#### TTW
# CR
python ttH-multilepton/ttH_plots.py no cr_3j_extr_synch syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_ttWctrl_SR
# Fake (removing TT cut)
python ttH-multilepton/ttH_plots.py no cr_3j_extr_synch_appl syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_ttWctrl_Fake

# Flip (inverting same-sign cut)
python ttH-multilepton/ttH_plots.py no cr_3j_extr_synch_flip syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_ttwctrl_Flip


### 4l 
python ttH-multilepton/ttH_plots.py no 4l_extr_synch syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_4l_SR


python ttH-multilepton/ttH_plots.py no 4l_extr_synch_relax syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_4l_Fake


#### ZZ 
python ttH-multilepton/ttH_plots.py no cr_zz_extr_synch syncTree_{cname}.root -p ttH | bash
rootcp syncTree_TTHnobb_pow.root:t synchTuples.root:syncTree_ZZctrl_SR
