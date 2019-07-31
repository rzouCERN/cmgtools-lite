## results presented on 19th July ##
## ------------------------------ ##
#for i in TTH_amc_2018 TTWToLNu_2018 TTZToLLNuNu_2018 TTTo2L2Nu_2018_part1 TTTo2L2Nu_2018_part2 TTTo2L2Nu_2018_part3 TTToSemiLeptonic_2018 TTJets_amc_SingleLeptFromTbar_2018
#do
#python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
#done

## extending to full mca after presenting on 19th July ##
## --------------------------------------------------- ##
#for i in THQ_2017 THW TTWW_2018 WZTo3LNu_2018 WZTo3LNu_ext_2018 W1JetsToLNu_2018 W2JetsToLNu_2018 W3JetsToLNu_2018 W4JetsToLNu_2018 ZZTo4L_2018+ZZTo4L_ext_2018_part1 ZZTo4L_ext_2018_part2 ZZTo4L_ext_2018_part3 ZZTo4L_ext_2018_part4 DYJetsToLL_M10to50_2018 DYJetsToLL_M50_amc_2018 GluGluToContinToZZTo2e2mu_2018 GluGluToContinToZZTo2e2tau_2018 GluGluToContinToZZTo2mu2tau_2018 GluGluToContinToZZTo4e_2018 GluGluToContinToZZTo4mu_2018 GluGluToContinToZZTo4tau_2018 GGHZZ4L VHToNonbb_2018 TTTT_2018 TTTW_2018 TTWH_2018 tZq_ll_2018_part1 tZq_ll_2018_part2 WW_DPS_2018 WpWpJJ_2017 WWW_2018 WWZ_2018 WZG_2018 WZZ_2018 ZZZ_2018 WJetsToLNu_2018 T_tch_2018 Tbar_tch_2018 T_sch_2018 TW_2018 TbarW_2018 WWTo2L2Nu_2018 TbarW_2018 TW_2018 WWTo2L2Nu_2018 TTGJets_2018 WGToLNuG_2018 ZGTo2LG_2018 TGJets_lep
#do
#python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
#done

## the below samples were missing from the above extension to all mca ##
## ------------------------------------------------------------------ ##

for i in ZZTo4L_2018 ZZTo4L_ext_2018_part1 ZZTo4L_ext_2018_part2 ZZTo4L_ext_2018_part3 ZZTo4L_ext_2018_part4 DYJetsToLL_M50_LO_2018_part1 DYJetsToLL_M50_LO_2018_part2
do
python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
done
