## results presented on 19th July ##
## ------------------------------ ##
#for i in TTH_amc_2017 TTW_LO_2017 TTW_LO_ext1_2017 TTZToLLNuNu_2017 TTZToLLNuNu_m1to10_2017 TTTo2L2Nu_2017_part1 TTTo2L2Nu_2017_part2 TTTo2L2Nu_2017_part3 TTToSemiLeptonic_2017 TTToSemiLeptonic_PSweights_2017 
#do 
#python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
#done

## below is ~66 files
## extending to full mca after presenting on 19th July ##
## --------------------------------------------------- ##
#for i in THQ_2017 THW_2017 TTWW_2017 WZTo3LNu_2017 W1JetsToLNu_2017 W2JetsToLNu_2017 W3JetsToLNu_2017 W4JetsToLNu_2017 ZZTo4L_2017 ZZTo4L_ext1_2017_part1 ZZTo4L_ext1_2017_part2 ZZTo4L_ext1_2017_part3 DYJetsToLL_M10to50_2017 DYJetsToLL_M50_2017_part1 DYJetsToLL_M50_2017_part2 DYJetsToLL_M50_2017_part3 DYJetsToLL_M50_2017_part4 DYJetsToLL_M50_2017_part5 GluGluToContinToZZTo2e2mu_2017 GluGluToContinToZZTo2e2mu_ext1_2017 GluGluToContinToZZTo2e2tau_2017 GluGluToContinToZZTo2e2tau_ext1_2017 GluGluToContinToZZTo2mu2tau_2017 GluGluToContinToZZTo2mu2tau_ext1_2017 GluGluToContinToZZTo4e_2017 GluGluToContinToZZTo4e_ext1_2017 GluGluToContinToZZTo4mu_2017 GluGluToContinToZZTo4mu_ext1_2017 GluGluToContinToZZTo4tau_2017 GGHZZ4L_2017 VHToNonbb_2017 TTTT_2017 TTTW_2017 TTWH_2017 tZq_ll_2017 WW_DPS_2017 WpWpJJ_2017 WWW_2017 WWZ_2017 WZG_2017 WZZ_2017 ZZZ_2017 WJetsToLNu_2017 T_tch_2017_2017 Tbar_tch_2017_2017 T_sch_2017_2017 TW_NoFullyHadronic_2017 TbarW_NoFullyHadronic_2017 WWTo2L2Nu_2017 TTGJets_2017 WGToLNuG_amcatanlo_ext WGToLNuG_amcatnlo_ext2_skimPromptGamma ZGTo2LG_ext W1JetsToLNu_2017 W2JetsToLNu_2017 W3JetsToLNu_2017 W4JetsToLNu_2017 DYJetsToLL_M50_0J_part1 DYJetsToLL_M50_0J_part2 DYJetsToLL_M50_1J_2017_part1 DYJetsToLL_M50_1J_2017_part2 DYJetsToLL_M50_1J_2017_part3 DYJetsToLL_M50_2J_2017_part1 DYJetsToLL_M50_2J_2017_part2 DYJetsToLL_M50_2J_2017_part3
#do 
#python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
#done

## the below samples were missing from the above extension to all mca ##
## ------------------------------------------------------------------ ##

for i in TTW_LO_2017 TTW_LO_ext1_2017 DYJetsToLL_M4to50_HT100to200_2017 DYJetsToLL_M4to50_HT100to200_ext1_2017 DYJetsToLL_M4to50_HT200to400_2017 DYJetsToLL_M4to50_HT200to400_ext1_2017 DYJetsToLL_M4to50_HT400to600_2017 DYJetsToLL_M4to50_HT400to600_ext1_2017 DYJetsToLL_M4to50_HT600toInf_2017 DYJetsToLL_M4to50_HT600toInf_ext1_2017 TGJets_lep_2017 
do 
python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
done
