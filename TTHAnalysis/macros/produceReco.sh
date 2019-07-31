## results presented on 19th July ## 
## ------------------------------ ##
#for i in TTHnobb_2016 TTWToLNu_2016 TTZToLLNuNu_ext_2016 TTZToLLNuNu_ext3_2016 TTZToLLNuNu_m1to10_2016 TTJets_SingleLeptonFromT_2016 TTJets_SingleLeptonFromTbar_2016 TTJets_DiLepton_2016 TTJets_DiLepton_ext_2016
#do
#python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
#done

## extending to full mca  after presenting on 19th July ## 
## ---------------------------------------------------- ##
for i in TTH THQ_2016 THW_2016 TTWW_2016 WZTo3LNu_2016 WJetsToLNu_LO_2016 WJetsToLNu_LO_ext1_2016 ZZTo4L_2016 GluGluToContinToZZTo2e2mu_2016 GluGluToContinToZZTo2e2tau_2016 GluGluToContinToZZTo2mu2tau_2016 GluGluToContinToZZTo4e_2016 GluGluToContinToZZTo4mu_2016 GluGluToContinToZZTo4tau_2016 DYJetsToLL_M10to50_LO_2016 DYJetsToLL_M50_LO_ext2_2016_part1 DYJetsToLL_M50_LO_ext2_2016_part2 DYJetsToLL_M50_LO_ext_2016 WpWpJJ_2016 WWDouble WWW_2016 WWZ_2016 WZZ_2016 ZZZ_2016 TTTT_2016 tZq_ll_2016 tWll T_tch_2016 Tbar_tch_2016 T_sch_2016 TW_2016 TbarW_2016 WWTo2L2Nu_2016 TTGJets_2016 TTGJets_ext1_2016 WGToLNuG_amcatanlo_ext1_2016 WGToLNuG_amcatanlo_ext2_2016 WGToLNuG_amcatanlo_ext3_2016 ZGTo2LG_2016 TGJets TGJets_ext_skimPromptGamma
do
python prepareEventVariablesFriendTree.py -I CMGTools.TTHAnalysis.tools.higgsRecoTTH -m higgsRecoTTH -t nanoAODskim -F sf/t {P}1_triggerDecision/evVarFriend_{cname}.root -F sf/t {P}2_countJets/evVarFriend_{cname}.root -F sf/t {P}2_countJets_index/evVarFriend_{cname}.root -F sf/t {P}3_lepVars/evVarFriend_{cname}.root -F sf/t {P}4_bweight_new/evVarFriend_{cname}.root -F sf/t {P}5_MET/evVarFriend_{cname}.root -F sf/t {P}5_taus/evVarFriend_{cname}.root -F sf/t {P}6_eventVars/evVarFriend_{cname}.root -d $i /nfs/user/pvischia/tth/ttH_EasterProduction_2lss_3l/ . -j 4
done




























