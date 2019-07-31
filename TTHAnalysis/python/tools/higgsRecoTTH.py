from CMGTools.TTHAnalysis.treeReAnalyzer import Collection
import ROOT, itertools

class HiggsRecoTTH:
    def __init__(self,label="_Recl",cut_BDT_rTT_score = 0.0, cuts_mW_had = (50.,110.), cuts_mH_vis = (90.,130.), btagDeepCSVveto = 0.4941):
        self.label = label
        self.branches = ["Hreco_%s"%x for x in ["minDRlj","visHmass","Wmass","lepIdx","j1Idx","j2Idx","dilepinvmass","pTHvis","nHiggs","nHiggs_status_prompt","nHiggs_status_isHP","nHiggs_status_fromHP","pTHgen"]]
        self.cut_BDT_rTT_score = cut_BDT_rTT_score
        self.cuts_mW_had = cuts_mW_had
        self.cuts_mH_vis = cuts_mH_vis
        self.btagDeepCSVveto = btagDeepCSVveto
       
    def listBranches(self):
        return self.branches
    def __call__(self,event):

        score = 1

        candidates=[]
	statusFlagsMap = {
	  'isPrompt'                           : 0,  # any decay product NOT coming from hadron, muon or tau decay
	  'isDecayedLeptonHadron'              : 1,  # a particle coming from hadron, muon, or tau decay
						     # (does not include resonance decays like W,Z,Higgs,top,etc)
						     # equivalent to status 2 in the current HepMC standard
	  'isTauDecayProduct'                  : 2,  # a direct or indirect tau decay product
	  'isPromptTauDecayProduct'            : 3,  # a direct or indirect decay product of a prompt tau
	  'isDirectTauDecayProduct'            : 4,  # a direct tau decay product
	  'isDirectPromptTauDecayProduct'      : 5,  # a direct decay product from a prompt tau
	  'isDirectHadronDecayProduct'         : 6,  # a direct decay product from a hadron
	  'isHardProcess'                      : 7,  # part of the hard process
	  'fromHardProcess'                    : 8,  # the direct descendant of a hard process particle of the same pdg id
	  'isHardProcessTauDecayProduct'       : 9,  # a direct or indirect decay product of a tau from the hard process
	  'isDirectHardProcessTauDecayProduct' : 10, # a direct decay product of a tau from the hard process
	  'fromHardProcessBeforeFSR'           : 11, # the direct descendant of a hard process particle of the same pdg id
					     # for outgoing particles the kinematics are those before QCD or QED FSR
	  'isFirstCopy'                        : 12, # the first copy of the particle in the chain with the same pdg id
	  'isLastCopy'                         : 13, # the last copy of the particle in the chain with the same pdg id
						     # (and therefore is more likely, but not guaranteed,
						     # to carry the final physical momentum)
	  'isLastCopyBeforeFSR'                : 14, # the last copy of the particle in the chain with the same pdg id
						     # before QED or QCD FSR (and therefore is more likely,
						     # but not guaranteed, to carry the momentum after ISR;
						     # only really makes sense for outgoing particles
	}
 
        genpar = Collection(event, "GenPart", "nGenPart")
        nHiggs = 0
        nHiggs_status = [0]*257
        pTHgen = 0
        for part in genpar:
            #if not part.pdgId==25 : continue
            if part.pdgId == 25:
               nHiggs +=1
               if not part.statusFlags &(1 << statusFlagsMap['isHardProcess']): continue  
               #print ("it is a hard Higgs")
               nHiggs_status[part.statusFlags &(1 << statusFlagsMap['isHardProcess'])] +=1
            #nHiggs_status[part.statusFlags &(1 << statusFlagsMap['isPrompt'])] +=1
            #nHiggs_status[part.statusFlags &(1 << statusFlagsMap['fromHardProcess'])] +=1
               pTHgen=part.p4().Pt()
               #print (part.pdgId) 
            

        if score>self.cut_BDT_rTT_score:

            nFO = getattr(event,"nLepFO"+self.label)
            lepsFO = Collection(event,"LepFO", "nLepFO")
            jets = [x for x in Collection(event,"JetSel"+self.label,"nJetSel"+self.label)]

            # Top-tagger to be used with new ntuples
            #j1top = getattr(event,"BDT_resolvedTopTagger_j1") 
            #j2top = getattr(event,"BDT_resolvedTopTagger_j2")
            #j3top = getattr(event,"BDT_resolvedTopTagger_j3")
            #jetsNoTopNoB = [j for i,j in enumerate(jets) if i not in [j1top,j2top,j3top] and j.btagDeepCSV<self.btagDeepCSVveto]

            jetsNoTopNoB = [j for i,j in enumerate(jets) if j.btagDeepB<self.btagDeepCSVveto]
	    for _lep,lep in [(ix,x.p4()) for ix,x in enumerate(lepsFO)]:
	       if not len(lepsFO)==2 : continue
	       if not lepsFO[0].charge==lepsFO[1].charge: continue 
      	       for _j1,_j2,j1,j2 in [(jets.index(x1),jets.index(x2),x1.p4(),x2.p4()) for x1,x2 in itertools.combinations(jetsNoTopNoB,2)]:
		   W = j1+j2
	           mW = W.M()
	     	   if mW<self.cuts_mW_had[0] or mW>self.cuts_mW_had[1]: continue
	       	   Wconstr = ROOT.TLorentzVector()
	   	   Wconstr.SetPtEtaPhiM(W.Pt(),W.Eta(),W.Phi(),80.4)
	 	   Hvisconstr = lep+Wconstr
	  	   mHvisconstr = Hvisconstr.M()
	      	   pTHvisconstr= Hvisconstr.Pt()
	  	   if mHvisconstr<self.cuts_mH_vis[0] or mHvisconstr>self.cuts_mH_vis[1]: continue
	     	   mindR = min(lep.DeltaR(j1),lep.DeltaR(j2))
		   candidates.append((mindR,mHvisconstr,mW,_lep,_j1,_j2,pTHvisconstr))
                                                             

        best = min(candidates) if len(candidates) else None

        ret = {}
        ret["Hreco_minDRlj"] = best[0] if best else -99
        ret["Hreco_visHmass"] = best[1] if best else -99
        ret["Hreco_Wmass"] = best[2] if best else -99
        ret["Hreco_lepIdx"] = best[3] if best else -99
        ret["Hreco_j1Idx"] = best[4] if best else -99
        ret["Hreco_j2Idx"] = best[5] if best else -99
        ret["Hreco_dilepinvmass"]= (lepsFO[0].p4()+lepsFO[best[3]].p4()).M() if best else -99 
        ret["Hreco_pTHvis"]=best[6] if best else -99
        ret["Hreco_nHiggs"]= nHiggs
        ret["Hreco_nHiggs_status_prompt"]= nHiggs_status[1]
        ret["Hreco_nHiggs_status_isHP"]= nHiggs_status[128]
        ret["Hreco_nHiggs_status_fromHP"]= nHiggs_status[256]
        ret["Hreco_pTHgen"]= pTHgen 
        return ret

            
MODULES=[]


MODULES.append( ('higgsRecoTTH', lambda : HiggsRecoTTH(label="") ) )
