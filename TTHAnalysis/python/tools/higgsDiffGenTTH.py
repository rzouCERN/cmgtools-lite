from __future__ import division
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from CMGTools.TTHAnalysis.tools.nanoAOD.friendVariableProducerTools import declareOutput, writeOutput
from CMGTools.TTHAnalysis.treeReAnalyzer import Collection as CMGCollection
import ROOT, itertools
from math import *
import sys
import CMGTools.TTHAnalysis.tools.higgsDiffUtils as diffUtils 
from PhysicsTools.Heppy.physicsobjects.Jet import _btagWPs as HiggsRecoTTHbtagwps

class HiggsDiffGenTTH(Module):
    def __init__(self, label = 'Hreco_', debug=False):
        self.debug = debug
        self.label = label
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Particle counters
        self.out.branch('%snHiggses'%self.label         , 'I')
        self.out.branch('%snTfromhardprocess'%self.label, 'I')
        self.out.branch('%snHadT'%self.label            , 'I')
        self.out.branch('%snLepT'%self.label            , 'I')
        self.out.branch('%snWFromHadT'%self.label       , 'I')
        self.out.branch('%snWFromLepT'%self.label       , 'I')
        self.out.branch('%snWFromH'%self.label          , 'I')
        self.out.branch('%snHadWFromH'%self.label       , 'I')
        self.out.branch('%snWFromT'%self.label          , 'I')
        self.out.branch('%snQFromW'%self.label          , 'I')
        self.out.branch('%snGenLep'%self.label          , 'I')
        self.out.branch('%snLFromW'%self.label          , 'I')
        self.out.branch('%snNuFromWFromH'%self.label    , 'I')
        self.out.branch('%snNuFromWFromT'%self.label    , 'I')
        self.out.branch('%snQFromWFromH'%self.label     , 'I')
        self.out.branch('%snQFromWFromT'%self.label     , 'I')
        self.out.branch('%snLFromWFromH'%self.label     , 'I')
        self.out.branch('%snLFromWFromT'%self.label     , 'I')
        self.out.branch('%snTauFromH'%self.label        , 'I')
        self.out.branch('%snhFromTauFromH'%self.label   , 'I')
        self.out.branch('%snLFromTauFromH'%self.label   , 'I')
        self.out.branch('%snNuFromTauFromH'%self.label  , 'I')
        self.out.branch('%snTNuFromTauFromH'%self.label , 'I')
        self.out.branch('%snZFromH'%self.label          , 'I')
        self.out.branch('%snQFromZFromH'%self.label     , 'I')
        self.out.branch('%snLFromZFromH'%self.label     , 'I')
        self.out.branch('%snNuFromZFromH'%self.label    , 'I')
        self.out.branch('%snTauFromWFromH'%self.label   , 'I')
        self.out.branch('%snLFromTauFromWFromH'%self.label  , 'I')
        self.out.branch('%snNuFromTauFromWFromH'%self.label , 'I')
        self.out.branch('%snhFromTauFromWFromH'%self.label  , 'I')

        # Although we expect a maximum of 2 objects per array, we allow for 4 of them to be stored, for safety and later checks
        for suffix in ["_pt", "_eta", "_phi", "_mass"]:
            
            self.out.branch('%sHiggses%s'%(self.label,suffix)          , 'F', 4)
            self.out.branch('%sTfromhardprocess%s'%(self.label,suffix) , 'F', 4)
            self.out.branch('%sHadT%s'%(self.label,suffix)             , 'F', 4)
            self.out.branch('%sLepT%s'%(self.label,suffix)             , 'F', 4)
            self.out.branch('%sWFromHadT%s'%(self.label,suffix)        , 'F', 4)
            self.out.branch('%sWFromLepT%s'%(self.label,suffix)        , 'F', 4)
            self.out.branch('%sWFromH%s'%(self.label,suffix)           , 'F', 4)
            self.out.branch('%sHadWFromH%s'%(self.label,suffix)        , 'F', 4)
            self.out.branch('%sWFromT%s'%(self.label,suffix)           , 'F', 4)
            self.out.branch('%sQFromW%s'%(self.label,suffix)           , 'F', 4)
            self.out.branch('%sGenLep%s'%(self.label,suffix)           , 'F', 4)
            self.out.branch('%sLFromW%s'%(self.label,suffix)           , 'F', 4)
            self.out.branch('%sNuFromWFromH%s'%(self.label,suffix)     , 'F', 4)
            self.out.branch('%sNuFromWFromT%s'%(self.label,suffix)     , 'F', 4)
            self.out.branch('%sQFromWFromH%s'%(self.label,suffix)      , 'F', 4)
            self.out.branch('%sQFromWFromT%s'%(self.label,suffix)      , 'F', 4)
            self.out.branch('%sLFromWFromH%s'%(self.label,suffix)      , 'F', 4)
            self.out.branch('%sLFromWFromT%s'%(self.label,suffix)      , 'F', 4)
            self.out.branch('%sTauFromH%s'%(self.label,suffix)         , 'F', 4)
            self.out.branch('%shFromTauFromH%s'%(self.label,suffix)    , 'F', 4)
            self.out.branch('%sLFromTauFromH%s'%(self.label,suffix)    , 'F', 4)
            self.out.branch('%sNuFromTauFromH%s'%(self.label,suffix)   , 'F', 4)
            self.out.branch('%sTNuFromTauFromH%s'%(self.label,suffix)  , 'F', 4)
            self.out.branch('%sZFromH%s'%(self.label,suffix)           , 'F', 4)
            self.out.branch('%sQFromZFromH%s'%(self.label,suffix)      , 'F', 4)
            self.out.branch('%sLFromZFromH%s'%(self.label,suffix)      , 'F', 4)
            self.out.branch('%sNuFromZFromH%s'%(self.label,suffix)     , 'F', 4)
            self.out.branch('%sTauFromWFromH%s'%(self.label,suffix)    , 'F', 4)
            self.out.branch('%sLFromTauFromWFromH%s'%(self.label,suffix)   , 'F', 4)
            self.out.branch('%sNuFromTauFromWFromH%s'%(self.label,suffix)  , 'F', 4)
            self.out.branch('%shFromTauFromWFromH%s'%(self.label,suffix)   , 'F', 4)

        # Some precomputed quantities of interest

        self.out.branch('%spTHgen'%self.label            , 'F')
        self.out.branch('%spTtgen'%self.label            , 'F', 4)
        self.out.branch('%sinv_mass_q1_q2'%self.label    , 'F')
        self.out.branch('%sdelR_partonsFromH'%self.label , 'F')
        self.out.branch('%spt_partonsFromH'%self.label   , 'F')
        self.out.branch('%squark1pT'%self.label          , 'F')
        self.out.branch('%squark2pT'%self.label          , 'F')
        self.out.branch('%squark1Eta'%self.label         , 'F')
        self.out.branch('%squark2Eta'%self.label         , 'F')
        self.out.branch('%squark1Flavour'%self.label     , 'F')
        self.out.branch('%squark2Flavour'%self.label     , 'F')
        self.out.branch('%sdelR_H_q1l'%self.label        , 'F')
        self.out.branch('%sdelR_H_q2l'%self.label        , 'F')
        self.out.branch('%spTTrueGen'%self.label         , 'F')
        self.out.branch('%sMTrueGen'%self.label          , 'F')
        self.out.branch('%spTTrueGenPlusNu'%self.label   , 'F')
        
    def analyze(self, event):

        # Input collections and maps
        genpar = Collection(event,"GenPart","nGenPart") 
        
        Higgses          = []
        Tfromhardprocess = []
        HadT             = []
        LepT             = []
        WFromHadT        = []
        WFromLepT        = []
        WFromH           = []
        HadWFromH        = []
        WFromT           = []
        QFromW           = []
        GenLep           = []
        LFromW           = []
        NuFromWFromH     = []
        NuFromWFromT     = []
        QFromWFromH      = []
        QFromWFromT      = []
        LFromWFromH      = []
        LFromWFromT      = []

        TauFromH         = []
        hFromTauFromH    = []
        LFromTauFromH    = []
        NuFromTauFromH   = []
        TNuFromTauFromH  = []

        ZFromH           = []
        QFromZFromH      = []
        LFromZFromH      = []
        NuFromZFromH     = []

        TauFromWFromH       = []
        LFromTauFromWFromH  = []
        NuFromTauFromWFromH = []
        hFromTauFromWFromH  = []

        for part in genpar:
            # higgs
            if part.pdgId == 25 and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess']):
                Higgses.append(part)

            # tops
            if abs(part.pdgId) == 6 and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess']):
                Tfromhardprocess.append(part)

            # W from higgs
            if (abs(part.pdgId) == 24 and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])
                and part.genPartIdxMother >= 0 and  genpar[part.genPartIdxMother].pdgId == 25 ):
                WFromH.append(part)
                if self.debug: print "it is a hard W coming from a Higgs"
                
            # W from tops 
            if (abs(part.pdgId) == 24 and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])
                and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 6):
                WFromT.append(part)
                if self.debug: print "it is a hard W coming from a top"

            # W decays to quarks
            if (abs(part.pdgId) in [1,2,3,4,5,6] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24): 
                QFromW.append(part)
                if self.debug: print "it is a quark coming from a W"

            # gen leptons
            if (abs(part.pdgId) in [11,13] and part.status == 1 and part.statusFlags &(1 << diffUtils.statusFlagsMap['isLastCopy']) 
                and not part.statusFlags &(1 << diffUtils.statusFlagsMap['isDirectHadronDecayProduct'])):
                if part.statusFlags &(1 << diffUtils.statusFlagsMap['isPrompt']) or part.statusFlags &(1 << diffUtils.statusFlagsMap['isDirectPromptTauDecayProduct']):
                    GenLep.append(part)
                    if self.debug: print "it is a prompt lepton"

            # gen leptons from W
            if (abs(part.pdgId) in [11,13] and part.status == 1 and part.statusFlags &(1 << diffUtils.statusFlagsMap['isLastCopy']) and not part.statusFlags &(1 << diffUtils.statusFlagsMap['isDirectHadronDecayProduct'])):
                if part.statusFlags &(1 << diffUtils.statusFlagsMap['isPrompt']) or part.statusFlags &(1 << diffUtils.statusFlagsMap['isDirectPromptTauDecayProduct']):
                    if part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24: 
                        LFromW.append(part)
                        if self.debug: print "it is a prompt lepton"
                        
            # neutrino from W from H
            if abs(part.pdgId) in [12, 14]:
                if self.debug: print "it is a neutrino"
                if part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24: 
                    if self.debug: print "the mother of this neutrino is W+ or W-"
                    if abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId) == 25:
                        NuFromWFromH.append(part)
                        if self.debug: print "the mother of this W is a Higgs"
                        
            # neutrino from W from top
            if abs(part.pdgId) in [12, 14, 16]:
                if self.debug: print "it is a neutrino"
                if part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24:
                    if self.debug: print "the mother of this neutrino is W+ or W-"
                    mpart = genpar[part.genPartIdxMother]
                    while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 24): mpart = genpar[mpart.genPartIdxMother]
                    if abs(genpar[mpart.genPartIdxMother].pdgId) == 6:
                        NuFromWFromT.append(part)
                        if mpart not in WFromLepT: WFromLepT.append(mpart)
                        if genpar[mpart.genPartIdxMother] not in LepT: LepT.append(genpar[mpart.genPartIdxMother])
                        if self.debug: print "the mother of this W is a top"
        
            # quarks from W from H
            if (abs(part.pdgId) in [1,2,3,4,5,6] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24
                     and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                if self.debug: print "it is a quark coming from a hard W"
                if (genpar[genpar[part.genPartIdxMother].genPartIdxMother].genPartIdxMother >= 0 
                    and genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId == 25
                    and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                    QFromWFromH.append(part)
                    if genpar[part.genPartIdxMother] not in HadWFromH: HadWFromH.append(genpar[part.genPartIdxMother])
                    if self.debug: print "the mother of this hard W is a hard Higgs"

            # quarks from W from T 
            if (abs(part.pdgId) in [1,2,3,4,5,6] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24): 
                if self.debug: print "it is a quark coming from a hard W"
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 24): mpart = genpar[mpart.genPartIdxMother]
                if abs(genpar[mpart.genPartIdxMother].pdgId) == 6:
                    QFromWFromT.append(part)
                    if mpart not in WFromHadT: WFromHadT.append(mpart)
                    if genpar[mpart.genPartIdxMother] not in HadT: HadT.append(genpar[mpart.genPartIdxMother])
                    if self.debug: print "the mother of this hard W is a hard top"

            # leptons (excl. taus) from W from H 
            if (abs(part.pdgId) in [11,13] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24 
                     and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                if self.debug: print "it is a lepton coming from a hard W"
                if (genpar[genpar[part.genPartIdxMother].genPartIdxMother].genPartIdxMother >= 0 
                    and genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId == 25 
                    and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                    LFromWFromH.append(part)
                    if self.debug: print "the mother of this hard W is a hard Higgs"
                    
            # leptons (excl. taus) from W from top
            if (abs(part.pdgId) in [11,13] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24
                     and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                if self.debug: print "it is a lepton coming from a hard W"
                if (genpar[genpar[part.genPartIdxMother].genPartIdxMother].genPartIdxMother >= 0 
                    and abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId) == 6 
                    and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                    LFromWFromT.append(part)
                    if self.debug: print "the mother of this W is a hard top"
        
            # Tau from H
            if (abs(part.pdgId) in [15] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==25):
                TauFromH.append(part)
                if self.debug: print "it is a hard tau coming from a higgs"

            # pion / kaon from tau from H
            if (abs(part.pdgId) in [111,211,311,321] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==15):
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 15): mpart = genpar[mpart.genPartIdxMother]
                if (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 25):
                    hFromTauFromH.append(part)
                    if self.debug: print "it is a pion or kaon coming from a hard tau from a higgs"

            # Lepton from tau from H
            if (abs(part.pdgId) in [11,13] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==15):
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 15): mpart = genpar[mpart.genPartIdxMother]
                if (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 25):
                    LFromTauFromH.append(part)
                    if self.debug: print "it is a lepton coming from a hard tau from a higgs"

            # lnu from tau from H
            if (abs(part.pdgId) in [12,14] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==15):
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 15): mpart = genpar[mpart.genPartIdxMother]
                if (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 25):
                    NuFromTauFromH.append(part)
                    if self.debug: print "it is a lepton neutrino coming from a hard tau from a higgs"

            # tnu from tau from H
            if (abs(part.pdgId) in [16] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==15):
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 15): mpart = genpar[mpart.genPartIdxMother]
                if (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 25):
                    TNuFromTauFromH.append(part)
                    if self.debug: print "it is a tau neutrino coming from a hard tau from a higgs"

            # Z from H
            if (abs(part.pdgId) in [23] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==25):
                ZFromH.append(part)
                if self.debug: print "it is a hard Z coming from a higgs"

            # Q from Z from H
            if (abs(part.pdgId) in [1,2,3,4,5,6] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==23):
                if (genpar[part.genPartIdxMother].genPartIdxMother >= 0 and abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId)==25):
                    QFromZFromH.append(part)
                    if self.debug: print "it is a quark coming from a hard Z from a higgs"

            # L from Z from H
            if (abs(part.pdgId) in [11,13] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==23):
                if (genpar[part.genPartIdxMother].genPartIdxMother >= 0 and abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId)==25):
                    LFromZFromH.append(part)
                    if self.debug: print "it is a lepton coming from a hard Z from a higgs"

            # Nu from Z from H
            if (abs(part.pdgId) in [12,14] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==23):
                if (genpar[part.genPartIdxMother].genPartIdxMother >= 0 and abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId)==25):
                    NuFromZFromH.append(part)
                    if self.debug: print "it is a lepton neutrino coming from a hard Z from a higgs"

            # Tau from W from H
            if (abs(part.pdgId) in [15] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24):
                if (genpar[part.genPartIdxMother].genPartIdxMother >= 0 and abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId)==25):
                    TauFromWFromH.append(part)
                    if self.debug: print "it is a tau coming from a hard W from a higgs"

            # L from tau from W from H
            if (abs(part.pdgId) in [11,13] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 15):
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 15): mpart = genpar[mpart.genPartIdxMother]
                if (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 24
                        and genpar[mpart.genPartIdxMother].genPartIdxMother >= 0
                        and abs(genpar[genpar[mpart.genPartIdxMother].genPartIdxMother].pdgId) == 25):
                    LFromTauFromWFromH.append(part)
                    if self.debug: print "it is a lepton coming from a tau from a hard W from a higgs"

            # Nu from tau from W from H
            if (abs(part.pdgId) in [12,14] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 15):
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 15): mpart = genpar[mpart.genPartIdxMother]
                if (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 24
                        and genpar[mpart.genPartIdxMother].genPartIdxMother >= 0
                        and abs(genpar[genpar[mpart.genPartIdxMother].genPartIdxMother].pdgId) == 25):
                    NuFromTauFromWFromH.append(part)
                    if self.debug: print "it is a lepton neutrino coming from a tau from a hard W from a higgs"

            # pion / kaon from tau from W from H
            if (abs(part.pdgId) in [111,211,311,321] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 15):
                mpart = genpar[part.genPartIdxMother]
                while (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 15): mpart = genpar[mpart.genPartIdxMother]
                if (mpart.genPartIdxMother >= 0 and abs(genpar[mpart.genPartIdxMother].pdgId) == 24
                        and genpar[mpart.genPartIdxMother].genPartIdxMother >= 0
                        and abs(genpar[genpar[mpart.genPartIdxMother].genPartIdxMother].pdgId) == 25):
                    hFromTauFromWFromH.append(part)
                    if self.debug: print "it is a pion or kaon coming from a tau from a hard W from a higgs"

            if self.debug:
                # taus
                if (abs(part.pdgId) in [15] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24):
                    if(genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                        tauFromW.append(part)
                        if self.debug: print "it is a tau coming from a hard W"
                    
                    # leptons from top as recommended by gen particle producer
                if (abs(part.pdgId) in [11,13] 
                    and part.statusFlags &(1 << diffUtils.statusFlagsMap['isPrompt'])
                    and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])
                    and part.statusFlags &(1 << diffUtils.statusFlagsMap['isFirstCopy'])
                    and not part.statusFlags &(1 << diffUtils.statusFlagsMap['isDirectHadronDecayProduct'])):
                    if self.debug: print "it should be a  lepton coming from top"
                    # LFromWFromT.append(part) # why to the same collection? No.
        # End loop on genparticles

        if self.debug:
            print (" >> in this event: "  
                   + " \n higgses              = " + str(len(Higgses)) 
                   + " \n W from Higgs         = " + str(len(WFromH))
                   + " \n hard tops            = " + str(len(Tfromhardprocess))
                   + " \n W from tops          = " + str(len(WFromT))
                   + " \n leptons              = " + str(len(GenLep)) 
                   + " \n QFromWFromH          = " + str(len(QFromWFromH))
                   + " \n QFromWFromT          = " + str(len(QFromWFromT)) 
                   + " \n LFromWFromH          = " + str(len(LFromWFromH)) 
                   + " \n LFromWFromT          = " + str(len(LFromWFromT)) 
                   + " \n NuFromWFromH         = " + str(len(NuFromWFromH)) 
                   + " \n NuFromWFromT         = " + str(len(NuFromWFromT))
                   + " \n quarks from W's      = " + str(len(QFromW))
                   + " \n leptons from W's     = " + str(len(LFromW))
                   #   + " \n taus from W's        = " + str(len(tauFromW))
                   + " \n <<")


        # Fill branches for particles

        self.out.fillBranch('%snHiggses'%self.label         , len(Higgses))
        self.out.fillBranch('%snTfromhardprocess'%self.label, len(Tfromhardprocess))
        self.out.fillBranch('%snHadT'%self.label            , len(HadT))
        self.out.fillBranch('%snLepT'%self.label            , len(LepT))
        self.out.fillBranch('%snWFromHadT'%self.label       , len(WFromHadT))
        self.out.fillBranch('%snWFromLepT'%self.label       , len(WFromLepT))
        self.out.fillBranch('%snWFromH'%self.label          , len(WFromH))
        self.out.fillBranch('%snHadWFromH'%self.label       , len(HadWFromH))
        self.out.fillBranch('%snWFromT'%self.label          , len(WFromT))
        self.out.fillBranch('%snQFromW'%self.label          , len(QFromW))
        self.out.fillBranch('%snGenLep'%self.label          , len(GenLep))
        self.out.fillBranch('%snLFromW'%self.label          , len(LFromW))
        self.out.fillBranch('%snNuFromWFromH'%self.label    , len(NuFromWFromH))
        self.out.fillBranch('%snNuFromWFromT'%self.label    , len(NuFromWFromT))
        self.out.fillBranch('%snQFromWFromH'%self.label     , len(QFromWFromH))
        self.out.fillBranch('%snQFromWFromT'%self.label     , len(QFromWFromT))
        self.out.fillBranch('%snLFromWFromH'%self.label     , len(LFromWFromH))
        self.out.fillBranch('%snLFromWFromT'%self.label     , len(LFromWFromT)) 
        self.out.fillBranch('%snTauFromH'%self.label        , len(TauFromH))
        self.out.fillBranch('%snhFromTauFromH'%self.label   , len(hFromTauFromH))
        self.out.fillBranch('%snLFromTauFromH'%self.label   , len(LFromTauFromH))
        self.out.fillBranch('%snNuFromTauFromH'%self.label  , len(NuFromTauFromH))
        self.out.fillBranch('%snTNuFromTauFromH'%self.label , len(TNuFromTauFromH))
        self.out.fillBranch('%snZFromH'%self.label          , len(ZFromH))
        self.out.fillBranch('%snQFromZFromH'%self.label     , len(QFromZFromH))
        self.out.fillBranch('%snLFromZFromH'%self.label     , len(LFromZFromH))
        self.out.fillBranch('%snNuFromZFromH'%self.label    , len(NuFromZFromH))
        self.out.fillBranch('%snTauFromWFromH'%self.label   , len(TauFromWFromH))
        self.out.fillBranch('%snLFromTauFromWFromH'%self.label  , len(LFromTauFromWFromH))
        self.out.fillBranch('%snNuFromTauFromWFromH'%self.label , len(NuFromTauFromWFromH))
        self.out.fillBranch('%snhFromTauFromWFromH'%self.label  , len(hFromTauFromWFromH))

        self.out.fillBranch('%sHiggses_pt'%self.label               , [ Higgses[i].p4().Pt()               if len(Higgses) > i             else -99 for i in range(4) ]) 
        self.out.fillBranch('%sTfromhardprocess_pt'%self.label      , [ Tfromhardprocess[i].p4().Pt()      if len(Tfromhardprocess) > i    else -99 for i in range(4) ]) 
        self.out.fillBranch('%sHadT_pt'%self.label                  , [ HadT[i].p4().Pt()                  if len(HadT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sLepT_pt'%self.label                  , [ LepT[i].p4().Pt()                  if len(LepT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromHadT_pt'%self.label             , [ WFromHadT[i].p4().Pt()             if len(WFromHadT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromLepT_pt'%self.label             , [ WFromLepT[i].p4().Pt()             if len(WFromLepT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromH_pt'%self.label                , [ WFromH[i].p4().Pt()                if len(WFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sHadWFromH_pt'%self.label             , [ HadWFromH[i].p4().Pt()             if len(HadWFromH) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromT_pt'%self.label                , [ WFromT[i].p4().Pt()                if len(WFromT) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromW_pt'%self.label                , [ QFromW[i].p4().Pt()                if len(QFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sGenLep_pt'%self.label                , [ GenLep[i].p4().Pt()                if len(GenLep) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromW_pt'%self.label                , [ LFromW[i].p4().Pt()                if len(LFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromH_pt'%self.label          , [ NuFromWFromH[i].p4().Pt()          if len(NuFromWFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromT_pt'%self.label          , [ NuFromWFromT[i].p4().Pt()          if len(NuFromWFromT) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromH_pt'%self.label           , [ QFromWFromH[i].p4().Pt()           if len(QFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromT_pt'%self.label           , [ QFromWFromT[i].p4().Pt()           if len(QFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromH_pt'%self.label           , [ LFromWFromH[i].p4().Pt()           if len(LFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromT_pt'%self.label           , [ LFromWFromT[i].p4().Pt()           if len(LFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromH_pt'%self.label              , [ TauFromH[i].p4().Pt()              if len(TauFromH) > i            else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromH_pt'%self.label         , [ hFromTauFromH[i].p4().Pt()         if len(hFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromH_pt'%self.label         , [ LFromTauFromH[i].p4().Pt()         if len(LFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromH_pt'%self.label        , [ NuFromTauFromH[i].p4().Pt()        if len(NuFromTauFromH) > i      else -99 for i in range(4) ])
        self.out.fillBranch('%sTNuFromTauFromH_pt'%self.label       , [ TNuFromTauFromH[i].p4().Pt()       if len(TNuFromTauFromH) > i     else -99 for i in range(4) ])
        self.out.fillBranch('%sZFromH_pt'%self.label                , [ ZFromH[i].p4().Pt()                if len(ZFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromZFromH_pt'%self.label           , [ QFromZFromH[i].p4().Pt()           if len(QFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromZFromH_pt'%self.label           , [ LFromZFromH[i].p4().Pt()           if len(LFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromZFromH_pt'%self.label          , [ NuFromZFromH[i].p4().Pt()          if len(NuFromZFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromWFromH_pt'%self.label         , [ TauFromWFromH[i].p4().Pt()         if len(TauFromWFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromWFromH_pt'%self.label    , [ LFromTauFromWFromH[i].p4().Pt()    if len(LFromTauFromWFromH) > i  else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromWFromH_pt'%self.label   , [ NuFromTauFromWFromH[i].p4().Pt()   if len(NuFromTauFromWFromH) > i else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromWFromH_pt'%self.label    , [ hFromTauFromWFromH[i].p4().Pt()    if len(hFromTauFromWFromH) > i  else -99 for i in range(4) ])

        self.out.fillBranch('%sHiggses_eta'%self.label              , [ Higgses[i].p4().Eta()              if len(Higgses) > i             else -99 for i in range(4) ]) 
        self.out.fillBranch('%sTfromhardprocess_eta'%self.label     , [ Tfromhardprocess[i].p4().Eta()     if len(Tfromhardprocess) > i    else -99 for i in range(4) ]) 
        self.out.fillBranch('%sHadT_eta'%self.label                 , [ HadT[i].p4().Eta()              if len(HadT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sLepT_eta'%self.label                 , [ LepT[i].p4().Eta()                 if len(LepT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromHadT_eta'%self.label            , [ WFromHadT[i].p4().Eta()            if len(WFromHadT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromLepT_eta'%self.label            , [ WFromLepT[i].p4().Eta()            if len(WFromLepT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromH_eta'%self.label               , [ WFromH[i].p4().Eta()               if len(WFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sHadWFromH_eta'%self.label            , [ HadWFromH[i].p4().Eta()            if len(HadWFromH) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromT_eta'%self.label               , [ WFromT[i].p4().Eta()               if len(WFromT) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromW_eta'%self.label               , [ QFromW[i].p4().Eta()               if len(QFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sGenLep_eta'%self.label               , [ GenLep[i].p4().Eta()               if len(GenLep) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromW_eta'%self.label               , [ LFromW[i].p4().Eta()               if len(LFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromH_eta'%self.label         , [ NuFromWFromH[i].p4().Eta()         if len(NuFromWFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromT_eta'%self.label         , [ NuFromWFromT[i].p4().Eta()         if len(NuFromWFromT) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromH_eta'%self.label          , [ QFromWFromH[i].p4().Eta()          if len(QFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromT_eta'%self.label          , [ QFromWFromT[i].p4().Eta()          if len(QFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromH_eta'%self.label          , [ LFromWFromH[i].p4().Eta()          if len(LFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromT_eta'%self.label          , [ LFromWFromT[i].p4().Eta()          if len(LFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromH_eta'%self.label             , [ TauFromH[i].p4().Eta()             if len(TauFromH) > i            else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromH_eta'%self.label        , [ hFromTauFromH[i].p4().Eta()        if len(hFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromH_eta'%self.label        , [ LFromTauFromH[i].p4().Eta()        if len(LFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromH_eta'%self.label       , [ NuFromTauFromH[i].p4().Eta()       if len(NuFromTauFromH) > i      else -99 for i in range(4) ])
        self.out.fillBranch('%sTNuFromTauFromH_eta'%self.label      , [ TNuFromTauFromH[i].p4().Eta()      if len(TNuFromTauFromH) > i     else -99 for i in range(4) ])
        self.out.fillBranch('%sZFromH_eta'%self.label               , [ ZFromH[i].p4().Eta()               if len(ZFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromZFromH_eta'%self.label          , [ QFromZFromH[i].p4().Eta()          if len(QFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromZFromH_eta'%self.label          , [ LFromZFromH[i].p4().Eta()          if len(LFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromZFromH_eta'%self.label         , [ NuFromZFromH[i].p4().Eta()         if len(NuFromZFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromWFromH_eta'%self.label        , [ TauFromWFromH[i].p4().Eta()        if len(TauFromWFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromWFromH_eta'%self.label   , [ LFromTauFromWFromH[i].p4().Eta()   if len(LFromTauFromWFromH) > i  else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromWFromH_eta'%self.label  , [ NuFromTauFromWFromH[i].p4().Eta()  if len(NuFromTauFromWFromH) > i else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromWFromH_eta'%self.label   , [ hFromTauFromWFromH[i].p4().Eta()   if len(hFromTauFromWFromH) > i  else -99 for i in range(4) ])

        self.out.fillBranch('%sHiggses_phi'%self.label              , [ Higgses[i].p4().Phi()              if len(Higgses) > i             else -99 for i in range(4) ]) 
        self.out.fillBranch('%sTfromhardprocess_phi'%self.label     , [ Tfromhardprocess[i].p4().Phi()     if len(Tfromhardprocess) > i    else -99 for i in range(4) ]) 
        self.out.fillBranch('%sHadT_phi'%self.label                 , [ HadT[i].p4().Phi()              if len(HadT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sLepT_phi'%self.label                 , [ LepT[i].p4().Phi()                 if len(LepT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromHadT_phi'%self.label            , [ WFromHadT[i].p4().Phi()            if len(WFromHadT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromLepT_phi'%self.label            , [ WFromLepT[i].p4().Phi()            if len(WFromLepT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromH_phi'%self.label               , [ WFromH[i].p4().Phi()               if len(WFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sHadWFromH_phi'%self.label            , [ HadWFromH[i].p4().Phi()            if len(HadWFromH) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromT_phi'%self.label               , [ WFromT[i].p4().Phi()               if len(WFromT) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromW_phi'%self.label               , [ QFromW[i].p4().Phi()               if len(QFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sGenLep_phi'%self.label               , [ GenLep[i].p4().Phi()               if len(GenLep) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromW_phi'%self.label               , [ LFromW[i].p4().Phi()               if len(LFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromH_phi'%self.label         , [ NuFromWFromH[i].p4().Phi()         if len(NuFromWFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromT_phi'%self.label         , [ NuFromWFromT[i].p4().Phi()         if len(NuFromWFromT) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromH_phi'%self.label          , [ QFromWFromH[i].p4().Phi()          if len(QFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromT_phi'%self.label          , [ QFromWFromT[i].p4().Phi()          if len(QFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromH_phi'%self.label          , [ LFromWFromH[i].p4().Phi()          if len(LFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromT_phi'%self.label          , [ LFromWFromT[i].p4().Phi()          if len(LFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromH_phi'%self.label             , [ TauFromH[i].p4().Phi()             if len(TauFromH) > i            else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromH_phi'%self.label        , [ hFromTauFromH[i].p4().Phi()        if len(hFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromH_phi'%self.label        , [ LFromTauFromH[i].p4().Phi()        if len(LFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromH_phi'%self.label       , [ NuFromTauFromH[i].p4().Phi()       if len(NuFromTauFromH) > i      else -99 for i in range(4) ])
        self.out.fillBranch('%sTNuFromTauFromH_phi'%self.label      , [ TNuFromTauFromH[i].p4().Phi()      if len(TNuFromTauFromH) > i     else -99 for i in range(4) ])
        self.out.fillBranch('%sZFromH_phi'%self.label               , [ ZFromH[i].p4().Phi()               if len(ZFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromZFromH_phi'%self.label          , [ QFromZFromH[i].p4().Phi()          if len(QFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromZFromH_phi'%self.label          , [ LFromZFromH[i].p4().Phi()          if len(LFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromZFromH_phi'%self.label         , [ NuFromZFromH[i].p4().Phi()         if len(NuFromZFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromWFromH_phi'%self.label        , [ TauFromWFromH[i].p4().Phi()        if len(TauFromWFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromWFromH_phi'%self.label   , [ LFromTauFromWFromH[i].p4().Phi()   if len(LFromTauFromWFromH) > i  else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromWFromH_phi'%self.label  , [ NuFromTauFromWFromH[i].p4().Phi()  if len(NuFromTauFromWFromH) > i else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromWFromH_phi'%self.label   , [ hFromTauFromWFromH[i].p4().Phi()   if len(hFromTauFromWFromH) > i  else -99 for i in range(4) ])

        self.out.fillBranch('%sHiggses_mass'%self.label             , [ Higgses[i].p4().M()                if len(Higgses) > i             else -99 for i in range(4) ]) 
        self.out.fillBranch('%sTfromhardprocess_mass'%self.label    , [ Tfromhardprocess[i].p4().M()       if len(Tfromhardprocess) > i    else -99 for i in range(4) ]) 
        self.out.fillBranch('%sHadT_mass'%self.label                , [ HadT[i].p4().M()                if len(HadT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sLepT_mass'%self.label                , [ LepT[i].p4().M()                   if len(LepT) > i                else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromHadT_mass'%self.label           , [ WFromHadT[i].p4().M()              if len(WFromHadT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromLepT_mass'%self.label           , [ WFromLepT[i].p4().M()              if len(WFromLepT) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromH_mass'%self.label              , [ WFromH[i].p4().M()                 if len(WFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sHadWFromH_mass'%self.label           , [ HadWFromH[i].p4().M()              if len(HadWFromH) > i           else -99 for i in range(4) ])
        self.out.fillBranch('%sWFromT_mass'%self.label              , [ WFromT[i].p4().M()                 if len(WFromT) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromW_mass'%self.label              , [ QFromW[i].p4().M()                 if len(QFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sGenLep_mass'%self.label              , [ GenLep[i].p4().M()                 if len(GenLep) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromW_mass'%self.label              , [ LFromW[i].p4().M()                 if len(LFromW) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromH_mass'%self.label        , [ NuFromWFromH[i].p4().M()           if len(NuFromWFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromWFromT_mass'%self.label        , [ NuFromWFromT[i].p4().M()           if len(NuFromWFromT) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromH_mass'%self.label         , [ QFromWFromH[i].p4().M()            if len(QFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromWFromT_mass'%self.label         , [ QFromWFromT[i].p4().M()            if len(QFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromH_mass'%self.label         , [ LFromWFromH[i].p4().M()            if len(LFromWFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromWFromT_mass'%self.label         , [ LFromWFromT[i].p4().M()            if len(LFromWFromT) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromH_mass'%self.label            , [ TauFromH[i].p4().M()               if len(TauFromH) > i            else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromH_mass'%self.label       , [ hFromTauFromH[i].p4().M()          if len(hFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromH_mass'%self.label       , [ LFromTauFromH[i].p4().M()          if len(LFromTauFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromH_mass'%self.label      , [ NuFromTauFromH[i].p4().M()         if len(NuFromTauFromH) > i      else -99 for i in range(4) ])
        self.out.fillBranch('%sTNuFromTauFromH_mass'%self.label     , [ TNuFromTauFromH[i].p4().M()        if len(TNuFromTauFromH) > i     else -99 for i in range(4) ])
        self.out.fillBranch('%sZFromH_mass'%self.label              , [ ZFromH[i].p4().M()                 if len(ZFromH) > i              else -99 for i in range(4) ])
        self.out.fillBranch('%sQFromZFromH_mass'%self.label         , [ QFromZFromH[i].p4().M()            if len(QFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromZFromH_mass'%self.label         , [ LFromZFromH[i].p4().M()            if len(LFromZFromH) > i         else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromZFromH_mass'%self.label        , [ NuFromZFromH[i].p4().M()           if len(NuFromZFromH) > i        else -99 for i in range(4) ])
        self.out.fillBranch('%sTauFromWFromH_mass'%self.label       , [ TauFromWFromH[i].p4().M()          if len(TauFromWFromH) > i       else -99 for i in range(4) ])
        self.out.fillBranch('%sLFromTauFromWFromH_mass'%self.label  , [ LFromTauFromWFromH[i].p4().M()     if len(LFromTauFromWFromH) > i  else -99 for i in range(4) ])
        self.out.fillBranch('%sNuFromTauFromWFromH_mass'%self.label , [ NuFromTauFromWFromH[i].p4().M()    if len(NuFromTauFromWFromH) > i else -99 for i in range(4) ])
        self.out.fillBranch('%shFromTauFromWFromH_mass'%self.label  , [ hFromTauFromWFromH[i].p4().M()     if len(hFromTauFromWFromH) > i  else -99 for i in range(4) ])

        # Fill branches for some precomputed variables
        self.out.fillBranch('%spTHgen'%self.label, Higgses[0].p4().Pt() if len(Higgses)==1 else -99)
        self.out.fillBranch('%spTtgen'%self.label, [ Tfromhardprocess[i].p4().Pt() if len(Tfromhardprocess) > i else -99 for i in range(4) ]) 
        self.out.fillBranch('%sinv_mass_q1_q2'%self.label, (QFromWFromH[0].p4()+QFromWFromH[1].p4()).M() if len(QFromWFromH)==2 else -99)
        self.out.fillBranch('%sdelR_partonsFromH'%self.label, QFromWFromH[0].p4().DeltaR(QFromWFromH[1].p4()) if len(QFromWFromH)==2 else -99)
        self.out.fillBranch('%spt_partonsFromH'%self.label, (QFromWFromH[0].p4()+QFromWFromH[1].p4()).Pt() if len(QFromWFromH)==2 else -99)

        self.out.fillBranch('%squark1pT'%self.label, QFromWFromH[0].p4().Pt() if len(QFromWFromH)==2 else -99)
        self.out.fillBranch('%squark2pT'%self.label, QFromWFromH[1].p4().Pt() if len(QFromWFromH)==2 else -99)
        self.out.fillBranch('%squark1Eta'%self.label, QFromWFromH[0].p4().Eta() if len(QFromWFromH)==2 else -99)
        self.out.fillBranch('%squark2Eta'%self.label, QFromWFromH[1].p4().Eta() if len(QFromWFromH)==2 else -99)
        q1pdgid = -99
        q2pdgid = -99
        if len(QFromWFromH)==2:
#           q1pdgid = QFromWFromH[0].pdgId if QFromWFromH[0].pdgId == 4 or QFromWFromH[0].pdgId == 5 else 0
#           q2pdgid = QFromWFromH[1].pdgId if QFromWFromH[1].pdgId == 4 or QFromWFromH[1].pdgId == 5 else 0
            q1pdgid = QFromWFromH[0].pdgId
            q2pdgid = QFromWFromH[1].pdgId
        self.out.fillBranch('%squark1Flavour'%self.label, q1pdgid)
        self.out.fillBranch('%squark2Flavour'%self.label, q2pdgid)

        if len(LFromWFromH)==1:
            self.out.fillBranch('%sdelR_H_q1l'%self.label, LFromWFromH[0].p4().DeltaR(QFromWFromH[0].p4()) if len(QFromWFromH)==2 else -99)
            self.out.fillBranch('%sdelR_H_q2l'%self.label, LFromWFromH[0].p4().DeltaR(QFromWFromH[1].p4()) if len(QFromWFromH)==2 else -99)
        
        else:
            # Different error code to highlight different cause
            self.out.fillBranch('%sdelR_H_q1l'%self.label, -88)
            self.out.fillBranch('%sdelR_H_q2l'%self.label, -88)

        # Build pt of Higgs decay particles
        if len(LFromWFromH)==1 and len(QFromWFromH)==2:
            trueVisibleHiggs = LFromWFromH[0].p4() + QFromWFromH[0].p4() + QFromWFromH[1].p4()
            self.out.fillBranch('%spTTrueGen'%self.label, trueVisibleHiggs.Pt())
            self.out.fillBranch('%sMTrueGen'%self.label, trueVisibleHiggs.M())
            if len(NuFromWFromH)==1:
                trueFullHiggs=trueVisibleHiggs+NuFromWFromH[0].p4()
                self.out.fillBranch('%spTTrueGenPlusNu'%self.label, trueFullHiggs.Pt())
            else:
                # Different error code to highligh different cause
                self.out.fillBranch('%spTTrueGenPlusNu'%self.label, -88)
        else:
            self.out.fillBranch('%spTTrueGen'%self.label, -99)
            self.out.fillBranch('%sMTrueGen'%self.label, -99)
            self.out.fillBranch('%spTTrueGenPlusNu'%self.label, -99)
            


        # Necessary.
        return True



higgsDiffGenTTH = lambda : HiggsDiffGenTTH()
