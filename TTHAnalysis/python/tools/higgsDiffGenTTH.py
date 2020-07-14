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
        self.out.branch('%snWFromH'%self.label          , 'I')
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
        self.out.branch('%snPiFromTauFromH'%self.label  , 'I')
        self.out.branch('%snKaFromTauFromH'%self.label  , 'I')
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
            
            self.out.branch('%sHiggses%s'%(self.label,suffix)          , 'F', 4, '%snHiggses'%self.label         ) 
            self.out.branch('%sTfromhardprocess%s'%(self.label,suffix) , 'F', 4, '%snTfromhardprocess'%self.label) 
            self.out.branch('%sWFromH%s'%(self.label,suffix)           , 'F', 4, '%snWFromH'%self.label          ) 
            self.out.branch('%sWFromT%s'%(self.label,suffix)           , 'F', 4, '%snWFromT'%self.label          ) 
            self.out.branch('%sQFromW%s'%(self.label,suffix)           , 'F', 4, '%snQFromW'%self.label          ) 
            self.out.branch('%sGenLep%s'%(self.label,suffix)           , 'F', 4, '%snGenLep'%self.label          ) 
            self.out.branch('%sLFromW%s'%(self.label,suffix)           , 'F', 4, '%snLFromW'%self.label          ) 
            self.out.branch('%sNuFromWFromH%s'%(self.label,suffix)     , 'F', 4, '%snNuFromWFromH'%self.label    ) 
            self.out.branch('%sNuFromWFromT%s'%(self.label,suffix)     , 'F', 4, '%snNuFromWFromT'%self.label    ) 
            self.out.branch('%sQFromWFromH%s'%(self.label,suffix)      , 'F', 4, '%snQFromWFromH'%self.label     ) 
            self.out.branch('%sQFromWFromT%s'%(self.label,suffix)      , 'F', 4, '%snQFromWFromT'%self.label     ) 
            self.out.branch('%sLFromWFromH%s'%(self.label,suffix)      , 'F', 4, '%snLFromWFromH'%self.label     ) 
            self.out.branch('%sLFromWFromT%s'%(self.label,suffix)      , 'F', 4, '%snLFromWFromT'%self.label     ) 
            self.out.branch('%sTauFromH%s'%(self.label,suffix)          , 'F', 4, '%snTauFromH'%self.label        )
            self.out.branch('%sPiFromTauFromH%s'%(self.label,suffix)    , 'F', 4, '%snPiFromTauFromH'%self.label  )
            self.out.branch('%sKaFromTauFromH%s'%(self.label,suffix)    , 'F', 4, '%snKaFromTauFromH'%self.label  )
            self.out.branch('%sLFromTauFromH%s'%(self.label,suffix)     , 'F', 4, '%snLFromTauFromH'%self.label   )
            self.out.branch('%sNuFromTauFromH%s'%(self.label,suffix)    , 'F', 4, '%snNuFromTauFromH'%self.label  )
            self.out.branch('%sTNuFromTauFromH%s'%(self.label,suffix)   , 'F', 4, '%snTNuFromTauFromH'%self.label )
            self.out.branch('%sZFromH%s'%(self.label,suffix)            , 'F', 4, '%snZFromH'%self.label          )
            self.out.branch('%sQFromZFromH%s'%(self.label,suffix)       , 'F', 4, '%snQFromZFromH'%self.label     )
            self.out.branch('%sLFromZFromH%s'%(self.label,suffix)       , 'F', 4, '%snLFromZFromH'%self.label     )
            self.out.branch('%sNuFromZFromH%s'%(self.label,suffix)      , 'F', 4, '%snNuFromZFromH'%self.label    )
            self.out.branch('%sTauFromWFromH%s'%(self.label,suffix)     , 'F', 4, '%snTauFromWFromH'%self.label   )
            self.out.branch('%sLFromTauFromWFromH%s'%(self.label,suffix)   , 'F', 4, '%snLFromTauFromWFromH'%self.label   )
            self.out.branch('%sNuFromTauFromWFromH%s'%(self.label,suffix)  , 'F', 4, '%snNuFromTauFromWFromH'%self.label  )
            self.out.branch('%shFromTauFromWFromH%s'%(self.label,suffix)   , 'F', 4, '%snhFromTauFromWFromH'%self.label   )
        # Some precomputed quantities of interest

        self.out.branch('%spTHgen'%self.label            , 'F')
        self.out.branch('%spTtgen'%self.label            , 'F', 4, '%snTfromhardprocess'%self.label)
        self.out.branch('%sinv_mass_q1_q2'%self.label    , 'F')
        self.out.branch('%sdelR_partonsFromH'%self.label , 'F')
        self.out.branch('%squark1pT'%self.label          , 'F')
        self.out.branch('%squark2pT'%self.label          , 'F')
        self.out.branch('%sdelR_H_q1l'%self.label        , 'F')
        self.out.branch('%sdelR_H_q2l'%self.label        , 'F')
        self.out.branch('%spTTrueGen'%self.label         , 'F')
        self.out.branch('%spTTrueGenPlusNu'%self.label   , 'F')
        
    def analyze(self, event):
        # Input collections and maps
        genpar = Collection(event,"GenPart","nGenPart") 
        
        Higgses          = []
        Tfromhardprocess = []
        WFromH           = []
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
        PiFromTauFromH   = []
        KaFromTauFromH   = []
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
            if abs(part.pdgId) in [12, 14]:
                if self.debug: print "it is a neutrino"
                if part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24: 
                    if self.debug: print "the mother of this neutrino is W+ or W-"
                    if abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId) == 6:
                        NuFromWFromT.append(part)
                        if self.debug: print "the mother of this W is a top"
        
            # quarks from W from H
            if (abs(part.pdgId) in [1,2,3,4,5,6] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24
                     and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                if self.debug: print "it is a quark coming from a hard W"
                if (genpar[genpar[part.genPartIdxMother].genPartIdxMother].genPartIdxMother >= 0 
                    and genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId == 25
                    and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                    QFromWFromH.append(part)
                    if self.debug: print "the mother of this hard W is a hard Higgs"

            # quarks from W from T 
            if (abs(part.pdgId) in [1,2,3,4,5,6] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId) == 24 
                     and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                if self.debug: print "it is a quark coming from a hard W"
                if (genpar[genpar[part.genPartIdxMother].genPartIdxMother].genPartIdxMother >= 0 and abs(genpar[genpar[part.genPartIdxMother].genPartIdxMother].pdgId) == 6
                        and genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                    QFromWFromT.append(part)
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
                if(genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
                    TauFromH.append(part)
                    if self.debug: print "it is a hard tau coming from a higgs"

            # Pion from tau from H
            if (abs(part.pdgId) in [111,211] and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcessTauDecayProduct'])):
                PiFromTauFromH.append(part)
                if self.debug: print "it is a pion coming from a hard tau from a higgs"

            # Kaon from tau from H
            if (abs(part.pdgId) in [311,321] and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcessTauDecayProduct'])):
                KaFromTauFromH.append(part)
                if self.debug: print "it is a kaon coming from a hard tau from a higgs"

            # Lepton from tau from H
            if (abs(part.pdgId) in [11,13] and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcessTauDecayProduct'])):
                LFromTauFromH.append(part)
                if self.debug: print "it is a lepton coming from a hard tau from a higgs"

            # lnu from tau from H
            if (abs(part.pdgId) in [12,14] and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcessTauDecayProduct'])):
                NuFromTauFromH.append(part)
                if self.debug: print "it is a lepton neutrino coming from a hard tau from a higgs"

            # tnu from tau from H
            if (abs(part.pdgId) in [16] and part.statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcessTauDecayProduct'])):
                TNuFromTauFromH.append(part)
                if self.debug: print "it is a tau neutrino coming from a hard tau from a higgs"

            # Z from H
            if (abs(part.pdgId) in [23] and part.genPartIdxMother >= 0 and abs(genpar[part.genPartIdxMother].pdgId)==25):
                if(genpar[part.genPartIdxMother].statusFlags &(1 << diffUtils.statusFlagsMap['isHardProcess'])):
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
        self.out.fillBranch('%snWFromH'%self.label          , len(WFromH))
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
        self.out.fillBranch('%snPiFromTauFromH'%self.label  , len(PiFromTauFromH))
        self.out.fillBranch('%snKaFromTauFromH'%self.label  , len(KaFromTauFromH))
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

        self.out.fillBranch('%sHiggses_pt'%self.label          , [ part.p4().Pt() for part in Higgses         ]) 
        self.out.fillBranch('%sTfromhardprocess_pt'%self.label , [ part.p4().Pt() for part in Tfromhardprocess]) 
        self.out.fillBranch('%sWFromH_pt'%self.label           , [ part.p4().Pt() for part in WFromH          ]) 
        self.out.fillBranch('%sWFromT_pt'%self.label           , [ part.p4().Pt() for part in WFromT          ]) 
        self.out.fillBranch('%sQFromW_pt'%self.label           , [ part.p4().Pt() for part in QFromW          ]) 
        self.out.fillBranch('%sGenLep_pt'%self.label           , [ part.p4().Pt() for part in GenLep          ]) 
        self.out.fillBranch('%sLFromW_pt'%self.label           , [ part.p4().Pt() for part in LFromW          ]) 
        self.out.fillBranch('%sNuFromWFromH_pt'%self.label     , [ part.p4().Pt() for part in NuFromWFromH    ]) 
        self.out.fillBranch('%sNuFromWFromT_pt'%self.label     , [ part.p4().Pt() for part in NuFromWFromT    ]) 
        self.out.fillBranch('%sQFromWFromH_pt'%self.label      , [ part.p4().Pt() for part in QFromWFromH     ]) 
        self.out.fillBranch('%sQFromWFromT_pt'%self.label      , [ part.p4().Pt() for part in QFromWFromT     ]) 
        self.out.fillBranch('%sLFromWFromH_pt'%self.label      , [ part.p4().Pt() for part in LFromWFromH     ]) 
        self.out.fillBranch('%sLFromWFromT_pt'%self.label      , [ part.p4().Pt() for part in LFromWFromT     ]) 
        self.out.fillBranch('%sTauFromH_pt'%self.label         , [ part.p4().Pt() for part in TauFromH        ])
        self.out.fillBranch('%sPiFromTauFromH_pt'%self.label   , [ part.p4().Pt() for part in PiFromTauFromH  ])
        self.out.fillBranch('%sKaFromTauFromH_pt'%self.label   , [ part.p4().Pt() for part in KaFromTauFromH  ])
        self.out.fillBranch('%sLFromTauFromH_pt'%self.label    , [ part.p4().Pt() for part in LFromTauFromH   ])
        self.out.fillBranch('%sNuFromTauFromH_pt'%self.label   , [ part.p4().Pt() for part in NuFromTauFromH  ])
        self.out.fillBranch('%sTNuFromTauFromH_pt'%self.label  , [ part.p4().Pt() for part in TNuFromTauFromH ])
        self.out.fillBranch('%sZFromH_pt'%self.label           , [ part.p4().Pt() for part in ZFromH          ])
        self.out.fillBranch('%sQFromZFromH_pt'%self.label      , [ part.p4().Pt() for part in QFromZFromH     ])
        self.out.fillBranch('%sLFromZFromH_pt'%self.label      , [ part.p4().Pt() for part in LFromZFromH     ])
        self.out.fillBranch('%sNuFromZFromH_pt'%self.label     , [ part.p4().Pt() for part in NuFromZFromH    ])
        self.out.fillBranch('%sTauFromWFromH_pt'%self.label    , [ part.p4().Pt() for part in TauFromWFromH   ])
        self.out.fillBranch('%sLFromTauFromWFromH_pt'%self.label   , [ part.p4().Pt() for part in LFromTauFromWFromH  ])
        self.out.fillBranch('%sNuFromTauFromWFromH_pt'%self.label  , [ part.p4().Pt() for part in NuFromTauFromWFromH ])
        self.out.fillBranch('%shFromTauFromWFromH_pt'%self.label   , [ part.p4().Pt() for part in hFromTauFromWFromH  ])

        self.out.fillBranch('%sHiggses_eta'%self.label          , [ part.p4().Eta() for part in Higgses         ]) 
        self.out.fillBranch('%sTfromhardprocess_eta'%self.label , [ part.p4().Eta() for part in Tfromhardprocess]) 
        self.out.fillBranch('%sWFromH_eta'%self.label           , [ part.p4().Eta() for part in WFromH          ]) 
        self.out.fillBranch('%sWFromT_eta'%self.label           , [ part.p4().Eta() for part in WFromT          ]) 
        self.out.fillBranch('%sQFromW_eta'%self.label           , [ part.p4().Eta() for part in QFromW          ]) 
        self.out.fillBranch('%sGenLep_eta'%self.label           , [ part.p4().Eta() for part in GenLep          ]) 
        self.out.fillBranch('%sLFromW_eta'%self.label           , [ part.p4().Eta() for part in LFromW          ]) 
        self.out.fillBranch('%sNuFromWFromH_eta'%self.label     , [ part.p4().Eta() for part in NuFromWFromH    ]) 
        self.out.fillBranch('%sNuFromWFromT_eta'%self.label     , [ part.p4().Eta() for part in NuFromWFromT    ]) 
        self.out.fillBranch('%sQFromWFromH_eta'%self.label      , [ part.p4().Eta() for part in QFromWFromH     ]) 
        self.out.fillBranch('%sQFromWFromT_eta'%self.label      , [ part.p4().Eta() for part in QFromWFromT     ]) 
        self.out.fillBranch('%sLFromWFromH_eta'%self.label      , [ part.p4().Eta() for part in LFromWFromH     ]) 
        self.out.fillBranch('%sLFromWFromT_eta'%self.label      , [ part.p4().Eta() for part in LFromWFromT     ]) 
        self.out.fillBranch('%sTauFromH_eta'%self.label         , [ part.p4().Eta() for part in TauFromH        ])
        self.out.fillBranch('%sPiFromTauFromH_eta'%self.label   , [ part.p4().Eta() for part in PiFromTauFromH  ])
        self.out.fillBranch('%sKaFromTauFromH_eta'%self.label   , [ part.p4().Eta() for part in KaFromTauFromH  ])
        self.out.fillBranch('%sLFromTauFromH_eta'%self.label    , [ part.p4().Eta() for part in LFromTauFromH   ])
        self.out.fillBranch('%sNuFromTauFromH_eta'%self.label   , [ part.p4().Eta() for part in NuFromTauFromH  ])
        self.out.fillBranch('%sTNuFromTauFromH_eta'%self.label  , [ part.p4().Eta() for part in TNuFromTauFromH ])
        self.out.fillBranch('%sZFromH_eta'%self.label           , [ part.p4().Eta() for part in ZFromH          ])
        self.out.fillBranch('%sQFromZFromH_eta'%self.label      , [ part.p4().Eta() for part in QFromZFromH     ])
        self.out.fillBranch('%sLFromZFromH_eta'%self.label      , [ part.p4().Eta() for part in LFromZFromH     ])
        self.out.fillBranch('%sNuFromZFromH_eta'%self.label     , [ part.p4().Eta() for part in NuFromZFromH    ])
        self.out.fillBranch('%sTauFromWFromH_eta'%self.label    , [ part.p4().Eta() for part in TauFromWFromH   ])
        self.out.fillBranch('%sLFromTauFromWFromH_eta'%self.label   , [ part.p4().Eta() for part in LFromTauFromWFromH  ])
        self.out.fillBranch('%sNuFromTauFromWFromH_eta'%self.label  , [ part.p4().Eta() for part in NuFromTauFromWFromH ])
        self.out.fillBranch('%shFromTauFromWFromH_eta'%self.label   , [ part.p4().Eta() for part in hFromTauFromWFromH  ])

        self.out.fillBranch('%sHiggses_phi'%self.label          , [ part.p4().Phi() for part in Higgses         ]) 
        self.out.fillBranch('%sTfromhardprocess_phi'%self.label , [ part.p4().Phi() for part in Tfromhardprocess]) 
        self.out.fillBranch('%sWFromH_phi'%self.label           , [ part.p4().Phi() for part in WFromH          ]) 
        self.out.fillBranch('%sWFromT_phi'%self.label           , [ part.p4().Phi() for part in WFromT          ]) 
        self.out.fillBranch('%sQFromW_phi'%self.label           , [ part.p4().Phi() for part in QFromW          ]) 
        self.out.fillBranch('%sGenLep_phi'%self.label           , [ part.p4().Phi() for part in GenLep          ]) 
        self.out.fillBranch('%sLFromW_phi'%self.label           , [ part.p4().Phi() for part in LFromW          ]) 
        self.out.fillBranch('%sNuFromWFromH_phi'%self.label     , [ part.p4().Phi() for part in NuFromWFromH    ]) 
        self.out.fillBranch('%sNuFromWFromT_phi'%self.label     , [ part.p4().Phi() for part in NuFromWFromT    ]) 
        self.out.fillBranch('%sQFromWFromH_phi'%self.label      , [ part.p4().Phi() for part in QFromWFromH     ]) 
        self.out.fillBranch('%sQFromWFromT_phi'%self.label      , [ part.p4().Phi() for part in QFromWFromT     ]) 
        self.out.fillBranch('%sLFromWFromH_phi'%self.label      , [ part.p4().Phi() for part in LFromWFromH     ]) 
        self.out.fillBranch('%sLFromWFromT_phi'%self.label      , [ part.p4().Phi() for part in LFromWFromT     ]) 
        self.out.fillBranch('%sTauFromH_phi'%self.label         , [ part.p4().Phi() for part in TauFromH        ])
        self.out.fillBranch('%sPiFromTauFromH_phi'%self.label   , [ part.p4().Phi() for part in PiFromTauFromH  ])
        self.out.fillBranch('%sKaFromTauFromH_phi'%self.label   , [ part.p4().Phi() for part in KaFromTauFromH  ])
        self.out.fillBranch('%sLFromTauFromH_phi'%self.label    , [ part.p4().Phi() for part in LFromTauFromH   ])
        self.out.fillBranch('%sNuFromTauFromH_phi'%self.label   , [ part.p4().Phi() for part in NuFromTauFromH  ])
        self.out.fillBranch('%sTNuFromTauFromH_phi'%self.label  , [ part.p4().Phi() for part in TNuFromTauFromH ])
        self.out.fillBranch('%sZFromH_phi'%self.label           , [ part.p4().Phi() for part in ZFromH          ])
        self.out.fillBranch('%sQFromZFromH_phi'%self.label      , [ part.p4().Phi() for part in QFromZFromH     ])
        self.out.fillBranch('%sLFromZFromH_phi'%self.label      , [ part.p4().Phi() for part in LFromZFromH     ])
        self.out.fillBranch('%sNuFromZFromH_phi'%self.label     , [ part.p4().Phi() for part in NuFromZFromH    ])
        self.out.fillBranch('%sTauFromWFromH_phi'%self.label    , [ part.p4().Phi() for part in TauFromWFromH   ])
        self.out.fillBranch('%sLFromTauFromWFromH_phi'%self.label   , [ part.p4().Phi() for part in LFromTauFromWFromH  ])
        self.out.fillBranch('%sNuFromTauFromWFromH_phi'%self.label  , [ part.p4().Phi() for part in NuFromTauFromWFromH ])
        self.out.fillBranch('%shFromTauFromWFromH_phi'%self.label   , [ part.p4().Phi() for part in hFromTauFromWFromH  ])

        self.out.fillBranch('%sHiggses_mass'%self.label          , [ part.p4().M() for part in Higgses         ]) 
        self.out.fillBranch('%sTfromhardprocess_mass'%self.label , [ part.p4().M() for part in Tfromhardprocess]) 
        self.out.fillBranch('%sWFromH_mass'%self.label           , [ part.p4().M() for part in WFromH          ]) 
        self.out.fillBranch('%sWFromT_mass'%self.label           , [ part.p4().M() for part in WFromT          ]) 
        self.out.fillBranch('%sQFromW_mass'%self.label           , [ part.p4().M() for part in QFromW          ]) 
        self.out.fillBranch('%sGenLep_mass'%self.label           , [ part.p4().M() for part in GenLep          ]) 
        self.out.fillBranch('%sLFromW_mass'%self.label           , [ part.p4().M() for part in LFromW          ]) 
        self.out.fillBranch('%sNuFromWFromH_mass'%self.label     , [ part.p4().M() for part in NuFromWFromH    ]) 
        self.out.fillBranch('%sNuFromWFromT_mass'%self.label     , [ part.p4().M() for part in NuFromWFromT    ]) 
        self.out.fillBranch('%sQFromWFromH_mass'%self.label      , [ part.p4().M() for part in QFromWFromH     ]) 
        self.out.fillBranch('%sQFromWFromT_mass'%self.label      , [ part.p4().M() for part in QFromWFromT     ]) 
        self.out.fillBranch('%sLFromWFromH_mass'%self.label      , [ part.p4().M() for part in LFromWFromH     ]) 
        self.out.fillBranch('%sLFromWFromT_mass'%self.label      , [ part.p4().M() for part in LFromWFromT     ]) 
        self.out.fillBranch('%sTauFromH_mass'%self.label         , [ part.p4().M() for part in TauFromH        ])
        self.out.fillBranch('%sPiFromTauFromH_mass'%self.label   , [ part.p4().M() for part in PiFromTauFromH  ])
        self.out.fillBranch('%sKaFromTauFromH_mass'%self.label   , [ part.p4().M() for part in KaFromTauFromH  ])
        self.out.fillBranch('%sLFromTauFromH_mass'%self.label    , [ part.p4().M() for part in LFromTauFromH   ])
        self.out.fillBranch('%sNuFromTauFromH_mass'%self.label   , [ part.p4().M() for part in NuFromTauFromH  ])
        self.out.fillBranch('%sTNuFromTauFromH_mass'%self.label  , [ part.p4().M() for part in TNuFromTauFromH ])
        self.out.fillBranch('%sZFromH_mass'%self.label           , [ part.p4().M() for part in ZFromH          ])
        self.out.fillBranch('%sQFromZFromH_mass'%self.label      , [ part.p4().M() for part in QFromZFromH     ])
        self.out.fillBranch('%sLFromZFromH_mass'%self.label      , [ part.p4().M() for part in LFromZFromH     ])
        self.out.fillBranch('%sNuFromZFromH_mass'%self.label     , [ part.p4().M() for part in NuFromZFromH    ])
        self.out.fillBranch('%sTauFromWFromH_mass'%self.label    , [ part.p4().M() for part in TauFromWFromH   ])
        self.out.fillBranch('%sLFromTauFromWFromH_mass'%self.label   , [ part.p4().M() for part in LFromTauFromWFromH  ])
        self.out.fillBranch('%sNuFromTauFromWFromH_mass'%self.label  , [ part.p4().M() for part in NuFromTauFromWFromH ])
        self.out.fillBranch('%shFromTauFromWFromH_mass'%self.label   , [ part.p4().M() for part in hFromTauFromWFromH  ])

        # Fill branches for some precomputed variables
        self.out.fillBranch('%spTHgen'%self.label, Higgses[0].p4().Pt() if len(Higgses)==1 else -99)
        self.out.fillBranch('%spTtgen'%self.label, [part.p4().Pt() for part in Tfromhardprocess] ) 
        self.out.fillBranch('%sinv_mass_q1_q2'%self.label, (QFromWFromH[0].p4()+QFromWFromH[1].p4()).M() if len(QFromWFromH)==2 else -99)
        self.out.fillBranch('%sdelR_partonsFromH'%self.label, QFromWFromH[0].p4().DeltaR(QFromWFromH[1].p4()) if len(QFromWFromH)==2 else -99)

        self.out.fillBranch('%squark1pT'%self.label, QFromWFromH[0].p4().Pt() if len(QFromWFromH)==2 else -99)
        self.out.fillBranch('%squark2pT'%self.label, QFromWFromH[1].p4().Pt() if len(QFromWFromH)==2 else -99)

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
            if len(NuFromWFromH)==1:
                trueFullHiggs=trueVisibleHiggs+NuFromWFromH[0].p4()
                self.out.fillBranch('%spTTrueGenPlusNu'%self.label, trueFullHiggs.Pt())
            else:
                # Different error code to highligh different cause
                self.out.fillBranch('%spTTrueGenPlusNu'%self.label, -88)
        else:
            self.out.fillBranch('%spTTrueGen'%self.label, -99)
            self.out.fillBranch('%spTTrueGenPlusNu'%self.label, -99)
            


        # Necessary.
        return True



higgsDiffGenTTH = lambda : HiggsDiffGenTTH()
