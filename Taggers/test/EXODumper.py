#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager


process = cms.Process("EXOTagDumper")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
process.load("flashgg.Systematics.flashggMuonSystematics_cfi")
process.load("flashgg.Systematics.flashggElectronSystematics_cfi")
process.load("flashgg.Systematics.flashggJetSystematics_cfi")
process.flashggDiPhotonSystematics.src='flashggDiPhotons'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 2500 )



process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_387.root"
                                "/store/group/phys_higgs/cmshgg/musella/flashgg/EXOSpring16_v1_p4/diphotons_80_v1/RSGravToGG_kMpl-001_M-750_TuneCUEP8M1_13TeV-pythia8/EXOSpring16_v1_p4-diphotons_80_v1-v0-RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/160527_170715/0000/diphotonsMicroAOD_1.root"
                                ))

import re
testString = str(process.source.fileNames)

print "\n\n"
fName = re.findall(r'\'([^]]*)\'', testString)[0]
if "DoubleEG" in fName:
    print "Use flashggSelectedElectrons"
    electronString = "flashggSelectedElectrons"
else:
    print "Use flashggElectrons"
    electronString = "flashggElectrons"
print "\n\n"

process.TFileService = cms.Service( "TFileService",
                                    fileName = cms.string("EXOTagsDump.root"),
                                    closeFileFast = cms.untracked.bool(True) )


import flashgg.Taggers.dumperConfigTools as cfgTools
from  flashgg.Taggers.tagsDumpers_cfi import createTagDumper


process.flashggUnpackedJets = cms.EDProducer( "FlashggVectorVectorJetUnpacker",
                                              JetsTag = cms.InputTag("flashggFinalJets"),
                                              NCollections = cms.uint32(8))

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
process.flashggEXOTag = cms.EDProducer("FlashggEXOTagProducer",
                                inputTagJets= UnpackedJetCollectionVInputTag,
                                ElectronTag= cms.InputTag(electronString),
                                DiPhotonTag     = cms.InputTag("flashggDiPhotons"),
                                rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll')
                                )



process.exoTagDumper = createTagDumper("EXOTag")
process.exoTagDumper.dumpTrees = True
process.exoTagDumper.dumpHistos = False
process.exoTagDumper.dumpWorkspace = False

diphoton_vars = [

        "eventID := getEventNumber()",

        "diphoton_Mass := getDiphotonMass()",
        "diphoton_NConv := getDiphotonNConv()",
        "diphoton_PullCov := getDiphotonPullConv()",

        "diphoton_CutsPass := getDiphotonCutsPass()",

        "diphoton_LeadPt := getDiphotonLeadPt()",
        "diphoton_SubleadPt := getDiphotonSubleadPt()",
        "diphoton_LeadEta := getDiphotonLeadEta()",
        "diphoton_SubleadEta := getDiphotonSubleadEta()",

        "diphoton_LeadR9 := getDiphotonLeadR9()",
        "diphoton_SubleadR9 := getDiphotonSubleadR9()",
        "diphoton_LeadEtaSC := getDiphotonLeadEtaSC()",
        "diphoton_SubleadEtaSC := getDiphotonSubleadEtaSC()",
        "diphoton_LeadPhiSC := getDiphotonLeadPhiSC()",
        "diphoton_SubleadPhiSC := getDiphotonSubleadPhiSC()",
        "diphoton_Category := getDiphotonCategory()",

        "diphoton_LeadCHI := getDiphotonLeadCHI()",
        "diphoton_SubleadCHI := getDiphotonSubleadCHI()",
        "diphoton_LeadEGPhoIso := getDiphotonLeadEGPhoIso()",
        "diphoton_SubleadEGPhoIso := getDiphotonSubleadEGPhoIso()",
        "diphoton_LeadF5x5SigmaIetaIeta := getDiphotonLeadF5x5SigmaIetaIeta()",
        "diphoton_SubleadF5x5SigmaIetaIeta := getDiphotonSubleadF5x5SigmaIetaIeta()",
        "diphoton_LeadFull5x5R9 := getDiphotonLeadFull5x5R9()",
        "diphoton_SubleadFull5x5R9 := getDiphotonSubleadFull5x5R9()",
        "diphoton_LeadHadronicOverEM := getDiphotonLeadHadronicOverEM()",
        "diphoton_SubleadHadronicOverEM := getDiphotonSubleadHadronicOverEM()",

        "diphoton_LeadIsSaturated := getDiphotonLeadIsSaturated()",
        "diphoton_SubleadIsSaturated := getDiphotonSubleadIsSaturated()",
        "diphoton_LeadPassElectronVeto := getDiphotonLeadPassElectronVeto()",
        "diphoton_SubleadPassElectronVeto := getDiphotonSubleadPassElectronVeto()"

        ]

jet_vars = [

        "jets_multiplicites_All := getJetMultiplicities_All()",
        "jets_multiplicites_EGT20 := getJetMultiplicities_EGT20()",
        "jets_multiplicites_EGT30 := getJetMultiplicities_EGT30()",
        "jets_multiplicites_EGT40 := getJetMultiplicities_EGT40()",

        "dijet_LeadPt := getDijetLeadPt()",
        "dijet_SubleadPt := getDijetSubleadPt()",
        "dijet_LeadEta := getDijetLeadEta()",
        "dijet_SubleadEta := getDijetSubleadEta()",
        "dijet_Mass := getDijetMass()",
        "dijet_DeltaEta := getDijetDeltaEta()",
        "dijet_Zeppenfeld := getDijetZeppenfeld()",
        "dijet_DeltaPhi_jj := getDijetDeltaPhi_jj()",
        "dijet_DeltaPhi_ggjj := getDijetDeltaPhi_ggjj()"

        ]

electron_vars = [

 
        "electrons_multiplicites_EGT35 := getElectronMultiplicity_EGT35()",
        "electrons_multiplicites_EGT75 := getElectronMultiplicity_EGT75()",

        "dielectron_LeadPt := getDielectronLeadPt()",
        "dielectron_SubleadPt := getDielectronSubleadPt()",
        "dielectron_LeadEta := getDielectronLeadEta()",
        "dielectron_SubleadEta := getDielectronSubleadEta()",
        "dielectron_Mass := getDielectronMass()",
        "dielectron_DeltaEta := getDielectronDeltaEta()",
        "dielectron_Zeppenfeld := getDielectronZeppenfeld()",
        "dielectron_DeltaPhi_ee := getDielectronDeltaPhi_ee()",
        "dielectron_DeltaPhi_ggee := getDielectronDeltaPhi_ggee()"

        ]
        
all_var = diphoton_vars + jet_vars + electron_vars

cfgTools.addCategories(process.exoTagDumper,
                        [("test","getDiphotonMass()>0",0)],
                        variables = all_var,
                        histograms = []
                        )

process.exoTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"

from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",1.e+4)
customize(process)

process.p1 = cms.Path( process.flashggUnpackedJets+process.flashggEXOTag+process.exoTagDumper )

print process.p1
