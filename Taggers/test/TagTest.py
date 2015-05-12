import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggTag")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

#Setting max events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#Old MC Samples
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_1.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_10.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_11.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_12.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_2.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_3.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_4.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_5.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_6.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_7.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_8.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_155949/0000/myMicroAODOutputFile_9.root"))

#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_1.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_10.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_11.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_2.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_3.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_4.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_5.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_6.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_7.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_8.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6/HggPhys14-Phys14MicroAODV2-v2-Phys14DR-PU20bx25_PHYS14_25_V1-v2/150213_054817/0000/myMicroAODOutputFile_9.root"))

#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/TTbarH_HToGG_M-125_13TeV_amcatnlo-pythia8-tauola/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_160101/0000/myMicroAODOutputFile_1.root"))

#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/GluGluToHToGG_M-125_13TeV-powheg-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150210_160020/0000/myMicroAODOutputFile_1.root"))

#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/VBF_HToGG_M-125_13TeV-powheg-pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_160130/0000/myMicroAODOutputFile_1.root"))

#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV2/WH_ZH_HToGG_M-125_13TeV_pythia6/HggPhys14-Phys14MicroAODV2-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150210_160158/0000/myMicroAODOutputFile_1.root"))

#New (15/4/15) MC sample
#ggH
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/GluGluToHToGG_M-125_13TeV-powheg-pythia6/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/150414_122024/0000/myMicroAODOutputFile_1.root"))
#ttH
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/TTbarH_HToGG_M-125_13TeV_amcatnlo-pythia8-tauola/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150413_121421/0000/myMicroAODOutputFile_1.root"))
#VBF
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/VBF_HToGG_M-125_13TeV-powheg-pythia6/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150413_121454/0000/myMicroAODOutputFile_1.root"))
#VH
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/WH_ZH_HToGG_M-125_13TeV_pythia6/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU20bx25_PHYS14_25_V1-v1/150413_121525/0000/myMicroAODOutputFile_1.root"))

#Drell-Yan
#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/DYJetsToLL_M-50_13TeV-madgraph-pythia8/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU4bx50_PHYS14_25_V1-v1/150413_121142/0000/myMicroAODOutputFile_1.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/DYJetsToLL_M-50_13TeV-madgraph-pythia8/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU4bx50_PHYS14_25_V1-v1/150413_121142/0000/myMicroAODOutputFile_10.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/DYJetsToLL_M-50_13TeV-madgraph-pythia8/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU4bx50_PHYS14_25_V1-v1/150413_121142/0000/myMicroAODOutputFile_11.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/DYJetsToLL_M-50_13TeV-madgraph-pythia8/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU4bx50_PHYS14_25_V1-v1/150413_121142/0000/myMicroAODOutputFile_2.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/DYJetsToLL_M-50_13TeV-madgraph-pythia8/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU4bx50_PHYS14_25_V1-v1/150413_121142/0000/myMicroAODOutputFile_6.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/DYJetsToLL_M-50_13TeV-madgraph-pythia8/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU4bx50_PHYS14_25_V1-v1/150413_121142/0000/myMicroAODOutputFile_8.root","/store/group/phys_higgs/cmshgg/sethzenz/flashgg/HggPhys14/Phys14MicroAODV3/DYJetsToLL_M-50_13TeV-madgraph-pythia8/HggPhys14-Phys14MicroAODV3-v0-Phys14DR-PU4bx50_PHYS14_25_V1-v1/150413_121142/0000/myMicroAODOutputFile_9.root"))



process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/Taggers/flashggTagTester_cfi")

from flashgg.Taggers.flashggTagOutputCommands_cff import tagDefaultOutputCommand

process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myTagOutputFile.root'),
                               outputCommands = tagDefaultOutputCommand			       
                               )

process.p = cms.Path(process.flashggTagSequence*process.flashggTagTester)
process.e = cms.EndPath(process.out)
