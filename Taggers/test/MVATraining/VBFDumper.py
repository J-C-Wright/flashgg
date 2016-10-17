#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager
import sys

process = cms.Process('VBFDumper')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(500)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #"/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/160707_142718/0000/myMicroAODOutputFile_100.root"
                                "/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/VBFHToGG_M-125_13TeV_powheg_pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/160707_150107/0000/myMicroAODOutputFile_10.root"
                                )
                            )

process.TFileService = cms.Service( "TFileService",
                                    fileName = cms.string("VBFDumperTree.root"),
                                    closeFileFast = cms.untracked.bool(True) )

from flashgg.Taggers.flashggTagOutputCommands_cff import tagDefaultOutputCommand
import flashgg.Taggers.dumperConfigTools as cfgTools
from flashgg.Taggers.tagsDumpers_cfi import createTagDumper

process.flashggUnpackedJets = cms.EDProducer( "FlashggVectorVectorJetUnpacker",
                                              JetsTag = cms.InputTag("flashggFinalJets"),
                                              NCollections = cms.uint32(8) )

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
                                    

process.load('flashgg.Taggers.flashggTagSequence_cfi')
process.load('flashgg.Taggers.flashggTagTester_cfi')

process.vbfTagDumper = createTagDumper("VBFTag")
process.vbfTagDumper.dumpTrees = True
process.vbfTagDumper.dumpHistos    = False
process.vbfTagDumper.dumpWorkspace = False


import flashgg.Taggers.VBFTagVariables as var

all_var = var.variable_exploration

cfgTools.addCategories(process.vbfTagDumper,
                        [("test","diPhoton.mass>0",0)],
                        variables = all_var,
                        histograms = []
                        )

process.vbfTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"

from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents" ,-1)    # max-number of events
customize.setDefault("targetLumi",2.11e+3) # define integrated lumi
customize(process)

process.p1 = cms.Path(
                    process.flashggTagSequence+
                    process.vbfTagDumper
                    )

print process.p1



