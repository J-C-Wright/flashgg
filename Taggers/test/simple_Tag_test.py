import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggTag")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )

# Uncomment the following if you notice you have a memory leak
# This is a lightweight tool to digg further
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#                                        ignoreTotal = cms.untracked.int32(1),
#                                        monitorPssAndPrivate = cms.untracked.bool(True)
#                                       )

#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring("file:myMicroAODOutputFile.root"))
#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-50ns/Spring15BetaV4/VBFHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-50ns-Spring15BetaV4-v0-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/150812_181849/0000/myMicroAODOutputFile_18.root"))
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-25ns/Spring15BetaV5/VBFHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-25ns-Spring15BetaV5-v0-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150922_093659/0000/myMicroAODOutputFile_11.root"))


process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/Taggers/flashggTagTester_cfi")

# For debugging
switchToPreselected = False
switchToFinal = False
switchToPuppi = False
assert(not switchToPreselected or not switchToFinal)

if switchToPreselected:
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggDiPhotons"),cms.InputTag("flashggPreselectedDiPhotons"))

if switchToFinal:
    from flashgg.MicroAOD.flashggFinalEGamma_cfi import flashggFinalEGamma
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.flashggTagSequence,cms.InputTag("flashggDiPhotons"),cms.InputTag("flashggFinalEGamma",flashggFinalEGamma.DiPhotonCollectionName.value()))

if switchToPuppi:
    process.flashggUnpackedJets.JetsTag = cms.InputTag("flashggFinalPuppiJets")

from flashgg.Taggers.flashggTagOutputCommands_cff import tagDefaultOutputCommand

process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myTagOutputFile.root'),
                               outputCommands = tagDefaultOutputCommand			       
                               )

process.p = cms.Path(process.flashggTagSequence*process.flashggTagTester)

process.e = cms.EndPath(process.out)
