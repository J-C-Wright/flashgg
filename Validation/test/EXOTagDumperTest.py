# ================================================
#   MicroAOD and Jet Tree Producer
#   Y. Haddad 01/2015
# ================================================

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
process = cms.Process("FLASHggEXOValidation")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("flashgg.Systematics.flashggDiPhotonSystematics_cfi")
process.load("flashgg.Systematics.flashggMuonSystematics_cfi")
process.load("flashgg.Systematics.flashggElectronSystematics_cfi")
process.load("flashgg.Systematics.flashggJetSystematics_cfi")
process.flashggDiPhotonSystematics.src='flashggDiPhotons'

# +++++ the number of processed events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

# +++++ the source file
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
                                                        #"/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall15DR76-1_3_0-25ns/1_3_0/DoubleEG/RunIIFall15DR76-1_3_0-25ns-1_3_0-v0-Run2015D-16Dec2015-v2/160116_222234/0000/myMicroAODOutputFile_720.root"
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_393.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_392.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_391.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_390.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_389.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_388.root",
														"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_387.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_386.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_385.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_384.root",
														#"/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_383.root"
														)) 

#process.MessageLogger.cerr.threshold = 'ERROR'

process.TFileService = cms.Service("TFileService",fileName  = cms.string("exoValidationTrees_test.root"))

# process.flashggEXOValidationTreeMaker = cms.EDAnalyzer('FlashggEXOValidationTreeMaker',
#                                                        GenParticleTag           = cms.untracked.InputTag('prunedGenParticles'),
#                                                        JetTagDz                 = cms.InputTag("flashggJets"),
#                                                        StringTag		= cms.string("PF"),
#                                                    )
#

#JetCollectionVInputTagPFCHS = cms.VInputTag()
#JetCollectionVInputTagPUPPI = cms.VInputTag()
#for i in range(0,5):
#    JetCollectionVInputTagPFCHS.append(cms.InputTag('flashggPFCHSJets' + str(i)))
    #JetCollectionVInputTagPUPPI.append(cms.InputTag('flashggPUPPIJets' + str(i)))

process.flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
                                             JetsTag = cms.InputTag("flashggFinalJets"),
                                             NCollections = cms.uint32(8))
	 
hltPaths = ["HLT_DoublePhoton85*","HLT_Photon250_NoHE*","HLT_DoublePhoton60*"]
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.HLTPaths = hltPaths

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
process.flashggEXOValidationTreeMaker = cms.EDAnalyzer('FlashggEXOValidationTreeMaker',
                                                            inputTagJets= UnpackedJetCollectionVInputTag,
                                                            ElectronTag= cms.InputTag("flashggSelectedElectrons"),
#                                                            ElectronTag= cms.InputTag("flashggElectrons"),
                                                            DiPhotonTag     = cms.InputTag("flashggDiPhotonSystematics"),
                                                            rhoFixedGridCollection = cms.InputTag('fixedGridRhoAll')
                                                            )


from flashgg.Taggers.tagDumpers_cfi import createTagDumper
import flashgg.Taggers.dumperConfigTools as cfgTools

process.exoTagDumper = createTagDumper("EXOTag")
process.exoTagDumper.dumpTrees = True

diphotonVariables = ["eventNumber := EXOTag().getEventNumber()",
                     "diphotonMass := EXOTag().getDiphotonMass()",
                     "diphotonLeadPt := EXOTag().getDiphotonLeadPt()",
                     "diphotonSubleadPt := EXOTag().getDiphotonSubleadPt()",
                     "diphotonLeadEta := EXOTag().getDiphotonLeadEta()",
                     "diphotonSubleadEta := EXOTag().getDiphotonSubleadEta()"]

cfgTools.addCategories( process.exoTagDumper,
                        [
                            ("EXODiphoton","EXOTag().getDiphotonMass() > 0") #Placeholder
                        ],
                        variables = diphotonVariables
                        histograms = [];
    )

process.vbfTagDumper.nameTemplate = "$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL"

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.p = cms.Path( 
    process.flashggUnpackedJets+
    process.hltHighLevel*
    process.flashggDiPhotonSystematics*
    process.flashggEXOValidationTreeMaker 
    )
process.e = cms.EndPath()

process.load('flashgg.Systematics.escales.escale76X_16DecRereco_2015')
#useEGMTools(process)

 # customization for job splitting, lumi weighting, etc.
from diphotons.MetaData.JobConfig import customize
customize.setDefault("maxEvents",3000)
customize.setDefault("processId",'Data')
customize(process)
from flashgg.Systematics.SystematicsCustomize import *
#if customize.processId == "Data":	
customizePhotonSystematicsForData(process)
#customizeSystematicsForSignal(process)


