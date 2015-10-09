import FWCore.ParameterSet.Config as cms
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

flashggTagTester = cms.EDAnalyzer('FlashggTagTestAnalyzer',
                                  TagSorter = cms.InputTag('flashggTagSorter'),
                                  inputTagJets= UnpackedJetCollectionVInputTag,
                                  DiPhotonTag = cms.InputTag('flashggDiPhotons'),
                                  )
