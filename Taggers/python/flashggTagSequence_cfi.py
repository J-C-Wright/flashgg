import FWCore.ParameterSet.Config as cms
from flashgg.Taggers.flashggDiPhotonMVA_cfi import flashggDiPhotonMVA
from flashgg.Taggers.flashggVBFMVA_cff import flashggVBFMVA,flashggVBFMVANew,flashggVBFDiPhoDiJetMVA, flashggVBFDiPhoDiJetMVANew, flashggVBFJetFilter
from flashgg.Taggers.flashggTags_cff import *
from flashgg.Taggers.flashggTagSorter_cfi import flashggTagSorter, flashggTagSorterNew

flashggTagSequence = cms.Sequence(flashggDiPhotonMVA
																	* flashggVBFJetFilter
                                  * flashggVBFMVA
                                  * flashggVBFMVANew
                                  * flashggVBFDiPhoDiJetMVA
                                  * flashggVBFDiPhoDiJetMVANew
                                  * (flashggUntaggedCategory
                                     + flashggVBFTag
                                     + flashggVBFTagNew
                                     + flashggTTHLeptonicTag
                                     + flashggTTHHadronicTag
                                     + flashggVHLooseTag
                                     + flashggVHTightTag
                                     + flashggVHHadronicTag)
                                  * flashggTagSorter
                                  * flashggTagSorterNew
                                  )

