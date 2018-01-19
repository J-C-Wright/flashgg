import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
import json
import sys


process = cms.Process("FLASHggTag")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(250) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 500 )

#sample_name = 
#file_path = "root://gfe02.grid.hep.ph.ic.ac.uk:1095//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_5-25ns_Moriond17/2_4_5/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16-2_4_5-25ns_Moriond17-2_4_5-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/171115_091126/0000/myMicroAODOutputFile_103.root"

args = sys.argv

job_json = args[2].split('=')[-1]
out_path = args[3].split('=')[-1]
run_number = int(args[4].split('=')[-1])
print job_json
print out_path
print run_number

with open(job_json,'r') as jf:
    job_info = json.loads(jf.read()) 

run_info = job_info['run_list'][run_number]
print run_info
file_path = run_info['path'].encode('ascii','ignore')
file_out = (out_path+'/%s_%s'%(run_info['sample'],run_info['path'].split('/')[-1])).encode('ascii','ignore')
print ''
print file_out

process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(file_path))

#file_list = ['file:ForJack_VBF_1k_10Nov.root','file:myMicroAODOutputFile_103.root']
#process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(file_list))



process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/Taggers/flashggDJINNTreeMaker_cfi")

process.flashggDJINNTreeMaker.lumiWeight = cms.double(0.5)
process.flashggDJINNTreeMaker.xs = cms.double(0.5)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(file_out),
                                   closeFileFast = cms.untracked.bool(True))

process.p = cms.Path(process.flashggTagSequence*process.flashggDJINNTreeMaker)

