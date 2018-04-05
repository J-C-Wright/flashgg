import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
import json
import sys
import os





process = cms.Process("FLASHggTag")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 500 )




'''
args = sys.argv

job_json = args[2].split('=')[-1]
out_path = args[3].split('=')[-1]
run_number = int(args[4].split('=')[-1])
print '\033[94m' + "Job JSON: " + '\033[0m',job_json
print '\033[94m' + "Out path: " + '\033[0m',out_path
print '\033[94m' + "File number: " + '\033[0m',run_number

with open(job_json,'r') as jf:
    job_info = json.loads(jf.read()) 

run_info = job_info['run_list'][run_number]
print '\033[94m' + 'Sample: ' + '\033[0m',run_info['sample']
print '\033[94m' + 'In path: ' + '\033[0m',run_info['path']
print '\033[94m' + 'N events: ' + '\033[0m',run_info['n_events']

file_path = run_info['path'].encode('ascii','ignore')
file_out = (out_path+'/%s_%s'%(run_info['sample'],run_info['path'].split('/')[-1])).encode('ascii','ignore')
print '\033[94m' + "Out file: " + '\033[0m',file_out

print '\033[94m' + 'XS etc.:' + '\033[0m',
with open('MetaData/data/cross_sections.json','r') as jf:
    xs_info = json.loads(jf.read())[run_info['sample']]
print xs_info
'''
xs_info = {}



file_path = "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_6-25ns_Moriond17/2_4_6/VBFHToGG_M125_13TeV_amcatnlo_pythia8/RunIISummer16-2_4_6-25ns_Moriond17-2_4_6-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/180207_080416/0000/myMicroAODOutputFile_5.root"
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(file_path))
process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/Taggers/flashggDJINNTreeMaker_cfi")


#Stuff for weighting
if 'xs' in xs_info.keys():
    xs = xs_info['xs']
else:
    xs = 1.0

if 'br' in xs_info.keys():
    br = xs_info['br']
else:
    br = 1.0

if 'kf' in xs_info.keys():
    kf = xs_info['kf']
else:
    kf = 1.0
print xs,br,kf,xs*br*kf

#Luminosity and XR*BR*KF weights
#process.flashggDJINNTreeMaker.lumiWeight = cms.double(1.0)
#process.flashggDJINNTreeMaker.xs = cms.double(xs*br*kf)

#RMS cut input
process.flashggDJINNTreeMaker.rmsforwardCut = cms.double(3.0)

#Pujid inputs
process.flashggDJINNTreeMaker.pujidWpPtBin1 = cms.vdouble([0.69, -0.35, -0.26, -0.21])
process.flashggDJINNTreeMaker.pujidWpPtBin2 = cms.vdouble([0.86, -0.1 , -0.05, -0.01])
process.flashggDJINNTreeMaker.pujidWpPtBin3 = cms.vdouble([0.95,  0.28,  0.31,  0.28])

#process.flashggDJINNTreeMaker.pujidWpPtBin1 = cms.vdouble([])
#process.flashggDJINNTreeMaker.pujidWpPtBin2 = cms.vdouble([])
#process.flashggDJINNTreeMaker.pujidWpPtBin3 = cms.vdouble([])

#Jet ID inputs
process.flashggDJINNTreeMaker.JetIDLevel = cms.string("Tight")
process.flashggDJINNTreeMaker.useJetID = cms.bool(True)

#Inputs for pileup weighting
process.flashggDJINNTreeMaker.dataPu = cms.vdouble( 239000.0, 838000.0, 2310000.0, 3120000.0, 4480000.0,
                                                    6000000.0, 7000000.0, 12900000.0, 35300000.0, 78700000.0,
                                                    177000000.0, 360000000.0, 603000000.0, 877000000.0, 1170000000.0,
                                                    1490000000.0, 1760000000.0, 1940000000.0, 2050000000.0, 2100000000.0,
                                                    2130000000.0, 2150000000.0, 2130000000.0, 2060000000.0, 1960000000.0,
                                                    1840000000.0, 1700000000.0, 1550000000.0, 1400000000.0, 1240000000.0,
                                                    1090000000.0, 937000000.0, 792000000.0, 657000000.0, 534000000.0,
                                                    427000000.0, 335000000.0, 258000000.0, 194000000.0, 142000000.0,
                                                    101000000.0, 69000000.0, 45500000.0, 28800000.0, 17500000.0,
                                                    10200000.0, 5640000.0, 2990000.0, 1510000.0, 732000.0,
                                                    340000.0, 153000.0, 67400.0, 30500.0, 15200.0,
                                                    8980.0, 6500.0, 5430.0, 4890.0, 4520.0,
                                                    4210.0, 3910.0, 3610.0, 3320.0, 3030.0,
                                                    2750.0, 2470.0, 2210.0, 1970.0, 1740.0,
                                                    1520.0, 1320.0, 1140.0, 983.0, 839.0 )

process.flashggDJINNTreeMaker.mcPu = cms.vdouble( 1.78653e-05, 2.56602e-05, 5.27857e-05, 8.88954e-05, 0.000109362,
                                                  0.000140973, 0.000240998, 0.00071209, 0.00130121, 0.00245255,
                                                  0.00502589, 0.00919534, 0.0146697, 0.0204126, 0.0267586,
                                                  0.0337697, 0.0401478, 0.0450159, 0.0490577, 0.0524855,
                                                  0.0548159, 0.0559937, 0.0554468, 0.0537687, 0.0512055,
                                                  0.0476713, 0.0435312, 0.0393107, 0.0349812, 0.0307413,
                                                  0.0272425, 0.0237115, 0.0208329, 0.0182459, 0.0160712,
                                                  0.0142498, 0.012804, 0.011571, 0.010547, 0.00959489,
                                                  0.00891718, 0.00829292, 0.0076195, 0.0069806, 0.0062025,
                                                  0.00546581, 0.00484127, 0.00407168, 0.00337681, 0.00269893,
                                                  0.00212473, 0.00160208, 0.00117884, 0.000859662, 0.000569085,
                                                  0.000365431, 0.000243565, 0.00015688, 9.88128e-05, 6.53783e-05,
                                                  3.73924e-05, 2.61382e-05, 2.0307e-05, 1.73032e-05, 1.435e-05,
                                                  1.36486e-05, 1.35555e-05, 1.37491e-05, 1.34255e-05, 1.33987e-05,
                                                  1.34061e-05, 1.34211e-05, 1.34177e-05, 1.32959e-05, 1.33287e-05 )

process.flashggDJINNTreeMaker.puBins = cms.vdouble( 0.0, 1.0, 2.0, 3.0, 4.0,
                                                    5.0, 6.0, 7.0, 8.0, 9.0,
                                                    10.0, 11.0, 12.0, 13.0, 14.0,
                                                    15.0, 16.0, 17.0, 18.0, 19.0,
                                                    20.0, 21.0, 22.0, 23.0, 24.0,
                                                    25.0, 26.0, 27.0, 28.0, 29.0,
                                                    30.0, 31.0, 32.0, 33.0, 34.0,
                                                    35.0, 36.0, 37.0, 38.0, 39.0,
                                                    40.0, 41.0, 42.0, 43.0, 44.0,
                                                    45.0, 46.0, 47.0, 48.0, 49.0,
                                                    50.0, 51.0, 52.0, 53.0, 54.0,
                                                    55.0, 56.0, 57.0, 58.0, 59.0,
                                                    60.0, 61.0, 62.0, 63.0, 64.0,
                                                    65.0, 66.0, 67.0, 68.0, 69.0,
                                                    70.0, 71.0, 72.0, 73.0, 74.0 )

process.flashggDJINNTreeMaker.puInfo = cms.InputTag("slimmedAddPileupInfo")

process.flashggDJINNTreeMaker.puReWeight = cms.bool(True)

process.flashggDJINNTreeMaker.rho = cms.InputTag("fixedGridRhoAll")

process.flashggDJINNTreeMaker.useTruePu = cms.bool(True)

process.flashggDJINNTreeMaker.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")

process.flashggDJINNTreeMaker.isData = cms.bool(False)

process.flashggDJINNTreeMaker.dijet_BDT_XML = cms.FileInPath("flashgg/Taggers/data/sklearn_training_moriond17_v8.xml")
process.flashggDJINNTreeMaker.combined_BDT_XML = cms.FileInPath("flashgg/Taggers/data/sklearn_combined_moriond17_v4.xml")
process.flashggDJINNTreeMaker.BDTMethod = cms.string("BDTG")

process.flashggDJINNTreeMaker.reweighGGHforNNLOPS = cms.untracked.bool(False)
process.flashggDJINNTreeMaker.NNLOPSWeightFile = cms.FileInPath("flashgg/Taggers/data/NNLOPS_reweight.root")




from flashgg.Taggers.globalVariables_cff import globalVariables
globalVariables.puReWeight = cms.bool(True)
globalVariables.dataPu = process.flashggDJINNTreeMaker.dataPu
globalVariables.puBins = process.flashggDJINNTreeMaker.puBins
globalVariables.mcPu = process.flashggDJINNTreeMaker.mcPu


process.flashggDJINNTreeMaker.globalVariables = globalVariables

process.flashggDJINNTreeMaker.calcPdfWeights = cms.bool(False)

#Tree output stuff
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('output.root'),
                                   closeFileFast = cms.untracked.bool(True))

# import flashgg customization
from flashgg.MetaData.JobConfig import customize
# set default options if needed
customize.setDefault("maxEvents",100)
customize.setDefault("targetLumi",10e+3)
# call the customization
customize(process)

print globalVariables
#Run
process.p = cms.Path(process.flashggTagSequence*process.flashggDJINNTreeMaker)

