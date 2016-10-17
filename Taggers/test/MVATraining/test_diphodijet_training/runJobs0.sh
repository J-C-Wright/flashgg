#!/bin/bash

WD=$PWD
echo
echo
echo
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
scram project CMSSW CMSSW_8_0_8_patch1
cd CMSSW_8_0_8_patch1
tar zxf /afs/cern.ch/work/j/jwright/private/VBF/CMSSW_8_0_8_patch1/src/flashgg/Taggers/test/MVATraining/test_diphodijet_training/sandbox.tgz
scram b
eval $(scram runtime -sh)
cd $WD
mkdir test_diphodijet_training
echo "ls $X509_USER_PROXY"
ls $X509_USER_PROXY
cmsRun /afs/cern.ch/work/j/jwright/private/VBF/CMSSW_8_0_8_patch1/src/flashgg/Taggers/test/MVATraining/test_diphodijet_training/VBFDiPhoDiJetMVA_Training.py maxEvents=100 campaign=RunIISpring16DR80X-2_2_0-25ns_DYWithPDF processIdMap=/afs/cern.ch/work/j/jwright/private/VBF/CMSSW_8_0_8_patch1/src/flashgg/Taggers/test/MVATraining/test_diphodijet_training/config.json dataset=/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 outputFile=test_diphodijet_training/output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root nJobs=1 jobId=0
retval=$?
cd test_diphodijet_training
echo
echo
echo "Job finished with exit code $?"
echo "Files in ouput folder"
ls -ltr
if [[ $retval == 0 ]]; then
    errors=""
    for file in $(find -name '*.root' -or -name '*.xml'); do
        cp -pv $file /afs/cern.ch/work/j/jwright/private/VBF/CMSSW_8_0_8_patch1/src/flashgg/Taggers/test/MVATraining/test_diphodijet_training
        if [[ $? != 0 ]]; then
            errors="$errors $file($?)"
        fi
    done
    if [[ -n "$errors" ]]; then
       echo "Errors while staging files"
       echo "$errors"
       exit -2
    fi
fi

exit $retval

