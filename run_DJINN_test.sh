#!/bin/bash

jsonfile=/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/flashgg/run_samples_signals_and_data.json

today=`date +%F`
outdir=/vols/cms/jwright/DJINN_configtest/test_job_config-${today}
fggRunJobs.py --load ${jsonfile} -d ${outdir} \
      -x cmsRun Taggers/test/djinn_treemaker.py maxEvents=-1 \
      -q hepmedium.q --no-use-tarball useAAA=1 targetLumi=35.90e+3 \
      -n 200

