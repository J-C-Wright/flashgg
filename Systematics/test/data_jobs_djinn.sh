# NB this command is specific to the configuration on lxplus and is not gaurenteed elsewhere
#outdir="/afs/cern.ch/work/s/sethzenz/ws/" # can't set absolute path on lsf because we're expecting to stage
queue="hepmedium.q"
useAAA=1
today=`date +%F`
outdir=/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-${today}
LM=${CMSSW_BASE}/src/flashgg/MetaData/work/jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt

fggRunJobs.py --load data_jobs.json -d ${outdir} -x cmsRun workspaceStd.py maxEvents=-1 -n 1000 -q ${queue} -D -P useAAA=$useAAA doFiducial=False tthTagsOnly=False lumiMask=${LM}

