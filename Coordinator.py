
import json
import os

def make_job_filelists(file_list,num_jobs):
    events_per_job = num_events/num_jobs
    jobs = []
    i = 0
    exit_state = False
    while not exit_state:
        count = 0
        files_in_job = []
        while count < events_per_job:

            if i == len(file_list):
                break

            files_in_job.append(file_list[i])
            count += file_list[i]['n_events']
            i += 1
        jobs.append(files_in_job)
        if i==len(file_list):
            exit_state = True

    return jobs

base_dir = '/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/'

catalogue_path = '/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/flashgg/MetaData//data/RunIISummer16-2_4_5-25ns_Moriond17/'
catalogue_names = ['datasets.json','datasets_1.json','datasets_2.json','datasets_3.json','datasets_4.json','datasets_5.json']

num_jobs = 1000
num_files = 0
num_events = 0

print "Dataset summary:"
file_list = []
for name in catalogue_names:

    with open(catalogue_path+name,'r') as cat:
        catalogue = json.loads(cat.read())

    for key,value in catalogue.iteritems():

        for f in value['files']:
            if f['bad']: continue
    
            file_info = {}
            file_info['sample'] = key.split('/')[1]
            file_info['path'] = f['name']
            file_info['n_events'] = f['nevents']
            file_list.append(file_info)

            num_events += f['nevents']
            num_files += 1

        print '\t',key.split('/')[1]

print "\tTotal number of files: %d"%num_files
print "\tTotal number of events: %d"%num_events

jobs = make_job_filelists(file_list,num_jobs)
script_paths = []
for i,run_list in enumerate(jobs):
    job_info = {}
    job_info['name'] = 'job_%d'%i
    job_info['run_list'] = run_list
    with open('jobs/'+job_info['name']+'.json','w') as jf:
        jf.write(json.dumps(job_info))

    script = ''
    script += '#!/bin/sh' + '\n'
    script += 'export X509_USER_PROXY=/home/hep/jw3914/x509up_u634344' + '\n'
    script += 'cd %s\n'%base_dir
    script += 'eval $(scram runtime -sh)' + '\n'
    script += 'cd flashgg/' + '\n'
    script += 'ls $X509_USER_PROXY' + '\n'

    for j in range(len(run_list)):
        script += 'cmsRun Taggers/test/djinn_treemaker.py jobJson=jobs/job_%d.json outPath=/vols/cms/jwright/DJINN_Trees/ runNumber=%d'%(i,j)
        script += '\n'
    script_path = 'jobs/scripts/'+job_info['name']+'.sh'

    with open(script_path,'w') as jf:
        jf.write(script)
    os.chmod(script_path, 0755)
    script_paths.append(script_path)

#for path in script_paths:
#    qsub_string =  'qsub -q hep.q -l h_rt=6:0:0 -l h_vmem=6G ./%s'%path
#    print qsub_string 
#    qsub_message = os.popen(qsub_string).read()
#    print qsub_message




