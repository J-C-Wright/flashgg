
import json
import os

trees_dirs = [
                #'/vols/cms/jwright/DJINN_configtest/train_job_config-2018-04-03/'
                '/vols/cms/jwright/DJINN_configtest/train_job_config-2018-04-04/'
              ]

out_dir = '/vols/cms/jwright/DJINN_Trees/systematics_hadded/train/'

keys = []
for trees_dir in trees_dirs:
    for file in os.listdir(trees_dir):
            if file.endswith(".root"):
                name = '_'.join(file.split('_')[:-1])
                if name not in keys:
                    keys.append(trees_dir+name)

keys = list(set(sorted(keys)))


for sample in keys:

    print '------'
    
    command = 'hadd -f '
    command += '%s%s.root '%(out_dir+'/',sample.split('/')[-1])
    command += '%s*.root'%(sample)
    print command

    hadd_message = os.popen(command).read()
    print hadd_message

