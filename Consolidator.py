
import json
import os

trees_dirs = [#'/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-2018-03-15/',
              #'/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-tth-2018-03-16/',
              #'/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-wh-2018-03-16/',
              #'/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-tth-2018-03-19/']
              '/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-2018-03-20/']

out_dir = '/vols/cms/jwright/DJINN_Trees/systematics_hadded/'

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


