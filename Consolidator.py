
import json
import os

'''
base_dir = '/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/'

#catalogue_path = '/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/flashgg/MetaData//data/RunIISummer16-2_4_5-25ns_Moriond17/'
#catalogue_names = ['datasets.json','datasets_1.json','datasets_2.json','datasets_3.json','datasets_4.json','datasets_5.json']
catalogue_path = '/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/flashgg/MetaData/data/RunIISummer16-2_4_6-25ns_Moriond17/'
catalogue_names = ['datasets.json','datasets_1.json','datasets_2.json','datasets_3.json','datasets_4.json']

keys = []
for name in catalogue_names:

    with open(catalogue_path+name,'r') as cat:
        catalogue = json.loads(cat.read())

    keys += [key.split('/')[1] for key in catalogue.keys()]

keys = [
        'VBFHToGG_M125_13TeV_amcatnlo_pythia8',
        'GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8',
        'ZHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8'
        ]
keys = [
        'WHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8'
        ]
#keys = [
#        'ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8'
#        ]
print keys
#trees_dir = '/vols/cms/jwright/DJINN_Trees/'
#trees_dir = '/vols/cms/jwright/DJINN_configtest/test_job_config-2018-03-05/'
#trees_dir = '/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-tth-2018-03-14/'
'''




trees_dirs = ['/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-2018-03-15/',
             '/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-tth-2018-03-16/',
             '/vols/cms/jwright/DJINN_configtest/test_sig_syst_config-wh-2018-03-16/']

out_dir = '/vols/cms/jwright/DJINN_Trees/systematics_hadded/'

keys = []
for trees_dir in trees_dirs:
    for file in os.listdir(trees_dir):
            if file.endswith(".root"):
                name = '_'.join(file.split('_')[:-1])
                if name not in keys:
                    keys.append(trees_dir+name)

keys = list(set(sorted(keys)))
#for name in keys:
#    print name

for sample in keys:

    print '------'
    
    command = 'hadd -f '
    command += '%s%s.root '%(out_dir+'/',sample.split('/')[-1])
    command += '%s*.root'%(sample)
    print command

    hadd_message = os.popen(command).read()
    print hadd_message
'''
'''


