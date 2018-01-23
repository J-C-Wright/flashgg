
import json
import os

base_dir = '/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/'

catalogue_path = '/home/hep/jw3914/Work/DJINN_Data/CMSSW_8_0_28/src/flashgg/MetaData//data/RunIISummer16-2_4_5-25ns_Moriond17/'
catalogue_names = ['datasets.json','datasets_1.json','datasets_2.json','datasets_3.json','datasets_4.json','datasets_5.json']

keys = []
for name in catalogue_names:

    with open(catalogue_path+name,'r') as cat:
        catalogue = json.loads(cat.read())

    keys += [key.split('/')[1] for key in catalogue.keys()]

trees_dir = '/vols/cms/jwright/DJINN_Trees/'
keys = sorted(list(set(keys)))

for sample in keys:

    print '------'
    print sample
    
    command = 'hadd -f '
    command += '%s%s.root '%(trees_dir+'hadded/',sample)
    command += '%s%s*.root'%(trees_dir,sample)
    print command

    hadd_message = os.popen(command).read()
    print hadd_message




