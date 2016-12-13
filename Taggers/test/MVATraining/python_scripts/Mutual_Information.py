import ROOT
import os
import sys
import subprocess
import stat
from os import popen
from math import fabs,sqrt
import numpy as np
from shutil import copyfile
from optparse import OptionParser
from VariableStudies import VariableStudies as vs
import scipy.stats
import scipy.optimize


def get_options():
    parser = OptionParser()

    parser.add_option('--treesFile',
                      dest='treesFile',default='',
                      help = '''
                      file containing info on the trees
                      must be in data/
                      ''',
                      metavar='FILE')
    parser.add_option('--varsFile',
                      dest='varsFile',default='',
                      help = '''
                      file containing var info
                      must be in data/
                      ''',
                      metavar = 'FILE')
    parser.add_option('--cut',
                      dest='cut_name',default='vbf_ps',
                      help='''
                      Specify cuts to apply:
                      vbf_ps, dipho_2j, vbf_ps_3j, or dipho_3j
                      ''')
    parser.add_option('--signalRegion',
                      dest='signal_region',action='store_true',
                      default=False,help='''
                      Only look at events with diphoton mass
                      between 115 and 135
                      ''')
    parser.add_option('--sample',
                      dest='sample',default='',
                      help='''
                      specify which sample to use: vbf, ggH,
                      gjet,dpbox,qcd,allbkg
                      ''')

    return parser.parse_args()




if __name__ == '__main__':

    (opt,args) = get_options()

    cut_name = opt.cut_name
    varsInfoFile = opt.varsFile
    treesFile = opt.treesFile
    signal_region = opt.signal_region

    if opt.sample == 'qcd':
        num_slices = 20
    else:
        num_slices = 100

    outpath = '/afs/cern.ch/user/j/jwright/www/'+opt.sample+'_corr_'+cut_name+'_'+varsInfoFile.split('.')[0]+'/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    copyfile('/afs/cern.ch/user/j/jwright/www/index.php',outpath+'index.php')

    tableName = opt.sample+'_corr_'+cut_name+'_'+varsInfoFile.split('.')[0]

    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    vbf_info, ggh_info, gjet_info, dbox_info, qcd_info, _ = vs.getInfoFromFile(name=treesFile)

    vbf_tfiles = []
    ggh_tfiles = []
    gjet_tfiles = []
    dbox_tfiles = []
    qcd_tfiles = []

    for filename in vbf_info[0]:
        vbf_tfiles.append(ROOT.TFile.Open(filename))
    for filename in ggh_info[0]:
        ggh_tfiles.append(ROOT.TFile.Open(filename))
    for filename in gjet_info[0]:
        gjet_tfiles.append(ROOT.TFile.Open(filename))
    for filename in dbox_info[0]:
        dbox_tfiles.append(ROOT.TFile.Open(filename))
    for filename in qcd_info[0]:
        qcd_tfiles.append(ROOT.TFile.Open(filename))

    vbf_trees = vs.getTrees(tfiles=vbf_tfiles,treenames=vbf_info[1])
    ggh_trees = vs.getTrees(tfiles=ggh_tfiles,treenames=ggh_info[1])
    gjet_trees = vs.getTrees(tfiles=gjet_tfiles,treenames=gjet_info[1])
    dbox_trees = vs.getTrees(tfiles=dbox_tfiles,treenames=dbox_info[1])
    qcd_trees = vs.getTrees(tfiles=qcd_tfiles,treenames=qcd_info[1])

    histos_info = vs.getHistosInfoFromFile(name=varsInfoFile)

    cuts = []
    if cut_name == 'vbf_ps':
        cuts = [vs.dph_ps_cut, vs.vbf_ps_cut]
    elif cut_name == 'dipho_2j':
        cuts = [vs.dph_ps_cut]
    elif cut_name == 'vbf_ps_3j':
        cuts = [vs.dph_ps_cut, vs.vbf_ps_cut, vs.jet3_cut]
    elif cut_name == 'dipho_3j':
        cuts = [vs.dph_ps_cut, vs.jet3_cut]

    if signal_region:
        cuts.append(vs.signal_region)

    cut = '&&'.join(cuts)

    trees = []
    if opt.sample == 'vbf':
        trees = vbf_trees
    elif opt.sample == 'ggH':
        trees = ggh_trees
    elif opt.sample == 'gjet':
        trees = gjet_trees
    elif opt.sample == 'dpbox':
        trees = dbox_trees
    elif opt.sample == 'qcd':
        trees = qcd_trees
    elif opt.sample == 'allbkg':
        trees = ggh_trees + gjet_trees + dbox_trees + qcd_trees

    x_info = histos_info[0]
    y_info = histos_info[1]

    c1 = ROOT.TCanvas('c1')
    
    test_hist = ROOT.TH2F(x_info[0]+"_vs_"+y_info[0],"",int(x_info[1]),float(x_info[2]),float(x_info[3]),
                                                        int(y_info[1]),float(y_info[2]),float(y_info[3]))
    vbf_trees[0].Draw("%s:%s"%(x_info[0],y_info[0]))
    c1.Print('test.pdf')

    #for tree in vbf_trees:
        #test_hist = ROOT.TH2F(x_info[0]+"_vs_"+y_info[0],"",int(x_info[1]),float(x_info[2]),float(x_info[3]),
        #                                                    int(y_info[1]),float(y_info[2]),float(y_info[3]))


















