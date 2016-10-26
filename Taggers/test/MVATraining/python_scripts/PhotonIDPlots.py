import ROOT
import os
import subprocess
import stat
from os import popen
from math import fabs
from shutil import copyfile
from optparse import OptionParser
from VariableStudies import VariableStudies as vs
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

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
                      dest='cut_name',default='none',
                      help='''
                      Specify cuts to apply:
                      vbf_ps, dipho_2j, vbf_ps_3j, dipho_3j, or none
                      ''')
    parser.add_option('--signalRegion',
                      dest='signal_region',action='store_true',
                      default=False,help='''
                      Only look at events with diphoton mass
                      between 115 and 135
                      ''')

    return parser.parse_args()




if __name__ == '__main__':

    (opt,args) = get_options()

    cut_name = opt.cut_name
    varsInfoFile = opt.varsFile
    treesFile = opt.treesFile
    signal_region = opt.signal_region


    outpath = '/afs/cern.ch/user/j/jwright/www/photon_id_study_'+cut_name+'_'+varsInfoFile.split('.')[0]+'/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    copyfile('/afs/cern.ch/user/j/jwright/www/index.php',outpath+'index.php')

    tableName = 'Correlation_'+cut_name+'_'+varsInfoFile.split('.')[0]

    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    ROOT.gStyle.SetOptStat(0)

    cuts = []
    if cut_name == 'vbf_ps':
        cuts = [vs.dph_ps_cut, vs.vbf_ps_cut]
    elif cut_name == 'dipho_2j':
        cuts = [vs.dph_ps_cut]
    elif cut_name == 'vbf_ps_3j':
        cuts = [vs.dph_ps_cut, vs.vbf_ps_cut, vs.jet3_cut]
    elif cut_name == 'dipho_3j':
        cuts = [vs.dph_ps_cut, vs.jet3_cut]
    elif cut_name == 'none':
        cuts = ['1']

    if signal_region:
        cuts.append(vs.signal_region)

    cut = '&&'.join(cuts)


    print "Table and plots can be found at:"
    print "http://jwright.web.cern.ch/jwright/"+outpath.split('www/')[-1]


    vbf_info, ggh_info, gjet_info, dbox_info, qcd_info, _ = vs.getInfoFromFile(name=treesFile)

    vbf_tfiles = []
    ggh_tfiles = []
    gjet_tfiles = []
    dbox_tfiles = []
    qcd_tfiles = []

    for filepath in vbf_info[0]:
        vbf_tfiles.append(ROOT.TFile.Open(filepath))
    for filepath in ggh_info[0]:
        ggh_tfiles.append(ROOT.TFile.Open(filepath))
    for filepath in gjet_info[0]:
        gjet_tfiles.append(ROOT.TFile.Open(filepath))
    for filepath in dbox_info[0]:
        dbox_tfiles.append(ROOT.TFile.Open(filepath))
    for filepath in qcd_info[0]:
        qcd_tfiles.append(ROOT.TFile.Open(filepath))

    vbf_trees = vs.getTrees(tfiles=vbf_tfiles,treenames=vbf_info[1])
    ggh_trees = vs.getTrees(tfiles=ggh_tfiles,treenames=ggh_info[1])
    gjet_trees = vs.getTrees(tfiles=gjet_tfiles,treenames=gjet_info[1])
    dbox_trees = vs.getTrees(tfiles=dbox_tfiles,treenames=dbox_info[1])
    qcd_trees = vs.getTrees(tfiles=qcd_tfiles,treenames=qcd_info[1])

    all_bkg_trees = gjet_trees+dbox_trees+qcd_trees

    histos_info = vs.getHistosInfoFromFile(name=varsInfoFile)

    print '\nvbf: ',
    for tree in vbf_trees:
        print tree.GetEntries(),
    print '\nggh: ',
    for tree in ggh_trees:
        print tree.GetEntries(),
    print '\ngjet: ',
    for tree in gjet_trees:
        print tree.GetEntries(),
    print '\ndbox: ',
    for tree in dbox_trees:
        print tree.GetEntries(),
    print '\nqcd: ',
    for tree in qcd_trees:
        print tree.GetEntries(),
    print
    for histo_info in histos_info:
        print histo_info

    tree_names = ['vbf','ggH','GJet','dpbox','QCD','all_bkg']
    trees_list = [vbf_trees,ggh_trees,gjet_trees,dbox_trees,qcd_trees,all_bkg_trees]

    x_info = histos_info[0]
    y_info = histos_info[1]

    c1 = ROOT.TCanvas('c1',"",500,500)
    c1.SetFixedAspectRatio()

    #2D Histograms
    histograms=[]
    for tree_name,trees in zip(tree_names,trees_list):
    
        histogram = vs.get2DHistogram(trees=trees,x_info=x_info,y_info=y_info,cut=cut)
        histogram.Scale(1.0/histogram.Integral())
        histograms.append(histogram)
        histogram.Draw('colz')
        c1.Print(outpath+tree_name+'.pdf')
        c1.Print(outpath+tree_name+'.png')

    #Mass plots
    id_region_cuts = ['(dipho_leadIDMVA > -0.9 && dipho_subleadIDMVA > -0.9)',
                      '((dipho_leadIDMVA < -0.9 && dipho_subleadIDMVA > -0.9) || (dipho_leadIDMVA > -0.9 && dipho_subleadIDMVA < -0.9))',
                      '(dipho_leadIDMVA < -0.9 && dipho_subleadIDMVA < -0.9)']
    id_region_names = ['gg','fg','ff']

    for tree_name,trees in zip(tree_names,trees_list):

        temp_hists = []
        for id_region_name,id_region_cut in zip(id_region_names,id_region_cuts):
            id_cut = '&&'.join([cut,id_region_cut])

            histogram = vs.get1DHistogram(trees=trees,x_info=histos_info[2],cut=id_cut)
            if histogram.Integral() > 0:
                histogram.Scale(1.0/histogram.Integral())
            temp_hists.append(histogram)

        
        temp_hists[0].Draw()
        temp_hists[0].SetLineColor(ROOT.kBlack)
        temp_hists[1].Draw('same')
        temp_hists[1].SetLineColor(ROOT.kBlue)
        temp_hists[2].Draw('same')
        temp_hists[2].SetLineColor(ROOT.kRed)

        leg = ROOT.TLegend(0.75,0.75,0.9,0.9)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)

        for i in range(3):
            leg.AddEntry(temp_hists[i],id_region_names[i],"L")
        leg.Draw('same')

        c1.Print(outpath+tree_name+'_mgg.pdf')
        c1.Print(outpath+tree_name+'_mgg.png')



