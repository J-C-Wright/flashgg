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

    return parser.parse_args()




if __name__ == '__main__':

    (opt,args) = get_options()

    cut_name = opt.cut_name
    varsInfoFile = opt.varsFile
    treesFile = opt.treesFile
    signal_region = opt.signal_region


    outpath = '/afs/cern.ch/user/j/jwright/www/correlations_'+varsInfoFile.split('.')[0]+'/'
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

    all_bkg_trees = ggh_trees+gjet_trees+dbox_trees+qcd_trees

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

#histos_info = histos_info[:5]

num_of_vars = len(histos_info)
spearmans = []
pearsons = []
for i in range(num_of_vars):
    spearmans.append([0]*num_of_vars)
    pearsons.append([0]*num_of_vars)

sample_size = 100000
bins_factor = 2

c1 = ROOT.TCanvas('c1',"",500,500)
c1.SetFixedAspectRatio()

sample_names = ['ggH']
tree_lists = [ggh_trees]

for sample_name,tree_list in zip(sample_names,tree_lists):
    for i,var1 in enumerate(histos_info):
        spearmans[i][i] = 1
        pearsons[i][i] = 1
        for j,var2 in enumerate(histos_info):

            if i >= j: continue
            print var1[0]+'_vs_'+var2[0],

            var1_sample = vs.getSampleOfBranch(trees=tree_list,var_info=var1,cut=cut,sample_size=sample_size)
            var2_sample = vs.getSampleOfBranch(trees=tree_list,var_info=var2,cut=cut,sample_size=sample_size)

            spearman_result = scipy.stats.spearmanr(a=var1_sample,b=var2_sample)
            pearson_result = scipy.stats.pearsonr(x=var1_sample,y=var2_sample)
            spearmans[i][j] = spearman_result[0]
            spearmans[j][i] = spearman_result[0]
            pearsons[i][j] = pearson_result[0]
            pearsons[j][i] = pearson_result[0]

            print spearman_result[0],pearson_result[0]


            hist = ROOT.TH2F(var1[0]+'_vs_'+var2[0],'',int(var1[1])*bins_factor,float(var1[2]),float(var1[3]),
                                                       int(var2[1])*bins_factor,float(var2[2]),float(var2[3]))
            hist.SetDirectory(0)
            for v1,v2 in zip(var1_sample,var2_sample):
                hist.Fill(v1,v2)

            profile = hist.ProfileX()
            profile.SetLineColor(ROOT.kRed)
            profile.SetDirectory(0)

            hist.Draw('colz')
            profile.Draw('same')

            c1.Print(outpath+var1[0]+'_vs_'+var2[0]+'_'+sample_name+'.pdf')
            c1.Print(outpath+var1[0]+'_vs_'+var2[0]+'_'+sample_name+'.png')

            profile.Draw()
            c1.Print(outpath+var1[0]+'_vs_'+var2[0]+'_'+sample_name+'_profile.pdf')
            c1.Print(outpath+var1[0]+'_vs_'+var2[0]+'_'+sample_name+'_profile.png')

            c1.Clear()

    np_spearmans = np.array(spearmans)
    np_spearmans = (np_spearmans+1)/2.0
    np_pearsons = np.array(pearsons)
    np_pearsons = (np_pearsons+1)/2.0

    #Spearman's r heatmap
    fig,ax = plt.subplots()
    heatmap = ax.pcolor(np_spearmans,cmap=plt.cm.seismic,alpha=1.0)
    fit = plt.gcf()
    fig.set_size_inches(10,10)
    ax.set_frame_on(False)

    labels = ['$'+row[4]+'$' for row in histos_info]

    ax.set_xticks(np.arange(len(labels))+0.5,minor=False)
    ax.set_yticks(np.arange(len(labels))+0.5,minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(labels,minor=False)
    ax.set_yticklabels(labels,minor=False)
    plt.xticks(rotation=90)

    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick10n = False
        t.tick20n = False
    for t in ax.yaxis.get_major_ticks():
        t.tick10n = False
        t.tick20n = False

    plt.savefig(outpath+sample_name+'_spearmans_heatmap.pdf',format='pdf',pad_inches=1.0)
    plt.savefig(outpath+sample_name+'_spearmans_heatmap.png',format='png',pad_inches=1.0)

    #Pearson's r heatmap
    heatmap = ax.pcolor(np_pearsons,cmap=plt.cm.seismic,alpha=1.0)
    plt.savefig(outpath+sample_name+'_pearsons_heatmap.pdf',format='pdf',pad_inches=1.0)
    plt.savefig(outpath+sample_name+'_pearsons_heatmap.png',format='png',pad_inches=1.0)



