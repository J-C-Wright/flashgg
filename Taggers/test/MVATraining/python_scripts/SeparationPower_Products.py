import ROOT
import os
import subprocess
import stat
from os import popen
from math import fabs
from shutil import copyfile
from optparse import OptionParser
from VariableStudies import VariableStudies as vs


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


    outpath = '/afs/cern.ch/user/j/jwright/www/vbf_sep_products_'+cut_name+'_'+varsInfoFile.split('.')[0]+'/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    copyfile('/afs/cern.ch/user/j/jwright/www/index.php',outpath+'index.php')

    tableName = 'SepPower_products_'+cut_name+'_'+varsInfoFile.split('.')[0]

    ROOT.gROOT.SetBatch(ROOT.kTRUE)


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


    print "Cut is: ",cut



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


    c1 = ROOT.TCanvas('c1')

    #VBF vs 
    bkg_channel_trees = [all_bkg_trees,ggh_trees,gjet_trees,dbox_trees,qcd_trees]
    vbf_vs = []
    bkg_names = ['All','ggH','GJet','DP-Box','QCD']

    opString = '%s*%s'

    variables = []

    for i in range(len(histos_info)):
        for j in range(len(histos_info)):

            if i > j: continue

            var1_info = histos_info[i]
            var2_info = histos_info[j]
            #if var1_info[0] == var2_info[0]: continue
            newVar = opString % (var1_info[0],var2_info[0])
            print newVar
            variables.append(opString % (var1_info[-1],var2_info[-1]))


    for bkg_trees,bkg_name in zip(bkg_channel_trees,bkg_names):

        print "Processing vbf vs ",bkg_name

        temp_vbf_vs = []

        for i in range(len(histos_info)):
            for j in range(len(histos_info)):

                var1_info = histos_info[i]
                var2_info = histos_info[j]

                if i > j: continue
                newVar = opString % (var1_info[0],var2_info[0])
                
                print "--> Processing ",newVar

                signal_histo,bkg_histo = vs.getCombinationVarSBHistoPair(signal_trees=vbf_trees,bkg_trees=bkg_trees,
                                                                         var1_info=var1_info,var2_info=var2_info,
                                                                         opString=opString,cut=cut)

                sep_power =  vs.separationPower(signal_histo,bkg_histo)
                temp_vbf_vs.append(sep_power)

                signal_histo.SetLineColor(ROOT.kRed)
                bkg_histo.SetLineColor(ROOT.kBlue)

                label_newVar = opString % (var1_info[-1],var2_info[-1])

                signal_histo.SetXTitle(label_newVar.replace('\\','#'))
                bkg_histo.SetXTitle(label_newVar.replace('\\','#'))

                signal_max = signal_histo.GetBinContent(signal_histo.GetMaximumBin())
                bkg_max = bkg_histo.GetBinContent(bkg_histo.GetMaximumBin())
                if signal_max > bkg_max:
                    signal_histo.Draw()
                    bkg_histo.Draw("same")
                else:
                    bkg_histo.Draw()
                    signal_histo.Draw("same")

                c1.Print(outpath+'vbf_vs_'+bkg_name+'_'+var1_info[0]+'_'+var2_info[0]+'.pdf')
                c1.Print(outpath+'vbf_vs_'+bkg_name+'_'+var1_info[0]+'_'+var2_info[0]+'.png')
                c1.Clear()

                print 'Separation power is ',sep_power

        vbf_vs.append(temp_vbf_vs)







    all_sort = sorted(range(len(vbf_vs[0])),key=lambda x:vbf_vs[0][x])
    all_sort.reverse()

    sorted_values = []
    sorted_variables = []
    for index in all_sort:
        row = []
        for vbf_vs_x in vbf_vs:
            row.append(vbf_vs_x[index])
        sorted_values.append(row)
        sorted_variables.append(variables[index])

    vs.makeTablePdf(column_labels=bkg_names,values=sorted_values,variables=sorted_variables,path=outpath,name=tableName)

    print "Table and plots can be found at:"
    print "http://jwright.web.cern.ch/jwright/"+outpath.split('www/')[-1]







