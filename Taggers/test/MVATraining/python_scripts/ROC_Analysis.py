import ROOT
import os
import subprocess
import stat
from os import popen
from math import fabs
import numpy as np
from optparse import OptionParser
from VariableStudies import VariableStudies as vs
from shutil import copyfile

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



    outpath = '/afs/cern.ch/user/j/jwright/www/vbf_ROC_'+cut_name+'_'+varsInfoFile.split('.')[0]+'/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    copyfile('/afs/cern.ch/user/j/jwright/www/index.php',outpath+'index.php')

    tableName = 'ROC_'+cut_name+'_'+varsInfoFile.split('.')[0]
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

    bkg_channel_trees = [all_bkg_trees,ggh_trees,gjet_trees,dbox_trees,qcd_trees]
    vbf_vs_graphs = []
    vbf_vs_aucs = []
    bkg_names = ['All','ggH','GJet','DP-Box','QCD']


    c1 = ROOT.TCanvas('c1',"",500,500)
    c1.SetFixedAspectRatio()
    for bkg_trees,bkg_name in zip(bkg_channel_trees,bkg_names):

        print "Processing vbs vs ",bkg_name

        temp_vbf_vs_auc = []
        temp_vbf_vs_graph = []

        for histo_info in histos_info:

            varname = histo_info[0]
            varname = varname.replace('/','_div_')
            varname = varname.replace('*','_prod_')


            print "--> Processing ",varname


            signal_histo,bkg_histo = vs.getSBHistoPair(
                                                    signal_trees=vbf_trees,bkg_trees=bkg_trees,
                                                    histo_info=histo_info,
                                                    cut=cut)

            x,y = vs.ROCCurve(signal=signal_histo,bkg=bkg_histo)
            x_array = np.asarray(x)
            y_array = np.asarray(y)


            ROC_graph = ROOT.TGraph(len(x_array),x_array,y_array)
            #ROC_graph.Draw()
            #c1.Print(outpath+'roc_vbf_vs_'+bkg_name+'_'+varname+'_ROC.pdf')
            #c1.Print(outpath+'roc_vbf_vs_'+bkg_name+'_'+varname+'_ROC.png')
            #c1.Clear()

            area = vs.areaUnderROC(x_array,y_array)

            temp_vbf_vs_auc.append(area)
            temp_vbf_vs_graph.append(ROC_graph)

        vbf_vs_aucs.append(temp_vbf_vs_auc)
        vbf_vs_graphs.append(temp_vbf_vs_graph)







    colours = [ROOT.kRed,ROOT.kBlack,ROOT.kBlue,ROOT.kGreen,ROOT.kMagenta]
    names = [row[0] for row in histos_info]
    
    for i,name in enumerate(names):

        graphs = [row[i] for row in vbf_vs_graphs]
     
        varname = name
        varname = varname.replace('/','_div_')
        varname = varname.replace('*','_prod_')

        leg = ROOT.TLegend(0.7,0.1,0.9,0.48)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)

        graphs[0].Draw()
        for graph,colour,bkg_name in zip(graphs,colours,bkg_names):
            graph.SetLineColor(colour)
            graph.Draw("same")
            leg.AddEntry(graph,bkg_name,"L")

        leg.Draw('same')

        c1.Print(outpath+'comp_roc_vbf_vs_bkg_'+varname+'.pdf')
        c1.Print(outpath+'comp_roc_vbf_vs_bkg_'+varname+'.png')

        c1.Clear()


    variables = [row[4] for row in histos_info]
    highlights =  [row[5] for row in histos_info]
    
    all_sort = sorted(range(len(vbf_vs_aucs[0])),key=lambda x:vbf_vs_aucs[0][x])
    all_sort.reverse()

    sorted_values = []
    sorted_variables = []
    sorted_highlights = []
    for index in all_sort:
        row = []
        for vbf_vs_x in vbf_vs_aucs:
            row.append(vbf_vs_x[index])
        sorted_values.append(row)
        sorted_variables.append(variables[index])
        sorted_highlights.append(highlights[index])
        print index, row

    vs.makeTablePdf(column_labels=bkg_names,values=sorted_values,variables=sorted_variables,path=outpath,name=tableName,highlights=sorted_highlights)

    print "Table and plots can be found at:"
    print "http://jwright.web.cern.ch/jwright/"+outpath.split('www/')[-1]




