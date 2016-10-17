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

    if signal_region:
        num_slices = 20
    else:
        num_slices = 50

    outpath = '/afs/cern.ch/user/j/jwright/www/vbf_corr_'+cut_name+'_'+varsInfoFile.split('.')[0]+'/'
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    copyfile('/afs/cern.ch/user/j/jwright/www/index.php',outpath+'index.php')

    tableName = 'VBFCorr_'+cut_name+'_'+varsInfoFile.split('.')[0]

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



    c1 = ROOT.TCanvas('c1')

    corr_coeffs = []
    std_devs = []

    for i in range(1,len(histos_info)):

        histo_info=histos_info[i]
        print histo_info[0]
        
        limits = vs.getSliceRanges(trees=vbf_trees,x_info=histo_info,num_slices=num_slices,cut=cut)
        means,mids,errs = vs.getSliceHistos(trees = vbf_trees,x_info=histo_info,y_info=histos_info[0],slices=limits,cut=cut)

        x = np.asarray(mids)
        y = np.asarray(means)


        graph1 = ROOT.TGraph(len(x),x,y)
        graph1.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
        graph1.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))

        graph2 = ROOT.TGraph(len(x),x,y-errs)
        graph2.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
        graph2.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
        graph2.SetLineColor(ROOT.kRed)

        graph3 = ROOT.TGraph(len(x),x,y+errs)
        graph3.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
        graph3.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
        graph3.SetLineColor(ROOT.kRed)

        graph1.Draw()
        graph2.Draw('same')
        graph3.Draw('same')

        graph1.GetYaxis().SetRangeUser(123.0,127.0)

        graph1.Draw()
        graph2.Draw('same')
        graph3.Draw('same')

        c1.Update()

        c1.Print(outpath+'vbf_corr_'+histo_info[0].replace('/','_div_')+'_vs_'+histos_info[0][0]+'.pdf')
        c1.Print(outpath+'vbf_corr_'+histo_info[0].replace('/','_div_')+'_vs_'+histos_info[0][0]+'.png')

        c1.Clear()

        coeff =  np.corrcoef(x=x,y=y)
        std_dev = np.std(y)

        std_devs.append(std_dev)
        corr_coeffs.append(coeff[0][1])






    #Background correlations
    bkg_trees_list = [ggh_trees,gjet_trees,dbox_trees,qcd_trees]
    bkg_names = ['ggH','GJet','DP-Box','QCD']

    #all together
    bkg_trees = ggh_trees#+gjet_trees+dbox_trees+qcd_trees
    corr_coeffs_allbkg = []
    std_devs_allbkg = []
    each_correlation = []
    each_stddev = []

    for bkg_trees,bkg_name in zip(bkg_trees_list,bkg_names):
        print bkg_trees[0].GetName()

        temp_correlation = []
        temp_stddev = []

        for i in range(1,len(histos_info)):

            histo_info=histos_info[i]

            print histo_info[0],

            limits_bkg = vs.getSliceRanges(trees=bkg_trees,x_info=histo_info,num_slices=num_slices,cut=cut)
            means_bkg,mids_bkg,err_bkg = vs.getSliceHistos(trees=bkg_trees,x_info=histo_info,y_info=histos_info[0],slices=limits_bkg,cut=cut)

            limits_vbf = vs.getSliceRanges(trees=vbf_trees,x_info=histo_info,num_slices=num_slices,cut=cut)
            means_vbf,mids_vbf,err_vbf = vs.getSliceHistos(trees = vbf_trees,x_info=histo_info,y_info=histos_info[0],slices=limits_vbf,cut=cut)

            x_bkg = np.asarray(mids_bkg)
            y_bkg = np.asarray(means_bkg)
            y_err_bkg = np.asarray(err_bkg)
            x_vbf = np.asarray(mids_vbf)
            y_vbf = np.asarray(means_vbf)
            y_err_vbf = np.asarray(err_vbf)


            graph_vbf_1 = ROOT.TGraph(len(x_vbf),x_vbf,y_vbf)
            graph_vbf_1.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
            graph_vbf_1.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
            graph_vbf_2 = ROOT.TGraph(len(x_vbf),x_vbf,y_vbf-y_err_vbf)
            graph_vbf_2.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
            graph_vbf_2.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
            graph_vbf_2.SetLineColor(ROOT.kRed)
            graph_vbf_3 = ROOT.TGraph(len(x_vbf),x_vbf,y_vbf+y_err_vbf)
            graph_vbf_3.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
            graph_vbf_3.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
            graph_vbf_3.SetLineColor(ROOT.kRed)

            graph_bkg_1 = ROOT.TGraph(len(x_bkg),x_bkg,y_bkg)
            graph_bkg_1.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
            graph_bkg_1.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
            graph_bkg_1.SetLineColor(ROOT.kBlue)
            graph_bkg_2 = ROOT.TGraph(len(x_bkg),x_bkg,y_bkg-y_err_bkg)
            graph_bkg_2.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
            graph_bkg_2.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
            graph_bkg_2.SetLineColor(ROOT.kGreen)
            graph_bkg_3 = ROOT.TGraph(len(x_bkg),x_bkg,y_bkg+y_err_bkg)
            graph_bkg_3.GetXaxis().SetTitle(histo_info[4].replace('\\','#'))
            graph_bkg_3.GetYaxis().SetTitle(histos_info[0][4].replace('\\','#'))
            graph_bkg_3.SetLineColor(ROOT.kGreen)


            graph_vbf_1.Draw()
            graph_vbf_2.Draw('same')
            graph_vbf_3.Draw('same')
            graph_bkg_1.Draw('same')
            graph_bkg_2.Draw('same')
            graph_bkg_3.Draw('same')
            graph_vbf_1.GetYaxis().SetRangeUser(123.0,127.0)
            graph_vbf_1.Draw()
            graph_vbf_2.Draw('same')
            graph_vbf_3.Draw('same')
            graph_bkg_1.Draw('same')
            graph_bkg_2.Draw('same')
            graph_bkg_3.Draw('same')

            c1.Update()
            c1.Print(outpath+'vbf_vs_'+bkg_name+'_'+histo_info[0].replace('/','_div_')+'_vs_'+histos_info[0][0]+'.pdf')
            c1.Print(outpath+'vbf_vs_'+bkg_name+'_'+histo_info[0].replace('/','_div_')+'_vs_'+histos_info[0][0]+'.png')

            c1.Clear()

            coeff =  np.corrcoef(x=x_bkg,y=y_bkg)
            std_dev = np.std(y_bkg)

            print coeff[0][1],std_dev

            temp_stddev.append(std_dev)
            temp_correlation.append(coeff[0][1])

        each_correlation.append(temp_correlation)
        each_stddev.append(temp_stddev)









    variables = [row[4] for row in histos_info[1:]]
    highlights = [row[5] for row in histos_info[1:]]

    all_sort = sorted(range(len(corr_coeffs)),key=lambda x:corr_coeffs[x])
    all_sort.reverse()

    print all_sort

    sorted_variables = []
    sorted_coeffs = []
    sorted_std_devs = []
    sorted_highlights = []

    sorted_values = []

    column_labels = ['vbf r','vbf $\\sigma$']
    for name in bkg_names:
        column_labels.append(name+' r')
        column_labels.append(name+' $\\sigma$')

    for index in all_sort:

        row = []
        row.append(corr_coeffs[index])
        row.append(std_devs[index])
        for name,correlation,stddev in zip(bkg_names,each_correlation,each_stddev):
            row.append(correlation[index])
            row.append(stddev[index])

        sorted_values.append(row)
        sorted_variables.append(variables[index])
        sorted_highlights.append(highlights[index])

    vs.makeTablePdf(column_labels=column_labels,values=sorted_values,variables=sorted_variables,path=outpath,name=tableName,orientation='landscape',highlights=sorted_highlights)

    print "Table and plots can be found at:"
    print "http://jwright.web.cern.ch/jwright/"+outpath.split('www/')[-1]





