import ROOT
import os
import subprocess
import stat
from os import popen
from math import fabs

def getInfoFromFile(path = 'data/',name=''):

    if name == '': raise Exception,"Need filename"

    tree_path = ''
    vbf_files = []
    ggh_files = []
    gjet_files = []
    dbox_files = []
    qcd_files = []
    data_files = []
    vbf_treenames = []
    ggh_treenames = []
    gjet_treenames = []
    dbox_treenames = []
    qcd_treenames = []
    data_treenames = []

    with open(path+name) as f:

        lines = f.read().split('\n')
        for line in lines:
            if line == '' or line[0] == '#':continue

            tag,string,_ = line.split()
            if tag == 'folder':
                tree_path = string

        for line in lines:
            if line == '' or line[0] == '#':continue

            tag,string,treename = line.split()

            if  tag == 'data':
                data_files.append(tree_path+string)
                data_treenames.append(treename)
            elif tag == 'gjet_mc':
                gjet_files.append(tree_path+string)
                gjet_treenames.append(treename)
            elif tag == 'dbox_mc':
                dbox_files.append(tree_path+string)
                dbox_treenames.append(treename)
            elif tag == 'qcd_mc':
                qcd_files.append(tree_path+string)
                qcd_treenames.append(treename)
            elif tag == 'ggh_mc':
                ggh_files.append(tree_path+string)
                ggh_treenames.append(treename)
            elif tag == 'vbf_mc':
                vbf_files.append(tree_path+string)
                vbf_treenames.append(treename)

    f.close()    

    return [[vbf_files,vbf_treenames],
            [ggh_files,ggh_treenames],
            [gjet_files,gjet_treenames],
            [dbox_files,dbox_treenames],
            [qcd_files,qcd_treenames],
            [data_files,data_treenames]]

def getTrees(tfiles=[],treenames=[],path_in_tfile = 'vbfTagDumper/trees/'):
    
    trees = []
    for tfile,treename in zip(tfiles,treenames):
        tree = tfile.Get(path_in_tfile+treename)
        trees.append(tree)

    return trees

def getHistosInfoFromFile(path='data/',name=''):

    if name == '': raise Exception,"Need filename"

    histos_info = []
    with open(path+name) as f:
        for line in f.read().split('\n'):
            if line == '' or line[0] == '#':continue
            histos_info.append(line.split())

    f.close()
    return histos_info

def getSBHistoPair(signal_trees=[],bkg_trees=[],histo_info=[],cut='',normalize=True,useWeights=True):

    if useWeights:
        cut = 'weight*'+cut
    
    signal_histo = ROOT.TH1F(histo_info[0]+'_signal','',int(histo_info[1]),float(histo_info[2]),float(histo_info[3]))
    for tree in signal_trees:
        tree.Draw(histo_info[0]+'>>var('+histo_info[1]+','+histo_info[2]+','+histo_info[3]+')',cut)
        hist = ROOT.gROOT.FindObject('var')
        hist.SetDirectory(0)
        signal_histo.Add(hist)
    signal_histo.SetDirectory(0)

    bkg_histo = ROOT.TH1F(histo_info[0]+'_bkg','',int(histo_info[1]),float(histo_info[2]),float(histo_info[3]))
    for tree in bkg_trees:
        tree.Draw(histo_info[0]+'>>var('+histo_info[1]+','+histo_info[2]+','+histo_info[3]+')',cut)
        hist = ROOT.gROOT.FindObject('var')
        hist.SetDirectory(0)
        bkg_histo.Add(hist)
    bkg_histo.SetDirectory(0)

    if normalize:
        signal_histo.Scale(1/signal_histo.Integral())
        bkg_histo.Scale(1/bkg_histo.Integral())

    return [signal_histo,bkg_histo]

def separationPower(signal=None,bkg=None):

    nbins = signal.GetXaxis().GetNbins()
    separation = 0
    for i in range(1,nbins):
        separation += 0.5*fabs(signal.GetBinContent(i)-bkg.GetBinContent(i))
    return separation


        







ROOT.gROOT.SetBatch(ROOT.kTRUE)

vbf_info, ggh_info, gjet_info, dbox_info, qcd_info, _ = getInfoFromFile(name='YacineTrees.dat')

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

vbf_trees = getTrees(tfiles=vbf_tfiles,treenames=vbf_info[1])
ggh_trees = getTrees(tfiles=ggh_tfiles,treenames=ggh_info[1])
gjet_trees = getTrees(tfiles=gjet_tfiles,treenames=gjet_info[1])
dbox_trees = getTrees(tfiles=dbox_tfiles,treenames=dbox_info[1])
qcd_trees = getTrees(tfiles=qcd_tfiles,treenames=qcd_info[1])

histos_info = getHistosInfoFromFile(name='histograms_info.dat')

ps_cut = '((dijet_Mjj>250)&&(dijet_LeadJPt>30)&&(dijet_SubJPt>20)&&(leadPho_PToM>1/3.0)&&(sublPho_PToM>1/4.0))'




#VBF vs all the others
vbf_vs_all = []
for histo_info in histos_info:

    signal_histo,bkg_histo = getSBHistoPair(signal_trees=vbf_trees,bkg_trees=gjet_trees+dbox_trees+qcd_trees+ggh_trees,
                                            histo_info=histo_info,cut=ps_cut)

    sep_power =  separationPower(signal_histo,bkg_histo)
    vbf_vs_all.append(sep_power)

    c1 = ROOT.TCanvas('c1')
    signal_histo.SetLineColor(ROOT.kRed)
    bkg_histo.SetLineColor(ROOT.kBlue)
    signal_histo.Draw()
    bkg_histo.Draw("same")
    c1.Print('plots/vbf_vs_all_'+histo_info[0]+'.pdf')

for row in vbf_vs_all:
    print row

#VBF vs the other non-higgs backgrounds
vbf_vs_nonhiggs = []
for histo_info in histos_info:

    signal_histo,bkg_histo = getSBHistoPair(signal_trees=vbf_trees,bkg_trees=gjet_trees+dbox_trees+qcd_trees,
                                            histo_info=histo_info,cut=ps_cut)

    sep_power =  separationPower(signal_histo,bkg_histo)
    vbf_vs_all.append(sep_power)

    c1 = ROOT.TCanvas('c1')
    signal_histo.SetLineColor(ROOT.kRed)
    bkg_histo.SetLineColor(ROOT.kBlue)
    signal_histo.Draw()
    bkg_histo.Draw("same")
    c1.Print('plots/vbf_vs_nonhiggs_'+histo_info[0]+'.pdf')

for row in vbf_vs_nonhiggs:
    print row


#VBF vs the individual background samples
bkg_channel_trees = [ggh_trees,gjet_trees,dbox_trees,qcd_trees]
vbf_vs_each = []
bkg_names = ['ggH','GJet','Box','QCD']
for bkg_trees,bkg_name in zip(bkg_channel_trees,bkg_names):

    temp_vbf_vs = []
    for histo_info in histos_info:

        signal_histo,bkg_histo = getSBHistoPair(signal_trees=vbf_trees,bkg_trees=bkg_trees,
                                                histo_info=histo_info,cut=ps_cut)

        sep_power =  separationPower(signal_histo,bkg_histo)
        temp_vbf_vs.append(sep_power)

        c1 = ROOT.TCanvas('c1')
        signal_histo.SetLineColor(ROOT.kRed)
        bkg_histo.SetLineColor(ROOT.kBlue)
        signal_histo.Draw()
        bkg_histo.Draw("same")
        c1.Print('plots/vbf_vs_'+bkg_name+'_'+histo_info[0]+'.pdf')

    vbf_vs_each.append(temp_vbf_vs)
    for row in temp_vbf_vs:
        print row

columns = ['All','Non H']+bkg_names
numbers = zip(vbf_vs_all,vbf_vs_nonhiggs,vbf_vs_each[0],vbf_vs_each[1],vbf_vs_each[2],vbf_vs_each[3])




print columns
for row in numbers:
    print row







