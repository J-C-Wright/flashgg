import ROOT
import os
import subprocess
import stat
from os import popen
from math import fabs
import numpy as np

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

def makeTableContents(column_labels = [], values = [],variables=[]):

    outString = 'Variable'
    for label in column_labels:
        outString += ' & '+label
    outString += ' \\\\ \\hline\n'

    for row,variable in zip(values,variables):
        outString += '$'+variable+'$'
        for value in row:
            outString += ' & %5.3f' % (value,)
        outString += ' \\\\ \n'

    return outString

def makeTablePdf(contents=[]):

    latex = '\\documentclass{article}\n'
    latex += '\\usepackage{graphicx}\n'

    latex += '\\begin{document}\n'
    latex += '\\begin{center}\n'
    latex += '\\begin{tabular}{| l || c | c | c | c | c | c |}\n'

    latex += '\\hline\n'

    latex += contents
    latex += '\\hline\n'
    latex += '\\end{tabular}\n'
    latex += '\\end{center}\n'
    latex += '\\end{document}\n'

    latexFile = open('ROC_Table.tex','w')
    latexFile.write(latex)
    latexFile.close()

    command = 'pdflatex ROC_Table.tex'
    output = popen(command).read()
    print output

def separationPower(signal=None,bkg=None):

    nbins = signal.GetXaxis().GetNbins()
    separation = 0
    for i in range(1,nbins):
        separation += 0.5*fabs(signal.GetBinContent(i)-bkg.GetBinContent(i))
    return separation

def ROCCurve(signal=None,bkg=None):

    nbins = signal.GetXaxis().GetNbins()
    x = [1]
    y = [1]
    for i in range(1,nbins):
        TP = signal.Integral(i,nbins-1)
        FP = bkg.Integral(i,nbins-1)
        if TP+FP > 0.0:
            x.append(FP)
            y.append(TP)
    x.append(0)
    y.append(0)

    x.reverse()
    y.reverse()

    x = x
    y = y

    return [x,y]
       
def areaUnderROC(x = [],y = []):

    area = 0
    for i in range(len(x)-1):
        
        d_area = 0.5*fabs((x[i+1]**2 - x[i]**2) - (x[i+1]-x[i])*(y[i+1]+y[i]))
        area += d_area

    return area







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

dph_ps_cut = '((dipho_mass > 100 && dipho_mass < 180)&&(dijet_LeadJPt>0)&&(dijet_SubJPt>0)&&(leadPho_PToM>1/2.0)&&(sublPho_PToM>1/4.0))'
vbf_ps_cut = '((dijet_Mjj>250)&&(dijet_LeadJPt>30)&&(dijet_SubJPt>20)&&(fabs(dijet_leadEta) < 4.7 && fabs(dijet_subleadEta) < 4.7))'

#cut = '&&'.join([dph_ps_cut, vbf_ps_cut])
cut = '&&'.join([dph_ps_cut])



#VBF vs all the others
vbf_vs_all_auc = []
vbf_vs_all_graph = []
for histo_info in histos_info:

    signal_histo,bkg_histo = getSBHistoPair(signal_trees=vbf_trees,bkg_trees=gjet_trees+dbox_trees+qcd_trees+ggh_trees,
                                            histo_info=histo_info,cut=cut)
    x,y = ROCCurve(signal=signal_histo,bkg=bkg_histo)


    x_array = np.asarray(x)
    y_array = np.asarray(y)

    c1 = ROOT.TCanvas('c1',"",500,500)
    c1.SetFixedAspectRatio()

    ROC_graph = ROOT.TGraph(len(x_array),x_array,y_array)
    ROC_graph.Draw()
    c1.Print('plots/ROCs/'+histo_info[0]+'_ROC.pdf')

    area = areaUnderROC(x_array,y_array)
    vbf_vs_all_auc.append(area)
    vbf_vs_all_graph.append(ROC_graph)




#VBF vs the individual background samples
bkg_channel_trees = [ggh_trees,gjet_trees,dbox_trees,qcd_trees]
vbf_vs_each_auc = []
vbf_vs_each_graph = []
bkg_names = ['ggH','GJet','DP-Box','QCD']
for bkg_trees,bkg_name in zip(bkg_channel_trees,bkg_names):

    temp_vbf_vs_auc = []
    temp_vbf_vs_graph = []
    for histo_info in histos_info:

        signal_histo,bkg_histo = getSBHistoPair(signal_trees=vbf_trees,bkg_trees=bkg_trees,
                                                histo_info=histo_info,cut=cut)

        x,y = ROCCurve(signal=signal_histo,bkg=bkg_histo)
        x_array = np.asarray(x)
        y_array = np.asarray(y)

        c1 = ROOT.TCanvas('c1',"",500,500)
        c1.SetFixedAspectRatio()

        ROC_graph = ROOT.TGraph(len(x_array),x_array,y_array)
        ROC_graph.Draw()
        c1.Print('plots/ROCs/'+bkg_name+'_'+histo_info[0]+'_ROC.pdf')
        c1.Print('plots/ROCs/'+bkg_name+'_'+histo_info[0]+'_ROC.png')

        area = areaUnderROC(x_array,y_array)

        temp_vbf_vs_auc.append(area)
        temp_vbf_vs_graph.append(ROC_graph)

    vbf_vs_each_auc.append(temp_vbf_vs_auc)
    vbf_vs_each_graph.append(temp_vbf_vs_graph)

names = [row[0] for row in histos_info]
colours = [ROOT.kRed,ROOT.kBlack,ROOT.kBlue,ROOT.kGreen,ROOT.kMagenta]
for i,row in enumerate(zip(vbf_vs_all_graph,vbf_vs_each_graph[0],vbf_vs_each_graph[1],vbf_vs_each_graph[2],vbf_vs_each_graph[3])):

    c1 = ROOT.TCanvas('c1',"",500,500)
    row[0].Draw()
    for entry,colour in zip(row,colours):
        entry.SetLineColor(colour)
        entry.Draw("same")

    c1.Print('plots/ROCs/Comp_'+names[i]+'_ROC.pdf')
    c1.Print('plots/ROCs/Comp_'+names[i]+'_ROC.png')

all_sort = sorted(range(len(vbf_vs_all_auc)),key=lambda x:vbf_vs_all_auc[x])
all_sort.reverse()
variables = [row[-1] for row in histos_info]

sorted_values = []
sorted_variables = []
for index in all_sort:
    row = [vbf_vs_all_auc[index],vbf_vs_each_auc[0][index],vbf_vs_each_auc[1][index],vbf_vs_each_auc[2][index],vbf_vs_each_auc[3][index]]
    sorted_values.append(row)
    sorted_variables.append(variables[index])
    print index, row


contents = makeTableContents(column_labels = ['All']+bkg_names, 
                            values = sorted_values,
                            variables = sorted_variables)

makeTablePdf(contents=contents)





