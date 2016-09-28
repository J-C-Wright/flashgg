import ROOT
import os
import sys
import subprocess
import stat
from os import popen
from math import fabs,sqrt
import numpy as np

def getInfoFromFile(path = 'data/',name=''):

    if name == '': raise Exception,'Need filename'

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

    if name == '': raise Exception,'Need filename'

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

def makeTablePdf(contents=[],filename=''):

    latex = '\\documentclass{article}\n'
    latex += '\\usepackage{graphicx}\n'

    latex += '\\begin{document}\n'
    latex += '\\begin{center}\n'
    latex += '\\begin{tabular}{| l || c | c |}\n'

    latex += '\\hline\n'

    latex += contents
    latex += '\\hline\n'
    latex += '\\end{tabular}\n'
    latex += '\\end{center}\n'
    latex += '\\end{document}\n'

    latexFile = open(filename+'.tex','w')
    latexFile.write(latex)
    latexFile.close()

    command = 'pdflatex '+filename
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
        area += 0.5*(x[i+1]-x[i])*(y[i+1]+y[i])

    return area



def getSliceRanges(trees=None,x_info=None,cut=None,num_slices=None):

    bins = int(x_info[1])*2
    var_hist = ROOT.TH1F('temp_mgg','',bins,float(x_info[2]),float(x_info[3]))
    for tree in trees:
        for tree in trees:
            tree.Draw(x_info[0]+'>>var('+str(bins)+','+x_info[2]+','+x_info[3]+')',cut)
            hist = ROOT.gROOT.FindObject('var')
            hist.SetDirectory(0)
            var_hist.Add(hist)
        var_hist.SetDirectory(0)

    total =  var_hist.Integral()
    slice_area = total/float(num_slices)
    low = 1
    hi = 2
    limits =[var_hist.GetBinLowEdge(1)]
    for i in range(num_slices):
        temp_area = 0

        count = 0
        while temp_area < slice_area:

            temp_area = var_hist.Integral(low,hi)
            hi += 1
            count += 1

            if hi > bins: 
                break
        limits.append(var_hist.GetBinLowEdge(hi))
        low = hi+1
        hi += 2
    limits.append(var_hist.GetBinLowEdge(bins+1))

    return limits







def getSliceHistos(trees=None,x_info=None,y_info=None,cut=None,slices=None):

    sliceMeans = []
    sliceMids = []
    sliceErrs = []
    
    for i in range(len(slices)-1):

        temp_slice_hist = ROOT.TH1F('temp_mgg','',int(y_info[1]),float(y_info[2]),float(y_info[3]))

        slice_cut = '('+x_info[0]+'>'+str(slices[i])+' && '+x_info[0]+'<'+str(slices[i+1])+')&&'+cut
        
        sliceMid = (slices[i] + slices[i+1])*0.5

        for tree in trees:
            tree.Draw(y_info[0]+'>>var('+y_info[1]+','+y_info[2]+','+y_info[3]+')',slice_cut)
            hist = ROOT.gROOT.FindObject('var')
            hist.SetDirectory(0)
            temp_slice_hist.Add(hist)
        temp_slice_hist.SetDirectory(0)

        if temp_slice_hist.GetMean() > 0:
            sliceMeans.append(temp_slice_hist.GetMean())
            sliceMids.append(sliceMid)
            sliceErrs.append(temp_slice_hist.GetMeanError())

    return [sliceMeans,sliceMids,sliceErrs]





ROOT.gROOT.SetBatch(ROOT.kTRUE)

vbf_info, ggh_info, gjet_info, dbox_info, qcd_info, _ = getInfoFromFile(name='YacineTrees.dat')

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

vbf_trees = getTrees(tfiles=vbf_tfiles,treenames=vbf_info[1])
ggh_trees = getTrees(tfiles=ggh_tfiles,treenames=ggh_info[1])
gjet_trees = getTrees(tfiles=gjet_tfiles,treenames=gjet_info[1])
dbox_trees = getTrees(tfiles=dbox_tfiles,treenames=dbox_info[1])
qcd_trees = getTrees(tfiles=qcd_tfiles,treenames=qcd_info[1])

histos_info = getHistosInfoFromFile(name='histograms_info_corr.dat')

dph_ps_cut = '((dipho_mass > 100 && dipho_mass < 180)&&(dijet_LeadJPt>0)&&(dijet_SubJPt>0)&&(leadPho_PToM>1/2.0)&&(sublPho_PToM>1/4.0))'
signal_region = '(dipho_mass > 115 && dipho_mass < 135)'
vbf_ps_cut = '((dijet_Mjj>250)&&(dijet_LeadJPt>30)&&(dijet_SubJPt>20)&&(fabs(dijet_leadEta) < 4.7 && fabs(dijet_subleadEta) < 4.7))'

#cut = '&&'.join([dph_ps_cut,signal_region])
cut = '&&'.join([dph_ps_cut, vbf_ps_cut, signal_region])



corr_coeffs = []
std_devs = []


for i in range(1,len(histos_info)):

    histo_info=histos_info[i]

    print histo_info[0]
    
    limits = getSliceRanges(trees=vbf_trees,x_info=histo_info,num_slices=100,cut=cut)
    means,mids,errs = getSliceHistos(trees = vbf_trees,x_info=histo_info,y_info=histos_info[0],slices=limits,cut=cut)

    x = np.asarray(mids)
    y = np.asarray(means)


    c1 = ROOT.TCanvas('c1')
    graph1 = ROOT.TGraph(len(x),x,y)
    graph1.GetXaxis().SetTitle(histo_info[-1].replace('\\','#'))
    graph1.GetYaxis().SetTitle(histos_info[0][-1].replace('\\','#'))

    graph2 = ROOT.TGraph(len(x),x,y-errs)
    graph2.GetXaxis().SetTitle(histo_info[-1].replace('\\','#'))
    graph2.GetYaxis().SetTitle(histos_info[0][-1].replace('\\','#'))
    graph2.SetLineColor(ROOT.kRed)

    graph3 = ROOT.TGraph(len(x),x,y+errs)
    graph3.GetXaxis().SetTitle(histo_info[-1].replace('\\','#'))
    graph3.GetYaxis().SetTitle(histos_info[0][-1].replace('\\','#'))
    graph3.SetLineColor(ROOT.kRed)

    graph1.Draw()
    graph2.Draw('same')
    graph3.Draw('same')

    graph1.GetYaxis().SetRangeUser(124.0,125.0)

    graph1.Draw()
    graph2.Draw('same')
    graph3.Draw('same')

    c1.Update()

    c1.Print('plots/correlation/'+histo_info[0]+'_vs_'+histos_info[0][0]+'.pdf')
    c1.Print('plots/correlation/'+histo_info[0]+'_vs_'+histos_info[0][0]+'.png')

    coeff =  np.corrcoef(x=x,y=y)
    std_dev = np.std(y)

    std_devs.append(std_dev)
    corr_coeffs.append(coeff[0][1])

variables = [row[-1] for row in histos_info[1:]]

all_sort = sorted(range(len(corr_coeffs)),key=lambda x:corr_coeffs[x])
all_sort.reverse()

sorted_variables = []
sorted_coeffs = []
sorted_std_devs = []

sorted_values = []

for index in all_sort:
    row = [corr_coeffs[index],std_devs[index]]
    sorted_values.append(row)
    sorted_variables.append(variables[index])

contents = makeTableContents(column_labels = ['r','$\\sigma$'],values=sorted_values,variables=sorted_variables)
makeTablePdf(contents=contents,filename='correlationTable')







