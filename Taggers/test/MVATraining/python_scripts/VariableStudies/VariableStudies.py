import ROOT
import os
import subprocess
import stat
from os import popen
from math import fabs

def getListOfBranchesFromTree(tree=None,outputPath='data/',outputFile='Branches.dat'):

    array = tree.GetListOfBranches()
    output = ''
    for i in range(array.GetEntries()):
        output += array.At(i).GetName() + '\n'

    outputFile = open(outputPath+outputFile,'w')
    outputFile.write(output)
    outputFile.close()


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
        if signal_histo.Integral() > 0 and bkg_histo.Integral() > 0:
            signal_histo.Scale(1/signal_histo.Integral())
            bkg_histo.Scale(1/bkg_histo.Integral())

    return [signal_histo,bkg_histo]

def separationPower(signal=None,bkg=None):

    nbins = signal.GetXaxis().GetNbins()
    separation = 0
    for i in range(1,nbins):
        separation += 0.5*fabs(signal.GetBinContent(i)-bkg.GetBinContent(i))
    return separation

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

    latexFile = open('separationTable.tex','w')
    latexFile.write(latex)
    latexFile.close()

    command = 'pdflatex separationTable.tex'
    output = popen(command).read()
    print output

    command = 'convert -density 300  separationTable.pdf -quality 100 separationTable.png'
    output = popen(command).read()
    print output
    


       



