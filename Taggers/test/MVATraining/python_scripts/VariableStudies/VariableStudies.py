import ROOT
import os
import subprocess
import stat
from os import popen
from math import fabs
from shutil import copyfile

dph_ps_cut = '((dipho_mass > 100 && dipho_mass < 180)&&(dijet_LeadJPt>0)&&(dijet_SubJPt>0)&&(leadPho_PToM>1/2.0)&&(sublPho_PToM>1/4.0))'
vbf_ps_cut = '((dijet_Mjj>250)&&(dijet_LeadJPt>30)&&(dijet_SubJPt>20)&&(fabs(dijet_leadEta) < 4.7 && fabs(dijet_subleadEta) < 4.7))'
jet3_cut = '(jet3_pt > 0)'
signal_region = '(dipho_mass > 115 && dipho_mass < 135)'


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

def getTrees(tfiles=None,treenames=None,path_in_tfile = 'vbfTagDumper/trees/'):
    
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

def getSBHistoPair(signal_trees=None,bkg_trees=None,histo_info=None,cut='',normalize=True,useWeights=True):

    if useWeights:
        cut = 'weight*'+cut
    
    signal_histo = ROOT.TH1F(histo_info[0].replace('/','div')+'_signal','',int(histo_info[1]),float(histo_info[2]),float(histo_info[3]))
    for tree in signal_trees:
        tree.Draw(histo_info[0]+'>>var('+histo_info[1]+','+histo_info[2]+','+histo_info[3]+')',cut)
        hist = ROOT.gROOT.FindObject('var')
        hist.SetDirectory(0)
        signal_histo.Add(hist)
    signal_histo.SetDirectory(0)

    bkg_histo = ROOT.TH1F(histo_info[0].replace('/','div')+'_bkg','',int(histo_info[1]),float(histo_info[2]),float(histo_info[3]))
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

def getCombinationVarSBHistoPair(signal_trees=None,bkg_trees=None,vars_info=None,opString='',cut='',normalize=True,useWeights=True,whiten=False):

    #Get index of each input
    unique_var_names = []
    for var in vars_info:
        if var[0] not in unique_var_names: unique_var_names.append(var[0])

    indices = []
    for var in vars_info:
        for index,name in enumerate(unique_var_names):
            if var[0] == name: indices.append(index)

    #Build expression
    xyz = ['x','y','z']
    form_inputs = []
    for index in indices:
        form_inputs.append(xyz[index])

    if not whiten:
        expression_xy = opString % tuple(form_inputs)
        expression = opString % tuple([row[0] for row in vars_info])
    else:
        whitened_vars_xy = []
        whitened_vars = []
        for var,var_info in zip(form_inputs,vars_info):
            whitened_vars_xy.append('('+var+'-'+var_info[2]+')/('+var_info[3]+'-'+var_info[2]+')')
            whitened_vars.append('('+var_info[0]+'-'+var_info[2]+')/('+var_info[3]+'-'+var_info[2]+')')
        expression_xy = opString % tuple(whitened_vars_xy)
        expression = opString % tuple(whitened_vars)

    n_inputs = len(unique_var_names)
    if (n_inputs) == 1:
        formula = ROOT.TF1('formula',expression_xy)
    elif (n_inputs) == 2:
        formula = ROOT.TF2('formula',expression_xy)
    elif (n_inputs) == 3:
        formula = ROOT.TF2('formula',expression_xy)
    else:
        raise Exception, 'Too many input variables (>3)'

    #Get bins
    bins = int(vars_info[0][1])
    for var in vars_info:
        if int(var[1]) < bins: bins = int(var[1]) < bins

    #Build limits
    limits = []
    unique_indices = []
    for index in indices:
        if index not in unique_indices: unique_indices.append(index)

    for index in unique_indices:
        limits.append([float(vars_info[index][2]),float(vars_info[index][3])])

    lim_inds = [0]*len(unique_indices)
    incr = 0
    varmin = 9999
    varmax = -9999

    permutation = 0
    incr = 0
    while sum(lim_inds) < len(lim_inds):

        testvalues = []
        for limit,index in zip(limits,lim_inds):
            testvalues.append(limit[index])

        testExpr = 'formula('+','.join(['%s']*n_inputs)+')'
        testExpr.replace(',)',')')
        test = eval(testExpr % tuple(testvalues))

        if test < varmin: varmin = test
        if test > varmax: varmax = test

        lim_inds = lim_inds[1:]+lim_inds[:1]
        permutation += 1

        if permutation == len(lim_inds):
            lim_inds[incr] = 1
            incr += 1
            permutation = 0

    print '\t\t',expression_xy,': ',bins,varmin,varmax

    if useWeights:
        cut = 'weight*'+cut

    varmin_str = str(varmin)
    varmax_str = str(varmax)
    
    signal_histo = ROOT.TH1F(expression+'_signal','',bins,float(varmin_str),float(varmax_str))
    for tree in signal_trees:
        tree.Draw(expression+'>>var('+str(bins)+','+varmin_str+','+varmax_str+')',cut)
        hist = ROOT.gROOT.FindObject('var')
        hist.SetDirectory(0)
        signal_histo.Add(hist)
    signal_histo.SetDirectory(0)

    bkg_histo = ROOT.TH1F(expression+'_bkg','',bins,float(varmin_str),float(varmax_str))
    for tree in bkg_trees:
        tree.Draw(expression+'>>var('+str(bins)+','+varmin_str+','+varmax_str+')',cut)
        hist = ROOT.gROOT.FindObject('var')
        hist.SetDirectory(0)
        bkg_histo.Add(hist)
    bkg_histo.SetDirectory(0)

    if normalize:
        if signal_histo.Integral() > 0 and bkg_histo.Integral() > 0:
            signal_histo.Scale(1/signal_histo.Integral())
            bkg_histo.Scale(1/bkg_histo.Integral())
        else:
            print 'WARNING: empty histo ', signal_histo.Integral(), bkg_histo.Integral()

    return [signal_histo,bkg_histo]

   



def separationPower(signal=None,bkg=None):

    nbins = signal.GetXaxis().GetNbins()
    separation = 0
    for i in range(1,nbins):
        separation += 0.5*fabs(signal.GetBinContent(i)-bkg.GetBinContent(i))
    return separation

def makeTablePdf(column_labels=None,values=None,variables=None,path=None,name=None,orientation='portrait',highlights=None):

    latex = '\\documentclass[a4paper,'+orientation+']{article}\n'
    latex += '\\usepackage{graphicx}\n'
    latex += '\\usepackage{color,colortbl}\n'
    latex += '\\definecolor{Salmon}{rgb}{1,0.63,0.48}\n'

    latex += '\\begin{document}\n'
    latex += '\\begin{center}\n'
    latex += '\\begin{tabular}{| l ||'
    latex += ' c |'*len(values[0])
    latex += '}\n'

    latex += '\\hline\n'
    latex += "Variable"
    for label in column_labels:
        latex += ' & '+label
    latex += ' \\\\ \\hline\n'

    for row,variable,highlight in zip(values,variables,highlights):

        if highlight == 'True': latex += '\\rowcolor{Salmon}\n'
        latex += '$'+variable+'$'
        for value in row:
            latex += ' & %5.3f' % (value)
        latex += ' \\\\ \n'


    latex += '\\hline\n'
    latex += '\\end{tabular}\n'
    latex += '\\end{center}\n'
    latex += '\\end{document}\n'

    latexFile = open(name+'.tex','w')
    latexFile.write(latex)
    latexFile.close()

    command = 'pdflatex '+name+'.tex'
    output = popen(command).read()
    print output

    command = 'convert -density 300  '+name+'.pdf -quality 100 '+name+'.png'
    output = popen(command).read()
    print output
    
    copyfile(name+'.pdf',path+name+'.pdf')
    copyfile(name+'.png',path+name+'.png')

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
       
def areaUnderROC(x = None,y = None):

    area = 0
    for i in range(len(x)-1):
        
        d_area = 0.5*fabs((x[i+1]**2 - x[i]**2) - (x[i+1]-x[i])*(y[i+1]+y[i]))
        area += d_area

    return area


def getSliceRanges(trees=None,x_info=None,cut=None,num_slices=None):

    bins = int(x_info[1])*2
    var_hist = ROOT.TH1F('temp_mgg','',bins,float(x_info[2]),float(x_info[3]))
    for tree in trees:
        for tree in trees:
            tree.Draw(x_info[0]+'>>var('+str(bins)+','+x_info[2]+','+x_info[3]+')','weight*('+cut+')')
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
        if var_hist.GetBinLowEdge(hi)>var_hist.GetBinLowEdge(bins+1):
            break
        limits.append(var_hist.GetBinLowEdge(hi))
        low = hi+1
        hi += 2
    #limits.append(var_hist.GetBinLowEdge(bins+1))

    return limits

def getSliceHistos(trees=None,x_info=None,y_info=None,cut=None,slices=None):

    sliceMeans = []
    sliceMids = []
    sliceErrs = []
    
    for i in range(len(slices)-1):

        temp_slice_hist = ROOT.TH1F('temp_mgg','',int(y_info[1]),float(y_info[2]),float(y_info[3]))

        slice_cut = 'weight*('+'('+x_info[0]+'>'+str(slices[i])+' && '+x_info[0]+'<'+str(slices[i+1])+')&&'+cut+')'
        
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








