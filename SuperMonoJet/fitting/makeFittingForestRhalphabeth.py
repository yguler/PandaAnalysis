#!/usr/bin/env python
from re import sub
from math import *
from array import array
from sys import argv,exit
from os import path,getenv
import os
import ROOT as r
from glob import glob
import argparse
parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--ddt',metavar='ddt',type=bool,default=False)
parser.add_argument('--input',metavar='input',type=str,default=getenv('PANDA_FLATDIR'))

args = parser.parse_args()
nddt = args.ddt
region = args.region
#region = out_region.split('_')[0]
if region=='test':
    is_test = True 
    region = 'signal'
else:
    is_test = False

argv=[]
import PandaAnalysis.Flat.fitting_forest as forest 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.MonoX.MonoXSelection as sel

basedir = args.input
lumi = 35900

def f(x):
    return basedir + x + '.root'

def addN2DDT(rootfile):
        #load the map
	trans = r.TFile(basedir+'/DDT.root')
	f56 = trans.Get("DDT_5by6")
	f53 = trans.Get("DDT_5by3")
        #open the file
        output = rootfile+"_tmp"
	fi = r.TFile(rootfile)
        fo = r.TFile(output,"RECREATE")
        for tkey in fi.GetListOfKeys():
            key=tkey.GetName()
            t = fi.Get("%s"%key)
            clonet = t.CloneTree(0)
            n = t.GetEntries()	
            ndd56 = array('f', [-99])
            clonet.Branch("n2ddt56", ndd56, 'n2ddt56'+'/F')
            ndd53 = array('f', [-99])
            clonet.Branch("n2ddt53", ndd53, 'n2ddt53'+'/F')
            for j in range(0,int(n)):
                t.GetEntry(j)
                pt = t.fj1Pt
                mass = t.fj1MSD_corr
                rho = log(mass*mass/pt/pt)
                jtN2b1sd_8 = t.n2

                rind6 = f56.GetXaxis().FindBin(rho)
                pind6 = f56.GetYaxis().FindBin(pt)
                    
                rind3 = f53.GetXaxis().FindBin(rho)
                pind3 = f53.GetYaxis().FindBin(pt)		
		
                if rho >  f56.GetXaxis().GetBinUpEdge( f56.GetXaxis().GetNbins() ) :
                    rind6 = f56.GetXaxis().GetNbins()
                if rho <  f56.GetXaxis().GetBinLowEdge( 1 ) :
                    rind6 = 1
                if pt >  f56.GetYaxis().GetBinUpEdge( f56.GetYaxis().GetNbins() ) :
                    pind6 = f56.GetYaxis().GetNbins()
                if pt < f56.GetYaxis().GetBinLowEdge( 1 ) :
                    pind6 = 1

                if rho > f53.GetXaxis().GetBinUpEdge( f53.GetXaxis().GetNbins() ) :
                    rind3 = f53.GetXaxis().GetNbins()
                if rho <  f53.GetXaxis().GetBinLowEdge( 1 ) :
                    rind3 = 1
                if pt > f53.GetYaxis().GetBinUpEdge( f53.GetYaxis().GetNbins() ) :
                    pind3 = f53.GetYaxis().GetNbins()
                if pt < f53.GetYaxis().GetBinLowEdge( 1 ) :
                    pind3 = 1
                
                ndd56[0] = jtN2b1sd_8 - f56.GetBinContent(rind6,pind6)
                ndd53[0] = jtN2b1sd_8 - f53.GetBinContent(rind3,pind3)

                clonet.Fill()
            fo.Write()
            
        fo.Close()
        fi.Close()
	trans.Close()
        os.system("mv -f %s %s" % (output, rootfile))


def shift_btags(additional=None):
    shifted_weights = {}
    #if not any([x in region for x in ['signal','top','w']]):
    #    return shifted_weights 
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_btag']:
            shiftedlabel = ''
            if 'sj' in cent:
                shiftedlabel += 'sj'
            if 'B' in shift:
                shiftedlabel += 'btag'
            else:
                shiftedlabel += 'mistag'
            if 'Up' in shift:
                shiftedlabel += 'Up'
            else:
                shiftedlabel += 'Down'
            weight = sel.weights[region+'_'+cent+shift]%lumi
            if additional:
                weight = tTIMES(weight,additional)
            shifted_weights[shiftedlabel] = weight
    return shifted_weights


#vmap definition
vmap = {}
mc_vmap = {'genBosonPt':'genBosonPt'}
if region in ['signal','test']:
    u,uphi, = ('puppimet','dphipuppimet')
elif 'pho' in region:
    u,uphi = ('puppiUAmag','dphipuppiUA')
elif 'wmn'or 'wen' or 'te' or 'tm' in region:
    u,uphi = ('puppiUWmag','dphipuppiUW')
elif 'zee'  or 'zmn' in region:
    u,uphi = ('puppiUZmag','dphipuppiUZ')
vmap['met'] = 'min(%s,999.9999)'%u 
vmap['fj1Pt'] = 'fj1Pt'
vmap['fj1MSD_corr'] = 'fj1MSD_corr'
vmap['fj1ECFN_2_3_10'] = 'fj1ECFN_2_3_10'
vmap['fj1ECFN_1_2_10'] = 'fj1ECFN_1_2_10'
vmap['n2'] = 'fj1ECFN_2_3_10/pow(fj1ECFN_1_2_10,2.00)'
#vmap['ndd56'] = 'n2ddt56'
#vmap['ndd53'] = 'n2ddt53'

#weights
weights = {'nominal' : sel.weights[region]%lumi}
#btag shifts to be added here, maybe

ps ={'ZtoNuNu','ZJets','WJets','SingleTop','Diboson','QCD','SingleElectron','TTbar','MET'}
#ps ={'ZtoNuNu','ZJets','WJets','SingleTop','Diboson','QCD','SingleElectron','TTbar','MET','GJets', 'SinglePhoton'}

forestDir = basedir + '/rhalphabeth/'
os.system('mkdir -p %s/%s'%(forestDir,region))

for rfiles in os.listdir(basedir):
	if '.root' in rfiles: rfile = rfiles.split('.')[0]
	if rfile not in ps: continue
	if  rfile == 'SingleElectron' and ('zee' in region or 'te' in region or 'wen' in region):	
	 	print rfile,0
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , "Events" , is_data=True, extra_cut=sel.eleTrigger )
		factory.run(forestDir+'/%s/fittingForest_%s.root'%(region,rfile))
	if rfile == 'MET' and ('zmm' in region or 'tm' in region or 'wmn' in region or 'signal' in region):
		print rfile,1
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , "Events" , is_data=True, extra_cut=sel.metTrigger )
                factory.run(forestDir+'/%s/fittingForest_%s.root'%(region,rfile))
	if rfile == 'SinglePhoton' and  'pho' in region:
		print rfile,2
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , "Events" , is_data=True)
                factory.run(forestDir+'/%s/fittingForest_%s.root'%(region,rfile))
	if rfile != 'SingleElectron' and rfile != 'SinglePhoton' and rfile != 'MET':
		print rfile,3
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , "Events" )
                factory.run(forestDir+'/%s/fittingForest_%s.root'%(region,rfile))
        else: continue
	
pd = {'ZpA0h_med-600_dm-300','ZpA0h_med-800_dm-300','ZpA0h_med-1000_dm-300','ZpA0h_med-1200_dm-300','ZpA0h_med-1400_dm-300','ZpA0h_med-1700_dm-300','ZpA0h_med-2000_dm-300','ZpA0h_med-2500_dm-300','ZpBaryonic_med-10_dm-1','ZpBaryonic_med-10_dm-10','ZpBaryonic_med-10_dm-50','ZpBaryonic_med-10_dm-150','ZpBaryonic_med-10_dm-500','ZpBaryonic_med-15_dm-10','ZpBaryonic_med-20_dm-1','ZpBaryonic_med-50_dm-1','ZpBaryonic_med-50_dm-10','ZpBaryonic_med-50_dm-50','ZpBaryonic_med-95_dm-50','ZpBaryonic_med-100_dm-1','ZpBaryonic_med-100_dm-10','ZpBaryonic_med-200_dm-1','ZpBaryonic_med-200_dm-50','ZpBaryonic_med-200_dm-150','ZpBaryonic_med-295_dm-150','ZpBaryonic_med-300_dm-1','ZpBaryonic_med-300_dm-50','ZpBaryonic_med-500_dm-1','ZpBaryonic_med-500_dm-150','ZpBaryonic_med-995_dm-500','ZpBaryonic_med-1000_dm-1','ZpBaryonic_med-1000_dm-150','ZpBaryonic_med-1000_dm-1000','ZpBaryonic_med-1995_dm-1000','ZpBaryonic_med-2000_dm-1','ZpBaryonic_med-2000_dm-500','ZpBaryonic_med-10000_dm-1','ZpBaryonic_med-10000_dm-150','ZpBaryonic_med-10000_dm-500','ZpBaryonic_med-10000_dm-10','ZpBaryonic_med-10000_dm-50'}

if region in ['signal_scalar','signal_vector','signal_thq','signal_stdm','signal','signal_fail']:

    for rfiles in os.listdir(basedir):
        if '.root' in rfiles: rfile = rfiles.split('.')[0]
        if rfile not in pd: continue
        #print rfile
        factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                                       cut = sel.cuts[region],
                                       variables = vmap,
                                       mc_variables = mc_vmap,
                                       mc_weights = weights)
        factory.add_process(f(rfile) , "Events" )
        factory.run(forestDir+'/%s/fittingForest_%s.root'%(region,rfile))


#Computing n2ddt variables in the ntuples stored inside the basedir directory
if nddt:
    for roo in os.listdir(forestDir+"/"+region):
        if '.root' in roo:
                print roo
                addN2DDT(forestDir+'/'+region+'/%s'%roo)

