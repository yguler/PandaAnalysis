#
import os
import ROOT
from ROOT import *
from array import array
import math
from math import *
import sys
import pdb

def ComputeDDT(name, point, nPtBins, nRhoBins, H):
	DDT = TH2F(name, "", nRhoBins, -6, -1.5, nPtBins, 200, 800)
	DDT.SetStats(0)
	nXb = H.GetXaxis().GetNbins()
	nYb = H.GetYaxis().GetNbins()
	for x in range(nXb):
		for y in range(nYb):
			proj = H.ProjectionZ("H3"+str(x)+str(y),x+1,x+1,y+1,y+1)
			print str(x+1) + "," + str(y+1) + ":    "+ str(proj.Integral())
			p = array('d', [point])
			q = array('d', [0.0]*len(p))
			proj.GetQuantiles(len(p), q, p)
			DDT.SetBinContent( x+1, y+1, q[0] );
	return DDT

def DisplayDDT(DDT, Title, SaveName):
	CMSLABL = TLatex()
	CMSLABL.SetNDC()
	CMSLABL.SetTextSize(0.045)
	PRELABL = TLatex()
	PRELABL.SetNDC()
	PRELABL.SetTextSize(0.04)
	THILABL = TLatex()
	THILABL.SetNDC()
	THILABL.SetTextSize(0.045)

	C = TCanvas("TempCanvas", "Title", 605, 495)
	C.cd()
	DDT.SetStats(0)
	DDT.GetXaxis().SetTitle("jet #rho")
	DDT.GetXaxis().SetTitleSize(0.045)
	DDT.GetYaxis().SetTitle("jet p_{T}")
	DDT.GetYaxis().SetTitleSize(0.045)
	DDT.GetYaxis().SetTitleOffset(1.145)
	DDT.Draw("COLZ")
	CMSLABL.DrawLatex(0.1465,0.85,"CMS")
	THILABL.DrawLatex(0.81,0.91,"#bf{13 TeV}")
	PRELABL.DrawLatex(0.1465,0.812,"#bf{#it{Simulation Preliminary}}")

	C.Print("MAP_"+SaveName+".gif")

H3 = TH3F("H3", "", 9, -6, -1.5, 12, 200, 800, 500, 0, 0.5)
H3.SetStats(0)

Bkgs =["../../zprimebits/zprimegamma/GJetsmvaEVv3.root"]
for B in Bkgs:
	F = TFile(B)
	T = F.Get("Events")
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		T.GetEntry(j)
		weight = T.puWeight*T.scale1fb*T.kfactor
		PT = T.AK8Puppijet0_pt
		if PT>0.: 
			preRHO = T.AK8Puppijet0_msd*T.AK8Puppijet0_msd/T.AK8Puppijet0_pt/T.AK8Puppijet0_pt
		else:
			continue
		if preRHO > 0.:
			RHO = math.log(preRHO)
			if PT > 200 and RHO < -1.5 and RHO > -6.0 and T.AK8Puppijet0_N2sdb1 > 0.:
				H3.Fill(RHO, PT, T.AK8Puppijet0_N2sdb1, weight)

DDT_5by6 = ComputeDDT("DDT", 0.05, 12, 9, H3)
DisplayDDT(DDT_5by6, "DDT vals at 5%", "DDT")

Fout = TFile("PhotonDDTs.root", "recreate")
Fout.cd()
DDT_5by6.Write()
Fout.Close()
for bks,B in Bkgs_tags:

	H3[bks] = TH3F("H3_%s"%(bks), "H3_%s"%(bks), 9, -6, -1.5, 12, 200, 800, 500, 0, 0.5)
	
for bks,B in Bkgs_tags:

	print 'Starting with '+bks+'-------------------------- :)'

	H3[bks].SetStats(0)
	F = TFile(B)
	if "test" in B:
		tree = "Diboson_test"
	T = F.Get("%s"%tree)
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
                if(j % (1 * n/100) == 0):
                        sys.stdout.write("\r[" + "="*int(20*j/n) + " " + str(round(100.*j/n,0)) + "% done")
                        sys.stdout.flush()
		T.GetEntry(j)
		#weight = T.puWeight*T.scale1fb*T.kfactor*T.kfactorNLO
		weight = 1
		PT = T.fj1Pt
		preRHO = T.fj1MSD_corr*T.fj1MSD_corr/T.fj1Pt/T.fj1Pt
		if preRHO > 0. and T.fj1ECFN_1_2_10 != 0.:
			RHO = math.log(preRHO)
			jtN2b1sd_8 = T.fj1ECFN_2_3_10/math.pow(T.fj1ECFN_1_2_10,2.00)
			if PT > 200 and RHO < -1.5 and RHO > -6.0 and jtN2b1sd_8 > 0.:
				H3[bks].Fill(RHO, PT, jtN2b1sd_8, weight)
DDT_5by6={}
DDT_5by3={}
Fout = TFile("DDTs.root", "recreate")
Fout.cd()
for key in H3:
	DDT_5by6[key]=ComputeDDT('DDT_5by6_%s'%(key), 0.05, 12, 9, H3[key])
	DDT_5by6[key].Write()
	DDT_5by3[key]=ComputeDDT('DDT_5by3_%s'%(key), 0.05, 3, 9, H3[key])
	DDT_5by3[key].Write()
Fout.Close()
