#!/usr/bin/env python
import sys
import os
import ROOT 
import math
import shutil
from array import array
import warnings
from ROOT import TFile
from ROOT import TDirectory
ROOT.gROOT.SetBatch(True)
from optparse import OptionParser

#usage: ./write_regression_systematic.py path

#os.mkdir(path+'/sys')
argv = sys.argv
parser = OptionParser()
parser.add_option("-P", "--path", dest="path", default="", 
                      help="path to samples")
#parser.add_option("-S", "--samples", dest="names", default="", 
#                      help="samples you want to run on")
#parser.add_option("-C", "--config", dest="config", default=[], action="append",
#                      help="configuration defining the plots to make")
(opts, args) = parser.parse_args(argv)

from myutils import BetterConfigParser, ParseInfo

ROOT.gSystem.Load("libPhysicsToolsKinFitter.so")
KinFitNameSpace=('/gpfs/ddn/cms/user/cvernier/KinFit/CMSSW_5_3_3_patch2/src/PhysicsTools/KinFitter/test/kinFit4b_h.so')
ROOT.gSystem.Load(KinFitNameSpace)
from ROOT import H4b
pathIN = opts.path
pathOUT = opts.path+'out/' 
job ='Selected_flatTree.root'
#job = 'flatTree_VBF-Powheg125_preselect_hard_tmva.root'
jobdir = 'Hbb'
jobtree = 'events'
print 'INput samples:\t%s'%pathIN
print 'OUTput samples:\t%s'%pathOUT


def deltaPhi(phi1, phi2): 
    result = phi1 - phi2
    while (result > math.pi): result -= 2*math.pi
    while (result <= -math.pi): result += 2*math.pi
    return result


def deltaR(eta1,phi1,eta2,phi2):
      deta = eta1 - eta2;
      dphi = deltaPhi(phi1, phi2);
      return math.sqrt(deta*deta + dphi*dphi)

#def Xmass():

ROOT.gROOT.ProcessLine(
        "struct H {\
        float         chi;\
        float         chi2;\
	float         mass;\
        float         pt;\
        float         eta;\
        float         phi;\
        float         dR;\
        float         dPhi;\
        float         dEta;\
        } ;"
    )
    
    
input = ROOT.TFile.Open(pathIN+job,'read')
output = ROOT.TFile.Open(pathOUT+job,'recreate')

input.cd()
#dir = input.Get(jobdir)
#tree = dir.Get(jobtree)
tree = input.Get(jobtree)
nEntries = tree.GetEntries()
print nEntries	
H = ROOT.H()
tree.SetBranchStatus('H',0)
output.cd()
newtree = tree.CloneTree(0)

#hJ0 = ROOT.TLorentzVector()
#hJ1 = ROOT.TLorentzVector()
#hJ2 = ROOT.TLorentzVector()
#hJ3 = ROOT.TLorentzVector()
#hJ4 = ROOT.TLorentzVector()
#hJ5 = ROOT.TLorentzVector()
#H_ = ROOT.TLorentzVector()
hJet_e = array('f',[0]*20) 
hJet_pt = array('f',[0]*20)
hJet_eta = array('f',[0]*20)
hJet_phi = array('f',[0]*20)
hJet_PU = array('f',[0]*20)
hJet_btag = array('f',[0]*20)
hJet_reg = array('f',[0]*20)
h_njets = ROOT.TH1F('njet','njets', 10, 0, 10)
h_eta = ROOT.TH1F('eta','eta',40, -5, 5)	
h_btag = ROOT.TH1F('btag','btag', 20, -1,1)
h_pt = ROOT.TH1F('pt','pt', 20, 0., 100.)
#chi = ROOT.Double()
cont = 0
newtree.Branch( 'H', H , 'chi/F:chi2/F:mass/F:pt/F:eta/F:phi/F:dR/F:dPhi/F:dEta/F' )
for entry in range(0,nEntries):
	 
	     Bid = []	
	     index =[]	
	     hJ0 = ROOT.TLorentzVector()
	     hJ1 = ROOT.TLorentzVector()
             hJ2 = ROOT.TLorentzVector()
             hJ3 = ROOT.TLorentzVector()
             hJ4 = ROOT.TLorentzVector()
             hJ5 = ROOT.TLorentzVector()
             chi = ROOT.Double()

	     tree.GetEntry(entry)
             Bid.append(tree.b1)
	     Bid.append(tree.b2)
	     Bid.append(tree.q1)
	     Bid.append(tree.q2)			
	     #print tree.nJets
	     for i in range(0,tree.nJets):
	       hJet_e[i] = tree.jetMass[i]
               hJet_pt[i] = tree.jetPt[i]
               hJet_eta[i] = tree.jetEta[i]
               hJet_phi[i] = tree.jetPhi[i]
	       hJet_PU[i] = tree.jetPuIdL[i]	
	       hJet_btag[i] = tree.jetBtag[i]	
	       hJet_reg[i] = tree.jetReg[i]     	
	       if i not in Bid and len(Bid)<5 : Bid.append(i)	
#	       print tree.jetBtag[i]
	   #    BTag.insert(i,tree.jetBtag[i])
           #    BTagORD.insert(i,tree.jetBtag[i])	
          #    for i in range(0,5):
	#	BTagORD.sort()
	#	BTagORD.reverse()
		#print '%s, %s, %s'%(BTagORD[i],BTag.index(BTagORD[i]),BTag[i])
	#	Bid.insert(i, BTag.index(BTagORD[i]))
#		if BTagORD[i] <0 : 
#			BTag.insert(i, -2.*i)
#			BTag.remove(BTagORD[i])
	     #print BTag
	     #print BTagORD	
#	     print Bid	
	     #print Bid.index(4)	
             hJ0.SetPtEtaPhiM(hJet_pt[Bid[0]],hJet_eta[Bid[0]],hJet_phi[Bid[0]],hJet_e[Bid[0]])
             hJ1.SetPtEtaPhiM(hJet_pt[Bid[1]],hJet_eta[Bid[1]],hJet_phi[Bid[1]],hJet_e[Bid[1]])
	     hJ0= hJ0*tree.jetReg[Bid[0]]	
	     hJ1= hJ1*tree.jetReg[Bid[1]]	
             #hJ0.SetPtEtaPhiM(hJet_pt[Bid[0]],hJet_eta[Bid[0]],hJet_phi[Bid[0]],hJet_e[Bid[0]])
             #hJ1.SetPtEtaPhiM(hJet_pt[Bid[1]],hJet_eta[Bid[1]],hJet_phi[Bid[1]],hJet_e[Bid[1]])
	     #if deltaPhi(hJ0.Phi(),hJ1.Phi())>2 : continue     
             hJ2.SetPtEtaPhiM(hJet_pt[Bid[2]],hJet_eta[Bid[2]],hJet_phi[Bid[2]],hJet_e[Bid[2]])
             hJ3.SetPtEtaPhiM(hJet_pt[Bid[3]],hJet_eta[Bid[3]],hJet_phi[Bid[3]],hJet_e[Bid[3]])
             #if(tree.nJets >4): hJ4.SetPtEtaPhiM(hJet_pt[Bid[4]],hJet_eta[Bid[4]],hJet_phi[Bid[4]],hJet_e[Bid[4]])

	     Px = hJ2.Px()+hJ3.Px()#-tree.met*math.cos(tree.metPhi)
	     Py = hJ2.Py()+hJ3.Py()#-tree.met*math.sin(tree.metPhi)	
	     for i in range(0,20) : index.insert(i,i)
	     #print index
	     #print Bid
	     for i in range(0,4)  : index.remove(Bid[i])
	     #print index
	     njet=0		
	     for i in range(0,8):
	      if hJet_PU[index[i]]>0. and hJet_pt[index[i]]>10 and abs( hJet_eta[index[i]])<2.5 and njet <5 :
		#and hJet_pt[index[i]]<50: 
		hJ4.SetPtEtaPhiM(hJet_pt[index[i]],hJet_eta[index[i]],hJet_phi[index[i]],hJet_e[index[i]])
		h_btag.Fill(hJet_btag[index[i]])
		h_eta.Fill(hJet_eta[index[i]])
		h_pt.Fill(hJet_pt[index[i]])
		if hJet_btag[index[i]]>0.244 : hJ4 = hJ4*hJet_reg[index[i]]	
		Px = Px + hJ4.Px() 
		Py = Py + hJ4.Py()
		njet= njet+1
	     h_njets.Fill(njet)
		
	     #METett = ROOT.TTreeFormula('METett','met',tree)
             #METphii = ROOT.TTreeFormula('METphii','metPhi',tree)	
	     #METet = METett.EvalInstance()
	     #METphi = METphii.EvalInstance()	
	     #metX = -METet*math.cos(METphi)
             #metY = -METet*math.sin(METphi)
	     H.eta = abs(hJ0.Eta()-hJ1.Eta())
	     H.mass = (hJ0+hJ1).M()
		
	     if H.mass-tree.mbbReg>0.1 :
				cont = cont+1	
			        print Bid
             			#print '2 %s,%s'%(Bid.index(0),Bid.index(1))
             			#print '3 btag %s, %s, %s, %s: '%(tree.jetBtag[Bid.index(0)],tree.jetBtag[Bid.index(1)],tree.jetBtag[Bid.index(2)],tree.jetBtag[Bid.index(3)])     
             			#print '4 btag %s, %s, %s, %s: '%(tree.jetBtag[0],tree.jetBtag[1],tree.jetBtag[2],tree.jetBtag[3]) 
			#	print '5 bindx %s, %s, %s, %s, %s: '%(tree.btagIdx_[0],tree.btagIdx_[1],tree.btagIdx_[2],tree.btagIdx_[3],tree.btagIdx_[4])
	#		print 'my mass %s, tree mass %s '%(H.mass,tree.mbb)	
	#		print 'deta : %s'%H.eta
	#		print 'b ordered %s,%s eta ordered: %s,%s '%(Bid.index(0),Bid.index(1),Etaid.index(0),Etaid.index(1))
             H.pt = (hJ0+hJ1).Pt()

             H.phi = (hJ0+hJ1).Phi()
             H.dR = hJ0.DeltaR(hJ1)
             H.dPhi = hJ0.DeltaPhi(hJ1)
             H.dEta = abs(hJ0.Eta()-hJ1.Eta())	
	     H.chi= (H4b.Xchi2(hJ0,hJ1,-Px,-Py,chi,0.)).M()    
	     H.chi2 = chi 
	     #if chi>300 : H.chi = H.mass		
	    # H_.SetPtEtaPhiE(H.pt,H.eta,H.phi,(hJ0+hJ1).E())	
	     #del H			
	     del Bid		
             newtree.Fill()
h_eta.Write()
h_pt.Write()
h_njets.Write()
h_btag.Write()
print cont
print 'Exit loop'
newtree.AutoSave()
print 'Save'
output.Close()
print 'Close'
