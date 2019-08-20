#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>

Int_t nPart;
Int_t bPartons30;
Int_t nJets;
Int_t nJets30;
Int_t partID[200];

Float_t jet_pt[200];
Float_t jet_eta[200];
Float_t jet_phi[200];
Float_t MCeta[200];
Float_t MCphi[200];
Float_t MCpt[200];
Float_t jet_bdisc[200];
Float_t Ht;

TH1I *bjets = new TH1I("nBJets","number of bjets per event (only 1:1 match jets, with 'loose' cut)",8,0,8);

void sixjetmatch()
{	
   Int_t bmatchind[200];
   Int_t jetmatchind[200];
   int matchindex;
   int bmatch;
   int b_jet;
   int onetoone;
   int repeats;
   int skipmatches[500];
   bool skipped;

/***************************
Be sure to match mass point
***************************/   
   TFile infile("BprimeAnalysis_650.root");
   TFile ofile("jetmatch/6Jets_BprimeJetMatching_650.root","RECREATE");
/***************************
***************************/

   TTree *btree = infile.GetObjectChecked("BprimeTree","TTree");

   btree->SetBranchAddress("nPart", &nPart);
   btree->SetBranchAddress("bPartons30", &bPartons30);
   btree->SetBranchAddress("nJets",&nJets);
   btree->SetBranchAddress("nJets30", &nJets30);
   btree->SetBranchAddress("partID[nPart]",partID);
   btree->SetBranchAddress("jet_pt[nJets]",jet_pt);
   btree->SetBranchAddress("jet_eta[nJets]",jet_eta);
   btree->SetBranchAddress("jet_phi[nJets]",jet_phi);
   btree->SetBranchAddress("MCeta[nPart]",MCeta);
   btree->SetBranchAddress("MCphi[nPart]",MCphi);
   btree->SetBranchAddress("MCpt[nPart]",MCpt);
   btree->SetBranchAddress("jet_bdisc[nJets]",jet_bdisc);
   btree->SetBranchAddress("Ht",&Ht);

   for(int j = 0; btree->GetEntry(j) > 0; ++j)
   {	btree->GetEntry(j);
 	if((bPartons30 == 6) && ((nJets30 == 6) || (nJets30 == 7)) && (Ht > 1000))//Require 6 b parts and 6 or 7 jets (w/ cut)
	{   matchindex = 0;
	    b_jet = 0;
	    repeats = 0;
	    skipped = false;
	    for(int i = 0; i < nPart; ++i)		//Loop over all particles
	    { 	   if((fabs(partID[i]) == 5) && (MCpt[i] > 30) && fabs((MCeta[i]) < 2.4))//If they are b partons (eta & pt cut)
         	   {	bmatch = 0;		//Counts how many times this b parton is matched to a jet
			for(int k = 0; k < nJets; ++k)			//Loop over all jets
			{   if((jet_pt[k] > 30) && (fabs(jet_eta[k]) < 2.4))	//pT and eta requirement
			    {	if(matcher(jet_eta[k],jet_phi[k],MCeta[i],MCphi[i]) < 0.1)//Do jet matching 
			    	{   bmatchind[matchindex] = i;	
			    	    jetmatchind[matchindex] = k;
			    	    ++matchindex;
				    ++bmatch;	//Counts how many times this b parton is matched to a jet
//				    if(jet_bdisc[k] >= 0.244)	//Count as bjet if bdisc makes cutoff (loose = 0.244)
//				    {	++b_jet;
//				    }
			    	}
			    }
			}
		   }
	    }
//	    bjets->Fill(b_jet);	//Histo number of bjets

/***The following code takes repeat matches and keeps the one with the best deltaR***/
	    onetoone = matchindex;
	    for(int k = 0; k < (matchindex-1); ++k)
	    {	for(int l = (k+1); l < matchindex; ++l)
		{   if(jetmatchind[k] == jetmatchind[l])
		    {	if(matcher(jet_eta[jetmatchind[k]],jet_phi[jetmatchind[k]],MCeta[bmatchind[k]],MCphi[bmatchind[k]]) > matcher(jet_eta[jetmatchind[l]],jet_phi[jetmatchind[l]],MCeta[bmatchind[l]],MCphi[bmatchind[l]]))
			{   jetmatchind[k] = jetmatchind[l];   
			}
			for(int p = l; p < matchindex-1; ++p)
			{   jetmatchind[p] = jetmatchind[p+1];
			}
			--onetoone;
		    }
		}
	    }
/***End repeat elimination***/
	    if(onetoone == 6)
	    {   for(int k = 0; k < 6; ++k)
		{   if(jet_bdisc[jetmatchind[k]] >= 0.244)
		    {	++b_jet;
		    }
		}
	    }
	    if(b_jet != 0)
	    {	bjets->Fill(b_jet);
	    }
	}
   }

//Write histograms and root file
   bjets->Write();
   ofile.Write();
}

float matcher(float jeta, float jphi, float beta, float bphi) //Calculates and returns DeltaR
{
   float pi = 4*atan(1.0);
   float dPhi = fabs(bphi - jphi);
   float dEta = fabs(beta - jeta);

   if((2*pi - dPhi) < dPhi) dPhi = ((2*pi) - dPhi); //Take minimum distance in phi (periodic)

   return sqrt(pow(dEta,2)+pow(dPhi,2));
}
 
