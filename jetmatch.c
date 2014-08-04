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

TH1F *deltaR = new TH1F("deltaR_values","deltaR values",500,0,4.5);
TH1F *bdisc = new TH1F("bdisc_onematch","b disc values of jets in 1:1 matches",1000,0,1);
TH1F *bdisc_all = new TH1F("all_bdiscs","b disc values of all matched jets",1000,0,1);
TH1F *nomatch_pt = new TH1F("nomatch_pt","pt of b partons with no matched jet",100,15,725);
TH1F *nomatch_eta = new TH1F("nomatch_eta","|eta| values of b partons with no match",100,0,3);
TH1I *jetsperparton = new TH1I("jetsperparton","jets matched per parton",4,0,4);

void jetmatch()
{	
   Int_t bmatchind[200];
   Int_t jetmatchind[200];
   int matchindex;
   int bmatch;
   int repeats;
   int skipmatches[500];
   bool skipped;

/***************************
Be sure to match mass point
***************************/   
   TFile infile("BprimeAnalysis_650.root");
   TFile ofile("jetmatch/BprimeJetMatching_650.root","RECREATE");
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
 	if((bPartons30 == 6) && ((nJets30 == 6) || (nJets30 == 7)) && (Ht > 1000))	//Require 6 b partons and 6 or 7 jets w/ eta phi cut
	{    matchindex = 0;
	     repeats = 0;
	     skipped = false;
	     for(int i = 0; i < nPart; ++i)		//Loop over all particles
	     { 	   if((fabs(partID[i]) == 5) && (MCpt[i] > 30) && fabs((MCeta[i]) < 2.4))//If they are b partons with eta and pt cut
         	   {	bmatch = 0;		//Counts how many times this b parton is matched to a jet
			for(int k = 0; k < nJets; ++k)			//Loop over all jets
			{   if((jet_pt[k] > 30) && (fabs(jet_eta[k]) < 2.4))	//pT and eta requirement
			    {	deltaR->Fill(matcher(jet_eta[k],jet_phi[k],MCeta[i],MCphi[i]));	//Histo deltaR values
				if(matcher(jet_eta[k],jet_phi[k],MCeta[i],MCphi[i]) < 0.1)//Do jet matching 
			    	{   bmatchind[matchindex] = i;	
			    	    jetmatchind[matchindex] = k;
			    	    ++matchindex;
				    ++bmatch;	//Counts how many times this b parton is matched to a jet
				    bdisc_all->Fill(jet_bdisc[k]);	//Record bdisc of every jet matched to b parton
			    	}
			    }
			}
			if(bmatch == 0)	//Histo pt and |eta| of b's with no matched jet
			{	nomatch_pt->Fill(MCpt[i]);
				nomatch_eta->Fill(fabs(MCeta[i]));
			}
			jetsperparton->Fill(bmatch);
		   }
	     }

/***The following code eliminates repeated jet matches, rather than keeping the one with best deltaR***/

	    for(int k = 0; k < (matchindex-1); ++k)
	    {	for(int l = (k+1); l < matchindex; ++l)
		{   if(jetmatchind[k] == jetmatchind[l])
		    {	skipmatches[repeats] = k;	//Create an array of indices to be skipped in jetmatchind[]
			skipmatches[repeats+1] = l;
			++repeats;
			++repeats;
		    }
		}
	    }
	    for(int k = 0; k < matchindex; ++k)
	    {	for(int l = 0; l < repeats; ++l)	//If index matches an index in array of indices to skip
		{   if(k == skipmatches[l])		//Assign it as "skipped"
		    {	skipped = true;
		    }
		}
		if(skipped == false)			//If it is not assigned as "skipped", histo its bdisc
		{   bdisc->Fill(jet_bdisc[jetmatchind[k]]);	//Interested in bdisc of 1:1 jet matches, hence repeat elimination
		}
	    }

/***End repeat elimination***/

	}
   }

//Write histograms and root file
   deltaR->Write();
   bdisc->Write();
   bdisc_all->Write();
   nomatch_pt->Write();
   nomatch_eta->Write();
   jetsperparton->Write();
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
 
