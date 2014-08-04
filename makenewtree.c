/*********************************
**********************************
** Created by Brian Francisco 	**
**	on 07-11-2014	            	**
**			                      	**
** Loops over file	          	**
** "OmniTree_BprimeBprimeTo..."	**
** and creates new, simplified	**
** tree for Bprime analysis	    **
**			                      	**
**********************************
*********************************/

#include <TFile.h>
#include <TTree.h>
#include <math.h>
#include <TROOT.h>

//Declare variables to be filled by old tree and to fill new tree

Int_t bPartons;		//b partons per event
Int_t bPartons30;	//b partons with pT>30, |eta|<2.4
Int_t nGenPart;		//Particles per event
Int_t nPart;
Int_t nPFJets;		//Jets per event
Int_t nJets;
Int_t nJets30;		//Jets with pT>30, |eta|<2.4
Int_t nHiggs;		//Higgs per event
Int_t partID[50];	//Particle IDs
Int_t pdgID[50];

Float_t MCpx[50];	//px, py and pz of all particles
Float_t MCpy[50];	//used to calculate pT, theta, eta and phi	
Float_t MCpz[50];	
Float_t MCpt[50];	//pT of all particles
Float_t b_pt[10];	//pT of b partons
Float_t Higgs_pt[2];	//pT of Higgs
Float_t MCeta[50];	//eta value (not abs val)
Float_t b_eta[10];	//eta for b partons
Float_t Higgs_eta[2];	//eta for Higgs
Float_t MCphi[50];	//phi
Float_t b_phi[10];	//b parton phi
Float_t Higgs_phi[2];	//Higgs phi
Float_t Ht;		//Ht of event
Int_t MCmotherind[50];	//Index of particle's mother
Int_t Motherindex[50];
Int_t b_mothers[10];
Int_t Higgs_mother[2];
Float_t jet_PF_pt[45];		//Jet's pT
Float_t jet_PF_eta[45];		//eta
Float_t jet_PF_phi[45];		//phi
Float_t jet_PF_px[45];
Float_t jet_PF_py[45];
Float_t bdiscCSV_PF[45];	//b discriminate value
Float_t jet_pt[45];
Float_t jet_eta[45];
Float_t jet_phi[45];
Float_t jet_bdisc[45];

void makenewtree()
{  
// Define input and output files  
   TFile infile("OmniTree_BprimeBprimeToBHBHinc_M-650_TuneZ2star_8TeV-madgraph_tree.root");
   TFile ofile("BprimeAnalysis_650.root","RECREATE");

// Point to old tree, declare new tree
   TTree *oldtree = infile.GetObjectChecked("EvTree","TTree");	
   TTree *newtree = new TTree("BprimeTree","Smaller Tree for B prime Analysis");

// Create new tree's branches
   newtree->Branch("nPart",&nPart,"nPart/I");
   newtree->Branch("bPartons",&bPartons,"bPartons/I");
   newtree->Branch("bPartons30",&bPartons30,"bPartons30/I");
   newtree->Branch("nJets",&nJets,"nJets/I");
   newtree->Branch("nJets30",&nJets30,"nJets30/I");
   newtree->Branch("nHiggs",&nHiggs,"nHiggs/I");
   newtree->Branch("partID[nPart]",partID,"partID[nPart]/I");
   newtree->Branch("Motherindex[nPart]",Motherindex,"Motherindex[nPart]/I");
   newtree->Branch("b_mothers[bPartons]",b_mothers,"b_mothers[bPartons]/I");
   newtree->Branch("Higgs_mother[nHiggs]",Higgs_mother,"Higgs_mother[nHiggs]/I");
   newtree->Branch("MCpt[nPart]",MCpt,"MCpt[nPart]/F");
   newtree->Branch("b_pt[bPartons]",b_pt,"b_pt[bPartons]/F");
   newtree->Branch("Higgs_pt[nHiggs]",Higgs_pt,"Higgs_pt[nHiggs]/F");
   newtree->Branch("MCeta[nPart]",MCeta,"MCeta[nPart]/F");
   newtree->Branch("b_eta[bPartons]",b_eta,"b_eta[bPartons]/F");
   newtree->Branch("Higgs_eta[nHiggs]",Higgs_eta,"Higgs_eta[nHiggs]/F");
   newtree->Branch("MCphi[nPart]",MCphi,"MCphi[nPart]/F");
   newtree->Branch("b_phi[bPartons]",b_phi,"b_phi[bPartons]/F");
   newtree->Branch("Higgs_phi[nHiggs]",Higgs_phi,"Higgs_phi[nHiggs]/F");
   newtree->Branch("jet_pt[nJets]",jet_pt,"jet_pt[nJets]/F");
   newtree->Branch("jet_eta[nJets]",jet_eta,"jet_eta[nJets]/F");
   newtree->Branch("jet_phi[nJets]",jet_phi,"jet_phi[nJets]/F");
   newtree->Branch("jet_bdisc[nJets]",jet_bdisc,"jet_bdisc[nJets]/F");
   newtree->Branch("Ht",&Ht,"Ht/F");

// Address old tree's branches
   oldtree->SetBranchAddress("nGenPart", &nGenPart);
   oldtree->SetBranchAddress("nPFJets", &nPFJets);
   oldtree->SetBranchAddress("MCpx[nGenPart]",MCpx);
   oldtree->SetBranchAddress("MCpy[nGenPart]",MCpy);
   oldtree->SetBranchAddress("MCpz[nGenPart]",MCpz);
   oldtree->SetBranchAddress("pdgID[nGenPart]",pdgID);
   oldtree->SetBranchAddress("MCmotherind[nGenPart]",MCmotherind);
   oldtree->SetBranchAddress("jet_PF_pt[nPFJets]",jet_PF_pt);
   oldtree->SetBranchAddress("jet_PF_eta[nPFJets]",jet_PF_eta);
   oldtree->SetBranchAddress("jet_PF_phi[nPFJets]",jet_PF_phi);
   oldtree->SetBranchAddress("jet_PF_px[nPFJets]",jet_PF_px);
   oldtree->SetBranchAddress("jet_PF_py[nPFJets]",jet_PF_py);
   oldtree->SetBranchAddress("bdiscCSV_PF[nPFJets]",bdiscCSV_PF);

// Loop over each event
   for(int i = 0; oldtree->GetEntry(i) > 0; ++i)
   {	oldtree->GetEntry(i);
	bPartons = 0;
	bPartons30 = 0;
	nHiggs = 0;
	Ht = 0;	
	nPart = nGenPart;
	nJets = nPFJets;
	nJets30 = 0;

//     	Loop over each particle	
	for(int j = 0; j < nGenPart; ++j)
	{	partID[j] = pdgID[j];					//copies particle ID array
		
		MCpt[j] = sqrt(pow(MCpx[j],2)+pow(MCpy[j],2));		//calculate pT
	  float p_mag = sqrt(pow(MCpx[j],2)+pow(MCpy[j],2)+pow(MCpz[j],2));	//magnitude of 3p vector
		MCeta[j] = .5*log((p_mag + MCpz[j])/(p_mag - MCpz[j]));			//calculate eta
		MCphi[j] = atan2(MCpy[j],MCpx[j]);			//calculate phi
		
		Motherindex[j] = MCmotherind[j];			//copies mother index array

		if(fabs(pdgID[j]) == 5)	
		{	b_pt[bPartons] = MCpt[j];		//pt, eta and phi of b partons
			b_eta[bPartons] = MCeta[j];
			b_phi[bPartons] = MCphi[j];
			b_mothers[bPartons] = pdgID[MCmotherind[j]];	//stores ID of mother, rather than index
			++bPartons;		//count b partons
			if((MCpt[j]>30)&&(fabs(MCeta[j])<2.4)) ++bPartons30;
		}
		if(fabs(pdgID[j]) == 25) 
		{	Higgs_pt[nHiggs] = MCpt[j];		//pt, eta and phi of Higgs
			Higgs_eta[nHiggs] = MCeta[j];
			Higgs_phi[nHiggs] = MCphi[j];
			Higgs_mother[nHiggs] = pdgID[MCmotherind[j]];	//stores ID of mother, rather than index
			++nHiggs;		//count Higgs
		}
	}

//	Loop over each jet
	for(int j = 0; j < nPFJets; ++j)
	{	jet_pt[j] = jet_PF_pt[j]; 	//Copy values from old tree
		jet_eta[j] = jet_PF_eta[j];
		jet_phi[j] = jet_PF_phi[j];
		jet_bdisc[j] = bdiscCSV_PF[j];	

		if((jet_PF_pt[j] > 30.0)&&(fabs(jet_PF_eta[j])<2.4))	
		{	Ht += jet_PF_pt[j];	//Ht is sum of pT of jets with pT > 30 GeV
			++nJets30;
		}
	}
	
//	Fill new tree
	newtree->Fill();
   }

// Write new tree to output file
   newtree->Print();
   newtree->Write();

// Write output file
   ofile.Write();
}
