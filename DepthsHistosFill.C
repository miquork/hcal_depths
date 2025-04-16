// Purpose: calibrate HCAL depths using tag-and-probe with dijet events.

#define DepthsHistosFill_cxx
#include "DepthsHistosFill.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TMath.h>

#include <iostream>

// UTILITIES
double DELTAPHI(double phi1, double phi2)
{
  double dphi = fabs(phi1 - phi2);
  return (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}
double DELTAR(double phi1, double phi2, double eta1, double eta2)
{
  return sqrt(pow(DELTAPHI(phi1, phi2), 2) + pow(eta1 - eta2, 2));
}


void DepthsHistosFill::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DepthsHistosFill.C
//      root> DepthsHistosFill t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();

   fChain->SetBranchStatus("*",0);

   fChain->SetBranchStatus("Flag_METFilters",1);
   fChain->SetBranchStatus("PuppiMET_pt",1);
   fChain->SetBranchStatus("PuppiMET_sumEt",1);
   
   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_jetId",1);

   fChain->SetBranchStatus("nJetPFCand",1);
   fChain->SetBranchStatus("JetPFCand_pfCandIdx",1);
   fChain->SetBranchStatus("JetPFCand_jetIdx",1);

   fChain->SetBranchStatus("PFCand_hcalDepthEnergyFraction1",1);
   fChain->SetBranchStatus("PFCand_hcalDepthEnergyFraction2",1);
   fChain->SetBranchStatus("PFCand_hcalDepthEnergyFraction3",1);
   fChain->SetBranchStatus("PFCand_hcalDepthEnergyFraction4",1);
   fChain->SetBranchStatus("PFCand_hcalDepthEnergyFraction5",1);
   fChain->SetBranchStatus("PFCand_hcalDepthEnergyFraction6",1);
   fChain->SetBranchStatus("PFCand_hcalDepthEnergyFraction7",1);
   //fChain->SetBranchStatus("PFCand_ecalEnergy",1);
   fChain->SetBranchStatus("PFCand_rawEcalEnergy",1);
   //fChain->SetBranchStatus("PFCand_hcalEnergy",1);
   fChain->SetBranchStatus("PFCand_rawHcalEnergy",1);
   //fChain->SetBranchStatus("PFCand_hoEnergy",1);
   fChain->SetBranchStatus("PFCand_rawHoEnergy",1);
   
   TFile *fout = new TFile("DepthsHistosFill.root","RECREATE");

   // Regular L2Relative eta binning
   double vx[] =
     {-5.191,
      -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
      -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
      -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
      -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435,
      -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
      0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
      1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
      2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
      4.363, 4.538, 4.716, 4.889, 5.191};
   const int nx = sizeof(vx) / sizeof(vx[0]) - 1;

   // Control for checking if HT_2000 is enough for pT>1100 GeV
   // Inclusive jets pT binning
   double vpti[] =
     {1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
      97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
      507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
      1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
      2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
      4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
   double npti = sizeof(vpti) / sizeof(vpti[0]) - 1;
   
   TH1D *heta = new TH1D("heta",";#eta_{probe} (GeV);N_{tag}",nx,vx);

   TProfile *peta = new TProfile("peta",";#eta_{probe};p_{T,probe}/p_{T,tag}",
				 nx,vx);

   const int ndepth = 1 + 7 + 1; // ECAL, 7xHCAL, HO
   TProfile* vpeta[ndepth];
   for (int depth = 0; depth != ndepth; ++depth) {
     vpeta[depth] = new TProfile(Form("peta_depth%d",depth),
				 Form(";#eta_{probe};p_{T,probe} / p_{T,tag} "
				      " #times f_{depth%d}",depth),nx,vx);
   } // depth


   TH1D *hpt = new TH1D("hpt",";p_{T,tag} (GeV);N_{tag}",npti,vpti);
   
   TProfile *ppt = new TProfile("ppt",";p_{T,tag};p_{T,probe}/p_{T,tag}",
				npti,vpti);
   
   TProfile* vppt[ndepth];
   for (int depth = 0; depth != ndepth; ++depth) {
     vppt[depth] = new TProfile(Form("ppt_depth%d",depth),
				Form(";p_{T,tag};p_{T,probe} / p_{T,tag} "
				     " #times f_{depth%d}",depth),npti,vpti);
   } // depth

   
   cout << "Processing DepthsHistosFill.C over " << nentries << " events"
	<< endl << flush;
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (jentry%10000==0) cout << "." << flush;
      
      // Both choices for tag and probe jets
      for (int itag = 0; itag != 2; ++itag) {
	int iprobe = (itag==0 ? 1 : 0);
      
	// Basic event selection cuts
	bool passEvent = (Flag_METFilters==1 &&
			  PuppiMET_pt < 0.3*PuppiMET_sumEt);
	bool passTnP = (nJet>=2 &&
			fabs(Jet_eta[itag])<1.3 && //fabs(Jet_pt[itag])>1200. &&
			DELTAPHI(Jet_phi[itag],Jet_phi[iprobe])>2.7 &&
			Jet_jetId[itag]>=4 && Jet_jetId[iprobe]>=4 &&
			(nJet<3 || Jet_pt[2]<30 || Jet_pt[2]<0.3*Jet_pt[itag]));

	if (passEvent && passTnP) {

	  // Calculate depth fractions for tag jet
	  double fdepth[ndepth];
	  for (int i = 0; i != ndepth; ++i) fdepth[i] = 0;

	  for (int ilink = 0; ilink != nJetPFCand; ++ilink) {

	    if (JetPFCand_jetIdx[ilink]==iprobe) {
	      int ipf = JetPFCand_pfCandIdx[ilink];

	      fdepth[0] += PFCand_rawEcalEnergy[ipf];
	      double rawhcal = PFCand_rawHcalEnergy[ipf];
	      fdepth[1] += rawhcal * PFCand_hcalDepthEnergyFraction1[ipf];
	      fdepth[2] += rawhcal * PFCand_hcalDepthEnergyFraction2[ipf];
	      fdepth[3] += rawhcal * PFCand_hcalDepthEnergyFraction3[ipf];
	      fdepth[4] += rawhcal * PFCand_hcalDepthEnergyFraction4[ipf];
	      fdepth[5] += rawhcal * PFCand_hcalDepthEnergyFraction5[ipf];
	      fdepth[6] += rawhcal * PFCand_hcalDepthEnergyFraction6[ipf];
	      fdepth[7] += rawhcal * PFCand_hcalDepthEnergyFraction7[ipf];
	      fdepth[8] += PFCand_rawHoEnergy[ipf];
	    } // ipf in iprobe
	  } // for ilink

	  // Calculate total raw energy
	  double sumraw(0);
	  for (int i = 0; i != ndepth; ++i) sumraw += fdepth[i];
	  
	  // Normalize fractions by total energy
	  if (sumraw!=0) for (int i = 0; i != ndepth; ++i) fdepth[i] /= sumraw;

	  double DB = Jet_pt[iprobe] / Jet_pt[itag];
	  if (Jet_pt[itag]>1200.)        {
	    heta->Fill(Jet_eta[iprobe]);
	    peta->Fill(Jet_eta[iprobe], DB);
	  }
	  if (fabs(Jet_eta[iprobe])<1.3) {
	    hpt->Fill(Jet_pt[itag]);
	    ppt->Fill(Jet_pt[itag], DB);
	  }
	  
	  for (int depth = 0; depth != ndepth; ++depth) {

	    if (Jet_pt[itag]>1200.) 
	      vpeta[depth]->Fill(Jet_eta[iprobe], DB * fdepth[depth]);
	    if (fabs(Jet_eta[iprobe])<1.3) 
	      vppt[depth]->Fill(Jet_pt[itag], DB * fdepth[depth]);
	  }
	}
      } // for itag
   }

   fout->Write();
   fout->Close();
}
