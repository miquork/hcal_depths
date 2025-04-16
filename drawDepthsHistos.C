// Purpose: draw comparison between baseline and new gains (vs MC)
#include "TFile.h"
#include "TProfile.h"

#include "tdrstyle_mod22.C"

// HB should be zero for deptha 4-7, but is not exactlyx, so zero explicitly
void resetHB(TH1D *h) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (fabs(h->GetXaxis()->GetBinCenter(i))<1.305) {
      h->SetBinContent(i, 0.);
      h->SetBinError(i, 0.);
    }
  }
} // resetHB

void drawDepthsHistos() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  TFile *fb = new TFile("DepthsHistosFill_Baseline_v3.root","READ");
  assert(fb && !fb->IsZombie());

  TFile *fu = new TFile("DepthsHistosFill_UpdatedCorrs_v3.root","READ");
  assert(fu && !fu->IsZombie());

  TFile *fm = new TFile("DepthsHistosFill_HT2000_v3.root","READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();

  double eps = 1e-4;
  TH1D *h = tdrHist("h","p_{T,probe} / p_{T,tag}"
		    " #times f_{depth,probe}",0+eps,0.68,
		    "#eta",-2.322,2.322);
  TH1D *hd = tdrHist("hd","Data / MC",0,8,//0.6,4.8,//2.4,
		     "#eta",-2.322,2.322);
  extraText = "Private";
  lumi_136TeV = "2024Iv2 JetMET1 (private re-reco)";
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);

  c1->cd(1);
  h->GetYaxis()->SetTitleOffset(1.00);//1.05);
  
  const int ndepth = 1 + 7 + 1; // ECAL, HCAL, HO
  TLegend *leg1 = tdrLeg(0.33,0.53-3*0.05,0.58,0.53);
  leg1->SetTextSize(0.045);
  TLegend *leg2 = tdrLeg(0.50,0.59-ndepth*0.035,0.80,0.59);
  leg2->SetTextSize(0.035);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  tex->DrawLatex(0.33,0.85,"Tag p_{T}>1200 GeV, |#eta|<1.3, #Delta#phi>2.7, #alpha<0.3");

  int color[ndepth] =
    {kBlue+2,  // ECAL
     kBlue, kOrange+1, kGreen+1, kRed, // HB
     kYellow+1, kOrange+3, kGray+1, // HE
     kGray+2}; // HO
  int marker[ndepth] =
    {kFullDiamond, // ECAL
     kOpenTriangleDown, kOpenSquare, kOpenCircle, kOpenTriangleUp, //HB
     kOpenStar, kOpenDiamond, kOpenCross, // HE
     kFullStar}; // HO
  double size[ndepth] =
    {1.0,
     0.7,0.7,0.7,0.7,
     0.7,0.7,0.7,
     1.0};
  
  for (int depth = 0; depth != ndepth; ++depth) {
    
    TProfile *pb = (TProfile*)fb->Get(Form("peta_depth%d",depth)); assert(pb);
    TProfile *pu = (TProfile*)fu->Get(Form("peta_depth%d",depth)); assert(pu);
    TProfile *pm = (TProfile*)fm->Get(Form("peta_depth%d",depth)); assert(pm);

    TH1D *hr = pu->ProjectionX(Form("hreta_depth%d",depth));
    hr->Divide(pb);

    TH1D *hrb = pb->ProjectionX(Form("hretab_depth%d",depth));
    hrb->Divide(pm);
    if (depth>4 && depth!=8) resetHB(hrb);
    
    TH1D *hru = pu->ProjectionX(Form("hretau_depth%d",depth));
    hru->Divide(pm);
    if (depth>4 && depth!=8) resetHB(hru);
    
    c1->cd(1);
    tdrDraw(pm,"HIST",kNone,color[depth],kSolid,-1,kNone,0,size[depth]);
    tdrDraw(pb,"HIST",kNone,color[depth],kDotted,-1,kNone,0,size[depth]);
    tdrDraw(pu,"Pz",marker[depth],color[depth],kSolid,-1,kNone,0,size[depth]);

    if (depth==2) {
      leg1->AddEntry(pu,"New","PLE");
      leg1->AddEntry(pb,"Base","L");
      leg1->AddEntry(pm,"MC","L");
    }
    if (depth==0)      leg2->AddEntry(pu,"ECAL","PLE");
    else if (depth<5)  leg2->AddEntry(pu,Form("HB/HE depth %d",depth),"PLE");
    else if (depth<8)  leg2->AddEntry(pu,Form("HE depth %d",depth),"PL");
    else if (depth==8) leg2->AddEntry(pu,"HO","PL");

    c1->cd(2);
    //tdrDraw(hr,depth<5 ? "Pz" : "HISTPz",
    //	    marker[depth],color[depth],kSolid,-1,kNone,0,0.7);
    tdrDraw(hrb,"HIST",kNone,color[depth],kSolid,-1,kNone,0,size[depth]);
    tdrDraw(hru,depth<5||depth==8 ? "Pz" : "HISTPz",
	    marker[depth],color[depth],kSolid,-1,kNone,0,size[depth]);
  } // for depth

  c1->cd(1);
  gPad->RedrawAxis();

  c1->cd(2);
  gPad->RedrawAxis();
  
  c1->SaveAs("pdf/drawDepthsHistos.pdf");



  double eps2 = 1e-1;
  TH1D *h_2 = tdrHist("h_2","p_{T,probe} / p_{T,tag}"
		      " #times f_{depth,probe}",0,0.72,
		      "p_{T,tag} (GeV)",1000-eps2,3000);
  TH1D *hd_2 = tdrHist("hd_2","Data / MC",0.6,4.8,//2.4,
		       "p_{T,tag} (GeV)",1000-eps2,3000);
  //extraText = "Private";
  //lumi_136TeV = "2024Iv2 JetMET1";
  TCanvas *c2 = tdrDiCanvas("c2",h_2,hd_2,8,11);

  c2->cd(1);
  gPad->SetLogx();
  h_2->GetYaxis()->SetTitleOffset(1.05);
  
  TLegend *leg2_1 = tdrLeg(0.33,0.50-3*0.05,0.58,0.50);
  leg2_1->SetTextSize(0.045);
  TLegend *leg2_2 = tdrLeg(0.50,0.55-(ndepth-3)*0.05,0.80,0.55);
  leg2_2->SetTextSize(0.045);

  tex->DrawLatex(0.33,0.85,"|#eta_{tag,probe}|<1.3, #Delta#phi>2.7, #alpha<0.3");

  c2->cd(2);
  gPad->SetLogx();

  for (int depth = 0; depth != ndepth; ++depth) {
    
    TProfile *pb = (TProfile*)fb->Get(Form("ppt_depth%d",depth)); assert(pb);
    TProfile *pu = (TProfile*)fu->Get(Form("ppt_depth%d",depth)); assert(pu);
    TProfile *pm = (TProfile*)fm->Get(Form("ppt_depth%d",depth)); assert(pm);

    TH1D *hr = pu->ProjectionX(Form("hrpt_depth%d",depth));
    hr->Divide(pb);

    TH1D *hrb = pb->ProjectionX(Form("hrptb_depth%d",depth));
    hrb->Divide(pm);
    if (depth>4 && depth!=8) resetHB(hrb);
    
    TH1D *hru = pu->ProjectionX(Form("hrptu_depth%d",depth));
    hru->Divide(pm);
    if (depth>4 && depth!=8) resetHB(hru);
    
    c2->cd(1);
    if (depth>=5 && depth<=7) continue;
    tdrDraw(pm,"HIST",kNone,color[depth],kSolid,-1,kNone,0,size[depth]);
    tdrDraw(pb,"HIST",kNone,color[depth],kDotted,-1,kNone,0,size[depth]);
    tdrDraw(pu,"Pz",marker[depth],color[depth],kSolid,-1,kNone,0,size[depth]);

    if (depth==2) {
      leg2_1->AddEntry(pu,"New","PLE");
      leg2_1->AddEntry(pb,"Base","L");
      leg2_1->AddEntry(pm,"MC","L");
    }
    if (depth==0)      leg2_2->AddEntry(pu,"ECAL","PLE");
    else if (depth<5)  leg2_2->AddEntry(pu,Form("HB depth %d",depth),"PLE");
    //else if (depth<8)  leg2_2->AddEntry(pu,Form("HE depth %d",depth),"PL");
    else if (depth==8) leg2_2->AddEntry(pu,"HO","PL");

    c2->cd(2);
    //tdrDraw(hr,depth<5 ? "Pz" : "HISTPz",
    //	    marker[depth],color[depth],kSolid,-1,kNone,0,0.7);
    tdrDraw(hrb,"HIST",kNone,color[depth],kSolid,-1,kNone,0,size[depth]);
    tdrDraw(hru,depth<5||depth==8 ? "Pz" : "HISTPz",
	    marker[depth],color[depth],kSolid,-1,kNone,0,size[depth]);
  } // for depth

  c2->cd(1);
  gPad->RedrawAxis();

  c2->cd(2);
  gPad->RedrawAxis();
  
  c2->SaveAs("pdf/drawDepthsHistos_vsPtTag.pdf");
  
} // drawDepthsHistos()
