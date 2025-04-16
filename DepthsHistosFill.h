//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  4 09:03:22 2025 by ROOT version 6.28/06
// from TTree Events/Events
// found on file: ../data/dijet/data/NanoMerged_Data24Iv2_JetMET1_TeVJet_Baseline_RAWToPFNANO_v0p1.root
//////////////////////////////////////////////////////////

#ifndef DepthsHistosFill_h
#define DepthsHistosFill_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DepthsHistosFill {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          bunchCrossing;
   UInt_t          orbitNumber;
   Short_t         BeamSpot_type;
   Float_t         BeamSpot_sigmaZ;
   Float_t         BeamSpot_sigmaZError;
   Float_t         BeamSpot_z;
   Float_t         BeamSpot_zError;
   Int_t           nboostedTau;
   UChar_t         boostedTau_idAntiEle2018[3];   //[nboostedTau]
   UChar_t         boostedTau_idAntiMu[3];   //[nboostedTau]
   UChar_t         boostedTau_idMVAnewDM2017v2[3];   //[nboostedTau]
   UChar_t         boostedTau_idMVAoldDM2017v2[3];   //[nboostedTau]
   Short_t         boostedTau_jetIdx[3];   //[nboostedTau]
   Short_t         boostedTau_rawAntiEleCat2018[3];   //[nboostedTau]
   Int_t           boostedTau_charge[3];   //[nboostedTau]
   Int_t           boostedTau_decayMode[3];   //[nboostedTau]
   Float_t         boostedTau_chargedIso[3];   //[nboostedTau]
   Float_t         boostedTau_eta[3];   //[nboostedTau]
   Float_t         boostedTau_leadTkDeltaEta[3];   //[nboostedTau]
   Float_t         boostedTau_leadTkDeltaPhi[3];   //[nboostedTau]
   Float_t         boostedTau_leadTkPtOverTauPt[3];   //[nboostedTau]
   Float_t         boostedTau_mass[3];   //[nboostedTau]
   Float_t         boostedTau_neutralIso[3];   //[nboostedTau]
   Float_t         boostedTau_phi[3];   //[nboostedTau]
   Float_t         boostedTau_photonsOutsideSignalCone[3];   //[nboostedTau]
   Float_t         boostedTau_pt[3];   //[nboostedTau]
   Float_t         boostedTau_puCorr[3];   //[nboostedTau]
   Float_t         boostedTau_rawAntiEle2018[3];   //[nboostedTau]
   Float_t         boostedTau_rawIso[3];   //[nboostedTau]
   Float_t         boostedTau_rawIsodR03[3];   //[nboostedTau]
   Float_t         boostedTau_rawMVAnewDM2017v2[3];   //[nboostedTau]
   Float_t         boostedTau_rawMVAoldDM2017v2[3];   //[nboostedTau]
   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Int_t           nCorrT1METJet;
   Float_t         CorrT1METJet_EmEF[46];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_area[46];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_eta[46];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_muonSubtrFactor[46];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_phi[46];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_rawPt[46];   //[nCorrT1METJet]
   Float_t         DeepMETResolutionTune_phi;
   Float_t         DeepMETResolutionTune_pt;
   Float_t         DeepMETResponseTune_phi;
   Float_t         DeepMETResponseTune_pt;
   Int_t           nElectron;
   Bool_t          Electron_convVeto[11];   //[nElectron]
   UChar_t         Electron_cutBased[11];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[11];   //[nElectron]
   Bool_t          Electron_isEB[11];   //[nElectron]
   Bool_t          Electron_isEcalDriven[11];   //[nElectron]
   Bool_t          Electron_isPFcand[11];   //[nElectron]
   UChar_t         Electron_jetNDauCharged[11];   //[nElectron]
   UChar_t         Electron_lostHits[11];   //[nElectron]
   Bool_t          Electron_mvaIso_WP80[11];   //[nElectron]
   Bool_t          Electron_mvaIso_WP90[11];   //[nElectron]
   Bool_t          Electron_mvaIso_WPHZZ[11];   //[nElectron]
   Bool_t          Electron_mvaNoIso_WP80[11];   //[nElectron]
   Bool_t          Electron_mvaNoIso_WP90[11];   //[nElectron]
   UChar_t         Electron_seedGain[11];   //[nElectron]
   UChar_t         Electron_tightCharge[11];   //[nElectron]
   Short_t         Electron_jetIdx[11];   //[nElectron]
   Short_t         Electron_photonIdx[11];   //[nElectron]
   Short_t         Electron_seediEtaOriX[11];   //[nElectron]
   Short_t         Electron_seediPhiOriY[11];   //[nElectron]
   Short_t         Electron_svIdx[11];   //[nElectron]
   Short_t         Electron_fsrPhotonIdx[11];   //[nElectron]
   Int_t           Electron_charge[11];   //[nElectron]
   Int_t           Electron_pdgId[11];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[11];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmapHEEP[11];   //[nElectron]
   Float_t         Electron_PreshowerEnergy[11];   //[nElectron]
   Float_t         Electron_deltaEtaSC[11];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[11];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[11];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[11];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[11];   //[nElectron]
   Float_t         Electron_dxy[11];   //[nElectron]
   Float_t         Electron_dxyErr[11];   //[nElectron]
   Float_t         Electron_dz[11];   //[nElectron]
   Float_t         Electron_dzErr[11];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[11];   //[nElectron]
   Float_t         Electron_ecalEnergy[11];   //[nElectron]
   Float_t         Electron_ecalEnergyError[11];   //[nElectron]
   Float_t         Electron_energyErr[11];   //[nElectron]
   Float_t         Electron_eta[11];   //[nElectron]
   Float_t         Electron_fbrem[11];   //[nElectron]
   Float_t         Electron_gsfTrketaMode[11];   //[nElectron]
   Float_t         Electron_gsfTrkpMode[11];   //[nElectron]
   Float_t         Electron_gsfTrkpModeErr[11];   //[nElectron]
   Float_t         Electron_gsfTrkphiMode[11];   //[nElectron]
   Float_t         Electron_hoe[11];   //[nElectron]
   Float_t         Electron_ip3d[11];   //[nElectron]
   Float_t         Electron_jetPtRelv2[11];   //[nElectron]
   Float_t         Electron_jetRelIso[11];   //[nElectron]
   Float_t         Electron_mass[11];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[11];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[11];   //[nElectron]
   Float_t         Electron_mvaHZZIso[11];   //[nElectron]
   Float_t         Electron_mvaIso[11];   //[nElectron]
   Float_t         Electron_mvaNoIso[11];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[11];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[11];   //[nElectron]
   Float_t         Electron_phi[11];   //[nElectron]
   Float_t         Electron_pt[11];   //[nElectron]
   Float_t         Electron_r9[11];   //[nElectron]
   Float_t         Electron_rawEnergy[11];   //[nElectron]
   Float_t         Electron_scEtOverPt[11];   //[nElectron]
   Float_t         Electron_sieie[11];   //[nElectron]
   Float_t         Electron_sip3d[11];   //[nElectron]
   Float_t         Electron_promptMVA[11];   //[nElectron]
   Int_t           nFatJet;
   UChar_t         FatJet_jetId[6];   //[nFatJet]
   UChar_t         FatJet_nConstituents[6];   //[nFatJet]
   Short_t         FatJet_chMultiplicity[6];   //[nFatJet]
   Short_t         FatJet_neMultiplicity[6];   //[nFatJet]
   Short_t         FatJet_subJetIdx1[6];   //[nFatJet]
   Short_t         FatJet_subJetIdx2[6];   //[nFatJet]
   Short_t         FatJet_electronIdx3SJ[6];   //[nFatJet]
   Short_t         FatJet_muonIdx3SJ[6];   //[nFatJet]
   Float_t         FatJet_area[6];   //[nFatJet]
   Float_t         FatJet_chEmEF[6];   //[nFatJet]
   Float_t         FatJet_chHEF[6];   //[nFatJet]
   Float_t         FatJet_eta[6];   //[nFatJet]
   Float_t         FatJet_mass[6];   //[nFatJet]
   Float_t         FatJet_msoftdrop[6];   //[nFatJet]
   Float_t         FatJet_muEF[6];   //[nFatJet]
   Float_t         FatJet_n2b1[6];   //[nFatJet]
   Float_t         FatJet_n3b1[6];   //[nFatJet]
   Float_t         FatJet_neEmEF[6];   //[nFatJet]
   Float_t         FatJet_neHEF[6];   //[nFatJet]
   Float_t         FatJet_particleNetLegacy_QCD[6];   //[nFatJet]
   Float_t         FatJet_particleNetLegacy_Xbb[6];   //[nFatJet]
   Float_t         FatJet_particleNetLegacy_Xcc[6];   //[nFatJet]
   Float_t         FatJet_particleNetLegacy_Xqq[6];   //[nFatJet]
   Float_t         FatJet_particleNetLegacy_mass[6];   //[nFatJet]
   Float_t         FatJet_particleNetWithMass_H4qvsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNetWithMass_HbbvsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNetWithMass_HccvsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNetWithMass_QCD[6];   //[nFatJet]
   Float_t         FatJet_particleNetWithMass_TvsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNetWithMass_WvsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNetWithMass_ZvsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_QCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_QCD0HF[6];   //[nFatJet]
   Float_t         FatJet_particleNet_QCD1HF[6];   //[nFatJet]
   Float_t         FatJet_particleNet_QCD2HF[6];   //[nFatJet]
   Float_t         FatJet_particleNet_XbbVsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_XccVsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_XggVsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_XqqVsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_XteVsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_XtmVsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_XttVsQCD[6];   //[nFatJet]
   Float_t         FatJet_particleNet_massCorr[6];   //[nFatJet]
   Float_t         FatJet_phi[6];   //[nFatJet]
   Float_t         FatJet_pt[6];   //[nFatJet]
   Float_t         FatJet_rawFactor[6];   //[nFatJet]
   Float_t         FatJet_tau1[6];   //[nFatJet]
   Float_t         FatJet_tau2[6];   //[nFatJet]
   Float_t         FatJet_tau3[6];   //[nFatJet]
   Float_t         FatJet_tau4[6];   //[nFatJet]
   Float_t         FatJet_lsf3[6];   //[nFatJet]
  static const int nMaxJetPFCand = 10000;
   Int_t           nJetPFCand;
   Int_t           JetPFCand_pfCandIdx[nMaxJetPFCand];   //[nJetPFCand]
   Int_t           JetPFCand_jetIdx[nMaxJetPFCand];   //[nJetPFCand]
   Int_t           nFatJetPFCand;
   Int_t           FatJetPFCand_pfCandIdx[467];   //[nFatJetPFCand]
   Int_t           FatJetPFCand_jetIdx[467];   //[nFatJetPFCand]
   Int_t           nFsrPhoton;
   Short_t         FsrPhoton_electronIdx[4];   //[nFsrPhoton]
   Short_t         FsrPhoton_muonIdx[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_dROverEt2[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_eta[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_phi[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_pt[4];   //[nFsrPhoton]
   Float_t         FsrPhoton_relIso03[4];   //[nFsrPhoton]
   Int_t           nIsoTrack;
   Bool_t          IsoTrack_isHighPurityTrack[23];   //[nIsoTrack]
   Bool_t          IsoTrack_isPFcand[23];   //[nIsoTrack]
   Bool_t          IsoTrack_isFromLostTrack[23];   //[nIsoTrack]
   Short_t         IsoTrack_charge[23];   //[nIsoTrack]
   Short_t         IsoTrack_fromPV[23];   //[nIsoTrack]
   Int_t           IsoTrack_pdgId[23];   //[nIsoTrack]
   Float_t         IsoTrack_dxy[23];   //[nIsoTrack]
   Float_t         IsoTrack_dz[23];   //[nIsoTrack]
   Float_t         IsoTrack_eta[23];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_all[23];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_chg[23];   //[nIsoTrack]
   Float_t         IsoTrack_phi[23];   //[nIsoTrack]
   Float_t         IsoTrack_pt[23];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_all[23];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_chg[23];   //[nIsoTrack]
  static const int nMaxJet = 500;
   Int_t           nJet;
   UChar_t         Jet_chMultiplicity[nMaxJet];   //[nJet]
   UChar_t         Jet_jetId[nMaxJet];   //[nJet]
   UChar_t         Jet_nConstituents[nMaxJet];   //[nJet]
   UChar_t         Jet_nElectrons[nMaxJet];   //[nJet]
   UChar_t         Jet_nMuons[nMaxJet];   //[nJet]
   UChar_t         Jet_nSVs[nMaxJet];   //[nJet]
   UChar_t         Jet_neMultiplicity[nMaxJet];   //[nJet]
   Short_t         Jet_electronIdx1[nMaxJet];   //[nJet]
   Short_t         Jet_electronIdx2[nMaxJet];   //[nJet]
   Short_t         Jet_muonIdx1[nMaxJet];   //[nJet]
   Short_t         Jet_muonIdx2[nMaxJet];   //[nJet]
   Short_t         Jet_svIdx1[nMaxJet];   //[nJet]
   Short_t         Jet_svIdx2[nMaxJet];   //[nJet]
   Int_t           Jet_hfadjacentEtaStripsSize[nMaxJet];   //[nJet]
   Int_t           Jet_hfcentralEtaStripSize[nMaxJet];   //[nJet]
   Float_t         Jet_PNetRegPtRawCorr[nMaxJet];   //[nJet]
   Float_t         Jet_PNetRegPtRawCorrNeutrino[nMaxJet];   //[nJet]
   Float_t         Jet_PNetRegPtRawRes[nMaxJet];   //[nJet]
   Float_t         Jet_UParTAK4RegPtRawCorr[nMaxJet];   //[nJet]
   Float_t         Jet_UParTAK4RegPtRawCorrNeutrino[nMaxJet];   //[nJet]
   Float_t         Jet_UParTAK4RegPtRawRes[nMaxJet];   //[nJet]
   Float_t         Jet_area[nMaxJet];   //[nJet]
   Float_t         Jet_btagDeepFlavB[nMaxJet];   //[nJet]
   Float_t         Jet_btagDeepFlavCvB[nMaxJet];   //[nJet]
   Float_t         Jet_btagDeepFlavCvL[nMaxJet];   //[nJet]
   Float_t         Jet_btagDeepFlavQG[nMaxJet];   //[nJet]
   Float_t         Jet_btagPNetB[nMaxJet];   //[nJet]
   Float_t         Jet_btagPNetCvB[nMaxJet];   //[nJet]
   Float_t         Jet_btagPNetCvL[nMaxJet];   //[nJet]
   Float_t         Jet_btagPNetCvNotB[nMaxJet];   //[nJet]
   Float_t         Jet_btagPNetQvG[nMaxJet];   //[nJet]
   Float_t         Jet_btagPNetTauVJet[nMaxJet];   //[nJet]
   Float_t         Jet_btagUParTAK4B[nMaxJet];   //[nJet]
   Float_t         Jet_btagUParTAK4CvB[nMaxJet];   //[nJet]
   Float_t         Jet_btagUParTAK4CvL[nMaxJet];   //[nJet]
   Float_t         Jet_btagUParTAK4CvNotB[nMaxJet];   //[nJet]
   Float_t         Jet_btagUParTAK4QvG[nMaxJet];   //[nJet]
   Float_t         Jet_btagUParTAK4TauVJet[nMaxJet];   //[nJet]
   Float_t         Jet_chEmEF[nMaxJet];   //[nJet]
   Float_t         Jet_chHEF[nMaxJet];   //[nJet]
   Float_t         Jet_eta[nMaxJet];   //[nJet]
   Float_t         Jet_hfEmEF[nMaxJet];   //[nJet]
   Float_t         Jet_hfHEF[nMaxJet];   //[nJet]
   Float_t         Jet_hfsigmaEtaEta[nMaxJet];   //[nJet]
   Float_t         Jet_hfsigmaPhiPhi[nMaxJet];   //[nJet]
   Float_t         Jet_mass[nMaxJet];   //[nJet]
   Float_t         Jet_muEF[nMaxJet];   //[nJet]
   Float_t         Jet_muonSubtrFactor[nMaxJet];   //[nJet]
   Float_t         Jet_neEmEF[nMaxJet];   //[nJet]
   Float_t         Jet_neHEF[nMaxJet];   //[nJet]
   Float_t         Jet_phi[nMaxJet];   //[nJet]
   Float_t         Jet_pt[nMaxJet];   //[nJet]
   Float_t         Jet_rawFactor[nMaxJet];   //[nJet]
   Int_t           nLowPtElectron;
   Bool_t          LowPtElectron_convVeto[9];   //[nLowPtElectron]
   UChar_t         LowPtElectron_convWP[9];   //[nLowPtElectron]
   UChar_t         LowPtElectron_lostHits[9];   //[nLowPtElectron]
   Short_t         LowPtElectron_electronIdx[9];   //[nLowPtElectron]
   Short_t         LowPtElectron_photonIdx[9];   //[nLowPtElectron]
   Int_t           LowPtElectron_charge[9];   //[nLowPtElectron]
   Int_t           LowPtElectron_pdgId[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_ID[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_convVtxRadius[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_deltaEtaSC[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dxy[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dxyErr[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dz[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_dzErr[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_eInvMinusPInv[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_energyErr[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_eta[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_hoe[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_mass[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_miniPFRelIso_all[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_miniPFRelIso_chg[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_phi[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_pt[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_ptbiased[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_r9[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_scEtOverPt[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_sieie[9];   //[nLowPtElectron]
   Float_t         LowPtElectron_unbiased[9];   //[nLowPtElectron]
   Int_t           nProton_multiRP;
   UChar_t         Proton_multiRP_arm[2];   //[nProton_multiRP]
   Float_t         Proton_multiRP_t[2];   //[nProton_multiRP]
   Float_t         Proton_multiRP_thetaX[2];   //[nProton_multiRP]
   Float_t         Proton_multiRP_thetaY[2];   //[nProton_multiRP]
   Float_t         Proton_multiRP_time[2];   //[nProton_multiRP]
   Float_t         Proton_multiRP_timeUnc[2];   //[nProton_multiRP]
   Float_t         Proton_multiRP_xi[2];   //[nProton_multiRP]
   Int_t           nMuon;
   UChar_t         Muon_highPtId[27];   //[nMuon]
   Bool_t          Muon_highPurity[27];   //[nMuon]
   Bool_t          Muon_inTimeMuon[27];   //[nMuon]
   Bool_t          Muon_isGlobal[27];   //[nMuon]
   Bool_t          Muon_isPFcand[27];   //[nMuon]
   Bool_t          Muon_isStandalone[27];   //[nMuon]
   Bool_t          Muon_isTracker[27];   //[nMuon]
   UChar_t         Muon_jetNDauCharged[27];   //[nMuon]
   Bool_t          Muon_looseId[27];   //[nMuon]
   Bool_t          Muon_mediumId[27];   //[nMuon]
   Bool_t          Muon_mediumPromptId[27];   //[nMuon]
   UChar_t         Muon_miniIsoId[27];   //[nMuon]
   UChar_t         Muon_multiIsoId[27];   //[nMuon]
   UChar_t         Muon_mvaMuID_WP[27];   //[nMuon]
   UChar_t         Muon_nStations[27];   //[nMuon]
   UChar_t         Muon_nTrackerLayers[27];   //[nMuon]
   UChar_t         Muon_pfIsoId[27];   //[nMuon]
   UChar_t         Muon_puppiIsoId[27];   //[nMuon]
   Bool_t          Muon_softId[27];   //[nMuon]
   Bool_t          Muon_softMvaId[27];   //[nMuon]
   UChar_t         Muon_tightCharge[27];   //[nMuon]
   Bool_t          Muon_tightId[27];   //[nMuon]
   UChar_t         Muon_tkIsoId[27];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[27];   //[nMuon]
   Short_t         Muon_jetIdx[27];   //[nMuon]
   Short_t         Muon_svIdx[27];   //[nMuon]
   Short_t         Muon_fsrPhotonIdx[27];   //[nMuon]
   Int_t           Muon_charge[27];   //[nMuon]
   Int_t           Muon_pdgId[27];   //[nMuon]
   Float_t         Muon_dxy[27];   //[nMuon]
   Float_t         Muon_dxyErr[27];   //[nMuon]
   Float_t         Muon_dxybs[27];   //[nMuon]
   Float_t         Muon_dz[27];   //[nMuon]
   Float_t         Muon_dzErr[27];   //[nMuon]
   Float_t         Muon_eta[27];   //[nMuon]
   Float_t         Muon_ip3d[27];   //[nMuon]
   Float_t         Muon_jetPtRelv2[27];   //[nMuon]
   Float_t         Muon_jetRelIso[27];   //[nMuon]
   Float_t         Muon_mass[27];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[27];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[27];   //[nMuon]
   Float_t         Muon_mvaMuID[27];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[27];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[27];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[27];   //[nMuon]
   Float_t         Muon_phi[27];   //[nMuon]
   Float_t         Muon_pt[27];   //[nMuon]
   Float_t         Muon_ptErr[27];   //[nMuon]
   Float_t         Muon_segmentComp[27];   //[nMuon]
   Float_t         Muon_sip3d[27];   //[nMuon]
   Float_t         Muon_softMva[27];   //[nMuon]
   Float_t         Muon_softMvaRun3[27];   //[nMuon]
   Float_t         Muon_tkRelIso[27];   //[nMuon]
   Float_t         Muon_tunepRelPt[27];   //[nMuon]
   Float_t         Muon_bsConstrainedChi2[27];   //[nMuon]
   Float_t         Muon_bsConstrainedPt[27];   //[nMuon]
   Float_t         Muon_bsConstrainedPtErr[27];   //[nMuon]
   Float_t         Muon_mvaLowPt[27];   //[nMuon]
   Float_t         Muon_promptMVA[27];   //[nMuon]
   Int_t           nPFCandToPFClusterHCAL;
   Int_t           PFCandToPFClusterHCAL_PFCandIdx[528];   //[nPFCandToPFClusterHCAL]
   Int_t           PFCandToPFClusterHCAL_PFClusterHCALIdx[528];   //[nPFCandToPFClusterHCAL]
   Int_t           nPFClusterECAL;
   Int_t           PFClusterECAL_layer[261];   //[nPFClusterECAL]
   Int_t           PFClusterECAL_nhits[261];   //[nPFClusterECAL]
   UInt_t          PFClusterECAL_seedhit_detId[261];   //[nPFClusterECAL]
   Float_t         PFClusterECAL_pt[261];   //[nPFClusterECAL]
   Float_t         PFClusterECAL_energy[261];   //[nPFClusterECAL]
   Float_t         PFClusterECAL_correctedEnergy[261];   //[nPFClusterECAL]
   Float_t         PFClusterECAL_eta[261];   //[nPFClusterECAL]
   Float_t         PFClusterECAL_phi[261];   //[nPFClusterECAL]
   Int_t           nPFClusterHCAL;
   Int_t           PFClusterHCAL_layer[238];   //[nPFClusterHCAL]
   Int_t           PFClusterHCAL_nhits[238];   //[nPFClusterHCAL]
   UInt_t          PFClusterHCAL_seedhit_detId[238];   //[nPFClusterHCAL]
   Float_t         PFClusterHCAL_pt[238];   //[nPFClusterHCAL]
   Float_t         PFClusterHCAL_energy[238];   //[nPFClusterHCAL]
   Float_t         PFClusterHCAL_eta[238];   //[nPFClusterHCAL]
   Float_t         PFClusterHCAL_phi[238];   //[nPFClusterHCAL]
   Int_t           nPFClusterHCALToPFRecHitHBHE;
   Int_t           PFClusterHCALToPFRecHitHBHE_PFClusterHCALIdx[15748];   //[nPFClusterHCALToPFRecHitHBHE]
   Int_t           PFClusterHCALToPFRecHitHBHE_PFRecHitHBHEIdx[15748];   //[nPFClusterHCALToPFRecHitHBHE]
   Float_t         PFClusterHCALToPFRecHitHBHE_fraction[15748];   //[nPFClusterHCALToPFRecHitHBHE]
   Int_t           nPFClusterPS;
   Int_t           PFClusterPS_layer[388];   //[nPFClusterPS]
   Int_t           PFClusterPS_nhits[388];   //[nPFClusterPS]
   UInt_t          PFClusterPS_seedhit_detId[388];   //[nPFClusterPS]
   Float_t         PFClusterPS_pt[388];   //[nPFClusterPS]
   Float_t         PFClusterPS_energy[388];   //[nPFClusterPS]
   Float_t         PFClusterPS_correctedEnergy[388];   //[nPFClusterPS]
   Float_t         PFClusterPS_eta[388];   //[nPFClusterPS]
   Float_t         PFClusterPS_phi[388];   //[nPFClusterPS]
   Int_t           nPFRecHitHBHE;
   Int_t           PFRecHitHBHE_ieta[3004];   //[nPFRecHitHBHE]
   Int_t           PFRecHitHBHE_iphi[3004];   //[nPFRecHitHBHE]
   Int_t           PFRecHitHBHE_depth[3004];   //[nPFRecHitHBHE]
   UInt_t          PFRecHitHBHE_detId[3004];   //[nPFRecHitHBHE]
   Float_t         PFRecHitHBHE_energy[3004];   //[nPFRecHitHBHE]
  static const int nMaxPFCand = 10000;
   Int_t           nPFCand;
   Bool_t          PFCand_isIsolatedChargedHadron[nMaxPFCand];   //[nPFCand]
   Bool_t          PFCand_trkHighPurity[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_charge[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_fromPV0[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_lostInnerHits[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_lostOuterHits[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_nHits[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_nPixelHits[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_pdgId[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_pvAssocQuality[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_trkAlgo[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_trkQuality[nMaxPFCand];   //[nPFCand]
   Int_t           PFCand_vertexRef[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_d0[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_d0Err[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_dz[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_dzErr[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_energy[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_eta[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_mass[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_phi[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_pt[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_puppiWeight[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_puppiWeightNoLep[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_trkChi2[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_trkEta[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_trkP[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_trkPhi[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_trkPt[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_vtxChi2[nMaxPFCand];   //[nPFCand]
   Float_t         PFMET_covXX;
   Float_t         PFMET_covXY;
   Float_t         PFMET_covYY;
   Float_t         PFMET_phi;
   Float_t         PFMET_phiUnclusteredDown;
   Float_t         PFMET_phiUnclusteredUp;
   Float_t         PFMET_pt;
   Float_t         PFMET_ptUnclusteredDown;
   Float_t         PFMET_ptUnclusteredUp;
   Float_t         PFMET_significance;
   Float_t         PFMET_sumEt;
   Float_t         PFMET_sumPtUnclustered;
   Int_t           nPhoton;
   UChar_t         Photon_cutBased[140];   //[nPhoton]
   Bool_t          Photon_electronVeto[140];   //[nPhoton]
   Bool_t          Photon_hasConversionTracks[140];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[140];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[140];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[140];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[140];   //[nPhoton]
   Bool_t          Photon_pixelSeed[140];   //[nPhoton]
   UChar_t         Photon_seedGain[140];   //[nPhoton]
   Short_t         Photon_electronIdx[140];   //[nPhoton]
   Short_t         Photon_jetIdx[140];   //[nPhoton]
   Short_t         Photon_seediEtaOriX[140];   //[nPhoton]
   Short_t         Photon_seediPhiOriY[140];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[140];   //[nPhoton]
   Float_t         Photon_ecalPFClusterIso[140];   //[nPhoton]
   Float_t         Photon_energyErr[140];   //[nPhoton]
   Float_t         Photon_energyRaw[140];   //[nPhoton]
   Float_t         Photon_esEffSigmaRR[140];   //[nPhoton]
   Float_t         Photon_esEnergyOverRawE[140];   //[nPhoton]
   Float_t         Photon_eta[140];   //[nPhoton]
   Float_t         Photon_etaWidth[140];   //[nPhoton]
   Float_t         Photon_haloTaggerMVAVal[140];   //[nPhoton]
   Float_t         Photon_hcalPFClusterIso[140];   //[nPhoton]
   Float_t         Photon_hoe[140];   //[nPhoton]
   Float_t         Photon_hoe_PUcorr[140];   //[nPhoton]
   Float_t         Photon_hoe_Tower[140];   //[nPhoton]
   Float_t         Photon_mvaID[140];   //[nPhoton]
   Float_t         Photon_pfChargedIso[140];   //[nPhoton]
   Float_t         Photon_pfChargedIsoPFPV[140];   //[nPhoton]
   Float_t         Photon_pfChargedIsoWorstVtx[140];   //[nPhoton]
   Float_t         Photon_pfPhoIso03[140];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all_quadratic[140];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg_quadratic[140];   //[nPhoton]
   Float_t         Photon_phi[140];   //[nPhoton]
   Float_t         Photon_phiWidth[140];   //[nPhoton]
   Float_t         Photon_pt[140];   //[nPhoton]
   Float_t         Photon_r9[140];   //[nPhoton]
   Float_t         Photon_s4[140];   //[nPhoton]
   Float_t         Photon_sieie[140];   //[nPhoton]
   Float_t         Photon_sieip[140];   //[nPhoton]
   Float_t         Photon_sipip[140];   //[nPhoton]
   Float_t         Photon_superclusterEta[140];   //[nPhoton]
   Float_t         Photon_trkSumPtHollowConeDR03[140];   //[nPhoton]
   Float_t         Photon_trkSumPtSolidConeDR04[140];   //[nPhoton]
   Float_t         Photon_x_calo[140];   //[nPhoton]
   Float_t         Photon_y_calo[140];   //[nPhoton]
   Float_t         Photon_z_calo[140];   //[nPhoton]
   Int_t           nPPSLocalTrack;
   Int_t           PPSLocalTrack_multiRPProtonIdx[25];   //[nPPSLocalTrack]
   Int_t           PPSLocalTrack_singleRPProtonIdx[25];   //[nPPSLocalTrack]
   Int_t           PPSLocalTrack_decRPId[25];   //[nPPSLocalTrack]
   Int_t           PPSLocalTrack_rpType[25];   //[nPPSLocalTrack]
   Float_t         PPSLocalTrack_x[25];   //[nPPSLocalTrack]
   Float_t         PPSLocalTrack_y[25];   //[nPPSLocalTrack]
   Float_t         PPSLocalTrack_time[25];   //[nPPSLocalTrack]
   Float_t         PPSLocalTrack_timeUnc[25];   //[nPPSLocalTrack]
   Float_t         PuppiMET_covXX;
   Float_t         PuppiMET_covXY;
   Float_t         PuppiMET_covYY;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_phiUnclusteredDown;
   Float_t         PuppiMET_phiUnclusteredUp;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_ptUnclusteredDown;
   Float_t         PuppiMET_ptUnclusteredUp;
   Float_t         PuppiMET_significance;
   Float_t         PuppiMET_sumEt;
   Float_t         PuppiMET_sumPtUnclustered;
   Int_t           nPVBS;
   Float_t         PVBS_chi2[1];   //[nPVBS]
   Float_t         PVBS_cov00[1];   //[nPVBS]
   Float_t         PVBS_cov10[1];   //[nPVBS]
   Float_t         PVBS_cov11[1];   //[nPVBS]
   Float_t         PVBS_cov20[1];   //[nPVBS]
   Float_t         PVBS_cov21[1];   //[nPVBS]
   Float_t         PVBS_cov22[1];   //[nPVBS]
   Float_t         PVBS_x[1];   //[nPVBS]
   Float_t         PVBS_y[1];   //[nPVBS]
   Float_t         PVBS_z[1];   //[nPVBS]
   Float_t         RawPFMET_phi;
   Float_t         RawPFMET_pt;
   Float_t         RawPFMET_sumEt;
   Float_t         RawPuppiMET_phi;
   Float_t         RawPuppiMET_pt;
   Float_t         RawPuppiMET_sumEt;
   Float_t         Rho_fixedGridRhoAll;
   Float_t         Rho_fixedGridRhoFastjetAll;
   Float_t         Rho_fixedGridRhoFastjetCentral;
   Float_t         Rho_fixedGridRhoFastjetCentralCalo;
   Float_t         Rho_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         Rho_fixedGridRhoFastjetCentralNeutral;
   Int_t           nSoftActivityJet;
   Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;
   Int_t           nProton_singleRP;
   Short_t         Proton_singleRP_decRPId[25];   //[nProton_singleRP]
   Float_t         Proton_singleRP_thetaY[25];   //[nProton_singleRP]
   Float_t         Proton_singleRP_xi[25];   //[nProton_singleRP]
   Int_t           nSubJet;
   Float_t         SubJet_area[14];   //[nSubJet]
   Float_t         SubJet_btagDeepB[14];   //[nSubJet]
   Float_t         SubJet_eta[14];   //[nSubJet]
   Float_t         SubJet_mass[14];   //[nSubJet]
   Float_t         SubJet_n2b1[14];   //[nSubJet]
   Float_t         SubJet_n3b1[14];   //[nSubJet]
   Float_t         SubJet_phi[14];   //[nSubJet]
   Float_t         SubJet_pt[14];   //[nSubJet]
   Float_t         SubJet_rawFactor[14];   //[nSubJet]
   Float_t         SubJet_tau1[14];   //[nSubJet]
   Float_t         SubJet_tau2[14];   //[nSubJet]
   Float_t         SubJet_tau3[14];   //[nSubJet]
   Float_t         SubJet_tau4[14];   //[nSubJet]
   Int_t           nTauProd;
   Short_t         TauProd_tauIdx[24];   //[nTauProd]
   Int_t           TauProd_pdgId[24];   //[nTauProd]
   Float_t         TauProd_eta[24];   //[nTauProd]
   Float_t         TauProd_phi[24];   //[nTauProd]
   Float_t         TauProd_pt[24];   //[nTauProd]
   Int_t           nTau;
   UChar_t         Tau_decayMode[6];   //[nTau]
   Bool_t          Tau_idAntiEleDeadECal[6];   //[nTau]
   UChar_t         Tau_idAntiMu[6];   //[nTau]
   Bool_t          Tau_idDecayModeNewDMs[6];   //[nTau]
   Bool_t          Tau_idDecayModeOldDMs[6];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSe[6];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSjet[6];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSmu[6];   //[nTau]
   UChar_t         Tau_idDeepTau2018v2p5VSe[6];   //[nTau]
   UChar_t         Tau_idDeepTau2018v2p5VSjet[6];   //[nTau]
   UChar_t         Tau_idDeepTau2018v2p5VSmu[6];   //[nTau]
   UChar_t         Tau_nSVs[6];   //[nTau]
   Short_t         Tau_charge[6];   //[nTau]
   Short_t         Tau_decayModePNet[6];   //[nTau]
   Short_t         Tau_decayModeUParT[6];   //[nTau]
   Short_t         Tau_eleIdx[6];   //[nTau]
   Short_t         Tau_jetIdx[6];   //[nTau]
   Short_t         Tau_muIdx[6];   //[nTau]
   Short_t         Tau_svIdx1[6];   //[nTau]
   Short_t         Tau_svIdx2[6];   //[nTau]
   Float_t         Tau_chargedIso[6];   //[nTau]
   Float_t         Tau_dxy[6];   //[nTau]
   Float_t         Tau_dz[6];   //[nTau]
   Float_t         Tau_eta[6];   //[nTau]
   Float_t         Tau_leadTkDeltaEta[6];   //[nTau]
   Float_t         Tau_leadTkDeltaPhi[6];   //[nTau]
   Float_t         Tau_leadTkPtOverTauPt[6];   //[nTau]
   Float_t         Tau_mass[6];   //[nTau]
   Float_t         Tau_neutralIso[6];   //[nTau]
   Float_t         Tau_phi[6];   //[nTau]
   Float_t         Tau_photonsOutsideSignalCone[6];   //[nTau]
   Float_t         Tau_probDM0PNet[6];   //[nTau]
   Float_t         Tau_probDM0UParT[6];   //[nTau]
   Float_t         Tau_probDM10PNet[6];   //[nTau]
   Float_t         Tau_probDM10UParT[6];   //[nTau]
   Float_t         Tau_probDM11PNet[6];   //[nTau]
   Float_t         Tau_probDM11UParT[6];   //[nTau]
   Float_t         Tau_probDM1PNet[6];   //[nTau]
   Float_t         Tau_probDM1UParT[6];   //[nTau]
   Float_t         Tau_probDM2PNet[6];   //[nTau]
   Float_t         Tau_probDM2UParT[6];   //[nTau]
   Float_t         Tau_pt[6];   //[nTau]
   Float_t         Tau_ptCorrPNet[6];   //[nTau]
   Float_t         Tau_ptCorrUParT[6];   //[nTau]
   Float_t         Tau_puCorr[6];   //[nTau]
   Float_t         Tau_qConfPNet[6];   //[nTau]
   Float_t         Tau_qConfUParT[6];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSe[6];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSjet[6];   //[nTau]
   Float_t         Tau_rawDeepTau2017v2p1VSmu[6];   //[nTau]
   Float_t         Tau_rawDeepTau2018v2p5VSe[6];   //[nTau]
   Float_t         Tau_rawDeepTau2018v2p5VSjet[6];   //[nTau]
   Float_t         Tau_rawDeepTau2018v2p5VSmu[6];   //[nTau]
   Float_t         Tau_rawIso[6];   //[nTau]
   Float_t         Tau_rawIsodR03[6];   //[nTau]
   Float_t         Tau_rawPNetVSe[6];   //[nTau]
   Float_t         Tau_rawPNetVSjet[6];   //[nTau]
   Float_t         Tau_rawPNetVSmu[6];   //[nTau]
   Float_t         Tau_rawUParTVSe[6];   //[nTau]
   Float_t         Tau_rawUParTVSjet[6];   //[nTau]
   Float_t         Tau_rawUParTVSmu[6];   //[nTau]
   Int_t           nTrigObj;
   Short_t         TrigObj_l1charge[110];   //[nTrigObj]
   UShort_t        TrigObj_id[110];   //[nTrigObj]
   Int_t           TrigObj_l1iso[110];   //[nTrigObj]
   Int_t           TrigObj_filterBits[110];   //[nTrigObj]
   Float_t         TrigObj_pt[110];   //[nTrigObj]
   Float_t         TrigObj_eta[110];   //[nTrigObj]
   Float_t         TrigObj_phi[110];   //[nTrigObj]
   Float_t         TrigObj_l1pt[110];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[110];   //[nTrigObj]
   Float_t         TrigObj_l2pt[110];   //[nTrigObj]
   Float_t         TrkMET_phi;
   Float_t         TrkMET_pt;
   Float_t         TrkMET_sumEt;
   Int_t           nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         OtherPV_score[3];   //[nOtherPV]
   UChar_t         PV_npvs;
   UChar_t         PV_npvsGood;
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Float_t         PV_sumpt2;
   Float_t         PV_sumpx;
   Float_t         PV_sumpy;
   Int_t           nSV;
   Short_t         SV_charge[24];   //[nSV]
   Float_t         SV_dlen[24];   //[nSV]
   Float_t         SV_dlenSig[24];   //[nSV]
   Float_t         SV_dxy[24];   //[nSV]
   Float_t         SV_dxySig[24];   //[nSV]
   Float_t         SV_pAngle[24];   //[nSV]
   Float_t         Electron_IPx[11];   //[nElectron]
   Float_t         Electron_IPy[11];   //[nElectron]
   Float_t         Electron_IPz[11];   //[nElectron]
   Float_t         Electron_ipLengthSig[11];   //[nElectron]
   Float_t         Muon_IPx[27];   //[nMuon]
   Float_t         Muon_IPy[27];   //[nMuon]
   Float_t         Muon_IPz[27];   //[nMuon]
   Float_t         Muon_ipLengthSig[27];   //[nMuon]
   Int_t           PFCand_fromPVvertexRef[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_blockIdx[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFTrack[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFClusterHCAL[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFClusterECAL[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFClusterPS[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFTrackInBlock[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFClusterHCALInBlock[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFClusterECALInBlock[nMaxPFCand];   //[nPFCand]
   UInt_t          PFCand_nPFClusterPSInBlock[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalDepthEnergyFraction1[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalDepthEnergyFraction2[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalDepthEnergyFraction3[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalDepthEnergyFraction4[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalDepthEnergyFraction5[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalDepthEnergyFraction6[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalDepthEnergyFraction7[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_ecalEnergy[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_rawEcalEnergy[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hcalEnergy[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_rawHcalEnergy[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_hoEnergy[nMaxPFCand];   //[nPFCand]
   Float_t         PFCand_rawHoEnergy[nMaxPFCand];   //[nPFCand]
   UChar_t         SV_ntracks[24];   //[nSV]
   Float_t         SV_chi2[24];   //[nSV]
   Float_t         SV_eta[24];   //[nSV]
   Float_t         SV_mass[24];   //[nSV]
   Float_t         SV_ndof[24];   //[nSV]
   Float_t         SV_phi[24];   //[nSV]
   Float_t         SV_pt[24];   //[nSV]
   Float_t         SV_x[24];   //[nSV]
   Float_t         SV_y[24];   //[nSV]
   Float_t         SV_z[24];   //[nSV]
   Bool_t          Tau_hasRefitSV[6];   //[nTau]
   Float_t         Tau_IPx[6];   //[nTau]
   Float_t         Tau_IPy[6];   //[nTau]
   Float_t         Tau_IPz[6];   //[nTau]
   Float_t         Tau_ipLengthSig[6];   //[nTau]
   Float_t         Tau_refitSVchi2[6];   //[nTau]
   Float_t         Tau_refitSVcov00[6];   //[nTau]
   Float_t         Tau_refitSVcov10[6];   //[nTau]
   Float_t         Tau_refitSVcov11[6];   //[nTau]
   Float_t         Tau_refitSVcov20[6];   //[nTau]
   Float_t         Tau_refitSVcov21[6];   //[nTau]
   Float_t         Tau_refitSVcov22[6];   //[nTau]
   Float_t         Tau_refitSVx[6];   //[nTau]
   Float_t         Tau_refitSVy[6];   //[nTau]
   Float_t         Tau_refitSVz[6];   //[nTau]
   Bool_t          L1Reco_step;
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadPFMuonDzFilter;
   Bool_t          Flag_hfNoisyHitsFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   Bool_t          L1_AXO_Loose;
   Bool_t          L1_AXO_Nominal;
   Bool_t          L1_AXO_Tight;
   Bool_t          L1_AXO_VLoose;
   Bool_t          L1_AXO_VTight;
   Bool_t          L1_AlwaysTrue;
   Bool_t          L1_BPTX_AND_Ref1_VME;
   Bool_t          L1_BPTX_AND_Ref3_VME;
   Bool_t          L1_BPTX_AND_Ref4_VME;
   Bool_t          L1_BPTX_BeamGas_B1_VME;
   Bool_t          L1_BPTX_BeamGas_B2_VME;
   Bool_t          L1_BPTX_BeamGas_Ref1_VME;
   Bool_t          L1_BPTX_BeamGas_Ref2_VME;
   Bool_t          L1_BPTX_NotOR_VME;
   Bool_t          L1_BPTX_OR_Ref3_VME;
   Bool_t          L1_BPTX_OR_Ref4_VME;
   Bool_t          L1_BPTX_RefAND_VME;
   Bool_t          L1_BptxMinus;
   Bool_t          L1_BptxOR;
   Bool_t          L1_BptxPlus;
   Bool_t          L1_BptxXOR;
   Bool_t          L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          L1_CICADA_Loose;
   Bool_t          L1_CICADA_Medium;
   Bool_t          L1_CICADA_Tight;
   Bool_t          L1_CICADA_VLoose;
   Bool_t          L1_CICADA_VTight;
   Bool_t          L1_DoubleEG11_er1p2_dR_Max0p6;
   Bool_t          L1_DoubleEG8er2p5_HTT280er;
   Bool_t          L1_DoubleEG8er2p5_HTT300er;
   Bool_t          L1_DoubleEG8er2p5_HTT320er;
   Bool_t          L1_DoubleEG_15_10_er2p5;
   Bool_t          L1_DoubleEG_20_10_er2p5;
   Bool_t          L1_DoubleEG_22_10_er2p5;
   Bool_t          L1_DoubleEG_25_12_er2p5;
   Bool_t          L1_DoubleEG_25_14_er2p5;
   Bool_t          L1_DoubleEG_27_14_er2p5;
   Bool_t          L1_DoubleEG_LooseIso16_LooseIso12_er1p5;
   Bool_t          L1_DoubleEG_LooseIso18_LooseIso12_er1p5;
   Bool_t          L1_DoubleEG_LooseIso20_LooseIso12_er1p5;
   Bool_t          L1_DoubleEG_LooseIso22_12_er2p5;
   Bool_t          L1_DoubleEG_LooseIso22_LooseIso12_er1p5;
   Bool_t          L1_DoubleEG_LooseIso25_12_er2p5;
   Bool_t          L1_DoubleEG_LooseIso25_LooseIso12_er1p5;
   Bool_t          L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5;
   Bool_t          L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5;
   Bool_t          L1_DoubleIsoTau28er2p1;
   Bool_t          L1_DoubleIsoTau30er2p1;
   Bool_t          L1_DoubleIsoTau32er2p1;
   Bool_t          L1_DoubleIsoTau32er2p1_Mass_Max80;
   Bool_t          L1_DoubleIsoTau34er2p1;
   Bool_t          L1_DoubleIsoTau35er2p1;
   Bool_t          L1_DoubleIsoTau36er2p1;
   Bool_t          L1_DoubleJet100er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet100er2p5;
   Bool_t          L1_DoubleJet112er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet120er2p5;
   Bool_t          L1_DoubleJet120er2p5_Mu3_dR_Max0p8;
   Bool_t          L1_DoubleJet150er2p5;
   Bool_t          L1_DoubleJet16er2p5_Mu3_dR_Max0p4;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5;
   Bool_t          L1_DoubleJet35er2p5_Mu3_dR_Max0p4;
   Bool_t          L1_DoubleJet40er2p5;
   Bool_t          L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5;
   Bool_t          L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2;
   Bool_t          L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5;
   Bool_t          L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2;
   Bool_t          L1_DoubleJet60er2p5_Mu3_dR_Max0p4;
   Bool_t          L1_DoubleJet80er2p5_Mu3_dR_Max0p4;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min800;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min850;
   Bool_t          L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50;
   Bool_t          L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50;
   Bool_t          L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65;
   Bool_t          L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65;
   Bool_t          L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ;
   Bool_t          L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ;
   Bool_t          L1_DoubleLLPJet40;
   Bool_t          L1_DoubleLooseIsoEG22er2p1;
   Bool_t          L1_DoubleLooseIsoEG24er2p1;
   Bool_t          L1_DoubleMu0;
   Bool_t          L1_DoubleMu0_Mass_Min1;
   Bool_t          L1_DoubleMu0_OQ;
   Bool_t          L1_DoubleMu0_SQ;
   Bool_t          L1_DoubleMu0_SQ_OS;
   Bool_t          L1_DoubleMu0_Upt15_Upt7;
   Bool_t          L1_DoubleMu0_Upt5_Upt5;
   Bool_t          L1_DoubleMu0_Upt6_IP_Min1_Upt4;
   Bool_t          L1_DoubleMu0_Upt6_SQ_er2p0;
   Bool_t          L1_DoubleMu0_Upt7_SQ_er2p0;
   Bool_t          L1_DoubleMu0_Upt8_SQ_er2p0;
   Bool_t          L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8;
   Bool_t          L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p5_SQ;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p5_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2;
   Bool_t          L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5;
   Bool_t          L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6;
   Bool_t          L1_DoubleMu0er2p0_SQ_dEta_Max1p5;
   Bool_t          L1_DoubleMu0er2p0_SQ_dEta_Max1p6;
   Bool_t          L1_DoubleMu18er2p1_SQ;
   Bool_t          L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20;
   Bool_t          L1_DoubleMu3_SQ_ETMHF30_HTT60er;
   Bool_t          L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t          L1_DoubleMu3_SQ_ETMHF40_HTT60er;
   Bool_t          L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_HTT60er;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5;
   Bool_t          L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5;
   Bool_t          L1_DoubleMu3_SQ_HTT220er;
   Bool_t          L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6;
   Bool_t          L1_DoubleMu4_SQ_EG9er2p5;
   Bool_t          L1_DoubleMu4_SQ_OS;
   Bool_t          L1_DoubleMu4_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6;
   Bool_t          L1_DoubleMu4p5_SQ_OS;
   Bool_t          L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7;
   Bool_t          L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20;
   Bool_t          L1_DoubleMu5_SQ_EG9er2p5;
   Bool_t          L1_DoubleMu5_SQ_OS_dR_Max1p6;
   Bool_t          L1_DoubleMu6_Upt6_SQ_er2p0;
   Bool_t          L1_DoubleMu7_Upt7_SQ_er2p0;
   Bool_t          L1_DoubleMu8_SQ;
   Bool_t          L1_DoubleMu8_Upt8_SQ_er2p0;
   Bool_t          L1_DoubleMu9_SQ;
   Bool_t          L1_DoubleMu_12_5;
   Bool_t          L1_DoubleMu_15_5_SQ;
   Bool_t          L1_DoubleMu_15_7;
   Bool_t          L1_DoubleMu_15_7_Mass_Min1;
   Bool_t          L1_DoubleMu_15_7_SQ;
   Bool_t          L1_DoubleTau70er2p1;
   Bool_t          L1_ETM120;
   Bool_t          L1_ETM150;
   Bool_t          L1_ETMHF100;
   Bool_t          L1_ETMHF100_HTT60er;
   Bool_t          L1_ETMHF110;
   Bool_t          L1_ETMHF110_HTT60er;
   Bool_t          L1_ETMHF120;
   Bool_t          L1_ETMHF120_HTT60er;
   Bool_t          L1_ETMHF130;
   Bool_t          L1_ETMHF130_HTT60er;
   Bool_t          L1_ETMHF140;
   Bool_t          L1_ETMHF150;
   Bool_t          L1_ETMHF70;
   Bool_t          L1_ETMHF70_HTT60er;
   Bool_t          L1_ETMHF80;
   Bool_t          L1_ETMHF80_HTT60er;
   Bool_t          L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1;
   Bool_t          L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6;
   Bool_t          L1_ETMHF90;
   Bool_t          L1_ETMHF90_HTT60er;
   Bool_t          L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1;
   Bool_t          L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6;
   Bool_t          L1_ETT2000;
   Bool_t          L1_FirstBunchAfterTrain;
   Bool_t          L1_FirstBunchBeforeTrain;
   Bool_t          L1_FirstBunchInTrain;
   Bool_t          L1_FirstCollisionInOrbit;
   Bool_t          L1_FirstCollisionInTrain;
   Bool_t          L1_HCAL_LaserMon_Trig;
   Bool_t          L1_HCAL_LaserMon_Veto;
   Bool_t          L1_HTMHF100;
   Bool_t          L1_HTMHF120;
   Bool_t          L1_HTMHF125;
   Bool_t          L1_HTMHF130;
   Bool_t          L1_HTMHF150;
   Bool_t          L1_HTT120_SingleLLPJet40;
   Bool_t          L1_HTT120er;
   Bool_t          L1_HTT160_SingleLLPJet50;
   Bool_t          L1_HTT160er;
   Bool_t          L1_HTT200_SingleLLPJet60;
   Bool_t          L1_HTT200er;
   Bool_t          L1_HTT240_SingleLLPJet70;
   Bool_t          L1_HTT255er;
   Bool_t          L1_HTT280er;
   Bool_t          L1_HTT280er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT320er;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p5;
   Bool_t          L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;
   Bool_t          L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;
   Bool_t          L1_HTT360er;
   Bool_t          L1_HTT400er;
   Bool_t          L1_HTT450er;
   Bool_t          L1_IsoEG32er2p5_Mt40;
   Bool_t          L1_IsoTau52er2p1_QuadJet36er2p5;
   Bool_t          L1_IsolatedBunch;
   Bool_t          L1_LastBunchInTrain;
   Bool_t          L1_LastCollisionInTrain;
   Bool_t          L1_LooseIsoEG14er2p5_HTT200er;
   Bool_t          L1_LooseIsoEG16er2p5_HTT200er;
   Bool_t          L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG26er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          L1_LooseIsoEG28er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;
   Bool_t          L1_MinimumBiasHF0;
   Bool_t          L1_MinimumBiasHF0_AND_BptxAND;
   Bool_t          L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;
   Bool_t          L1_Mu12_HTT150er;
   Bool_t          L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;
   Bool_t          L1_Mu14_HTT150er;
   Bool_t          L1_Mu18er2p1_Tau24er2p1;
   Bool_t          L1_Mu18er2p1_Tau26er2p1;
   Bool_t          L1_Mu18er2p1_Tau26er2p1_Jet55;
   Bool_t          L1_Mu18er2p1_Tau26er2p1_Jet70;
   Bool_t          L1_Mu20_EG10er2p5;
   Bool_t          L1_Mu22er2p1_IsoTau30er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau32er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau34er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau40er2p1;
   Bool_t          L1_Mu22er2p1_Tau70er2p1;
   Bool_t          L1_Mu3_Jet120er2p5_dR_Max0p4;
   Bool_t          L1_Mu3_Jet16er2p5_dR_Max0p4;
   Bool_t          L1_Mu3_Jet30er2p5;
   Bool_t          L1_Mu3_Jet60er2p5_dR_Max0p4;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF30;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF40;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF50;
   Bool_t          L1_Mu5_EG23er2p5;
   Bool_t          L1_Mu5_LooseIsoEG20er2p5;
   Bool_t          L1_Mu6_DoubleEG10er2p5;
   Bool_t          L1_Mu6_DoubleEG12er2p5;
   Bool_t          L1_Mu6_DoubleEG15er2p5;
   Bool_t          L1_Mu6_DoubleEG17er2p5;
   Bool_t          L1_Mu6_HTT240er;
   Bool_t          L1_Mu6_HTT250er;
   Bool_t          L1_Mu7_EG20er2p5;
   Bool_t          L1_Mu7_EG23er2p5;
   Bool_t          L1_Mu7_LooseIsoEG20er2p5;
   Bool_t          L1_Mu7_LooseIsoEG23er2p5;
   Bool_t          L1_NotBptxOR;
   Bool_t          L1_QuadJet60er2p5;
   Bool_t          L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0;
   Bool_t          L1_QuadMu0;
   Bool_t          L1_QuadMu0_OQ;
   Bool_t          L1_QuadMu0_SQ;
   Bool_t          L1_SecondBunchInTrain;
   Bool_t          L1_SecondLastBunchInTrain;
   Bool_t          L1_SingleEG10er2p5;
   Bool_t          L1_SingleEG15er2p5;
   Bool_t          L1_SingleEG26er2p5;
   Bool_t          L1_SingleEG28_FWD2p5;
   Bool_t          L1_SingleEG28er1p5;
   Bool_t          L1_SingleEG28er2p1;
   Bool_t          L1_SingleEG28er2p5;
   Bool_t          L1_SingleEG34er2p5;
   Bool_t          L1_SingleEG36er2p5;
   Bool_t          L1_SingleEG38er2p5;
   Bool_t          L1_SingleEG40er2p5;
   Bool_t          L1_SingleEG42er2p5;
   Bool_t          L1_SingleEG45er2p5;
   Bool_t          L1_SingleEG50;
   Bool_t          L1_SingleEG60;
   Bool_t          L1_SingleEG8er2p5;
   Bool_t          L1_SingleIsoEG24er2p1;
   Bool_t          L1_SingleIsoEG26er2p1;
   Bool_t          L1_SingleIsoEG26er2p5;
   Bool_t          L1_SingleIsoEG28_FWD2p5;
   Bool_t          L1_SingleIsoEG28er1p5;
   Bool_t          L1_SingleIsoEG28er2p1;
   Bool_t          L1_SingleIsoEG28er2p5;
   Bool_t          L1_SingleIsoEG30er2p1;
   Bool_t          L1_SingleIsoEG30er2p5;
   Bool_t          L1_SingleIsoEG32er2p1;
   Bool_t          L1_SingleIsoEG32er2p5;
   Bool_t          L1_SingleIsoEG34er2p5;
   Bool_t          L1_SingleJet10erHE;
   Bool_t          L1_SingleJet120;
   Bool_t          L1_SingleJet120_FWD2p5;
   Bool_t          L1_SingleJet120_FWD3p0;
   Bool_t          L1_SingleJet120er1p3;
   Bool_t          L1_SingleJet120er2p5;
   Bool_t          L1_SingleJet12erHE;
   Bool_t          L1_SingleJet140er2p5;
   Bool_t          L1_SingleJet160er2p5;
   Bool_t          L1_SingleJet180;
   Bool_t          L1_SingleJet180er2p5;
   Bool_t          L1_SingleJet200;
   Bool_t          L1_SingleJet20er2p5_NotBptxOR;
   Bool_t          L1_SingleJet20er2p5_NotBptxOR_3BX;
   Bool_t          L1_SingleJet35;
   Bool_t          L1_SingleJet35_FWD2p5;
   Bool_t          L1_SingleJet35_FWD3p0;
   Bool_t          L1_SingleJet35er1p3;
   Bool_t          L1_SingleJet35er2p5;
   Bool_t          L1_SingleJet43er2p5_NotBptxOR_3BX;
   Bool_t          L1_SingleJet46er2p5_NotBptxOR_3BX;
   Bool_t          L1_SingleJet60;
   Bool_t          L1_SingleJet60_FWD2p5;
   Bool_t          L1_SingleJet8erHE;
   Bool_t          L1_SingleJet90;
   Bool_t          L1_SingleJet90_FWD2p5;
   Bool_t          L1_SingleLooseIsoEG26er1p5;
   Bool_t          L1_SingleLooseIsoEG26er2p5;
   Bool_t          L1_SingleLooseIsoEG28_FWD2p5;
   Bool_t          L1_SingleLooseIsoEG28er1p5;
   Bool_t          L1_SingleLooseIsoEG28er2p1;
   Bool_t          L1_SingleLooseIsoEG28er2p5;
   Bool_t          L1_SingleLooseIsoEG30er1p5;
   Bool_t          L1_SingleLooseIsoEG30er2p5;
   Bool_t          L1_SingleMu0_BMTF;
   Bool_t          L1_SingleMu0_DQ;
   Bool_t          L1_SingleMu0_EMTF;
   Bool_t          L1_SingleMu0_OMTF;
   Bool_t          L1_SingleMu0_SQ13_BMTF;
   Bool_t          L1_SingleMu0_SQ14_BMTF;
   Bool_t          L1_SingleMu0_SQ15_BMTF;
   Bool_t          L1_SingleMu0_Upt10;
   Bool_t          L1_SingleMu0_Upt10_BMTF;
   Bool_t          L1_SingleMu0_Upt10_EMTF;
   Bool_t          L1_SingleMu0_Upt10_OMTF;
   Bool_t          L1_SingleMu0_Upt10_SQ14_BMTF;
   Bool_t          L1_SingleMu0_Upt15_SQ14_BMTF;
   Bool_t          L1_SingleMu0_Upt20_SQ14_BMTF;
   Bool_t          L1_SingleMu0_Upt25_SQ14_BMTF;
   Bool_t          L1_SingleMu10_SQ14_BMTF;
   Bool_t          L1_SingleMu11_SQ14_BMTF;
   Bool_t          L1_SingleMu12_DQ_BMTF;
   Bool_t          L1_SingleMu12_DQ_EMTF;
   Bool_t          L1_SingleMu12_DQ_OMTF;
   Bool_t          L1_SingleMu15_DQ;
   Bool_t          L1_SingleMu18;
   Bool_t          L1_SingleMu20;
   Bool_t          L1_SingleMu22;
   Bool_t          L1_SingleMu22_BMTF;
   Bool_t          L1_SingleMu22_BMTF_NEG;
   Bool_t          L1_SingleMu22_BMTF_POS;
   Bool_t          L1_SingleMu22_DQ;
   Bool_t          L1_SingleMu22_EMTF;
   Bool_t          L1_SingleMu22_EMTF_NEG;
   Bool_t          L1_SingleMu22_EMTF_POS;
   Bool_t          L1_SingleMu22_OMTF;
   Bool_t          L1_SingleMu22_OMTF_NEG;
   Bool_t          L1_SingleMu22_OMTF_POS;
   Bool_t          L1_SingleMu22_OQ;
   Bool_t          L1_SingleMu25;
   Bool_t          L1_SingleMu3;
   Bool_t          L1_SingleMu5;
   Bool_t          L1_SingleMu5_SQ14_BMTF;
   Bool_t          L1_SingleMu6_SQ14_BMTF;
   Bool_t          L1_SingleMu7;
   Bool_t          L1_SingleMu7_DQ;
   Bool_t          L1_SingleMu7_SQ14_BMTF;
   Bool_t          L1_SingleMu8_SQ14_BMTF;
   Bool_t          L1_SingleMu9_SQ14_BMTF;
   Bool_t          L1_SingleMuCosmics;
   Bool_t          L1_SingleMuCosmics_BMTF;
   Bool_t          L1_SingleMuCosmics_EMTF;
   Bool_t          L1_SingleMuCosmics_OMTF;
   Bool_t          L1_SingleMuOpen;
   Bool_t          L1_SingleMuOpen_BMTF;
   Bool_t          L1_SingleMuOpen_EMTF;
   Bool_t          L1_SingleMuOpen_NotBptxOR;
   Bool_t          L1_SingleMuOpen_OMTF;
   Bool_t          L1_SingleMuOpen_er1p1_NotBptxOR_3BX;
   Bool_t          L1_SingleMuOpen_er1p4_NotBptxOR_3BX;
   Bool_t          L1_SingleMuShower_Nominal;
   Bool_t          L1_SingleMuShower_Tight;
   Bool_t          L1_SingleTau120er2p1;
   Bool_t          L1_SingleTau130er2p1;
   Bool_t          L1_TOTEM_1;
   Bool_t          L1_TOTEM_2;
   Bool_t          L1_TOTEM_3;
   Bool_t          L1_TOTEM_4;
   Bool_t          L1_TripleEG_18_17_8_er2p5;
   Bool_t          L1_TripleEG_18_18_12_er2p5;
   Bool_t          L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5;
   Bool_t          L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5;
   Bool_t          L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5;
   Bool_t          L1_TripleMu0;
   Bool_t          L1_TripleMu0_OQ;
   Bool_t          L1_TripleMu0_SQ;
   Bool_t          L1_TripleMu3;
   Bool_t          L1_TripleMu3_SQ;
   Bool_t          L1_TripleMu_3SQ_2p5SQ_0;
   Bool_t          L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12;
   Bool_t          L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12;
   Bool_t          L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5_3_3;
   Bool_t          L1_TripleMu_5_3_3_SQ;
   Bool_t          L1_TripleMu_5_3p5_2p5;
   Bool_t          L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_5_3;
   Bool_t          L1_TwoMuShower_Loose;
   Bool_t          L1_UnpairedBunchBptxMinus;
   Bool_t          L1_UnpairedBunchBptxPlus;
   Bool_t          L1_ZeroBias;
   Bool_t          L1_ZeroBias_copy;
   Bool_t          L1_UnprefireableEvent_TriggerRules;
   Bool_t          L1_UnprefireableEvent_FirstBxInTrain;
   Bool_t          L1_FinalOR_BXmin1;
   Bool_t          L1_FinalOR_BXmin2;
   Bool_t          L1Reco_step_pRECO;
   Bool_t          Flag_HBHENoiseFilter_pRECO;
   Bool_t          Flag_HBHENoiseIsoFilter_pRECO;
   Bool_t          Flag_CSCTightHaloFilter_pRECO;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO;
   Bool_t          Flag_CSCTightHalo2015Filter_pRECO;
   Bool_t          Flag_globalTightHalo2016Filter_pRECO;
   Bool_t          Flag_globalSuperTightHalo2016Filter_pRECO;
   Bool_t          Flag_HcalStripHaloFilter_pRECO;
   Bool_t          Flag_hcalLaserEventFilter_pRECO;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter_pRECO;
   Bool_t          Flag_ecalBadCalibFilter_pRECO;
   Bool_t          Flag_goodVertices_pRECO;
   Bool_t          Flag_eeBadScFilter_pRECO;
   Bool_t          Flag_ecalLaserCorrFilter_pRECO;
   Bool_t          Flag_trkPOGFilters_pRECO;
   Bool_t          Flag_chargedHadronTrackResolutionFilter_pRECO;
   Bool_t          Flag_muonBadTrackFilter_pRECO;
   Bool_t          Flag_BadChargedCandidateFilter_pRECO;
   Bool_t          Flag_BadPFMuonFilter_pRECO;
   Bool_t          Flag_BadPFMuonDzFilter_pRECO;
   Bool_t          Flag_hfNoisyHitsFilter_pRECO;
   Bool_t          Flag_BadChargedCandidateSummer16Filter_pRECO;
   Bool_t          Flag_BadPFMuonSummer16Filter_pRECO;
   Bool_t          Flag_trkPOG_manystripclus53X_pRECO;
   Bool_t          Flag_trkPOG_toomanystripclus53X_pRECO;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters_pRECO;
   Bool_t          Flag_METFilters_pRECO;
   Bool_t          L1_AXO_Loose_pRECO;
   Bool_t          L1_AXO_Nominal_pRECO;
   Bool_t          L1_AXO_Tight_pRECO;
   Bool_t          L1_AXO_VLoose_pRECO;
   Bool_t          L1_AXO_VTight_pRECO;
   Bool_t          L1_AlwaysTrue_pRECO;
   Bool_t          L1_BPTX_AND_Ref1_VME_pRECO;
   Bool_t          L1_BPTX_AND_Ref3_VME_pRECO;
   Bool_t          L1_BPTX_AND_Ref4_VME_pRECO;
   Bool_t          L1_BPTX_BeamGas_B1_VME_pRECO;
   Bool_t          L1_BPTX_BeamGas_B2_VME_pRECO;
   Bool_t          L1_BPTX_BeamGas_Ref1_VME_pRECO;
   Bool_t          L1_BPTX_BeamGas_Ref2_VME_pRECO;
   Bool_t          L1_BPTX_NotOR_VME_pRECO;
   Bool_t          L1_BPTX_OR_Ref3_VME_pRECO;
   Bool_t          L1_BPTX_OR_Ref4_VME_pRECO;
   Bool_t          L1_BPTX_RefAND_VME_pRECO;
   Bool_t          L1_BptxMinus_pRECO;
   Bool_t          L1_BptxOR_pRECO;
   Bool_t          L1_BptxPlus_pRECO;
   Bool_t          L1_BptxXOR_pRECO;
   Bool_t          L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_pRECO;
   Bool_t          L1_CICADA_Loose_pRECO;
   Bool_t          L1_CICADA_Medium_pRECO;
   Bool_t          L1_CICADA_Tight_pRECO;
   Bool_t          L1_CICADA_VLoose_pRECO;
   Bool_t          L1_CICADA_VTight_pRECO;
   Bool_t          L1_DoubleEG11_er1p2_dR_Max0p6_pRECO;
   Bool_t          L1_DoubleEG8er2p5_HTT280er_pRECO;
   Bool_t          L1_DoubleEG8er2p5_HTT300er_pRECO;
   Bool_t          L1_DoubleEG8er2p5_HTT320er_pRECO;
   Bool_t          L1_DoubleEG_15_10_er2p5_pRECO;
   Bool_t          L1_DoubleEG_20_10_er2p5_pRECO;
   Bool_t          L1_DoubleEG_22_10_er2p5_pRECO;
   Bool_t          L1_DoubleEG_25_12_er2p5_pRECO;
   Bool_t          L1_DoubleEG_25_14_er2p5_pRECO;
   Bool_t          L1_DoubleEG_27_14_er2p5_pRECO;
   Bool_t          L1_DoubleEG_LooseIso16_LooseIso12_er1p5_pRECO;
   Bool_t          L1_DoubleEG_LooseIso18_LooseIso12_er1p5_pRECO;
   Bool_t          L1_DoubleEG_LooseIso20_LooseIso12_er1p5_pRECO;
   Bool_t          L1_DoubleEG_LooseIso22_12_er2p5_pRECO;
   Bool_t          L1_DoubleEG_LooseIso22_LooseIso12_er1p5_pRECO;
   Bool_t          L1_DoubleEG_LooseIso25_12_er2p5_pRECO;
   Bool_t          L1_DoubleEG_LooseIso25_LooseIso12_er1p5_pRECO;
   Bool_t          L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5_pRECO;
   Bool_t          L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5_pRECO;
   Bool_t          L1_DoubleIsoTau28er2p1_pRECO;
   Bool_t          L1_DoubleIsoTau30er2p1_pRECO;
   Bool_t          L1_DoubleIsoTau32er2p1_pRECO;
   Bool_t          L1_DoubleIsoTau32er2p1_Mass_Max80_pRECO;
   Bool_t          L1_DoubleIsoTau34er2p1_pRECO;
   Bool_t          L1_DoubleIsoTau35er2p1_pRECO;
   Bool_t          L1_DoubleIsoTau36er2p1_pRECO;
   Bool_t          L1_DoubleJet100er2p3_dEta_Max1p6_pRECO;
   Bool_t          L1_DoubleJet100er2p5_pRECO;
   Bool_t          L1_DoubleJet112er2p3_dEta_Max1p6_pRECO;
   Bool_t          L1_DoubleJet120er2p5_pRECO;
   Bool_t          L1_DoubleJet120er2p5_Mu3_dR_Max0p8_pRECO;
   Bool_t          L1_DoubleJet150er2p5_pRECO;
   Bool_t          L1_DoubleJet16er2p5_Mu3_dR_Max0p4_pRECO;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5_pRECO;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5_pRECO;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5_pRECO;
   Bool_t          L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5_pRECO;
   Bool_t          L1_DoubleJet35er2p5_Mu3_dR_Max0p4_pRECO;
   Bool_t          L1_DoubleJet40er2p5_pRECO;
   Bool_t          L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5_pRECO;
   Bool_t          L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO;
   Bool_t          L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5_pRECO;
   Bool_t          L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO;
   Bool_t          L1_DoubleJet60er2p5_Mu3_dR_Max0p4_pRECO;
   Bool_t          L1_DoubleJet80er2p5_Mu3_dR_Max0p4_pRECO;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min800_pRECO;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min850_pRECO;
   Bool_t          L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50_pRECO;
   Bool_t          L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50_pRECO;
   Bool_t          L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65_pRECO;
   Bool_t          L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65_pRECO;
   Bool_t          L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ_pRECO;
   Bool_t          L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ_pRECO;
   Bool_t          L1_DoubleLLPJet40_pRECO;
   Bool_t          L1_DoubleLooseIsoEG22er2p1_pRECO;
   Bool_t          L1_DoubleLooseIsoEG24er2p1_pRECO;
   Bool_t          L1_DoubleMu0_pRECO;
   Bool_t          L1_DoubleMu0_Mass_Min1_pRECO;
   Bool_t          L1_DoubleMu0_OQ_pRECO;
   Bool_t          L1_DoubleMu0_SQ_pRECO;
   Bool_t          L1_DoubleMu0_SQ_OS_pRECO;
   Bool_t          L1_DoubleMu0_Upt15_Upt7_pRECO;
   Bool_t          L1_DoubleMu0_Upt5_Upt5_pRECO;
   Bool_t          L1_DoubleMu0_Upt6_IP_Min1_Upt4_pRECO;
   Bool_t          L1_DoubleMu0_Upt6_SQ_er2p0_pRECO;
   Bool_t          L1_DoubleMu0_Upt7_SQ_er2p0_pRECO;
   Bool_t          L1_DoubleMu0_Upt8_SQ_er2p0_pRECO;
   Bool_t          L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8_pRECO;
   Bool_t          L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6_pRECO;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2_pRECO;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_pRECO;
   Bool_t          L1_DoubleMu0er1p5_SQ_pRECO;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_pRECO;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2_pRECO;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_pRECO;
   Bool_t          L1_DoubleMu0er1p5_SQ_dR_Max1p4_pRECO;
   Bool_t          L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2_pRECO;
   Bool_t          L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5_pRECO;
   Bool_t          L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6_pRECO;
   Bool_t          L1_DoubleMu0er2p0_SQ_dEta_Max1p5_pRECO;
   Bool_t          L1_DoubleMu0er2p0_SQ_dEta_Max1p6_pRECO;
   Bool_t          L1_DoubleMu18er2p1_SQ_pRECO;
   Bool_t          L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF30_HTT60er_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF40_HTT60er_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_HTT60er_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5_pRECO;
   Bool_t          L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5_pRECO;
   Bool_t          L1_DoubleMu3_SQ_HTT220er_pRECO;
   Bool_t          L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6_pRECO;
   Bool_t          L1_DoubleMu4_SQ_EG9er2p5_pRECO;
   Bool_t          L1_DoubleMu4_SQ_OS_pRECO;
   Bool_t          L1_DoubleMu4_SQ_OS_dR_Max1p2_pRECO;
   Bool_t          L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6_pRECO;
   Bool_t          L1_DoubleMu4p5_SQ_OS_pRECO;
   Bool_t          L1_DoubleMu4p5_SQ_OS_dR_Max1p2_pRECO;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_pRECO;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18_pRECO;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_pRECO;
   Bool_t          L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20_pRECO;
   Bool_t          L1_DoubleMu5_SQ_EG9er2p5_pRECO;
   Bool_t          L1_DoubleMu5_SQ_OS_dR_Max1p6_pRECO;
   Bool_t          L1_DoubleMu6_Upt6_SQ_er2p0_pRECO;
   Bool_t          L1_DoubleMu7_Upt7_SQ_er2p0_pRECO;
   Bool_t          L1_DoubleMu8_SQ_pRECO;
   Bool_t          L1_DoubleMu8_Upt8_SQ_er2p0_pRECO;
   Bool_t          L1_DoubleMu9_SQ_pRECO;
   Bool_t          L1_DoubleMu_12_5_pRECO;
   Bool_t          L1_DoubleMu_15_5_SQ_pRECO;
   Bool_t          L1_DoubleMu_15_7_pRECO;
   Bool_t          L1_DoubleMu_15_7_Mass_Min1_pRECO;
   Bool_t          L1_DoubleMu_15_7_SQ_pRECO;
   Bool_t          L1_DoubleTau70er2p1_pRECO;
   Bool_t          L1_ETM120_pRECO;
   Bool_t          L1_ETM150_pRECO;
   Bool_t          L1_ETMHF100_pRECO;
   Bool_t          L1_ETMHF100_HTT60er_pRECO;
   Bool_t          L1_ETMHF110_pRECO;
   Bool_t          L1_ETMHF110_HTT60er_pRECO;
   Bool_t          L1_ETMHF120_pRECO;
   Bool_t          L1_ETMHF120_HTT60er_pRECO;
   Bool_t          L1_ETMHF130_pRECO;
   Bool_t          L1_ETMHF130_HTT60er_pRECO;
   Bool_t          L1_ETMHF140_pRECO;
   Bool_t          L1_ETMHF150_pRECO;
   Bool_t          L1_ETMHF70_pRECO;
   Bool_t          L1_ETMHF70_HTT60er_pRECO;
   Bool_t          L1_ETMHF80_pRECO;
   Bool_t          L1_ETMHF80_HTT60er_pRECO;
   Bool_t          L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1_pRECO;
   Bool_t          L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6_pRECO;
   Bool_t          L1_ETMHF90_pRECO;
   Bool_t          L1_ETMHF90_HTT60er_pRECO;
   Bool_t          L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1_pRECO;
   Bool_t          L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6_pRECO;
   Bool_t          L1_ETT2000_pRECO;
   Bool_t          L1_FirstBunchAfterTrain_pRECO;
   Bool_t          L1_FirstBunchBeforeTrain_pRECO;
   Bool_t          L1_FirstBunchInTrain_pRECO;
   Bool_t          L1_FirstCollisionInOrbit_pRECO;
   Bool_t          L1_FirstCollisionInTrain_pRECO;
   Bool_t          L1_HCAL_LaserMon_Trig_pRECO;
   Bool_t          L1_HCAL_LaserMon_Veto_pRECO;
   Bool_t          L1_HTMHF100_pRECO;
   Bool_t          L1_HTMHF120_pRECO;
   Bool_t          L1_HTMHF125_pRECO;
   Bool_t          L1_HTMHF130_pRECO;
   Bool_t          L1_HTMHF150_pRECO;
   Bool_t          L1_HTT120_SingleLLPJet40_pRECO;
   Bool_t          L1_HTT120er_pRECO;
   Bool_t          L1_HTT160_SingleLLPJet50_pRECO;
   Bool_t          L1_HTT160er_pRECO;
   Bool_t          L1_HTT200_SingleLLPJet60_pRECO;
   Bool_t          L1_HTT200er_pRECO;
   Bool_t          L1_HTT240_SingleLLPJet70_pRECO;
   Bool_t          L1_HTT255er_pRECO;
   Bool_t          L1_HTT280er_pRECO;
   Bool_t          L1_HTT280er_QuadJet_70_55_40_35_er2p5_pRECO;
   Bool_t          L1_HTT320er_pRECO;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p5_pRECO;
   Bool_t          L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3_pRECO;
   Bool_t          L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3_pRECO;
   Bool_t          L1_HTT360er_pRECO;
   Bool_t          L1_HTT400er_pRECO;
   Bool_t          L1_HTT450er_pRECO;
   Bool_t          L1_IsoEG32er2p5_Mt40_pRECO;
   Bool_t          L1_IsoTau52er2p1_QuadJet36er2p5_pRECO;
   Bool_t          L1_IsolatedBunch_pRECO;
   Bool_t          L1_LastBunchInTrain_pRECO;
   Bool_t          L1_LastCollisionInTrain_pRECO;
   Bool_t          L1_LooseIsoEG14er2p5_HTT200er_pRECO;
   Bool_t          L1_LooseIsoEG16er2p5_HTT200er_pRECO;
   Bool_t          L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3_pRECO;
   Bool_t          L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3_pRECO;
   Bool_t          L1_LooseIsoEG24er2p1_HTT100er_pRECO;
   Bool_t          L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3_pRECO;
   Bool_t          L1_LooseIsoEG26er2p1_HTT100er_pRECO;
   Bool_t          L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_pRECO;
   Bool_t          L1_LooseIsoEG28er2p1_HTT100er_pRECO;
   Bool_t          L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_pRECO;
   Bool_t          L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_pRECO;
   Bool_t          L1_MinimumBiasHF0_pRECO;
   Bool_t          L1_MinimumBiasHF0_AND_BptxAND_pRECO;
   Bool_t          L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6_pRECO;
   Bool_t          L1_Mu12_HTT150er_pRECO;
   Bool_t          L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6_pRECO;
   Bool_t          L1_Mu14_HTT150er_pRECO;
   Bool_t          L1_Mu18er2p1_Tau24er2p1_pRECO;
   Bool_t          L1_Mu18er2p1_Tau26er2p1_pRECO;
   Bool_t          L1_Mu18er2p1_Tau26er2p1_Jet55_pRECO;
   Bool_t          L1_Mu18er2p1_Tau26er2p1_Jet70_pRECO;
   Bool_t          L1_Mu20_EG10er2p5_pRECO;
   Bool_t          L1_Mu22er2p1_IsoTau30er2p1_pRECO;
   Bool_t          L1_Mu22er2p1_IsoTau32er2p1_pRECO;
   Bool_t          L1_Mu22er2p1_IsoTau34er2p1_pRECO;
   Bool_t          L1_Mu22er2p1_IsoTau40er2p1_pRECO;
   Bool_t          L1_Mu22er2p1_Tau70er2p1_pRECO;
   Bool_t          L1_Mu3_Jet120er2p5_dR_Max0p4_pRECO;
   Bool_t          L1_Mu3_Jet16er2p5_dR_Max0p4_pRECO;
   Bool_t          L1_Mu3_Jet30er2p5_pRECO;
   Bool_t          L1_Mu3_Jet60er2p5_dR_Max0p4_pRECO;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF30_pRECO;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF40_pRECO;
   Bool_t          L1_Mu3er1p5_Jet100er2p5_ETMHF50_pRECO;
   Bool_t          L1_Mu5_EG23er2p5_pRECO;
   Bool_t          L1_Mu5_LooseIsoEG20er2p5_pRECO;
   Bool_t          L1_Mu6_DoubleEG10er2p5_pRECO;
   Bool_t          L1_Mu6_DoubleEG12er2p5_pRECO;
   Bool_t          L1_Mu6_DoubleEG15er2p5_pRECO;
   Bool_t          L1_Mu6_DoubleEG17er2p5_pRECO;
   Bool_t          L1_Mu6_HTT240er_pRECO;
   Bool_t          L1_Mu6_HTT250er_pRECO;
   Bool_t          L1_Mu7_EG20er2p5_pRECO;
   Bool_t          L1_Mu7_EG23er2p5_pRECO;
   Bool_t          L1_Mu7_LooseIsoEG20er2p5_pRECO;
   Bool_t          L1_Mu7_LooseIsoEG23er2p5_pRECO;
   Bool_t          L1_NotBptxOR_pRECO;
   Bool_t          L1_QuadJet60er2p5_pRECO;
   Bool_t          L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0_pRECO;
   Bool_t          L1_QuadMu0_pRECO;
   Bool_t          L1_QuadMu0_OQ_pRECO;
   Bool_t          L1_QuadMu0_SQ_pRECO;
   Bool_t          L1_SecondBunchInTrain_pRECO;
   Bool_t          L1_SecondLastBunchInTrain_pRECO;
   Bool_t          L1_SingleEG10er2p5_pRECO;
   Bool_t          L1_SingleEG15er2p5_pRECO;
   Bool_t          L1_SingleEG26er2p5_pRECO;
   Bool_t          L1_SingleEG28_FWD2p5_pRECO;
   Bool_t          L1_SingleEG28er1p5_pRECO;
   Bool_t          L1_SingleEG28er2p1_pRECO;
   Bool_t          L1_SingleEG28er2p5_pRECO;
   Bool_t          L1_SingleEG34er2p5_pRECO;
   Bool_t          L1_SingleEG36er2p5_pRECO;
   Bool_t          L1_SingleEG38er2p5_pRECO;
   Bool_t          L1_SingleEG40er2p5_pRECO;
   Bool_t          L1_SingleEG42er2p5_pRECO;
   Bool_t          L1_SingleEG45er2p5_pRECO;
   Bool_t          L1_SingleEG50_pRECO;
   Bool_t          L1_SingleEG60_pRECO;
   Bool_t          L1_SingleEG8er2p5_pRECO;
   Bool_t          L1_SingleIsoEG24er2p1_pRECO;
   Bool_t          L1_SingleIsoEG26er2p1_pRECO;
   Bool_t          L1_SingleIsoEG26er2p5_pRECO;
   Bool_t          L1_SingleIsoEG28_FWD2p5_pRECO;
   Bool_t          L1_SingleIsoEG28er1p5_pRECO;
   Bool_t          L1_SingleIsoEG28er2p1_pRECO;
   Bool_t          L1_SingleIsoEG28er2p5_pRECO;
   Bool_t          L1_SingleIsoEG30er2p1_pRECO;
   Bool_t          L1_SingleIsoEG30er2p5_pRECO;
   Bool_t          L1_SingleIsoEG32er2p1_pRECO;
   Bool_t          L1_SingleIsoEG32er2p5_pRECO;
   Bool_t          L1_SingleIsoEG34er2p5_pRECO;
   Bool_t          L1_SingleJet10erHE_pRECO;
   Bool_t          L1_SingleJet120_pRECO;
   Bool_t          L1_SingleJet120_FWD2p5_pRECO;
   Bool_t          L1_SingleJet120_FWD3p0_pRECO;
   Bool_t          L1_SingleJet120er1p3_pRECO;
   Bool_t          L1_SingleJet120er2p5_pRECO;
   Bool_t          L1_SingleJet12erHE_pRECO;
   Bool_t          L1_SingleJet140er2p5_pRECO;
   Bool_t          L1_SingleJet160er2p5_pRECO;
   Bool_t          L1_SingleJet180_pRECO;
   Bool_t          L1_SingleJet180er2p5_pRECO;
   Bool_t          L1_SingleJet200_pRECO;
   Bool_t          L1_SingleJet20er2p5_NotBptxOR_pRECO;
   Bool_t          L1_SingleJet20er2p5_NotBptxOR_3BX_pRECO;
   Bool_t          L1_SingleJet35_pRECO;
   Bool_t          L1_SingleJet35_FWD2p5_pRECO;
   Bool_t          L1_SingleJet35_FWD3p0_pRECO;
   Bool_t          L1_SingleJet35er1p3_pRECO;
   Bool_t          L1_SingleJet35er2p5_pRECO;
   Bool_t          L1_SingleJet43er2p5_NotBptxOR_3BX_pRECO;
   Bool_t          L1_SingleJet46er2p5_NotBptxOR_3BX_pRECO;
   Bool_t          L1_SingleJet60_pRECO;
   Bool_t          L1_SingleJet60_FWD2p5_pRECO;
   Bool_t          L1_SingleJet8erHE_pRECO;
   Bool_t          L1_SingleJet90_pRECO;
   Bool_t          L1_SingleJet90_FWD2p5_pRECO;
   Bool_t          L1_SingleLooseIsoEG26er1p5_pRECO;
   Bool_t          L1_SingleLooseIsoEG26er2p5_pRECO;
   Bool_t          L1_SingleLooseIsoEG28_FWD2p5_pRECO;
   Bool_t          L1_SingleLooseIsoEG28er1p5_pRECO;
   Bool_t          L1_SingleLooseIsoEG28er2p1_pRECO;
   Bool_t          L1_SingleLooseIsoEG28er2p5_pRECO;
   Bool_t          L1_SingleLooseIsoEG30er1p5_pRECO;
   Bool_t          L1_SingleLooseIsoEG30er2p5_pRECO;
   Bool_t          L1_SingleMu0_BMTF_pRECO;
   Bool_t          L1_SingleMu0_DQ_pRECO;
   Bool_t          L1_SingleMu0_EMTF_pRECO;
   Bool_t          L1_SingleMu0_OMTF_pRECO;
   Bool_t          L1_SingleMu0_SQ13_BMTF_pRECO;
   Bool_t          L1_SingleMu0_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu0_SQ15_BMTF_pRECO;
   Bool_t          L1_SingleMu0_Upt10_pRECO;
   Bool_t          L1_SingleMu0_Upt10_BMTF_pRECO;
   Bool_t          L1_SingleMu0_Upt10_EMTF_pRECO;
   Bool_t          L1_SingleMu0_Upt10_OMTF_pRECO;
   Bool_t          L1_SingleMu0_Upt10_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu0_Upt15_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu0_Upt20_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu0_Upt25_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu10_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu11_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu12_DQ_BMTF_pRECO;
   Bool_t          L1_SingleMu12_DQ_EMTF_pRECO;
   Bool_t          L1_SingleMu12_DQ_OMTF_pRECO;
   Bool_t          L1_SingleMu15_DQ_pRECO;
   Bool_t          L1_SingleMu18_pRECO;
   Bool_t          L1_SingleMu20_pRECO;
   Bool_t          L1_SingleMu22_pRECO;
   Bool_t          L1_SingleMu22_BMTF_pRECO;
   Bool_t          L1_SingleMu22_BMTF_NEG_pRECO;
   Bool_t          L1_SingleMu22_BMTF_POS_pRECO;
   Bool_t          L1_SingleMu22_DQ_pRECO;
   Bool_t          L1_SingleMu22_EMTF_pRECO;
   Bool_t          L1_SingleMu22_EMTF_NEG_pRECO;
   Bool_t          L1_SingleMu22_EMTF_POS_pRECO;
   Bool_t          L1_SingleMu22_OMTF_pRECO;
   Bool_t          L1_SingleMu22_OMTF_NEG_pRECO;
   Bool_t          L1_SingleMu22_OMTF_POS_pRECO;
   Bool_t          L1_SingleMu22_OQ_pRECO;
   Bool_t          L1_SingleMu25_pRECO;
   Bool_t          L1_SingleMu3_pRECO;
   Bool_t          L1_SingleMu5_pRECO;
   Bool_t          L1_SingleMu5_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu6_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu7_pRECO;
   Bool_t          L1_SingleMu7_DQ_pRECO;
   Bool_t          L1_SingleMu7_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu8_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMu9_SQ14_BMTF_pRECO;
   Bool_t          L1_SingleMuCosmics_pRECO;
   Bool_t          L1_SingleMuCosmics_BMTF_pRECO;
   Bool_t          L1_SingleMuCosmics_EMTF_pRECO;
   Bool_t          L1_SingleMuCosmics_OMTF_pRECO;
   Bool_t          L1_SingleMuOpen_pRECO;
   Bool_t          L1_SingleMuOpen_BMTF_pRECO;
   Bool_t          L1_SingleMuOpen_EMTF_pRECO;
   Bool_t          L1_SingleMuOpen_NotBptxOR_pRECO;
   Bool_t          L1_SingleMuOpen_OMTF_pRECO;
   Bool_t          L1_SingleMuOpen_er1p1_NotBptxOR_3BX_pRECO;
   Bool_t          L1_SingleMuOpen_er1p4_NotBptxOR_3BX_pRECO;
   Bool_t          L1_SingleMuShower_Nominal_pRECO;
   Bool_t          L1_SingleMuShower_Tight_pRECO;
   Bool_t          L1_SingleTau120er2p1_pRECO;
   Bool_t          L1_SingleTau130er2p1_pRECO;
   Bool_t          L1_TOTEM_1_pRECO;
   Bool_t          L1_TOTEM_2_pRECO;
   Bool_t          L1_TOTEM_3_pRECO;
   Bool_t          L1_TOTEM_4_pRECO;
   Bool_t          L1_TripleEG_18_17_8_er2p5_pRECO;
   Bool_t          L1_TripleEG_18_18_12_er2p5_pRECO;
   Bool_t          L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5_pRECO;
   Bool_t          L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5_pRECO;
   Bool_t          L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5_pRECO;
   Bool_t          L1_TripleMu0_pRECO;
   Bool_t          L1_TripleMu0_OQ_pRECO;
   Bool_t          L1_TripleMu0_SQ_pRECO;
   Bool_t          L1_TripleMu3_pRECO;
   Bool_t          L1_TripleMu3_SQ_pRECO;
   Bool_t          L1_TripleMu_3SQ_2p5SQ_0_pRECO;
   Bool_t          L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12_pRECO;
   Bool_t          L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12_pRECO;
   Bool_t          L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12_pRECO;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ_pRECO;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO;
   Bool_t          L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO;
   Bool_t          L1_TripleMu_5_3_3_pRECO;
   Bool_t          L1_TripleMu_5_3_3_SQ_pRECO;
   Bool_t          L1_TripleMu_5_3p5_2p5_pRECO;
   Bool_t          L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO;
   Bool_t          L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO;
   Bool_t          L1_TripleMu_5_5_3_pRECO;
   Bool_t          L1_TwoMuShower_Loose_pRECO;
   Bool_t          L1_UnpairedBunchBptxMinus_pRECO;
   Bool_t          L1_UnpairedBunchBptxPlus_pRECO;
   Bool_t          L1_ZeroBias_pRECO;
   Bool_t          L1_ZeroBias_copy_pRECO;
   Bool_t          L1_UnprefireableEvent_TriggerRules_pRECO;
   Bool_t          L1_UnprefireableEvent_FirstBxInTrain_pRECO;
   Bool_t          L1_FinalOR_BXmin1_pRECO;
   Bool_t          L1_FinalOR_BXmin2_pRECO;
   Bool_t          HLTriggerFirstPath;
   Bool_t          DST_PFScouting_DoubleMuon;
   Bool_t          DST_PFScouting_DoubleEG;
   Bool_t          DST_PFScouting_JetHT;
   Bool_t          DST_PFScouting_DatasetMuon;
   Bool_t          DST_PFScouting_AXOVLoose;
   Bool_t          DST_PFScouting_AXOLoose;
   Bool_t          DST_PFScouting_AXONominal;
   Bool_t          DST_PFScouting_AXOTight;
   Bool_t          DST_PFScouting_AXOVTight;
   Bool_t          DST_PFScouting_CICADAVLoose;
   Bool_t          DST_PFScouting_CICADALoose;
   Bool_t          DST_PFScouting_CICADAMedium;
   Bool_t          DST_PFScouting_CICADATight;
   Bool_t          DST_PFScouting_CICADAVTight;
   Bool_t          DST_PFScouting_SingleMuon;
   Bool_t          DST_PFScouting_SinglePhotonEB;
   Bool_t          DST_PFScouting_ZeroBias;
   Bool_t          HLT_EphemeralPhysics;
   Bool_t          HLT_TestPhysics;
   Bool_t          HLT_EphemeralZeroBias;
   Bool_t          HLT_SpecialZeroBias;
   Bool_t          HLT_EcalCalibration;
   Bool_t          HLT_HcalCalibration;
   Bool_t          HLT_HcalNZS;
   Bool_t          HLT_HcalPhiSym;
   Bool_t          HLT_Random;
   Bool_t          HLT_Physics;
   Bool_t          HLT_ZeroBias;
   Bool_t          HLT_ZeroBias_Alignment;
   Bool_t          HLT_ZeroBias_Beamspot;
   Bool_t          HLT_ZeroBias_IsolatedBunches;
   Bool_t          HLT_ZeroBias_FirstBXAfterTrain;
   Bool_t          HLT_ZeroBias_FirstCollisionAfterAbortGap;
   Bool_t          HLT_ZeroBias_FirstCollisionInTrain;
   Bool_t          HLT_ZeroBias_LastCollisionInTrain;
   Bool_t          HLT_HT300_Beamspot;
   Bool_t          HLT_IsoTrackHB;
   Bool_t          HLT_IsoTrackHE;
   Bool_t          HLT_PFJet40_GPUvsCPU;
   Bool_t          HLT_AK8PFJet380_SoftDropMass30;
   Bool_t          HLT_AK8PFJet400_SoftDropMass30;
   Bool_t          HLT_AK8PFJet425_SoftDropMass30;
   Bool_t          HLT_AK8PFJet450_SoftDropMass30;
   Bool_t          HLT_AK8DiPFJet250_250_SoftDropMass40;
   Bool_t          HLT_AK8DiPFJet250_250_SoftDropMass50;
   Bool_t          HLT_AK8DiPFJet260_260_SoftDropMass30;
   Bool_t          HLT_AK8DiPFJet260_260_SoftDropMass40;
   Bool_t          HLT_AK8DiPFJet270_270_SoftDropMass30;
   Bool_t          HLT_AK8DiPFJet280_280_SoftDropMass30;
   Bool_t          HLT_AK8DiPFJet290_290_SoftDropMass30;
   Bool_t          HLT_CaloJet500_NoJetID;
   Bool_t          HLT_CaloJet550_NoJetID;
   Bool_t          HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;
   Bool_t          HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;
   Bool_t          HLT_Trimuon5_3p5_2_Upsilon_Muon;
   Bool_t          HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;
   Bool_t          HLT_DoubleEle25_CaloIdL_MW;
   Bool_t          HLT_DoubleEle27_CaloIdL_MW;
   Bool_t          HLT_DoubleEle33_CaloIdL_MW;
   Bool_t          HLT_DoubleEle24_eta2p1_WPTight_Gsf;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;
   Bool_t          HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Mu27_Ele37_CaloIdL_MW;
   Bool_t          HLT_Mu37_Ele27_CaloIdL_MW;
   Bool_t          HLT_Mu37_TkMu27;
   Bool_t          HLT_DoubleMu4_3_Bs;
   Bool_t          HLT_DoubleMu4_3_Jpsi;
   Bool_t          HLT_DoubleMu4_3_LowMass;
   Bool_t          HLT_DoubleMu4_LowMass_Displaced;
   Bool_t          HLT_Mu0_L1DoubleMu;
   Bool_t          HLT_Mu4_L1DoubleMu;
   Bool_t          HLT_DoubleMu2_Jpsi_LowPt;
   Bool_t          HLT_DoubleMu4_3_Photon4_BsToMMG;
   Bool_t          HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu;
   Bool_t          HLT_DoubleMu3_TkMu_DsTau3Mu;
   Bool_t          HLT_DoubleMu4_Mass3p8_DZ_PFHT350;
   Bool_t          HLT_DoubleMu4_MuMuTrk_Displaced;
   Bool_t          HLT_Mu3_PFJet40;
   Bool_t          HLT_Mu7p5_L2Mu2_Jpsi;
   Bool_t          HLT_Mu7p5_L2Mu2_Upsilon;
   Bool_t          HLT_Mu3_L1SingleMu5orSingleMu7;
   Bool_t          HLT_Mu0_Barrel;
   Bool_t          HLT_Mu0_Barrel_L1HP6;
   Bool_t          HLT_Mu0_Barrel_L1HP7;
   Bool_t          HLT_Mu0_Barrel_L1HP8;
   Bool_t          HLT_Mu0_Barrel_L1HP9;
   Bool_t          HLT_Mu0_Barrel_L1HP10;
   Bool_t          HLT_Mu0_Barrel_L1HP11;
   Bool_t          HLT_Mu0_Barrel_L1HP6_IP6;
   Bool_t          HLT_Mu6_Barrel_L1HP7_IP6;
   Bool_t          HLT_Mu7_Barrel_L1HP8_IP6;
   Bool_t          HLT_Mu8_Barrel_L1HP9_IP6;
   Bool_t          HLT_Mu9_Barrel_L1HP10_IP6;
   Bool_t          HLT_Mu10_Barrel_L1HP11_IP6;
   Bool_t          HLT_DoublePhoton33_CaloIdL;
   Bool_t          HLT_DoublePhoton70;
   Bool_t          HLT_DoublePhoton85;
   Bool_t          HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;
   Bool_t          HLT_Ele30_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele35_WPTight_Gsf;
   Bool_t          HLT_Ele38_WPTight_Gsf;
   Bool_t          HLT_Ele40_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   Bool_t          HLT_IsoMu20;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX;
   Bool_t          HLT_UncorrectedJetE30_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE60_NoBPTX3BX;
   Bool_t          HLT_UncorrectedJetE70_NoBPTX3BX;
   Bool_t          HLT_L1SingleMuCosmics;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX3BX;
   Bool_t          HLT_L2Mu10_NoVertex_NoBPTX;
   Bool_t          HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
   Bool_t          HLT_L2Mu23NoVtx_2Cha;
   Bool_t          HLT_L2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t          HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t          HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;
   Bool_t          HLT_DoubleL2Mu50;
   Bool_t          HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;
   Bool_t          HLT_DoubleL2Mu23NoVtx_2Cha;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha;
   Bool_t          HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30;
   Bool_t          HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;
   Bool_t          HLT_Mu30_TkMu0_Psi;
   Bool_t          HLT_Mu30_TkMu0_Upsilon;
   Bool_t          HLT_Mu25_TkMu0_Phi;
   Bool_t          HLT_Mu15;
   Bool_t          HLT_Mu20;
   Bool_t          HLT_Mu27;
   Bool_t          HLT_Mu50;
   Bool_t          HLT_Mu55;
   Bool_t          HLT_CascadeMu100;
   Bool_t          HLT_HighPtTkMu100;
   Bool_t          HLT_DiPFJetAve40;
   Bool_t          HLT_DiPFJetAve60;
   Bool_t          HLT_DiPFJetAve80;
   Bool_t          HLT_DiPFJetAve140;
   Bool_t          HLT_DiPFJetAve200;
   Bool_t          HLT_DiPFJetAve260;
   Bool_t          HLT_DiPFJetAve320;
   Bool_t          HLT_DiPFJetAve400;
   Bool_t          HLT_DiPFJetAve500;
   Bool_t          HLT_DiPFJetAve60_HFJEC;
   Bool_t          HLT_DiPFJetAve80_HFJEC;
   Bool_t          HLT_DiPFJetAve100_HFJEC;
   Bool_t          HLT_DiPFJetAve160_HFJEC;
   Bool_t          HLT_DiPFJetAve220_HFJEC;
   Bool_t          HLT_DiPFJetAve260_HFJEC;
   Bool_t          HLT_DiPFJetAve300_HFJEC;
   Bool_t          HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP;
   Bool_t          HLT_AK8PFJet40;
   Bool_t          HLT_AK8PFJet60;
   Bool_t          HLT_AK8PFJet80;
   Bool_t          HLT_AK8PFJet140;
   Bool_t          HLT_AK8PFJet200;
   Bool_t          HLT_AK8PFJet260;
   Bool_t          HLT_AK8PFJet320;
   Bool_t          HLT_AK8PFJet400;
   Bool_t          HLT_AK8PFJet450;
   Bool_t          HLT_AK8PFJet500;
   Bool_t          HLT_AK8PFJet550;
   Bool_t          HLT_PFJet40;
   Bool_t          HLT_PFJet40_ZeroBiasCopy;
   Bool_t          HLT_PFJet60;
   Bool_t          HLT_PFJet80;
   Bool_t          HLT_PFJet80_L1SingleJet60;
   Bool_t          HLT_PFJet110;
   Bool_t          HLT_PFJet140;
   Bool_t          HLT_PFJet200;
   Bool_t          HLT_PFJet260;
   Bool_t          HLT_PFJet320;
   Bool_t          HLT_PFJet400;
   Bool_t          HLT_PFJet450;
   Bool_t          HLT_PFJet500;
   Bool_t          HLT_PFJet550;
   Bool_t          HLT_PFJetFwd40;
   Bool_t          HLT_PFJetFwd60;
   Bool_t          HLT_PFJetFwd80;
   Bool_t          HLT_PFJetFwd140;
   Bool_t          HLT_PFJetFwd200;
   Bool_t          HLT_PFJetFwd260;
   Bool_t          HLT_PFJetFwd320;
   Bool_t          HLT_PFJetFwd400;
   Bool_t          HLT_PFJetFwd450;
   Bool_t          HLT_PFJetFwd500;
   Bool_t          HLT_AK8PFJetFwd40;
   Bool_t          HLT_AK8PFJetFwd60;
   Bool_t          HLT_AK8PFJetFwd80;
   Bool_t          HLT_AK8PFJetFwd140;
   Bool_t          HLT_AK8PFJetFwd200;
   Bool_t          HLT_AK8PFJetFwd260;
   Bool_t          HLT_AK8PFJetFwd320;
   Bool_t          HLT_AK8PFJetFwd400;
   Bool_t          HLT_AK8PFJetFwd450;
   Bool_t          HLT_AK8PFJetFwd500;
   Bool_t          HLT_PFHT180;
   Bool_t          HLT_PFHT250;
   Bool_t          HLT_PFHT370;
   Bool_t          HLT_PFHT430;
   Bool_t          HLT_PFHT510;
   Bool_t          HLT_PFHT590;
   Bool_t          HLT_PFHT680;
   Bool_t          HLT_PFHT780;
   Bool_t          HLT_PFHT890;
   Bool_t          HLT_PFHT1050;
   Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_PFHT500_PFMET110_PFMHT110_IDTight;
   Bool_t          HLT_PFHT700_PFMET85_PFMHT85_IDTight;
   Bool_t          HLT_PFHT800_PFMET75_PFMHT75_IDTight;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight;
   Bool_t          HLT_PFMET130_PFMHT130_IDTight;
   Bool_t          HLT_PFMET140_PFMHT140_IDTight;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_PFHT60;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
   Bool_t          HLT_PFMETTypeOne140_PFMHT140_IDTight;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
   Bool_t          HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;
   Bool_t          HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;
   Bool_t          HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;
   Bool_t          HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;
   Bool_t          HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;
   Bool_t          HLT_L1ETMHadSeeds;
   Bool_t          HLT_CaloMHT90;
   Bool_t          HLT_CaloMET90_NotCleaned;
   Bool_t          HLT_CaloMET350_NotCleaned;
   Bool_t          HLT_PFMET200_NotCleaned;
   Bool_t          HLT_PFMET250_NotCleaned;
   Bool_t          HLT_PFMET300_NotCleaned;
   Bool_t          HLT_PFMET200_BeamHaloCleaned;
   Bool_t          HLT_PFMETTypeOne200_BeamHaloCleaned;
   Bool_t          HLT_MET105_IsoTrk50;
   Bool_t          HLT_MET120_IsoTrk50;
   Bool_t          HLT_Mu12eta2p3;
   Bool_t          HLT_Mu12eta2p3_PFJet40;
   Bool_t          HLT_Mu12_DoublePFJets40_PNetBTag_0p11;
   Bool_t          HLT_Mu12_DoublePFJets100_PNetBTag_0p11;
   Bool_t          HLT_Mu12_DoublePFJets200_PNetBTag_0p11;
   Bool_t          HLT_Mu12_DoublePFJets350_PNetBTag_0p11;
   Bool_t          HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11;
   Bool_t          HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11;
   Bool_t          HLT_DoublePFJets40_PNetBTag_0p11;
   Bool_t          HLT_DoublePFJets100_PNetBTag_0p11;
   Bool_t          HLT_DoublePFJets200_PNetBTag_0p11;
   Bool_t          HLT_DoublePFJets350_PNetBTag_0p11;
   Bool_t          HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11;
   Bool_t          HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11;
   Bool_t          HLT_Photon300_NoHE;
   Bool_t          HLT_Mu8_TrkIsoVVL;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;
   Bool_t          HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu17_TrkIsoVVL;
   Bool_t          HLT_Mu19_TrkIsoVVL;
   Bool_t          HLT_BTagMu_AK4DiJet20_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet40_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet70_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet110_Mu5;
   Bool_t          HLT_BTagMu_AK4DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK4Jet300_Mu5;
   Bool_t          HLT_BTagMu_AK8DiJet170_Mu5;
   Bool_t          HLT_BTagMu_AK8Jet170_DoubleMu5;
   Bool_t          HLT_BTagMu_AK8Jet300_Mu5;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53;
   Bool_t          HLT_Photon33;
   Bool_t          HLT_Photon50;
   Bool_t          HLT_Photon75;
   Bool_t          HLT_Photon90;
   Bool_t          HLT_Photon120;
   Bool_t          HLT_Photon150;
   Bool_t          HLT_Photon175;
   Bool_t          HLT_Photon200;
   Bool_t          HLT_Photon45EB;
   Bool_t          HLT_Photon40EB;
   Bool_t          HLT_Photon50EB;
   Bool_t          HLT_Photon30EB_TightID_TightIso;
   Bool_t          HLT_Photon40EB_TightID_TightIso;
   Bool_t          HLT_Photon45EB_TightID_TightIso;
   Bool_t          HLT_Photon50EB_TightID_TightIso_PFJet30;
   Bool_t          HLT_Photon50EB_TightID_TightIso_CaloJet30;
   Bool_t          HLT_Photon50EB_TightID_TightIso_AK8PFJet30;
   Bool_t          HLT_Photon50EB_TightID_TightIso_AK8CaloJet30;
   Bool_t          HLT_Photon50EB_TightID_TightIso;
   Bool_t          HLT_Photon55EB_TightID_TightIso;
   Bool_t          HLT_Photon75EB_TightID_TightIso;
   Bool_t          HLT_Photon90EB_TightID_TightIso;
   Bool_t          HLT_Photon110EB_TightID_TightIso;
   Bool_t          HLT_Photon110EB_TightID_TightIso_PFJet30;
   Bool_t          HLT_Photon110EB_TightID_TightIso_CaloJet30;
   Bool_t          HLT_Photon110EB_TightID_TightIso_AK8PFJet30;
   Bool_t          HLT_Photon110EB_TightID_TightIso_AK8CaloJet30;
   Bool_t          HLT_Photon100EBHE10;
   Bool_t          HLT_Photon50_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon90_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon120_R9Id90_HE10_IsoM;
   Bool_t          HLT_Photon165_R9Id90_HE10_IsoM;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;
   Bool_t          HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;
   Bool_t          HLT_DiphotonMVA14p25_Mass90;
   Bool_t          HLT_DiphotonMVA14p25_Tight_Mass90;
   Bool_t          HLT_Photon35_TwoProngs35;
   Bool_t          HLT_IsoMu24_TwoProngs35;
   Bool_t          HLT_Dimuon0_Jpsi_L1_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_NoOS;
   Bool_t          HLT_Dimuon0_Jpsi;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing;
   Bool_t          HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;
   Bool_t          HLT_Dimuon0_Jpsi3p5_Muon2;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0;
   Bool_t          HLT_Dimuon0_Upsilon_L1_4p5er2p0M;
   Bool_t          HLT_Dimuon0_Upsilon_NoVertexing;
   Bool_t          HLT_Dimuon0_LowMass_L1_0er1p5;
   Bool_t          HLT_Dimuon0_LowMass;
   Bool_t          HLT_Dimuon0_LowMass_L1_4;
   Bool_t          HLT_Dimuon0_LowMass_L1_TM530;
   Bool_t          HLT_Dimuon0_Upsilon_Muon_NoL1Mass;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8_DZ;
   Bool_t          HLT_TripleMu_10_5_5_DZ;
   Bool_t          HLT_TripleMu_12_10_5;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;
   Bool_t          HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;
   Bool_t          HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
   Bool_t          HLT_DoubleMu3_DZ_PFMET70_PFMHT70;
   Bool_t          HLT_DoubleMu3_DZ_PFMET90_PFMHT90;
   Bool_t          HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;
   Bool_t          HLT_DoubleMu4_3_LowMass_SS;
   Bool_t          HLT_DoubleMu4_Jpsi_Displaced;
   Bool_t          HLT_DoubleMu4_Jpsi_NoVertexing;
   Bool_t          HLT_DoubleMu4_JpsiTrkTrk_Displaced;
   Bool_t          HLT_DoubleMu4_JpsiTrk_Bc;
   Bool_t          HLT_DoubleMu43NoFiltersNoVtx;
   Bool_t          HLT_DoubleMu48NoFiltersNoVtx;
   Bool_t          HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;
   Bool_t          HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;
   Bool_t          HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;
   Bool_t          HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;
   Bool_t          HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;
   Bool_t          HLT_Ele28_eta2p1_WPTight_Gsf_HT150;
   Bool_t          HLT_Ele28_HighEta_SC20_Mass55;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele50_IsoVVVL_PFHT450;
   Bool_t          HLT_Ele15_IsoVVVL_PFHT600;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450_PFMET50;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu50_IsoVVVL_PFHT450;
   Bool_t          HLT_Mu15_IsoVVVL_PFHT600;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;
   Bool_t          HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;
   Bool_t          HLT_Dimuon10_Upsilon_y1p4;
   Bool_t          HLT_Dimuon12_Upsilon_y1p4;
   Bool_t          HLT_Dimuon14_Phi_Barrel_Seagulls;
   Bool_t          HLT_Dimuon25_Jpsi;
   Bool_t          HLT_Dimuon14_PsiPrime;
   Bool_t          HLT_Dimuon14_PsiPrime_noCorrL1;
   Bool_t          HLT_Dimuon18_PsiPrime;
   Bool_t          HLT_Dimuon18_PsiPrime_noCorrL1;
   Bool_t          HLT_Dimuon24_Upsilon_noCorrL1;
   Bool_t          HLT_Dimuon24_Phi_noCorrL1;
   Bool_t          HLT_Dimuon25_Jpsi_noCorrL1;
   Bool_t          HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
   Bool_t          HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
   Bool_t          HLT_DoubleIsoMu20_eta2p1;
   Bool_t          HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;
   Bool_t          HLT_Mu8;
   Bool_t          HLT_Mu17;
   Bool_t          HLT_Mu19;
   Bool_t          HLT_Mu17_Photon30_IsoCaloId;
   Bool_t          HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53;
   Bool_t          HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
   Bool_t          HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_Ele135_CaloIdVT_GsfTrkIdT;
   Bool_t          HLT_PFHT330PT30_QuadPFJet_75_60_45_40;
   Bool_t          HLT_PFHT400_SixPFJet32;
   Bool_t          HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50;
   Bool_t          HLT_PFHT450_SixPFJet36;
   Bool_t          HLT_PFHT450_SixPFJet36_PNetBTag0p35;
   Bool_t          HLT_PFHT400_FivePFJet_120_120_60_30_30;
   Bool_t          HLT_PFHT350;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380;
   Bool_t          HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400;
   Bool_t          HLT_ECALHT800;
   Bool_t          HLT_DiSC30_18_EIso_AND_HE_Mass70;
   Bool_t          HLT_Photon20_HoverELoose;
   Bool_t          HLT_Photon30_HoverELoose;
   Bool_t          HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          HLT_CDC_L2cosmic_10_er1p0;
   Bool_t          HLT_CDC_L2cosmic_5p5_er1p0;
   Bool_t          HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
   Bool_t          HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1;
   Bool_t          HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t          HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;
   Bool_t          HLT_Mu18_Mu9_SameSign;
   Bool_t          HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;
   Bool_t          HLT_DoubleMu3_DCA_PFMET50_PFMHT60;
   Bool_t          HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0;
   Bool_t          HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA;
   Bool_t          HLT_TripleMu_5_3_3_Mass3p8_DCA;
   Bool_t          HLT_QuadPFJet103_88_75_15;
   Bool_t          HLT_QuadPFJet105_88_76_15;
   Bool_t          HLT_QuadPFJet111_90_80_15;
   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;
   Bool_t          HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;
   Bool_t          HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1;
   Bool_t          HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1;
   Bool_t          HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1;
   Bool_t          HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1;
   Bool_t          HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5;
   Bool_t          HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3;
   Bool_t          HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0;
   Bool_t          HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6;
   Bool_t          HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3;
   Bool_t          HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2;
   Bool_t          HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1;
   Bool_t          HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2;
   Bool_t          HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1;
   Bool_t          HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2;
   Bool_t          HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5;
   Bool_t          HLT_PFHT250_QuadPFJet25;
   Bool_t          HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55;
   Bool_t          HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50;
   Bool_t          HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55;
   Bool_t          HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50;
   Bool_t          HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50;
   Bool_t          HLT_PFHT280_QuadPFJet30;
   Bool_t          HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55;
   Bool_t          HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60;
   Bool_t          HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60;
   Bool_t          HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50;
   Bool_t          HLT_QuadPFJet100_88_70_30;
   Bool_t          HLT_QuadPFJet105_88_75_30;
   Bool_t          HLT_QuadPFJet111_90_80_30;
   Bool_t          HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight;
   Bool_t          HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight;
   Bool_t          HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight;
   Bool_t          HLT_AK8PFJet220_SoftDropMass40;
   Bool_t          HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50;
   Bool_t          HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53;
   Bool_t          HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55;
   Bool_t          HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60;
   Bool_t          HLT_AK8PFJet230_SoftDropMass40;
   Bool_t          HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;
   Bool_t          HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10;
   Bool_t          HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03;
   Bool_t          HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05;
   Bool_t          HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06;
   Bool_t          HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10;
   Bool_t          HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03;
   Bool_t          HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05;
   Bool_t          HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06;
   Bool_t          HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10;
   Bool_t          HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03;
   Bool_t          HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05;
   Bool_t          HLT_AK8PFJet275_Nch45;
   Bool_t          HLT_AK8PFJet275_Nch40;
   Bool_t          HLT_IsoMu50_AK8PFJet220_SoftDropMass40;
   Bool_t          HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06;
   Bool_t          HLT_IsoMu50_AK8PFJet230_SoftDropMass40;
   Bool_t          HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06;
   Bool_t          HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06;
   Bool_t          HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20;
   Bool_t          HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60;
   Bool_t          HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75;
   Bool_t          HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring;
   Bool_t          HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PFHT250;
   Bool_t          HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25;
   Bool_t          HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50;
   Bool_t          HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50;
   Bool_t          HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3;
   Bool_t          HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3;
   Bool_t          HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60;
   Bool_t          HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1;
   Bool_t          HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1;
   Bool_t          HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3;
   Bool_t          HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3;
   Bool_t          HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3;
   Bool_t          HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm;
   Bool_t          HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm;
   Bool_t          HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm;
   Bool_t          HLT_L2Mu10NoVtx_2Cha;
   Bool_t          HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_L3Mu10NoVtx;
   Bool_t          HLT_L3Mu10NoVtx_DxyMin0p01cm;
   Bool_t          HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm;
   Bool_t          HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm;
   Bool_t          HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_L2Mu10NoVtx_2Cha_CosmicSeed;
   Bool_t          HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm;
   Bool_t          HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm;
   Bool_t          HLT_L3Mu30NoVtx_DxyMin0p01cm;
   Bool_t          HLT_L3Mu50NoVtx_DxyMin0p01cm;
   Bool_t          HLT_L3dTksMu10_NoVtx_DxyMin0p01cm;
   Bool_t          HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId;
   Bool_t          HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24;
   Bool_t          HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34;
   Bool_t          HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1;
   Bool_t          HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1;
   Bool_t          HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy;
   Bool_t          HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1;
   Bool_t          HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive;
   Bool_t          HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive;
   Bool_t          HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive;
   Bool_t          HLT_HT350_DelayedJet40_SingleDelay3nsInclusive;
   Bool_t          HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive;
   Bool_t          HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay2nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1nsTrackless;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless;
   Bool_t          HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless;
   Bool_t          HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive;
   Bool_t          HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless;
   Bool_t          HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless;
   Bool_t          HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless;
   Bool_t          HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless;
   Bool_t          HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless;
   Bool_t          HLT_L1Mu6HT240;
   Bool_t          HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose;
   Bool_t          HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5;
   Bool_t          HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose;
   Bool_t          HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5;
   Bool_t          HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose;
   Bool_t          HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5;
   Bool_t          HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5;
   Bool_t          HLT_HT350;
   Bool_t          HLT_HT425;
   Bool_t          HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT400_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT430_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT550_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT650_DisplacedDijet60_Inclusive;
   Bool_t          HLT_CaloMET60_DTCluster50;
   Bool_t          HLT_CaloMET60_DTClusterNoMB1S50;
   Bool_t          HLT_L1MET_DTCluster50;
   Bool_t          HLT_L1MET_DTClusterNoMB1S50;
   Bool_t          HLT_CscCluster_Loose;
   Bool_t          HLT_CscCluster_Medium;
   Bool_t          HLT_CscCluster_Tight;
   Bool_t          HLT_CscCluster50_Photon20Unseeded;
   Bool_t          HLT_CscCluster50_Photon30Unseeded;
   Bool_t          HLT_CscCluster100_Ele5;
   Bool_t          HLT_CscCluster100_Mu5;
   Bool_t          HLT_CscCluster100_PNetTauhPFJet10_Loose;
   Bool_t          HLT_DoubleCscCluster75;
   Bool_t          HLT_IsoTrk200_L1SingleMuShower;
   Bool_t          HLT_IsoTrk400_L1SingleMuShower;
   Bool_t          HLT_DoubleCscCluster100;
   Bool_t          HLT_L1CSCShower_DTCluster50;
   Bool_t          HLT_L1CSCShower_DTCluster75;
   Bool_t          HLT_PFMET105_IsoTrk50;
   Bool_t          HLT_L1SingleLLPJet;
   Bool_t          HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack;
   Bool_t          HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;
   Bool_t          HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;
   Bool_t          HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive;
   Bool_t          HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless;
   Bool_t          HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive;
   Bool_t          HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless;
   Bool_t          HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive;
   Bool_t          HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7;
   Bool_t          HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8;
   Bool_t          HLT_DiPhoton10Time1ns;
   Bool_t          HLT_DiPhoton10Time1p2ns;
   Bool_t          HLT_DiPhoton10Time1p4ns;
   Bool_t          HLT_DiPhoton10Time1p6ns;
   Bool_t          HLT_DiPhoton10Time1p8ns;
   Bool_t          HLT_DiPhoton10Time2ns;
   Bool_t          HLT_DiPhoton10_CaloIdL;
   Bool_t          HLT_DoubleEle6p5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle8_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle10_eta1p22_mMax6;
   Bool_t          HLT_SingleEle8;
   Bool_t          HLT_SingleEle8_SingleEGL1;
   Bool_t          HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t          HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;
   Bool_t          HLT_Mu50_L1SingleMuShower;
   Bool_t          HLT_IsoMu24_OneProng32;
   Bool_t          HLT_Photon32_OneProng32_M50To105;
   Bool_t          HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng;
   Bool_t          HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1;
   Bool_t          HLT_VBF_DiPFJet125_45_Mjj1050;
   Bool_t          HLT_VBF_DiPFJet125_45_Mjj1200;
   Bool_t          HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60;
   Bool_t          HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60;
   Bool_t          HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85;
   Bool_t          HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85;
   Bool_t          HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL;
   Bool_t          HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL;
   Bool_t          HLT_VBF_DiPFJet50_Mjj650_Photon22;
   Bool_t          HLT_VBF_DiPFJet50_Mjj750_Photon22;
   Bool_t          HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf;
   Bool_t          HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf;
   Bool_t          HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1;
   Bool_t          HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1;
   Bool_t          HLT_VBF_DoublePNetTauhPFJet20_eta2p2;
   Bool_t          HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3;
   Bool_t          HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3;
   Bool_t          HLT_PFJet200_TimeLtNeg2p5ns;
   Bool_t          HLT_PFJet200_TimeGt2p5ns;
   Bool_t          HLT_Photon50_TimeLtNeg2p5ns;
   Bool_t          HLT_Photon50_TimeGt2p5ns;
   Bool_t          HLT_PPSMaxTracksPerArm1;
   Bool_t          HLT_PPSMaxTracksPerRP4;
   Bool_t          HLT_PPSRandom;
   Bool_t          HLT_L1AXOVTight;
   Bool_t          HLTriggerFinalPath;
   Bool_t          Dataset_ScoutingPFMonitor;
   Bool_t          Dataset_ScoutingPFRun3;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!
   TBranch        *b_BeamSpot_type;   //!
   TBranch        *b_BeamSpot_sigmaZ;   //!
   TBranch        *b_BeamSpot_sigmaZError;   //!
   TBranch        *b_BeamSpot_z;   //!
   TBranch        *b_BeamSpot_zError;   //!
   TBranch        *b_nboostedTau;   //!
   TBranch        *b_boostedTau_idAntiEle2018;   //!
   TBranch        *b_boostedTau_idAntiMu;   //!
   TBranch        *b_boostedTau_idMVAnewDM2017v2;   //!
   TBranch        *b_boostedTau_idMVAoldDM2017v2;   //!
   TBranch        *b_boostedTau_jetIdx;   //!
   TBranch        *b_boostedTau_rawAntiEleCat2018;   //!
   TBranch        *b_boostedTau_charge;   //!
   TBranch        *b_boostedTau_decayMode;   //!
   TBranch        *b_boostedTau_chargedIso;   //!
   TBranch        *b_boostedTau_eta;   //!
   TBranch        *b_boostedTau_leadTkDeltaEta;   //!
   TBranch        *b_boostedTau_leadTkDeltaPhi;   //!
   TBranch        *b_boostedTau_leadTkPtOverTauPt;   //!
   TBranch        *b_boostedTau_mass;   //!
   TBranch        *b_boostedTau_neutralIso;   //!
   TBranch        *b_boostedTau_phi;   //!
   TBranch        *b_boostedTau_photonsOutsideSignalCone;   //!
   TBranch        *b_boostedTau_pt;   //!
   TBranch        *b_boostedTau_puCorr;   //!
   TBranch        *b_boostedTau_rawAntiEle2018;   //!
   TBranch        *b_boostedTau_rawIso;   //!
   TBranch        *b_boostedTau_rawIsodR03;   //!
   TBranch        *b_boostedTau_rawMVAnewDM2017v2;   //!
   TBranch        *b_boostedTau_rawMVAoldDM2017v2;   //!
   TBranch        *b_CaloMET_phi;   //!
   TBranch        *b_CaloMET_pt;   //!
   TBranch        *b_CaloMET_sumEt;   //!
   TBranch        *b_nCorrT1METJet;   //!
   TBranch        *b_CorrT1METJet_EmEF;   //!
   TBranch        *b_CorrT1METJet_area;   //!
   TBranch        *b_CorrT1METJet_eta;   //!
   TBranch        *b_CorrT1METJet_muonSubtrFactor;   //!
   TBranch        *b_CorrT1METJet_phi;   //!
   TBranch        *b_CorrT1METJet_rawPt;   //!
   TBranch        *b_DeepMETResolutionTune_phi;   //!
   TBranch        *b_DeepMETResolutionTune_pt;   //!
   TBranch        *b_DeepMETResponseTune_phi;   //!
   TBranch        *b_DeepMETResponseTune_pt;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isEB;   //!
   TBranch        *b_Electron_isEcalDriven;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_jetNDauCharged;   //!
   TBranch        *b_Electron_lostHits;   //!
   TBranch        *b_Electron_mvaIso_WP80;   //!
   TBranch        *b_Electron_mvaIso_WP90;   //!
   TBranch        *b_Electron_mvaIso_WPHZZ;   //!
   TBranch        *b_Electron_mvaNoIso_WP80;   //!
   TBranch        *b_Electron_mvaNoIso_WP90;   //!
   TBranch        *b_Electron_seedGain;   //!
   TBranch        *b_Electron_tightCharge;   //!
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_seediEtaOriX;   //!
   TBranch        *b_Electron_seediPhiOriY;   //!
   TBranch        *b_Electron_svIdx;   //!
   TBranch        *b_Electron_fsrPhotonIdx;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_Electron_vidNestedWPBitmapHEEP;   //!
   TBranch        *b_Electron_PreshowerEnergy;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_ecalEnergy;   //!
   TBranch        *b_Electron_ecalEnergyError;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_fbrem;   //!
   TBranch        *b_Electron_gsfTrketaMode;   //!
   TBranch        *b_Electron_gsfTrkpMode;   //!
   TBranch        *b_Electron_gsfTrkpModeErr;   //!
   TBranch        *b_Electron_gsfTrkphiMode;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetPtRelv2;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_Electron_mvaHZZIso;   //!
   TBranch        *b_Electron_mvaIso;   //!
   TBranch        *b_Electron_mvaNoIso;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_rawEnergy;   //!
   TBranch        *b_Electron_scEtOverPt;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_promptMVA;   //!
   TBranch        *b_nFatJet;   //!
   TBranch        *b_FatJet_jetId;   //!
   TBranch        *b_FatJet_nConstituents;   //!
   TBranch        *b_FatJet_chMultiplicity;   //!
   TBranch        *b_FatJet_neMultiplicity;   //!
   TBranch        *b_FatJet_subJetIdx1;   //!
   TBranch        *b_FatJet_subJetIdx2;   //!
   TBranch        *b_FatJet_electronIdx3SJ;   //!
   TBranch        *b_FatJet_muonIdx3SJ;   //!
   TBranch        *b_FatJet_area;   //!
   TBranch        *b_FatJet_chEmEF;   //!
   TBranch        *b_FatJet_chHEF;   //!
   TBranch        *b_FatJet_eta;   //!
   TBranch        *b_FatJet_mass;   //!
   TBranch        *b_FatJet_msoftdrop;   //!
   TBranch        *b_FatJet_muEF;   //!
   TBranch        *b_FatJet_n2b1;   //!
   TBranch        *b_FatJet_n3b1;   //!
   TBranch        *b_FatJet_neEmEF;   //!
   TBranch        *b_FatJet_neHEF;   //!
   TBranch        *b_FatJet_particleNetLegacy_QCD;   //!
   TBranch        *b_FatJet_particleNetLegacy_Xbb;   //!
   TBranch        *b_FatJet_particleNetLegacy_Xcc;   //!
   TBranch        *b_FatJet_particleNetLegacy_Xqq;   //!
   TBranch        *b_FatJet_particleNetLegacy_mass;   //!
   TBranch        *b_FatJet_particleNetWithMass_H4qvsQCD;   //!
   TBranch        *b_FatJet_particleNetWithMass_HbbvsQCD;   //!
   TBranch        *b_FatJet_particleNetWithMass_HccvsQCD;   //!
   TBranch        *b_FatJet_particleNetWithMass_QCD;   //!
   TBranch        *b_FatJet_particleNetWithMass_TvsQCD;   //!
   TBranch        *b_FatJet_particleNetWithMass_WvsQCD;   //!
   TBranch        *b_FatJet_particleNetWithMass_ZvsQCD;   //!
   TBranch        *b_FatJet_particleNet_QCD;   //!
   TBranch        *b_FatJet_particleNet_QCD0HF;   //!
   TBranch        *b_FatJet_particleNet_QCD1HF;   //!
   TBranch        *b_FatJet_particleNet_QCD2HF;   //!
   TBranch        *b_FatJet_particleNet_XbbVsQCD;   //!
   TBranch        *b_FatJet_particleNet_XccVsQCD;   //!
   TBranch        *b_FatJet_particleNet_XggVsQCD;   //!
   TBranch        *b_FatJet_particleNet_XqqVsQCD;   //!
   TBranch        *b_FatJet_particleNet_XteVsQCD;   //!
   TBranch        *b_FatJet_particleNet_XtmVsQCD;   //!
   TBranch        *b_FatJet_particleNet_XttVsQCD;   //!
   TBranch        *b_FatJet_particleNet_massCorr;   //!
   TBranch        *b_FatJet_phi;   //!
   TBranch        *b_FatJet_pt;   //!
   TBranch        *b_FatJet_rawFactor;   //!
   TBranch        *b_FatJet_tau1;   //!
   TBranch        *b_FatJet_tau2;   //!
   TBranch        *b_FatJet_tau3;   //!
   TBranch        *b_FatJet_tau4;   //!
   TBranch        *b_FatJet_lsf3;   //!
   TBranch        *b_nJetPFCand;   //!
   TBranch        *b_JetPFCand_pfCandIdx;   //!
   TBranch        *b_JetPFCand_jetIdx;   //!
   TBranch        *b_nFatJetPFCand;   //!
   TBranch        *b_FatJetPFCand_pfCandIdx;   //!
   TBranch        *b_FatJetPFCand_jetIdx;   //!
   TBranch        *b_nFsrPhoton;   //!
   TBranch        *b_FsrPhoton_electronIdx;   //!
   TBranch        *b_FsrPhoton_muonIdx;   //!
   TBranch        *b_FsrPhoton_dROverEt2;   //!
   TBranch        *b_FsrPhoton_eta;   //!
   TBranch        *b_FsrPhoton_phi;   //!
   TBranch        *b_FsrPhoton_pt;   //!
   TBranch        *b_FsrPhoton_relIso03;   //!
   TBranch        *b_nIsoTrack;   //!
   TBranch        *b_IsoTrack_isHighPurityTrack;   //!
   TBranch        *b_IsoTrack_isPFcand;   //!
   TBranch        *b_IsoTrack_isFromLostTrack;   //!
   TBranch        *b_IsoTrack_charge;   //!
   TBranch        *b_IsoTrack_fromPV;   //!
   TBranch        *b_IsoTrack_pdgId;   //!
   TBranch        *b_IsoTrack_dxy;   //!
   TBranch        *b_IsoTrack_dz;   //!
   TBranch        *b_IsoTrack_eta;   //!
   TBranch        *b_IsoTrack_pfRelIso03_all;   //!
   TBranch        *b_IsoTrack_pfRelIso03_chg;   //!
   TBranch        *b_IsoTrack_phi;   //!
   TBranch        *b_IsoTrack_pt;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_all;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_chg;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_chMultiplicity;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_nSVs;   //!
   TBranch        *b_Jet_neMultiplicity;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_svIdx1;   //!
   TBranch        *b_Jet_svIdx2;   //!
   TBranch        *b_Jet_hfadjacentEtaStripsSize;   //!
   TBranch        *b_Jet_hfcentralEtaStripSize;   //!
   TBranch        *b_Jet_PNetRegPtRawCorr;   //!
   TBranch        *b_Jet_PNetRegPtRawCorrNeutrino;   //!
   TBranch        *b_Jet_PNetRegPtRawRes;   //!
   TBranch        *b_Jet_UParTAK4RegPtRawCorr;   //!
   TBranch        *b_Jet_UParTAK4RegPtRawCorrNeutrino;   //!
   TBranch        *b_Jet_UParTAK4RegPtRawRes;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_btagDeepFlavCvB;   //!
   TBranch        *b_Jet_btagDeepFlavCvL;   //!
   TBranch        *b_Jet_btagDeepFlavQG;   //!
   TBranch        *b_Jet_btagPNetB;   //!
   TBranch        *b_Jet_btagPNetCvB;   //!
   TBranch        *b_Jet_btagPNetCvL;   //!
   TBranch        *b_Jet_btagPNetCvNotB;   //!
   TBranch        *b_Jet_btagPNetQvG;   //!
   TBranch        *b_Jet_btagPNetTauVJet;   //!
   TBranch        *b_Jet_btagUParTAK4B;   //!
   TBranch        *b_Jet_btagUParTAK4CvB;   //!
   TBranch        *b_Jet_btagUParTAK4CvL;   //!
   TBranch        *b_Jet_btagUParTAK4CvNotB;   //!
   TBranch        *b_Jet_btagUParTAK4QvG;   //!
   TBranch        *b_Jet_btagUParTAK4TauVJet;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_hfEmEF;   //!
   TBranch        *b_Jet_hfHEF;   //!
   TBranch        *b_Jet_hfsigmaEtaEta;   //!
   TBranch        *b_Jet_hfsigmaPhiPhi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_nLowPtElectron;   //!
   TBranch        *b_LowPtElectron_convVeto;   //!
   TBranch        *b_LowPtElectron_convWP;   //!
   TBranch        *b_LowPtElectron_lostHits;   //!
   TBranch        *b_LowPtElectron_electronIdx;   //!
   TBranch        *b_LowPtElectron_photonIdx;   //!
   TBranch        *b_LowPtElectron_charge;   //!
   TBranch        *b_LowPtElectron_pdgId;   //!
   TBranch        *b_LowPtElectron_ID;   //!
   TBranch        *b_LowPtElectron_convVtxRadius;   //!
   TBranch        *b_LowPtElectron_deltaEtaSC;   //!
   TBranch        *b_LowPtElectron_dxy;   //!
   TBranch        *b_LowPtElectron_dxyErr;   //!
   TBranch        *b_LowPtElectron_dz;   //!
   TBranch        *b_LowPtElectron_dzErr;   //!
   TBranch        *b_LowPtElectron_eInvMinusPInv;   //!
   TBranch        *b_LowPtElectron_energyErr;   //!
   TBranch        *b_LowPtElectron_eta;   //!
   TBranch        *b_LowPtElectron_hoe;   //!
   TBranch        *b_LowPtElectron_mass;   //!
   TBranch        *b_LowPtElectron_miniPFRelIso_all;   //!
   TBranch        *b_LowPtElectron_miniPFRelIso_chg;   //!
   TBranch        *b_LowPtElectron_phi;   //!
   TBranch        *b_LowPtElectron_pt;   //!
   TBranch        *b_LowPtElectron_ptbiased;   //!
   TBranch        *b_LowPtElectron_r9;   //!
   TBranch        *b_LowPtElectron_scEtOverPt;   //!
   TBranch        *b_LowPtElectron_sieie;   //!
   TBranch        *b_LowPtElectron_unbiased;   //!
   TBranch        *b_nProton_multiRP;   //!
   TBranch        *b_Proton_multiRP_arm;   //!
   TBranch        *b_Proton_multiRP_t;   //!
   TBranch        *b_Proton_multiRP_thetaX;   //!
   TBranch        *b_Proton_multiRP_thetaY;   //!
   TBranch        *b_Proton_multiRP_time;   //!
   TBranch        *b_Proton_multiRP_timeUnc;   //!
   TBranch        *b_Proton_multiRP_xi;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_highPurity;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isStandalone;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_jetNDauCharged;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaMuID_WP;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_puppiIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_svIdx;   //!
   TBranch        *b_Muon_fsrPhotonIdx;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dxybs;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetPtRelv2;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_mvaMuID;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_softMva;   //!
   TBranch        *b_Muon_softMvaRun3;   //!
   TBranch        *b_Muon_tkRelIso;   //!
   TBranch        *b_Muon_tunepRelPt;   //!
   TBranch        *b_Muon_bsConstrainedChi2;   //!
   TBranch        *b_Muon_bsConstrainedPt;   //!
   TBranch        *b_Muon_bsConstrainedPtErr;   //!
   TBranch        *b_Muon_mvaLowPt;   //!
   TBranch        *b_Muon_promptMVA;   //!
   TBranch        *b_nPFCandToPFClusterHCAL;   //!
   TBranch        *b_PFCandToPFClusterHCAL_PFCandIdx;   //!
   TBranch        *b_PFCandToPFClusterHCAL_PFClusterHCALIdx;   //!
   TBranch        *b_nPFClusterECAL;   //!
   TBranch        *b_PFClusterECAL_layer;   //!
   TBranch        *b_PFClusterECAL_nhits;   //!
   TBranch        *b_PFClusterECAL_seedhit_detId;   //!
   TBranch        *b_PFClusterECAL_pt;   //!
   TBranch        *b_PFClusterECAL_energy;   //!
   TBranch        *b_PFClusterECAL_correctedEnergy;   //!
   TBranch        *b_PFClusterECAL_eta;   //!
   TBranch        *b_PFClusterECAL_phi;   //!
   TBranch        *b_nPFClusterHCAL;   //!
   TBranch        *b_PFClusterHCAL_layer;   //!
   TBranch        *b_PFClusterHCAL_nhits;   //!
   TBranch        *b_PFClusterHCAL_seedhit_detId;   //!
   TBranch        *b_PFClusterHCAL_pt;   //!
   TBranch        *b_PFClusterHCAL_energy;   //!
   TBranch        *b_PFClusterHCAL_eta;   //!
   TBranch        *b_PFClusterHCAL_phi;   //!
   TBranch        *b_nPFClusterHCALToPFRecHitHBHE;   //!
   TBranch        *b_PFClusterHCALToPFRecHitHBHE_PFClusterHCALIdx;   //!
   TBranch        *b_PFClusterHCALToPFRecHitHBHE_PFRecHitHBHEIdx;   //!
   TBranch        *b_PFClusterHCALToPFRecHitHBHE_fraction;   //!
   TBranch        *b_nPFClusterPS;   //!
   TBranch        *b_PFClusterPS_layer;   //!
   TBranch        *b_PFClusterPS_nhits;   //!
   TBranch        *b_PFClusterPS_seedhit_detId;   //!
   TBranch        *b_PFClusterPS_pt;   //!
   TBranch        *b_PFClusterPS_energy;   //!
   TBranch        *b_PFClusterPS_correctedEnergy;   //!
   TBranch        *b_PFClusterPS_eta;   //!
   TBranch        *b_PFClusterPS_phi;   //!
   TBranch        *b_nPFRecHitHBHE;   //!
   TBranch        *b_PFRecHitHBHE_ieta;   //!
   TBranch        *b_PFRecHitHBHE_iphi;   //!
   TBranch        *b_PFRecHitHBHE_depth;   //!
   TBranch        *b_PFRecHitHBHE_detId;   //!
   TBranch        *b_PFRecHitHBHE_energy;   //!
   TBranch        *b_nPFCand;   //!
   TBranch        *b_PFCand_isIsolatedChargedHadron;   //!
   TBranch        *b_PFCand_trkHighPurity;   //!
   TBranch        *b_PFCand_charge;   //!
   TBranch        *b_PFCand_fromPV0;   //!
   TBranch        *b_PFCand_lostInnerHits;   //!
   TBranch        *b_PFCand_lostOuterHits;   //!
   TBranch        *b_PFCand_nHits;   //!
   TBranch        *b_PFCand_nPixelHits;   //!
   TBranch        *b_PFCand_pdgId;   //!
   TBranch        *b_PFCand_pvAssocQuality;   //!
   TBranch        *b_PFCand_trkAlgo;   //!
   TBranch        *b_PFCand_trkQuality;   //!
   TBranch        *b_PFCand_vertexRef;   //!
   TBranch        *b_PFCand_d0;   //!
   TBranch        *b_PFCand_d0Err;   //!
   TBranch        *b_PFCand_dz;   //!
   TBranch        *b_PFCand_dzErr;   //!
   TBranch        *b_PFCand_energy;   //!
   TBranch        *b_PFCand_eta;   //!
   TBranch        *b_PFCand_mass;   //!
   TBranch        *b_PFCand_phi;   //!
   TBranch        *b_PFCand_pt;   //!
   TBranch        *b_PFCand_puppiWeight;   //!
   TBranch        *b_PFCand_puppiWeightNoLep;   //!
   TBranch        *b_PFCand_trkChi2;   //!
   TBranch        *b_PFCand_trkEta;   //!
   TBranch        *b_PFCand_trkP;   //!
   TBranch        *b_PFCand_trkPhi;   //!
   TBranch        *b_PFCand_trkPt;   //!
   TBranch        *b_PFCand_vtxChi2;   //!
   TBranch        *b_PFMET_covXX;   //!
   TBranch        *b_PFMET_covXY;   //!
   TBranch        *b_PFMET_covYY;   //!
   TBranch        *b_PFMET_phi;   //!
   TBranch        *b_PFMET_phiUnclusteredDown;   //!
   TBranch        *b_PFMET_phiUnclusteredUp;   //!
   TBranch        *b_PFMET_pt;   //!
   TBranch        *b_PFMET_ptUnclusteredDown;   //!
   TBranch        *b_PFMET_ptUnclusteredUp;   //!
   TBranch        *b_PFMET_significance;   //!
   TBranch        *b_PFMET_sumEt;   //!
   TBranch        *b_PFMET_sumPtUnclustered;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_cutBased;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_hasConversionTracks;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_seedGain;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
   TBranch        *b_Photon_seediEtaOriX;   //!
   TBranch        *b_Photon_seediPhiOriY;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_ecalPFClusterIso;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_energyRaw;   //!
   TBranch        *b_Photon_esEffSigmaRR;   //!
   TBranch        *b_Photon_esEnergyOverRawE;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_etaWidth;   //!
   TBranch        *b_Photon_haloTaggerMVAVal;   //!
   TBranch        *b_Photon_hcalPFClusterIso;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_hoe_PUcorr;   //!
   TBranch        *b_Photon_hoe_Tower;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_pfChargedIso;   //!
   TBranch        *b_Photon_pfChargedIsoPFPV;   //!
   TBranch        *b_Photon_pfChargedIsoWorstVtx;   //!
   TBranch        *b_Photon_pfPhoIso03;   //!
   TBranch        *b_Photon_pfRelIso03_all_quadratic;   //!
   TBranch        *b_Photon_pfRelIso03_chg_quadratic;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_phiWidth;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_s4;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_sieip;   //!
   TBranch        *b_Photon_sipip;   //!
   TBranch        *b_Photon_superclusterEta;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_x_calo;   //!
   TBranch        *b_Photon_y_calo;   //!
   TBranch        *b_Photon_z_calo;   //!
   TBranch        *b_nPPSLocalTrack;   //!
   TBranch        *b_PPSLocalTrack_multiRPProtonIdx;   //!
   TBranch        *b_PPSLocalTrack_singleRPProtonIdx;   //!
   TBranch        *b_PPSLocalTrack_decRPId;   //!
   TBranch        *b_PPSLocalTrack_rpType;   //!
   TBranch        *b_PPSLocalTrack_x;   //!
   TBranch        *b_PPSLocalTrack_y;   //!
   TBranch        *b_PPSLocalTrack_time;   //!
   TBranch        *b_PPSLocalTrack_timeUnc;   //!
   TBranch        *b_PuppiMET_covXX;   //!
   TBranch        *b_PuppiMET_covXY;   //!
   TBranch        *b_PuppiMET_covYY;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_phiUnclusteredDown;   //!
   TBranch        *b_PuppiMET_phiUnclusteredUp;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_ptUnclusteredDown;   //!
   TBranch        *b_PuppiMET_ptUnclusteredUp;   //!
   TBranch        *b_PuppiMET_significance;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_PuppiMET_sumPtUnclustered;   //!
   TBranch        *b_nPVBS;   //!
   TBranch        *b_PVBS_chi2;   //!
   TBranch        *b_PVBS_cov00;   //!
   TBranch        *b_PVBS_cov10;   //!
   TBranch        *b_PVBS_cov11;   //!
   TBranch        *b_PVBS_cov20;   //!
   TBranch        *b_PVBS_cov21;   //!
   TBranch        *b_PVBS_cov22;   //!
   TBranch        *b_PVBS_x;   //!
   TBranch        *b_PVBS_y;   //!
   TBranch        *b_PVBS_z;   //!
   TBranch        *b_RawPFMET_phi;   //!
   TBranch        *b_RawPFMET_pt;   //!
   TBranch        *b_RawPFMET_sumEt;   //!
   TBranch        *b_RawPuppiMET_phi;   //!
   TBranch        *b_RawPuppiMET_pt;   //!
   TBranch        *b_RawPuppiMET_sumEt;   //!
   TBranch        *b_Rho_fixedGridRhoAll;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetAll;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_Rho_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nSoftActivityJet;   //!
   TBranch        *b_SoftActivityJet_eta;   //!
   TBranch        *b_SoftActivityJet_phi;   //!
   TBranch        *b_SoftActivityJet_pt;   //!
   TBranch        *b_SoftActivityJetNjets10;   //!
   TBranch        *b_SoftActivityJetNjets2;   //!
   TBranch        *b_SoftActivityJetNjets5;   //!
   TBranch        *b_SoftActivityJetHT;   //!
   TBranch        *b_SoftActivityJetHT10;   //!
   TBranch        *b_SoftActivityJetHT2;   //!
   TBranch        *b_SoftActivityJetHT5;   //!
   TBranch        *b_nProton_singleRP;   //!
   TBranch        *b_Proton_singleRP_decRPId;   //!
   TBranch        *b_Proton_singleRP_thetaY;   //!
   TBranch        *b_Proton_singleRP_xi;   //!
   TBranch        *b_nSubJet;   //!
   TBranch        *b_SubJet_area;   //!
   TBranch        *b_SubJet_btagDeepB;   //!
   TBranch        *b_SubJet_eta;   //!
   TBranch        *b_SubJet_mass;   //!
   TBranch        *b_SubJet_n2b1;   //!
   TBranch        *b_SubJet_n3b1;   //!
   TBranch        *b_SubJet_phi;   //!
   TBranch        *b_SubJet_pt;   //!
   TBranch        *b_SubJet_rawFactor;   //!
   TBranch        *b_SubJet_tau1;   //!
   TBranch        *b_SubJet_tau2;   //!
   TBranch        *b_SubJet_tau3;   //!
   TBranch        *b_SubJet_tau4;   //!
   TBranch        *b_nTauProd;   //!
   TBranch        *b_TauProd_tauIdx;   //!
   TBranch        *b_TauProd_pdgId;   //!
   TBranch        *b_TauProd_eta;   //!
   TBranch        *b_TauProd_phi;   //!
   TBranch        *b_TauProd_pt;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_idAntiEleDeadECal;   //!
   TBranch        *b_Tau_idAntiMu;   //!
   TBranch        *b_Tau_idDecayModeNewDMs;   //!
   TBranch        *b_Tau_idDecayModeOldDMs;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_idDeepTau2017v2p1VSmu;   //!
   TBranch        *b_Tau_idDeepTau2018v2p5VSe;   //!
   TBranch        *b_Tau_idDeepTau2018v2p5VSjet;   //!
   TBranch        *b_Tau_idDeepTau2018v2p5VSmu;   //!
   TBranch        *b_Tau_nSVs;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayModePNet;   //!
   TBranch        *b_Tau_decayModeUParT;   //!
   TBranch        *b_Tau_eleIdx;   //!
   TBranch        *b_Tau_jetIdx;   //!
   TBranch        *b_Tau_muIdx;   //!
   TBranch        *b_Tau_svIdx1;   //!
   TBranch        *b_Tau_svIdx2;   //!
   TBranch        *b_Tau_chargedIso;   //!
   TBranch        *b_Tau_dxy;   //!
   TBranch        *b_Tau_dz;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_leadTkDeltaEta;   //!
   TBranch        *b_Tau_leadTkDeltaPhi;   //!
   TBranch        *b_Tau_leadTkPtOverTauPt;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_neutralIso;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_photonsOutsideSignalCone;   //!
   TBranch        *b_Tau_probDM0PNet;   //!
   TBranch        *b_Tau_probDM0UParT;   //!
   TBranch        *b_Tau_probDM10PNet;   //!
   TBranch        *b_Tau_probDM10UParT;   //!
   TBranch        *b_Tau_probDM11PNet;   //!
   TBranch        *b_Tau_probDM11UParT;   //!
   TBranch        *b_Tau_probDM1PNet;   //!
   TBranch        *b_Tau_probDM1UParT;   //!
   TBranch        *b_Tau_probDM2PNet;   //!
   TBranch        *b_Tau_probDM2UParT;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_ptCorrPNet;   //!
   TBranch        *b_Tau_ptCorrUParT;   //!
   TBranch        *b_Tau_puCorr;   //!
   TBranch        *b_Tau_qConfPNet;   //!
   TBranch        *b_Tau_qConfUParT;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSe;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSjet;   //!
   TBranch        *b_Tau_rawDeepTau2017v2p1VSmu;   //!
   TBranch        *b_Tau_rawDeepTau2018v2p5VSe;   //!
   TBranch        *b_Tau_rawDeepTau2018v2p5VSjet;   //!
   TBranch        *b_Tau_rawDeepTau2018v2p5VSmu;   //!
   TBranch        *b_Tau_rawIso;   //!
   TBranch        *b_Tau_rawIsodR03;   //!
   TBranch        *b_Tau_rawPNetVSe;   //!
   TBranch        *b_Tau_rawPNetVSjet;   //!
   TBranch        *b_Tau_rawPNetVSmu;   //!
   TBranch        *b_Tau_rawUParTVSe;   //!
   TBranch        *b_Tau_rawUParTVSjet;   //!
   TBranch        *b_Tau_rawUParTVSmu;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrkMET_phi;   //!
   TBranch        *b_TrkMET_pt;   //!
   TBranch        *b_TrkMET_sumEt;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_OtherPV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_sumpt2;   //!
   TBranch        *b_PV_sumpx;   //!
   TBranch        *b_PV_sumpy;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_charge;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_dxy;   //!
   TBranch        *b_SV_dxySig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_Electron_IPx;   //!
   TBranch        *b_Electron_IPy;   //!
   TBranch        *b_Electron_IPz;   //!
   TBranch        *b_Electron_ipLengthSig;   //!
   TBranch        *b_Muon_IPx;   //!
   TBranch        *b_Muon_IPy;   //!
   TBranch        *b_Muon_IPz;   //!
   TBranch        *b_Muon_ipLengthSig;   //!
   TBranch        *b_PFCand_fromPVvertexRef;   //!
   TBranch        *b_PFCand_blockIdx;   //!
   TBranch        *b_PFCand_nPFTrack;   //!
   TBranch        *b_PFCand_nPFClusterHCAL;   //!
   TBranch        *b_PFCand_nPFClusterECAL;   //!
   TBranch        *b_PFCand_nPFClusterPS;   //!
   TBranch        *b_PFCand_nPFTrackInBlock;   //!
   TBranch        *b_PFCand_nPFClusterHCALInBlock;   //!
   TBranch        *b_PFCand_nPFClusterECALInBlock;   //!
   TBranch        *b_PFCand_nPFClusterPSInBlock;   //!
   TBranch        *b_PFCand_hcalDepthEnergyFraction1;   //!
   TBranch        *b_PFCand_hcalDepthEnergyFraction2;   //!
   TBranch        *b_PFCand_hcalDepthEnergyFraction3;   //!
   TBranch        *b_PFCand_hcalDepthEnergyFraction4;   //!
   TBranch        *b_PFCand_hcalDepthEnergyFraction5;   //!
   TBranch        *b_PFCand_hcalDepthEnergyFraction6;   //!
   TBranch        *b_PFCand_hcalDepthEnergyFraction7;   //!
   TBranch        *b_PFCand_ecalEnergy;   //!
   TBranch        *b_PFCand_rawEcalEnergy;   //!
   TBranch        *b_PFCand_hcalEnergy;   //!
   TBranch        *b_PFCand_rawHcalEnergy;   //!
   TBranch        *b_PFCand_hoEnergy;   //!
   TBranch        *b_PFCand_rawHoEnergy;   //!
   TBranch        *b_SV_ntracks;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_Tau_hasRefitSV;   //!
   TBranch        *b_Tau_IPx;   //!
   TBranch        *b_Tau_IPy;   //!
   TBranch        *b_Tau_IPz;   //!
   TBranch        *b_Tau_ipLengthSig;   //!
   TBranch        *b_Tau_refitSVchi2;   //!
   TBranch        *b_Tau_refitSVcov00;   //!
   TBranch        *b_Tau_refitSVcov10;   //!
   TBranch        *b_Tau_refitSVcov11;   //!
   TBranch        *b_Tau_refitSVcov20;   //!
   TBranch        *b_Tau_refitSVcov21;   //!
   TBranch        *b_Tau_refitSVcov22;   //!
   TBranch        *b_Tau_refitSVx;   //!
   TBranch        *b_Tau_refitSVy;   //!
   TBranch        *b_Tau_refitSVz;   //!
   TBranch        *b_L1Reco_step;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_L1_AXO_Loose;   //!
   TBranch        *b_L1_AXO_Nominal;   //!
   TBranch        *b_L1_AXO_Tight;   //!
   TBranch        *b_L1_AXO_VLoose;   //!
   TBranch        *b_L1_AXO_VTight;   //!
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_BPTX_AND_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B2_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref2_VME;   //!
   TBranch        *b_L1_BPTX_NotOR_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_RefAND_VME;   //!
   TBranch        *b_L1_BptxMinus;   //!
   TBranch        *b_L1_BptxOR;   //!
   TBranch        *b_L1_BptxPlus;   //!
   TBranch        *b_L1_BptxXOR;   //!
   TBranch        *b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_L1_CICADA_Loose;   //!
   TBranch        *b_L1_CICADA_Medium;   //!
   TBranch        *b_L1_CICADA_Tight;   //!
   TBranch        *b_L1_CICADA_VLoose;   //!
   TBranch        *b_L1_CICADA_VTight;   //!
   TBranch        *b_L1_DoubleEG11_er1p2_dR_Max0p6;   //!
   TBranch        *b_L1_DoubleEG8er2p5_HTT280er;   //!
   TBranch        *b_L1_DoubleEG8er2p5_HTT300er;   //!
   TBranch        *b_L1_DoubleEG8er2p5_HTT320er;   //!
   TBranch        *b_L1_DoubleEG_15_10_er2p5;   //!
   TBranch        *b_L1_DoubleEG_20_10_er2p5;   //!
   TBranch        *b_L1_DoubleEG_22_10_er2p5;   //!
   TBranch        *b_L1_DoubleEG_25_12_er2p5;   //!
   TBranch        *b_L1_DoubleEG_25_14_er2p5;   //!
   TBranch        *b_L1_DoubleEG_27_14_er2p5;   //!
   TBranch        *b_L1_DoubleEG_LooseIso16_LooseIso12_er1p5;   //!
   TBranch        *b_L1_DoubleEG_LooseIso18_LooseIso12_er1p5;   //!
   TBranch        *b_L1_DoubleEG_LooseIso20_LooseIso12_er1p5;   //!
   TBranch        *b_L1_DoubleEG_LooseIso22_12_er2p5;   //!
   TBranch        *b_L1_DoubleEG_LooseIso22_LooseIso12_er1p5;   //!
   TBranch        *b_L1_DoubleEG_LooseIso25_12_er2p5;   //!
   TBranch        *b_L1_DoubleEG_LooseIso25_LooseIso12_er1p5;   //!
   TBranch        *b_L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5;   //!
   TBranch        *b_L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5;   //!
   TBranch        *b_L1_DoubleIsoTau28er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau30er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau32er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau32er2p1_Mass_Max80;   //!
   TBranch        *b_L1_DoubleIsoTau34er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau35er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau36er2p1;   //!
   TBranch        *b_L1_DoubleJet100er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet100er2p5;   //!
   TBranch        *b_L1_DoubleJet112er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet120er2p5;   //!
   TBranch        *b_L1_DoubleJet120er2p5_Mu3_dR_Max0p8;   //!
   TBranch        *b_L1_DoubleJet150er2p5;   //!
   TBranch        *b_L1_DoubleJet16er2p5_Mu3_dR_Max0p4;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet35er2p5_Mu3_dR_Max0p4;   //!
   TBranch        *b_L1_DoubleJet40er2p5;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2;   //!
   TBranch        *b_L1_DoubleJet60er2p5_Mu3_dR_Max0p4;   //!
   TBranch        *b_L1_DoubleJet80er2p5_Mu3_dR_Max0p4;   //!
   TBranch        *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min800;   //!
   TBranch        *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min850;   //!
   TBranch        *b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50;   //!
   TBranch        *b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50;   //!
   TBranch        *b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65;   //!
   TBranch        *b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65;   //!
   TBranch        *b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ;   //!
   TBranch        *b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ;   //!
   TBranch        *b_L1_DoubleLLPJet40;   //!
   TBranch        *b_L1_DoubleLooseIsoEG22er2p1;   //!
   TBranch        *b_L1_DoubleLooseIsoEG24er2p1;   //!
   TBranch        *b_L1_DoubleMu0;   //!
   TBranch        *b_L1_DoubleMu0_Mass_Min1;   //!
   TBranch        *b_L1_DoubleMu0_OQ;   //!
   TBranch        *b_L1_DoubleMu0_SQ;   //!
   TBranch        *b_L1_DoubleMu0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0_Upt15_Upt7;   //!
   TBranch        *b_L1_DoubleMu0_Upt5_Upt5;   //!
   TBranch        *b_L1_DoubleMu0_Upt6_IP_Min1_Upt4;   //!
   TBranch        *b_L1_DoubleMu0_Upt6_SQ_er2p0;   //!
   TBranch        *b_L1_DoubleMu0_Upt7_SQ_er2p0;   //!
   TBranch        *b_L1_DoubleMu0_Upt8_SQ_er2p0;   //!
   TBranch        *b_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8;   //!
   TBranch        *b_L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2;   //!
   TBranch        *b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleMu18er2p1_SQ;   //!
   TBranch        *b_L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF30_HTT60er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF40_HTT60er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_HTT60er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT220er;   //!
   TBranch        *b_L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6;   //!
   TBranch        *b_L1_DoubleMu4_SQ_EG9er2p5;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7;   //!
   TBranch        *b_L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20;   //!
   TBranch        *b_L1_DoubleMu5_SQ_EG9er2p5;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS_dR_Max1p6;   //!
   TBranch        *b_L1_DoubleMu6_Upt6_SQ_er2p0;   //!
   TBranch        *b_L1_DoubleMu7_Upt7_SQ_er2p0;   //!
   TBranch        *b_L1_DoubleMu8_SQ;   //!
   TBranch        *b_L1_DoubleMu8_Upt8_SQ_er2p0;   //!
   TBranch        *b_L1_DoubleMu9_SQ;   //!
   TBranch        *b_L1_DoubleMu_12_5;   //!
   TBranch        *b_L1_DoubleMu_15_5_SQ;   //!
   TBranch        *b_L1_DoubleMu_15_7;   //!
   TBranch        *b_L1_DoubleMu_15_7_Mass_Min1;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ;   //!
   TBranch        *b_L1_DoubleTau70er2p1;   //!
   TBranch        *b_L1_ETM120;   //!
   TBranch        *b_L1_ETM150;   //!
   TBranch        *b_L1_ETMHF100;   //!
   TBranch        *b_L1_ETMHF100_HTT60er;   //!
   TBranch        *b_L1_ETMHF110;   //!
   TBranch        *b_L1_ETMHF110_HTT60er;   //!
   TBranch        *b_L1_ETMHF120;   //!
   TBranch        *b_L1_ETMHF120_HTT60er;   //!
   TBranch        *b_L1_ETMHF130;   //!
   TBranch        *b_L1_ETMHF130_HTT60er;   //!
   TBranch        *b_L1_ETMHF140;   //!
   TBranch        *b_L1_ETMHF150;   //!
   TBranch        *b_L1_ETMHF70;   //!
   TBranch        *b_L1_ETMHF70_HTT60er;   //!
   TBranch        *b_L1_ETMHF80;   //!
   TBranch        *b_L1_ETMHF80_HTT60er;   //!
   TBranch        *b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1;   //!
   TBranch        *b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6;   //!
   TBranch        *b_L1_ETMHF90;   //!
   TBranch        *b_L1_ETMHF90_HTT60er;   //!
   TBranch        *b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1;   //!
   TBranch        *b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6;   //!
   TBranch        *b_L1_ETT2000;   //!
   TBranch        *b_L1_FirstBunchAfterTrain;   //!
   TBranch        *b_L1_FirstBunchBeforeTrain;   //!
   TBranch        *b_L1_FirstBunchInTrain;   //!
   TBranch        *b_L1_FirstCollisionInOrbit;   //!
   TBranch        *b_L1_FirstCollisionInTrain;   //!
   TBranch        *b_L1_HCAL_LaserMon_Trig;   //!
   TBranch        *b_L1_HCAL_LaserMon_Veto;   //!
   TBranch        *b_L1_HTMHF100;   //!
   TBranch        *b_L1_HTMHF120;   //!
   TBranch        *b_L1_HTMHF125;   //!
   TBranch        *b_L1_HTMHF130;   //!
   TBranch        *b_L1_HTMHF150;   //!
   TBranch        *b_L1_HTT120_SingleLLPJet40;   //!
   TBranch        *b_L1_HTT120er;   //!
   TBranch        *b_L1_HTT160_SingleLLPJet50;   //!
   TBranch        *b_L1_HTT160er;   //!
   TBranch        *b_L1_HTT200_SingleLLPJet60;   //!
   TBranch        *b_L1_HTT200er;   //!
   TBranch        *b_L1_HTT240_SingleLLPJet70;   //!
   TBranch        *b_L1_HTT255er;   //!
   TBranch        *b_L1_HTT280er;   //!
   TBranch        *b_L1_HTT280er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT320er;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p5;   //!
   TBranch        *b_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;   //!
   TBranch        *b_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;   //!
   TBranch        *b_L1_HTT360er;   //!
   TBranch        *b_L1_HTT400er;   //!
   TBranch        *b_L1_HTT450er;   //!
   TBranch        *b_L1_IsoEG32er2p5_Mt40;   //!
   TBranch        *b_L1_IsoTau52er2p1_QuadJet36er2p5;   //!
   TBranch        *b_L1_IsolatedBunch;   //!
   TBranch        *b_L1_LastBunchInTrain;   //!
   TBranch        *b_L1_LastCollisionInTrain;   //!
   TBranch        *b_L1_LooseIsoEG14er2p5_HTT200er;   //!
   TBranch        *b_L1_LooseIsoEG16er2p5_HTT200er;   //!
   TBranch        *b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3;   //!
   TBranch        *b_L1_MinimumBiasHF0;   //!
   TBranch        *b_L1_MinimumBiasHF0_AND_BptxAND;   //!
   TBranch        *b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu12_HTT150er;   //!
   TBranch        *b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu14_HTT150er;   //!
   TBranch        *b_L1_Mu18er2p1_Tau24er2p1;   //!
   TBranch        *b_L1_Mu18er2p1_Tau26er2p1;   //!
   TBranch        *b_L1_Mu18er2p1_Tau26er2p1_Jet55;   //!
   TBranch        *b_L1_Mu18er2p1_Tau26er2p1_Jet70;   //!
   TBranch        *b_L1_Mu20_EG10er2p5;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau30er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau32er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau34er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau40er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_Tau70er2p1;   //!
   TBranch        *b_L1_Mu3_Jet120er2p5_dR_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet16er2p5_dR_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet30er2p5;   //!
   TBranch        *b_L1_Mu3_Jet60er2p5_dR_Max0p4;   //!
   TBranch        *b_L1_Mu3er1p5_Jet100er2p5_ETMHF30;   //!
   TBranch        *b_L1_Mu3er1p5_Jet100er2p5_ETMHF40;   //!
   TBranch        *b_L1_Mu3er1p5_Jet100er2p5_ETMHF50;   //!
   TBranch        *b_L1_Mu5_EG23er2p5;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG20er2p5;   //!
   TBranch        *b_L1_Mu6_DoubleEG10er2p5;   //!
   TBranch        *b_L1_Mu6_DoubleEG12er2p5;   //!
   TBranch        *b_L1_Mu6_DoubleEG15er2p5;   //!
   TBranch        *b_L1_Mu6_DoubleEG17er2p5;   //!
   TBranch        *b_L1_Mu6_HTT240er;   //!
   TBranch        *b_L1_Mu6_HTT250er;   //!
   TBranch        *b_L1_Mu7_EG20er2p5;   //!
   TBranch        *b_L1_Mu7_EG23er2p5;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG20er2p5;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG23er2p5;   //!
   TBranch        *b_L1_NotBptxOR;   //!
   TBranch        *b_L1_QuadJet60er2p5;   //!
   TBranch        *b_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0;   //!
   TBranch        *b_L1_QuadMu0;   //!
   TBranch        *b_L1_QuadMu0_OQ;   //!
   TBranch        *b_L1_QuadMu0_SQ;   //!
   TBranch        *b_L1_SecondBunchInTrain;   //!
   TBranch        *b_L1_SecondLastBunchInTrain;   //!
   TBranch        *b_L1_SingleEG10er2p5;   //!
   TBranch        *b_L1_SingleEG15er2p5;   //!
   TBranch        *b_L1_SingleEG26er2p5;   //!
   TBranch        *b_L1_SingleEG28_FWD2p5;   //!
   TBranch        *b_L1_SingleEG28er1p5;   //!
   TBranch        *b_L1_SingleEG28er2p1;   //!
   TBranch        *b_L1_SingleEG28er2p5;   //!
   TBranch        *b_L1_SingleEG34er2p5;   //!
   TBranch        *b_L1_SingleEG36er2p5;   //!
   TBranch        *b_L1_SingleEG38er2p5;   //!
   TBranch        *b_L1_SingleEG40er2p5;   //!
   TBranch        *b_L1_SingleEG42er2p5;   //!
   TBranch        *b_L1_SingleEG45er2p5;   //!
   TBranch        *b_L1_SingleEG50;   //!
   TBranch        *b_L1_SingleEG60;   //!
   TBranch        *b_L1_SingleEG8er2p5;   //!
   TBranch        *b_L1_SingleIsoEG24er2p1;   //!
   TBranch        *b_L1_SingleIsoEG26er2p1;   //!
   TBranch        *b_L1_SingleIsoEG26er2p5;   //!
   TBranch        *b_L1_SingleIsoEG28_FWD2p5;   //!
   TBranch        *b_L1_SingleIsoEG28er1p5;   //!
   TBranch        *b_L1_SingleIsoEG28er2p1;   //!
   TBranch        *b_L1_SingleIsoEG28er2p5;   //!
   TBranch        *b_L1_SingleIsoEG30er2p1;   //!
   TBranch        *b_L1_SingleIsoEG30er2p5;   //!
   TBranch        *b_L1_SingleIsoEG32er2p1;   //!
   TBranch        *b_L1_SingleIsoEG32er2p5;   //!
   TBranch        *b_L1_SingleIsoEG34er2p5;   //!
   TBranch        *b_L1_SingleJet10erHE;   //!
   TBranch        *b_L1_SingleJet120;   //!
   TBranch        *b_L1_SingleJet120_FWD2p5;   //!
   TBranch        *b_L1_SingleJet120_FWD3p0;   //!
   TBranch        *b_L1_SingleJet120er1p3;   //!
   TBranch        *b_L1_SingleJet120er2p5;   //!
   TBranch        *b_L1_SingleJet12erHE;   //!
   TBranch        *b_L1_SingleJet140er2p5;   //!
   TBranch        *b_L1_SingleJet160er2p5;   //!
   TBranch        *b_L1_SingleJet180;   //!
   TBranch        *b_L1_SingleJet180er2p5;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleJet20er2p5_NotBptxOR;   //!
   TBranch        *b_L1_SingleJet20er2p5_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet35;   //!
   TBranch        *b_L1_SingleJet35_FWD2p5;   //!
   TBranch        *b_L1_SingleJet35_FWD3p0;   //!
   TBranch        *b_L1_SingleJet35er1p3;   //!
   TBranch        *b_L1_SingleJet35er2p5;   //!
   TBranch        *b_L1_SingleJet43er2p5_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet46er2p5_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet60;   //!
   TBranch        *b_L1_SingleJet60_FWD2p5;   //!
   TBranch        *b_L1_SingleJet8erHE;   //!
   TBranch        *b_L1_SingleJet90;   //!
   TBranch        *b_L1_SingleJet90_FWD2p5;   //!
   TBranch        *b_L1_SingleLooseIsoEG26er1p5;   //!
   TBranch        *b_L1_SingleLooseIsoEG26er2p5;   //!
   TBranch        *b_L1_SingleLooseIsoEG28_FWD2p5;   //!
   TBranch        *b_L1_SingleLooseIsoEG28er1p5;   //!
   TBranch        *b_L1_SingleLooseIsoEG28er2p1;   //!
   TBranch        *b_L1_SingleLooseIsoEG28er2p5;   //!
   TBranch        *b_L1_SingleLooseIsoEG30er1p5;   //!
   TBranch        *b_L1_SingleLooseIsoEG30er2p5;   //!
   TBranch        *b_L1_SingleMu0_BMTF;   //!
   TBranch        *b_L1_SingleMu0_DQ;   //!
   TBranch        *b_L1_SingleMu0_EMTF;   //!
   TBranch        *b_L1_SingleMu0_OMTF;   //!
   TBranch        *b_L1_SingleMu0_SQ13_BMTF;   //!
   TBranch        *b_L1_SingleMu0_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu0_SQ15_BMTF;   //!
   TBranch        *b_L1_SingleMu0_Upt10;   //!
   TBranch        *b_L1_SingleMu0_Upt10_BMTF;   //!
   TBranch        *b_L1_SingleMu0_Upt10_EMTF;   //!
   TBranch        *b_L1_SingleMu0_Upt10_OMTF;   //!
   TBranch        *b_L1_SingleMu0_Upt10_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu0_Upt15_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu0_Upt20_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu0_Upt25_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu10_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu11_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu12_DQ_BMTF;   //!
   TBranch        *b_L1_SingleMu12_DQ_EMTF;   //!
   TBranch        *b_L1_SingleMu12_DQ_OMTF;   //!
   TBranch        *b_L1_SingleMu15_DQ;   //!
   TBranch        *b_L1_SingleMu18;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu22;   //!
   TBranch        *b_L1_SingleMu22_BMTF;   //!
   TBranch        *b_L1_SingleMu22_BMTF_NEG;   //!
   TBranch        *b_L1_SingleMu22_BMTF_POS;   //!
   TBranch        *b_L1_SingleMu22_DQ;   //!
   TBranch        *b_L1_SingleMu22_EMTF;   //!
   TBranch        *b_L1_SingleMu22_EMTF_NEG;   //!
   TBranch        *b_L1_SingleMu22_EMTF_POS;   //!
   TBranch        *b_L1_SingleMu22_OMTF;   //!
   TBranch        *b_L1_SingleMu22_OMTF_NEG;   //!
   TBranch        *b_L1_SingleMu22_OMTF_POS;   //!
   TBranch        *b_L1_SingleMu22_OQ;   //!
   TBranch        *b_L1_SingleMu25;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu5;   //!
   TBranch        *b_L1_SingleMu5_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu6_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMu7_DQ;   //!
   TBranch        *b_L1_SingleMu7_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu8_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMu9_SQ14_BMTF;   //!
   TBranch        *b_L1_SingleMuCosmics;   //!
   TBranch        *b_L1_SingleMuCosmics_BMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_EMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_OMTF;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_BMTF;   //!
   TBranch        *b_L1_SingleMuOpen_EMTF;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_L1_SingleMuOpen_OMTF;   //!
   TBranch        *b_L1_SingleMuOpen_er1p1_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleMuOpen_er1p4_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleMuShower_Nominal;   //!
   TBranch        *b_L1_SingleMuShower_Tight;   //!
   TBranch        *b_L1_SingleTau120er2p1;   //!
   TBranch        *b_L1_SingleTau130er2p1;   //!
   TBranch        *b_L1_TOTEM_1;   //!
   TBranch        *b_L1_TOTEM_2;   //!
   TBranch        *b_L1_TOTEM_3;   //!
   TBranch        *b_L1_TOTEM_4;   //!
   TBranch        *b_L1_TripleEG_18_17_8_er2p5;   //!
   TBranch        *b_L1_TripleEG_18_18_12_er2p5;   //!
   TBranch        *b_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5;   //!
   TBranch        *b_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5;   //!
   TBranch        *b_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5;   //!
   TBranch        *b_L1_TripleMu0;   //!
   TBranch        *b_L1_TripleMu0_OQ;   //!
   TBranch        *b_L1_TripleMu0_SQ;   //!
   TBranch        *b_L1_TripleMu3;   //!
   TBranch        *b_L1_TripleMu3_SQ;   //!
   TBranch        *b_L1_TripleMu_3SQ_2p5SQ_0;   //!
   TBranch        *b_L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12;   //!
   TBranch        *b_L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12;   //!
   TBranch        *b_L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5_3_3;   //!
   TBranch        *b_L1_TripleMu_5_3_3_SQ;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_5_3;   //!
   TBranch        *b_L1_TwoMuShower_Loose;   //!
   TBranch        *b_L1_UnpairedBunchBptxMinus;   //!
   TBranch        *b_L1_UnpairedBunchBptxPlus;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_copy;   //!
   TBranch        *b_L1_UnprefireableEvent_TriggerRules;   //!
   TBranch        *b_L1_UnprefireableEvent_FirstBxInTrain;   //!
   TBranch        *b_L1_FinalOR_BXmin1;   //!
   TBranch        *b_L1_FinalOR_BXmin2;   //!
   TBranch        *b_L1Reco_step_pRECO;   //!
   TBranch        *b_Flag_HBHENoiseFilter_pRECO;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter_pRECO;   //!
   TBranch        *b_Flag_CSCTightHaloFilter_pRECO;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter_pRECO;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter_pRECO;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter_pRECO;   //!
   TBranch        *b_Flag_HcalStripHaloFilter_pRECO;   //!
   TBranch        *b_Flag_hcalLaserEventFilter_pRECO;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter_pRECO;   //!
   TBranch        *b_Flag_ecalBadCalibFilter_pRECO;   //!
   TBranch        *b_Flag_goodVertices_pRECO;   //!
   TBranch        *b_Flag_eeBadScFilter_pRECO;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter_pRECO;   //!
   TBranch        *b_Flag_trkPOGFilters_pRECO;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter_pRECO;   //!
   TBranch        *b_Flag_muonBadTrackFilter_pRECO;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter_pRECO;   //!
   TBranch        *b_Flag_BadPFMuonFilter_pRECO;   //!
   TBranch        *b_Flag_BadPFMuonDzFilter_pRECO;   //!
   TBranch        *b_Flag_hfNoisyHitsFilter_pRECO;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter_pRECO;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter_pRECO;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X_pRECO;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X_pRECO;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters_pRECO;   //!
   TBranch        *b_Flag_METFilters_pRECO;   //!
   TBranch        *b_L1_AXO_Loose_pRECO;   //!
   TBranch        *b_L1_AXO_Nominal_pRECO;   //!
   TBranch        *b_L1_AXO_Tight_pRECO;   //!
   TBranch        *b_L1_AXO_VLoose_pRECO;   //!
   TBranch        *b_L1_AXO_VTight_pRECO;   //!
   TBranch        *b_L1_AlwaysTrue_pRECO;   //!
   TBranch        *b_L1_BPTX_AND_Ref1_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_AND_Ref3_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_AND_Ref4_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_BeamGas_B1_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_BeamGas_B2_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref1_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref2_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_NotOR_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_OR_Ref3_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_OR_Ref4_VME_pRECO;   //!
   TBranch        *b_L1_BPTX_RefAND_VME_pRECO;   //!
   TBranch        *b_L1_BptxMinus_pRECO;   //!
   TBranch        *b_L1_BptxOR_pRECO;   //!
   TBranch        *b_L1_BptxPlus_pRECO;   //!
   TBranch        *b_L1_BptxXOR_pRECO;   //!
   TBranch        *b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_pRECO;   //!
   TBranch        *b_L1_CICADA_Loose_pRECO;   //!
   TBranch        *b_L1_CICADA_Medium_pRECO;   //!
   TBranch        *b_L1_CICADA_Tight_pRECO;   //!
   TBranch        *b_L1_CICADA_VLoose_pRECO;   //!
   TBranch        *b_L1_CICADA_VTight_pRECO;   //!
   TBranch        *b_L1_DoubleEG11_er1p2_dR_Max0p6_pRECO;   //!
   TBranch        *b_L1_DoubleEG8er2p5_HTT280er_pRECO;   //!
   TBranch        *b_L1_DoubleEG8er2p5_HTT300er_pRECO;   //!
   TBranch        *b_L1_DoubleEG8er2p5_HTT320er_pRECO;   //!
   TBranch        *b_L1_DoubleEG_15_10_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_20_10_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_22_10_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_25_12_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_25_14_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_27_14_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_LooseIso16_LooseIso12_er1p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_LooseIso18_LooseIso12_er1p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_LooseIso20_LooseIso12_er1p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_LooseIso22_12_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_LooseIso22_LooseIso12_er1p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_LooseIso25_12_er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleEG_LooseIso25_LooseIso12_er1p5_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau28er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau30er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau32er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau32er2p1_Mass_Max80_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau34er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau35er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleIsoTau36er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleJet100er2p3_dEta_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleJet100er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet112er2p3_dEta_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleJet120er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet120er2p5_Mu3_dR_Max0p8_pRECO;   //!
   TBranch        *b_L1_DoubleJet150er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet16er2p5_Mu3_dR_Max0p4_pRECO;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet35er2p5_Mu3_dR_Max0p4_pRECO;   //!
   TBranch        *b_L1_DoubleJet40er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5_pRECO;   //!
   TBranch        *b_L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO;   //!
   TBranch        *b_L1_DoubleJet60er2p5_Mu3_dR_Max0p4_pRECO;   //!
   TBranch        *b_L1_DoubleJet80er2p5_Mu3_dR_Max0p4_pRECO;   //!
   TBranch        *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min800_pRECO;   //!
   TBranch        *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min850_pRECO;   //!
   TBranch        *b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50_pRECO;   //!
   TBranch        *b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50_pRECO;   //!
   TBranch        *b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65_pRECO;   //!
   TBranch        *b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65_pRECO;   //!
   TBranch        *b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ_pRECO;   //!
   TBranch        *b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ_pRECO;   //!
   TBranch        *b_L1_DoubleLLPJet40_pRECO;   //!
   TBranch        *b_L1_DoubleLooseIsoEG22er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleLooseIsoEG24er2p1_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_Mass_Min1_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_OQ_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_SQ_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_SQ_OS_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_Upt15_Upt7_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_Upt5_Upt5_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_Upt6_IP_Min1_Upt4_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_Upt6_SQ_er2p0_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_Upt7_SQ_er2p0_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_Upt8_SQ_er2p0_pRECO;   //!
   TBranch        *b_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_dR_Max1p4_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_dEta_Max1p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu0er2p0_SQ_dEta_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleMu18er2p1_SQ_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF30_HTT60er_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF40_HTT60er_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_HTT60er_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT220er_pRECO;   //!
   TBranch        *b_L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleMu4_SQ_EG9er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS_pRECO;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS_dR_Max1p2_pRECO;   //!
   TBranch        *b_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS_pRECO;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_pRECO;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_pRECO;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18_pRECO;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_pRECO;   //!
   TBranch        *b_L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20_pRECO;   //!
   TBranch        *b_L1_DoubleMu5_SQ_EG9er2p5_pRECO;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS_dR_Max1p6_pRECO;   //!
   TBranch        *b_L1_DoubleMu6_Upt6_SQ_er2p0_pRECO;   //!
   TBranch        *b_L1_DoubleMu7_Upt7_SQ_er2p0_pRECO;   //!
   TBranch        *b_L1_DoubleMu8_SQ_pRECO;   //!
   TBranch        *b_L1_DoubleMu8_Upt8_SQ_er2p0_pRECO;   //!
   TBranch        *b_L1_DoubleMu9_SQ_pRECO;   //!
   TBranch        *b_L1_DoubleMu_12_5_pRECO;   //!
   TBranch        *b_L1_DoubleMu_15_5_SQ_pRECO;   //!
   TBranch        *b_L1_DoubleMu_15_7_pRECO;   //!
   TBranch        *b_L1_DoubleMu_15_7_Mass_Min1_pRECO;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ_pRECO;   //!
   TBranch        *b_L1_DoubleTau70er2p1_pRECO;   //!
   TBranch        *b_L1_ETM120_pRECO;   //!
   TBranch        *b_L1_ETM150_pRECO;   //!
   TBranch        *b_L1_ETMHF100_pRECO;   //!
   TBranch        *b_L1_ETMHF100_HTT60er_pRECO;   //!
   TBranch        *b_L1_ETMHF110_pRECO;   //!
   TBranch        *b_L1_ETMHF110_HTT60er_pRECO;   //!
   TBranch        *b_L1_ETMHF120_pRECO;   //!
   TBranch        *b_L1_ETMHF120_HTT60er_pRECO;   //!
   TBranch        *b_L1_ETMHF130_pRECO;   //!
   TBranch        *b_L1_ETMHF130_HTT60er_pRECO;   //!
   TBranch        *b_L1_ETMHF140_pRECO;   //!
   TBranch        *b_L1_ETMHF150_pRECO;   //!
   TBranch        *b_L1_ETMHF70_pRECO;   //!
   TBranch        *b_L1_ETMHF70_HTT60er_pRECO;   //!
   TBranch        *b_L1_ETMHF80_pRECO;   //!
   TBranch        *b_L1_ETMHF80_HTT60er_pRECO;   //!
   TBranch        *b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1_pRECO;   //!
   TBranch        *b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6_pRECO;   //!
   TBranch        *b_L1_ETMHF90_pRECO;   //!
   TBranch        *b_L1_ETMHF90_HTT60er_pRECO;   //!
   TBranch        *b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1_pRECO;   //!
   TBranch        *b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6_pRECO;   //!
   TBranch        *b_L1_ETT2000_pRECO;   //!
   TBranch        *b_L1_FirstBunchAfterTrain_pRECO;   //!
   TBranch        *b_L1_FirstBunchBeforeTrain_pRECO;   //!
   TBranch        *b_L1_FirstBunchInTrain_pRECO;   //!
   TBranch        *b_L1_FirstCollisionInOrbit_pRECO;   //!
   TBranch        *b_L1_FirstCollisionInTrain_pRECO;   //!
   TBranch        *b_L1_HCAL_LaserMon_Trig_pRECO;   //!
   TBranch        *b_L1_HCAL_LaserMon_Veto_pRECO;   //!
   TBranch        *b_L1_HTMHF100_pRECO;   //!
   TBranch        *b_L1_HTMHF120_pRECO;   //!
   TBranch        *b_L1_HTMHF125_pRECO;   //!
   TBranch        *b_L1_HTMHF130_pRECO;   //!
   TBranch        *b_L1_HTMHF150_pRECO;   //!
   TBranch        *b_L1_HTT120_SingleLLPJet40_pRECO;   //!
   TBranch        *b_L1_HTT120er_pRECO;   //!
   TBranch        *b_L1_HTT160_SingleLLPJet50_pRECO;   //!
   TBranch        *b_L1_HTT160er_pRECO;   //!
   TBranch        *b_L1_HTT200_SingleLLPJet60_pRECO;   //!
   TBranch        *b_L1_HTT200er_pRECO;   //!
   TBranch        *b_L1_HTT240_SingleLLPJet70_pRECO;   //!
   TBranch        *b_L1_HTT255er_pRECO;   //!
   TBranch        *b_L1_HTT280er_pRECO;   //!
   TBranch        *b_L1_HTT280er_QuadJet_70_55_40_35_er2p5_pRECO;   //!
   TBranch        *b_L1_HTT320er_pRECO;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p5_pRECO;   //!
   TBranch        *b_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3_pRECO;   //!
   TBranch        *b_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3_pRECO;   //!
   TBranch        *b_L1_HTT360er_pRECO;   //!
   TBranch        *b_L1_HTT400er_pRECO;   //!
   TBranch        *b_L1_HTT450er_pRECO;   //!
   TBranch        *b_L1_IsoEG32er2p5_Mt40_pRECO;   //!
   TBranch        *b_L1_IsoTau52er2p1_QuadJet36er2p5_pRECO;   //!
   TBranch        *b_L1_IsolatedBunch_pRECO;   //!
   TBranch        *b_L1_LastBunchInTrain_pRECO;   //!
   TBranch        *b_L1_LastCollisionInTrain_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG14er2p5_HTT200er_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG16er2p5_HTT200er_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_HTT100er_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_HTT100er_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_HTT100er_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_pRECO;   //!
   TBranch        *b_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_pRECO;   //!
   TBranch        *b_L1_MinimumBiasHF0_pRECO;   //!
   TBranch        *b_L1_MinimumBiasHF0_AND_BptxAND_pRECO;   //!
   TBranch        *b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6_pRECO;   //!
   TBranch        *b_L1_Mu12_HTT150er_pRECO;   //!
   TBranch        *b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6_pRECO;   //!
   TBranch        *b_L1_Mu14_HTT150er_pRECO;   //!
   TBranch        *b_L1_Mu18er2p1_Tau24er2p1_pRECO;   //!
   TBranch        *b_L1_Mu18er2p1_Tau26er2p1_pRECO;   //!
   TBranch        *b_L1_Mu18er2p1_Tau26er2p1_Jet55_pRECO;   //!
   TBranch        *b_L1_Mu18er2p1_Tau26er2p1_Jet70_pRECO;   //!
   TBranch        *b_L1_Mu20_EG10er2p5_pRECO;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau30er2p1_pRECO;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau32er2p1_pRECO;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau34er2p1_pRECO;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau40er2p1_pRECO;   //!
   TBranch        *b_L1_Mu22er2p1_Tau70er2p1_pRECO;   //!
   TBranch        *b_L1_Mu3_Jet120er2p5_dR_Max0p4_pRECO;   //!
   TBranch        *b_L1_Mu3_Jet16er2p5_dR_Max0p4_pRECO;   //!
   TBranch        *b_L1_Mu3_Jet30er2p5_pRECO;   //!
   TBranch        *b_L1_Mu3_Jet60er2p5_dR_Max0p4_pRECO;   //!
   TBranch        *b_L1_Mu3er1p5_Jet100er2p5_ETMHF30_pRECO;   //!
   TBranch        *b_L1_Mu3er1p5_Jet100er2p5_ETMHF40_pRECO;   //!
   TBranch        *b_L1_Mu3er1p5_Jet100er2p5_ETMHF50_pRECO;   //!
   TBranch        *b_L1_Mu5_EG23er2p5_pRECO;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG20er2p5_pRECO;   //!
   TBranch        *b_L1_Mu6_DoubleEG10er2p5_pRECO;   //!
   TBranch        *b_L1_Mu6_DoubleEG12er2p5_pRECO;   //!
   TBranch        *b_L1_Mu6_DoubleEG15er2p5_pRECO;   //!
   TBranch        *b_L1_Mu6_DoubleEG17er2p5_pRECO;   //!
   TBranch        *b_L1_Mu6_HTT240er_pRECO;   //!
   TBranch        *b_L1_Mu6_HTT250er_pRECO;   //!
   TBranch        *b_L1_Mu7_EG20er2p5_pRECO;   //!
   TBranch        *b_L1_Mu7_EG23er2p5_pRECO;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG20er2p5_pRECO;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG23er2p5_pRECO;   //!
   TBranch        *b_L1_NotBptxOR_pRECO;   //!
   TBranch        *b_L1_QuadJet60er2p5_pRECO;   //!
   TBranch        *b_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0_pRECO;   //!
   TBranch        *b_L1_QuadMu0_pRECO;   //!
   TBranch        *b_L1_QuadMu0_OQ_pRECO;   //!
   TBranch        *b_L1_QuadMu0_SQ_pRECO;   //!
   TBranch        *b_L1_SecondBunchInTrain_pRECO;   //!
   TBranch        *b_L1_SecondLastBunchInTrain_pRECO;   //!
   TBranch        *b_L1_SingleEG10er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG15er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG26er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG28_FWD2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG28er1p5_pRECO;   //!
   TBranch        *b_L1_SingleEG28er2p1_pRECO;   //!
   TBranch        *b_L1_SingleEG28er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG34er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG36er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG38er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG40er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG42er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG45er2p5_pRECO;   //!
   TBranch        *b_L1_SingleEG50_pRECO;   //!
   TBranch        *b_L1_SingleEG60_pRECO;   //!
   TBranch        *b_L1_SingleEG8er2p5_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG24er2p1_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG26er2p1_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG26er2p5_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG28_FWD2p5_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG28er1p5_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG28er2p1_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG28er2p5_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG30er2p1_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG30er2p5_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG32er2p1_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG32er2p5_pRECO;   //!
   TBranch        *b_L1_SingleIsoEG34er2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet10erHE_pRECO;   //!
   TBranch        *b_L1_SingleJet120_pRECO;   //!
   TBranch        *b_L1_SingleJet120_FWD2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet120_FWD3p0_pRECO;   //!
   TBranch        *b_L1_SingleJet120er1p3_pRECO;   //!
   TBranch        *b_L1_SingleJet120er2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet12erHE_pRECO;   //!
   TBranch        *b_L1_SingleJet140er2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet160er2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet180_pRECO;   //!
   TBranch        *b_L1_SingleJet180er2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet200_pRECO;   //!
   TBranch        *b_L1_SingleJet20er2p5_NotBptxOR_pRECO;   //!
   TBranch        *b_L1_SingleJet20er2p5_NotBptxOR_3BX_pRECO;   //!
   TBranch        *b_L1_SingleJet35_pRECO;   //!
   TBranch        *b_L1_SingleJet35_FWD2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet35_FWD3p0_pRECO;   //!
   TBranch        *b_L1_SingleJet35er1p3_pRECO;   //!
   TBranch        *b_L1_SingleJet35er2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet43er2p5_NotBptxOR_3BX_pRECO;   //!
   TBranch        *b_L1_SingleJet46er2p5_NotBptxOR_3BX_pRECO;   //!
   TBranch        *b_L1_SingleJet60_pRECO;   //!
   TBranch        *b_L1_SingleJet60_FWD2p5_pRECO;   //!
   TBranch        *b_L1_SingleJet8erHE_pRECO;   //!
   TBranch        *b_L1_SingleJet90_pRECO;   //!
   TBranch        *b_L1_SingleJet90_FWD2p5_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG26er1p5_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG26er2p5_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG28_FWD2p5_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG28er1p5_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG28er2p1_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG28er2p5_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG30er1p5_pRECO;   //!
   TBranch        *b_L1_SingleLooseIsoEG30er2p5_pRECO;   //!
   TBranch        *b_L1_SingleMu0_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_DQ_pRECO;   //!
   TBranch        *b_L1_SingleMu0_EMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_OMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_SQ13_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_SQ15_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt10_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt10_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt10_EMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt10_OMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt10_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt15_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt20_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu0_Upt25_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu10_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu11_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu12_DQ_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu12_DQ_EMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu12_DQ_OMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu15_DQ_pRECO;   //!
   TBranch        *b_L1_SingleMu18_pRECO;   //!
   TBranch        *b_L1_SingleMu20_pRECO;   //!
   TBranch        *b_L1_SingleMu22_pRECO;   //!
   TBranch        *b_L1_SingleMu22_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu22_BMTF_NEG_pRECO;   //!
   TBranch        *b_L1_SingleMu22_BMTF_POS_pRECO;   //!
   TBranch        *b_L1_SingleMu22_DQ_pRECO;   //!
   TBranch        *b_L1_SingleMu22_EMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu22_EMTF_NEG_pRECO;   //!
   TBranch        *b_L1_SingleMu22_EMTF_POS_pRECO;   //!
   TBranch        *b_L1_SingleMu22_OMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu22_OMTF_NEG_pRECO;   //!
   TBranch        *b_L1_SingleMu22_OMTF_POS_pRECO;   //!
   TBranch        *b_L1_SingleMu22_OQ_pRECO;   //!
   TBranch        *b_L1_SingleMu25_pRECO;   //!
   TBranch        *b_L1_SingleMu3_pRECO;   //!
   TBranch        *b_L1_SingleMu5_pRECO;   //!
   TBranch        *b_L1_SingleMu5_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu6_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu7_pRECO;   //!
   TBranch        *b_L1_SingleMu7_DQ_pRECO;   //!
   TBranch        *b_L1_SingleMu7_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu8_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMu9_SQ14_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMuCosmics_pRECO;   //!
   TBranch        *b_L1_SingleMuCosmics_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMuCosmics_EMTF_pRECO;   //!
   TBranch        *b_L1_SingleMuCosmics_OMTF_pRECO;   //!
   TBranch        *b_L1_SingleMuOpen_pRECO;   //!
   TBranch        *b_L1_SingleMuOpen_BMTF_pRECO;   //!
   TBranch        *b_L1_SingleMuOpen_EMTF_pRECO;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR_pRECO;   //!
   TBranch        *b_L1_SingleMuOpen_OMTF_pRECO;   //!
   TBranch        *b_L1_SingleMuOpen_er1p1_NotBptxOR_3BX_pRECO;   //!
   TBranch        *b_L1_SingleMuOpen_er1p4_NotBptxOR_3BX_pRECO;   //!
   TBranch        *b_L1_SingleMuShower_Nominal_pRECO;   //!
   TBranch        *b_L1_SingleMuShower_Tight_pRECO;   //!
   TBranch        *b_L1_SingleTau120er2p1_pRECO;   //!
   TBranch        *b_L1_SingleTau130er2p1_pRECO;   //!
   TBranch        *b_L1_TOTEM_1_pRECO;   //!
   TBranch        *b_L1_TOTEM_2_pRECO;   //!
   TBranch        *b_L1_TOTEM_3_pRECO;   //!
   TBranch        *b_L1_TOTEM_4_pRECO;   //!
   TBranch        *b_L1_TripleEG_18_17_8_er2p5_pRECO;   //!
   TBranch        *b_L1_TripleEG_18_18_12_er2p5_pRECO;   //!
   TBranch        *b_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5_pRECO;   //!
   TBranch        *b_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5_pRECO;   //!
   TBranch        *b_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5_pRECO;   //!
   TBranch        *b_L1_TripleMu0_pRECO;   //!
   TBranch        *b_L1_TripleMu0_OQ_pRECO;   //!
   TBranch        *b_L1_TripleMu0_SQ_pRECO;   //!
   TBranch        *b_L1_TripleMu3_pRECO;   //!
   TBranch        *b_L1_TripleMu3_SQ_pRECO;   //!
   TBranch        *b_L1_TripleMu_3SQ_2p5SQ_0_pRECO;   //!
   TBranch        *b_L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12_pRECO;   //!
   TBranch        *b_L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12_pRECO;   //!
   TBranch        *b_L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12_pRECO;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ_pRECO;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO;   //!
   TBranch        *b_L1_TripleMu_5_3_3_pRECO;   //!
   TBranch        *b_L1_TripleMu_5_3_3_SQ_pRECO;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5_pRECO;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO;   //!
   TBranch        *b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO;   //!
   TBranch        *b_L1_TripleMu_5_5_3_pRECO;   //!
   TBranch        *b_L1_TwoMuShower_Loose_pRECO;   //!
   TBranch        *b_L1_UnpairedBunchBptxMinus_pRECO;   //!
   TBranch        *b_L1_UnpairedBunchBptxPlus_pRECO;   //!
   TBranch        *b_L1_ZeroBias_pRECO;   //!
   TBranch        *b_L1_ZeroBias_copy_pRECO;   //!
   TBranch        *b_L1_UnprefireableEvent_TriggerRules_pRECO;   //!
   TBranch        *b_L1_UnprefireableEvent_FirstBxInTrain_pRECO;   //!
   TBranch        *b_L1_FinalOR_BXmin1_pRECO;   //!
   TBranch        *b_L1_FinalOR_BXmin2_pRECO;   //!
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_DST_PFScouting_DoubleMuon;   //!
   TBranch        *b_DST_PFScouting_DoubleEG;   //!
   TBranch        *b_DST_PFScouting_JetHT;   //!
   TBranch        *b_DST_PFScouting_DatasetMuon;   //!
   TBranch        *b_DST_PFScouting_AXOVLoose;   //!
   TBranch        *b_DST_PFScouting_AXOLoose;   //!
   TBranch        *b_DST_PFScouting_AXONominal;   //!
   TBranch        *b_DST_PFScouting_AXOTight;   //!
   TBranch        *b_DST_PFScouting_AXOVTight;   //!
   TBranch        *b_DST_PFScouting_CICADAVLoose;   //!
   TBranch        *b_DST_PFScouting_CICADALoose;   //!
   TBranch        *b_DST_PFScouting_CICADAMedium;   //!
   TBranch        *b_DST_PFScouting_CICADATight;   //!
   TBranch        *b_DST_PFScouting_CICADAVTight;   //!
   TBranch        *b_DST_PFScouting_SingleMuon;   //!
   TBranch        *b_DST_PFScouting_SinglePhotonEB;   //!
   TBranch        *b_DST_PFScouting_ZeroBias;   //!
   TBranch        *b_HLT_EphemeralPhysics;   //!
   TBranch        *b_HLT_TestPhysics;   //!
   TBranch        *b_HLT_EphemeralZeroBias;   //!
   TBranch        *b_HLT_SpecialZeroBias;   //!
   TBranch        *b_HLT_EcalCalibration;   //!
   TBranch        *b_HLT_HcalCalibration;   //!
   TBranch        *b_HLT_HcalNZS;   //!
   TBranch        *b_HLT_HcalPhiSym;   //!
   TBranch        *b_HLT_Random;   //!
   TBranch        *b_HLT_Physics;   //!
   TBranch        *b_HLT_ZeroBias;   //!
   TBranch        *b_HLT_ZeroBias_Alignment;   //!
   TBranch        *b_HLT_ZeroBias_Beamspot;   //!
   TBranch        *b_HLT_ZeroBias_IsolatedBunches;   //!
   TBranch        *b_HLT_ZeroBias_FirstBXAfterTrain;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionAfterAbortGap;   //!
   TBranch        *b_HLT_ZeroBias_FirstCollisionInTrain;   //!
   TBranch        *b_HLT_ZeroBias_LastCollisionInTrain;   //!
   TBranch        *b_HLT_HT300_Beamspot;   //!
   TBranch        *b_HLT_IsoTrackHB;   //!
   TBranch        *b_HLT_IsoTrackHE;   //!
   TBranch        *b_HLT_PFJet40_GPUvsCPU;   //!
   TBranch        *b_HLT_AK8PFJet380_SoftDropMass30;   //!
   TBranch        *b_HLT_AK8PFJet400_SoftDropMass30;   //!
   TBranch        *b_HLT_AK8PFJet425_SoftDropMass30;   //!
   TBranch        *b_HLT_AK8PFJet450_SoftDropMass30;   //!
   TBranch        *b_HLT_AK8DiPFJet250_250_SoftDropMass40;   //!
   TBranch        *b_HLT_AK8DiPFJet250_250_SoftDropMass50;   //!
   TBranch        *b_HLT_AK8DiPFJet260_260_SoftDropMass30;   //!
   TBranch        *b_HLT_AK8DiPFJet260_260_SoftDropMass40;   //!
   TBranch        *b_HLT_AK8DiPFJet270_270_SoftDropMass30;   //!
   TBranch        *b_HLT_AK8DiPFJet280_280_SoftDropMass30;   //!
   TBranch        *b_HLT_AK8DiPFJet290_290_SoftDropMass30;   //!
   TBranch        *b_HLT_CaloJet500_NoJetID;   //!
   TBranch        *b_HLT_CaloJet550_NoJetID;   //!
   TBranch        *b_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon;   //!
   TBranch        *b_HLT_Trimuon5_3p5_2_Upsilon_Muon;   //!
   TBranch        *b_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon;   //!
   TBranch        *b_HLT_DoubleEle25_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle27_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle33_CaloIdL_MW;   //!
   TBranch        *b_HLT_DoubleEle24_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350;   //!
   TBranch        *b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Mu27_Ele37_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_Ele27_CaloIdL_MW;   //!
   TBranch        *b_HLT_Mu37_TkMu27;   //!
   TBranch        *b_HLT_DoubleMu4_3_Bs;   //!
   TBranch        *b_HLT_DoubleMu4_3_Jpsi;   //!
   TBranch        *b_HLT_DoubleMu4_3_LowMass;   //!
   TBranch        *b_HLT_DoubleMu4_LowMass_Displaced;   //!
   TBranch        *b_HLT_Mu0_L1DoubleMu;   //!
   TBranch        *b_HLT_Mu4_L1DoubleMu;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_LowPt;   //!
   TBranch        *b_HLT_DoubleMu4_3_Photon4_BsToMMG;   //!
   TBranch        *b_HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG;   //!
   TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu;   //!
   TBranch        *b_HLT_DoubleMu3_TkMu_DsTau3Mu;   //!
   TBranch        *b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350;   //!
   TBranch        *b_HLT_DoubleMu4_MuMuTrk_Displaced;   //!
   TBranch        *b_HLT_Mu3_PFJet40;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Jpsi;   //!
   TBranch        *b_HLT_Mu7p5_L2Mu2_Upsilon;   //!
   TBranch        *b_HLT_Mu3_L1SingleMu5orSingleMu7;   //!
   TBranch        *b_HLT_Mu0_Barrel;   //!
   TBranch        *b_HLT_Mu0_Barrel_L1HP6;   //!
   TBranch        *b_HLT_Mu0_Barrel_L1HP7;   //!
   TBranch        *b_HLT_Mu0_Barrel_L1HP8;   //!
   TBranch        *b_HLT_Mu0_Barrel_L1HP9;   //!
   TBranch        *b_HLT_Mu0_Barrel_L1HP10;   //!
   TBranch        *b_HLT_Mu0_Barrel_L1HP11;   //!
   TBranch        *b_HLT_Mu0_Barrel_L1HP6_IP6;   //!
   TBranch        *b_HLT_Mu6_Barrel_L1HP7_IP6;   //!
   TBranch        *b_HLT_Mu7_Barrel_L1HP8_IP6;   //!
   TBranch        *b_HLT_Mu8_Barrel_L1HP9_IP6;   //!
   TBranch        *b_HLT_Mu9_Barrel_L1HP10_IP6;   //!
   TBranch        *b_HLT_Mu10_Barrel_L1HP11_IP6;   //!
   TBranch        *b_HLT_DoublePhoton33_CaloIdL;   //!
   TBranch        *b_HLT_DoublePhoton70;   //!
   TBranch        *b_HLT_DoublePhoton85;   //!
   TBranch        *b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG;   //!
   TBranch        *b_HLT_Ele30_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele38_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele40_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_L1DoubleEG;   //!
   TBranch        *b_HLT_IsoMu20;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX;   //!
   TBranch        *b_HLT_UncorrectedJetE30_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE60_NoBPTX3BX;   //!
   TBranch        *b_HLT_UncorrectedJetE70_NoBPTX3BX;   //!
   TBranch        *b_HLT_L1SingleMuCosmics;   //!
   TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu10_NoVertex_NoBPTX;   //!
   TBranch        *b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;   //!
   TBranch        *b_HLT_L2Mu23NoVtx_2Cha;   //!
   TBranch        *b_HLT_L2Mu23NoVtx_2Cha_CosmicSeed;   //!
   TBranch        *b_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4;   //!
   TBranch        *b_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4;   //!
   TBranch        *b_HLT_DoubleL2Mu50;   //!
   TBranch        *b_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;   //!
   TBranch        *b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed;   //!
   TBranch        *b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4;   //!
   TBranch        *b_HLT_DoubleL2Mu23NoVtx_2Cha;   //!
   TBranch        *b_HLT_DoubleL2Mu25NoVtx_2Cha;   //!
   TBranch        *b_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_Mu30_TkMu0_Psi;   //!
   TBranch        *b_HLT_Mu30_TkMu0_Upsilon;   //!
   TBranch        *b_HLT_Mu25_TkMu0_Phi;   //!
   TBranch        *b_HLT_Mu15;   //!
   TBranch        *b_HLT_Mu20;   //!
   TBranch        *b_HLT_Mu27;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_CascadeMu100;   //!
   TBranch        *b_HLT_HighPtTkMu100;   //!
   TBranch        *b_HLT_DiPFJetAve40;   //!
   TBranch        *b_HLT_DiPFJetAve60;   //!
   TBranch        *b_HLT_DiPFJetAve80;   //!
   TBranch        *b_HLT_DiPFJetAve140;   //!
   TBranch        *b_HLT_DiPFJetAve200;   //!
   TBranch        *b_HLT_DiPFJetAve260;   //!
   TBranch        *b_HLT_DiPFJetAve320;   //!
   TBranch        *b_HLT_DiPFJetAve400;   //!
   TBranch        *b_HLT_DiPFJetAve500;   //!
   TBranch        *b_HLT_DiPFJetAve60_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve80_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve100_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve160_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve220_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve260_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve300_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP;   //!
   TBranch        *b_HLT_AK8PFJet40;   //!
   TBranch        *b_HLT_AK8PFJet60;   //!
   TBranch        *b_HLT_AK8PFJet80;   //!
   TBranch        *b_HLT_AK8PFJet140;   //!
   TBranch        *b_HLT_AK8PFJet200;   //!
   TBranch        *b_HLT_AK8PFJet260;   //!
   TBranch        *b_HLT_AK8PFJet320;   //!
   TBranch        *b_HLT_AK8PFJet400;   //!
   TBranch        *b_HLT_AK8PFJet450;   //!
   TBranch        *b_HLT_AK8PFJet500;   //!
   TBranch        *b_HLT_AK8PFJet550;   //!
   TBranch        *b_HLT_PFJet40;   //!
   TBranch        *b_HLT_PFJet40_ZeroBiasCopy;   //!
   TBranch        *b_HLT_PFJet60;   //!
   TBranch        *b_HLT_PFJet80;   //!
   TBranch        *b_HLT_PFJet80_L1SingleJet60;   //!
   TBranch        *b_HLT_PFJet110;   //!
   TBranch        *b_HLT_PFJet140;   //!
   TBranch        *b_HLT_PFJet200;   //!
   TBranch        *b_HLT_PFJet260;   //!
   TBranch        *b_HLT_PFJet320;   //!
   TBranch        *b_HLT_PFJet400;   //!
   TBranch        *b_HLT_PFJet450;   //!
   TBranch        *b_HLT_PFJet500;   //!
   TBranch        *b_HLT_PFJet550;   //!
   TBranch        *b_HLT_PFJetFwd40;   //!
   TBranch        *b_HLT_PFJetFwd60;   //!
   TBranch        *b_HLT_PFJetFwd80;   //!
   TBranch        *b_HLT_PFJetFwd140;   //!
   TBranch        *b_HLT_PFJetFwd200;   //!
   TBranch        *b_HLT_PFJetFwd260;   //!
   TBranch        *b_HLT_PFJetFwd320;   //!
   TBranch        *b_HLT_PFJetFwd400;   //!
   TBranch        *b_HLT_PFJetFwd450;   //!
   TBranch        *b_HLT_PFJetFwd500;   //!
   TBranch        *b_HLT_AK8PFJetFwd40;   //!
   TBranch        *b_HLT_AK8PFJetFwd60;   //!
   TBranch        *b_HLT_AK8PFJetFwd80;   //!
   TBranch        *b_HLT_AK8PFJetFwd140;   //!
   TBranch        *b_HLT_AK8PFJetFwd200;   //!
   TBranch        *b_HLT_AK8PFJetFwd260;   //!
   TBranch        *b_HLT_AK8PFJetFwd320;   //!
   TBranch        *b_HLT_AK8PFJetFwd400;   //!
   TBranch        *b_HLT_AK8PFJetFwd450;   //!
   TBranch        *b_HLT_AK8PFJetFwd500;   //!
   TBranch        *b_HLT_PFHT180;   //!
   TBranch        *b_HLT_PFHT250;   //!
   TBranch        *b_HLT_PFHT370;   //!
   TBranch        *b_HLT_PFHT430;   //!
   TBranch        *b_HLT_PFHT510;   //!
   TBranch        *b_HLT_PFHT590;   //!
   TBranch        *b_HLT_PFHT680;   //!
   TBranch        *b_HLT_PFHT780;   //!
   TBranch        *b_HLT_PFHT890;   //!
   TBranch        *b_HLT_PFHT1050;   //!
   TBranch        *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_PFHT500_PFMET110_PFMHT110_IDTight;   //!
   TBranch        *b_HLT_PFHT700_PFMET85_PFMHT85_IDTight;   //!
   TBranch        *b_HLT_PFHT800_PFMET75_PFMHT75_IDTight;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight;   //!
   TBranch        *b_HLT_PFMET130_PFMHT130_IDTight;   //!
   TBranch        *b_HLT_PFMET140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;   //!
   TBranch        *b_HLT_PFMETTypeOne140_PFMHT140_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight;   //!
   TBranch        *b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF;   //!
   TBranch        *b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF;   //!
   TBranch        *b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF;   //!
   TBranch        *b_HLT_L1ETMHadSeeds;   //!
   TBranch        *b_HLT_CaloMHT90;   //!
   TBranch        *b_HLT_CaloMET90_NotCleaned;   //!
   TBranch        *b_HLT_CaloMET350_NotCleaned;   //!
   TBranch        *b_HLT_PFMET200_NotCleaned;   //!
   TBranch        *b_HLT_PFMET250_NotCleaned;   //!
   TBranch        *b_HLT_PFMET300_NotCleaned;   //!
   TBranch        *b_HLT_PFMET200_BeamHaloCleaned;   //!
   TBranch        *b_HLT_PFMETTypeOne200_BeamHaloCleaned;   //!
   TBranch        *b_HLT_MET105_IsoTrk50;   //!
   TBranch        *b_HLT_MET120_IsoTrk50;   //!
   TBranch        *b_HLT_Mu12eta2p3;   //!
   TBranch        *b_HLT_Mu12eta2p3_PFJet40;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets40_PNetBTag_0p11;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets100_PNetBTag_0p11;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets200_PNetBTag_0p11;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets350_PNetBTag_0p11;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11;   //!
   TBranch        *b_HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11;   //!
   TBranch        *b_HLT_DoublePFJets40_PNetBTag_0p11;   //!
   TBranch        *b_HLT_DoublePFJets100_PNetBTag_0p11;   //!
   TBranch        *b_HLT_DoublePFJets200_PNetBTag_0p11;   //!
   TBranch        *b_HLT_DoublePFJets350_PNetBTag_0p11;   //!
   TBranch        *b_HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11;   //!
   TBranch        *b_HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11;   //!
   TBranch        *b_HLT_Photon300_NoHE;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ;   //!
   TBranch        *b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu19_TrkIsoVVL;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet20_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet40_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet70_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet110_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4DiJet170_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK4Jet300_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK8DiJet170_Mu5;   //!
   TBranch        *b_HLT_BTagMu_AK8Jet170_DoubleMu5;   //!
   TBranch        *b_HLT_BTagMu_AK8Jet300_Mu5;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53;   //!
   TBranch        *b_HLT_Photon33;   //!
   TBranch        *b_HLT_Photon50;   //!
   TBranch        *b_HLT_Photon75;   //!
   TBranch        *b_HLT_Photon90;   //!
   TBranch        *b_HLT_Photon120;   //!
   TBranch        *b_HLT_Photon150;   //!
   TBranch        *b_HLT_Photon175;   //!
   TBranch        *b_HLT_Photon200;   //!
   TBranch        *b_HLT_Photon45EB;   //!
   TBranch        *b_HLT_Photon40EB;   //!
   TBranch        *b_HLT_Photon50EB;   //!
   TBranch        *b_HLT_Photon30EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon40EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon45EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon50EB_TightID_TightIso_PFJet30;   //!
   TBranch        *b_HLT_Photon50EB_TightID_TightIso_CaloJet30;   //!
   TBranch        *b_HLT_Photon50EB_TightID_TightIso_AK8PFJet30;   //!
   TBranch        *b_HLT_Photon50EB_TightID_TightIso_AK8CaloJet30;   //!
   TBranch        *b_HLT_Photon50EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon55EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon75EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon90EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon110EB_TightID_TightIso;   //!
   TBranch        *b_HLT_Photon110EB_TightID_TightIso_PFJet30;   //!
   TBranch        *b_HLT_Photon110EB_TightID_TightIso_CaloJet30;   //!
   TBranch        *b_HLT_Photon110EB_TightID_TightIso_AK8PFJet30;   //!
   TBranch        *b_HLT_Photon110EB_TightID_TightIso_AK8CaloJet30;   //!
   TBranch        *b_HLT_Photon100EBHE10;   //!
   TBranch        *b_HLT_Photon50_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon90_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon120_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Photon165_R9Id90_HE10_IsoM;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90;   //!
   TBranch        *b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95;   //!
   TBranch        *b_HLT_DiphotonMVA14p25_Mass90;   //!
   TBranch        *b_HLT_DiphotonMVA14p25_Tight_Mass90;   //!
   TBranch        *b_HLT_Photon35_TwoProngs35;   //!
   TBranch        *b_HLT_IsoMu24_TwoProngs35;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_L1_NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R;   //!
   TBranch        *b_HLT_Dimuon0_Jpsi3p5_Muon2;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_NoVertexing;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_0er1p5;   //!
   TBranch        *b_HLT_Dimuon0_LowMass;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_4;   //!
   TBranch        *b_HLT_Dimuon0_LowMass_L1_TM530;   //!
   TBranch        *b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass;   //!
   TBranch        *b_HLT_TripleMu_5_3_3_Mass3p8_DZ;   //!
   TBranch        *b_HLT_TripleMu_10_5_5_DZ;   //!
   TBranch        *b_HLT_TripleMu_12_10_5;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15;   //!
   TBranch        *b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70;   //!
   TBranch        *b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90;   //!
   TBranch        *b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass;   //!
   TBranch        *b_HLT_DoubleMu4_3_LowMass_SS;   //!
   TBranch        *b_HLT_DoubleMu4_Jpsi_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_Jpsi_NoVertexing;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrkTrk_Displaced;   //!
   TBranch        *b_HLT_DoubleMu4_JpsiTrk_Bc;   //!
   TBranch        *b_HLT_DoubleMu43NoFiltersNoVtx;   //!
   TBranch        *b_HLT_DoubleMu48NoFiltersNoVtx;   //!
   TBranch        *b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL;   //!
   TBranch        *b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL;   //!
   TBranch        *b_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL;   //!
   TBranch        *b_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL;   //!
   TBranch        *b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned;   //!
   TBranch        *b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150;   //!
   TBranch        *b_HLT_Ele28_HighEta_SC20_Mass55;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Ele50_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Mu50_IsoVVVL_PFHT450;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight;   //!
   TBranch        *b_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight;   //!
   TBranch        *b_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight;   //!
   TBranch        *b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight;   //!
   TBranch        *b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight;   //!
   TBranch        *b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight;   //!
   TBranch        *b_HLT_Dimuon10_Upsilon_y1p4;   //!
   TBranch        *b_HLT_Dimuon12_Upsilon_y1p4;   //!
   TBranch        *b_HLT_Dimuon14_Phi_Barrel_Seagulls;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi;   //!
   TBranch        *b_HLT_Dimuon14_PsiPrime;   //!
   TBranch        *b_HLT_Dimuon14_PsiPrime_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon18_PsiPrime;   //!
   TBranch        *b_HLT_Dimuon18_PsiPrime_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon24_Upsilon_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon24_Phi_noCorrL1;   //!
   TBranch        *b_HLT_Dimuon25_Jpsi_noCorrL1;   //!
   TBranch        *b_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;   //!
   TBranch        *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;   //!
   TBranch        *b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_DoubleIsoMu20_eta2p1;   //!
   TBranch        *b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx;   //!
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu17;   //!
   TBranch        *b_HLT_Mu19;   //!
   TBranch        *b_HLT_Mu17_Photon30_IsoCaloId;   //!
   TBranch        *b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;   //!
   TBranch        *b_HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele135_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40;   //!
   TBranch        *b_HLT_PFHT400_SixPFJet32;   //!
   TBranch        *b_HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50;   //!
   TBranch        *b_HLT_PFHT450_SixPFJet36;   //!
   TBranch        *b_HLT_PFHT450_SixPFJet36_PNetBTag0p35;   //!
   TBranch        *b_HLT_PFHT400_FivePFJet_120_120_60_30_30;   //!
   TBranch        *b_HLT_PFHT350;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380;   //!
   TBranch        *b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400;   //!
   TBranch        *b_HLT_ECALHT800;   //!
   TBranch        *b_HLT_DiSC30_18_EIso_AND_HE_Mass70;   //!
   TBranch        *b_HLT_Photon20_HoverELoose;   //!
   TBranch        *b_HLT_Photon30_HoverELoose;   //!
   TBranch        *b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_HLT_CDC_L2cosmic_10_er1p0;   //!
   TBranch        *b_HLT_CDC_L2cosmic_5p5_er1p0;   //!
   TBranch        *b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   //!
   TBranch        *b_HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1;   //!
   TBranch        *b_HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;   //!
   TBranch        *b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3;   //!
   TBranch        *b_HLT_Mu18_Mu9_SameSign;   //!
   TBranch        *b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05;   //!
   TBranch        *b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60;   //!
   TBranch        *b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0;   //!
   TBranch        *b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA;   //!
   TBranch        *b_HLT_TripleMu_5_3_3_Mass3p8_DCA;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15;   //!
   TBranch        *b_HLT_QuadPFJet105_88_76_15;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15;   //!
   TBranch        *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId;   //!
   TBranch        *b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55;   //!
   TBranch        *b_HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1;   //!
   TBranch        *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5;   //!
   TBranch        *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3;   //!
   TBranch        *b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0;   //!
   TBranch        *b_HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6;   //!
   TBranch        *b_HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5;   //!
   TBranch        *b_HLT_PFHT250_QuadPFJet25;   //!
   TBranch        *b_HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55;   //!
   TBranch        *b_HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50;   //!
   TBranch        *b_HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55;   //!
   TBranch        *b_HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50;   //!
   TBranch        *b_HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50;   //!
   TBranch        *b_HLT_PFHT280_QuadPFJet30;   //!
   TBranch        *b_HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55;   //!
   TBranch        *b_HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60;   //!
   TBranch        *b_HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60;   //!
   TBranch        *b_HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50;   //!
   TBranch        *b_HLT_QuadPFJet100_88_70_30;   //!
   TBranch        *b_HLT_QuadPFJet105_88_75_30;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_30;   //!
   TBranch        *b_HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight;   //!
   TBranch        *b_HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight;   //!
   TBranch        *b_HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55;   //!
   TBranch        *b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03;   //!
   TBranch        *b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03;   //!
   TBranch        *b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03;   //!
   TBranch        *b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05;   //!
   TBranch        *b_HLT_AK8PFJet275_Nch45;   //!
   TBranch        *b_HLT_AK8PFJet275_Nch40;   //!
   TBranch        *b_HLT_IsoMu50_AK8PFJet220_SoftDropMass40;   //!
   TBranch        *b_HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40;   //!
   TBranch        *b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06;   //!
   TBranch        *b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20;   //!
   TBranch        *b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60;   //!
   TBranch        *b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PFHT250;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50;   //!
   TBranch        *b_HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3;   //!
   TBranch        *b_HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3;   //!
   TBranch        *b_HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60;   //!
   TBranch        *b_HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1;   //!
   TBranch        *b_HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3;   //!
   TBranch        *b_HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3;   //!
   TBranch        *b_HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3;   //!
   TBranch        *b_HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_L2Mu10NoVtx_2Cha;   //!
   TBranch        *b_HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_L3Mu10NoVtx;   //!
   TBranch        *b_HLT_L3Mu10NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm;   //!
   TBranch        *b_HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm;   //!
   TBranch        *b_HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed;   //!
   TBranch        *b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm;   //!
   TBranch        *b_HLT_L3Mu30NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_L3Mu50NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_L3dTksMu10_NoVtx_DxyMin0p01cm;   //!
   TBranch        *b_HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId;   //!
   TBranch        *b_HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24;   //!
   TBranch        *b_HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1;   //!
   TBranch        *b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy;   //!
   TBranch        *b_HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1;   //!
   TBranch        *b_HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive;   //!
   TBranch        *b_HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive;   //!
   TBranch        *b_HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive;   //!
   TBranch        *b_HLT_HT350_DelayedJet40_SingleDelay3nsInclusive;   //!
   TBranch        *b_HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive;   //!
   TBranch        *b_HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay2nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1nsTrackless;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless;   //!
   TBranch        *b_HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless;   //!
   TBranch        *b_HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless;   //!
   TBranch        *b_HLT_L1Mu6HT240;   //!
   TBranch        *b_HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose;   //!
   TBranch        *b_HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5;   //!
   TBranch        *b_HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose;   //!
   TBranch        *b_HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5;   //!
   TBranch        *b_HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose;   //!
   TBranch        *b_HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5;   //!
   TBranch        *b_HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5;   //!
   TBranch        *b_HLT_HT350;   //!
   TBranch        *b_HLT_HT425;   //!
   TBranch        *b_HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT400_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT430_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT550_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_HT650_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_CaloMET60_DTCluster50;   //!
   TBranch        *b_HLT_CaloMET60_DTClusterNoMB1S50;   //!
   TBranch        *b_HLT_L1MET_DTCluster50;   //!
   TBranch        *b_HLT_L1MET_DTClusterNoMB1S50;   //!
   TBranch        *b_HLT_CscCluster_Loose;   //!
   TBranch        *b_HLT_CscCluster_Medium;   //!
   TBranch        *b_HLT_CscCluster_Tight;   //!
   TBranch        *b_HLT_CscCluster50_Photon20Unseeded;   //!
   TBranch        *b_HLT_CscCluster50_Photon30Unseeded;   //!
   TBranch        *b_HLT_CscCluster100_Ele5;   //!
   TBranch        *b_HLT_CscCluster100_Mu5;   //!
   TBranch        *b_HLT_CscCluster100_PNetTauhPFJet10_Loose;   //!
   TBranch        *b_HLT_DoubleCscCluster75;   //!
   TBranch        *b_HLT_IsoTrk200_L1SingleMuShower;   //!
   TBranch        *b_HLT_IsoTrk400_L1SingleMuShower;   //!
   TBranch        *b_HLT_DoubleCscCluster100;   //!
   TBranch        *b_HLT_L1CSCShower_DTCluster50;   //!
   TBranch        *b_HLT_L1CSCShower_DTCluster75;   //!
   TBranch        *b_HLT_PFMET105_IsoTrk50;   //!
   TBranch        *b_HLT_L1SingleLLPJet;   //!
   TBranch        *b_HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack;   //!
   TBranch        *b_HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5;   //!
   TBranch        *b_HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7;   //!
   TBranch        *b_HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8;   //!
   TBranch        *b_HLT_DiPhoton10Time1ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p2ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p4ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p6ns;   //!
   TBranch        *b_HLT_DiPhoton10Time1p8ns;   //!
   TBranch        *b_HLT_DiPhoton10Time2ns;   //!
   TBranch        *b_HLT_DiPhoton10_CaloIdL;   //!
   TBranch        *b_HLT_DoubleEle6p5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle8_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle10_eta1p22_mMax6;   //!
   TBranch        *b_HLT_SingleEle8;   //!
   TBranch        *b_HLT_SingleEle8_SingleEGL1;   //!
   TBranch        *b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT;   //!
   TBranch        *b_HLT_Mu50_L1SingleMuShower;   //!
   TBranch        *b_HLT_IsoMu24_OneProng32;   //!
   TBranch        *b_HLT_Photon32_OneProng32_M50To105;   //!
   TBranch        *b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1;   //!
   TBranch        *b_HLT_VBF_DiPFJet125_45_Mjj1050;   //!
   TBranch        *b_HLT_VBF_DiPFJet125_45_Mjj1200;   //!
   TBranch        *b_HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60;   //!
   TBranch        *b_HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60;   //!
   TBranch        *b_HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85;   //!
   TBranch        *b_HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85;   //!
   TBranch        *b_HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL;   //!
   TBranch        *b_HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL;   //!
   TBranch        *b_HLT_VBF_DiPFJet50_Mjj650_Photon22;   //!
   TBranch        *b_HLT_VBF_DiPFJet50_Mjj750_Photon22;   //!
   TBranch        *b_HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1;   //!
   TBranch        *b_HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1;   //!
   TBranch        *b_HLT_VBF_DoublePNetTauhPFJet20_eta2p2;   //!
   TBranch        *b_HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3;   //!
   TBranch        *b_HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3;   //!
   TBranch        *b_HLT_PFJet200_TimeLtNeg2p5ns;   //!
   TBranch        *b_HLT_PFJet200_TimeGt2p5ns;   //!
   TBranch        *b_HLT_Photon50_TimeLtNeg2p5ns;   //!
   TBranch        *b_HLT_Photon50_TimeGt2p5ns;   //!
   TBranch        *b_HLT_PPSMaxTracksPerArm1;   //!
   TBranch        *b_HLT_PPSMaxTracksPerRP4;   //!
   TBranch        *b_HLT_PPSRandom;   //!
   TBranch        *b_HLT_L1AXOVTight;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_Dataset_ScoutingPFMonitor;   //!
   TBranch        *b_Dataset_ScoutingPFRun3;   //!

   DepthsHistosFill(TTree *tree=0);
   virtual ~DepthsHistosFill();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DepthsHistosFill_cxx
DepthsHistosFill::DepthsHistosFill(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/dijet/data/NanoMerged_Data24Iv2_JetMET1_TeVJet_Baseline_RAWToPFNANO_v0p1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data/dijet/data/NanoMerged_Data24Iv2_JetMET1_TeVJet_Baseline_RAWToPFNANO_v0p1.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

DepthsHistosFill::~DepthsHistosFill()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DepthsHistosFill::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DepthsHistosFill::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DepthsHistosFill::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   /*
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("BeamSpot_type", &BeamSpot_type, &b_BeamSpot_type);
   fChain->SetBranchAddress("BeamSpot_sigmaZ", &BeamSpot_sigmaZ, &b_BeamSpot_sigmaZ);
   fChain->SetBranchAddress("BeamSpot_sigmaZError", &BeamSpot_sigmaZError, &b_BeamSpot_sigmaZError);
   fChain->SetBranchAddress("BeamSpot_z", &BeamSpot_z, &b_BeamSpot_z);
   fChain->SetBranchAddress("BeamSpot_zError", &BeamSpot_zError, &b_BeamSpot_zError);
   fChain->SetBranchAddress("nboostedTau", &nboostedTau, &b_nboostedTau);
   fChain->SetBranchAddress("boostedTau_idAntiEle2018", boostedTau_idAntiEle2018, &b_boostedTau_idAntiEle2018);
   fChain->SetBranchAddress("boostedTau_idAntiMu", boostedTau_idAntiMu, &b_boostedTau_idAntiMu);
   fChain->SetBranchAddress("boostedTau_idMVAnewDM2017v2", boostedTau_idMVAnewDM2017v2, &b_boostedTau_idMVAnewDM2017v2);
   fChain->SetBranchAddress("boostedTau_idMVAoldDM2017v2", boostedTau_idMVAoldDM2017v2, &b_boostedTau_idMVAoldDM2017v2);
   fChain->SetBranchAddress("boostedTau_jetIdx", boostedTau_jetIdx, &b_boostedTau_jetIdx);
   fChain->SetBranchAddress("boostedTau_rawAntiEleCat2018", boostedTau_rawAntiEleCat2018, &b_boostedTau_rawAntiEleCat2018);
   fChain->SetBranchAddress("boostedTau_charge", boostedTau_charge, &b_boostedTau_charge);
   fChain->SetBranchAddress("boostedTau_decayMode", boostedTau_decayMode, &b_boostedTau_decayMode);
   fChain->SetBranchAddress("boostedTau_chargedIso", boostedTau_chargedIso, &b_boostedTau_chargedIso);
   fChain->SetBranchAddress("boostedTau_eta", boostedTau_eta, &b_boostedTau_eta);
   fChain->SetBranchAddress("boostedTau_leadTkDeltaEta", boostedTau_leadTkDeltaEta, &b_boostedTau_leadTkDeltaEta);
   fChain->SetBranchAddress("boostedTau_leadTkDeltaPhi", boostedTau_leadTkDeltaPhi, &b_boostedTau_leadTkDeltaPhi);
   fChain->SetBranchAddress("boostedTau_leadTkPtOverTauPt", boostedTau_leadTkPtOverTauPt, &b_boostedTau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("boostedTau_mass", boostedTau_mass, &b_boostedTau_mass);
   fChain->SetBranchAddress("boostedTau_neutralIso", boostedTau_neutralIso, &b_boostedTau_neutralIso);
   fChain->SetBranchAddress("boostedTau_phi", boostedTau_phi, &b_boostedTau_phi);
   fChain->SetBranchAddress("boostedTau_photonsOutsideSignalCone", boostedTau_photonsOutsideSignalCone, &b_boostedTau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("boostedTau_pt", boostedTau_pt, &b_boostedTau_pt);
   fChain->SetBranchAddress("boostedTau_puCorr", boostedTau_puCorr, &b_boostedTau_puCorr);
   fChain->SetBranchAddress("boostedTau_rawAntiEle2018", boostedTau_rawAntiEle2018, &b_boostedTau_rawAntiEle2018);
   fChain->SetBranchAddress("boostedTau_rawIso", boostedTau_rawIso, &b_boostedTau_rawIso);
   fChain->SetBranchAddress("boostedTau_rawIsodR03", boostedTau_rawIsodR03, &b_boostedTau_rawIsodR03);
   fChain->SetBranchAddress("boostedTau_rawMVAnewDM2017v2", boostedTau_rawMVAnewDM2017v2, &b_boostedTau_rawMVAnewDM2017v2);
   fChain->SetBranchAddress("boostedTau_rawMVAoldDM2017v2", boostedTau_rawMVAoldDM2017v2, &b_boostedTau_rawMVAoldDM2017v2);
   fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   fChain->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet, &b_nCorrT1METJet);
   fChain->SetBranchAddress("CorrT1METJet_EmEF", CorrT1METJet_EmEF, &b_CorrT1METJet_EmEF);
   fChain->SetBranchAddress("CorrT1METJet_area", CorrT1METJet_area, &b_CorrT1METJet_area);
   fChain->SetBranchAddress("CorrT1METJet_eta", CorrT1METJet_eta, &b_CorrT1METJet_eta);
   fChain->SetBranchAddress("CorrT1METJet_muonSubtrFactor", CorrT1METJet_muonSubtrFactor, &b_CorrT1METJet_muonSubtrFactor);
   fChain->SetBranchAddress("CorrT1METJet_phi", CorrT1METJet_phi, &b_CorrT1METJet_phi);
   fChain->SetBranchAddress("CorrT1METJet_rawPt", CorrT1METJet_rawPt, &b_CorrT1METJet_rawPt);
   fChain->SetBranchAddress("DeepMETResolutionTune_phi", &DeepMETResolutionTune_phi, &b_DeepMETResolutionTune_phi);
   fChain->SetBranchAddress("DeepMETResolutionTune_pt", &DeepMETResolutionTune_pt, &b_DeepMETResolutionTune_pt);
   fChain->SetBranchAddress("DeepMETResponseTune_phi", &DeepMETResponseTune_phi, &b_DeepMETResponseTune_phi);
   fChain->SetBranchAddress("DeepMETResponseTune_pt", &DeepMETResponseTune_pt, &b_DeepMETResponseTune_pt);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isEB", Electron_isEB, &b_Electron_isEB);
   fChain->SetBranchAddress("Electron_isEcalDriven", Electron_isEcalDriven, &b_Electron_isEcalDriven);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_jetNDauCharged", Electron_jetNDauCharged, &b_Electron_jetNDauCharged);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaIso_WP80", Electron_mvaIso_WP80, &b_Electron_mvaIso_WP80);
   fChain->SetBranchAddress("Electron_mvaIso_WP90", Electron_mvaIso_WP90, &b_Electron_mvaIso_WP90);
   fChain->SetBranchAddress("Electron_mvaIso_WPHZZ", Electron_mvaIso_WPHZZ, &b_Electron_mvaIso_WPHZZ);
   fChain->SetBranchAddress("Electron_mvaNoIso_WP80", Electron_mvaNoIso_WP80, &b_Electron_mvaNoIso_WP80);
   fChain->SetBranchAddress("Electron_mvaNoIso_WP90", Electron_mvaNoIso_WP90, &b_Electron_mvaNoIso_WP90);
   fChain->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_seediEtaOriX", Electron_seediEtaOriX, &b_Electron_seediEtaOriX);
   fChain->SetBranchAddress("Electron_seediPhiOriY", Electron_seediPhiOriY, &b_Electron_seediPhiOriY);
   fChain->SetBranchAddress("Electron_svIdx", Electron_svIdx, &b_Electron_svIdx);
   fChain->SetBranchAddress("Electron_fsrPhotonIdx", Electron_fsrPhotonIdx, &b_Electron_fsrPhotonIdx);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmapHEEP", Electron_vidNestedWPBitmapHEEP, &b_Electron_vidNestedWPBitmapHEEP);
   fChain->SetBranchAddress("Electron_PreshowerEnergy", Electron_PreshowerEnergy, &b_Electron_PreshowerEnergy);
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_ecalEnergy", Electron_ecalEnergy, &b_Electron_ecalEnergy);
   fChain->SetBranchAddress("Electron_ecalEnergyError", Electron_ecalEnergyError, &b_Electron_ecalEnergyError);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_fbrem", Electron_fbrem, &b_Electron_fbrem);
   fChain->SetBranchAddress("Electron_gsfTrketaMode", Electron_gsfTrketaMode, &b_Electron_gsfTrketaMode);
   fChain->SetBranchAddress("Electron_gsfTrkpMode", Electron_gsfTrkpMode, &b_Electron_gsfTrkpMode);
   fChain->SetBranchAddress("Electron_gsfTrkpModeErr", Electron_gsfTrkpModeErr, &b_Electron_gsfTrkpModeErr);
   fChain->SetBranchAddress("Electron_gsfTrkphiMode", Electron_gsfTrkphiMode, &b_Electron_gsfTrkphiMode);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2, &b_Electron_jetPtRelv2);
   fChain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaHZZIso", Electron_mvaHZZIso, &b_Electron_mvaHZZIso);
   fChain->SetBranchAddress("Electron_mvaIso", Electron_mvaIso, &b_Electron_mvaIso);
   fChain->SetBranchAddress("Electron_mvaNoIso", Electron_mvaNoIso, &b_Electron_mvaNoIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_rawEnergy", Electron_rawEnergy, &b_Electron_rawEnergy);
   fChain->SetBranchAddress("Electron_scEtOverPt", Electron_scEtOverPt, &b_Electron_scEtOverPt);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_promptMVA", Electron_promptMVA, &b_Electron_promptMVA);
   fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   fChain->SetBranchAddress("FatJet_nConstituents", FatJet_nConstituents, &b_FatJet_nConstituents);
   fChain->SetBranchAddress("FatJet_chMultiplicity", FatJet_chMultiplicity, &b_FatJet_chMultiplicity);
   fChain->SetBranchAddress("FatJet_neMultiplicity", FatJet_neMultiplicity, &b_FatJet_neMultiplicity);
   fChain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   fChain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   fChain->SetBranchAddress("FatJet_electronIdx3SJ", FatJet_electronIdx3SJ, &b_FatJet_electronIdx3SJ);
   fChain->SetBranchAddress("FatJet_muonIdx3SJ", FatJet_muonIdx3SJ, &b_FatJet_muonIdx3SJ);
   fChain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   fChain->SetBranchAddress("FatJet_chEmEF", FatJet_chEmEF, &b_FatJet_chEmEF);
   fChain->SetBranchAddress("FatJet_chHEF", FatJet_chHEF, &b_FatJet_chHEF);
   fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   fChain->SetBranchAddress("FatJet_muEF", FatJet_muEF, &b_FatJet_muEF);
   fChain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   fChain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   fChain->SetBranchAddress("FatJet_neEmEF", FatJet_neEmEF, &b_FatJet_neEmEF);
   fChain->SetBranchAddress("FatJet_neHEF", FatJet_neHEF, &b_FatJet_neHEF);
   fChain->SetBranchAddress("FatJet_particleNetLegacy_QCD", FatJet_particleNetLegacy_QCD, &b_FatJet_particleNetLegacy_QCD);
   fChain->SetBranchAddress("FatJet_particleNetLegacy_Xbb", FatJet_particleNetLegacy_Xbb, &b_FatJet_particleNetLegacy_Xbb);
   fChain->SetBranchAddress("FatJet_particleNetLegacy_Xcc", FatJet_particleNetLegacy_Xcc, &b_FatJet_particleNetLegacy_Xcc);
   fChain->SetBranchAddress("FatJet_particleNetLegacy_Xqq", FatJet_particleNetLegacy_Xqq, &b_FatJet_particleNetLegacy_Xqq);
   fChain->SetBranchAddress("FatJet_particleNetLegacy_mass", FatJet_particleNetLegacy_mass, &b_FatJet_particleNetLegacy_mass);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_H4qvsQCD", FatJet_particleNetWithMass_H4qvsQCD, &b_FatJet_particleNetWithMass_H4qvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_HbbvsQCD", FatJet_particleNetWithMass_HbbvsQCD, &b_FatJet_particleNetWithMass_HbbvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_HccvsQCD", FatJet_particleNetWithMass_HccvsQCD, &b_FatJet_particleNetWithMass_HccvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_QCD", FatJet_particleNetWithMass_QCD, &b_FatJet_particleNetWithMass_QCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_TvsQCD", FatJet_particleNetWithMass_TvsQCD, &b_FatJet_particleNetWithMass_TvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_WvsQCD", FatJet_particleNetWithMass_WvsQCD, &b_FatJet_particleNetWithMass_WvsQCD);
   fChain->SetBranchAddress("FatJet_particleNetWithMass_ZvsQCD", FatJet_particleNetWithMass_ZvsQCD, &b_FatJet_particleNetWithMass_ZvsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_QCD", FatJet_particleNet_QCD, &b_FatJet_particleNet_QCD);
   fChain->SetBranchAddress("FatJet_particleNet_QCD0HF", FatJet_particleNet_QCD0HF, &b_FatJet_particleNet_QCD0HF);
   fChain->SetBranchAddress("FatJet_particleNet_QCD1HF", FatJet_particleNet_QCD1HF, &b_FatJet_particleNet_QCD1HF);
   fChain->SetBranchAddress("FatJet_particleNet_QCD2HF", FatJet_particleNet_QCD2HF, &b_FatJet_particleNet_QCD2HF);
   fChain->SetBranchAddress("FatJet_particleNet_XbbVsQCD", FatJet_particleNet_XbbVsQCD, &b_FatJet_particleNet_XbbVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XccVsQCD", FatJet_particleNet_XccVsQCD, &b_FatJet_particleNet_XccVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XggVsQCD", FatJet_particleNet_XggVsQCD, &b_FatJet_particleNet_XggVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XqqVsQCD", FatJet_particleNet_XqqVsQCD, &b_FatJet_particleNet_XqqVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XteVsQCD", FatJet_particleNet_XteVsQCD, &b_FatJet_particleNet_XteVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XtmVsQCD", FatJet_particleNet_XtmVsQCD, &b_FatJet_particleNet_XtmVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_XttVsQCD", FatJet_particleNet_XttVsQCD, &b_FatJet_particleNet_XttVsQCD);
   fChain->SetBranchAddress("FatJet_particleNet_massCorr", FatJet_particleNet_massCorr, &b_FatJet_particleNet_massCorr);
   fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   fChain->SetBranchAddress("FatJet_lsf3", FatJet_lsf3, &b_FatJet_lsf3);
   */
   fChain->SetBranchAddress("nJetPFCand", &nJetPFCand, &b_nJetPFCand);
   fChain->SetBranchAddress("JetPFCand_pfCandIdx", JetPFCand_pfCandIdx, &b_JetPFCand_pfCandIdx);
   fChain->SetBranchAddress("JetPFCand_jetIdx", JetPFCand_jetIdx, &b_JetPFCand_jetIdx);
   /*
   fChain->SetBranchAddress("nFatJetPFCand", &nFatJetPFCand, &b_nFatJetPFCand);
   fChain->SetBranchAddress("FatJetPFCand_pfCandIdx", FatJetPFCand_pfCandIdx, &b_FatJetPFCand_pfCandIdx);
   fChain->SetBranchAddress("FatJetPFCand_jetIdx", FatJetPFCand_jetIdx, &b_FatJetPFCand_jetIdx);
   fChain->SetBranchAddress("nFsrPhoton", &nFsrPhoton, &b_nFsrPhoton);
   fChain->SetBranchAddress("FsrPhoton_electronIdx", FsrPhoton_electronIdx, &b_FsrPhoton_electronIdx);
   fChain->SetBranchAddress("FsrPhoton_muonIdx", FsrPhoton_muonIdx, &b_FsrPhoton_muonIdx);
   fChain->SetBranchAddress("FsrPhoton_dROverEt2", FsrPhoton_dROverEt2, &b_FsrPhoton_dROverEt2);
   fChain->SetBranchAddress("FsrPhoton_eta", FsrPhoton_eta, &b_FsrPhoton_eta);
   fChain->SetBranchAddress("FsrPhoton_phi", FsrPhoton_phi, &b_FsrPhoton_phi);
   fChain->SetBranchAddress("FsrPhoton_pt", FsrPhoton_pt, &b_FsrPhoton_pt);
   fChain->SetBranchAddress("FsrPhoton_relIso03", FsrPhoton_relIso03, &b_FsrPhoton_relIso03);
   fChain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   fChain->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   fChain->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   fChain->SetBranchAddress("IsoTrack_charge", IsoTrack_charge, &b_IsoTrack_charge);
   fChain->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   fChain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   */
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_chMultiplicity", Jet_chMultiplicity, &b_Jet_chMultiplicity);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_nSVs", Jet_nSVs, &b_Jet_nSVs);
   fChain->SetBranchAddress("Jet_neMultiplicity", Jet_neMultiplicity, &b_Jet_neMultiplicity);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_svIdx1", Jet_svIdx1, &b_Jet_svIdx1);
   fChain->SetBranchAddress("Jet_svIdx2", Jet_svIdx2, &b_Jet_svIdx2);
   fChain->SetBranchAddress("Jet_hfadjacentEtaStripsSize", Jet_hfadjacentEtaStripsSize, &b_Jet_hfadjacentEtaStripsSize);
   fChain->SetBranchAddress("Jet_hfcentralEtaStripSize", Jet_hfcentralEtaStripSize, &b_Jet_hfcentralEtaStripSize);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorr", Jet_PNetRegPtRawCorr, &b_Jet_PNetRegPtRawCorr);
   fChain->SetBranchAddress("Jet_PNetRegPtRawCorrNeutrino", Jet_PNetRegPtRawCorrNeutrino, &b_Jet_PNetRegPtRawCorrNeutrino);
   fChain->SetBranchAddress("Jet_PNetRegPtRawRes", Jet_PNetRegPtRawRes, &b_Jet_PNetRegPtRawRes);
   fChain->SetBranchAddress("Jet_UParTAK4RegPtRawCorr", Jet_UParTAK4RegPtRawCorr, &b_Jet_UParTAK4RegPtRawCorr);
   fChain->SetBranchAddress("Jet_UParTAK4RegPtRawCorrNeutrino", Jet_UParTAK4RegPtRawCorrNeutrino, &b_Jet_UParTAK4RegPtRawCorrNeutrino);
   fChain->SetBranchAddress("Jet_UParTAK4RegPtRawRes", Jet_UParTAK4RegPtRawRes, &b_Jet_UParTAK4RegPtRawRes);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB, &b_Jet_btagDeepFlavCvB);
   fChain->SetBranchAddress("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL, &b_Jet_btagDeepFlavCvL);
   fChain->SetBranchAddress("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG, &b_Jet_btagDeepFlavQG);
   fChain->SetBranchAddress("Jet_btagPNetB", Jet_btagPNetB, &b_Jet_btagPNetB);
   fChain->SetBranchAddress("Jet_btagPNetCvB", Jet_btagPNetCvB, &b_Jet_btagPNetCvB);
   fChain->SetBranchAddress("Jet_btagPNetCvL", Jet_btagPNetCvL, &b_Jet_btagPNetCvL);
   fChain->SetBranchAddress("Jet_btagPNetCvNotB", Jet_btagPNetCvNotB, &b_Jet_btagPNetCvNotB);
   fChain->SetBranchAddress("Jet_btagPNetQvG", Jet_btagPNetQvG, &b_Jet_btagPNetQvG);
   fChain->SetBranchAddress("Jet_btagPNetTauVJet", Jet_btagPNetTauVJet, &b_Jet_btagPNetTauVJet);
   fChain->SetBranchAddress("Jet_btagUParTAK4B", Jet_btagUParTAK4B, &b_Jet_btagUParTAK4B);
   fChain->SetBranchAddress("Jet_btagUParTAK4CvB", Jet_btagUParTAK4CvB, &b_Jet_btagUParTAK4CvB);
   fChain->SetBranchAddress("Jet_btagUParTAK4CvL", Jet_btagUParTAK4CvL, &b_Jet_btagUParTAK4CvL);
   fChain->SetBranchAddress("Jet_btagUParTAK4CvNotB", Jet_btagUParTAK4CvNotB, &b_Jet_btagUParTAK4CvNotB);
   fChain->SetBranchAddress("Jet_btagUParTAK4QvG", Jet_btagUParTAK4QvG, &b_Jet_btagUParTAK4QvG);
   fChain->SetBranchAddress("Jet_btagUParTAK4TauVJet", Jet_btagUParTAK4TauVJet, &b_Jet_btagUParTAK4TauVJet);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_hfEmEF", Jet_hfEmEF, &b_Jet_hfEmEF);
   fChain->SetBranchAddress("Jet_hfHEF", Jet_hfHEF, &b_Jet_hfHEF);
   fChain->SetBranchAddress("Jet_hfsigmaEtaEta", Jet_hfsigmaEtaEta, &b_Jet_hfsigmaEtaEta);
   fChain->SetBranchAddress("Jet_hfsigmaPhiPhi", Jet_hfsigmaPhiPhi, &b_Jet_hfsigmaPhiPhi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   /*
   fChain->SetBranchAddress("nLowPtElectron", &nLowPtElectron, &b_nLowPtElectron);
   fChain->SetBranchAddress("LowPtElectron_convVeto", LowPtElectron_convVeto, &b_LowPtElectron_convVeto);
   fChain->SetBranchAddress("LowPtElectron_convWP", LowPtElectron_convWP, &b_LowPtElectron_convWP);
   fChain->SetBranchAddress("LowPtElectron_lostHits", LowPtElectron_lostHits, &b_LowPtElectron_lostHits);
   fChain->SetBranchAddress("LowPtElectron_electronIdx", LowPtElectron_electronIdx, &b_LowPtElectron_electronIdx);
   fChain->SetBranchAddress("LowPtElectron_photonIdx", LowPtElectron_photonIdx, &b_LowPtElectron_photonIdx);
   fChain->SetBranchAddress("LowPtElectron_charge", LowPtElectron_charge, &b_LowPtElectron_charge);
   fChain->SetBranchAddress("LowPtElectron_pdgId", LowPtElectron_pdgId, &b_LowPtElectron_pdgId);
   fChain->SetBranchAddress("LowPtElectron_ID", LowPtElectron_ID, &b_LowPtElectron_ID);
   fChain->SetBranchAddress("LowPtElectron_convVtxRadius", LowPtElectron_convVtxRadius, &b_LowPtElectron_convVtxRadius);
   fChain->SetBranchAddress("LowPtElectron_deltaEtaSC", LowPtElectron_deltaEtaSC, &b_LowPtElectron_deltaEtaSC);
   fChain->SetBranchAddress("LowPtElectron_dxy", LowPtElectron_dxy, &b_LowPtElectron_dxy);
   fChain->SetBranchAddress("LowPtElectron_dxyErr", LowPtElectron_dxyErr, &b_LowPtElectron_dxyErr);
   fChain->SetBranchAddress("LowPtElectron_dz", LowPtElectron_dz, &b_LowPtElectron_dz);
   fChain->SetBranchAddress("LowPtElectron_dzErr", LowPtElectron_dzErr, &b_LowPtElectron_dzErr);
   fChain->SetBranchAddress("LowPtElectron_eInvMinusPInv", LowPtElectron_eInvMinusPInv, &b_LowPtElectron_eInvMinusPInv);
   fChain->SetBranchAddress("LowPtElectron_energyErr", LowPtElectron_energyErr, &b_LowPtElectron_energyErr);
   fChain->SetBranchAddress("LowPtElectron_eta", LowPtElectron_eta, &b_LowPtElectron_eta);
   fChain->SetBranchAddress("LowPtElectron_hoe", LowPtElectron_hoe, &b_LowPtElectron_hoe);
   fChain->SetBranchAddress("LowPtElectron_mass", LowPtElectron_mass, &b_LowPtElectron_mass);
   fChain->SetBranchAddress("LowPtElectron_miniPFRelIso_all", LowPtElectron_miniPFRelIso_all, &b_LowPtElectron_miniPFRelIso_all);
   fChain->SetBranchAddress("LowPtElectron_miniPFRelIso_chg", LowPtElectron_miniPFRelIso_chg, &b_LowPtElectron_miniPFRelIso_chg);
   fChain->SetBranchAddress("LowPtElectron_phi", LowPtElectron_phi, &b_LowPtElectron_phi);
   fChain->SetBranchAddress("LowPtElectron_pt", LowPtElectron_pt, &b_LowPtElectron_pt);
   fChain->SetBranchAddress("LowPtElectron_ptbiased", LowPtElectron_ptbiased, &b_LowPtElectron_ptbiased);
   fChain->SetBranchAddress("LowPtElectron_r9", LowPtElectron_r9, &b_LowPtElectron_r9);
   fChain->SetBranchAddress("LowPtElectron_scEtOverPt", LowPtElectron_scEtOverPt, &b_LowPtElectron_scEtOverPt);
   fChain->SetBranchAddress("LowPtElectron_sieie", LowPtElectron_sieie, &b_LowPtElectron_sieie);
   fChain->SetBranchAddress("LowPtElectron_unbiased", LowPtElectron_unbiased, &b_LowPtElectron_unbiased);
   fChain->SetBranchAddress("nProton_multiRP", &nProton_multiRP, &b_nProton_multiRP);
   fChain->SetBranchAddress("Proton_multiRP_arm", Proton_multiRP_arm, &b_Proton_multiRP_arm);
   fChain->SetBranchAddress("Proton_multiRP_t", Proton_multiRP_t, &b_Proton_multiRP_t);
   fChain->SetBranchAddress("Proton_multiRP_thetaX", Proton_multiRP_thetaX, &b_Proton_multiRP_thetaX);
   fChain->SetBranchAddress("Proton_multiRP_thetaY", Proton_multiRP_thetaY, &b_Proton_multiRP_thetaY);
   fChain->SetBranchAddress("Proton_multiRP_time", Proton_multiRP_time, &b_Proton_multiRP_time);
   fChain->SetBranchAddress("Proton_multiRP_timeUnc", Proton_multiRP_timeUnc, &b_Proton_multiRP_timeUnc);
   fChain->SetBranchAddress("Proton_multiRP_xi", Proton_multiRP_xi, &b_Proton_multiRP_xi);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_highPurity", Muon_highPurity, &b_Muon_highPurity);
   fChain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isStandalone", Muon_isStandalone, &b_Muon_isStandalone);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_jetNDauCharged", Muon_jetNDauCharged, &b_Muon_jetNDauCharged);
   fChain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaMuID_WP", Muon_mvaMuID_WP, &b_Muon_mvaMuID_WP);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_puppiIsoId", Muon_puppiIsoId, &b_Muon_puppiIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_svIdx", Muon_svIdx, &b_Muon_svIdx);
   fChain->SetBranchAddress("Muon_fsrPhotonIdx", Muon_fsrPhotonIdx, &b_Muon_fsrPhotonIdx);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dxybs", Muon_dxybs, &b_Muon_dxybs);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   fChain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_mvaMuID", Muon_mvaMuID, &b_Muon_mvaMuID);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_softMva", Muon_softMva, &b_Muon_softMva);
   fChain->SetBranchAddress("Muon_softMvaRun3", Muon_softMvaRun3, &b_Muon_softMvaRun3);
   fChain->SetBranchAddress("Muon_tkRelIso", Muon_tkRelIso, &b_Muon_tkRelIso);
   fChain->SetBranchAddress("Muon_tunepRelPt", Muon_tunepRelPt, &b_Muon_tunepRelPt);
   fChain->SetBranchAddress("Muon_bsConstrainedChi2", Muon_bsConstrainedChi2, &b_Muon_bsConstrainedChi2);
   fChain->SetBranchAddress("Muon_bsConstrainedPt", Muon_bsConstrainedPt, &b_Muon_bsConstrainedPt);
   fChain->SetBranchAddress("Muon_bsConstrainedPtErr", Muon_bsConstrainedPtErr, &b_Muon_bsConstrainedPtErr);
   fChain->SetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt, &b_Muon_mvaLowPt);
   fChain->SetBranchAddress("Muon_promptMVA", Muon_promptMVA, &b_Muon_promptMVA);
   fChain->SetBranchAddress("nPFCandToPFClusterHCAL", &nPFCandToPFClusterHCAL, &b_nPFCandToPFClusterHCAL);
   fChain->SetBranchAddress("PFCandToPFClusterHCAL_PFCandIdx", PFCandToPFClusterHCAL_PFCandIdx, &b_PFCandToPFClusterHCAL_PFCandIdx);
   fChain->SetBranchAddress("PFCandToPFClusterHCAL_PFClusterHCALIdx", PFCandToPFClusterHCAL_PFClusterHCALIdx, &b_PFCandToPFClusterHCAL_PFClusterHCALIdx);
   fChain->SetBranchAddress("nPFClusterECAL", &nPFClusterECAL, &b_nPFClusterECAL);
   fChain->SetBranchAddress("PFClusterECAL_layer", PFClusterECAL_layer, &b_PFClusterECAL_layer);
   fChain->SetBranchAddress("PFClusterECAL_nhits", PFClusterECAL_nhits, &b_PFClusterECAL_nhits);
   fChain->SetBranchAddress("PFClusterECAL_seedhit_detId", PFClusterECAL_seedhit_detId, &b_PFClusterECAL_seedhit_detId);
   fChain->SetBranchAddress("PFClusterECAL_pt", PFClusterECAL_pt, &b_PFClusterECAL_pt);
   fChain->SetBranchAddress("PFClusterECAL_energy", PFClusterECAL_energy, &b_PFClusterECAL_energy);
   fChain->SetBranchAddress("PFClusterECAL_correctedEnergy", PFClusterECAL_correctedEnergy, &b_PFClusterECAL_correctedEnergy);
   fChain->SetBranchAddress("PFClusterECAL_eta", PFClusterECAL_eta, &b_PFClusterECAL_eta);
   fChain->SetBranchAddress("PFClusterECAL_phi", PFClusterECAL_phi, &b_PFClusterECAL_phi);
   fChain->SetBranchAddress("nPFClusterHCAL", &nPFClusterHCAL, &b_nPFClusterHCAL);
   fChain->SetBranchAddress("PFClusterHCAL_layer", PFClusterHCAL_layer, &b_PFClusterHCAL_layer);
   fChain->SetBranchAddress("PFClusterHCAL_nhits", PFClusterHCAL_nhits, &b_PFClusterHCAL_nhits);
   fChain->SetBranchAddress("PFClusterHCAL_seedhit_detId", PFClusterHCAL_seedhit_detId, &b_PFClusterHCAL_seedhit_detId);
   fChain->SetBranchAddress("PFClusterHCAL_pt", PFClusterHCAL_pt, &b_PFClusterHCAL_pt);
   fChain->SetBranchAddress("PFClusterHCAL_energy", PFClusterHCAL_energy, &b_PFClusterHCAL_energy);
   fChain->SetBranchAddress("PFClusterHCAL_eta", PFClusterHCAL_eta, &b_PFClusterHCAL_eta);
   fChain->SetBranchAddress("PFClusterHCAL_phi", PFClusterHCAL_phi, &b_PFClusterHCAL_phi);
   fChain->SetBranchAddress("nPFClusterHCALToPFRecHitHBHE", &nPFClusterHCALToPFRecHitHBHE, &b_nPFClusterHCALToPFRecHitHBHE);
   fChain->SetBranchAddress("PFClusterHCALToPFRecHitHBHE_PFClusterHCALIdx", PFClusterHCALToPFRecHitHBHE_PFClusterHCALIdx, &b_PFClusterHCALToPFRecHitHBHE_PFClusterHCALIdx);
   fChain->SetBranchAddress("PFClusterHCALToPFRecHitHBHE_PFRecHitHBHEIdx", PFClusterHCALToPFRecHitHBHE_PFRecHitHBHEIdx, &b_PFClusterHCALToPFRecHitHBHE_PFRecHitHBHEIdx);
   fChain->SetBranchAddress("PFClusterHCALToPFRecHitHBHE_fraction", PFClusterHCALToPFRecHitHBHE_fraction, &b_PFClusterHCALToPFRecHitHBHE_fraction);
   fChain->SetBranchAddress("nPFClusterPS", &nPFClusterPS, &b_nPFClusterPS);
   fChain->SetBranchAddress("PFClusterPS_layer", PFClusterPS_layer, &b_PFClusterPS_layer);
   fChain->SetBranchAddress("PFClusterPS_nhits", PFClusterPS_nhits, &b_PFClusterPS_nhits);
   fChain->SetBranchAddress("PFClusterPS_seedhit_detId", PFClusterPS_seedhit_detId, &b_PFClusterPS_seedhit_detId);
   fChain->SetBranchAddress("PFClusterPS_pt", PFClusterPS_pt, &b_PFClusterPS_pt);
   fChain->SetBranchAddress("PFClusterPS_energy", PFClusterPS_energy, &b_PFClusterPS_energy);
   fChain->SetBranchAddress("PFClusterPS_correctedEnergy", PFClusterPS_correctedEnergy, &b_PFClusterPS_correctedEnergy);
   fChain->SetBranchAddress("PFClusterPS_eta", PFClusterPS_eta, &b_PFClusterPS_eta);
   fChain->SetBranchAddress("PFClusterPS_phi", PFClusterPS_phi, &b_PFClusterPS_phi);
   fChain->SetBranchAddress("nPFRecHitHBHE", &nPFRecHitHBHE, &b_nPFRecHitHBHE);
   fChain->SetBranchAddress("PFRecHitHBHE_ieta", PFRecHitHBHE_ieta, &b_PFRecHitHBHE_ieta);
   fChain->SetBranchAddress("PFRecHitHBHE_iphi", PFRecHitHBHE_iphi, &b_PFRecHitHBHE_iphi);
   fChain->SetBranchAddress("PFRecHitHBHE_depth", PFRecHitHBHE_depth, &b_PFRecHitHBHE_depth);
   fChain->SetBranchAddress("PFRecHitHBHE_detId", PFRecHitHBHE_detId, &b_PFRecHitHBHE_detId);
   fChain->SetBranchAddress("PFRecHitHBHE_energy", PFRecHitHBHE_energy, &b_PFRecHitHBHE_energy);
   */
   fChain->SetBranchAddress("nPFCand", &nPFCand, &b_nPFCand);
   fChain->SetBranchAddress("PFCand_isIsolatedChargedHadron", PFCand_isIsolatedChargedHadron, &b_PFCand_isIsolatedChargedHadron);
   fChain->SetBranchAddress("PFCand_trkHighPurity", PFCand_trkHighPurity, &b_PFCand_trkHighPurity);
   fChain->SetBranchAddress("PFCand_charge", PFCand_charge, &b_PFCand_charge);
   fChain->SetBranchAddress("PFCand_fromPV0", PFCand_fromPV0, &b_PFCand_fromPV0);
   fChain->SetBranchAddress("PFCand_lostInnerHits", PFCand_lostInnerHits, &b_PFCand_lostInnerHits);
   fChain->SetBranchAddress("PFCand_lostOuterHits", PFCand_lostOuterHits, &b_PFCand_lostOuterHits);
   fChain->SetBranchAddress("PFCand_nHits", PFCand_nHits, &b_PFCand_nHits);
   fChain->SetBranchAddress("PFCand_nPixelHits", PFCand_nPixelHits, &b_PFCand_nPixelHits);
   fChain->SetBranchAddress("PFCand_pdgId", PFCand_pdgId, &b_PFCand_pdgId);
   fChain->SetBranchAddress("PFCand_pvAssocQuality", PFCand_pvAssocQuality, &b_PFCand_pvAssocQuality);
   fChain->SetBranchAddress("PFCand_trkAlgo", PFCand_trkAlgo, &b_PFCand_trkAlgo);
   fChain->SetBranchAddress("PFCand_trkQuality", PFCand_trkQuality, &b_PFCand_trkQuality);
   fChain->SetBranchAddress("PFCand_vertexRef", PFCand_vertexRef, &b_PFCand_vertexRef);
   fChain->SetBranchAddress("PFCand_d0", PFCand_d0, &b_PFCand_d0);
   fChain->SetBranchAddress("PFCand_d0Err", PFCand_d0Err, &b_PFCand_d0Err);
   fChain->SetBranchAddress("PFCand_dz", PFCand_dz, &b_PFCand_dz);
   fChain->SetBranchAddress("PFCand_dzErr", PFCand_dzErr, &b_PFCand_dzErr);
   fChain->SetBranchAddress("PFCand_energy", PFCand_energy, &b_PFCand_energy);
   fChain->SetBranchAddress("PFCand_eta", PFCand_eta, &b_PFCand_eta);
   fChain->SetBranchAddress("PFCand_mass", PFCand_mass, &b_PFCand_mass);
   fChain->SetBranchAddress("PFCand_phi", PFCand_phi, &b_PFCand_phi);
   fChain->SetBranchAddress("PFCand_pt", PFCand_pt, &b_PFCand_pt);
   fChain->SetBranchAddress("PFCand_puppiWeight", PFCand_puppiWeight, &b_PFCand_puppiWeight);
   fChain->SetBranchAddress("PFCand_puppiWeightNoLep", PFCand_puppiWeightNoLep, &b_PFCand_puppiWeightNoLep);
   fChain->SetBranchAddress("PFCand_trkChi2", PFCand_trkChi2, &b_PFCand_trkChi2);
   fChain->SetBranchAddress("PFCand_trkEta", PFCand_trkEta, &b_PFCand_trkEta);
   fChain->SetBranchAddress("PFCand_trkP", PFCand_trkP, &b_PFCand_trkP);
   fChain->SetBranchAddress("PFCand_trkPhi", PFCand_trkPhi, &b_PFCand_trkPhi);
   fChain->SetBranchAddress("PFCand_trkPt", PFCand_trkPt, &b_PFCand_trkPt);
   fChain->SetBranchAddress("PFCand_vtxChi2", PFCand_vtxChi2, &b_PFCand_vtxChi2);
   /*
   fChain->SetBranchAddress("PFMET_covXX", &PFMET_covXX, &b_PFMET_covXX);
   fChain->SetBranchAddress("PFMET_covXY", &PFMET_covXY, &b_PFMET_covXY);
   fChain->SetBranchAddress("PFMET_covYY", &PFMET_covYY, &b_PFMET_covYY);
   fChain->SetBranchAddress("PFMET_phi", &PFMET_phi, &b_PFMET_phi);
   fChain->SetBranchAddress("PFMET_phiUnclusteredDown", &PFMET_phiUnclusteredDown, &b_PFMET_phiUnclusteredDown);
   fChain->SetBranchAddress("PFMET_phiUnclusteredUp", &PFMET_phiUnclusteredUp, &b_PFMET_phiUnclusteredUp);
   fChain->SetBranchAddress("PFMET_pt", &PFMET_pt, &b_PFMET_pt);
   fChain->SetBranchAddress("PFMET_ptUnclusteredDown", &PFMET_ptUnclusteredDown, &b_PFMET_ptUnclusteredDown);
   fChain->SetBranchAddress("PFMET_ptUnclusteredUp", &PFMET_ptUnclusteredUp, &b_PFMET_ptUnclusteredUp);
   fChain->SetBranchAddress("PFMET_significance", &PFMET_significance, &b_PFMET_significance);
   fChain->SetBranchAddress("PFMET_sumEt", &PFMET_sumEt, &b_PFMET_sumEt);
   fChain->SetBranchAddress("PFMET_sumPtUnclustered", &PFMET_sumPtUnclustered, &b_PFMET_sumPtUnclustered);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_hasConversionTracks", Photon_hasConversionTracks, &b_Photon_hasConversionTracks);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_seediEtaOriX", Photon_seediEtaOriX, &b_Photon_seediEtaOriX);
   fChain->SetBranchAddress("Photon_seediPhiOriY", Photon_seediPhiOriY, &b_Photon_seediPhiOriY);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_ecalPFClusterIso", Photon_ecalPFClusterIso, &b_Photon_ecalPFClusterIso);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_energyRaw", Photon_energyRaw, &b_Photon_energyRaw);
   fChain->SetBranchAddress("Photon_esEffSigmaRR", Photon_esEffSigmaRR, &b_Photon_esEffSigmaRR);
   fChain->SetBranchAddress("Photon_esEnergyOverRawE", Photon_esEnergyOverRawE, &b_Photon_esEnergyOverRawE);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
   fChain->SetBranchAddress("Photon_haloTaggerMVAVal", Photon_haloTaggerMVAVal, &b_Photon_haloTaggerMVAVal);
   fChain->SetBranchAddress("Photon_hcalPFClusterIso", Photon_hcalPFClusterIso, &b_Photon_hcalPFClusterIso);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_hoe_PUcorr", Photon_hoe_PUcorr, &b_Photon_hoe_PUcorr);
   fChain->SetBranchAddress("Photon_hoe_Tower", Photon_hoe_Tower, &b_Photon_hoe_Tower);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_pfChargedIso", Photon_pfChargedIso, &b_Photon_pfChargedIso);
   fChain->SetBranchAddress("Photon_pfChargedIsoPFPV", Photon_pfChargedIsoPFPV, &b_Photon_pfChargedIsoPFPV);
   fChain->SetBranchAddress("Photon_pfChargedIsoWorstVtx", Photon_pfChargedIsoWorstVtx, &b_Photon_pfChargedIsoWorstVtx);
   fChain->SetBranchAddress("Photon_pfPhoIso03", Photon_pfPhoIso03, &b_Photon_pfPhoIso03);
   fChain->SetBranchAddress("Photon_pfRelIso03_all_quadratic", Photon_pfRelIso03_all_quadratic, &b_Photon_pfRelIso03_all_quadratic);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg_quadratic", Photon_pfRelIso03_chg_quadratic, &b_Photon_pfRelIso03_chg_quadratic);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_s4", Photon_s4, &b_Photon_s4);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_sieip", Photon_sieip, &b_Photon_sieip);
   fChain->SetBranchAddress("Photon_sipip", Photon_sipip, &b_Photon_sipip);
   fChain->SetBranchAddress("Photon_superclusterEta", Photon_superclusterEta, &b_Photon_superclusterEta);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_x_calo", Photon_x_calo, &b_Photon_x_calo);
   fChain->SetBranchAddress("Photon_y_calo", Photon_y_calo, &b_Photon_y_calo);
   fChain->SetBranchAddress("Photon_z_calo", Photon_z_calo, &b_Photon_z_calo);
   fChain->SetBranchAddress("nPPSLocalTrack", &nPPSLocalTrack, &b_nPPSLocalTrack);
   fChain->SetBranchAddress("PPSLocalTrack_multiRPProtonIdx", PPSLocalTrack_multiRPProtonIdx, &b_PPSLocalTrack_multiRPProtonIdx);
   fChain->SetBranchAddress("PPSLocalTrack_singleRPProtonIdx", PPSLocalTrack_singleRPProtonIdx, &b_PPSLocalTrack_singleRPProtonIdx);
   fChain->SetBranchAddress("PPSLocalTrack_decRPId", PPSLocalTrack_decRPId, &b_PPSLocalTrack_decRPId);
   fChain->SetBranchAddress("PPSLocalTrack_rpType", PPSLocalTrack_rpType, &b_PPSLocalTrack_rpType);
   fChain->SetBranchAddress("PPSLocalTrack_x", PPSLocalTrack_x, &b_PPSLocalTrack_x);
   fChain->SetBranchAddress("PPSLocalTrack_y", PPSLocalTrack_y, &b_PPSLocalTrack_y);
   fChain->SetBranchAddress("PPSLocalTrack_time", PPSLocalTrack_time, &b_PPSLocalTrack_time);
   fChain->SetBranchAddress("PPSLocalTrack_timeUnc", PPSLocalTrack_timeUnc, &b_PPSLocalTrack_timeUnc);
   */
   fChain->SetBranchAddress("PuppiMET_covXX", &PuppiMET_covXX, &b_PuppiMET_covXX);
   fChain->SetBranchAddress("PuppiMET_covXY", &PuppiMET_covXY, &b_PuppiMET_covXY);
   fChain->SetBranchAddress("PuppiMET_covYY", &PuppiMET_covYY, &b_PuppiMET_covYY);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredDown", &PuppiMET_phiUnclusteredDown, &b_PuppiMET_phiUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_phiUnclusteredUp", &PuppiMET_phiUnclusteredUp, &b_PuppiMET_phiUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredDown", &PuppiMET_ptUnclusteredDown, &b_PuppiMET_ptUnclusteredDown);
   fChain->SetBranchAddress("PuppiMET_ptUnclusteredUp", &PuppiMET_ptUnclusteredUp, &b_PuppiMET_ptUnclusteredUp);
   fChain->SetBranchAddress("PuppiMET_significance", &PuppiMET_significance, &b_PuppiMET_significance);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("PuppiMET_sumPtUnclustered", &PuppiMET_sumPtUnclustered, &b_PuppiMET_sumPtUnclustered);
   /*
   fChain->SetBranchAddress("nPVBS", &nPVBS, &b_nPVBS);
   fChain->SetBranchAddress("PVBS_chi2", PVBS_chi2, &b_PVBS_chi2);
   fChain->SetBranchAddress("PVBS_cov00", PVBS_cov00, &b_PVBS_cov00);
   fChain->SetBranchAddress("PVBS_cov10", PVBS_cov10, &b_PVBS_cov10);
   fChain->SetBranchAddress("PVBS_cov11", PVBS_cov11, &b_PVBS_cov11);
   fChain->SetBranchAddress("PVBS_cov20", PVBS_cov20, &b_PVBS_cov20);
   fChain->SetBranchAddress("PVBS_cov21", PVBS_cov21, &b_PVBS_cov21);
   fChain->SetBranchAddress("PVBS_cov22", PVBS_cov22, &b_PVBS_cov22);
   fChain->SetBranchAddress("PVBS_x", PVBS_x, &b_PVBS_x);
   fChain->SetBranchAddress("PVBS_y", PVBS_y, &b_PVBS_y);
   fChain->SetBranchAddress("PVBS_z", PVBS_z, &b_PVBS_z);
   fChain->SetBranchAddress("RawPFMET_phi", &RawPFMET_phi, &b_RawPFMET_phi);
   fChain->SetBranchAddress("RawPFMET_pt", &RawPFMET_pt, &b_RawPFMET_pt);
   fChain->SetBranchAddress("RawPFMET_sumEt", &RawPFMET_sumEt, &b_RawPFMET_sumEt);
   fChain->SetBranchAddress("RawPuppiMET_phi", &RawPuppiMET_phi, &b_RawPuppiMET_phi);
   fChain->SetBranchAddress("RawPuppiMET_pt", &RawPuppiMET_pt, &b_RawPuppiMET_pt);
   fChain->SetBranchAddress("RawPuppiMET_sumEt", &RawPuppiMET_sumEt, &b_RawPuppiMET_sumEt);
   fChain->SetBranchAddress("Rho_fixedGridRhoAll", &Rho_fixedGridRhoAll, &b_Rho_fixedGridRhoAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetAll", &Rho_fixedGridRhoFastjetAll, &b_Rho_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentral", &Rho_fixedGridRhoFastjetCentral, &b_Rho_fixedGridRhoFastjetCentral);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralCalo", &Rho_fixedGridRhoFastjetCentralCalo, &b_Rho_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralChargedPileUp", &Rho_fixedGridRhoFastjetCentralChargedPileUp, &b_Rho_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("Rho_fixedGridRhoFastjetCentralNeutral", &Rho_fixedGridRhoFastjetCentralNeutral, &b_Rho_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   fChain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   fChain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   fChain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   fChain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   fChain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   fChain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   fChain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   fChain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   fChain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   fChain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   fChain->SetBranchAddress("nProton_singleRP", &nProton_singleRP, &b_nProton_singleRP);
   fChain->SetBranchAddress("Proton_singleRP_decRPId", Proton_singleRP_decRPId, &b_Proton_singleRP_decRPId);
   fChain->SetBranchAddress("Proton_singleRP_thetaY", Proton_singleRP_thetaY, &b_Proton_singleRP_thetaY);
   fChain->SetBranchAddress("Proton_singleRP_xi", Proton_singleRP_xi, &b_Proton_singleRP_xi);
   fChain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   fChain->SetBranchAddress("SubJet_area", SubJet_area, &b_SubJet_area);
   fChain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   fChain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   fChain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   fChain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   fChain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   fChain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   fChain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   fChain->SetBranchAddress("SubJet_rawFactor", SubJet_rawFactor, &b_SubJet_rawFactor);
   fChain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   fChain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   fChain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   fChain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   fChain->SetBranchAddress("nTauProd", &nTauProd, &b_nTauProd);
   fChain->SetBranchAddress("TauProd_tauIdx", TauProd_tauIdx, &b_TauProd_tauIdx);
   fChain->SetBranchAddress("TauProd_pdgId", TauProd_pdgId, &b_TauProd_pdgId);
   fChain->SetBranchAddress("TauProd_eta", TauProd_eta, &b_TauProd_eta);
   fChain->SetBranchAddress("TauProd_phi", TauProd_phi, &b_TauProd_phi);
   fChain->SetBranchAddress("TauProd_pt", TauProd_pt, &b_TauProd_pt);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_idAntiEleDeadECal", Tau_idAntiEleDeadECal, &b_Tau_idAntiEleDeadECal);
   fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   fChain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("Tau_idDecayModeOldDMs", Tau_idDecayModeOldDMs, &b_Tau_idDecayModeOldDMs);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSe", Tau_idDeepTau2017v2p1VSe, &b_Tau_idDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSjet", Tau_idDeepTau2017v2p1VSjet, &b_Tau_idDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSmu", Tau_idDeepTau2017v2p1VSmu, &b_Tau_idDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_idDeepTau2018v2p5VSe", Tau_idDeepTau2018v2p5VSe, &b_Tau_idDeepTau2018v2p5VSe);
   fChain->SetBranchAddress("Tau_idDeepTau2018v2p5VSjet", Tau_idDeepTau2018v2p5VSjet, &b_Tau_idDeepTau2018v2p5VSjet);
   fChain->SetBranchAddress("Tau_idDeepTau2018v2p5VSmu", Tau_idDeepTau2018v2p5VSmu, &b_Tau_idDeepTau2018v2p5VSmu);
   fChain->SetBranchAddress("Tau_nSVs", Tau_nSVs, &b_Tau_nSVs);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayModePNet", Tau_decayModePNet, &b_Tau_decayModePNet);
   fChain->SetBranchAddress("Tau_decayModeUParT", Tau_decayModeUParT, &b_Tau_decayModeUParT);
   fChain->SetBranchAddress("Tau_eleIdx", Tau_eleIdx, &b_Tau_eleIdx);
   fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   fChain->SetBranchAddress("Tau_muIdx", Tau_muIdx, &b_Tau_muIdx);
   fChain->SetBranchAddress("Tau_svIdx1", Tau_svIdx1, &b_Tau_svIdx1);
   fChain->SetBranchAddress("Tau_svIdx2", Tau_svIdx2, &b_Tau_svIdx2);
   fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("Tau_probDM0PNet", Tau_probDM0PNet, &b_Tau_probDM0PNet);
   fChain->SetBranchAddress("Tau_probDM0UParT", Tau_probDM0UParT, &b_Tau_probDM0UParT);
   fChain->SetBranchAddress("Tau_probDM10PNet", Tau_probDM10PNet, &b_Tau_probDM10PNet);
   fChain->SetBranchAddress("Tau_probDM10UParT", Tau_probDM10UParT, &b_Tau_probDM10UParT);
   fChain->SetBranchAddress("Tau_probDM11PNet", Tau_probDM11PNet, &b_Tau_probDM11PNet);
   fChain->SetBranchAddress("Tau_probDM11UParT", Tau_probDM11UParT, &b_Tau_probDM11UParT);
   fChain->SetBranchAddress("Tau_probDM1PNet", Tau_probDM1PNet, &b_Tau_probDM1PNet);
   fChain->SetBranchAddress("Tau_probDM1UParT", Tau_probDM1UParT, &b_Tau_probDM1UParT);
   fChain->SetBranchAddress("Tau_probDM2PNet", Tau_probDM2PNet, &b_Tau_probDM2PNet);
   fChain->SetBranchAddress("Tau_probDM2UParT", Tau_probDM2UParT, &b_Tau_probDM2UParT);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_ptCorrPNet", Tau_ptCorrPNet, &b_Tau_ptCorrPNet);
   fChain->SetBranchAddress("Tau_ptCorrUParT", Tau_ptCorrUParT, &b_Tau_ptCorrUParT);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("Tau_qConfPNet", Tau_qConfPNet, &b_Tau_qConfPNet);
   fChain->SetBranchAddress("Tau_qConfUParT", Tau_qConfUParT, &b_Tau_qConfUParT);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSe", Tau_rawDeepTau2017v2p1VSe, &b_Tau_rawDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSjet", Tau_rawDeepTau2017v2p1VSjet, &b_Tau_rawDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_rawDeepTau2017v2p1VSmu", Tau_rawDeepTau2017v2p1VSmu, &b_Tau_rawDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_rawDeepTau2018v2p5VSe", Tau_rawDeepTau2018v2p5VSe, &b_Tau_rawDeepTau2018v2p5VSe);
   fChain->SetBranchAddress("Tau_rawDeepTau2018v2p5VSjet", Tau_rawDeepTau2018v2p5VSjet, &b_Tau_rawDeepTau2018v2p5VSjet);
   fChain->SetBranchAddress("Tau_rawDeepTau2018v2p5VSmu", Tau_rawDeepTau2018v2p5VSmu, &b_Tau_rawDeepTau2018v2p5VSmu);
   fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   fChain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   fChain->SetBranchAddress("Tau_rawPNetVSe", Tau_rawPNetVSe, &b_Tau_rawPNetVSe);
   fChain->SetBranchAddress("Tau_rawPNetVSjet", Tau_rawPNetVSjet, &b_Tau_rawPNetVSjet);
   fChain->SetBranchAddress("Tau_rawPNetVSmu", Tau_rawPNetVSmu, &b_Tau_rawPNetVSmu);
   fChain->SetBranchAddress("Tau_rawUParTVSe", Tau_rawUParTVSe, &b_Tau_rawUParTVSe);
   fChain->SetBranchAddress("Tau_rawUParTVSjet", Tau_rawUParTVSjet, &b_Tau_rawUParTVSjet);
   fChain->SetBranchAddress("Tau_rawUParTVSmu", Tau_rawUParTVSmu, &b_Tau_rawUParTVSmu);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrkMET_phi", &TrkMET_phi, &b_TrkMET_phi);
   fChain->SetBranchAddress("TrkMET_pt", &TrkMET_pt, &b_TrkMET_pt);
   fChain->SetBranchAddress("TrkMET_sumEt", &TrkMET_sumEt, &b_TrkMET_sumEt);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("OtherPV_score", OtherPV_score, &b_OtherPV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_sumpt2", &PV_sumpt2, &b_PV_sumpt2);
   fChain->SetBranchAddress("PV_sumpx", &PV_sumpx, &b_PV_sumpx);
   fChain->SetBranchAddress("PV_sumpy", &PV_sumpy, &b_PV_sumpy);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_charge", SV_charge, &b_SV_charge);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_dxy", SV_dxy, &b_SV_dxy);
   fChain->SetBranchAddress("SV_dxySig", SV_dxySig, &b_SV_dxySig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("Electron_IPx", Electron_IPx, &b_Electron_IPx);
   fChain->SetBranchAddress("Electron_IPy", Electron_IPy, &b_Electron_IPy);
   fChain->SetBranchAddress("Electron_IPz", Electron_IPz, &b_Electron_IPz);
   fChain->SetBranchAddress("Electron_ipLengthSig", Electron_ipLengthSig, &b_Electron_ipLengthSig);
   fChain->SetBranchAddress("Muon_IPx", Muon_IPx, &b_Muon_IPx);
   fChain->SetBranchAddress("Muon_IPy", Muon_IPy, &b_Muon_IPy);
   fChain->SetBranchAddress("Muon_IPz", Muon_IPz, &b_Muon_IPz);
   fChain->SetBranchAddress("Muon_ipLengthSig", Muon_ipLengthSig, &b_Muon_ipLengthSig);
   */
   fChain->SetBranchAddress("PFCand_fromPVvertexRef", PFCand_fromPVvertexRef, &b_PFCand_fromPVvertexRef);
   fChain->SetBranchAddress("PFCand_blockIdx", PFCand_blockIdx, &b_PFCand_blockIdx);
   fChain->SetBranchAddress("PFCand_nPFTrack", PFCand_nPFTrack, &b_PFCand_nPFTrack);
   fChain->SetBranchAddress("PFCand_nPFClusterHCAL", PFCand_nPFClusterHCAL, &b_PFCand_nPFClusterHCAL);
   fChain->SetBranchAddress("PFCand_nPFClusterECAL", PFCand_nPFClusterECAL, &b_PFCand_nPFClusterECAL);
   fChain->SetBranchAddress("PFCand_nPFClusterPS", PFCand_nPFClusterPS, &b_PFCand_nPFClusterPS);
   fChain->SetBranchAddress("PFCand_nPFTrackInBlock", PFCand_nPFTrackInBlock, &b_PFCand_nPFTrackInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterHCALInBlock", PFCand_nPFClusterHCALInBlock, &b_PFCand_nPFClusterHCALInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterECALInBlock", PFCand_nPFClusterECALInBlock, &b_PFCand_nPFClusterECALInBlock);
   fChain->SetBranchAddress("PFCand_nPFClusterPSInBlock", PFCand_nPFClusterPSInBlock, &b_PFCand_nPFClusterPSInBlock);
   fChain->SetBranchAddress("PFCand_hcalDepthEnergyFraction1", PFCand_hcalDepthEnergyFraction1, &b_PFCand_hcalDepthEnergyFraction1);
   fChain->SetBranchAddress("PFCand_hcalDepthEnergyFraction2", PFCand_hcalDepthEnergyFraction2, &b_PFCand_hcalDepthEnergyFraction2);
   fChain->SetBranchAddress("PFCand_hcalDepthEnergyFraction3", PFCand_hcalDepthEnergyFraction3, &b_PFCand_hcalDepthEnergyFraction3);
   fChain->SetBranchAddress("PFCand_hcalDepthEnergyFraction4", PFCand_hcalDepthEnergyFraction4, &b_PFCand_hcalDepthEnergyFraction4);
   fChain->SetBranchAddress("PFCand_hcalDepthEnergyFraction5", PFCand_hcalDepthEnergyFraction5, &b_PFCand_hcalDepthEnergyFraction5);
   fChain->SetBranchAddress("PFCand_hcalDepthEnergyFraction6", PFCand_hcalDepthEnergyFraction6, &b_PFCand_hcalDepthEnergyFraction6);
   fChain->SetBranchAddress("PFCand_hcalDepthEnergyFraction7", PFCand_hcalDepthEnergyFraction7, &b_PFCand_hcalDepthEnergyFraction7);
   fChain->SetBranchAddress("PFCand_ecalEnergy", PFCand_ecalEnergy, &b_PFCand_ecalEnergy);
   fChain->SetBranchAddress("PFCand_rawEcalEnergy", PFCand_rawEcalEnergy, &b_PFCand_rawEcalEnergy);
   fChain->SetBranchAddress("PFCand_hcalEnergy", PFCand_hcalEnergy, &b_PFCand_hcalEnergy);
   fChain->SetBranchAddress("PFCand_rawHcalEnergy", PFCand_rawHcalEnergy, &b_PFCand_rawHcalEnergy);
   fChain->SetBranchAddress("PFCand_hoEnergy", PFCand_hoEnergy, &b_PFCand_hoEnergy);
   fChain->SetBranchAddress("PFCand_rawHoEnergy", PFCand_rawHoEnergy, &b_PFCand_rawHoEnergy);
   /*
   fChain->SetBranchAddress("SV_ntracks", SV_ntracks, &b_SV_ntracks);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("Tau_hasRefitSV", Tau_hasRefitSV, &b_Tau_hasRefitSV);
   fChain->SetBranchAddress("Tau_IPx", Tau_IPx, &b_Tau_IPx);
   fChain->SetBranchAddress("Tau_IPy", Tau_IPy, &b_Tau_IPy);
   fChain->SetBranchAddress("Tau_IPz", Tau_IPz, &b_Tau_IPz);
   fChain->SetBranchAddress("Tau_ipLengthSig", Tau_ipLengthSig, &b_Tau_ipLengthSig);
   fChain->SetBranchAddress("Tau_refitSVchi2", Tau_refitSVchi2, &b_Tau_refitSVchi2);
   fChain->SetBranchAddress("Tau_refitSVcov00", Tau_refitSVcov00, &b_Tau_refitSVcov00);
   fChain->SetBranchAddress("Tau_refitSVcov10", Tau_refitSVcov10, &b_Tau_refitSVcov10);
   fChain->SetBranchAddress("Tau_refitSVcov11", Tau_refitSVcov11, &b_Tau_refitSVcov11);
   fChain->SetBranchAddress("Tau_refitSVcov20", Tau_refitSVcov20, &b_Tau_refitSVcov20);
   fChain->SetBranchAddress("Tau_refitSVcov21", Tau_refitSVcov21, &b_Tau_refitSVcov21);
   fChain->SetBranchAddress("Tau_refitSVcov22", Tau_refitSVcov22, &b_Tau_refitSVcov22);
   fChain->SetBranchAddress("Tau_refitSVx", Tau_refitSVx, &b_Tau_refitSVx);
   fChain->SetBranchAddress("Tau_refitSVy", Tau_refitSVy, &b_Tau_refitSVy);
   fChain->SetBranchAddress("Tau_refitSVz", Tau_refitSVz, &b_Tau_refitSVz);
   fChain->SetBranchAddress("L1Reco_step", &L1Reco_step, &b_L1Reco_step);
   */
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter, &b_Flag_hfNoisyHitsFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   /*
   fChain->SetBranchAddress("L1_AXO_Loose", &L1_AXO_Loose, &b_L1_AXO_Loose);
   fChain->SetBranchAddress("L1_AXO_Nominal", &L1_AXO_Nominal, &b_L1_AXO_Nominal);
   fChain->SetBranchAddress("L1_AXO_Tight", &L1_AXO_Tight, &b_L1_AXO_Tight);
   fChain->SetBranchAddress("L1_AXO_VLoose", &L1_AXO_VLoose, &b_L1_AXO_VLoose);
   fChain->SetBranchAddress("L1_AXO_VTight", &L1_AXO_VTight, &b_L1_AXO_VTight);
   fChain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME, &b_L1_BPTX_AND_Ref1_VME);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME, &b_L1_BPTX_AND_Ref3_VME);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME, &b_L1_BPTX_AND_Ref4_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME, &b_L1_BPTX_BeamGas_B1_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME, &b_L1_BPTX_BeamGas_B2_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME, &b_L1_BPTX_BeamGas_Ref1_VME);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME, &b_L1_BPTX_BeamGas_Ref2_VME);
   fChain->SetBranchAddress("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME, &b_L1_BPTX_NotOR_VME);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME, &b_L1_BPTX_OR_Ref3_VME);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME, &b_L1_BPTX_OR_Ref4_VME);
   fChain->SetBranchAddress("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME, &b_L1_BPTX_RefAND_VME);
   fChain->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus, &b_L1_BptxMinus);
   fChain->SetBranchAddress("L1_BptxOR", &L1_BptxOR, &b_L1_BptxOR);
   fChain->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus, &b_L1_BptxPlus);
   fChain->SetBranchAddress("L1_BptxXOR", &L1_BptxXOR, &b_L1_BptxXOR);
   fChain->SetBranchAddress("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("L1_CICADA_Loose", &L1_CICADA_Loose, &b_L1_CICADA_Loose);
   fChain->SetBranchAddress("L1_CICADA_Medium", &L1_CICADA_Medium, &b_L1_CICADA_Medium);
   fChain->SetBranchAddress("L1_CICADA_Tight", &L1_CICADA_Tight, &b_L1_CICADA_Tight);
   fChain->SetBranchAddress("L1_CICADA_VLoose", &L1_CICADA_VLoose, &b_L1_CICADA_VLoose);
   fChain->SetBranchAddress("L1_CICADA_VTight", &L1_CICADA_VTight, &b_L1_CICADA_VTight);
   fChain->SetBranchAddress("L1_DoubleEG11_er1p2_dR_Max0p6", &L1_DoubleEG11_er1p2_dR_Max0p6, &b_L1_DoubleEG11_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er, &b_L1_DoubleEG8er2p5_HTT280er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er, &b_L1_DoubleEG8er2p5_HTT300er);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er, &b_L1_DoubleEG8er2p5_HTT320er);
   fChain->SetBranchAddress("L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5, &b_L1_DoubleEG_15_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5, &b_L1_DoubleEG_20_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5, &b_L1_DoubleEG_22_10_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5, &b_L1_DoubleEG_25_12_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5, &b_L1_DoubleEG_25_14_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5, &b_L1_DoubleEG_27_14_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso16_LooseIso12_er1p5", &L1_DoubleEG_LooseIso16_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso16_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso18_LooseIso12_er1p5", &L1_DoubleEG_LooseIso18_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso18_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso20_LooseIso12_er1p5", &L1_DoubleEG_LooseIso20_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso20_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5, &b_L1_DoubleEG_LooseIso22_12_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso22_LooseIso12_er1p5", &L1_DoubleEG_LooseIso22_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso22_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5, &b_L1_DoubleEG_LooseIso25_12_er2p5);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso25_LooseIso12_er1p5", &L1_DoubleEG_LooseIso25_LooseIso12_er1p5, &b_L1_DoubleEG_LooseIso25_LooseIso12_er1p5);
   fChain->SetBranchAddress("L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5, &b_L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5);
   fChain->SetBranchAddress("L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5, &b_L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5);
   fChain->SetBranchAddress("L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1, &b_L1_DoubleIsoTau28er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1, &b_L1_DoubleIsoTau30er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1, &b_L1_DoubleIsoTau32er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau32er2p1_Mass_Max80", &L1_DoubleIsoTau32er2p1_Mass_Max80, &b_L1_DoubleIsoTau32er2p1_Mass_Max80);
   fChain->SetBranchAddress("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1, &b_L1_DoubleIsoTau34er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau35er2p1", &L1_DoubleIsoTau35er2p1, &b_L1_DoubleIsoTau35er2p1);
   fChain->SetBranchAddress("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1, &b_L1_DoubleIsoTau36er2p1);
   fChain->SetBranchAddress("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6, &b_L1_DoubleJet100er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5, &b_L1_DoubleJet100er2p5);
   fChain->SetBranchAddress("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6, &b_L1_DoubleJet112er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5, &b_L1_DoubleJet120er2p5);
   fChain->SetBranchAddress("L1_DoubleJet120er2p5_Mu3_dR_Max0p8", &L1_DoubleJet120er2p5_Mu3_dR_Max0p8, &b_L1_DoubleJet120er2p5_Mu3_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5, &b_L1_DoubleJet150er2p5);
   fChain->SetBranchAddress("L1_DoubleJet16er2p5_Mu3_dR_Max0p4", &L1_DoubleJet16er2p5_Mu3_dR_Max0p4, &b_L1_DoubleJet16er2p5_Mu3_dR_Max0p4);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5, &b_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleJet35er2p5_Mu3_dR_Max0p4", &L1_DoubleJet35er2p5_Mu3_dR_Max0p4, &b_L1_DoubleJet35er2p5_Mu3_dR_Max0p4);
   fChain->SetBranchAddress("L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5, &b_L1_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5", &L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5, &b_L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2", &L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2, &b_L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5", &L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5, &b_L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2", &L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2, &b_L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2);
   fChain->SetBranchAddress("L1_DoubleJet60er2p5_Mu3_dR_Max0p4", &L1_DoubleJet60er2p5_Mu3_dR_Max0p4, &b_L1_DoubleJet60er2p5_Mu3_dR_Max0p4);
   fChain->SetBranchAddress("L1_DoubleJet80er2p5_Mu3_dR_Max0p4", &L1_DoubleJet80er2p5_Mu3_dR_Max0p4, &b_L1_DoubleJet80er2p5_Mu3_dR_Max0p4);
   fChain->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min800", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min800, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min800);
   fChain->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min850", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min850, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min850);
   fChain->SetBranchAddress("L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50", &L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50, &b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50);
   fChain->SetBranchAddress("L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50", &L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50, &b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50);
   fChain->SetBranchAddress("L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65", &L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65, &b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65);
   fChain->SetBranchAddress("L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65", &L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65, &b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65);
   fChain->SetBranchAddress("L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ", &L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ, &b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ);
   fChain->SetBranchAddress("L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ", &L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ, &b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ);
   fChain->SetBranchAddress("L1_DoubleLLPJet40", &L1_DoubleLLPJet40, &b_L1_DoubleLLPJet40);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1, &b_L1_DoubleLooseIsoEG22er2p1);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1, &b_L1_DoubleLooseIsoEG24er2p1);
   fChain->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0, &b_L1_DoubleMu0);
   fChain->SetBranchAddress("L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1, &b_L1_DoubleMu0_Mass_Min1);
   fChain->SetBranchAddress("L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ, &b_L1_DoubleMu0_OQ);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ, &b_L1_DoubleMu0_SQ);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS, &b_L1_DoubleMu0_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt15_Upt7", &L1_DoubleMu0_Upt15_Upt7, &b_L1_DoubleMu0_Upt15_Upt7);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt5_Upt5", &L1_DoubleMu0_Upt5_Upt5, &b_L1_DoubleMu0_Upt5_Upt5);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt6_IP_Min1_Upt4", &L1_DoubleMu0_Upt6_IP_Min1_Upt4, &b_L1_DoubleMu0_Upt6_IP_Min1_Upt4);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt6_SQ_er2p0", &L1_DoubleMu0_Upt6_SQ_er2p0, &b_L1_DoubleMu0_Upt6_SQ_er2p0);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt7_SQ_er2p0", &L1_DoubleMu0_Upt7_SQ_er2p0, &b_L1_DoubleMu0_Upt7_SQ_er2p0);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt8_SQ_er2p0", &L1_DoubleMu0_Upt8_SQ_er2p0, &b_L1_DoubleMu0_Upt8_SQ_er2p0);
   fChain->SetBranchAddress("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8, &b_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6, &b_L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2", &L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2, &b_L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ, &b_L1_DoubleMu0er1p5_SQ);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS, &b_L1_DoubleMu0er1p5_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2", &L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2, &b_L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_dR_Max1p4);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_dEta_Max1p5, &b_L1_DoubleMu0er2p0_SQ_dEta_Max1p5);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_dEta_Max1p6, &b_L1_DoubleMu0er2p0_SQ_dEta_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu18er2p1_SQ", &L1_DoubleMu18er2p1_SQ, &b_L1_DoubleMu18er2p1_SQ);
   fChain->SetBranchAddress("L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20", &L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20, &b_L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF30_HTT60er", &L1_DoubleMu3_SQ_ETMHF30_HTT60er, &b_L1_DoubleMu3_SQ_ETMHF30_HTT60er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5, &b_L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_HTT60er", &L1_DoubleMu3_SQ_ETMHF40_HTT60er, &b_L1_DoubleMu3_SQ_ETMHF40_HTT60er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5, &b_L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er, &b_L1_DoubleMu3_SQ_ETMHF50_HTT60er);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5, &b_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er, &b_L1_DoubleMu3_SQ_HTT220er);
   fChain->SetBranchAddress("L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6", &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6, &b_L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5, &b_L1_DoubleMu4_SQ_EG9er2p5);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS, &b_L1_DoubleMu4_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6", &L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6, &b_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS, &b_L1_DoubleMu4p5_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS, &b_L1_DoubleMu4p5er2p0_SQ_OS);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7);
   fChain->SetBranchAddress("L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20", &L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20, &b_L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20);
   fChain->SetBranchAddress("L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5, &b_L1_DoubleMu5_SQ_EG9er2p5);
   fChain->SetBranchAddress("L1_DoubleMu5_SQ_OS_dR_Max1p6", &L1_DoubleMu5_SQ_OS_dR_Max1p6, &b_L1_DoubleMu5_SQ_OS_dR_Max1p6);
   fChain->SetBranchAddress("L1_DoubleMu6_Upt6_SQ_er2p0", &L1_DoubleMu6_Upt6_SQ_er2p0, &b_L1_DoubleMu6_Upt6_SQ_er2p0);
   fChain->SetBranchAddress("L1_DoubleMu7_Upt7_SQ_er2p0", &L1_DoubleMu7_Upt7_SQ_er2p0, &b_L1_DoubleMu7_Upt7_SQ_er2p0);
   fChain->SetBranchAddress("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ, &b_L1_DoubleMu8_SQ);
   fChain->SetBranchAddress("L1_DoubleMu8_Upt8_SQ_er2p0", &L1_DoubleMu8_Upt8_SQ_er2p0, &b_L1_DoubleMu8_Upt8_SQ_er2p0);
   fChain->SetBranchAddress("L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ, &b_L1_DoubleMu9_SQ);
   fChain->SetBranchAddress("L1_DoubleMu_12_5", &L1_DoubleMu_12_5, &b_L1_DoubleMu_12_5);
   fChain->SetBranchAddress("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ, &b_L1_DoubleMu_15_5_SQ);
   fChain->SetBranchAddress("L1_DoubleMu_15_7", &L1_DoubleMu_15_7, &b_L1_DoubleMu_15_7);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1, &b_L1_DoubleMu_15_7_Mass_Min1);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ, &b_L1_DoubleMu_15_7_SQ);
   fChain->SetBranchAddress("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1, &b_L1_DoubleTau70er2p1);
   fChain->SetBranchAddress("L1_ETM120", &L1_ETM120, &b_L1_ETM120);
   fChain->SetBranchAddress("L1_ETM150", &L1_ETM150, &b_L1_ETM150);
   fChain->SetBranchAddress("L1_ETMHF100", &L1_ETMHF100, &b_L1_ETMHF100);
   fChain->SetBranchAddress("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er, &b_L1_ETMHF100_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF110", &L1_ETMHF110, &b_L1_ETMHF110);
   fChain->SetBranchAddress("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er, &b_L1_ETMHF110_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF120", &L1_ETMHF120, &b_L1_ETMHF120);
   fChain->SetBranchAddress("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er, &b_L1_ETMHF120_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF130", &L1_ETMHF130, &b_L1_ETMHF130);
   fChain->SetBranchAddress("L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er, &b_L1_ETMHF130_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF140", &L1_ETMHF140, &b_L1_ETMHF140);
   fChain->SetBranchAddress("L1_ETMHF150", &L1_ETMHF150, &b_L1_ETMHF150);
   fChain->SetBranchAddress("L1_ETMHF70", &L1_ETMHF70, &b_L1_ETMHF70);
   fChain->SetBranchAddress("L1_ETMHF70_HTT60er", &L1_ETMHF70_HTT60er, &b_L1_ETMHF70_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF80", &L1_ETMHF80, &b_L1_ETMHF80);
   fChain->SetBranchAddress("L1_ETMHF80_HTT60er", &L1_ETMHF80_HTT60er, &b_L1_ETMHF80_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1, &b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1);
   fChain->SetBranchAddress("L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6, &b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6);
   fChain->SetBranchAddress("L1_ETMHF90", &L1_ETMHF90, &b_L1_ETMHF90);
   fChain->SetBranchAddress("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er, &b_L1_ETMHF90_HTT60er);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1, &b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6, &b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6);
   fChain->SetBranchAddress("L1_ETT2000", &L1_ETT2000, &b_L1_ETT2000);
   fChain->SetBranchAddress("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain, &b_L1_FirstBunchAfterTrain);
   fChain->SetBranchAddress("L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain, &b_L1_FirstBunchBeforeTrain);
   fChain->SetBranchAddress("L1_FirstBunchInTrain", &L1_FirstBunchInTrain, &b_L1_FirstBunchInTrain);
   fChain->SetBranchAddress("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit, &b_L1_FirstCollisionInOrbit);
   fChain->SetBranchAddress("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain, &b_L1_FirstCollisionInTrain);
   fChain->SetBranchAddress("L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig, &b_L1_HCAL_LaserMon_Trig);
   fChain->SetBranchAddress("L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto, &b_L1_HCAL_LaserMon_Veto);
   fChain->SetBranchAddress("L1_HTMHF100", &L1_HTMHF100, &b_L1_HTMHF100);
   fChain->SetBranchAddress("L1_HTMHF120", &L1_HTMHF120, &b_L1_HTMHF120);
   fChain->SetBranchAddress("L1_HTMHF125", &L1_HTMHF125, &b_L1_HTMHF125);
   fChain->SetBranchAddress("L1_HTMHF130", &L1_HTMHF130, &b_L1_HTMHF130);
   fChain->SetBranchAddress("L1_HTMHF150", &L1_HTMHF150, &b_L1_HTMHF150);
   fChain->SetBranchAddress("L1_HTT120_SingleLLPJet40", &L1_HTT120_SingleLLPJet40, &b_L1_HTT120_SingleLLPJet40);
   fChain->SetBranchAddress("L1_HTT120er", &L1_HTT120er, &b_L1_HTT120er);
   fChain->SetBranchAddress("L1_HTT160_SingleLLPJet50", &L1_HTT160_SingleLLPJet50, &b_L1_HTT160_SingleLLPJet50);
   fChain->SetBranchAddress("L1_HTT160er", &L1_HTT160er, &b_L1_HTT160er);
   fChain->SetBranchAddress("L1_HTT200_SingleLLPJet60", &L1_HTT200_SingleLLPJet60, &b_L1_HTT200_SingleLLPJet60);
   fChain->SetBranchAddress("L1_HTT200er", &L1_HTT200er, &b_L1_HTT200er);
   fChain->SetBranchAddress("L1_HTT240_SingleLLPJet70", &L1_HTT240_SingleLLPJet70, &b_L1_HTT240_SingleLLPJet70);
   fChain->SetBranchAddress("L1_HTT255er", &L1_HTT255er, &b_L1_HTT255er);
   fChain->SetBranchAddress("L1_HTT280er", &L1_HTT280er, &b_L1_HTT280er);
   fChain->SetBranchAddress("L1_HTT280er_QuadJet_70_55_40_35_er2p5", &L1_HTT280er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT280er_QuadJet_70_55_40_35_er2p5);
   fChain->SetBranchAddress("L1_HTT320er", &L1_HTT320er, &b_L1_HTT320er);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p5", &L1_HTT320er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p5);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3, &b_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3, &b_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3);
   fChain->SetBranchAddress("L1_HTT360er", &L1_HTT360er, &b_L1_HTT360er);
   fChain->SetBranchAddress("L1_HTT400er", &L1_HTT400er, &b_L1_HTT400er);
   fChain->SetBranchAddress("L1_HTT450er", &L1_HTT450er, &b_L1_HTT450er);
   fChain->SetBranchAddress("L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40, &b_L1_IsoEG32er2p5_Mt40);
   fChain->SetBranchAddress("L1_IsoTau52er2p1_QuadJet36er2p5", &L1_IsoTau52er2p1_QuadJet36er2p5, &b_L1_IsoTau52er2p1_QuadJet36er2p5);
   fChain->SetBranchAddress("L1_IsolatedBunch", &L1_IsolatedBunch, &b_L1_IsolatedBunch);
   fChain->SetBranchAddress("L1_LastBunchInTrain", &L1_LastBunchInTrain, &b_L1_LastBunchInTrain);
   fChain->SetBranchAddress("L1_LastCollisionInTrain", &L1_LastCollisionInTrain, &b_L1_LastCollisionInTrain);
   fChain->SetBranchAddress("L1_LooseIsoEG14er2p5_HTT200er", &L1_LooseIsoEG14er2p5_HTT200er, &b_L1_LooseIsoEG14er2p5_HTT200er);
   fChain->SetBranchAddress("L1_LooseIsoEG16er2p5_HTT200er", &L1_LooseIsoEG16er2p5_HTT200er, &b_L1_LooseIsoEG16er2p5_HTT200er);
   fChain->SetBranchAddress("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3, &b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3, &b_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er, &b_L1_LooseIsoEG24er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er, &b_L1_LooseIsoEG26er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3, &b_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er, &b_L1_LooseIsoEG28er2p1_HTT100er);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3, &b_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3, &b_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3);
   fChain->SetBranchAddress("L1_MinimumBiasHF0", &L1_MinimumBiasHF0, &b_L1_MinimumBiasHF0);
   fChain->SetBranchAddress("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND, &b_L1_MinimumBiasHF0_AND_BptxAND);
   fChain->SetBranchAddress("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6, &b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_Mu12_HTT150er", &L1_Mu12_HTT150er, &b_L1_Mu12_HTT150er);
   fChain->SetBranchAddress("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6, &b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6);
   fChain->SetBranchAddress("L1_Mu14_HTT150er", &L1_Mu14_HTT150er, &b_L1_Mu14_HTT150er);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1, &b_L1_Mu18er2p1_Tau24er2p1);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1, &b_L1_Mu18er2p1_Tau26er2p1);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1_Jet55", &L1_Mu18er2p1_Tau26er2p1_Jet55, &b_L1_Mu18er2p1_Tau26er2p1_Jet55);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1_Jet70", &L1_Mu18er2p1_Tau26er2p1_Jet70, &b_L1_Mu18er2p1_Tau26er2p1_Jet70);
   fChain->SetBranchAddress("L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5, &b_L1_Mu20_EG10er2p5);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1, &b_L1_Mu22er2p1_IsoTau30er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1, &b_L1_Mu22er2p1_IsoTau32er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1, &b_L1_Mu22er2p1_IsoTau34er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1, &b_L1_Mu22er2p1_IsoTau40er2p1);
   fChain->SetBranchAddress("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1, &b_L1_Mu22er2p1_Tau70er2p1);
   fChain->SetBranchAddress("L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4, &b_L1_Mu3_Jet120er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4, &b_L1_Mu3_Jet16er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5, &b_L1_Mu3_Jet30er2p5);
   fChain->SetBranchAddress("L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4, &b_L1_Mu3_Jet60er2p5_dR_Max0p4);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF30", &L1_Mu3er1p5_Jet100er2p5_ETMHF30, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF30);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF40);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF50);
   fChain->SetBranchAddress("L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5, &b_L1_Mu5_EG23er2p5);
   fChain->SetBranchAddress("L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5, &b_L1_Mu5_LooseIsoEG20er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5, &b_L1_Mu6_DoubleEG10er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5, &b_L1_Mu6_DoubleEG12er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5, &b_L1_Mu6_DoubleEG15er2p5);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5, &b_L1_Mu6_DoubleEG17er2p5);
   fChain->SetBranchAddress("L1_Mu6_HTT240er", &L1_Mu6_HTT240er, &b_L1_Mu6_HTT240er);
   fChain->SetBranchAddress("L1_Mu6_HTT250er", &L1_Mu6_HTT250er, &b_L1_Mu6_HTT250er);
   fChain->SetBranchAddress("L1_Mu7_EG20er2p5", &L1_Mu7_EG20er2p5, &b_L1_Mu7_EG20er2p5);
   fChain->SetBranchAddress("L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5, &b_L1_Mu7_EG23er2p5);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5, &b_L1_Mu7_LooseIsoEG20er2p5);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5, &b_L1_Mu7_LooseIsoEG23er2p5);
   fChain->SetBranchAddress("L1_NotBptxOR", &L1_NotBptxOR, &b_L1_NotBptxOR);
   fChain->SetBranchAddress("L1_QuadJet60er2p5", &L1_QuadJet60er2p5, &b_L1_QuadJet60er2p5);
   fChain->SetBranchAddress("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0, &b_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0);
   fChain->SetBranchAddress("L1_QuadMu0", &L1_QuadMu0, &b_L1_QuadMu0);
   fChain->SetBranchAddress("L1_QuadMu0_OQ", &L1_QuadMu0_OQ, &b_L1_QuadMu0_OQ);
   fChain->SetBranchAddress("L1_QuadMu0_SQ", &L1_QuadMu0_SQ, &b_L1_QuadMu0_SQ);
   fChain->SetBranchAddress("L1_SecondBunchInTrain", &L1_SecondBunchInTrain, &b_L1_SecondBunchInTrain);
   fChain->SetBranchAddress("L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain, &b_L1_SecondLastBunchInTrain);
   fChain->SetBranchAddress("L1_SingleEG10er2p5", &L1_SingleEG10er2p5, &b_L1_SingleEG10er2p5);
   fChain->SetBranchAddress("L1_SingleEG15er2p5", &L1_SingleEG15er2p5, &b_L1_SingleEG15er2p5);
   fChain->SetBranchAddress("L1_SingleEG26er2p5", &L1_SingleEG26er2p5, &b_L1_SingleEG26er2p5);
   fChain->SetBranchAddress("L1_SingleEG28_FWD2p5", &L1_SingleEG28_FWD2p5, &b_L1_SingleEG28_FWD2p5);
   fChain->SetBranchAddress("L1_SingleEG28er1p5", &L1_SingleEG28er1p5, &b_L1_SingleEG28er1p5);
   fChain->SetBranchAddress("L1_SingleEG28er2p1", &L1_SingleEG28er2p1, &b_L1_SingleEG28er2p1);
   fChain->SetBranchAddress("L1_SingleEG28er2p5", &L1_SingleEG28er2p5, &b_L1_SingleEG28er2p5);
   fChain->SetBranchAddress("L1_SingleEG34er2p5", &L1_SingleEG34er2p5, &b_L1_SingleEG34er2p5);
   fChain->SetBranchAddress("L1_SingleEG36er2p5", &L1_SingleEG36er2p5, &b_L1_SingleEG36er2p5);
   fChain->SetBranchAddress("L1_SingleEG38er2p5", &L1_SingleEG38er2p5, &b_L1_SingleEG38er2p5);
   fChain->SetBranchAddress("L1_SingleEG40er2p5", &L1_SingleEG40er2p5, &b_L1_SingleEG40er2p5);
   fChain->SetBranchAddress("L1_SingleEG42er2p5", &L1_SingleEG42er2p5, &b_L1_SingleEG42er2p5);
   fChain->SetBranchAddress("L1_SingleEG45er2p5", &L1_SingleEG45er2p5, &b_L1_SingleEG45er2p5);
   fChain->SetBranchAddress("L1_SingleEG50", &L1_SingleEG50, &b_L1_SingleEG50);
   fChain->SetBranchAddress("L1_SingleEG60", &L1_SingleEG60, &b_L1_SingleEG60);
   fChain->SetBranchAddress("L1_SingleEG8er2p5", &L1_SingleEG8er2p5, &b_L1_SingleEG8er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1, &b_L1_SingleIsoEG24er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1, &b_L1_SingleIsoEG26er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5, &b_L1_SingleIsoEG26er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG28_FWD2p5", &L1_SingleIsoEG28_FWD2p5, &b_L1_SingleIsoEG28_FWD2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5, &b_L1_SingleIsoEG28er1p5);
   fChain->SetBranchAddress("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1, &b_L1_SingleIsoEG28er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5, &b_L1_SingleIsoEG28er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1, &b_L1_SingleIsoEG30er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5, &b_L1_SingleIsoEG30er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1, &b_L1_SingleIsoEG32er2p1);
   fChain->SetBranchAddress("L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5, &b_L1_SingleIsoEG32er2p5);
   fChain->SetBranchAddress("L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5, &b_L1_SingleIsoEG34er2p5);
   fChain->SetBranchAddress("L1_SingleJet10erHE", &L1_SingleJet10erHE, &b_L1_SingleJet10erHE);
   fChain->SetBranchAddress("L1_SingleJet120", &L1_SingleJet120, &b_L1_SingleJet120);
   fChain->SetBranchAddress("L1_SingleJet120_FWD2p5", &L1_SingleJet120_FWD2p5, &b_L1_SingleJet120_FWD2p5);
   fChain->SetBranchAddress("L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0, &b_L1_SingleJet120_FWD3p0);
   fChain->SetBranchAddress("L1_SingleJet120er1p3", &L1_SingleJet120er1p3, &b_L1_SingleJet120er1p3);
   fChain->SetBranchAddress("L1_SingleJet120er2p5", &L1_SingleJet120er2p5, &b_L1_SingleJet120er2p5);
   fChain->SetBranchAddress("L1_SingleJet12erHE", &L1_SingleJet12erHE, &b_L1_SingleJet12erHE);
   fChain->SetBranchAddress("L1_SingleJet140er2p5", &L1_SingleJet140er2p5, &b_L1_SingleJet140er2p5);
   fChain->SetBranchAddress("L1_SingleJet160er2p5", &L1_SingleJet160er2p5, &b_L1_SingleJet160er2p5);
   fChain->SetBranchAddress("L1_SingleJet180", &L1_SingleJet180, &b_L1_SingleJet180);
   fChain->SetBranchAddress("L1_SingleJet180er2p5", &L1_SingleJet180er2p5, &b_L1_SingleJet180er2p5);
   fChain->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   fChain->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR, &b_L1_SingleJet20er2p5_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX, &b_L1_SingleJet20er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet35", &L1_SingleJet35, &b_L1_SingleJet35);
   fChain->SetBranchAddress("L1_SingleJet35_FWD2p5", &L1_SingleJet35_FWD2p5, &b_L1_SingleJet35_FWD2p5);
   fChain->SetBranchAddress("L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0, &b_L1_SingleJet35_FWD3p0);
   fChain->SetBranchAddress("L1_SingleJet35er1p3", &L1_SingleJet35er1p3, &b_L1_SingleJet35er1p3);
   fChain->SetBranchAddress("L1_SingleJet35er2p5", &L1_SingleJet35er2p5, &b_L1_SingleJet35er2p5);
   fChain->SetBranchAddress("L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX, &b_L1_SingleJet43er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX, &b_L1_SingleJet46er2p5_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleJet60", &L1_SingleJet60, &b_L1_SingleJet60);
   fChain->SetBranchAddress("L1_SingleJet60_FWD2p5", &L1_SingleJet60_FWD2p5, &b_L1_SingleJet60_FWD2p5);
   fChain->SetBranchAddress("L1_SingleJet8erHE", &L1_SingleJet8erHE, &b_L1_SingleJet8erHE);
   fChain->SetBranchAddress("L1_SingleJet90", &L1_SingleJet90, &b_L1_SingleJet90);
   fChain->SetBranchAddress("L1_SingleJet90_FWD2p5", &L1_SingleJet90_FWD2p5, &b_L1_SingleJet90_FWD2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG26er1p5", &L1_SingleLooseIsoEG26er1p5, &b_L1_SingleLooseIsoEG26er1p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG26er2p5", &L1_SingleLooseIsoEG26er2p5, &b_L1_SingleLooseIsoEG26er2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28_FWD2p5", &L1_SingleLooseIsoEG28_FWD2p5, &b_L1_SingleLooseIsoEG28_FWD2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5, &b_L1_SingleLooseIsoEG28er1p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er2p1", &L1_SingleLooseIsoEG28er2p1, &b_L1_SingleLooseIsoEG28er2p1);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er2p5", &L1_SingleLooseIsoEG28er2p5, &b_L1_SingleLooseIsoEG28er2p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5, &b_L1_SingleLooseIsoEG30er1p5);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG30er2p5", &L1_SingleLooseIsoEG30er2p5, &b_L1_SingleLooseIsoEG30er2p5);
   fChain->SetBranchAddress("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF, &b_L1_SingleMu0_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_DQ", &L1_SingleMu0_DQ, &b_L1_SingleMu0_DQ);
   fChain->SetBranchAddress("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF, &b_L1_SingleMu0_EMTF);
   fChain->SetBranchAddress("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF, &b_L1_SingleMu0_OMTF);
   fChain->SetBranchAddress("L1_SingleMu0_SQ13_BMTF", &L1_SingleMu0_SQ13_BMTF, &b_L1_SingleMu0_SQ13_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_SQ14_BMTF", &L1_SingleMu0_SQ14_BMTF, &b_L1_SingleMu0_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_SQ15_BMTF", &L1_SingleMu0_SQ15_BMTF, &b_L1_SingleMu0_SQ15_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10", &L1_SingleMu0_Upt10, &b_L1_SingleMu0_Upt10);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_BMTF", &L1_SingleMu0_Upt10_BMTF, &b_L1_SingleMu0_Upt10_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_EMTF", &L1_SingleMu0_Upt10_EMTF, &b_L1_SingleMu0_Upt10_EMTF);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_OMTF", &L1_SingleMu0_Upt10_OMTF, &b_L1_SingleMu0_Upt10_OMTF);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_SQ14_BMTF", &L1_SingleMu0_Upt10_SQ14_BMTF, &b_L1_SingleMu0_Upt10_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_Upt15_SQ14_BMTF", &L1_SingleMu0_Upt15_SQ14_BMTF, &b_L1_SingleMu0_Upt15_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_Upt20_SQ14_BMTF", &L1_SingleMu0_Upt20_SQ14_BMTF, &b_L1_SingleMu0_Upt20_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu0_Upt25_SQ14_BMTF", &L1_SingleMu0_Upt25_SQ14_BMTF, &b_L1_SingleMu0_Upt25_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu10_SQ14_BMTF", &L1_SingleMu10_SQ14_BMTF, &b_L1_SingleMu10_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu11_SQ14_BMTF", &L1_SingleMu11_SQ14_BMTF, &b_L1_SingleMu11_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF, &b_L1_SingleMu12_DQ_BMTF);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF, &b_L1_SingleMu12_DQ_EMTF);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF, &b_L1_SingleMu12_DQ_OMTF);
   fChain->SetBranchAddress("L1_SingleMu15_DQ", &L1_SingleMu15_DQ, &b_L1_SingleMu15_DQ);
   fChain->SetBranchAddress("L1_SingleMu18", &L1_SingleMu18, &b_L1_SingleMu18);
   fChain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   fChain->SetBranchAddress("L1_SingleMu22", &L1_SingleMu22, &b_L1_SingleMu22);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF, &b_L1_SingleMu22_BMTF);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF_NEG", &L1_SingleMu22_BMTF_NEG, &b_L1_SingleMu22_BMTF_NEG);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF_POS", &L1_SingleMu22_BMTF_POS, &b_L1_SingleMu22_BMTF_POS);
   fChain->SetBranchAddress("L1_SingleMu22_DQ", &L1_SingleMu22_DQ, &b_L1_SingleMu22_DQ);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF, &b_L1_SingleMu22_EMTF);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF_NEG", &L1_SingleMu22_EMTF_NEG, &b_L1_SingleMu22_EMTF_NEG);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF_POS", &L1_SingleMu22_EMTF_POS, &b_L1_SingleMu22_EMTF_POS);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF, &b_L1_SingleMu22_OMTF);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF_NEG", &L1_SingleMu22_OMTF_NEG, &b_L1_SingleMu22_OMTF_NEG);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF_POS", &L1_SingleMu22_OMTF_POS, &b_L1_SingleMu22_OMTF_POS);
   fChain->SetBranchAddress("L1_SingleMu22_OQ", &L1_SingleMu22_OQ, &b_L1_SingleMu22_OQ);
   fChain->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25, &b_L1_SingleMu25);
   fChain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   fChain->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5, &b_L1_SingleMu5);
   fChain->SetBranchAddress("L1_SingleMu5_SQ14_BMTF", &L1_SingleMu5_SQ14_BMTF, &b_L1_SingleMu5_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu6_SQ14_BMTF", &L1_SingleMu6_SQ14_BMTF, &b_L1_SingleMu6_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   fChain->SetBranchAddress("L1_SingleMu7_DQ", &L1_SingleMu7_DQ, &b_L1_SingleMu7_DQ);
   fChain->SetBranchAddress("L1_SingleMu7_SQ14_BMTF", &L1_SingleMu7_SQ14_BMTF, &b_L1_SingleMu7_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu8_SQ14_BMTF", &L1_SingleMu8_SQ14_BMTF, &b_L1_SingleMu8_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMu9_SQ14_BMTF", &L1_SingleMu9_SQ14_BMTF, &b_L1_SingleMu9_SQ14_BMTF);
   fChain->SetBranchAddress("L1_SingleMuCosmics", &L1_SingleMuCosmics, &b_L1_SingleMuCosmics);
   fChain->SetBranchAddress("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF, &b_L1_SingleMuCosmics_BMTF);
   fChain->SetBranchAddress("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF, &b_L1_SingleMuCosmics_EMTF);
   fChain->SetBranchAddress("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF, &b_L1_SingleMuCosmics_OMTF);
   fChain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   fChain->SetBranchAddress("L1_SingleMuOpen_BMTF", &L1_SingleMuOpen_BMTF, &b_L1_SingleMuOpen_BMTF);
   fChain->SetBranchAddress("L1_SingleMuOpen_EMTF", &L1_SingleMuOpen_EMTF, &b_L1_SingleMuOpen_EMTF);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR, &b_L1_SingleMuOpen_NotBptxOR);
   fChain->SetBranchAddress("L1_SingleMuOpen_OMTF", &L1_SingleMuOpen_OMTF, &b_L1_SingleMuOpen_OMTF);
   fChain->SetBranchAddress("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX, &b_L1_SingleMuOpen_er1p1_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX, &b_L1_SingleMuOpen_er1p4_NotBptxOR_3BX);
   fChain->SetBranchAddress("L1_SingleMuShower_Nominal", &L1_SingleMuShower_Nominal, &b_L1_SingleMuShower_Nominal);
   fChain->SetBranchAddress("L1_SingleMuShower_Tight", &L1_SingleMuShower_Tight, &b_L1_SingleMuShower_Tight);
   fChain->SetBranchAddress("L1_SingleTau120er2p1", &L1_SingleTau120er2p1, &b_L1_SingleTau120er2p1);
   fChain->SetBranchAddress("L1_SingleTau130er2p1", &L1_SingleTau130er2p1, &b_L1_SingleTau130er2p1);
   fChain->SetBranchAddress("L1_TOTEM_1", &L1_TOTEM_1, &b_L1_TOTEM_1);
   fChain->SetBranchAddress("L1_TOTEM_2", &L1_TOTEM_2, &b_L1_TOTEM_2);
   fChain->SetBranchAddress("L1_TOTEM_3", &L1_TOTEM_3, &b_L1_TOTEM_3);
   fChain->SetBranchAddress("L1_TOTEM_4", &L1_TOTEM_4, &b_L1_TOTEM_4);
   fChain->SetBranchAddress("L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5, &b_L1_TripleEG_18_17_8_er2p5);
   fChain->SetBranchAddress("L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5, &b_L1_TripleEG_18_18_12_er2p5);
   fChain->SetBranchAddress("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5, &b_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5);
   fChain->SetBranchAddress("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5, &b_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5);
   fChain->SetBranchAddress("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5, &b_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5);
   fChain->SetBranchAddress("L1_TripleMu0", &L1_TripleMu0, &b_L1_TripleMu0);
   fChain->SetBranchAddress("L1_TripleMu0_OQ", &L1_TripleMu0_OQ, &b_L1_TripleMu0_OQ);
   fChain->SetBranchAddress("L1_TripleMu0_SQ", &L1_TripleMu0_SQ, &b_L1_TripleMu0_SQ);
   fChain->SetBranchAddress("L1_TripleMu3", &L1_TripleMu3, &b_L1_TripleMu3);
   fChain->SetBranchAddress("L1_TripleMu3_SQ", &L1_TripleMu3_SQ, &b_L1_TripleMu3_SQ);
   fChain->SetBranchAddress("L1_TripleMu_3SQ_2p5SQ_0", &L1_TripleMu_3SQ_2p5SQ_0, &b_L1_TripleMu_3SQ_2p5SQ_0);
   fChain->SetBranchAddress("L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12", &L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12, &b_L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12);
   fChain->SetBranchAddress("L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12", &L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12, &b_L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12);
   fChain->SetBranchAddress("L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12", &L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12, &b_L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ, &b_L1_TripleMu_5SQ_3SQ_0OQ);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9);
   fChain->SetBranchAddress("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3, &b_L1_TripleMu_5_3_3);
   fChain->SetBranchAddress("L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ, &b_L1_TripleMu_5_3_3_SQ);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5, &b_L1_TripleMu_5_3p5_2p5);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   fChain->SetBranchAddress("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3, &b_L1_TripleMu_5_5_3);
   fChain->SetBranchAddress("L1_TwoMuShower_Loose", &L1_TwoMuShower_Loose, &b_L1_TwoMuShower_Loose);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus, &b_L1_UnpairedBunchBptxMinus);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus, &b_L1_UnpairedBunchBptxPlus);
   fChain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   fChain->SetBranchAddress("L1_ZeroBias_copy", &L1_ZeroBias_copy, &b_L1_ZeroBias_copy);
   fChain->SetBranchAddress("L1_UnprefireableEvent_TriggerRules", &L1_UnprefireableEvent_TriggerRules, &b_L1_UnprefireableEvent_TriggerRules);
   fChain->SetBranchAddress("L1_UnprefireableEvent_FirstBxInTrain", &L1_UnprefireableEvent_FirstBxInTrain, &b_L1_UnprefireableEvent_FirstBxInTrain);
   fChain->SetBranchAddress("L1_FinalOR_BXmin1", &L1_FinalOR_BXmin1, &b_L1_FinalOR_BXmin1);
   fChain->SetBranchAddress("L1_FinalOR_BXmin2", &L1_FinalOR_BXmin2, &b_L1_FinalOR_BXmin2);
   fChain->SetBranchAddress("L1Reco_step_pRECO", &L1Reco_step_pRECO, &b_L1Reco_step_pRECO);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter_pRECO", &Flag_HBHENoiseFilter_pRECO, &b_Flag_HBHENoiseFilter_pRECO);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter_pRECO", &Flag_HBHENoiseIsoFilter_pRECO, &b_Flag_HBHENoiseIsoFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter_pRECO", &Flag_CSCTightHaloFilter_pRECO, &b_Flag_CSCTightHaloFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO", &Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO, &b_Flag_CSCTightHaloTrkMuUnvetoFilter_pRECO);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter_pRECO", &Flag_CSCTightHalo2015Filter_pRECO, &b_Flag_CSCTightHalo2015Filter_pRECO);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter_pRECO", &Flag_globalTightHalo2016Filter_pRECO, &b_Flag_globalTightHalo2016Filter_pRECO);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter_pRECO", &Flag_globalSuperTightHalo2016Filter_pRECO, &b_Flag_globalSuperTightHalo2016Filter_pRECO);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter_pRECO", &Flag_HcalStripHaloFilter_pRECO, &b_Flag_HcalStripHaloFilter_pRECO);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter_pRECO", &Flag_hcalLaserEventFilter_pRECO, &b_Flag_hcalLaserEventFilter_pRECO);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO", &Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO, &b_Flag_EcalDeadCellTriggerPrimitiveFilter_pRECO);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter_pRECO", &Flag_EcalDeadCellBoundaryEnergyFilter_pRECO, &b_Flag_EcalDeadCellBoundaryEnergyFilter_pRECO);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter_pRECO", &Flag_ecalBadCalibFilter_pRECO, &b_Flag_ecalBadCalibFilter_pRECO);
   fChain->SetBranchAddress("Flag_goodVertices_pRECO", &Flag_goodVertices_pRECO, &b_Flag_goodVertices_pRECO);
   fChain->SetBranchAddress("Flag_eeBadScFilter_pRECO", &Flag_eeBadScFilter_pRECO, &b_Flag_eeBadScFilter_pRECO);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter_pRECO", &Flag_ecalLaserCorrFilter_pRECO, &b_Flag_ecalLaserCorrFilter_pRECO);
   fChain->SetBranchAddress("Flag_trkPOGFilters_pRECO", &Flag_trkPOGFilters_pRECO, &b_Flag_trkPOGFilters_pRECO);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter_pRECO", &Flag_chargedHadronTrackResolutionFilter_pRECO, &b_Flag_chargedHadronTrackResolutionFilter_pRECO);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter_pRECO", &Flag_muonBadTrackFilter_pRECO, &b_Flag_muonBadTrackFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter_pRECO", &Flag_BadChargedCandidateFilter_pRECO, &b_Flag_BadChargedCandidateFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter_pRECO", &Flag_BadPFMuonFilter_pRECO, &b_Flag_BadPFMuonFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonDzFilter_pRECO", &Flag_BadPFMuonDzFilter_pRECO, &b_Flag_BadPFMuonDzFilter_pRECO);
   fChain->SetBranchAddress("Flag_hfNoisyHitsFilter_pRECO", &Flag_hfNoisyHitsFilter_pRECO, &b_Flag_hfNoisyHitsFilter_pRECO);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter_pRECO", &Flag_BadChargedCandidateSummer16Filter_pRECO, &b_Flag_BadChargedCandidateSummer16Filter_pRECO);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter_pRECO", &Flag_BadPFMuonSummer16Filter_pRECO, &b_Flag_BadPFMuonSummer16Filter_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X_pRECO", &Flag_trkPOG_manystripclus53X_pRECO, &b_Flag_trkPOG_manystripclus53X_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X_pRECO", &Flag_trkPOG_toomanystripclus53X_pRECO, &b_Flag_trkPOG_toomanystripclus53X_pRECO);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters_pRECO", &Flag_trkPOG_logErrorTooManyClusters_pRECO, &b_Flag_trkPOG_logErrorTooManyClusters_pRECO);
   fChain->SetBranchAddress("Flag_METFilters_pRECO", &Flag_METFilters_pRECO, &b_Flag_METFilters_pRECO);
   fChain->SetBranchAddress("L1_AXO_Loose_pRECO", &L1_AXO_Loose_pRECO, &b_L1_AXO_Loose_pRECO);
   fChain->SetBranchAddress("L1_AXO_Nominal_pRECO", &L1_AXO_Nominal_pRECO, &b_L1_AXO_Nominal_pRECO);
   fChain->SetBranchAddress("L1_AXO_Tight_pRECO", &L1_AXO_Tight_pRECO, &b_L1_AXO_Tight_pRECO);
   fChain->SetBranchAddress("L1_AXO_VLoose_pRECO", &L1_AXO_VLoose_pRECO, &b_L1_AXO_VLoose_pRECO);
   fChain->SetBranchAddress("L1_AXO_VTight_pRECO", &L1_AXO_VTight_pRECO, &b_L1_AXO_VTight_pRECO);
   fChain->SetBranchAddress("L1_AlwaysTrue_pRECO", &L1_AlwaysTrue_pRECO, &b_L1_AlwaysTrue_pRECO);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref1_VME_pRECO", &L1_BPTX_AND_Ref1_VME_pRECO, &b_L1_BPTX_AND_Ref1_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref3_VME_pRECO", &L1_BPTX_AND_Ref3_VME_pRECO, &b_L1_BPTX_AND_Ref3_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_AND_Ref4_VME_pRECO", &L1_BPTX_AND_Ref4_VME_pRECO, &b_L1_BPTX_AND_Ref4_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B1_VME_pRECO", &L1_BPTX_BeamGas_B1_VME_pRECO, &b_L1_BPTX_BeamGas_B1_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_B2_VME_pRECO", &L1_BPTX_BeamGas_B2_VME_pRECO, &b_L1_BPTX_BeamGas_B2_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref1_VME_pRECO", &L1_BPTX_BeamGas_Ref1_VME_pRECO, &b_L1_BPTX_BeamGas_Ref1_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_BeamGas_Ref2_VME_pRECO", &L1_BPTX_BeamGas_Ref2_VME_pRECO, &b_L1_BPTX_BeamGas_Ref2_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_NotOR_VME_pRECO", &L1_BPTX_NotOR_VME_pRECO, &b_L1_BPTX_NotOR_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref3_VME_pRECO", &L1_BPTX_OR_Ref3_VME_pRECO, &b_L1_BPTX_OR_Ref3_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_OR_Ref4_VME_pRECO", &L1_BPTX_OR_Ref4_VME_pRECO, &b_L1_BPTX_OR_Ref4_VME_pRECO);
   fChain->SetBranchAddress("L1_BPTX_RefAND_VME_pRECO", &L1_BPTX_RefAND_VME_pRECO, &b_L1_BPTX_RefAND_VME_pRECO);
   fChain->SetBranchAddress("L1_BptxMinus_pRECO", &L1_BptxMinus_pRECO, &b_L1_BptxMinus_pRECO);
   fChain->SetBranchAddress("L1_BptxOR_pRECO", &L1_BptxOR_pRECO, &b_L1_BptxOR_pRECO);
   fChain->SetBranchAddress("L1_BptxPlus_pRECO", &L1_BptxPlus_pRECO, &b_L1_BptxPlus_pRECO);
   fChain->SetBranchAddress("L1_BptxXOR_pRECO", &L1_BptxXOR_pRECO, &b_L1_BptxXOR_pRECO);
   fChain->SetBranchAddress("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_pRECO", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_pRECO, &b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_pRECO);
   fChain->SetBranchAddress("L1_CICADA_Loose_pRECO", &L1_CICADA_Loose_pRECO, &b_L1_CICADA_Loose_pRECO);
   fChain->SetBranchAddress("L1_CICADA_Medium_pRECO", &L1_CICADA_Medium_pRECO, &b_L1_CICADA_Medium_pRECO);
   fChain->SetBranchAddress("L1_CICADA_Tight_pRECO", &L1_CICADA_Tight_pRECO, &b_L1_CICADA_Tight_pRECO);
   fChain->SetBranchAddress("L1_CICADA_VLoose_pRECO", &L1_CICADA_VLoose_pRECO, &b_L1_CICADA_VLoose_pRECO);
   fChain->SetBranchAddress("L1_CICADA_VTight_pRECO", &L1_CICADA_VTight_pRECO, &b_L1_CICADA_VTight_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG11_er1p2_dR_Max0p6_pRECO", &L1_DoubleEG11_er1p2_dR_Max0p6_pRECO, &b_L1_DoubleEG11_er1p2_dR_Max0p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT280er_pRECO", &L1_DoubleEG8er2p5_HTT280er_pRECO, &b_L1_DoubleEG8er2p5_HTT280er_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT300er_pRECO", &L1_DoubleEG8er2p5_HTT300er_pRECO, &b_L1_DoubleEG8er2p5_HTT300er_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG8er2p5_HTT320er_pRECO", &L1_DoubleEG8er2p5_HTT320er_pRECO, &b_L1_DoubleEG8er2p5_HTT320er_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_15_10_er2p5_pRECO", &L1_DoubleEG_15_10_er2p5_pRECO, &b_L1_DoubleEG_15_10_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_20_10_er2p5_pRECO", &L1_DoubleEG_20_10_er2p5_pRECO, &b_L1_DoubleEG_20_10_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_22_10_er2p5_pRECO", &L1_DoubleEG_22_10_er2p5_pRECO, &b_L1_DoubleEG_22_10_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_25_12_er2p5_pRECO", &L1_DoubleEG_25_12_er2p5_pRECO, &b_L1_DoubleEG_25_12_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_25_14_er2p5_pRECO", &L1_DoubleEG_25_14_er2p5_pRECO, &b_L1_DoubleEG_25_14_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_27_14_er2p5_pRECO", &L1_DoubleEG_27_14_er2p5_pRECO, &b_L1_DoubleEG_27_14_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso16_LooseIso12_er1p5_pRECO", &L1_DoubleEG_LooseIso16_LooseIso12_er1p5_pRECO, &b_L1_DoubleEG_LooseIso16_LooseIso12_er1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso18_LooseIso12_er1p5_pRECO", &L1_DoubleEG_LooseIso18_LooseIso12_er1p5_pRECO, &b_L1_DoubleEG_LooseIso18_LooseIso12_er1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso20_LooseIso12_er1p5_pRECO", &L1_DoubleEG_LooseIso20_LooseIso12_er1p5_pRECO, &b_L1_DoubleEG_LooseIso20_LooseIso12_er1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso22_12_er2p5_pRECO", &L1_DoubleEG_LooseIso22_12_er2p5_pRECO, &b_L1_DoubleEG_LooseIso22_12_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso22_LooseIso12_er1p5_pRECO", &L1_DoubleEG_LooseIso22_LooseIso12_er1p5_pRECO, &b_L1_DoubleEG_LooseIso22_LooseIso12_er1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso25_12_er2p5_pRECO", &L1_DoubleEG_LooseIso25_12_er2p5_pRECO, &b_L1_DoubleEG_LooseIso25_12_er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleEG_LooseIso25_LooseIso12_er1p5_pRECO", &L1_DoubleEG_LooseIso25_LooseIso12_er1p5_pRECO, &b_L1_DoubleEG_LooseIso25_LooseIso12_er1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5_pRECO", &L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5_pRECO, &b_L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5_pRECO", &L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5_pRECO, &b_L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau28er2p1_pRECO", &L1_DoubleIsoTau28er2p1_pRECO, &b_L1_DoubleIsoTau28er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau30er2p1_pRECO", &L1_DoubleIsoTau30er2p1_pRECO, &b_L1_DoubleIsoTau30er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau32er2p1_pRECO", &L1_DoubleIsoTau32er2p1_pRECO, &b_L1_DoubleIsoTau32er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau32er2p1_Mass_Max80_pRECO", &L1_DoubleIsoTau32er2p1_Mass_Max80_pRECO, &b_L1_DoubleIsoTau32er2p1_Mass_Max80_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau34er2p1_pRECO", &L1_DoubleIsoTau34er2p1_pRECO, &b_L1_DoubleIsoTau34er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau35er2p1_pRECO", &L1_DoubleIsoTau35er2p1_pRECO, &b_L1_DoubleIsoTau35er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleIsoTau36er2p1_pRECO", &L1_DoubleIsoTau36er2p1_pRECO, &b_L1_DoubleIsoTau36er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet100er2p3_dEta_Max1p6_pRECO", &L1_DoubleJet100er2p3_dEta_Max1p6_pRECO, &b_L1_DoubleJet100er2p3_dEta_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet100er2p5_pRECO", &L1_DoubleJet100er2p5_pRECO, &b_L1_DoubleJet100er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet112er2p3_dEta_Max1p6_pRECO", &L1_DoubleJet112er2p3_dEta_Max1p6_pRECO, &b_L1_DoubleJet112er2p3_dEta_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet120er2p5_pRECO", &L1_DoubleJet120er2p5_pRECO, &b_L1_DoubleJet120er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet120er2p5_Mu3_dR_Max0p8_pRECO", &L1_DoubleJet120er2p5_Mu3_dR_Max0p8_pRECO, &b_L1_DoubleJet120er2p5_Mu3_dR_Max0p8_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet150er2p5_pRECO", &L1_DoubleJet150er2p5_pRECO, &b_L1_DoubleJet150er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet16er2p5_Mu3_dR_Max0p4_pRECO", &L1_DoubleJet16er2p5_Mu3_dR_Max0p4_pRECO, &b_L1_DoubleJet16er2p5_Mu3_dR_Max0p4_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5_pRECO", &L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5_pRECO, &b_L1_DoubleJet30er2p5_Mass_Min225_dEta_Max1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5_pRECO", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5_pRECO, &b_L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5_pRECO", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5_pRECO, &b_L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5_pRECO", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5_pRECO, &b_L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet35er2p5_Mu3_dR_Max0p4_pRECO", &L1_DoubleJet35er2p5_Mu3_dR_Max0p4_pRECO, &b_L1_DoubleJet35er2p5_Mu3_dR_Max0p4_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet40er2p5_pRECO", &L1_DoubleJet40er2p5_pRECO, &b_L1_DoubleJet40er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5_pRECO", &L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5_pRECO, &b_L1_DoubleJet45_Mass_Min550_IsoTau45er2p1_RmOvlp_dR0p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO", &L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO, &b_L1_DoubleJet45_Mass_Min550_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5_pRECO", &L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5_pRECO, &b_L1_DoubleJet45_Mass_Min600_IsoTau45er2p1_RmOvlp_dR0p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO", &L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO, &b_L1_DoubleJet45_Mass_Min600_LooseIsoEG20er2p1_RmOvlp_dR0p2_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet60er2p5_Mu3_dR_Max0p4_pRECO", &L1_DoubleJet60er2p5_Mu3_dR_Max0p4_pRECO, &b_L1_DoubleJet60er2p5_Mu3_dR_Max0p4_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet80er2p5_Mu3_dR_Max0p4_pRECO", &L1_DoubleJet80er2p5_Mu3_dR_Max0p4_pRECO, &b_L1_DoubleJet80er2p5_Mu3_dR_Max0p4_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min800_pRECO", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min800_pRECO, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min800_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min850_pRECO", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min850_pRECO, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min850_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50_pRECO", &L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50_pRECO, &b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min600_DoubleJetCentral50_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50_pRECO", &L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50_pRECO, &b_L1_DoubleJet_65_35_DoubleJet35_Mass_Min650_DoubleJetCentral50_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65_pRECO", &L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65_pRECO, &b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min500_ETMHF65_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65_pRECO", &L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65_pRECO, &b_L1_DoubleJet_70_35_DoubleJet35_Mass_Min550_ETMHF65_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ_pRECO", &L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ_pRECO, &b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min600_Mu3OQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ_pRECO", &L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ_pRECO, &b_L1_DoubleJet_85_35_DoubleJet35_Mass_Min650_Mu3OQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleLLPJet40_pRECO", &L1_DoubleLLPJet40_pRECO, &b_L1_DoubleLLPJet40_pRECO);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG22er2p1_pRECO", &L1_DoubleLooseIsoEG22er2p1_pRECO, &b_L1_DoubleLooseIsoEG22er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleLooseIsoEG24er2p1_pRECO", &L1_DoubleLooseIsoEG24er2p1_pRECO, &b_L1_DoubleLooseIsoEG24er2p1_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_pRECO", &L1_DoubleMu0_pRECO, &b_L1_DoubleMu0_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_Mass_Min1_pRECO", &L1_DoubleMu0_Mass_Min1_pRECO, &b_L1_DoubleMu0_Mass_Min1_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_OQ_pRECO", &L1_DoubleMu0_OQ_pRECO, &b_L1_DoubleMu0_OQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ_pRECO", &L1_DoubleMu0_SQ_pRECO, &b_L1_DoubleMu0_SQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_SQ_OS_pRECO", &L1_DoubleMu0_SQ_OS_pRECO, &b_L1_DoubleMu0_SQ_OS_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt15_Upt7_pRECO", &L1_DoubleMu0_Upt15_Upt7_pRECO, &b_L1_DoubleMu0_Upt15_Upt7_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt5_Upt5_pRECO", &L1_DoubleMu0_Upt5_Upt5_pRECO, &b_L1_DoubleMu0_Upt5_Upt5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt6_IP_Min1_Upt4_pRECO", &L1_DoubleMu0_Upt6_IP_Min1_Upt4_pRECO, &b_L1_DoubleMu0_Upt6_IP_Min1_Upt4_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt6_SQ_er2p0_pRECO", &L1_DoubleMu0_Upt6_SQ_er2p0_pRECO, &b_L1_DoubleMu0_Upt6_SQ_er2p0_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt7_SQ_er2p0_pRECO", &L1_DoubleMu0_Upt7_SQ_er2p0_pRECO, &b_L1_DoubleMu0_Upt7_SQ_er2p0_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_Upt8_SQ_er2p0_pRECO", &L1_DoubleMu0_Upt8_SQ_er2p0_pRECO, &b_L1_DoubleMu0_Upt8_SQ_er2p0_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8_pRECO", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8_pRECO, &b_L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6_pRECO", &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6_pRECO, &b_L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2_pRECO", &L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2_pRECO, &b_L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_pRECO", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_pRECO, &b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_pRECO", &L1_DoubleMu0er1p5_SQ_pRECO, &b_L1_DoubleMu0er1p5_SQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_pRECO", &L1_DoubleMu0er1p5_SQ_OS_pRECO, &b_L1_DoubleMu0er1p5_SQ_OS_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2_pRECO", &L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2_pRECO, &b_L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_pRECO", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_pRECO, &b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_dR_Max1p4_pRECO", &L1_DoubleMu0er1p5_SQ_dR_Max1p4_pRECO, &b_L1_DoubleMu0er1p5_SQ_dR_Max1p4_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2_pRECO", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2_pRECO, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max0p3_dPhi_0p8to1p2_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5_pRECO", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5_pRECO, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6_pRECO", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6_pRECO, &b_L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dEta_Max1p5_pRECO", &L1_DoubleMu0er2p0_SQ_dEta_Max1p5_pRECO, &b_L1_DoubleMu0er2p0_SQ_dEta_Max1p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu0er2p0_SQ_dEta_Max1p6_pRECO", &L1_DoubleMu0er2p0_SQ_dEta_Max1p6_pRECO, &b_L1_DoubleMu0er2p0_SQ_dEta_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu18er2p1_SQ_pRECO", &L1_DoubleMu18er2p1_SQ_pRECO, &b_L1_DoubleMu18er2p1_SQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20_pRECO", &L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20_pRECO, &b_L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF30_HTT60er_pRECO", &L1_DoubleMu3_SQ_ETMHF30_HTT60er_pRECO, &b_L1_DoubleMu3_SQ_ETMHF30_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5_pRECO", &L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5_pRECO, &b_L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_HTT60er_pRECO", &L1_DoubleMu3_SQ_ETMHF40_HTT60er_pRECO, &b_L1_DoubleMu3_SQ_ETMHF40_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5_pRECO", &L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5_pRECO, &b_L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_HTT60er_pRECO", &L1_DoubleMu3_SQ_ETMHF50_HTT60er_pRECO, &b_L1_DoubleMu3_SQ_ETMHF50_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_pRECO", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_pRECO, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5_pRECO", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5_pRECO, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5_pRECO", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5_pRECO, &b_L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3_SQ_HTT220er_pRECO", &L1_DoubleMu3_SQ_HTT220er_pRECO, &b_L1_DoubleMu3_SQ_HTT220er_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6_pRECO", &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6_pRECO, &b_L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_EG9er2p5_pRECO", &L1_DoubleMu4_SQ_EG9er2p5_pRECO, &b_L1_DoubleMu4_SQ_EG9er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS_pRECO", &L1_DoubleMu4_SQ_OS_pRECO, &b_L1_DoubleMu4_SQ_OS_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4_SQ_OS_dR_Max1p2_pRECO", &L1_DoubleMu4_SQ_OS_dR_Max1p2_pRECO, &b_L1_DoubleMu4_SQ_OS_dR_Max1p2_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6_pRECO", &L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6_pRECO, &b_L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_pRECO", &L1_DoubleMu4p5_SQ_OS_pRECO, &b_L1_DoubleMu4p5_SQ_OS_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_dR_Max1p2_pRECO", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2_pRECO, &b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_pRECO", &L1_DoubleMu4p5er2p0_SQ_OS_pRECO, &b_L1_DoubleMu4p5er2p0_SQ_OS_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18_pRECO", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18_pRECO, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_pRECO", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_pRECO, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20_pRECO", &L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20_pRECO, &b_L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu5_SQ_EG9er2p5_pRECO", &L1_DoubleMu5_SQ_EG9er2p5_pRECO, &b_L1_DoubleMu5_SQ_EG9er2p5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu5_SQ_OS_dR_Max1p6_pRECO", &L1_DoubleMu5_SQ_OS_dR_Max1p6_pRECO, &b_L1_DoubleMu5_SQ_OS_dR_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu6_Upt6_SQ_er2p0_pRECO", &L1_DoubleMu6_Upt6_SQ_er2p0_pRECO, &b_L1_DoubleMu6_Upt6_SQ_er2p0_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu7_Upt7_SQ_er2p0_pRECO", &L1_DoubleMu7_Upt7_SQ_er2p0_pRECO, &b_L1_DoubleMu7_Upt7_SQ_er2p0_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu8_SQ_pRECO", &L1_DoubleMu8_SQ_pRECO, &b_L1_DoubleMu8_SQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu8_Upt8_SQ_er2p0_pRECO", &L1_DoubleMu8_Upt8_SQ_er2p0_pRECO, &b_L1_DoubleMu8_Upt8_SQ_er2p0_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu9_SQ_pRECO", &L1_DoubleMu9_SQ_pRECO, &b_L1_DoubleMu9_SQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu_12_5_pRECO", &L1_DoubleMu_12_5_pRECO, &b_L1_DoubleMu_12_5_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu_15_5_SQ_pRECO", &L1_DoubleMu_15_5_SQ_pRECO, &b_L1_DoubleMu_15_5_SQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_pRECO", &L1_DoubleMu_15_7_pRECO, &b_L1_DoubleMu_15_7_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_Mass_Min1_pRECO", &L1_DoubleMu_15_7_Mass_Min1_pRECO, &b_L1_DoubleMu_15_7_Mass_Min1_pRECO);
   fChain->SetBranchAddress("L1_DoubleMu_15_7_SQ_pRECO", &L1_DoubleMu_15_7_SQ_pRECO, &b_L1_DoubleMu_15_7_SQ_pRECO);
   fChain->SetBranchAddress("L1_DoubleTau70er2p1_pRECO", &L1_DoubleTau70er2p1_pRECO, &b_L1_DoubleTau70er2p1_pRECO);
   fChain->SetBranchAddress("L1_ETM120_pRECO", &L1_ETM120_pRECO, &b_L1_ETM120_pRECO);
   fChain->SetBranchAddress("L1_ETM150_pRECO", &L1_ETM150_pRECO, &b_L1_ETM150_pRECO);
   fChain->SetBranchAddress("L1_ETMHF100_pRECO", &L1_ETMHF100_pRECO, &b_L1_ETMHF100_pRECO);
   fChain->SetBranchAddress("L1_ETMHF100_HTT60er_pRECO", &L1_ETMHF100_HTT60er_pRECO, &b_L1_ETMHF100_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_ETMHF110_pRECO", &L1_ETMHF110_pRECO, &b_L1_ETMHF110_pRECO);
   fChain->SetBranchAddress("L1_ETMHF110_HTT60er_pRECO", &L1_ETMHF110_HTT60er_pRECO, &b_L1_ETMHF110_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_ETMHF120_pRECO", &L1_ETMHF120_pRECO, &b_L1_ETMHF120_pRECO);
   fChain->SetBranchAddress("L1_ETMHF120_HTT60er_pRECO", &L1_ETMHF120_HTT60er_pRECO, &b_L1_ETMHF120_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_ETMHF130_pRECO", &L1_ETMHF130_pRECO, &b_L1_ETMHF130_pRECO);
   fChain->SetBranchAddress("L1_ETMHF130_HTT60er_pRECO", &L1_ETMHF130_HTT60er_pRECO, &b_L1_ETMHF130_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_ETMHF140_pRECO", &L1_ETMHF140_pRECO, &b_L1_ETMHF140_pRECO);
   fChain->SetBranchAddress("L1_ETMHF150_pRECO", &L1_ETMHF150_pRECO, &b_L1_ETMHF150_pRECO);
   fChain->SetBranchAddress("L1_ETMHF70_pRECO", &L1_ETMHF70_pRECO, &b_L1_ETMHF70_pRECO);
   fChain->SetBranchAddress("L1_ETMHF70_HTT60er_pRECO", &L1_ETMHF70_HTT60er_pRECO, &b_L1_ETMHF70_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_ETMHF80_pRECO", &L1_ETMHF80_pRECO, &b_L1_ETMHF80_pRECO);
   fChain->SetBranchAddress("L1_ETMHF80_HTT60er_pRECO", &L1_ETMHF80_HTT60er_pRECO, &b_L1_ETMHF80_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1_pRECO", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1_pRECO, &b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1_pRECO);
   fChain->SetBranchAddress("L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6_pRECO", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6_pRECO, &b_L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6_pRECO);
   fChain->SetBranchAddress("L1_ETMHF90_pRECO", &L1_ETMHF90_pRECO, &b_L1_ETMHF90_pRECO);
   fChain->SetBranchAddress("L1_ETMHF90_HTT60er_pRECO", &L1_ETMHF90_HTT60er_pRECO, &b_L1_ETMHF90_HTT60er_pRECO);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1_pRECO", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1_pRECO, &b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1_pRECO);
   fChain->SetBranchAddress("L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6_pRECO", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6_pRECO, &b_L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6_pRECO);
   fChain->SetBranchAddress("L1_ETT2000_pRECO", &L1_ETT2000_pRECO, &b_L1_ETT2000_pRECO);
   fChain->SetBranchAddress("L1_FirstBunchAfterTrain_pRECO", &L1_FirstBunchAfterTrain_pRECO, &b_L1_FirstBunchAfterTrain_pRECO);
   fChain->SetBranchAddress("L1_FirstBunchBeforeTrain_pRECO", &L1_FirstBunchBeforeTrain_pRECO, &b_L1_FirstBunchBeforeTrain_pRECO);
   fChain->SetBranchAddress("L1_FirstBunchInTrain_pRECO", &L1_FirstBunchInTrain_pRECO, &b_L1_FirstBunchInTrain_pRECO);
   fChain->SetBranchAddress("L1_FirstCollisionInOrbit_pRECO", &L1_FirstCollisionInOrbit_pRECO, &b_L1_FirstCollisionInOrbit_pRECO);
   fChain->SetBranchAddress("L1_FirstCollisionInTrain_pRECO", &L1_FirstCollisionInTrain_pRECO, &b_L1_FirstCollisionInTrain_pRECO);
   fChain->SetBranchAddress("L1_HCAL_LaserMon_Trig_pRECO", &L1_HCAL_LaserMon_Trig_pRECO, &b_L1_HCAL_LaserMon_Trig_pRECO);
   fChain->SetBranchAddress("L1_HCAL_LaserMon_Veto_pRECO", &L1_HCAL_LaserMon_Veto_pRECO, &b_L1_HCAL_LaserMon_Veto_pRECO);
   fChain->SetBranchAddress("L1_HTMHF100_pRECO", &L1_HTMHF100_pRECO, &b_L1_HTMHF100_pRECO);
   fChain->SetBranchAddress("L1_HTMHF120_pRECO", &L1_HTMHF120_pRECO, &b_L1_HTMHF120_pRECO);
   fChain->SetBranchAddress("L1_HTMHF125_pRECO", &L1_HTMHF125_pRECO, &b_L1_HTMHF125_pRECO);
   fChain->SetBranchAddress("L1_HTMHF130_pRECO", &L1_HTMHF130_pRECO, &b_L1_HTMHF130_pRECO);
   fChain->SetBranchAddress("L1_HTMHF150_pRECO", &L1_HTMHF150_pRECO, &b_L1_HTMHF150_pRECO);
   fChain->SetBranchAddress("L1_HTT120_SingleLLPJet40_pRECO", &L1_HTT120_SingleLLPJet40_pRECO, &b_L1_HTT120_SingleLLPJet40_pRECO);
   fChain->SetBranchAddress("L1_HTT120er_pRECO", &L1_HTT120er_pRECO, &b_L1_HTT120er_pRECO);
   fChain->SetBranchAddress("L1_HTT160_SingleLLPJet50_pRECO", &L1_HTT160_SingleLLPJet50_pRECO, &b_L1_HTT160_SingleLLPJet50_pRECO);
   fChain->SetBranchAddress("L1_HTT160er_pRECO", &L1_HTT160er_pRECO, &b_L1_HTT160er_pRECO);
   fChain->SetBranchAddress("L1_HTT200_SingleLLPJet60_pRECO", &L1_HTT200_SingleLLPJet60_pRECO, &b_L1_HTT200_SingleLLPJet60_pRECO);
   fChain->SetBranchAddress("L1_HTT200er_pRECO", &L1_HTT200er_pRECO, &b_L1_HTT200er_pRECO);
   fChain->SetBranchAddress("L1_HTT240_SingleLLPJet70_pRECO", &L1_HTT240_SingleLLPJet70_pRECO, &b_L1_HTT240_SingleLLPJet70_pRECO);
   fChain->SetBranchAddress("L1_HTT255er_pRECO", &L1_HTT255er_pRECO, &b_L1_HTT255er_pRECO);
   fChain->SetBranchAddress("L1_HTT280er_pRECO", &L1_HTT280er_pRECO, &b_L1_HTT280er_pRECO);
   fChain->SetBranchAddress("L1_HTT280er_QuadJet_70_55_40_35_er2p5_pRECO", &L1_HTT280er_QuadJet_70_55_40_35_er2p5_pRECO, &b_L1_HTT280er_QuadJet_70_55_40_35_er2p5_pRECO);
   fChain->SetBranchAddress("L1_HTT320er_pRECO", &L1_HTT320er_pRECO, &b_L1_HTT320er_pRECO);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p5_pRECO", &L1_HTT320er_QuadJet_70_55_40_40_er2p5_pRECO, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p5_pRECO);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3_pRECO", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3_pRECO, &b_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3_pRECO);
   fChain->SetBranchAddress("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3_pRECO", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3_pRECO, &b_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3_pRECO);
   fChain->SetBranchAddress("L1_HTT360er_pRECO", &L1_HTT360er_pRECO, &b_L1_HTT360er_pRECO);
   fChain->SetBranchAddress("L1_HTT400er_pRECO", &L1_HTT400er_pRECO, &b_L1_HTT400er_pRECO);
   fChain->SetBranchAddress("L1_HTT450er_pRECO", &L1_HTT450er_pRECO, &b_L1_HTT450er_pRECO);
   fChain->SetBranchAddress("L1_IsoEG32er2p5_Mt40_pRECO", &L1_IsoEG32er2p5_Mt40_pRECO, &b_L1_IsoEG32er2p5_Mt40_pRECO);
   fChain->SetBranchAddress("L1_IsoTau52er2p1_QuadJet36er2p5_pRECO", &L1_IsoTau52er2p1_QuadJet36er2p5_pRECO, &b_L1_IsoTau52er2p1_QuadJet36er2p5_pRECO);
   fChain->SetBranchAddress("L1_IsolatedBunch_pRECO", &L1_IsolatedBunch_pRECO, &b_L1_IsolatedBunch_pRECO);
   fChain->SetBranchAddress("L1_LastBunchInTrain_pRECO", &L1_LastBunchInTrain_pRECO, &b_L1_LastBunchInTrain_pRECO);
   fChain->SetBranchAddress("L1_LastCollisionInTrain_pRECO", &L1_LastCollisionInTrain_pRECO, &b_L1_LastCollisionInTrain_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG14er2p5_HTT200er_pRECO", &L1_LooseIsoEG14er2p5_HTT200er_pRECO, &b_L1_LooseIsoEG14er2p5_HTT200er_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG16er2p5_HTT200er_pRECO", &L1_LooseIsoEG16er2p5_HTT200er_pRECO, &b_L1_LooseIsoEG16er2p5_HTT200er_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3_pRECO", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3_pRECO, &b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3_pRECO", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3_pRECO, &b_L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_HTT100er_pRECO", &L1_LooseIsoEG24er2p1_HTT100er_pRECO, &b_L1_LooseIsoEG24er2p1_HTT100er_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3_pRECO", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3_pRECO, &b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_HTT100er_pRECO", &L1_LooseIsoEG26er2p1_HTT100er_pRECO, &b_L1_LooseIsoEG26er2p1_HTT100er_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_pRECO", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_pRECO, &b_L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_HTT100er_pRECO", &L1_LooseIsoEG28er2p1_HTT100er_pRECO, &b_L1_LooseIsoEG28er2p1_HTT100er_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_pRECO", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_pRECO, &b_L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3_pRECO);
   fChain->SetBranchAddress("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_pRECO", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_pRECO, &b_L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3_pRECO);
   fChain->SetBranchAddress("L1_MinimumBiasHF0_pRECO", &L1_MinimumBiasHF0_pRECO, &b_L1_MinimumBiasHF0_pRECO);
   fChain->SetBranchAddress("L1_MinimumBiasHF0_AND_BptxAND_pRECO", &L1_MinimumBiasHF0_AND_BptxAND_pRECO, &b_L1_MinimumBiasHF0_AND_BptxAND_pRECO);
   fChain->SetBranchAddress("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6_pRECO", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6_pRECO, &b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_Mu12_HTT150er_pRECO", &L1_Mu12_HTT150er_pRECO, &b_L1_Mu12_HTT150er_pRECO);
   fChain->SetBranchAddress("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6_pRECO", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6_pRECO, &b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6_pRECO);
   fChain->SetBranchAddress("L1_Mu14_HTT150er_pRECO", &L1_Mu14_HTT150er_pRECO, &b_L1_Mu14_HTT150er_pRECO);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau24er2p1_pRECO", &L1_Mu18er2p1_Tau24er2p1_pRECO, &b_L1_Mu18er2p1_Tau24er2p1_pRECO);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1_pRECO", &L1_Mu18er2p1_Tau26er2p1_pRECO, &b_L1_Mu18er2p1_Tau26er2p1_pRECO);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1_Jet55_pRECO", &L1_Mu18er2p1_Tau26er2p1_Jet55_pRECO, &b_L1_Mu18er2p1_Tau26er2p1_Jet55_pRECO);
   fChain->SetBranchAddress("L1_Mu18er2p1_Tau26er2p1_Jet70_pRECO", &L1_Mu18er2p1_Tau26er2p1_Jet70_pRECO, &b_L1_Mu18er2p1_Tau26er2p1_Jet70_pRECO);
   fChain->SetBranchAddress("L1_Mu20_EG10er2p5_pRECO", &L1_Mu20_EG10er2p5_pRECO, &b_L1_Mu20_EG10er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau30er2p1_pRECO", &L1_Mu22er2p1_IsoTau30er2p1_pRECO, &b_L1_Mu22er2p1_IsoTau30er2p1_pRECO);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau32er2p1_pRECO", &L1_Mu22er2p1_IsoTau32er2p1_pRECO, &b_L1_Mu22er2p1_IsoTau32er2p1_pRECO);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau34er2p1_pRECO", &L1_Mu22er2p1_IsoTau34er2p1_pRECO, &b_L1_Mu22er2p1_IsoTau34er2p1_pRECO);
   fChain->SetBranchAddress("L1_Mu22er2p1_IsoTau40er2p1_pRECO", &L1_Mu22er2p1_IsoTau40er2p1_pRECO, &b_L1_Mu22er2p1_IsoTau40er2p1_pRECO);
   fChain->SetBranchAddress("L1_Mu22er2p1_Tau70er2p1_pRECO", &L1_Mu22er2p1_Tau70er2p1_pRECO, &b_L1_Mu22er2p1_Tau70er2p1_pRECO);
   fChain->SetBranchAddress("L1_Mu3_Jet120er2p5_dR_Max0p4_pRECO", &L1_Mu3_Jet120er2p5_dR_Max0p4_pRECO, &b_L1_Mu3_Jet120er2p5_dR_Max0p4_pRECO);
   fChain->SetBranchAddress("L1_Mu3_Jet16er2p5_dR_Max0p4_pRECO", &L1_Mu3_Jet16er2p5_dR_Max0p4_pRECO, &b_L1_Mu3_Jet16er2p5_dR_Max0p4_pRECO);
   fChain->SetBranchAddress("L1_Mu3_Jet30er2p5_pRECO", &L1_Mu3_Jet30er2p5_pRECO, &b_L1_Mu3_Jet30er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu3_Jet60er2p5_dR_Max0p4_pRECO", &L1_Mu3_Jet60er2p5_dR_Max0p4_pRECO, &b_L1_Mu3_Jet60er2p5_dR_Max0p4_pRECO);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF30_pRECO", &L1_Mu3er1p5_Jet100er2p5_ETMHF30_pRECO, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF30_pRECO);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF40_pRECO", &L1_Mu3er1p5_Jet100er2p5_ETMHF40_pRECO, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF40_pRECO);
   fChain->SetBranchAddress("L1_Mu3er1p5_Jet100er2p5_ETMHF50_pRECO", &L1_Mu3er1p5_Jet100er2p5_ETMHF50_pRECO, &b_L1_Mu3er1p5_Jet100er2p5_ETMHF50_pRECO);
   fChain->SetBranchAddress("L1_Mu5_EG23er2p5_pRECO", &L1_Mu5_EG23er2p5_pRECO, &b_L1_Mu5_EG23er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu5_LooseIsoEG20er2p5_pRECO", &L1_Mu5_LooseIsoEG20er2p5_pRECO, &b_L1_Mu5_LooseIsoEG20er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG10er2p5_pRECO", &L1_Mu6_DoubleEG10er2p5_pRECO, &b_L1_Mu6_DoubleEG10er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG12er2p5_pRECO", &L1_Mu6_DoubleEG12er2p5_pRECO, &b_L1_Mu6_DoubleEG12er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG15er2p5_pRECO", &L1_Mu6_DoubleEG15er2p5_pRECO, &b_L1_Mu6_DoubleEG15er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu6_DoubleEG17er2p5_pRECO", &L1_Mu6_DoubleEG17er2p5_pRECO, &b_L1_Mu6_DoubleEG17er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu6_HTT240er_pRECO", &L1_Mu6_HTT240er_pRECO, &b_L1_Mu6_HTT240er_pRECO);
   fChain->SetBranchAddress("L1_Mu6_HTT250er_pRECO", &L1_Mu6_HTT250er_pRECO, &b_L1_Mu6_HTT250er_pRECO);
   fChain->SetBranchAddress("L1_Mu7_EG20er2p5_pRECO", &L1_Mu7_EG20er2p5_pRECO, &b_L1_Mu7_EG20er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu7_EG23er2p5_pRECO", &L1_Mu7_EG23er2p5_pRECO, &b_L1_Mu7_EG23er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG20er2p5_pRECO", &L1_Mu7_LooseIsoEG20er2p5_pRECO, &b_L1_Mu7_LooseIsoEG20er2p5_pRECO);
   fChain->SetBranchAddress("L1_Mu7_LooseIsoEG23er2p5_pRECO", &L1_Mu7_LooseIsoEG23er2p5_pRECO, &b_L1_Mu7_LooseIsoEG23er2p5_pRECO);
   fChain->SetBranchAddress("L1_NotBptxOR_pRECO", &L1_NotBptxOR_pRECO, &b_L1_NotBptxOR_pRECO);
   fChain->SetBranchAddress("L1_QuadJet60er2p5_pRECO", &L1_QuadJet60er2p5_pRECO, &b_L1_QuadJet60er2p5_pRECO);
   fChain->SetBranchAddress("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0_pRECO", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0_pRECO, &b_L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0_pRECO);
   fChain->SetBranchAddress("L1_QuadMu0_pRECO", &L1_QuadMu0_pRECO, &b_L1_QuadMu0_pRECO);
   fChain->SetBranchAddress("L1_QuadMu0_OQ_pRECO", &L1_QuadMu0_OQ_pRECO, &b_L1_QuadMu0_OQ_pRECO);
   fChain->SetBranchAddress("L1_QuadMu0_SQ_pRECO", &L1_QuadMu0_SQ_pRECO, &b_L1_QuadMu0_SQ_pRECO);
   fChain->SetBranchAddress("L1_SecondBunchInTrain_pRECO", &L1_SecondBunchInTrain_pRECO, &b_L1_SecondBunchInTrain_pRECO);
   fChain->SetBranchAddress("L1_SecondLastBunchInTrain_pRECO", &L1_SecondLastBunchInTrain_pRECO, &b_L1_SecondLastBunchInTrain_pRECO);
   fChain->SetBranchAddress("L1_SingleEG10er2p5_pRECO", &L1_SingleEG10er2p5_pRECO, &b_L1_SingleEG10er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG15er2p5_pRECO", &L1_SingleEG15er2p5_pRECO, &b_L1_SingleEG15er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG26er2p5_pRECO", &L1_SingleEG26er2p5_pRECO, &b_L1_SingleEG26er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG28_FWD2p5_pRECO", &L1_SingleEG28_FWD2p5_pRECO, &b_L1_SingleEG28_FWD2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG28er1p5_pRECO", &L1_SingleEG28er1p5_pRECO, &b_L1_SingleEG28er1p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG28er2p1_pRECO", &L1_SingleEG28er2p1_pRECO, &b_L1_SingleEG28er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleEG28er2p5_pRECO", &L1_SingleEG28er2p5_pRECO, &b_L1_SingleEG28er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG34er2p5_pRECO", &L1_SingleEG34er2p5_pRECO, &b_L1_SingleEG34er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG36er2p5_pRECO", &L1_SingleEG36er2p5_pRECO, &b_L1_SingleEG36er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG38er2p5_pRECO", &L1_SingleEG38er2p5_pRECO, &b_L1_SingleEG38er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG40er2p5_pRECO", &L1_SingleEG40er2p5_pRECO, &b_L1_SingleEG40er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG42er2p5_pRECO", &L1_SingleEG42er2p5_pRECO, &b_L1_SingleEG42er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG45er2p5_pRECO", &L1_SingleEG45er2p5_pRECO, &b_L1_SingleEG45er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleEG50_pRECO", &L1_SingleEG50_pRECO, &b_L1_SingleEG50_pRECO);
   fChain->SetBranchAddress("L1_SingleEG60_pRECO", &L1_SingleEG60_pRECO, &b_L1_SingleEG60_pRECO);
   fChain->SetBranchAddress("L1_SingleEG8er2p5_pRECO", &L1_SingleEG8er2p5_pRECO, &b_L1_SingleEG8er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG24er2p1_pRECO", &L1_SingleIsoEG24er2p1_pRECO, &b_L1_SingleIsoEG24er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG26er2p1_pRECO", &L1_SingleIsoEG26er2p1_pRECO, &b_L1_SingleIsoEG26er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG26er2p5_pRECO", &L1_SingleIsoEG26er2p5_pRECO, &b_L1_SingleIsoEG26er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG28_FWD2p5_pRECO", &L1_SingleIsoEG28_FWD2p5_pRECO, &b_L1_SingleIsoEG28_FWD2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG28er1p5_pRECO", &L1_SingleIsoEG28er1p5_pRECO, &b_L1_SingleIsoEG28er1p5_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG28er2p1_pRECO", &L1_SingleIsoEG28er2p1_pRECO, &b_L1_SingleIsoEG28er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG28er2p5_pRECO", &L1_SingleIsoEG28er2p5_pRECO, &b_L1_SingleIsoEG28er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG30er2p1_pRECO", &L1_SingleIsoEG30er2p1_pRECO, &b_L1_SingleIsoEG30er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG30er2p5_pRECO", &L1_SingleIsoEG30er2p5_pRECO, &b_L1_SingleIsoEG30er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG32er2p1_pRECO", &L1_SingleIsoEG32er2p1_pRECO, &b_L1_SingleIsoEG32er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG32er2p5_pRECO", &L1_SingleIsoEG32er2p5_pRECO, &b_L1_SingleIsoEG32er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleIsoEG34er2p5_pRECO", &L1_SingleIsoEG34er2p5_pRECO, &b_L1_SingleIsoEG34er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet10erHE_pRECO", &L1_SingleJet10erHE_pRECO, &b_L1_SingleJet10erHE_pRECO);
   fChain->SetBranchAddress("L1_SingleJet120_pRECO", &L1_SingleJet120_pRECO, &b_L1_SingleJet120_pRECO);
   fChain->SetBranchAddress("L1_SingleJet120_FWD2p5_pRECO", &L1_SingleJet120_FWD2p5_pRECO, &b_L1_SingleJet120_FWD2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet120_FWD3p0_pRECO", &L1_SingleJet120_FWD3p0_pRECO, &b_L1_SingleJet120_FWD3p0_pRECO);
   fChain->SetBranchAddress("L1_SingleJet120er1p3_pRECO", &L1_SingleJet120er1p3_pRECO, &b_L1_SingleJet120er1p3_pRECO);
   fChain->SetBranchAddress("L1_SingleJet120er2p5_pRECO", &L1_SingleJet120er2p5_pRECO, &b_L1_SingleJet120er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet12erHE_pRECO", &L1_SingleJet12erHE_pRECO, &b_L1_SingleJet12erHE_pRECO);
   fChain->SetBranchAddress("L1_SingleJet140er2p5_pRECO", &L1_SingleJet140er2p5_pRECO, &b_L1_SingleJet140er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet160er2p5_pRECO", &L1_SingleJet160er2p5_pRECO, &b_L1_SingleJet160er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet180_pRECO", &L1_SingleJet180_pRECO, &b_L1_SingleJet180_pRECO);
   fChain->SetBranchAddress("L1_SingleJet180er2p5_pRECO", &L1_SingleJet180er2p5_pRECO, &b_L1_SingleJet180er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet200_pRECO", &L1_SingleJet200_pRECO, &b_L1_SingleJet200_pRECO);
   fChain->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR_pRECO", &L1_SingleJet20er2p5_NotBptxOR_pRECO, &b_L1_SingleJet20er2p5_NotBptxOR_pRECO);
   fChain->SetBranchAddress("L1_SingleJet20er2p5_NotBptxOR_3BX_pRECO", &L1_SingleJet20er2p5_NotBptxOR_3BX_pRECO, &b_L1_SingleJet20er2p5_NotBptxOR_3BX_pRECO);
   fChain->SetBranchAddress("L1_SingleJet35_pRECO", &L1_SingleJet35_pRECO, &b_L1_SingleJet35_pRECO);
   fChain->SetBranchAddress("L1_SingleJet35_FWD2p5_pRECO", &L1_SingleJet35_FWD2p5_pRECO, &b_L1_SingleJet35_FWD2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet35_FWD3p0_pRECO", &L1_SingleJet35_FWD3p0_pRECO, &b_L1_SingleJet35_FWD3p0_pRECO);
   fChain->SetBranchAddress("L1_SingleJet35er1p3_pRECO", &L1_SingleJet35er1p3_pRECO, &b_L1_SingleJet35er1p3_pRECO);
   fChain->SetBranchAddress("L1_SingleJet35er2p5_pRECO", &L1_SingleJet35er2p5_pRECO, &b_L1_SingleJet35er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet43er2p5_NotBptxOR_3BX_pRECO", &L1_SingleJet43er2p5_NotBptxOR_3BX_pRECO, &b_L1_SingleJet43er2p5_NotBptxOR_3BX_pRECO);
   fChain->SetBranchAddress("L1_SingleJet46er2p5_NotBptxOR_3BX_pRECO", &L1_SingleJet46er2p5_NotBptxOR_3BX_pRECO, &b_L1_SingleJet46er2p5_NotBptxOR_3BX_pRECO);
   fChain->SetBranchAddress("L1_SingleJet60_pRECO", &L1_SingleJet60_pRECO, &b_L1_SingleJet60_pRECO);
   fChain->SetBranchAddress("L1_SingleJet60_FWD2p5_pRECO", &L1_SingleJet60_FWD2p5_pRECO, &b_L1_SingleJet60_FWD2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleJet8erHE_pRECO", &L1_SingleJet8erHE_pRECO, &b_L1_SingleJet8erHE_pRECO);
   fChain->SetBranchAddress("L1_SingleJet90_pRECO", &L1_SingleJet90_pRECO, &b_L1_SingleJet90_pRECO);
   fChain->SetBranchAddress("L1_SingleJet90_FWD2p5_pRECO", &L1_SingleJet90_FWD2p5_pRECO, &b_L1_SingleJet90_FWD2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG26er1p5_pRECO", &L1_SingleLooseIsoEG26er1p5_pRECO, &b_L1_SingleLooseIsoEG26er1p5_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG26er2p5_pRECO", &L1_SingleLooseIsoEG26er2p5_pRECO, &b_L1_SingleLooseIsoEG26er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28_FWD2p5_pRECO", &L1_SingleLooseIsoEG28_FWD2p5_pRECO, &b_L1_SingleLooseIsoEG28_FWD2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er1p5_pRECO", &L1_SingleLooseIsoEG28er1p5_pRECO, &b_L1_SingleLooseIsoEG28er1p5_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er2p1_pRECO", &L1_SingleLooseIsoEG28er2p1_pRECO, &b_L1_SingleLooseIsoEG28er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG28er2p5_pRECO", &L1_SingleLooseIsoEG28er2p5_pRECO, &b_L1_SingleLooseIsoEG28er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG30er1p5_pRECO", &L1_SingleLooseIsoEG30er1p5_pRECO, &b_L1_SingleLooseIsoEG30er1p5_pRECO);
   fChain->SetBranchAddress("L1_SingleLooseIsoEG30er2p5_pRECO", &L1_SingleLooseIsoEG30er2p5_pRECO, &b_L1_SingleLooseIsoEG30er2p5_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_BMTF_pRECO", &L1_SingleMu0_BMTF_pRECO, &b_L1_SingleMu0_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_DQ_pRECO", &L1_SingleMu0_DQ_pRECO, &b_L1_SingleMu0_DQ_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_EMTF_pRECO", &L1_SingleMu0_EMTF_pRECO, &b_L1_SingleMu0_EMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_OMTF_pRECO", &L1_SingleMu0_OMTF_pRECO, &b_L1_SingleMu0_OMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_SQ13_BMTF_pRECO", &L1_SingleMu0_SQ13_BMTF_pRECO, &b_L1_SingleMu0_SQ13_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_SQ14_BMTF_pRECO", &L1_SingleMu0_SQ14_BMTF_pRECO, &b_L1_SingleMu0_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_SQ15_BMTF_pRECO", &L1_SingleMu0_SQ15_BMTF_pRECO, &b_L1_SingleMu0_SQ15_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_pRECO", &L1_SingleMu0_Upt10_pRECO, &b_L1_SingleMu0_Upt10_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_BMTF_pRECO", &L1_SingleMu0_Upt10_BMTF_pRECO, &b_L1_SingleMu0_Upt10_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_EMTF_pRECO", &L1_SingleMu0_Upt10_EMTF_pRECO, &b_L1_SingleMu0_Upt10_EMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_OMTF_pRECO", &L1_SingleMu0_Upt10_OMTF_pRECO, &b_L1_SingleMu0_Upt10_OMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt10_SQ14_BMTF_pRECO", &L1_SingleMu0_Upt10_SQ14_BMTF_pRECO, &b_L1_SingleMu0_Upt10_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt15_SQ14_BMTF_pRECO", &L1_SingleMu0_Upt15_SQ14_BMTF_pRECO, &b_L1_SingleMu0_Upt15_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt20_SQ14_BMTF_pRECO", &L1_SingleMu0_Upt20_SQ14_BMTF_pRECO, &b_L1_SingleMu0_Upt20_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu0_Upt25_SQ14_BMTF_pRECO", &L1_SingleMu0_Upt25_SQ14_BMTF_pRECO, &b_L1_SingleMu0_Upt25_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu10_SQ14_BMTF_pRECO", &L1_SingleMu10_SQ14_BMTF_pRECO, &b_L1_SingleMu10_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu11_SQ14_BMTF_pRECO", &L1_SingleMu11_SQ14_BMTF_pRECO, &b_L1_SingleMu11_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_BMTF_pRECO", &L1_SingleMu12_DQ_BMTF_pRECO, &b_L1_SingleMu12_DQ_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_EMTF_pRECO", &L1_SingleMu12_DQ_EMTF_pRECO, &b_L1_SingleMu12_DQ_EMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu12_DQ_OMTF_pRECO", &L1_SingleMu12_DQ_OMTF_pRECO, &b_L1_SingleMu12_DQ_OMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu15_DQ_pRECO", &L1_SingleMu15_DQ_pRECO, &b_L1_SingleMu15_DQ_pRECO);
   fChain->SetBranchAddress("L1_SingleMu18_pRECO", &L1_SingleMu18_pRECO, &b_L1_SingleMu18_pRECO);
   fChain->SetBranchAddress("L1_SingleMu20_pRECO", &L1_SingleMu20_pRECO, &b_L1_SingleMu20_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_pRECO", &L1_SingleMu22_pRECO, &b_L1_SingleMu22_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF_pRECO", &L1_SingleMu22_BMTF_pRECO, &b_L1_SingleMu22_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF_NEG_pRECO", &L1_SingleMu22_BMTF_NEG_pRECO, &b_L1_SingleMu22_BMTF_NEG_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_BMTF_POS_pRECO", &L1_SingleMu22_BMTF_POS_pRECO, &b_L1_SingleMu22_BMTF_POS_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_DQ_pRECO", &L1_SingleMu22_DQ_pRECO, &b_L1_SingleMu22_DQ_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF_pRECO", &L1_SingleMu22_EMTF_pRECO, &b_L1_SingleMu22_EMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF_NEG_pRECO", &L1_SingleMu22_EMTF_NEG_pRECO, &b_L1_SingleMu22_EMTF_NEG_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_EMTF_POS_pRECO", &L1_SingleMu22_EMTF_POS_pRECO, &b_L1_SingleMu22_EMTF_POS_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF_pRECO", &L1_SingleMu22_OMTF_pRECO, &b_L1_SingleMu22_OMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF_NEG_pRECO", &L1_SingleMu22_OMTF_NEG_pRECO, &b_L1_SingleMu22_OMTF_NEG_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_OMTF_POS_pRECO", &L1_SingleMu22_OMTF_POS_pRECO, &b_L1_SingleMu22_OMTF_POS_pRECO);
   fChain->SetBranchAddress("L1_SingleMu22_OQ_pRECO", &L1_SingleMu22_OQ_pRECO, &b_L1_SingleMu22_OQ_pRECO);
   fChain->SetBranchAddress("L1_SingleMu25_pRECO", &L1_SingleMu25_pRECO, &b_L1_SingleMu25_pRECO);
   fChain->SetBranchAddress("L1_SingleMu3_pRECO", &L1_SingleMu3_pRECO, &b_L1_SingleMu3_pRECO);
   fChain->SetBranchAddress("L1_SingleMu5_pRECO", &L1_SingleMu5_pRECO, &b_L1_SingleMu5_pRECO);
   fChain->SetBranchAddress("L1_SingleMu5_SQ14_BMTF_pRECO", &L1_SingleMu5_SQ14_BMTF_pRECO, &b_L1_SingleMu5_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu6_SQ14_BMTF_pRECO", &L1_SingleMu6_SQ14_BMTF_pRECO, &b_L1_SingleMu6_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu7_pRECO", &L1_SingleMu7_pRECO, &b_L1_SingleMu7_pRECO);
   fChain->SetBranchAddress("L1_SingleMu7_DQ_pRECO", &L1_SingleMu7_DQ_pRECO, &b_L1_SingleMu7_DQ_pRECO);
   fChain->SetBranchAddress("L1_SingleMu7_SQ14_BMTF_pRECO", &L1_SingleMu7_SQ14_BMTF_pRECO, &b_L1_SingleMu7_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu8_SQ14_BMTF_pRECO", &L1_SingleMu8_SQ14_BMTF_pRECO, &b_L1_SingleMu8_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMu9_SQ14_BMTF_pRECO", &L1_SingleMu9_SQ14_BMTF_pRECO, &b_L1_SingleMu9_SQ14_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMuCosmics_pRECO", &L1_SingleMuCosmics_pRECO, &b_L1_SingleMuCosmics_pRECO);
   fChain->SetBranchAddress("L1_SingleMuCosmics_BMTF_pRECO", &L1_SingleMuCosmics_BMTF_pRECO, &b_L1_SingleMuCosmics_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMuCosmics_EMTF_pRECO", &L1_SingleMuCosmics_EMTF_pRECO, &b_L1_SingleMuCosmics_EMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMuCosmics_OMTF_pRECO", &L1_SingleMuCosmics_OMTF_pRECO, &b_L1_SingleMuCosmics_OMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMuOpen_pRECO", &L1_SingleMuOpen_pRECO, &b_L1_SingleMuOpen_pRECO);
   fChain->SetBranchAddress("L1_SingleMuOpen_BMTF_pRECO", &L1_SingleMuOpen_BMTF_pRECO, &b_L1_SingleMuOpen_BMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMuOpen_EMTF_pRECO", &L1_SingleMuOpen_EMTF_pRECO, &b_L1_SingleMuOpen_EMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR_pRECO", &L1_SingleMuOpen_NotBptxOR_pRECO, &b_L1_SingleMuOpen_NotBptxOR_pRECO);
   fChain->SetBranchAddress("L1_SingleMuOpen_OMTF_pRECO", &L1_SingleMuOpen_OMTF_pRECO, &b_L1_SingleMuOpen_OMTF_pRECO);
   fChain->SetBranchAddress("L1_SingleMuOpen_er1p1_NotBptxOR_3BX_pRECO", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX_pRECO, &b_L1_SingleMuOpen_er1p1_NotBptxOR_3BX_pRECO);
   fChain->SetBranchAddress("L1_SingleMuOpen_er1p4_NotBptxOR_3BX_pRECO", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX_pRECO, &b_L1_SingleMuOpen_er1p4_NotBptxOR_3BX_pRECO);
   fChain->SetBranchAddress("L1_SingleMuShower_Nominal_pRECO", &L1_SingleMuShower_Nominal_pRECO, &b_L1_SingleMuShower_Nominal_pRECO);
   fChain->SetBranchAddress("L1_SingleMuShower_Tight_pRECO", &L1_SingleMuShower_Tight_pRECO, &b_L1_SingleMuShower_Tight_pRECO);
   fChain->SetBranchAddress("L1_SingleTau120er2p1_pRECO", &L1_SingleTau120er2p1_pRECO, &b_L1_SingleTau120er2p1_pRECO);
   fChain->SetBranchAddress("L1_SingleTau130er2p1_pRECO", &L1_SingleTau130er2p1_pRECO, &b_L1_SingleTau130er2p1_pRECO);
   fChain->SetBranchAddress("L1_TOTEM_1_pRECO", &L1_TOTEM_1_pRECO, &b_L1_TOTEM_1_pRECO);
   fChain->SetBranchAddress("L1_TOTEM_2_pRECO", &L1_TOTEM_2_pRECO, &b_L1_TOTEM_2_pRECO);
   fChain->SetBranchAddress("L1_TOTEM_3_pRECO", &L1_TOTEM_3_pRECO, &b_L1_TOTEM_3_pRECO);
   fChain->SetBranchAddress("L1_TOTEM_4_pRECO", &L1_TOTEM_4_pRECO, &b_L1_TOTEM_4_pRECO);
   fChain->SetBranchAddress("L1_TripleEG_18_17_8_er2p5_pRECO", &L1_TripleEG_18_17_8_er2p5_pRECO, &b_L1_TripleEG_18_17_8_er2p5_pRECO);
   fChain->SetBranchAddress("L1_TripleEG_18_18_12_er2p5_pRECO", &L1_TripleEG_18_18_12_er2p5_pRECO, &b_L1_TripleEG_18_18_12_er2p5_pRECO);
   fChain->SetBranchAddress("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5_pRECO", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5_pRECO, &b_L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5_pRECO);
   fChain->SetBranchAddress("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5_pRECO", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5_pRECO, &b_L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5_pRECO);
   fChain->SetBranchAddress("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5_pRECO", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5_pRECO, &b_L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5_pRECO);
   fChain->SetBranchAddress("L1_TripleMu0_pRECO", &L1_TripleMu0_pRECO, &b_L1_TripleMu0_pRECO);
   fChain->SetBranchAddress("L1_TripleMu0_OQ_pRECO", &L1_TripleMu0_OQ_pRECO, &b_L1_TripleMu0_OQ_pRECO);
   fChain->SetBranchAddress("L1_TripleMu0_SQ_pRECO", &L1_TripleMu0_SQ_pRECO, &b_L1_TripleMu0_SQ_pRECO);
   fChain->SetBranchAddress("L1_TripleMu3_pRECO", &L1_TripleMu3_pRECO, &b_L1_TripleMu3_pRECO);
   fChain->SetBranchAddress("L1_TripleMu3_SQ_pRECO", &L1_TripleMu3_SQ_pRECO, &b_L1_TripleMu3_SQ_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_3SQ_2p5SQ_0_pRECO", &L1_TripleMu_3SQ_2p5SQ_0_pRECO, &b_L1_TripleMu_3SQ_2p5SQ_0_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12_pRECO", &L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12_pRECO, &b_L1_TripleMu_3SQ_2p5SQ_0_Mass_Max12_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12_pRECO", &L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12_pRECO, &b_L1_TripleMu_3SQ_2p5SQ_0_OS_Mass_Max12_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12_pRECO", &L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12_pRECO, &b_L1_TripleMu_4SQ_2p5SQ_0_OS_Mass_Max12_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_pRECO", &L1_TripleMu_5SQ_3SQ_0OQ_pRECO, &b_L1_TripleMu_5SQ_3SQ_0OQ_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO, &b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO, &b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5_3_3_pRECO", &L1_TripleMu_5_3_3_pRECO, &b_L1_TripleMu_5_3_3_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5_3_3_SQ_pRECO", &L1_TripleMu_5_3_3_SQ_pRECO, &b_L1_TripleMu_5_3_3_SQ_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5_pRECO", &L1_TripleMu_5_3p5_2p5_pRECO, &b_L1_TripleMu_5_3p5_2p5_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO, &b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO, &b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17_pRECO);
   fChain->SetBranchAddress("L1_TripleMu_5_5_3_pRECO", &L1_TripleMu_5_5_3_pRECO, &b_L1_TripleMu_5_5_3_pRECO);
   fChain->SetBranchAddress("L1_TwoMuShower_Loose_pRECO", &L1_TwoMuShower_Loose_pRECO, &b_L1_TwoMuShower_Loose_pRECO);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxMinus_pRECO", &L1_UnpairedBunchBptxMinus_pRECO, &b_L1_UnpairedBunchBptxMinus_pRECO);
   fChain->SetBranchAddress("L1_UnpairedBunchBptxPlus_pRECO", &L1_UnpairedBunchBptxPlus_pRECO, &b_L1_UnpairedBunchBptxPlus_pRECO);
   fChain->SetBranchAddress("L1_ZeroBias_pRECO", &L1_ZeroBias_pRECO, &b_L1_ZeroBias_pRECO);
   fChain->SetBranchAddress("L1_ZeroBias_copy_pRECO", &L1_ZeroBias_copy_pRECO, &b_L1_ZeroBias_copy_pRECO);
   fChain->SetBranchAddress("L1_UnprefireableEvent_TriggerRules_pRECO", &L1_UnprefireableEvent_TriggerRules_pRECO, &b_L1_UnprefireableEvent_TriggerRules_pRECO);
   fChain->SetBranchAddress("L1_UnprefireableEvent_FirstBxInTrain_pRECO", &L1_UnprefireableEvent_FirstBxInTrain_pRECO, &b_L1_UnprefireableEvent_FirstBxInTrain_pRECO);
   fChain->SetBranchAddress("L1_FinalOR_BXmin1_pRECO", &L1_FinalOR_BXmin1_pRECO, &b_L1_FinalOR_BXmin1_pRECO);
   fChain->SetBranchAddress("L1_FinalOR_BXmin2_pRECO", &L1_FinalOR_BXmin2_pRECO, &b_L1_FinalOR_BXmin2_pRECO);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("DST_PFScouting_DoubleMuon", &DST_PFScouting_DoubleMuon, &b_DST_PFScouting_DoubleMuon);
   fChain->SetBranchAddress("DST_PFScouting_DoubleEG", &DST_PFScouting_DoubleEG, &b_DST_PFScouting_DoubleEG);
   fChain->SetBranchAddress("DST_PFScouting_JetHT", &DST_PFScouting_JetHT, &b_DST_PFScouting_JetHT);
   fChain->SetBranchAddress("DST_PFScouting_DatasetMuon", &DST_PFScouting_DatasetMuon, &b_DST_PFScouting_DatasetMuon);
   fChain->SetBranchAddress("DST_PFScouting_AXOVLoose", &DST_PFScouting_AXOVLoose, &b_DST_PFScouting_AXOVLoose);
   fChain->SetBranchAddress("DST_PFScouting_AXOLoose", &DST_PFScouting_AXOLoose, &b_DST_PFScouting_AXOLoose);
   fChain->SetBranchAddress("DST_PFScouting_AXONominal", &DST_PFScouting_AXONominal, &b_DST_PFScouting_AXONominal);
   fChain->SetBranchAddress("DST_PFScouting_AXOTight", &DST_PFScouting_AXOTight, &b_DST_PFScouting_AXOTight);
   fChain->SetBranchAddress("DST_PFScouting_AXOVTight", &DST_PFScouting_AXOVTight, &b_DST_PFScouting_AXOVTight);
   fChain->SetBranchAddress("DST_PFScouting_CICADAVLoose", &DST_PFScouting_CICADAVLoose, &b_DST_PFScouting_CICADAVLoose);
   fChain->SetBranchAddress("DST_PFScouting_CICADALoose", &DST_PFScouting_CICADALoose, &b_DST_PFScouting_CICADALoose);
   fChain->SetBranchAddress("DST_PFScouting_CICADAMedium", &DST_PFScouting_CICADAMedium, &b_DST_PFScouting_CICADAMedium);
   fChain->SetBranchAddress("DST_PFScouting_CICADATight", &DST_PFScouting_CICADATight, &b_DST_PFScouting_CICADATight);
   fChain->SetBranchAddress("DST_PFScouting_CICADAVTight", &DST_PFScouting_CICADAVTight, &b_DST_PFScouting_CICADAVTight);
   fChain->SetBranchAddress("DST_PFScouting_SingleMuon", &DST_PFScouting_SingleMuon, &b_DST_PFScouting_SingleMuon);
   fChain->SetBranchAddress("DST_PFScouting_SinglePhotonEB", &DST_PFScouting_SinglePhotonEB, &b_DST_PFScouting_SinglePhotonEB);
   fChain->SetBranchAddress("DST_PFScouting_ZeroBias", &DST_PFScouting_ZeroBias, &b_DST_PFScouting_ZeroBias);
   fChain->SetBranchAddress("HLT_EphemeralPhysics", &HLT_EphemeralPhysics, &b_HLT_EphemeralPhysics);
   fChain->SetBranchAddress("HLT_TestPhysics", &HLT_TestPhysics, &b_HLT_TestPhysics);
   fChain->SetBranchAddress("HLT_EphemeralZeroBias", &HLT_EphemeralZeroBias, &b_HLT_EphemeralZeroBias);
   fChain->SetBranchAddress("HLT_SpecialZeroBias", &HLT_SpecialZeroBias, &b_HLT_SpecialZeroBias);
   fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   fChain->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLT_ZeroBias_Alignment", &HLT_ZeroBias_Alignment, &b_HLT_ZeroBias_Alignment);
   fChain->SetBranchAddress("HLT_ZeroBias_Beamspot", &HLT_ZeroBias_Beamspot, &b_HLT_ZeroBias_Beamspot);
   fChain->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstBXAfterTrain", &HLT_ZeroBias_FirstBXAfterTrain, &b_HLT_ZeroBias_FirstBXAfterTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionInTrain", &HLT_ZeroBias_FirstCollisionInTrain, &b_HLT_ZeroBias_FirstCollisionInTrain);
   fChain->SetBranchAddress("HLT_ZeroBias_LastCollisionInTrain", &HLT_ZeroBias_LastCollisionInTrain, &b_HLT_ZeroBias_LastCollisionInTrain);
   fChain->SetBranchAddress("HLT_HT300_Beamspot", &HLT_HT300_Beamspot, &b_HLT_HT300_Beamspot);
   fChain->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
   fChain->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);
   fChain->SetBranchAddress("HLT_PFJet40_GPUvsCPU", &HLT_PFJet40_GPUvsCPU, &b_HLT_PFJet40_GPUvsCPU);
   fChain->SetBranchAddress("HLT_AK8PFJet380_SoftDropMass30", &HLT_AK8PFJet380_SoftDropMass30, &b_HLT_AK8PFJet380_SoftDropMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet400_SoftDropMass30", &HLT_AK8PFJet400_SoftDropMass30, &b_HLT_AK8PFJet400_SoftDropMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet425_SoftDropMass30", &HLT_AK8PFJet425_SoftDropMass30, &b_HLT_AK8PFJet425_SoftDropMass30);
   fChain->SetBranchAddress("HLT_AK8PFJet450_SoftDropMass30", &HLT_AK8PFJet450_SoftDropMass30, &b_HLT_AK8PFJet450_SoftDropMass30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_250_SoftDropMass40", &HLT_AK8DiPFJet250_250_SoftDropMass40, &b_HLT_AK8DiPFJet250_250_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_250_SoftDropMass50", &HLT_AK8DiPFJet250_250_SoftDropMass50, &b_HLT_AK8DiPFJet250_250_SoftDropMass50);
   fChain->SetBranchAddress("HLT_AK8DiPFJet260_260_SoftDropMass30", &HLT_AK8DiPFJet260_260_SoftDropMass30, &b_HLT_AK8DiPFJet260_260_SoftDropMass30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet260_260_SoftDropMass40", &HLT_AK8DiPFJet260_260_SoftDropMass40, &b_HLT_AK8DiPFJet260_260_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8DiPFJet270_270_SoftDropMass30", &HLT_AK8DiPFJet270_270_SoftDropMass30, &b_HLT_AK8DiPFJet270_270_SoftDropMass30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_280_SoftDropMass30", &HLT_AK8DiPFJet280_280_SoftDropMass30, &b_HLT_AK8DiPFJet280_280_SoftDropMass30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet290_290_SoftDropMass30", &HLT_AK8DiPFJet290_290_SoftDropMass30, &b_HLT_AK8DiPFJet290_290_SoftDropMass30);
   fChain->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
   fChain->SetBranchAddress("HLT_CaloJet550_NoJetID", &HLT_CaloJet550_NoJetID, &b_HLT_CaloJet550_NoJetID);
   fChain->SetBranchAddress("HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL", &HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL, &b_HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon", &HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon, &b_HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon);
   fChain->SetBranchAddress("HLT_Trimuon5_3p5_2_Upsilon_Muon", &HLT_Trimuon5_3p5_2_Upsilon_Muon, &b_HLT_Trimuon5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon", &HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon, &b_HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon);
   fChain->SetBranchAddress("HLT_DoubleEle25_CaloIdL_MW", &HLT_DoubleEle25_CaloIdL_MW, &b_HLT_DoubleEle25_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle27_CaloIdL_MW", &HLT_DoubleEle27_CaloIdL_MW, &b_HLT_DoubleEle27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle24_eta2p1_WPTight_Gsf", &HLT_DoubleEle24_eta2p1_WPTight_Gsf, &b_HLT_DoubleEle24_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_MW", &HLT_Mu27_Ele37_CaloIdL_MW, &b_HLT_Mu27_Ele37_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_MW", &HLT_Mu37_Ele27_CaloIdL_MW, &b_HLT_Mu37_Ele27_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_Mu37_TkMu27", &HLT_Mu37_TkMu27, &b_HLT_Mu37_TkMu27);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi", &HLT_DoubleMu4_3_Jpsi, &b_HLT_DoubleMu4_3_Jpsi);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_LowMass", &HLT_DoubleMu4_3_LowMass, &b_HLT_DoubleMu4_3_LowMass);
   fChain->SetBranchAddress("HLT_DoubleMu4_LowMass_Displaced", &HLT_DoubleMu4_LowMass_Displaced, &b_HLT_DoubleMu4_LowMass_Displaced);
   fChain->SetBranchAddress("HLT_Mu0_L1DoubleMu", &HLT_Mu0_L1DoubleMu, &b_HLT_Mu0_L1DoubleMu);
   fChain->SetBranchAddress("HLT_Mu4_L1DoubleMu", &HLT_Mu4_L1DoubleMu, &b_HLT_Mu4_L1DoubleMu);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_LowPt", &HLT_DoubleMu2_Jpsi_LowPt, &b_HLT_DoubleMu2_Jpsi_LowPt);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Photon4_BsToMMG", &HLT_DoubleMu4_3_Photon4_BsToMMG, &b_HLT_DoubleMu4_3_Photon4_BsToMMG);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG", &HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG, &b_HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu, &b_HLT_DoubleMu3_Trk_Tau3mu);
   fChain->SetBranchAddress("HLT_DoubleMu3_TkMu_DsTau3Mu", &HLT_DoubleMu3_TkMu_DsTau3Mu, &b_HLT_DoubleMu3_TkMu_DsTau3Mu);
   fChain->SetBranchAddress("HLT_DoubleMu4_Mass3p8_DZ_PFHT350", &HLT_DoubleMu4_Mass3p8_DZ_PFHT350, &b_HLT_DoubleMu4_Mass3p8_DZ_PFHT350);
   fChain->SetBranchAddress("HLT_DoubleMu4_MuMuTrk_Displaced", &HLT_DoubleMu4_MuMuTrk_Displaced, &b_HLT_DoubleMu4_MuMuTrk_Displaced);
   fChain->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40, &b_HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu3_L1SingleMu5orSingleMu7", &HLT_Mu3_L1SingleMu5orSingleMu7, &b_HLT_Mu3_L1SingleMu5orSingleMu7);
   fChain->SetBranchAddress("HLT_Mu0_Barrel", &HLT_Mu0_Barrel, &b_HLT_Mu0_Barrel);
   fChain->SetBranchAddress("HLT_Mu0_Barrel_L1HP6", &HLT_Mu0_Barrel_L1HP6, &b_HLT_Mu0_Barrel_L1HP6);
   fChain->SetBranchAddress("HLT_Mu0_Barrel_L1HP7", &HLT_Mu0_Barrel_L1HP7, &b_HLT_Mu0_Barrel_L1HP7);
   fChain->SetBranchAddress("HLT_Mu0_Barrel_L1HP8", &HLT_Mu0_Barrel_L1HP8, &b_HLT_Mu0_Barrel_L1HP8);
   fChain->SetBranchAddress("HLT_Mu0_Barrel_L1HP9", &HLT_Mu0_Barrel_L1HP9, &b_HLT_Mu0_Barrel_L1HP9);
   fChain->SetBranchAddress("HLT_Mu0_Barrel_L1HP10", &HLT_Mu0_Barrel_L1HP10, &b_HLT_Mu0_Barrel_L1HP10);
   fChain->SetBranchAddress("HLT_Mu0_Barrel_L1HP11", &HLT_Mu0_Barrel_L1HP11, &b_HLT_Mu0_Barrel_L1HP11);
   fChain->SetBranchAddress("HLT_Mu0_Barrel_L1HP6_IP6", &HLT_Mu0_Barrel_L1HP6_IP6, &b_HLT_Mu0_Barrel_L1HP6_IP6);
   fChain->SetBranchAddress("HLT_Mu6_Barrel_L1HP7_IP6", &HLT_Mu6_Barrel_L1HP7_IP6, &b_HLT_Mu6_Barrel_L1HP7_IP6);
   fChain->SetBranchAddress("HLT_Mu7_Barrel_L1HP8_IP6", &HLT_Mu7_Barrel_L1HP8_IP6, &b_HLT_Mu7_Barrel_L1HP8_IP6);
   fChain->SetBranchAddress("HLT_Mu8_Barrel_L1HP9_IP6", &HLT_Mu8_Barrel_L1HP9_IP6, &b_HLT_Mu8_Barrel_L1HP9_IP6);
   fChain->SetBranchAddress("HLT_Mu9_Barrel_L1HP10_IP6", &HLT_Mu9_Barrel_L1HP10_IP6, &b_HLT_Mu9_Barrel_L1HP10_IP6);
   fChain->SetBranchAddress("HLT_Mu10_Barrel_L1HP11_IP6", &HLT_Mu10_Barrel_L1HP11_IP6, &b_HLT_Mu10_Barrel_L1HP11_IP6);
   fChain->SetBranchAddress("HLT_DoublePhoton33_CaloIdL", &HLT_DoublePhoton33_CaloIdL, &b_HLT_DoublePhoton33_CaloIdL);
   fChain->SetBranchAddress("HLT_DoublePhoton70", &HLT_DoublePhoton70, &b_HLT_DoublePhoton70);
   fChain->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
   fChain->SetBranchAddress("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG, &b_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG);
   fChain->SetBranchAddress("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf, &b_HLT_Ele30_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf, &b_HLT_Ele35_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf, &b_HLT_Ele38_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf, &b_HLT_Ele40_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX", &HLT_UncorrectedJetE30_NoBPTX, &b_HLT_UncorrectedJetE30_NoBPTX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE30_NoBPTX3BX", &HLT_UncorrectedJetE30_NoBPTX3BX, &b_HLT_UncorrectedJetE30_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE60_NoBPTX3BX", &HLT_UncorrectedJetE60_NoBPTX3BX, &b_HLT_UncorrectedJetE60_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_UncorrectedJetE70_NoBPTX3BX", &HLT_UncorrectedJetE70_NoBPTX3BX, &b_HLT_UncorrectedJetE70_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L1SingleMuCosmics", &HLT_L1SingleMuCosmics, &b_HLT_L1SingleMuCosmics);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_HLT_L2Mu10_NoVertex_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX, &b_HLT_L2Mu10_NoVertex_NoBPTX);
   fChain->SetBranchAddress("HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha", &HLT_L2Mu23NoVtx_2Cha, &b_HLT_L2Mu23NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_L2Mu23NoVtx_2Cha_CosmicSeed", &HLT_L2Mu23NoVtx_2Cha_CosmicSeed, &b_HLT_L2Mu23NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4, &b_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4);
   fChain->SetBranchAddress("HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4, &b_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4);
   fChain->SetBranchAddress("HLT_DoubleL2Mu50", &HLT_DoubleL2Mu50, &b_HLT_DoubleL2Mu50);
   fChain->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed, &b_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed, &b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4, &b_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4);
   fChain->SetBranchAddress("HLT_DoubleL2Mu23NoVtx_2Cha", &HLT_DoubleL2Mu23NoVtx_2Cha, &b_HLT_DoubleL2Mu23NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha", &HLT_DoubleL2Mu25NoVtx_2Cha, &b_HLT_DoubleL2Mu25NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4", &HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4, &b_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet30);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_CaloJet30);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8PFJet30);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_AK8CaloJet30);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8", &HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8, &b_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Psi", &HLT_Mu30_TkMu0_Psi, &b_HLT_Mu30_TkMu0_Psi);
   fChain->SetBranchAddress("HLT_Mu30_TkMu0_Upsilon", &HLT_Mu30_TkMu0_Upsilon, &b_HLT_Mu30_TkMu0_Upsilon);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_Phi", &HLT_Mu25_TkMu0_Phi, &b_HLT_Mu25_TkMu0_Phi);
   fChain->SetBranchAddress("HLT_Mu15", &HLT_Mu15, &b_HLT_Mu15);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_CascadeMu100", &HLT_CascadeMu100, &b_HLT_CascadeMu100);
   fChain->SetBranchAddress("HLT_HighPtTkMu100", &HLT_HighPtTkMu100, &b_HLT_HighPtTkMu100);
   fChain->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
   fChain->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
   fChain->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
   fChain->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
   fChain->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
   fChain->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
   fChain->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
   fChain->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
   fChain->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve260_HFJEC", &HLT_DiPFJetAve260_HFJEC, &b_HLT_DiPFJetAve260_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP", &HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP, &b_HLT_DiPFJetAve180_PPSMatch_Xi0p3_QuadJet_Max2ProtPerRP);
   fChain->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
   fChain->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
   fChain->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
   fChain->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
   fChain->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
   fChain->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
   fChain->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
   fChain->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
   fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   fChain->SetBranchAddress("HLT_AK8PFJet550", &HLT_AK8PFJet550, &b_HLT_AK8PFJet550);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet40_ZeroBiasCopy", &HLT_PFJet40_ZeroBiasCopy, &b_HLT_PFJet40_ZeroBiasCopy);
   fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet80_L1SingleJet60", &HLT_PFJet80_L1SingleJet60, &b_HLT_PFJet80_L1SingleJet60);
   fChain->SetBranchAddress("HLT_PFJet110", &HLT_PFJet110, &b_HLT_PFJet110);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_PFJet550", &HLT_PFJet550, &b_HLT_PFJet550);
   fChain->SetBranchAddress("HLT_PFJetFwd40", &HLT_PFJetFwd40, &b_HLT_PFJetFwd40);
   fChain->SetBranchAddress("HLT_PFJetFwd60", &HLT_PFJetFwd60, &b_HLT_PFJetFwd60);
   fChain->SetBranchAddress("HLT_PFJetFwd80", &HLT_PFJetFwd80, &b_HLT_PFJetFwd80);
   fChain->SetBranchAddress("HLT_PFJetFwd140", &HLT_PFJetFwd140, &b_HLT_PFJetFwd140);
   fChain->SetBranchAddress("HLT_PFJetFwd200", &HLT_PFJetFwd200, &b_HLT_PFJetFwd200);
   fChain->SetBranchAddress("HLT_PFJetFwd260", &HLT_PFJetFwd260, &b_HLT_PFJetFwd260);
   fChain->SetBranchAddress("HLT_PFJetFwd320", &HLT_PFJetFwd320, &b_HLT_PFJetFwd320);
   fChain->SetBranchAddress("HLT_PFJetFwd400", &HLT_PFJetFwd400, &b_HLT_PFJetFwd400);
   fChain->SetBranchAddress("HLT_PFJetFwd450", &HLT_PFJetFwd450, &b_HLT_PFJetFwd450);
   fChain->SetBranchAddress("HLT_PFJetFwd500", &HLT_PFJetFwd500, &b_HLT_PFJetFwd500);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd40", &HLT_AK8PFJetFwd40, &b_HLT_AK8PFJetFwd40);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd60", &HLT_AK8PFJetFwd60, &b_HLT_AK8PFJetFwd60);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd80", &HLT_AK8PFJetFwd80, &b_HLT_AK8PFJetFwd80);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd140", &HLT_AK8PFJetFwd140, &b_HLT_AK8PFJetFwd140);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd200", &HLT_AK8PFJetFwd200, &b_HLT_AK8PFJetFwd200);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd260", &HLT_AK8PFJetFwd260, &b_HLT_AK8PFJetFwd260);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd320", &HLT_AK8PFJetFwd320, &b_HLT_AK8PFJetFwd320);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd400", &HLT_AK8PFJetFwd400, &b_HLT_AK8PFJetFwd400);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd450", &HLT_AK8PFJetFwd450, &b_HLT_AK8PFJetFwd450);
   fChain->SetBranchAddress("HLT_AK8PFJetFwd500", &HLT_AK8PFJetFwd500, &b_HLT_AK8PFJetFwd500);
   fChain->SetBranchAddress("HLT_PFHT180", &HLT_PFHT180, &b_HLT_PFHT180);
   fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   fChain->SetBranchAddress("HLT_PFHT370", &HLT_PFHT370, &b_HLT_PFHT370);
   fChain->SetBranchAddress("HLT_PFHT430", &HLT_PFHT430, &b_HLT_PFHT430);
   fChain->SetBranchAddress("HLT_PFHT510", &HLT_PFHT510, &b_HLT_PFHT510);
   fChain->SetBranchAddress("HLT_PFHT590", &HLT_PFHT590, &b_HLT_PFHT590);
   fChain->SetBranchAddress("HLT_PFHT680", &HLT_PFHT680, &b_HLT_PFHT680);
   fChain->SetBranchAddress("HLT_PFHT780", &HLT_PFHT780, &b_HLT_PFHT780);
   fChain->SetBranchAddress("HLT_PFHT890", &HLT_PFHT890, &b_HLT_PFHT890);
   fChain->SetBranchAddress("HLT_PFHT1050", &HLT_PFHT1050, &b_HLT_PFHT1050);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight", &HLT_PFHT500_PFMET100_PFMHT100_IDTight, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight", &HLT_PFHT500_PFMET110_PFMHT110_IDTight, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight", &HLT_PFHT700_PFMET85_PFMHT85_IDTight, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight", &HLT_PFHT800_PFMET75_PFMHT75_IDTight, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight, &b_HLT_PFMET130_PFMHT130_IDTight);
   fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight, &b_HLT_PFMET140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
   fChain->SetBranchAddress("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight, &b_HLT_PFMETTypeOne140_PFMHT140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF);
   fChain->SetBranchAddress("HLT_L1ETMHadSeeds", &HLT_L1ETMHadSeeds, &b_HLT_L1ETMHadSeeds);
   fChain->SetBranchAddress("HLT_CaloMHT90", &HLT_CaloMHT90, &b_HLT_CaloMHT90);
   fChain->SetBranchAddress("HLT_CaloMET90_NotCleaned", &HLT_CaloMET90_NotCleaned, &b_HLT_CaloMET90_NotCleaned);
   fChain->SetBranchAddress("HLT_CaloMET350_NotCleaned", &HLT_CaloMET350_NotCleaned, &b_HLT_CaloMET350_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned, &b_HLT_PFMET200_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET250_NotCleaned", &HLT_PFMET250_NotCleaned, &b_HLT_PFMET250_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET300_NotCleaned", &HLT_PFMET300_NotCleaned, &b_HLT_PFMET300_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET200_BeamHaloCleaned", &HLT_PFMET200_BeamHaloCleaned, &b_HLT_PFMET200_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMETTypeOne200_BeamHaloCleaned", &HLT_PFMETTypeOne200_BeamHaloCleaned, &b_HLT_PFMETTypeOne200_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_MET105_IsoTrk50", &HLT_MET105_IsoTrk50, &b_HLT_MET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_MET120_IsoTrk50", &HLT_MET120_IsoTrk50, &b_HLT_MET120_IsoTrk50);
   fChain->SetBranchAddress("HLT_Mu12eta2p3", &HLT_Mu12eta2p3, &b_HLT_Mu12eta2p3);
   fChain->SetBranchAddress("HLT_Mu12eta2p3_PFJet40", &HLT_Mu12eta2p3_PFJet40, &b_HLT_Mu12eta2p3_PFJet40);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40_PNetBTag_0p11", &HLT_Mu12_DoublePFJets40_PNetBTag_0p11, &b_HLT_Mu12_DoublePFJets40_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets100_PNetBTag_0p11", &HLT_Mu12_DoublePFJets100_PNetBTag_0p11, &b_HLT_Mu12_DoublePFJets100_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets200_PNetBTag_0p11", &HLT_Mu12_DoublePFJets200_PNetBTag_0p11, &b_HLT_Mu12_DoublePFJets200_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets350_PNetBTag_0p11", &HLT_Mu12_DoublePFJets350_PNetBTag_0p11, &b_HLT_Mu12_DoublePFJets350_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11", &HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11, &b_HLT_Mu12_DoublePFJets40MaxDeta1p6_PNet2BTag_0p11);
   fChain->SetBranchAddress("HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11", &HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11, &b_HLT_Mu12_DoublePFJets54MaxDeta1p6_PNet2BTag_0p11);
   fChain->SetBranchAddress("HLT_DoublePFJets40_PNetBTag_0p11", &HLT_DoublePFJets40_PNetBTag_0p11, &b_HLT_DoublePFJets40_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_DoublePFJets100_PNetBTag_0p11", &HLT_DoublePFJets100_PNetBTag_0p11, &b_HLT_DoublePFJets100_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_DoublePFJets200_PNetBTag_0p11", &HLT_DoublePFJets200_PNetBTag_0p11, &b_HLT_DoublePFJets200_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_DoublePFJets350_PNetBTag_0p11", &HLT_DoublePFJets350_PNetBTag_0p11, &b_HLT_DoublePFJets350_PNetBTag_0p11);
   fChain->SetBranchAddress("HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11", &HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11, &b_HLT_DoublePFJets116MaxDeta1p6_PNet2BTag_0p11);
   fChain->SetBranchAddress("HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11", &HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11, &b_HLT_DoublePFJets128MaxDeta1p6_PNet2BTag_0p11);
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu19_TrkIsoVVL", &HLT_Mu19_TrkIsoVVL, &b_HLT_Mu19_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet20_Mu5", &HLT_BTagMu_AK4DiJet20_Mu5, &b_HLT_BTagMu_AK4DiJet20_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet40_Mu5", &HLT_BTagMu_AK4DiJet40_Mu5, &b_HLT_BTagMu_AK4DiJet40_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet70_Mu5", &HLT_BTagMu_AK4DiJet70_Mu5, &b_HLT_BTagMu_AK4DiJet70_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet110_Mu5", &HLT_BTagMu_AK4DiJet110_Mu5, &b_HLT_BTagMu_AK4DiJet110_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4DiJet170_Mu5", &HLT_BTagMu_AK4DiJet170_Mu5, &b_HLT_BTagMu_AK4DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK4Jet300_Mu5", &HLT_BTagMu_AK4Jet300_Mu5, &b_HLT_BTagMu_AK4Jet300_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8DiJet170_Mu5", &HLT_BTagMu_AK8DiJet170_Mu5, &b_HLT_BTagMu_AK8DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8Jet170_DoubleMu5", &HLT_BTagMu_AK8Jet170_DoubleMu5, &b_HLT_BTagMu_AK8Jet170_DoubleMu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5, &b_HLT_BTagMu_AK8Jet300_Mu5);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53", &HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53, &b_HLT_Mu12_IsoVVL_PFHT150_PNetBTag0p53);
   fChain->SetBranchAddress("HLT_Photon33", &HLT_Photon33, &b_HLT_Photon33);
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_Photon200", &HLT_Photon200, &b_HLT_Photon200);
   fChain->SetBranchAddress("HLT_Photon45EB", &HLT_Photon45EB, &b_HLT_Photon45EB);
   fChain->SetBranchAddress("HLT_Photon40EB", &HLT_Photon40EB, &b_HLT_Photon40EB);
   fChain->SetBranchAddress("HLT_Photon50EB", &HLT_Photon50EB, &b_HLT_Photon50EB);
   fChain->SetBranchAddress("HLT_Photon30EB_TightID_TightIso", &HLT_Photon30EB_TightID_TightIso, &b_HLT_Photon30EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon40EB_TightID_TightIso", &HLT_Photon40EB_TightID_TightIso, &b_HLT_Photon40EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon45EB_TightID_TightIso", &HLT_Photon45EB_TightID_TightIso, &b_HLT_Photon45EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso_PFJet30", &HLT_Photon50EB_TightID_TightIso_PFJet30, &b_HLT_Photon50EB_TightID_TightIso_PFJet30);
   fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso_CaloJet30", &HLT_Photon50EB_TightID_TightIso_CaloJet30, &b_HLT_Photon50EB_TightID_TightIso_CaloJet30);
   fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso_AK8PFJet30", &HLT_Photon50EB_TightID_TightIso_AK8PFJet30, &b_HLT_Photon50EB_TightID_TightIso_AK8PFJet30);
   fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso_AK8CaloJet30", &HLT_Photon50EB_TightID_TightIso_AK8CaloJet30, &b_HLT_Photon50EB_TightID_TightIso_AK8CaloJet30);
   fChain->SetBranchAddress("HLT_Photon50EB_TightID_TightIso", &HLT_Photon50EB_TightID_TightIso, &b_HLT_Photon50EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon55EB_TightID_TightIso", &HLT_Photon55EB_TightID_TightIso, &b_HLT_Photon55EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon75EB_TightID_TightIso", &HLT_Photon75EB_TightID_TightIso, &b_HLT_Photon75EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon90EB_TightID_TightIso", &HLT_Photon90EB_TightID_TightIso, &b_HLT_Photon90EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso, &b_HLT_Photon110EB_TightID_TightIso);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso_PFJet30", &HLT_Photon110EB_TightID_TightIso_PFJet30, &b_HLT_Photon110EB_TightID_TightIso_PFJet30);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso_CaloJet30", &HLT_Photon110EB_TightID_TightIso_CaloJet30, &b_HLT_Photon110EB_TightID_TightIso_CaloJet30);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso_AK8PFJet30", &HLT_Photon110EB_TightID_TightIso_AK8PFJet30, &b_HLT_Photon110EB_TightID_TightIso_AK8PFJet30);
   fChain->SetBranchAddress("HLT_Photon110EB_TightID_TightIso_AK8CaloJet30", &HLT_Photon110EB_TightID_TightIso_AK8CaloJet30, &b_HLT_Photon110EB_TightID_TightIso_AK8CaloJet30);
   fChain->SetBranchAddress("HLT_Photon100EBHE10", &HLT_Photon100EBHE10, &b_HLT_Photon100EBHE10);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95);
   fChain->SetBranchAddress("HLT_DiphotonMVA14p25_Mass90", &HLT_DiphotonMVA14p25_Mass90, &b_HLT_DiphotonMVA14p25_Mass90);
   fChain->SetBranchAddress("HLT_DiphotonMVA14p25_Tight_Mass90", &HLT_DiphotonMVA14p25_Tight_Mass90, &b_HLT_DiphotonMVA14p25_Tight_Mass90);
   fChain->SetBranchAddress("HLT_Photon35_TwoProngs35", &HLT_Photon35_TwoProngs35, &b_HLT_Photon35_TwoProngs35);
   fChain->SetBranchAddress("HLT_IsoMu24_TwoProngs35", &HLT_IsoMu24_TwoProngs35, &b_HLT_IsoMu24_TwoProngs35);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_NoOS", &HLT_Dimuon0_Jpsi_L1_NoOS, &b_HLT_Dimuon0_Jpsi_L1_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_NoOS", &HLT_Dimuon0_Jpsi_NoVertexing_NoOS, &b_HLT_Dimuon0_Jpsi_NoVertexing_NoOS);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi", &HLT_Dimuon0_Jpsi, &b_HLT_Dimuon0_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing", &HLT_Dimuon0_Jpsi_NoVertexing, &b_HLT_Dimuon0_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R", &HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R, &b_HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi3p5_Muon2", &HLT_Dimuon0_Jpsi3p5_Muon2, &b_HLT_Dimuon0_Jpsi3p5_Muon2);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5", &HLT_Dimuon0_Upsilon_L1_4p5, &b_HLT_Dimuon0_Upsilon_L1_4p5);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0", &HLT_Dimuon0_Upsilon_L1_4p5er2p0, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_L1_4p5er2p0M", &HLT_Dimuon0_Upsilon_L1_4p5er2p0M, &b_HLT_Dimuon0_Upsilon_L1_4p5er2p0M);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_NoVertexing", &HLT_Dimuon0_Upsilon_NoVertexing, &b_HLT_Dimuon0_Upsilon_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_0er1p5", &HLT_Dimuon0_LowMass_L1_0er1p5, &b_HLT_Dimuon0_LowMass_L1_0er1p5);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass", &HLT_Dimuon0_LowMass, &b_HLT_Dimuon0_LowMass);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_4", &HLT_Dimuon0_LowMass_L1_4, &b_HLT_Dimuon0_LowMass_L1_4);
   fChain->SetBranchAddress("HLT_Dimuon0_LowMass_L1_TM530", &HLT_Dimuon0_LowMass_L1_TM530, &b_HLT_Dimuon0_LowMass_L1_TM530);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon_NoL1Mass", &HLT_Dimuon0_Upsilon_Muon_NoL1Mass, &b_HLT_Dimuon0_Upsilon_Muon_NoL1Mass);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DZ", &HLT_TripleMu_5_3_3_Mass3p8_DZ, &b_HLT_TripleMu_5_3_3_Mass3p8_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_10_5_5_DZ", &HLT_TripleMu_10_5_5_DZ, &b_HLT_TripleMu_10_5_5_DZ);
   fChain->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5, &b_HLT_TripleMu_12_10_5);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15);
   fChain->SetBranchAddress("HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1", &HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1, &b_HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET50_PFMHT60", &HLT_DoubleMu3_DZ_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DZ_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET70_PFMHT70", &HLT_DoubleMu3_DZ_PFMET70_PFMHT70, &b_HLT_DoubleMu3_DZ_PFMET70_PFMHT70);
   fChain->SetBranchAddress("HLT_DoubleMu3_DZ_PFMET90_PFMHT90", &HLT_DoubleMu3_DZ_PFMET90_PFMHT90, &b_HLT_DoubleMu3_DZ_PFMET90_PFMHT90);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", &HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass, &b_HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_LowMass_SS", &HLT_DoubleMu4_3_LowMass_SS, &b_HLT_DoubleMu4_3_LowMass_SS);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_Displaced", &HLT_DoubleMu4_Jpsi_Displaced, &b_HLT_DoubleMu4_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_Jpsi_NoVertexing", &HLT_DoubleMu4_Jpsi_NoVertexing, &b_HLT_DoubleMu4_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrkTrk_Displaced", &HLT_DoubleMu4_JpsiTrkTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrkTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Bc", &HLT_DoubleMu4_JpsiTrk_Bc, &b_HLT_DoubleMu4_JpsiTrk_Bc);
   fChain->SetBranchAddress("HLT_DoubleMu43NoFiltersNoVtx", &HLT_DoubleMu43NoFiltersNoVtx, &b_HLT_DoubleMu43NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_DoubleMu48NoFiltersNoVtx", &HLT_DoubleMu48NoFiltersNoVtx, &b_HLT_DoubleMu48NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL, &b_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL", &HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL, &b_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL, &b_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL", &HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL, &b_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL);
   fChain->SetBranchAddress("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned, &b_HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned);
   fChain->SetBranchAddress("HLT_Ele28_eta2p1_WPTight_Gsf_HT150", &HLT_Ele28_eta2p1_WPTight_Gsf_HT150, &b_HLT_Ele28_eta2p1_WPTight_Gsf_HT150);
   fChain->SetBranchAddress("HLT_Ele28_HighEta_SC20_Mass55", &HLT_Ele28_HighEta_SC20_Mass55, &b_HLT_Ele28_HighEta_SC20_Mass55);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450_PFMET50", &HLT_Ele15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT450", &HLT_Ele15_IsoVVVL_PFHT450, &b_HLT_Ele15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450, &b_HLT_Ele50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450_PFMET50", &HLT_Mu15_IsoVVVL_PFHT450_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT450_PFMET50);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT450", &HLT_Mu15_IsoVVVL_PFHT450, &b_HLT_Mu15_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT450", &HLT_Mu50_IsoVVVL_PFHT450, &b_HLT_Mu50_IsoVVVL_PFHT450);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight);
   fChain->SetBranchAddress("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight, &b_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight);
   fChain->SetBranchAddress("HLT_Dimuon10_Upsilon_y1p4", &HLT_Dimuon10_Upsilon_y1p4, &b_HLT_Dimuon10_Upsilon_y1p4);
   fChain->SetBranchAddress("HLT_Dimuon12_Upsilon_y1p4", &HLT_Dimuon12_Upsilon_y1p4, &b_HLT_Dimuon12_Upsilon_y1p4);
   fChain->SetBranchAddress("HLT_Dimuon14_Phi_Barrel_Seagulls", &HLT_Dimuon14_Phi_Barrel_Seagulls, &b_HLT_Dimuon14_Phi_Barrel_Seagulls);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon14_PsiPrime", &HLT_Dimuon14_PsiPrime, &b_HLT_Dimuon14_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon14_PsiPrime_noCorrL1", &HLT_Dimuon14_PsiPrime_noCorrL1, &b_HLT_Dimuon14_PsiPrime_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime", &HLT_Dimuon18_PsiPrime, &b_HLT_Dimuon18_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon18_PsiPrime_noCorrL1", &HLT_Dimuon18_PsiPrime_noCorrL1, &b_HLT_Dimuon18_PsiPrime_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Upsilon_noCorrL1", &HLT_Dimuon24_Upsilon_noCorrL1, &b_HLT_Dimuon24_Upsilon_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon24_Phi_noCorrL1", &HLT_Dimuon24_Phi_noCorrL1, &b_HLT_Dimuon24_Phi_noCorrL1);
   fChain->SetBranchAddress("HLT_Dimuon25_Jpsi_noCorrL1", &HLT_Dimuon25_Jpsi_noCorrL1, &b_HLT_Dimuon25_Jpsi_noCorrL1);
   fChain->SetBranchAddress("HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8", &HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8, &b_HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_DoubleIsoMu20_eta2p1", &HLT_DoubleIsoMu20_eta2p1, &b_HLT_DoubleIsoMu20_eta2p1);
   fChain->SetBranchAddress("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx, &b_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu19", &HLT_Mu19, &b_HLT_Mu19);
   fChain->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId, &b_HLT_Mu17_Photon30_IsoCaloId);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53", &HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53, &b_HLT_Ele14_eta2p5_IsoVVVL_Gsf_PFHT200_PNetBTag0p53);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT, &b_HLT_Ele135_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40);
   fChain->SetBranchAddress("HLT_PFHT400_SixPFJet32", &HLT_PFHT400_SixPFJet32, &b_HLT_PFHT400_SixPFJet32);
   fChain->SetBranchAddress("HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50", &HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50, &b_HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50);
   fChain->SetBranchAddress("HLT_PFHT450_SixPFJet36", &HLT_PFHT450_SixPFJet36, &b_HLT_PFHT450_SixPFJet36);
   fChain->SetBranchAddress("HLT_PFHT450_SixPFJet36_PNetBTag0p35", &HLT_PFHT450_SixPFJet36_PNetBTag0p35, &b_HLT_PFHT450_SixPFJet36_PNetBTag0p35);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_120_120_60_30_30", &HLT_PFHT400_FivePFJet_120_120_60_30_30, &b_HLT_PFHT400_FivePFJet_120_120_60_30_30);
   fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400);
   fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   fChain->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   fChain->SetBranchAddress("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose, &b_HLT_Photon20_HoverELoose);
   fChain->SetBranchAddress("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose, &b_HLT_Photon30_HoverELoose);
   fChain->SetBranchAddress("HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   fChain->SetBranchAddress("HLT_CDC_L2cosmic_10_er1p0", &HLT_CDC_L2cosmic_10_er1p0, &b_HLT_CDC_L2cosmic_10_er1p0);
   fChain->SetBranchAddress("HLT_CDC_L2cosmic_5p5_er1p0", &HLT_CDC_L2cosmic_5p5_er1p0, &b_HLT_CDC_L2cosmic_5p5_er1p0);
   fChain->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1", &HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1, &b_HLT_VBF_DoubleMediumDeepTauPFTauHPS20_eta2p1);
   fChain->SetBranchAddress("HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3, &b_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3);
   fChain->SetBranchAddress("HLT_Mu18_Mu9_SameSign", &HLT_Mu18_Mu9_SameSign, &b_HLT_Mu18_Mu9_SameSign);
   fChain->SetBranchAddress("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05", &HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05, &b_HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05);
   fChain->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60, &b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60);
   fChain->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0, &b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0);
   fChain->SetBranchAddress("HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA", &HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA, &b_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_Mass2p0_noDCA);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3_Mass3p8_DCA", &HLT_TripleMu_5_3_3_Mass3p8_DCA, &b_HLT_TripleMu_5_3_3_Mass3p8_DCA);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15", &HLT_QuadPFJet103_88_75_15, &b_HLT_QuadPFJet103_88_75_15);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15", &HLT_QuadPFJet105_88_76_15, &b_HLT_QuadPFJet105_88_76_15);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15", &HLT_QuadPFJet111_90_80_15, &b_HLT_QuadPFJet111_90_80_15);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId, &b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55", &HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55, &b_HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55);
   fChain->SetBranchAddress("HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1", &HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1, &b_HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Tight_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Medium_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1", &HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1, &b_HLT_Ele24_eta2p1_WPTight_Gsf_PNetTauhPFJet30_Loose_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1", &HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1, &b_HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1", &HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1, &b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1", &HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1, &b_HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1);
   fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5);
   fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_4p3);
   fChain->SetBranchAddress("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0", &HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0, &b_HLT_PFHT330PT30_QuadPFJet_75_60_45_40_PNet3BTag_2p0);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6", &HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6, &b_HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_5p6);
   fChain->SetBranchAddress("HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3", &HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3, &b_HLT_PFHT400_FivePFJet_120_120_60_30_30_PNet2BTag_4p3);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2", &HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2, &b_HLT_QuadPFJet103_88_75_15_PNetBTag_0p4_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1", &HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1, &b_HLT_QuadPFJet103_88_75_15_PNet2BTag_0p4_0p12_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2", &HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2, &b_HLT_QuadPFJet105_88_76_15_PNetBTag_0p4_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1", &HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1, &b_HLT_QuadPFJet105_88_76_15_PNet2BTag_0p4_0p12_VBF1);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2", &HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2, &b_HLT_QuadPFJet111_90_80_15_PNetBTag_0p4_VBF2);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1", &HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1, &b_HLT_QuadPFJet111_90_80_15_PNet2BTag_0p4_0p12_VBF1);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5);
   fChain->SetBranchAddress("HLT_PFHT250_QuadPFJet25", &HLT_PFHT250_QuadPFJet25, &b_HLT_PFHT250_QuadPFJet25);
   fChain->SetBranchAddress("HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55", &HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55, &b_HLT_PFHT250_QuadPFJet25_PNet2BTagMean0p55);
   fChain->SetBranchAddress("HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50", &HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50, &b_HLT_PFHT250_QuadPFJet25_PNet1BTag0p20_PNet1Tauh0p50);
   fChain->SetBranchAddress("HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55", &HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55, &b_HLT_PFHT250_QuadPFJet30_PNet2BTagMean0p55);
   fChain->SetBranchAddress("HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50", &HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50, &b_HLT_PFHT250_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50);
   fChain->SetBranchAddress("HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50", &HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50, &b_HLT_PFHT280_QuadPFJet30_PNet1BTag0p20_PNet1Tauh0p50);
   fChain->SetBranchAddress("HLT_PFHT280_QuadPFJet30", &HLT_PFHT280_QuadPFJet30, &b_HLT_PFHT280_QuadPFJet30);
   fChain->SetBranchAddress("HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55", &HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55, &b_HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55);
   fChain->SetBranchAddress("HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60", &HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60, &b_HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60);
   fChain->SetBranchAddress("HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60", &HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60, &b_HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60);
   fChain->SetBranchAddress("HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70", &HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70, &b_HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT280_QuadPFJet30_PNet2BTagMean0p55);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PNet2BTagMean0p50);
   fChain->SetBranchAddress("HLT_QuadPFJet100_88_70_30", &HLT_QuadPFJet100_88_70_30, &b_HLT_QuadPFJet100_88_70_30);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_75_30", &HLT_QuadPFJet105_88_75_30, &b_HLT_QuadPFJet105_88_75_30);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_30", &HLT_QuadPFJet111_90_80_30, &b_HLT_QuadPFJet111_90_80_30);
   fChain->SetBranchAddress("HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight", &HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight, &b_HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight);
   fChain->SetBranchAddress("HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight", &HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight, &b_HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight);
   fChain->SetBranchAddress("HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight", &HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight, &b_HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40", &HLT_AK8PFJet220_SoftDropMass40, &b_HLT_AK8PFJet220_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p50);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p53);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p55);
   fChain->SetBranchAddress("HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60", &HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60, &b_HLT_AK8PFJet220_SoftDropMass40_PNetBB0p06_DoubleAK4PFJet60_30_PNet2BTagMean0p60);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40", &HLT_AK8PFJet230_SoftDropMass40, &b_HLT_AK8PFJet230_SoftDropMass40);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06, &b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10", &HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10, &b_HLT_AK8PFJet230_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03", &HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03, &b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03);
   fChain->SetBranchAddress("HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05", &HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05, &b_HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p05);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06, &b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10", &HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10, &b_HLT_AK8PFJet250_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03", &HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03, &b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p03);
   fChain->SetBranchAddress("HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05", &HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05, &b_HLT_AK8PFJet250_SoftDropMass40_PNetTauTau0p05);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06", &HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06, &b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10", &HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10, &b_HLT_AK8PFJet275_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03", &HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03, &b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p03);
   fChain->SetBranchAddress("HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05", &HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05, &b_HLT_AK8PFJet275_SoftDropMass40_PNetTauTau0p05);
   fChain->SetBranchAddress("HLT_AK8PFJet275_Nch45", &HLT_AK8PFJet275_Nch45, &b_HLT_AK8PFJet275_Nch45);
   fChain->SetBranchAddress("HLT_AK8PFJet275_Nch40", &HLT_AK8PFJet275_Nch40, &b_HLT_AK8PFJet275_Nch40);
   fChain->SetBranchAddress("HLT_IsoMu50_AK8PFJet220_SoftDropMass40", &HLT_IsoMu50_AK8PFJet220_SoftDropMass40, &b_HLT_IsoMu50_AK8PFJet220_SoftDropMass40);
   fChain->SetBranchAddress("HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06", &HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06, &b_HLT_IsoMu50_AK8PFJet220_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40, &b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40);
   fChain->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06, &b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10", &HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10, &b_HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p06);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_DoubleAK4PFJet60_30_PNet2BTagMean0p50);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet2BTagMean0p55);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFHT250_QuadPFJet25_PNet1BTag0p20);
   fChain->SetBranchAddress("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60", &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60, &b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60);
   fChain->SetBranchAddress("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75", &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75, &b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_L2NN_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60", &HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet60);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75", &HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1_PFJet75);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring", &HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Tight_eta2p3_CrossL1_ETau_Monitoring);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring", &HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_eta2p3_CrossL1_ETau_Monitoring);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring", &HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Loose_eta2p3_CrossL1_ETau_Monitoring);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1", &HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1, &b_HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Tight_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1", &HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1, &b_HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1", &HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1, &b_HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Loose_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet20_eta2p2_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet45_L2NN_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PFHT250", &HLT_IsoMu24_eta2p1_PFHT250, &b_HLT_IsoMu24_eta2p1_PFHT250);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25", &HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25, &b_HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50", &HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50, &b_HLT_IsoMu24_eta2p1_PFHT250_QuadPFJet25_PNet1Tauh0p50);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50", &HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50, &b_HLT_IsoMu24_eta2p1_SinglePFJet25_PNet1Tauh0p50);
   fChain->SetBranchAddress("HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3", &HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3, &b_HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3);
   fChain->SetBranchAddress("HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3", &HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3, &b_HLT_DoublePNetTauhPFJet30_Tight_L2NN_eta2p3);
   fChain->SetBranchAddress("HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60", &HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60, &b_HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet60);
   fChain->SetBranchAddress("HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75", &HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75, &b_HLT_DoublePNetTauhPFJet26_L2NN_eta2p3_PFJet75);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Loose_L2NN_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Medium_L2NN_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1", &HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1, &b_HLT_IsoMu24_eta2p1_PNetTauhPFJet130_Tight_L2NN_eta2p3_CrossL1);
   fChain->SetBranchAddress("HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3", &HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3, &b_HLT_SinglePNetTauhPFJet130_Loose_L2NN_eta2p3);
   fChain->SetBranchAddress("HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3", &HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3, &b_HLT_SinglePNetTauhPFJet130_Medium_L2NN_eta2p3);
   fChain->SetBranchAddress("HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3", &HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3, &b_HLT_SinglePNetTauhPFJet130_Tight_L2NN_eta2p3);
   fChain->SetBranchAddress("HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu12NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu14NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3Mu18_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3Mu20_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha", &HLT_L2Mu10NoVtx_2Cha, &b_HLT_L2Mu10NoVtx_2Cha);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm", &HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm, &b_HLT_L2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_L3Mu10NoVtx", &HLT_L3Mu10NoVtx, &b_HLT_L3Mu10NoVtx);
   fChain->SetBranchAddress("HLT_L3Mu10NoVtx_DxyMin0p01cm", &HLT_L3Mu10NoVtx_DxyMin0p01cm, &b_HLT_L3Mu10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm", &HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm, &b_HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm", &HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm, &b_HLT_DoubleL2Mu_L3Mu18NoVtx_VetoL3Mu0DxyMax0p1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm", &HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm, &b_HLT_DoubleL2Mu12NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha_CosmicSeed", &HLT_L2Mu10NoVtx_2Cha_CosmicSeed, &b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed);
   fChain->SetBranchAddress("HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm", &HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm, &b_HLT_L2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm", &HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm, &b_HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm", &HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm, &b_HLT_L2Mu50NoVtx_3Cha_CosmicSeed_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm", &HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm, &b_HLT_L2Mu50NoVtx_3Cha_VetoL3Mu0DxyMax1cm);
   fChain->SetBranchAddress("HLT_L3Mu30NoVtx_DxyMin0p01cm", &HLT_L3Mu30NoVtx_DxyMin0p01cm, &b_HLT_L3Mu30NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_L3Mu50NoVtx_DxyMin0p01cm", &HLT_L3Mu50NoVtx_DxyMin0p01cm, &b_HLT_L3Mu50NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_L3dTksMu10_NoVtx_DxyMin0p01cm", &HLT_L3dTksMu10_NoVtx_DxyMin0p01cm, &b_HLT_L3dTksMu10_NoVtx_DxyMin0p01cm);
   fChain->SetBranchAddress("HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId", &HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId, &b_HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId);
   fChain->SetBranchAddress("HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24", &HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24, &b_HLT_DisplacedMu24_MediumChargedIsoDisplacedPFTauHPS24);
   fChain->SetBranchAddress("HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34", &HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34, &b_HLT_Photon34_R9Id90_CaloIdL_IsoL_DisplacedIdL_MediumChargedIsoDisplacedPFTauHPS34);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1", &HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1, &b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1", &HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1, &b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1);
   fChain->SetBranchAddress("HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy", &HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy, &b_HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_noDxy);
   fChain->SetBranchAddress("HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1", &HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1, &b_HLT_IsoMu27_MediumChargedIsoDisplacedPFTauHPS24_eta2p1_SingleL1);
   fChain->SetBranchAddress("HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive", &HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive, &b_HLT_HT350_DelayedJet40_SingleDelay1p5To3p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive", &HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive, &b_HLT_HT350_DelayedJet40_SingleDelay1p6To3p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive", &HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive, &b_HLT_HT350_DelayedJet40_SingleDelay1p75To3p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT350_DelayedJet40_SingleDelay3nsInclusive", &HLT_HT350_DelayedJet40_SingleDelay3nsInclusive, &b_HLT_HT350_DelayedJet40_SingleDelay3nsInclusive);
   fChain->SetBranchAddress("HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive", &HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive, &b_HLT_HT350_DelayedJet40_SingleDelay3p25nsInclusive);
   fChain->SetBranchAddress("HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive", &HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive, &b_HLT_HT350_DelayedJet40_SingleDelay3p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay1nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay1nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay1To1p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay1p1To1p6nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay1p25To1p75nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay2nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay2nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay2nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive", &HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive, &b_HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless", &HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless, &b_HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1nsTrackless", &HLT_HT430_DelayedJet40_SingleDelay1nsTrackless, &b_HLT_HT430_DelayedJet40_SingleDelay1nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless", &HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless, &b_HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless", &HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless, &b_HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive", &HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive, &b_HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive", &HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive, &b_HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive", &HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive, &b_HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive", &HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive, &b_HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless", &HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless, &b_HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless", &HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless, &b_HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless);
   fChain->SetBranchAddress("HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless", &HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless, &b_HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive", &HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive, &b_HLT_L1Tau_DelayedJet40_SingleDelay2p5To4nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive", &HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive, &b_HLT_L1Tau_DelayedJet40_SingleDelay2p6To4nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive", &HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive, &b_HLT_L1Tau_DelayedJet40_SingleDelay2p75To4nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive", &HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive, &b_HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive", &HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive, &b_HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive", &HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive, &b_HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless", &HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless, &b_HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless", &HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless, &b_HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless", &HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless, &b_HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive", &HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive, &b_HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive", &HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive, &b_HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive", &HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive, &b_HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive", &HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive, &b_HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless", &HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless, &b_HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless", &HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless, &b_HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless", &HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless, &b_HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless);
   fChain->SetBranchAddress("HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless", &HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless, &b_HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless);
   fChain->SetBranchAddress("HLT_L1Mu6HT240", &HLT_L1Mu6HT240, &b_HLT_L1Mu6HT240);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose", &HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose, &b_HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5", &HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5, &b_HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose", &HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose, &b_HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5", &HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5, &b_HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose", &HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose, &b_HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5", &HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5, &b_HLT_Mu6HT240_DisplacedDijet45_Inclusive0PtrkShortSig5);
   fChain->SetBranchAddress("HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5", &HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5, &b_HLT_Mu6HT240_DisplacedDijet50_Inclusive0PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT350", &HLT_HT350, &b_HLT_HT350);
   fChain->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
   fChain->SetBranchAddress("HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT360_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5", &HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5, &b_HLT_HT360_DisplacedDijet45_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT390_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5", &HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5, &b_HLT_HT390_DisplacedDijet45_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT390eta2p0_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT430_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack", &HLT_HT400_DisplacedDijet40_DisplacedTrack, &b_HLT_HT400_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack", &HLT_HT430_DisplacedDijet40_DisplacedTrack, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive", &HLT_HT550_DisplacedDijet60_Inclusive, &b_HLT_HT550_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive", &HLT_HT650_DisplacedDijet60_Inclusive, &b_HLT_HT650_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_CaloMET60_DTCluster50", &HLT_CaloMET60_DTCluster50, &b_HLT_CaloMET60_DTCluster50);
   fChain->SetBranchAddress("HLT_CaloMET60_DTClusterNoMB1S50", &HLT_CaloMET60_DTClusterNoMB1S50, &b_HLT_CaloMET60_DTClusterNoMB1S50);
   fChain->SetBranchAddress("HLT_L1MET_DTCluster50", &HLT_L1MET_DTCluster50, &b_HLT_L1MET_DTCluster50);
   fChain->SetBranchAddress("HLT_L1MET_DTClusterNoMB1S50", &HLT_L1MET_DTClusterNoMB1S50, &b_HLT_L1MET_DTClusterNoMB1S50);
   fChain->SetBranchAddress("HLT_CscCluster_Loose", &HLT_CscCluster_Loose, &b_HLT_CscCluster_Loose);
   fChain->SetBranchAddress("HLT_CscCluster_Medium", &HLT_CscCluster_Medium, &b_HLT_CscCluster_Medium);
   fChain->SetBranchAddress("HLT_CscCluster_Tight", &HLT_CscCluster_Tight, &b_HLT_CscCluster_Tight);
   fChain->SetBranchAddress("HLT_CscCluster50_Photon20Unseeded", &HLT_CscCluster50_Photon20Unseeded, &b_HLT_CscCluster50_Photon20Unseeded);
   fChain->SetBranchAddress("HLT_CscCluster50_Photon30Unseeded", &HLT_CscCluster50_Photon30Unseeded, &b_HLT_CscCluster50_Photon30Unseeded);
   fChain->SetBranchAddress("HLT_CscCluster100_Ele5", &HLT_CscCluster100_Ele5, &b_HLT_CscCluster100_Ele5);
   fChain->SetBranchAddress("HLT_CscCluster100_Mu5", &HLT_CscCluster100_Mu5, &b_HLT_CscCluster100_Mu5);
   fChain->SetBranchAddress("HLT_CscCluster100_PNetTauhPFJet10_Loose", &HLT_CscCluster100_PNetTauhPFJet10_Loose, &b_HLT_CscCluster100_PNetTauhPFJet10_Loose);
   fChain->SetBranchAddress("HLT_DoubleCscCluster75", &HLT_DoubleCscCluster75, &b_HLT_DoubleCscCluster75);
   fChain->SetBranchAddress("HLT_IsoTrk200_L1SingleMuShower", &HLT_IsoTrk200_L1SingleMuShower, &b_HLT_IsoTrk200_L1SingleMuShower);
   fChain->SetBranchAddress("HLT_IsoTrk400_L1SingleMuShower", &HLT_IsoTrk400_L1SingleMuShower, &b_HLT_IsoTrk400_L1SingleMuShower);
   fChain->SetBranchAddress("HLT_DoubleCscCluster100", &HLT_DoubleCscCluster100, &b_HLT_DoubleCscCluster100);
   fChain->SetBranchAddress("HLT_L1CSCShower_DTCluster50", &HLT_L1CSCShower_DTCluster50, &b_HLT_L1CSCShower_DTCluster50);
   fChain->SetBranchAddress("HLT_L1CSCShower_DTCluster75", &HLT_L1CSCShower_DTCluster75, &b_HLT_L1CSCShower_DTCluster75);
   fChain->SetBranchAddress("HLT_PFMET105_IsoTrk50", &HLT_PFMET105_IsoTrk50, &b_HLT_PFMET105_IsoTrk50);
   fChain->SetBranchAddress("HLT_L1SingleLLPJet", &HLT_L1SingleLLPJet, &b_HLT_L1SingleLLPJet);
   fChain->SetBranchAddress("HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack", &HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack, &b_HLT_HT170_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack", &HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack", &HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet60_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack", &HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack, &b_HLT_HT270_L1SingleLLPJet_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT200_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT240_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5", &HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5, &b_HLT_HT280_L1SingleLLPJet_DisplacedDijet40_Inclusive1PtrkShortSig5);
   fChain->SetBranchAddress("HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive", &HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive, &b_HLT_HT320_L1SingleLLPJet_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive", &HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive, &b_HLT_HT420_L1SingleLLPJet_DisplacedDijet60_Inclusive);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless", &HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay1nsTrackless);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive", &HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_SingleDelay2nsInclusive);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless", &HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive", &HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive, &b_HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay1nsInclusive);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7", &HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7, &b_HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p7);
   fChain->SetBranchAddress("HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8", &HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8, &b_HLT_HT200_L1SingleLLPJet_PFJet60_NeutralHadronFrac0p8);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1ns", &HLT_DiPhoton10Time1ns, &b_HLT_DiPhoton10Time1ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p2ns", &HLT_DiPhoton10Time1p2ns, &b_HLT_DiPhoton10Time1p2ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p4ns", &HLT_DiPhoton10Time1p4ns, &b_HLT_DiPhoton10Time1p4ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p6ns", &HLT_DiPhoton10Time1p6ns, &b_HLT_DiPhoton10Time1p6ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time1p8ns", &HLT_DiPhoton10Time1p8ns, &b_HLT_DiPhoton10Time1p8ns);
   fChain->SetBranchAddress("HLT_DiPhoton10Time2ns", &HLT_DiPhoton10Time2ns, &b_HLT_DiPhoton10Time2ns);
   fChain->SetBranchAddress("HLT_DiPhoton10_CaloIdL", &HLT_DiPhoton10_CaloIdL, &b_HLT_DiPhoton10_CaloIdL);
   fChain->SetBranchAddress("HLT_DoubleEle6p5_eta1p22_mMax6", &HLT_DoubleEle6p5_eta1p22_mMax6, &b_HLT_DoubleEle6p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle8_eta1p22_mMax6", &HLT_DoubleEle8_eta1p22_mMax6, &b_HLT_DoubleEle8_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle10_eta1p22_mMax6", &HLT_DoubleEle10_eta1p22_mMax6, &b_HLT_DoubleEle10_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_SingleEle8", &HLT_SingleEle8, &b_HLT_SingleEle8);
   fChain->SetBranchAddress("HLT_SingleEle8_SingleEGL1", &HLT_SingleEle8_SingleEGL1, &b_HLT_SingleEle8_SingleEGL1);
   fChain->SetBranchAddress("HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT", &HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT, &b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton20_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT", &HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT, &b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HE_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton22_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton24_14_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT", &HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT, &b_HLT_Diphoton24_16_eta1p5_R9IdL_AND_HET_AND_IsoTCaloIdT);
   fChain->SetBranchAddress("HLT_Mu50_L1SingleMuShower", &HLT_Mu50_L1SingleMuShower, &b_HLT_Mu50_L1SingleMuShower);
   fChain->SetBranchAddress("HLT_IsoMu24_OneProng32", &HLT_IsoMu24_OneProng32, &b_HLT_IsoMu24_OneProng32);
   fChain->SetBranchAddress("HLT_Photon32_OneProng32_M50To105", &HLT_Photon32_OneProng32_M50To105, &b_HLT_Photon32_OneProng32_M50To105);
   fChain->SetBranchAddress("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng", &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng, &b_HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1", &HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1, &b_HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_CrossL1);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet125_45_Mjj1050", &HLT_VBF_DiPFJet125_45_Mjj1050, &b_HLT_VBF_DiPFJet125_45_Mjj1050);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet125_45_Mjj1200", &HLT_VBF_DiPFJet125_45_Mjj1200, &b_HLT_VBF_DiPFJet125_45_Mjj1200);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60", &HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60, &b_HLT_VBF_DiPFJet75_45_Mjj800_DiPFJet60);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60", &HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60, &b_HLT_VBF_DiPFJet75_45_Mjj850_DiPFJet60);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85", &HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85, &b_HLT_VBF_DiPFJet80_45_Mjj650_PFMETNoMu85);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85", &HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85, &b_HLT_VBF_DiPFJet80_45_Mjj750_PFMETNoMu85);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL", &HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL, &b_HLT_VBF_DiPFJet95_45_Mjj750_Mu3_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL", &HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL, &b_HLT_VBF_DiPFJet95_45_Mjj850_Mu3_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet50_Mjj650_Photon22", &HLT_VBF_DiPFJet50_Mjj650_Photon22, &b_HLT_VBF_DiPFJet50_Mjj650_Photon22);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet50_Mjj750_Photon22", &HLT_VBF_DiPFJet50_Mjj750_Photon22, &b_HLT_VBF_DiPFJet50_Mjj750_Photon22);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf", &HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf, &b_HLT_VBF_DiPFJet50_Mjj600_Ele22_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf", &HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf, &b_HLT_VBF_DiPFJet50_Mjj650_Ele22_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1", &HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1, &b_HLT_VBF_DiPFJet45_Mjj650_MediumDeepTauPFTauHPS45_L2NN_eta2p1);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1", &HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1, &b_HLT_VBF_DiPFJet45_Mjj750_MediumDeepTauPFTauHPS45_L2NN_eta2p1);
   fChain->SetBranchAddress("HLT_VBF_DoublePNetTauhPFJet20_eta2p2", &HLT_VBF_DoublePNetTauhPFJet20_eta2p2, &b_HLT_VBF_DoublePNetTauhPFJet20_eta2p2);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3", &HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3, &b_HLT_VBF_DiPFJet45_Mjj650_PNetTauhPFJet45_L2NN_eta2p3);
   fChain->SetBranchAddress("HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3", &HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3, &b_HLT_VBF_DiPFJet45_Mjj750_PNetTauhPFJet45_L2NN_eta2p3);
   fChain->SetBranchAddress("HLT_PFJet200_TimeLtNeg2p5ns", &HLT_PFJet200_TimeLtNeg2p5ns, &b_HLT_PFJet200_TimeLtNeg2p5ns);
   fChain->SetBranchAddress("HLT_PFJet200_TimeGt2p5ns", &HLT_PFJet200_TimeGt2p5ns, &b_HLT_PFJet200_TimeGt2p5ns);
   fChain->SetBranchAddress("HLT_Photon50_TimeLtNeg2p5ns", &HLT_Photon50_TimeLtNeg2p5ns, &b_HLT_Photon50_TimeLtNeg2p5ns);
   fChain->SetBranchAddress("HLT_Photon50_TimeGt2p5ns", &HLT_Photon50_TimeGt2p5ns, &b_HLT_Photon50_TimeGt2p5ns);
   fChain->SetBranchAddress("HLT_PPSMaxTracksPerArm1", &HLT_PPSMaxTracksPerArm1, &b_HLT_PPSMaxTracksPerArm1);
   fChain->SetBranchAddress("HLT_PPSMaxTracksPerRP4", &HLT_PPSMaxTracksPerRP4, &b_HLT_PPSMaxTracksPerRP4);
   fChain->SetBranchAddress("HLT_PPSRandom", &HLT_PPSRandom, &b_HLT_PPSRandom);
   fChain->SetBranchAddress("HLT_L1AXOVTight", &HLT_L1AXOVTight, &b_HLT_L1AXOVTight);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("Dataset_ScoutingPFMonitor", &Dataset_ScoutingPFMonitor, &b_Dataset_ScoutingPFMonitor);
   fChain->SetBranchAddress("Dataset_ScoutingPFRun3", &Dataset_ScoutingPFRun3, &b_Dataset_ScoutingPFRun3);
   */
   Notify();
}

Bool_t DepthsHistosFill::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DepthsHistosFill::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DepthsHistosFill::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DepthsHistosFill_cxx
