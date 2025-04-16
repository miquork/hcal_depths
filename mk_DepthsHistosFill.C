#include "DepthsHistosFill.h"

//#define GPU
#define LOCAL

// Compile these libraries into *.so first with root -l -b -q mk_compile.C
#ifdef LOCAL
// (works for 6.18.04?)
R__LOAD_LIBRARY(DepthsHistosFill_C)
#else
// (works for 6.26/10)
R__LOAD_LIBRARY(DepthsHistosFill_C.so)
#endif

void mk_DepthsHistosFill() {

  TChain *c = new TChain("Events");
  //c->AddFile("../data/dijet/data/NanoMerged_Data24Iv2_JetMET1_TeVJet_Baseline_RAWToPFNANO_v0p1.root");
  //c->AddFile("../data/dijet/data/NanoMerged_Data24Iv2_JetMET1_TeVJet_UpdatedGainsRespCorrs_RAWToPFNANO_v0p1.root");

  for (int i = 10; i != 50; ++i)
    c->AddFile(Form("../data/dijet/mc/QCD-4Jets_Bin-HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8_RunIII2024Summer24NanoAOD_Baseline_RAWToPFNANO_v0p1/NANOAODSIM_%d.root",i));

  DepthsHistosFill filler(c);
  filler.Loop();
}
