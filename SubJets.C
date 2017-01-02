#ifndef SubJets_cxx
#define SubJets_cxx
#include "SubJets.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void SubJets::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SubJets.C
//      Root > SubJets t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}


SubJets::SubJets(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/opt/sbg/cms/ui1_data1/vanhove/NTUPLES/CMSSW_7_6_3/MC/QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/Pt-80to120/JetTree_mc_FatJets_Subjets_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/opt/sbg/cms/ui1_data1/vanhove/NTUPLES/CMSSW_7_6_3/MC/QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/Pt-80to120/JetTree_mc_FatJets_Subjets_1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/opt/sbg/cms/ui1_data1/vanhove/NTUPLES/CMSSW_7_6_3/MC/QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/Pt-80to120/JetTree_mc_FatJets_Subjets_1.root:/btaganaFatJets");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

SubJets::~SubJets()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SubJets::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SubJets::LoadTree(Long64_t entry)
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

void SubJets::Init(TTree *tree)
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

   fChain->SetBranchAddress("FatJetInfo.nJet", &FatJetInfo_nJet, &b_FatJetInfo_nJet);
   fChain->SetBranchAddress("FatJetInfo.Jet_pt", FatJetInfo_Jet_pt, &b_FatJetInfo_Jet_pt);
   fChain->SetBranchAddress("FatJetInfo.Jet_genpt", FatJetInfo_Jet_genpt, &b_FatJetInfo_Jet_genpt);
   fChain->SetBranchAddress("FatJetInfo.Jet_residual", FatJetInfo_Jet_residual, &b_FatJetInfo_Jet_residual);
   fChain->SetBranchAddress("FatJetInfo.Jet_area", FatJetInfo_Jet_area, &b_FatJetInfo_Jet_area);
   fChain->SetBranchAddress("FatJetInfo.Jet_jes", FatJetInfo_Jet_jes, &b_FatJetInfo_Jet_jes);
   fChain->SetBranchAddress("FatJetInfo.Jet_eta", FatJetInfo_Jet_eta, &b_FatJetInfo_Jet_eta);
   fChain->SetBranchAddress("FatJetInfo.Jet_phi", FatJetInfo_Jet_phi, &b_FatJetInfo_Jet_phi);
   fChain->SetBranchAddress("FatJetInfo.Jet_mass", FatJetInfo_Jet_mass, &b_FatJetInfo_Jet_mass);
   fChain->SetBranchAddress("FatJetInfo.Jet_ntracks", FatJetInfo_Jet_ntracks, &b_FatJetInfo_Jet_ntracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_nseltracks", FatJetInfo_Jet_nseltracks, &b_FatJetInfo_Jet_nseltracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_flavour", FatJetInfo_Jet_flavour, &b_FatJetInfo_Jet_flavour);
   fChain->SetBranchAddress("FatJetInfo.Jet_nbHadrons", FatJetInfo_Jet_nbHadrons, &b_FatJetInfo_Jet_nbHadrons);
   fChain->SetBranchAddress("FatJetInfo.Jet_ncHadrons", FatJetInfo_Jet_ncHadrons, &b_FatJetInfo_Jet_ncHadrons);
   fChain->SetBranchAddress("FatJetInfo.Jet_Ip2N", FatJetInfo_Jet_Ip2N, &b_FatJetInfo_Jet_Ip2N);
   fChain->SetBranchAddress("FatJetInfo.Jet_Ip2P", FatJetInfo_Jet_Ip2P, &b_FatJetInfo_Jet_Ip2P);
   fChain->SetBranchAddress("FatJetInfo.Jet_Ip3N", FatJetInfo_Jet_Ip3N, &b_FatJetInfo_Jet_Ip3N);
   fChain->SetBranchAddress("FatJetInfo.Jet_Ip3P", FatJetInfo_Jet_Ip3P, &b_FatJetInfo_Jet_Ip3P);
   fChain->SetBranchAddress("FatJetInfo.Jet_ProbaN", FatJetInfo_Jet_ProbaN, &b_FatJetInfo_Jet_ProbaN);
   fChain->SetBranchAddress("FatJetInfo.Jet_ProbaP", FatJetInfo_Jet_ProbaP, &b_FatJetInfo_Jet_ProbaP);
   fChain->SetBranchAddress("FatJetInfo.Jet_Proba", FatJetInfo_Jet_Proba, &b_FatJetInfo_Jet_Proba);
   fChain->SetBranchAddress("FatJetInfo.Jet_BprobN", FatJetInfo_Jet_BprobN, &b_FatJetInfo_Jet_BprobN);
   fChain->SetBranchAddress("FatJetInfo.Jet_BprobP", FatJetInfo_Jet_BprobP, &b_FatJetInfo_Jet_BprobP);
   fChain->SetBranchAddress("FatJetInfo.Jet_Bprob", FatJetInfo_Jet_Bprob, &b_FatJetInfo_Jet_Bprob);
   fChain->SetBranchAddress("FatJetInfo.Jet_SvxN", FatJetInfo_Jet_SvxN, &b_FatJetInfo_Jet_SvxN);
   fChain->SetBranchAddress("FatJetInfo.Jet_Svx", FatJetInfo_Jet_Svx, &b_FatJetInfo_Jet_Svx);
   fChain->SetBranchAddress("FatJetInfo.Jet_SvxNHP", FatJetInfo_Jet_SvxNHP, &b_FatJetInfo_Jet_SvxNHP);
   fChain->SetBranchAddress("FatJetInfo.Jet_SvxHP", FatJetInfo_Jet_SvxHP, &b_FatJetInfo_Jet_SvxHP);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombSvxN", FatJetInfo_Jet_CombSvxN, &b_FatJetInfo_Jet_CombSvxN);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombSvxP", FatJetInfo_Jet_CombSvxP, &b_FatJetInfo_Jet_CombSvxP);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombSvx", FatJetInfo_Jet_CombSvx, &b_FatJetInfo_Jet_CombSvx);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombIVF", FatJetInfo_Jet_CombIVF, &b_FatJetInfo_Jet_CombIVF);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombIVF_P", FatJetInfo_Jet_CombIVF_P, &b_FatJetInfo_Jet_CombIVF_P);
   fChain->SetBranchAddress("FatJetInfo.Jet_CombIVF_N", FatJetInfo_Jet_CombIVF_N, &b_FatJetInfo_Jet_CombIVF_N);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftMuN", FatJetInfo_Jet_SoftMuN, &b_FatJetInfo_Jet_SoftMuN);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftMuP", FatJetInfo_Jet_SoftMuP, &b_FatJetInfo_Jet_SoftMuP);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftMu", FatJetInfo_Jet_SoftMu, &b_FatJetInfo_Jet_SoftMu);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftElN", FatJetInfo_Jet_SoftElN, &b_FatJetInfo_Jet_SoftElN);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftElP", FatJetInfo_Jet_SoftElP, &b_FatJetInfo_Jet_SoftElP);
   fChain->SetBranchAddress("FatJetInfo.Jet_SoftEl", FatJetInfo_Jet_SoftEl, &b_FatJetInfo_Jet_SoftEl);
   fChain->SetBranchAddress("FatJetInfo.Jet_DoubleSV", FatJetInfo_Jet_DoubleSV, &b_FatJetInfo_Jet_DoubleSV);
   fChain->SetBranchAddress("FatJetInfo.Jet_cMVA", FatJetInfo_Jet_cMVA, &b_FatJetInfo_Jet_cMVA);
   fChain->SetBranchAddress("FatJetInfo.Jet_cMVAv2", FatJetInfo_Jet_cMVAv2, &b_FatJetInfo_Jet_cMVAv2);
   fChain->SetBranchAddress("FatJetInfo.Jet_cMVAv2N", FatJetInfo_Jet_cMVAv2N, &b_FatJetInfo_Jet_cMVAv2N);
   fChain->SetBranchAddress("FatJetInfo.Jet_cMVAv2P", FatJetInfo_Jet_cMVAv2P, &b_FatJetInfo_Jet_cMVAv2P);
   fChain->SetBranchAddress("FatJetInfo.Jet_hist1", FatJetInfo_Jet_hist1, &b_FatJetInfo_Jet_hist1);
   fChain->SetBranchAddress("FatJetInfo.Jet_hist2", FatJetInfo_Jet_hist2, &b_FatJetInfo_Jet_hist2);
   fChain->SetBranchAddress("FatJetInfo.Jet_hist3", FatJetInfo_Jet_hist3, &b_FatJetInfo_Jet_hist3);
   fChain->SetBranchAddress("FatJetInfo.Jet_histJet", FatJetInfo_Jet_histJet, &b_FatJetInfo_Jet_histJet);
   fChain->SetBranchAddress("FatJetInfo.Jet_histSvx", FatJetInfo_Jet_histSvx, &b_FatJetInfo_Jet_histSvx);
   fChain->SetBranchAddress("FatJetInfo.Jet_SV_multi", FatJetInfo_Jet_SV_multi, &b_FatJetInfo_Jet_SV_multi);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSM", FatJetInfo_Jet_nSM, &b_FatJetInfo_Jet_nSM);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSE", FatJetInfo_Jet_nSE, &b_FatJetInfo_Jet_nSE);
   fChain->SetBranchAddress("FatJetInfo.Jet_looseID", FatJetInfo_Jet_looseID, &b_FatJetInfo_Jet_looseID);
   fChain->SetBranchAddress("FatJetInfo.Jet_tightID", FatJetInfo_Jet_tightID, &b_FatJetInfo_Jet_tightID);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSE", FatJetInfo_Jet_nFirstSE, &b_FatJetInfo_Jet_nFirstSE);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSE", FatJetInfo_Jet_nLastSE, &b_FatJetInfo_Jet_nLastSE);
   fChain->SetBranchAddress("FatJetInfo.nPFElectron", &FatJetInfo_nPFElectron, &b_FatJetInfo_nPFElectron);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_IdxJet", FatJetInfo_PFElectron_IdxJet, &b_FatJetInfo_PFElectron_IdxJet);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_pt", FatJetInfo_PFElectron_pt, &b_FatJetInfo_PFElectron_pt);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_eta", FatJetInfo_PFElectron_eta, &b_FatJetInfo_PFElectron_eta);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_phi", FatJetInfo_PFElectron_phi, &b_FatJetInfo_PFElectron_phi);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_ptrel", FatJetInfo_PFElectron_ptrel, &b_FatJetInfo_PFElectron_ptrel);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_deltaR", FatJetInfo_PFElectron_deltaR, &b_FatJetInfo_PFElectron_deltaR);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_ratio", FatJetInfo_PFElectron_ratio, &b_FatJetInfo_PFElectron_ratio);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_ratioRel", FatJetInfo_PFElectron_ratioRel, &b_FatJetInfo_PFElectron_ratioRel);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_IP", FatJetInfo_PFElectron_IP, &b_FatJetInfo_PFElectron_IP);
   fChain->SetBranchAddress("FatJetInfo.PFElectron_IP2D", FatJetInfo_PFElectron_IP2D, &b_FatJetInfo_PFElectron_IP2D);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSM", FatJetInfo_Jet_nFirstSM, &b_FatJetInfo_Jet_nFirstSM);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSM", FatJetInfo_Jet_nLastSM, &b_FatJetInfo_Jet_nLastSM);
   fChain->SetBranchAddress("FatJetInfo.nPFMuon", &FatJetInfo_nPFMuon, &b_FatJetInfo_nPFMuon);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IdxJet", FatJetInfo_PFMuon_IdxJet, &b_FatJetInfo_PFMuon_IdxJet);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nMuHit", FatJetInfo_PFMuon_nMuHit, &b_FatJetInfo_PFMuon_nMuHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nTkHit", FatJetInfo_PFMuon_nTkHit, &b_FatJetInfo_PFMuon_nTkHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nPixHit", FatJetInfo_PFMuon_nPixHit, &b_FatJetInfo_PFMuon_nPixHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nOutHit", FatJetInfo_PFMuon_nOutHit, &b_FatJetInfo_PFMuon_nOutHit);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nTkLwM", FatJetInfo_PFMuon_nTkLwM, &b_FatJetInfo_PFMuon_nTkLwM);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nPixLwM", FatJetInfo_PFMuon_nPixLwM, &b_FatJetInfo_PFMuon_nPixLwM);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_nMatched", FatJetInfo_PFMuon_nMatched, &b_FatJetInfo_PFMuon_nMatched);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_chi2", FatJetInfo_PFMuon_chi2, &b_FatJetInfo_PFMuon_chi2);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_chi2Tk", FatJetInfo_PFMuon_chi2Tk, &b_FatJetInfo_PFMuon_chi2Tk);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_isGlobal", FatJetInfo_PFMuon_isGlobal, &b_FatJetInfo_PFMuon_isGlobal);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_hist", FatJetInfo_PFMuon_hist, &b_FatJetInfo_PFMuon_hist);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_pt", FatJetInfo_PFMuon_pt, &b_FatJetInfo_PFMuon_pt);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_eta", FatJetInfo_PFMuon_eta, &b_FatJetInfo_PFMuon_eta);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_phi", FatJetInfo_PFMuon_phi, &b_FatJetInfo_PFMuon_phi);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_ptrel", FatJetInfo_PFMuon_ptrel, &b_FatJetInfo_PFMuon_ptrel);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_deltaR", FatJetInfo_PFMuon_deltaR, &b_FatJetInfo_PFMuon_deltaR);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_ratio", FatJetInfo_PFMuon_ratio, &b_FatJetInfo_PFMuon_ratio);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_ratioRel", FatJetInfo_PFMuon_ratioRel, &b_FatJetInfo_PFMuon_ratioRel);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IP", FatJetInfo_PFMuon_IP, &b_FatJetInfo_PFMuon_IP);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IP2D", FatJetInfo_PFMuon_IP2D, &b_FatJetInfo_PFMuon_IP2D);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IPsig", FatJetInfo_PFMuon_IPsig, &b_FatJetInfo_PFMuon_IPsig);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_IP2Dsig", FatJetInfo_PFMuon_IP2Dsig, &b_FatJetInfo_PFMuon_IP2Dsig);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_dz", FatJetInfo_PFMuon_dz, &b_FatJetInfo_PFMuon_dz);
   fChain->SetBranchAddress("FatJetInfo.PFMuon_GoodQuality", FatJetInfo_PFMuon_GoodQuality, &b_FatJetInfo_PFMuon_GoodQuality);
   fChain->SetBranchAddress("FatJetInfo.Jet_ptSoftDrop", FatJetInfo_Jet_ptSoftDrop, &b_FatJetInfo_Jet_ptSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_etaSoftDrop", FatJetInfo_Jet_etaSoftDrop, &b_FatJetInfo_Jet_etaSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_phiSoftDrop", FatJetInfo_Jet_phiSoftDrop, &b_FatJetInfo_Jet_phiSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_massSoftDrop", FatJetInfo_Jet_massSoftDrop, &b_FatJetInfo_Jet_massSoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_jecF0SoftDrop", FatJetInfo_Jet_jecF0SoftDrop, &b_FatJetInfo_Jet_jecF0SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_ptPruned", FatJetInfo_Jet_ptPruned, &b_FatJetInfo_Jet_ptPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_etaPruned", FatJetInfo_Jet_etaPruned, &b_FatJetInfo_Jet_etaPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_phiPruned", FatJetInfo_Jet_phiPruned, &b_FatJetInfo_Jet_phiPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_massPruned", FatJetInfo_Jet_massPruned, &b_FatJetInfo_Jet_massPruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_jecF0Pruned", FatJetInfo_Jet_jecF0Pruned, &b_FatJetInfo_Jet_jecF0Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1", FatJetInfo_Jet_tau1, &b_FatJetInfo_Jet_tau1);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2", FatJetInfo_Jet_tau2, &b_FatJetInfo_Jet_tau2);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis1_px", FatJetInfo_Jet_tauAxis1_px, &b_FatJetInfo_Jet_tauAxis1_px);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis1_py", FatJetInfo_Jet_tauAxis1_py, &b_FatJetInfo_Jet_tauAxis1_py);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis1_pz", FatJetInfo_Jet_tauAxis1_pz, &b_FatJetInfo_Jet_tauAxis1_pz);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis2_px", FatJetInfo_Jet_tauAxis2_px, &b_FatJetInfo_Jet_tauAxis2_px);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis2_py", FatJetInfo_Jet_tauAxis2_py, &b_FatJetInfo_Jet_tauAxis2_py);
   fChain->SetBranchAddress("FatJetInfo.Jet_tauAxis2_pz", FatJetInfo_Jet_tauAxis2_pz, &b_FatJetInfo_Jet_tauAxis2_pz);
   fChain->SetBranchAddress("FatJetInfo.Jet_z_ratio", FatJetInfo_Jet_z_ratio, &b_FatJetInfo_Jet_z_ratio);
   fChain->SetBranchAddress("FatJetInfo.Jet_nTracks_fat", FatJetInfo_Jet_nTracks_fat, &b_FatJetInfo_Jet_nTracks_fat);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSV_fat", FatJetInfo_Jet_nSV_fat, &b_FatJetInfo_Jet_nSV_fat);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackEtaRel_2", FatJetInfo_Jet_trackEtaRel_2, &b_FatJetInfo_Jet_trackEtaRel_2);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackEtaRel_1", FatJetInfo_Jet_trackEtaRel_1, &b_FatJetInfo_Jet_trackEtaRel_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_trackEtaRel_0", FatJetInfo_Jet_trackEtaRel_0, &b_FatJetInfo_Jet_trackEtaRel_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_trackEtaRel_0", FatJetInfo_Jet_tau1_trackEtaRel_0, &b_FatJetInfo_Jet_tau1_trackEtaRel_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_trackEtaRel_1", FatJetInfo_Jet_tau1_trackEtaRel_1, &b_FatJetInfo_Jet_tau1_trackEtaRel_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_trackEtaRel_2", FatJetInfo_Jet_tau1_trackEtaRel_2, &b_FatJetInfo_Jet_tau1_trackEtaRel_2);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_trackEtaRel_0", FatJetInfo_Jet_tau2_trackEtaRel_0, &b_FatJetInfo_Jet_tau2_trackEtaRel_0);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_trackEtaRel_1", FatJetInfo_Jet_tau2_trackEtaRel_1, &b_FatJetInfo_Jet_tau2_trackEtaRel_1);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_trackEtaRel_2", FatJetInfo_Jet_tau2_trackEtaRel_2, &b_FatJetInfo_Jet_tau2_trackEtaRel_2);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_nSecondaryVertices", FatJetInfo_Jet_tau1_nSecondaryVertices, &b_FatJetInfo_Jet_tau1_nSecondaryVertices);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_nSecondaryVertices", FatJetInfo_Jet_tau2_nSecondaryVertices, &b_FatJetInfo_Jet_tau2_nSecondaryVertices);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_flightDistance2dSig", FatJetInfo_Jet_tau1_flightDistance2dSig, &b_FatJetInfo_Jet_tau1_flightDistance2dSig);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_flightDistance2dSig", FatJetInfo_Jet_tau2_flightDistance2dSig, &b_FatJetInfo_Jet_tau2_flightDistance2dSig);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexDeltaR", FatJetInfo_Jet_tau1_vertexDeltaR, &b_FatJetInfo_Jet_tau1_vertexDeltaR);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexDeltaR", FatJetInfo_Jet_tau2_vertexDeltaR, &b_FatJetInfo_Jet_tau2_vertexDeltaR);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexEnergyRatio", FatJetInfo_Jet_tau1_vertexEnergyRatio, &b_FatJetInfo_Jet_tau1_vertexEnergyRatio);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexEnergyRatio", FatJetInfo_Jet_tau2_vertexEnergyRatio, &b_FatJetInfo_Jet_tau2_vertexEnergyRatio);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexMass", FatJetInfo_Jet_tau1_vertexMass, &b_FatJetInfo_Jet_tau1_vertexMass);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexMass", FatJetInfo_Jet_tau2_vertexMass, &b_FatJetInfo_Jet_tau2_vertexMass);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexMass_corrected", FatJetInfo_Jet_tau1_vertexMass_corrected, &b_FatJetInfo_Jet_tau1_vertexMass_corrected);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexMass_corrected", FatJetInfo_Jet_tau2_vertexMass_corrected, &b_FatJetInfo_Jet_tau2_vertexMass_corrected);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau1_vertexNTracks", FatJetInfo_Jet_tau1_vertexNTracks, &b_FatJetInfo_Jet_tau1_vertexNTracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_tau2_vertexNTracks", FatJetInfo_Jet_tau2_vertexNTracks, &b_FatJetInfo_Jet_tau2_vertexNTracks);
   fChain->SetBranchAddress("FatJetInfo.Jet_BDTG_SV", FatJetInfo_Jet_BDTG_SV, &b_FatJetInfo_Jet_BDTG_SV);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstTrkInc", FatJetInfo_Jet_nFirstTrkInc, &b_FatJetInfo_Jet_nFirstTrkInc);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastTrkInc", FatJetInfo_Jet_nLastTrkInc, &b_FatJetInfo_Jet_nLastTrkInc);
   fChain->SetBranchAddress("FatJetInfo.nTrkInc", &FatJetInfo_nTrkInc, &b_FatJetInfo_nTrkInc);
   fChain->SetBranchAddress("FatJetInfo.TrkInc_pt", FatJetInfo_TrkInc_pt, &b_FatJetInfo_TrkInc_pt);
   fChain->SetBranchAddress("FatJetInfo.TrkInc_eta", FatJetInfo_TrkInc_eta, &b_FatJetInfo_TrkInc_eta);
   fChain->SetBranchAddress("FatJetInfo.TrkInc_phi", FatJetInfo_TrkInc_phi, &b_FatJetInfo_TrkInc_phi);
   fChain->SetBranchAddress("FatJetInfo.TrkInc_ptrel", FatJetInfo_TrkInc_ptrel, &b_FatJetInfo_TrkInc_ptrel);
   fChain->SetBranchAddress("FatJetInfo.TrkInc_IPsig", FatJetInfo_TrkInc_IPsig, &b_FatJetInfo_TrkInc_IPsig);
   fChain->SetBranchAddress("FatJetInfo.TrkInc_IP", FatJetInfo_TrkInc_IP, &b_FatJetInfo_TrkInc_IP);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSubJets_SoftDrop", FatJetInfo_Jet_nSubJets_SoftDrop, &b_FatJetInfo_Jet_nSubJets_SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSJ_SoftDrop", FatJetInfo_Jet_nFirstSJ_SoftDrop, &b_FatJetInfo_Jet_nFirstSJ_SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSJ_SoftDrop", FatJetInfo_Jet_nLastSJ_SoftDrop, &b_FatJetInfo_Jet_nLastSJ_SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsharedtracks_SoftDrop", FatJetInfo_Jet_nsharedtracks_SoftDrop, &b_FatJetInfo_Jet_nsharedtracks_SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsubjettracks_SoftDrop", FatJetInfo_Jet_nsubjettracks_SoftDrop, &b_FatJetInfo_Jet_nsubjettracks_SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsharedsubjettracks_SoftDrop", FatJetInfo_Jet_nsharedsubjettracks_SoftDrop, &b_FatJetInfo_Jet_nsharedsubjettracks_SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.nSubJet_SoftDrop", &FatJetInfo_nSubJet_SoftDrop, &b_FatJetInfo_nSubJet_SoftDrop);
   fChain->SetBranchAddress("FatJetInfo.SubJetIdx_SoftDrop", FatJetInfo_SubJetIdx_SoftDrop, &b_FatJetInfo_SubJetIdx_SoftDrop);
   fChain->SetBranchAddress("SoftDropSubJetInfo.nJet", &SoftDropSubJetInfo_nJet, &b_SoftDropSubJetInfo_nJet);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_pt", SoftDropSubJetInfo_Jet_pt, &b_SoftDropSubJetInfo_Jet_pt);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_genpt", SoftDropSubJetInfo_Jet_genpt, &b_SoftDropSubJetInfo_Jet_genpt);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_residual", SoftDropSubJetInfo_Jet_residual, &b_SoftDropSubJetInfo_Jet_residual);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_area", SoftDropSubJetInfo_Jet_area, &b_SoftDropSubJetInfo_Jet_area);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_jes", SoftDropSubJetInfo_Jet_jes, &b_SoftDropSubJetInfo_Jet_jes);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_eta", SoftDropSubJetInfo_Jet_eta, &b_SoftDropSubJetInfo_Jet_eta);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_phi", SoftDropSubJetInfo_Jet_phi, &b_SoftDropSubJetInfo_Jet_phi);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_mass", SoftDropSubJetInfo_Jet_mass, &b_SoftDropSubJetInfo_Jet_mass);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_ntracks", SoftDropSubJetInfo_Jet_ntracks, &b_SoftDropSubJetInfo_Jet_ntracks);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nseltracks", SoftDropSubJetInfo_Jet_nseltracks, &b_SoftDropSubJetInfo_Jet_nseltracks);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_flavour", SoftDropSubJetInfo_Jet_flavour, &b_SoftDropSubJetInfo_Jet_flavour);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nbHadrons", SoftDropSubJetInfo_Jet_nbHadrons, &b_SoftDropSubJetInfo_Jet_nbHadrons);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_ncHadrons", SoftDropSubJetInfo_Jet_ncHadrons, &b_SoftDropSubJetInfo_Jet_ncHadrons);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_Ip2N", SoftDropSubJetInfo_Jet_Ip2N, &b_SoftDropSubJetInfo_Jet_Ip2N);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_Ip2P", SoftDropSubJetInfo_Jet_Ip2P, &b_SoftDropSubJetInfo_Jet_Ip2P);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_Ip3N", SoftDropSubJetInfo_Jet_Ip3N, &b_SoftDropSubJetInfo_Jet_Ip3N);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_Ip3P", SoftDropSubJetInfo_Jet_Ip3P, &b_SoftDropSubJetInfo_Jet_Ip3P);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_ProbaN", SoftDropSubJetInfo_Jet_ProbaN, &b_SoftDropSubJetInfo_Jet_ProbaN);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_ProbaP", SoftDropSubJetInfo_Jet_ProbaP, &b_SoftDropSubJetInfo_Jet_ProbaP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_Proba", SoftDropSubJetInfo_Jet_Proba, &b_SoftDropSubJetInfo_Jet_Proba);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_BprobN", SoftDropSubJetInfo_Jet_BprobN, &b_SoftDropSubJetInfo_Jet_BprobN);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_BprobP", SoftDropSubJetInfo_Jet_BprobP, &b_SoftDropSubJetInfo_Jet_BprobP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_Bprob", SoftDropSubJetInfo_Jet_Bprob, &b_SoftDropSubJetInfo_Jet_Bprob);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SvxN", SoftDropSubJetInfo_Jet_SvxN, &b_SoftDropSubJetInfo_Jet_SvxN);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_Svx", SoftDropSubJetInfo_Jet_Svx, &b_SoftDropSubJetInfo_Jet_Svx);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SvxNHP", SoftDropSubJetInfo_Jet_SvxNHP, &b_SoftDropSubJetInfo_Jet_SvxNHP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SvxHP", SoftDropSubJetInfo_Jet_SvxHP, &b_SoftDropSubJetInfo_Jet_SvxHP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_CombSvxN", SoftDropSubJetInfo_Jet_CombSvxN, &b_SoftDropSubJetInfo_Jet_CombSvxN);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_CombSvxP", SoftDropSubJetInfo_Jet_CombSvxP, &b_SoftDropSubJetInfo_Jet_CombSvxP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_CombSvx", SoftDropSubJetInfo_Jet_CombSvx, &b_SoftDropSubJetInfo_Jet_CombSvx);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_CombIVF", SoftDropSubJetInfo_Jet_CombIVF, &b_SoftDropSubJetInfo_Jet_CombIVF);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_CombIVF_P", SoftDropSubJetInfo_Jet_CombIVF_P, &b_SoftDropSubJetInfo_Jet_CombIVF_P);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_CombIVF_N", SoftDropSubJetInfo_Jet_CombIVF_N, &b_SoftDropSubJetInfo_Jet_CombIVF_N);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SoftMuN", SoftDropSubJetInfo_Jet_SoftMuN, &b_SoftDropSubJetInfo_Jet_SoftMuN);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SoftMuP", SoftDropSubJetInfo_Jet_SoftMuP, &b_SoftDropSubJetInfo_Jet_SoftMuP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SoftMu", SoftDropSubJetInfo_Jet_SoftMu, &b_SoftDropSubJetInfo_Jet_SoftMu);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SoftElN", SoftDropSubJetInfo_Jet_SoftElN, &b_SoftDropSubJetInfo_Jet_SoftElN);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SoftElP", SoftDropSubJetInfo_Jet_SoftElP, &b_SoftDropSubJetInfo_Jet_SoftElP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SoftEl", SoftDropSubJetInfo_Jet_SoftEl, &b_SoftDropSubJetInfo_Jet_SoftEl);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_DoubleSV", SoftDropSubJetInfo_Jet_DoubleSV, &b_SoftDropSubJetInfo_Jet_DoubleSV);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_cMVA", SoftDropSubJetInfo_Jet_cMVA, &b_SoftDropSubJetInfo_Jet_cMVA);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_cMVAv2", SoftDropSubJetInfo_Jet_cMVAv2, &b_SoftDropSubJetInfo_Jet_cMVAv2);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_cMVAv2N", SoftDropSubJetInfo_Jet_cMVAv2N, &b_SoftDropSubJetInfo_Jet_cMVAv2N);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_cMVAv2P", SoftDropSubJetInfo_Jet_cMVAv2P, &b_SoftDropSubJetInfo_Jet_cMVAv2P);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_hist1", SoftDropSubJetInfo_Jet_hist1, &b_SoftDropSubJetInfo_Jet_hist1);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_hist2", SoftDropSubJetInfo_Jet_hist2, &b_SoftDropSubJetInfo_Jet_hist2);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_hist3", SoftDropSubJetInfo_Jet_hist3, &b_SoftDropSubJetInfo_Jet_hist3);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_histJet", SoftDropSubJetInfo_Jet_histJet, &b_SoftDropSubJetInfo_Jet_histJet);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_histSvx", SoftDropSubJetInfo_Jet_histSvx, &b_SoftDropSubJetInfo_Jet_histSvx);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_SV_multi", SoftDropSubJetInfo_Jet_SV_multi, &b_SoftDropSubJetInfo_Jet_SV_multi);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nSM", SoftDropSubJetInfo_Jet_nSM, &b_SoftDropSubJetInfo_Jet_nSM);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nSE", SoftDropSubJetInfo_Jet_nSE, &b_SoftDropSubJetInfo_Jet_nSE);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_looseID", SoftDropSubJetInfo_Jet_looseID, &b_SoftDropSubJetInfo_Jet_looseID);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_tightID", SoftDropSubJetInfo_Jet_tightID, &b_SoftDropSubJetInfo_Jet_tightID);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nFirstSE", SoftDropSubJetInfo_Jet_nFirstSE, &b_SoftDropSubJetInfo_Jet_nFirstSE);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nLastSE", SoftDropSubJetInfo_Jet_nLastSE, &b_SoftDropSubJetInfo_Jet_nLastSE);
   fChain->SetBranchAddress("SoftDropSubJetInfo.nPFElectron", &SoftDropSubJetInfo_nPFElectron, &b_SoftDropSubJetInfo_nPFElectron);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_IdxJet", SoftDropSubJetInfo_PFElectron_IdxJet, &b_SoftDropSubJetInfo_PFElectron_IdxJet);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_pt", SoftDropSubJetInfo_PFElectron_pt, &b_SoftDropSubJetInfo_PFElectron_pt);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_eta", SoftDropSubJetInfo_PFElectron_eta, &b_SoftDropSubJetInfo_PFElectron_eta);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_phi", SoftDropSubJetInfo_PFElectron_phi, &b_SoftDropSubJetInfo_PFElectron_phi);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_ptrel", SoftDropSubJetInfo_PFElectron_ptrel, &b_SoftDropSubJetInfo_PFElectron_ptrel);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_deltaR", SoftDropSubJetInfo_PFElectron_deltaR, &b_SoftDropSubJetInfo_PFElectron_deltaR);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_ratio", SoftDropSubJetInfo_PFElectron_ratio, &b_SoftDropSubJetInfo_PFElectron_ratio);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_ratioRel", SoftDropSubJetInfo_PFElectron_ratioRel, &b_SoftDropSubJetInfo_PFElectron_ratioRel);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_IP", SoftDropSubJetInfo_PFElectron_IP, &b_SoftDropSubJetInfo_PFElectron_IP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFElectron_IP2D", SoftDropSubJetInfo_PFElectron_IP2D, &b_SoftDropSubJetInfo_PFElectron_IP2D);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nFirstSM", SoftDropSubJetInfo_Jet_nFirstSM, &b_SoftDropSubJetInfo_Jet_nFirstSM);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nLastSM", SoftDropSubJetInfo_Jet_nLastSM, &b_SoftDropSubJetInfo_Jet_nLastSM);
   fChain->SetBranchAddress("SoftDropSubJetInfo.nPFMuon", &SoftDropSubJetInfo_nPFMuon, &b_SoftDropSubJetInfo_nPFMuon);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_IdxJet", SoftDropSubJetInfo_PFMuon_IdxJet, &b_SoftDropSubJetInfo_PFMuon_IdxJet);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_nMuHit", SoftDropSubJetInfo_PFMuon_nMuHit, &b_SoftDropSubJetInfo_PFMuon_nMuHit);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_nTkHit", SoftDropSubJetInfo_PFMuon_nTkHit, &b_SoftDropSubJetInfo_PFMuon_nTkHit);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_nPixHit", SoftDropSubJetInfo_PFMuon_nPixHit, &b_SoftDropSubJetInfo_PFMuon_nPixHit);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_nOutHit", SoftDropSubJetInfo_PFMuon_nOutHit, &b_SoftDropSubJetInfo_PFMuon_nOutHit);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_nTkLwM", SoftDropSubJetInfo_PFMuon_nTkLwM, &b_SoftDropSubJetInfo_PFMuon_nTkLwM);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_nPixLwM", SoftDropSubJetInfo_PFMuon_nPixLwM, &b_SoftDropSubJetInfo_PFMuon_nPixLwM);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_nMatched", SoftDropSubJetInfo_PFMuon_nMatched, &b_SoftDropSubJetInfo_PFMuon_nMatched);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_chi2", SoftDropSubJetInfo_PFMuon_chi2, &b_SoftDropSubJetInfo_PFMuon_chi2);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_chi2Tk", SoftDropSubJetInfo_PFMuon_chi2Tk, &b_SoftDropSubJetInfo_PFMuon_chi2Tk);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_isGlobal", SoftDropSubJetInfo_PFMuon_isGlobal, &b_SoftDropSubJetInfo_PFMuon_isGlobal);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_hist", SoftDropSubJetInfo_PFMuon_hist, &b_SoftDropSubJetInfo_PFMuon_hist);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_pt", SoftDropSubJetInfo_PFMuon_pt, &b_SoftDropSubJetInfo_PFMuon_pt);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_eta", SoftDropSubJetInfo_PFMuon_eta, &b_SoftDropSubJetInfo_PFMuon_eta);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_phi", SoftDropSubJetInfo_PFMuon_phi, &b_SoftDropSubJetInfo_PFMuon_phi);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_ptrel", SoftDropSubJetInfo_PFMuon_ptrel, &b_SoftDropSubJetInfo_PFMuon_ptrel);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_deltaR", SoftDropSubJetInfo_PFMuon_deltaR, &b_SoftDropSubJetInfo_PFMuon_deltaR);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_ratio", SoftDropSubJetInfo_PFMuon_ratio, &b_SoftDropSubJetInfo_PFMuon_ratio);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_ratioRel", SoftDropSubJetInfo_PFMuon_ratioRel, &b_SoftDropSubJetInfo_PFMuon_ratioRel);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_IP", SoftDropSubJetInfo_PFMuon_IP, &b_SoftDropSubJetInfo_PFMuon_IP);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_IP2D", SoftDropSubJetInfo_PFMuon_IP2D, &b_SoftDropSubJetInfo_PFMuon_IP2D);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_IPsig", SoftDropSubJetInfo_PFMuon_IPsig, &b_SoftDropSubJetInfo_PFMuon_IPsig);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_IP2Dsig", SoftDropSubJetInfo_PFMuon_IP2Dsig, &b_SoftDropSubJetInfo_PFMuon_IP2Dsig);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_dz", SoftDropSubJetInfo_PFMuon_dz, &b_SoftDropSubJetInfo_PFMuon_dz);
   fChain->SetBranchAddress("SoftDropSubJetInfo.PFMuon_GoodQuality", SoftDropSubJetInfo_PFMuon_GoodQuality, &b_SoftDropSubJetInfo_PFMuon_GoodQuality);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_FatJetIdx", SoftDropSubJetInfo_Jet_FatJetIdx, &b_SoftDropSubJetInfo_Jet_FatJetIdx);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nFirstTrkInc", SoftDropSubJetInfo_Jet_nFirstTrkInc, &b_SoftDropSubJetInfo_Jet_nFirstTrkInc);
   fChain->SetBranchAddress("SoftDropSubJetInfo.Jet_nLastTrkInc", SoftDropSubJetInfo_Jet_nLastTrkInc, &b_SoftDropSubJetInfo_Jet_nLastTrkInc);
   fChain->SetBranchAddress("SoftDropSubJetInfo.nTrkInc", &SoftDropSubJetInfo_nTrkInc, &b_SoftDropSubJetInfo_nTrkInc);
   fChain->SetBranchAddress("SoftDropSubJetInfo.TrkInc_pt", SoftDropSubJetInfo_TrkInc_pt, &b_SoftDropSubJetInfo_TrkInc_pt);
   fChain->SetBranchAddress("SoftDropSubJetInfo.TrkInc_eta", SoftDropSubJetInfo_TrkInc_eta, &b_SoftDropSubJetInfo_TrkInc_eta);
   fChain->SetBranchAddress("SoftDropSubJetInfo.TrkInc_phi", SoftDropSubJetInfo_TrkInc_phi, &b_SoftDropSubJetInfo_TrkInc_phi);
   fChain->SetBranchAddress("SoftDropSubJetInfo.TrkInc_ptrel", SoftDropSubJetInfo_TrkInc_ptrel, &b_SoftDropSubJetInfo_TrkInc_ptrel);
   fChain->SetBranchAddress("SoftDropSubJetInfo.TrkInc_IPsig", SoftDropSubJetInfo_TrkInc_IPsig, &b_SoftDropSubJetInfo_TrkInc_IPsig);
   fChain->SetBranchAddress("SoftDropSubJetInfo.TrkInc_IP", SoftDropSubJetInfo_TrkInc_IP, &b_SoftDropSubJetInfo_TrkInc_IP);
   fChain->SetBranchAddress("FatJetInfo.Jet_nSubJets_Pruned", FatJetInfo_Jet_nSubJets_Pruned, &b_FatJetInfo_Jet_nSubJets_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nFirstSJ_Pruned", FatJetInfo_Jet_nFirstSJ_Pruned, &b_FatJetInfo_Jet_nFirstSJ_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nLastSJ_Pruned", FatJetInfo_Jet_nLastSJ_Pruned, &b_FatJetInfo_Jet_nLastSJ_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsharedtracks_Pruned", FatJetInfo_Jet_nsharedtracks_Pruned, &b_FatJetInfo_Jet_nsharedtracks_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsubjettracks_Pruned", FatJetInfo_Jet_nsubjettracks_Pruned, &b_FatJetInfo_Jet_nsubjettracks_Pruned);
   fChain->SetBranchAddress("FatJetInfo.Jet_nsharedsubjettracks_Pruned", FatJetInfo_Jet_nsharedsubjettracks_Pruned, &b_FatJetInfo_Jet_nsharedsubjettracks_Pruned);
   fChain->SetBranchAddress("FatJetInfo.nSubJet_Pruned", &FatJetInfo_nSubJet_Pruned, &b_FatJetInfo_nSubJet_Pruned);
   fChain->SetBranchAddress("FatJetInfo.SubJetIdx_Pruned", FatJetInfo_SubJetIdx_Pruned, &b_FatJetInfo_SubJetIdx_Pruned);
   fChain->SetBranchAddress("PrunedSubJetInfo.nJet", &PrunedSubJetInfo_nJet, &b_PrunedSubJetInfo_nJet);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_pt", PrunedSubJetInfo_Jet_pt, &b_PrunedSubJetInfo_Jet_pt);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_genpt", PrunedSubJetInfo_Jet_genpt, &b_PrunedSubJetInfo_Jet_genpt);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_residual", PrunedSubJetInfo_Jet_residual, &b_PrunedSubJetInfo_Jet_residual);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_area", PrunedSubJetInfo_Jet_area, &b_PrunedSubJetInfo_Jet_area);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_jes", PrunedSubJetInfo_Jet_jes, &b_PrunedSubJetInfo_Jet_jes);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_eta", PrunedSubJetInfo_Jet_eta, &b_PrunedSubJetInfo_Jet_eta);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_phi", PrunedSubJetInfo_Jet_phi, &b_PrunedSubJetInfo_Jet_phi);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_mass", PrunedSubJetInfo_Jet_mass, &b_PrunedSubJetInfo_Jet_mass);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ntracks", PrunedSubJetInfo_Jet_ntracks, &b_PrunedSubJetInfo_Jet_ntracks);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nseltracks", PrunedSubJetInfo_Jet_nseltracks, &b_PrunedSubJetInfo_Jet_nseltracks);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_flavour", PrunedSubJetInfo_Jet_flavour, &b_PrunedSubJetInfo_Jet_flavour);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nbHadrons", PrunedSubJetInfo_Jet_nbHadrons, &b_PrunedSubJetInfo_Jet_nbHadrons);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ncHadrons", PrunedSubJetInfo_Jet_ncHadrons, &b_PrunedSubJetInfo_Jet_ncHadrons);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Ip2N", PrunedSubJetInfo_Jet_Ip2N, &b_PrunedSubJetInfo_Jet_Ip2N);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Ip2P", PrunedSubJetInfo_Jet_Ip2P, &b_PrunedSubJetInfo_Jet_Ip2P);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Ip3N", PrunedSubJetInfo_Jet_Ip3N, &b_PrunedSubJetInfo_Jet_Ip3N);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Ip3P", PrunedSubJetInfo_Jet_Ip3P, &b_PrunedSubJetInfo_Jet_Ip3P);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ProbaN", PrunedSubJetInfo_Jet_ProbaN, &b_PrunedSubJetInfo_Jet_ProbaN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_ProbaP", PrunedSubJetInfo_Jet_ProbaP, &b_PrunedSubJetInfo_Jet_ProbaP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Proba", PrunedSubJetInfo_Jet_Proba, &b_PrunedSubJetInfo_Jet_Proba);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_BprobN", PrunedSubJetInfo_Jet_BprobN, &b_PrunedSubJetInfo_Jet_BprobN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_BprobP", PrunedSubJetInfo_Jet_BprobP, &b_PrunedSubJetInfo_Jet_BprobP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Bprob", PrunedSubJetInfo_Jet_Bprob, &b_PrunedSubJetInfo_Jet_Bprob);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SvxN", PrunedSubJetInfo_Jet_SvxN, &b_PrunedSubJetInfo_Jet_SvxN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_Svx", PrunedSubJetInfo_Jet_Svx, &b_PrunedSubJetInfo_Jet_Svx);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SvxNHP", PrunedSubJetInfo_Jet_SvxNHP, &b_PrunedSubJetInfo_Jet_SvxNHP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SvxHP", PrunedSubJetInfo_Jet_SvxHP, &b_PrunedSubJetInfo_Jet_SvxHP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombSvxN", PrunedSubJetInfo_Jet_CombSvxN, &b_PrunedSubJetInfo_Jet_CombSvxN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombSvxP", PrunedSubJetInfo_Jet_CombSvxP, &b_PrunedSubJetInfo_Jet_CombSvxP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombSvx", PrunedSubJetInfo_Jet_CombSvx, &b_PrunedSubJetInfo_Jet_CombSvx);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombIVF", PrunedSubJetInfo_Jet_CombIVF, &b_PrunedSubJetInfo_Jet_CombIVF);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombIVF_P", PrunedSubJetInfo_Jet_CombIVF_P, &b_PrunedSubJetInfo_Jet_CombIVF_P);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_CombIVF_N", PrunedSubJetInfo_Jet_CombIVF_N, &b_PrunedSubJetInfo_Jet_CombIVF_N);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftMuN", PrunedSubJetInfo_Jet_SoftMuN, &b_PrunedSubJetInfo_Jet_SoftMuN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftMuP", PrunedSubJetInfo_Jet_SoftMuP, &b_PrunedSubJetInfo_Jet_SoftMuP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftMu", PrunedSubJetInfo_Jet_SoftMu, &b_PrunedSubJetInfo_Jet_SoftMu);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftElN", PrunedSubJetInfo_Jet_SoftElN, &b_PrunedSubJetInfo_Jet_SoftElN);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftElP", PrunedSubJetInfo_Jet_SoftElP, &b_PrunedSubJetInfo_Jet_SoftElP);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SoftEl", PrunedSubJetInfo_Jet_SoftEl, &b_PrunedSubJetInfo_Jet_SoftEl);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_DoubleSV", PrunedSubJetInfo_Jet_DoubleSV, &b_PrunedSubJetInfo_Jet_DoubleSV);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_cMVA", PrunedSubJetInfo_Jet_cMVA, &b_PrunedSubJetInfo_Jet_cMVA);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_cMVAv2", PrunedSubJetInfo_Jet_cMVAv2, &b_PrunedSubJetInfo_Jet_cMVAv2);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_cMVAv2N", PrunedSubJetInfo_Jet_cMVAv2N, &b_PrunedSubJetInfo_Jet_cMVAv2N);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_cMVAv2P", PrunedSubJetInfo_Jet_cMVAv2P, &b_PrunedSubJetInfo_Jet_cMVAv2P);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_hist1", PrunedSubJetInfo_Jet_hist1, &b_PrunedSubJetInfo_Jet_hist1);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_hist2", PrunedSubJetInfo_Jet_hist2, &b_PrunedSubJetInfo_Jet_hist2);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_hist3", PrunedSubJetInfo_Jet_hist3, &b_PrunedSubJetInfo_Jet_hist3);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_histJet", PrunedSubJetInfo_Jet_histJet, &b_PrunedSubJetInfo_Jet_histJet);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_histSvx", PrunedSubJetInfo_Jet_histSvx, &b_PrunedSubJetInfo_Jet_histSvx);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_SV_multi", PrunedSubJetInfo_Jet_SV_multi, &b_PrunedSubJetInfo_Jet_SV_multi);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nSM", PrunedSubJetInfo_Jet_nSM, &b_PrunedSubJetInfo_Jet_nSM);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nSE", PrunedSubJetInfo_Jet_nSE, &b_PrunedSubJetInfo_Jet_nSE);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_looseID", PrunedSubJetInfo_Jet_looseID, &b_PrunedSubJetInfo_Jet_looseID);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_tightID", PrunedSubJetInfo_Jet_tightID, &b_PrunedSubJetInfo_Jet_tightID);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nFirstSE", PrunedSubJetInfo_Jet_nFirstSE, &b_PrunedSubJetInfo_Jet_nFirstSE);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nLastSE", PrunedSubJetInfo_Jet_nLastSE, &b_PrunedSubJetInfo_Jet_nLastSE);
   fChain->SetBranchAddress("PrunedSubJetInfo.nPFElectron", &PrunedSubJetInfo_nPFElectron, &b_PrunedSubJetInfo_nPFElectron);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_IdxJet", PrunedSubJetInfo_PFElectron_IdxJet, &b_PrunedSubJetInfo_PFElectron_IdxJet);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_pt", PrunedSubJetInfo_PFElectron_pt, &b_PrunedSubJetInfo_PFElectron_pt);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_eta", PrunedSubJetInfo_PFElectron_eta, &b_PrunedSubJetInfo_PFElectron_eta);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_phi", PrunedSubJetInfo_PFElectron_phi, &b_PrunedSubJetInfo_PFElectron_phi);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_ptrel", PrunedSubJetInfo_PFElectron_ptrel, &b_PrunedSubJetInfo_PFElectron_ptrel);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_deltaR", PrunedSubJetInfo_PFElectron_deltaR, &b_PrunedSubJetInfo_PFElectron_deltaR);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_ratio", PrunedSubJetInfo_PFElectron_ratio, &b_PrunedSubJetInfo_PFElectron_ratio);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_ratioRel", PrunedSubJetInfo_PFElectron_ratioRel, &b_PrunedSubJetInfo_PFElectron_ratioRel);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_IP", PrunedSubJetInfo_PFElectron_IP, &b_PrunedSubJetInfo_PFElectron_IP);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFElectron_IP2D", PrunedSubJetInfo_PFElectron_IP2D, &b_PrunedSubJetInfo_PFElectron_IP2D);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nFirstSM", PrunedSubJetInfo_Jet_nFirstSM, &b_PrunedSubJetInfo_Jet_nFirstSM);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nLastSM", PrunedSubJetInfo_Jet_nLastSM, &b_PrunedSubJetInfo_Jet_nLastSM);
   fChain->SetBranchAddress("PrunedSubJetInfo.nPFMuon", &PrunedSubJetInfo_nPFMuon, &b_PrunedSubJetInfo_nPFMuon);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_IdxJet", PrunedSubJetInfo_PFMuon_IdxJet, &b_PrunedSubJetInfo_PFMuon_IdxJet);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_nMuHit", PrunedSubJetInfo_PFMuon_nMuHit, &b_PrunedSubJetInfo_PFMuon_nMuHit);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_nTkHit", PrunedSubJetInfo_PFMuon_nTkHit, &b_PrunedSubJetInfo_PFMuon_nTkHit);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_nPixHit", PrunedSubJetInfo_PFMuon_nPixHit, &b_PrunedSubJetInfo_PFMuon_nPixHit);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_nOutHit", PrunedSubJetInfo_PFMuon_nOutHit, &b_PrunedSubJetInfo_PFMuon_nOutHit);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_nTkLwM", PrunedSubJetInfo_PFMuon_nTkLwM, &b_PrunedSubJetInfo_PFMuon_nTkLwM);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_nPixLwM", PrunedSubJetInfo_PFMuon_nPixLwM, &b_PrunedSubJetInfo_PFMuon_nPixLwM);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_nMatched", PrunedSubJetInfo_PFMuon_nMatched, &b_PrunedSubJetInfo_PFMuon_nMatched);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_chi2", PrunedSubJetInfo_PFMuon_chi2, &b_PrunedSubJetInfo_PFMuon_chi2);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_chi2Tk", PrunedSubJetInfo_PFMuon_chi2Tk, &b_PrunedSubJetInfo_PFMuon_chi2Tk);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_isGlobal", PrunedSubJetInfo_PFMuon_isGlobal, &b_PrunedSubJetInfo_PFMuon_isGlobal);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_hist", PrunedSubJetInfo_PFMuon_hist, &b_PrunedSubJetInfo_PFMuon_hist);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_pt", PrunedSubJetInfo_PFMuon_pt, &b_PrunedSubJetInfo_PFMuon_pt);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_eta", PrunedSubJetInfo_PFMuon_eta, &b_PrunedSubJetInfo_PFMuon_eta);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_phi", PrunedSubJetInfo_PFMuon_phi, &b_PrunedSubJetInfo_PFMuon_phi);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_ptrel", PrunedSubJetInfo_PFMuon_ptrel, &b_PrunedSubJetInfo_PFMuon_ptrel);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_deltaR", PrunedSubJetInfo_PFMuon_deltaR, &b_PrunedSubJetInfo_PFMuon_deltaR);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_ratio", PrunedSubJetInfo_PFMuon_ratio, &b_PrunedSubJetInfo_PFMuon_ratio);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_ratioRel", PrunedSubJetInfo_PFMuon_ratioRel, &b_PrunedSubJetInfo_PFMuon_ratioRel);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_IP", PrunedSubJetInfo_PFMuon_IP, &b_PrunedSubJetInfo_PFMuon_IP);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_IP2D", PrunedSubJetInfo_PFMuon_IP2D, &b_PrunedSubJetInfo_PFMuon_IP2D);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_IPsig", PrunedSubJetInfo_PFMuon_IPsig, &b_PrunedSubJetInfo_PFMuon_IPsig);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_IP2Dsig", PrunedSubJetInfo_PFMuon_IP2Dsig, &b_PrunedSubJetInfo_PFMuon_IP2Dsig);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_dz", PrunedSubJetInfo_PFMuon_dz, &b_PrunedSubJetInfo_PFMuon_dz);
   fChain->SetBranchAddress("PrunedSubJetInfo.PFMuon_GoodQuality", PrunedSubJetInfo_PFMuon_GoodQuality, &b_PrunedSubJetInfo_PFMuon_GoodQuality);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_FatJetIdx", PrunedSubJetInfo_Jet_FatJetIdx, &b_PrunedSubJetInfo_Jet_FatJetIdx);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nFirstTrkInc", PrunedSubJetInfo_Jet_nFirstTrkInc, &b_PrunedSubJetInfo_Jet_nFirstTrkInc);
   fChain->SetBranchAddress("PrunedSubJetInfo.Jet_nLastTrkInc", PrunedSubJetInfo_Jet_nLastTrkInc, &b_PrunedSubJetInfo_Jet_nLastTrkInc);
   fChain->SetBranchAddress("PrunedSubJetInfo.nTrkInc", &PrunedSubJetInfo_nTrkInc, &b_PrunedSubJetInfo_nTrkInc);
   fChain->SetBranchAddress("PrunedSubJetInfo.TrkInc_pt", PrunedSubJetInfo_TrkInc_pt, &b_PrunedSubJetInfo_TrkInc_pt);
   fChain->SetBranchAddress("PrunedSubJetInfo.TrkInc_eta", PrunedSubJetInfo_TrkInc_eta, &b_PrunedSubJetInfo_TrkInc_eta);
   fChain->SetBranchAddress("PrunedSubJetInfo.TrkInc_phi", PrunedSubJetInfo_TrkInc_phi, &b_PrunedSubJetInfo_TrkInc_phi);
   fChain->SetBranchAddress("PrunedSubJetInfo.TrkInc_ptrel", PrunedSubJetInfo_TrkInc_ptrel, &b_PrunedSubJetInfo_TrkInc_ptrel);
   fChain->SetBranchAddress("PrunedSubJetInfo.TrkInc_IPsig", PrunedSubJetInfo_TrkInc_IPsig, &b_PrunedSubJetInfo_TrkInc_IPsig);
   fChain->SetBranchAddress("PrunedSubJetInfo.TrkInc_IP", PrunedSubJetInfo_TrkInc_IP, &b_PrunedSubJetInfo_TrkInc_IP);
   Notify();
}

Bool_t SubJets::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SubJets::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SubJets::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SubJets_cxx

