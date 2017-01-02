//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb  2 16:43:05 2016 by ROOT version 5.34/18
// from TTree ttree/ttree
// found on file: /opt/sbg/cms/ui1_data1/vanhove/NTUPLES/CMSSW_7_6_3/MC/QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/Pt-80to120/JetTree_mc_FatJets_Subjets_1.root
//////////////////////////////////////////////////////////

#ifndef SubJets_h
#define SubJets_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define Def_FatJetInfo_nJet 1000
#define Def_FatJetInfo_nPFElectron 1000
#define Def_FatJetInfo_nPFMuon 1000
#define Def_FatJetInfo_nTrkInc 1000
#define Def_FatJetInfo_nSubJet_SoftDrop 1000
#define Def_FatJetInfo_nSubJet_Pruned 1000

#define Def_SoftDropSubJetInfo_nJet 1000
#define Def_SoftDropSubJetInfo_nPFElectron 1000
#define Def_SoftDropSubJetInfo_nPFMuon 1000
#define Def_SoftDropSubJetInfo_nTrkInc 1000

#define Def_PrunedSubJetInfo_nJet 1000
#define Def_PrunedSubJetInfo_nPFElectron 1000
#define Def_PrunedSubJetInfo_nPFMuon 1000
#define Def_PrunedSubJetInfo_nTrkInc 1000

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SubJets {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           FatJetInfo_nJet;
   Float_t         FatJetInfo_Jet_pt[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_genpt[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_residual[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_area[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_jes[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_eta[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_phi[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_mass[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_ntracks[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nseltracks[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_flavour[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nbHadrons[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_ncHadrons[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_Ip2N[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_Ip2P[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_Ip3N[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_Ip3P[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_ProbaN[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_ProbaP[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_Proba[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_BprobN[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_BprobP[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_Bprob[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SvxN[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_Svx[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SvxNHP[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SvxHP[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_CombSvxN[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_CombSvxP[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_CombSvx[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_CombIVF[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_CombIVF_P[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_CombIVF_N[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SoftMuN[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SoftMuP[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SoftMu[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SoftElN[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SoftElP[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_SoftEl[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_DoubleSV[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_cMVA[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_cMVAv2[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_cMVAv2N[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_cMVAv2P[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_hist1[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_hist2[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_hist3[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_histJet[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_histSvx[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_SV_multi[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nSM[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nSE[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_looseID[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_tightID[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nFirstSE[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nLastSE[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_nPFElectron;
   Int_t           FatJetInfo_PFElectron_IdxJet[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_pt[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_eta[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_phi[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_ptrel[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_deltaR[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_ratio[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_ratioRel[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_IP[Def_FatJetInfo_nPFElectron];
   Float_t         FatJetInfo_PFElectron_IP2D[Def_FatJetInfo_nPFElectron];
   Int_t           FatJetInfo_Jet_nFirstSM[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nLastSM[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_nPFMuon;
   Int_t           FatJetInfo_PFMuon_IdxJet[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_nMuHit[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_nTkHit[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_nPixHit[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_nOutHit[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_nTkLwM[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_nPixLwM[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_nMatched[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_chi2[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_chi2Tk[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_isGlobal[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_hist[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_pt[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_eta[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_phi[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_ptrel[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_deltaR[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_ratio[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_ratioRel[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_IP[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_IP2D[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_IPsig[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_IP2Dsig[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_PFMuon_dz[Def_FatJetInfo_nPFMuon];
   Int_t           FatJetInfo_PFMuon_GoodQuality[Def_FatJetInfo_nPFMuon];
   Float_t         FatJetInfo_Jet_ptSoftDrop[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_etaSoftDrop[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_phiSoftDrop[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_massSoftDrop[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_jecF0SoftDrop[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_ptPruned[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_etaPruned[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_phiPruned[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_massPruned[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_jecF0Pruned[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tauAxis1_px[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tauAxis1_py[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tauAxis1_pz[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tauAxis2_px[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tauAxis2_py[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tauAxis2_pz[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_z_ratio[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_nTracks_fat[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_nSV_fat[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_trackEtaRel_2[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_trackEtaRel_1[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_trackEtaRel_0[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_trackEtaRel_0[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_trackEtaRel_1[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_trackEtaRel_2[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_trackEtaRel_0[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_trackEtaRel_1[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_trackEtaRel_2[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_nSecondaryVertices[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_nSecondaryVertices[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_flightDistance2dSig[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_flightDistance2dSig[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_vertexDeltaR[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_vertexDeltaR[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_vertexEnergyRatio[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_vertexEnergyRatio[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_vertexMass[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_vertexMass[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_vertexMass_corrected[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_vertexMass_corrected[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau1_vertexNTracks[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_tau2_vertexNTracks[Def_FatJetInfo_nJet];
   Float_t         FatJetInfo_Jet_BDTG_SV[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nFirstTrkInc[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nLastTrkInc[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_nTrkInc;
   Float_t         FatJetInfo_TrkInc_pt[Def_FatJetInfo_nTrkInc];
   Float_t         FatJetInfo_TrkInc_eta[Def_FatJetInfo_nTrkInc];
   Float_t         FatJetInfo_TrkInc_phi[Def_FatJetInfo_nTrkInc];
   Float_t         FatJetInfo_TrkInc_ptrel[Def_FatJetInfo_nTrkInc];
   Float_t         FatJetInfo_TrkInc_IPsig[Def_FatJetInfo_nTrkInc];
   Float_t         FatJetInfo_TrkInc_IP[Def_FatJetInfo_nTrkInc];
   Int_t           FatJetInfo_Jet_nSubJets_SoftDrop[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nFirstSJ_SoftDrop[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nLastSJ_SoftDrop[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nsharedtracks_SoftDrop[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nsubjettracks_SoftDrop[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nsharedsubjettracks_SoftDrop[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_nSubJet_SoftDrop;
   Int_t           FatJetInfo_SubJetIdx_SoftDrop[Def_FatJetInfo_nSubJet_SoftDrop];
   Int_t           SoftDropSubJetInfo_nJet;
   Float_t         SoftDropSubJetInfo_Jet_pt[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_genpt[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_residual[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_area[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_jes[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_eta[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_phi[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_mass[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_ntracks[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nseltracks[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_flavour[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nbHadrons[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_ncHadrons[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_Ip2N[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_Ip2P[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_Ip3N[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_Ip3P[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_ProbaN[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_ProbaP[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_Proba[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_BprobN[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_BprobP[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_Bprob[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SvxN[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_Svx[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SvxNHP[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SvxHP[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_CombSvxN[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_CombSvxP[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_CombSvx[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_CombIVF[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_CombIVF_P[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_CombIVF_N[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SoftMuN[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SoftMuP[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SoftMu[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SoftElN[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SoftElP[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_SoftEl[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_DoubleSV[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_cMVA[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_cMVAv2[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_cMVAv2N[Def_SoftDropSubJetInfo_nJet];
   Float_t         SoftDropSubJetInfo_Jet_cMVAv2P[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_hist1[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_hist2[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_hist3[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_histJet[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_histSvx[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_SV_multi[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nSM[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nSE[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_looseID[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_tightID[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nFirstSE[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nLastSE[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_nPFElectron;
   Int_t           SoftDropSubJetInfo_PFElectron_IdxJet[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_pt[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_eta[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_phi[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_ptrel[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_deltaR[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_ratio[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_ratioRel[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_IP[Def_SoftDropSubJetInfo_nPFElectron];
   Float_t         SoftDropSubJetInfo_PFElectron_IP2D[Def_SoftDropSubJetInfo_nPFElectron];
   Int_t           SoftDropSubJetInfo_Jet_nFirstSM[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nLastSM[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_nPFMuon;
   Int_t           SoftDropSubJetInfo_PFMuon_IdxJet[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_nMuHit[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_nTkHit[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_nPixHit[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_nOutHit[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_nTkLwM[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_nPixLwM[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_nMatched[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_chi2[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_chi2Tk[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_isGlobal[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_hist[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_pt[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_eta[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_phi[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_ptrel[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_deltaR[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_ratio[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_ratioRel[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_IP[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_IP2D[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_IPsig[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_IP2Dsig[Def_SoftDropSubJetInfo_nPFMuon];
   Float_t         SoftDropSubJetInfo_PFMuon_dz[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_PFMuon_GoodQuality[Def_SoftDropSubJetInfo_nPFMuon];
   Int_t           SoftDropSubJetInfo_Jet_FatJetIdx[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nFirstTrkInc[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_Jet_nLastTrkInc[Def_SoftDropSubJetInfo_nJet];
   Int_t           SoftDropSubJetInfo_nTrkInc;
   Float_t         SoftDropSubJetInfo_TrkInc_pt[Def_SoftDropSubJetInfo_nTrkInc];
   Float_t         SoftDropSubJetInfo_TrkInc_eta[Def_SoftDropSubJetInfo_nTrkInc];
   Float_t         SoftDropSubJetInfo_TrkInc_phi[Def_SoftDropSubJetInfo_nTrkInc];
   Float_t         SoftDropSubJetInfo_TrkInc_ptrel[Def_SoftDropSubJetInfo_nTrkInc];
   Float_t         SoftDropSubJetInfo_TrkInc_IPsig[Def_SoftDropSubJetInfo_nTrkInc];
   Float_t         SoftDropSubJetInfo_TrkInc_IP[Def_SoftDropSubJetInfo_nTrkInc];
   Int_t           FatJetInfo_Jet_nSubJets_Pruned[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nFirstSJ_Pruned[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nLastSJ_Pruned[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nsharedtracks_Pruned[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nsubjettracks_Pruned[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_Jet_nsharedsubjettracks_Pruned[Def_FatJetInfo_nJet];
   Int_t           FatJetInfo_nSubJet_Pruned;
   Int_t           FatJetInfo_SubJetIdx_Pruned[Def_FatJetInfo_nSubJet_Pruned];
   Int_t           PrunedSubJetInfo_nJet;
   Float_t         PrunedSubJetInfo_Jet_pt[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_genpt[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_residual[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_area[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_jes[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_eta[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_phi[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_mass[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_ntracks[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nseltracks[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_flavour[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nbHadrons[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_ncHadrons[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_Ip2N[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_Ip2P[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_Ip3N[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_Ip3P[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_ProbaN[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_ProbaP[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_Proba[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_BprobN[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_BprobP[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_Bprob[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SvxN[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_Svx[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SvxNHP[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SvxHP[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_CombSvxN[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_CombSvxP[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_CombSvx[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_CombIVF[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_CombIVF_P[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_CombIVF_N[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SoftMuN[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SoftMuP[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SoftMu[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SoftElN[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SoftElP[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_SoftEl[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_DoubleSV[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_cMVA[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_cMVAv2[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_cMVAv2N[Def_PrunedSubJetInfo_nJet];
   Float_t         PrunedSubJetInfo_Jet_cMVAv2P[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_hist1[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_hist2[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_hist3[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_histJet[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_histSvx[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_SV_multi[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nSM[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nSE[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_looseID[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_tightID[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nFirstSE[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nLastSE[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_nPFElectron;
   Int_t           PrunedSubJetInfo_PFElectron_IdxJet[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_pt[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_eta[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_phi[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_ptrel[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_deltaR[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_ratio[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_ratioRel[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_IP[Def_PrunedSubJetInfo_nPFElectron];
   Float_t         PrunedSubJetInfo_PFElectron_IP2D[Def_PrunedSubJetInfo_nPFElectron];
   Int_t           PrunedSubJetInfo_Jet_nFirstSM[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nLastSM[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_nPFMuon;
   Int_t           PrunedSubJetInfo_PFMuon_IdxJet[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_nMuHit[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_nTkHit[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_nPixHit[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_nOutHit[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_nTkLwM[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_nPixLwM[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_nMatched[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_chi2[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_chi2Tk[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_isGlobal[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_hist[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_pt[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_eta[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_phi[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_ptrel[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_deltaR[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_ratio[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_ratioRel[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_IP[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_IP2D[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_IPsig[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_IP2Dsig[Def_PrunedSubJetInfo_nPFMuon];
   Float_t         PrunedSubJetInfo_PFMuon_dz[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_PFMuon_GoodQuality[Def_PrunedSubJetInfo_nPFMuon];
   Int_t           PrunedSubJetInfo_Jet_FatJetIdx[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nFirstTrkInc[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_Jet_nLastTrkInc[Def_PrunedSubJetInfo_nJet];
   Int_t           PrunedSubJetInfo_nTrkInc;
   Float_t         PrunedSubJetInfo_TrkInc_pt[Def_PrunedSubJetInfo_nTrkInc];
   Float_t         PrunedSubJetInfo_TrkInc_eta[Def_PrunedSubJetInfo_nTrkInc];
   Float_t         PrunedSubJetInfo_TrkInc_phi[Def_PrunedSubJetInfo_nTrkInc];
   Float_t         PrunedSubJetInfo_TrkInc_ptrel[Def_PrunedSubJetInfo_nTrkInc];
   Float_t         PrunedSubJetInfo_TrkInc_IPsig[Def_PrunedSubJetInfo_nTrkInc];
   Float_t         PrunedSubJetInfo_TrkInc_IP[Def_PrunedSubJetInfo_nTrkInc];

   // List of branches
   TBranch        *b_FatJetInfo_nJet;   //!
   TBranch        *b_FatJetInfo_Jet_pt;   //!
   TBranch        *b_FatJetInfo_Jet_genpt;   //!
   TBranch        *b_FatJetInfo_Jet_residual;   //!
   TBranch        *b_FatJetInfo_Jet_area;   //!
   TBranch        *b_FatJetInfo_Jet_jes;   //!
   TBranch        *b_FatJetInfo_Jet_eta;   //!
   TBranch        *b_FatJetInfo_Jet_phi;   //!
   TBranch        *b_FatJetInfo_Jet_mass;   //!
   TBranch        *b_FatJetInfo_Jet_ntracks;   //!
   TBranch        *b_FatJetInfo_Jet_nseltracks;   //!
   TBranch        *b_FatJetInfo_Jet_flavour;   //!
   TBranch        *b_FatJetInfo_Jet_nbHadrons;   //!
   TBranch        *b_FatJetInfo_Jet_ncHadrons;   //!
   TBranch        *b_FatJetInfo_Jet_Ip2N;   //!
   TBranch        *b_FatJetInfo_Jet_Ip2P;   //!
   TBranch        *b_FatJetInfo_Jet_Ip3N;   //!
   TBranch        *b_FatJetInfo_Jet_Ip3P;   //!
   TBranch        *b_FatJetInfo_Jet_ProbaN;   //!
   TBranch        *b_FatJetInfo_Jet_ProbaP;   //!
   TBranch        *b_FatJetInfo_Jet_Proba;   //!
   TBranch        *b_FatJetInfo_Jet_BprobN;   //!
   TBranch        *b_FatJetInfo_Jet_BprobP;   //!
   TBranch        *b_FatJetInfo_Jet_Bprob;   //!
   TBranch        *b_FatJetInfo_Jet_SvxN;   //!
   TBranch        *b_FatJetInfo_Jet_Svx;   //!
   TBranch        *b_FatJetInfo_Jet_SvxNHP;   //!
   TBranch        *b_FatJetInfo_Jet_SvxHP;   //!
   TBranch        *b_FatJetInfo_Jet_CombSvxN;   //!
   TBranch        *b_FatJetInfo_Jet_CombSvxP;   //!
   TBranch        *b_FatJetInfo_Jet_CombSvx;   //!
   TBranch        *b_FatJetInfo_Jet_CombIVF;   //!
   TBranch        *b_FatJetInfo_Jet_CombIVF_P;   //!
   TBranch        *b_FatJetInfo_Jet_CombIVF_N;   //!
   TBranch        *b_FatJetInfo_Jet_SoftMuN;   //!
   TBranch        *b_FatJetInfo_Jet_SoftMuP;   //!
   TBranch        *b_FatJetInfo_Jet_SoftMu;   //!
   TBranch        *b_FatJetInfo_Jet_SoftElN;   //!
   TBranch        *b_FatJetInfo_Jet_SoftElP;   //!
   TBranch        *b_FatJetInfo_Jet_SoftEl;   //!
   TBranch        *b_FatJetInfo_Jet_DoubleSV;   //!
   TBranch        *b_FatJetInfo_Jet_cMVA;   //!
   TBranch        *b_FatJetInfo_Jet_cMVAv2;   //!
   TBranch        *b_FatJetInfo_Jet_cMVAv2N;   //!
   TBranch        *b_FatJetInfo_Jet_cMVAv2P;   //!
   TBranch        *b_FatJetInfo_Jet_hist1;   //!
   TBranch        *b_FatJetInfo_Jet_hist2;   //!
   TBranch        *b_FatJetInfo_Jet_hist3;   //!
   TBranch        *b_FatJetInfo_Jet_histJet;   //!
   TBranch        *b_FatJetInfo_Jet_histSvx;   //!
   TBranch        *b_FatJetInfo_Jet_SV_multi;   //!
   TBranch        *b_FatJetInfo_Jet_nSM;   //!
   TBranch        *b_FatJetInfo_Jet_nSE;   //!
   TBranch        *b_FatJetInfo_Jet_looseID;   //!
   TBranch        *b_FatJetInfo_Jet_tightID;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSE;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSE;   //!
   TBranch        *b_FatJetInfo_nPFElectron;   //!
   TBranch        *b_FatJetInfo_PFElectron_IdxJet;   //!
   TBranch        *b_FatJetInfo_PFElectron_pt;   //!
   TBranch        *b_FatJetInfo_PFElectron_eta;   //!
   TBranch        *b_FatJetInfo_PFElectron_phi;   //!
   TBranch        *b_FatJetInfo_PFElectron_ptrel;   //!
   TBranch        *b_FatJetInfo_PFElectron_deltaR;   //!
   TBranch        *b_FatJetInfo_PFElectron_ratio;   //!
   TBranch        *b_FatJetInfo_PFElectron_ratioRel;   //!
   TBranch        *b_FatJetInfo_PFElectron_IP;   //!
   TBranch        *b_FatJetInfo_PFElectron_IP2D;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSM;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSM;   //!
   TBranch        *b_FatJetInfo_nPFMuon;   //!
   TBranch        *b_FatJetInfo_PFMuon_IdxJet;   //!
   TBranch        *b_FatJetInfo_PFMuon_nMuHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nTkHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nPixHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nOutHit;   //!
   TBranch        *b_FatJetInfo_PFMuon_nTkLwM;   //!
   TBranch        *b_FatJetInfo_PFMuon_nPixLwM;   //!
   TBranch        *b_FatJetInfo_PFMuon_nMatched;   //!
   TBranch        *b_FatJetInfo_PFMuon_chi2;   //!
   TBranch        *b_FatJetInfo_PFMuon_chi2Tk;   //!
   TBranch        *b_FatJetInfo_PFMuon_isGlobal;   //!
   TBranch        *b_FatJetInfo_PFMuon_hist;   //!
   TBranch        *b_FatJetInfo_PFMuon_pt;   //!
   TBranch        *b_FatJetInfo_PFMuon_eta;   //!
   TBranch        *b_FatJetInfo_PFMuon_phi;   //!
   TBranch        *b_FatJetInfo_PFMuon_ptrel;   //!
   TBranch        *b_FatJetInfo_PFMuon_deltaR;   //!
   TBranch        *b_FatJetInfo_PFMuon_ratio;   //!
   TBranch        *b_FatJetInfo_PFMuon_ratioRel;   //!
   TBranch        *b_FatJetInfo_PFMuon_IP;   //!
   TBranch        *b_FatJetInfo_PFMuon_IP2D;   //!
   TBranch        *b_FatJetInfo_PFMuon_IPsig;   //!
   TBranch        *b_FatJetInfo_PFMuon_IP2Dsig;   //!
   TBranch        *b_FatJetInfo_PFMuon_dz;   //!
   TBranch        *b_FatJetInfo_PFMuon_GoodQuality;   //!
   TBranch        *b_FatJetInfo_Jet_ptSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_etaSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_phiSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_massSoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_jecF0SoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_ptPruned;   //!
   TBranch        *b_FatJetInfo_Jet_etaPruned;   //!
   TBranch        *b_FatJetInfo_Jet_phiPruned;   //!
   TBranch        *b_FatJetInfo_Jet_massPruned;   //!
   TBranch        *b_FatJetInfo_Jet_jecF0Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_tau1;   //!
   TBranch        *b_FatJetInfo_Jet_tau2;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis1_px;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis1_py;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis1_pz;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis2_px;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis2_py;   //!
   TBranch        *b_FatJetInfo_Jet_tauAxis2_pz;   //!
   TBranch        *b_FatJetInfo_Jet_z_ratio;   //!
   TBranch        *b_FatJetInfo_Jet_nTracks_fat;   //!
   TBranch        *b_FatJetInfo_Jet_nSV_fat;   //!
   TBranch        *b_FatJetInfo_Jet_trackEtaRel_2;   //!
   TBranch        *b_FatJetInfo_Jet_trackEtaRel_1;   //!
   TBranch        *b_FatJetInfo_Jet_trackEtaRel_0;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_trackEtaRel_0;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_trackEtaRel_1;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_trackEtaRel_2;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_trackEtaRel_0;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_trackEtaRel_1;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_trackEtaRel_2;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_nSecondaryVertices;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_nSecondaryVertices;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_flightDistance2dSig;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_flightDistance2dSig;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexDeltaR;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexDeltaR;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexEnergyRatio;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexEnergyRatio;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexMass;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexMass;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexMass_corrected;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexMass_corrected;   //!
   TBranch        *b_FatJetInfo_Jet_tau1_vertexNTracks;   //!
   TBranch        *b_FatJetInfo_Jet_tau2_vertexNTracks;   //!
   TBranch        *b_FatJetInfo_Jet_BDTG_SV;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstTrkInc;   //!
   TBranch        *b_FatJetInfo_Jet_nLastTrkInc;   //!
   TBranch        *b_FatJetInfo_nTrkInc;   //!
   TBranch        *b_FatJetInfo_TrkInc_pt;   //!
   TBranch        *b_FatJetInfo_TrkInc_eta;   //!
   TBranch        *b_FatJetInfo_TrkInc_phi;   //!
   TBranch        *b_FatJetInfo_TrkInc_ptrel;   //!
   TBranch        *b_FatJetInfo_TrkInc_IPsig;   //!
   TBranch        *b_FatJetInfo_TrkInc_IP;   //!
   TBranch        *b_FatJetInfo_Jet_nSubJets_SoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSJ_SoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSJ_SoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_nsharedtracks_SoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_nsubjettracks_SoftDrop;   //!
   TBranch        *b_FatJetInfo_Jet_nsharedsubjettracks_SoftDrop;   //!
   TBranch        *b_FatJetInfo_nSubJet_SoftDrop;   //!
   TBranch        *b_FatJetInfo_SubJetIdx_SoftDrop;   //!
   TBranch        *b_SoftDropSubJetInfo_nJet;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_pt;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_genpt;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_residual;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_area;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_jes;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_eta;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_phi;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_mass;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_ntracks;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nseltracks;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_flavour;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nbHadrons;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_ncHadrons;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_Ip2N;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_Ip2P;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_Ip3N;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_Ip3P;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_ProbaN;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_ProbaP;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_Proba;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_BprobN;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_BprobP;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_Bprob;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SvxN;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_Svx;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SvxNHP;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SvxHP;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_CombSvxN;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_CombSvxP;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_CombSvx;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_CombIVF;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_CombIVF_P;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_CombIVF_N;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SoftMuN;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SoftMuP;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SoftMu;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SoftElN;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SoftElP;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SoftEl;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_DoubleSV;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_cMVA;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_cMVAv2;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_cMVAv2N;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_cMVAv2P;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_hist1;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_hist2;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_hist3;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_histJet;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_histSvx;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_SV_multi;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nSM;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nSE;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_looseID;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_tightID;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nFirstSE;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nLastSE;   //!
   TBranch        *b_SoftDropSubJetInfo_nPFElectron;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_IdxJet;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_pt;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_eta;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_phi;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_ptrel;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_deltaR;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_ratio;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_ratioRel;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_IP;   //!
   TBranch        *b_SoftDropSubJetInfo_PFElectron_IP2D;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nFirstSM;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nLastSM;   //!
   TBranch        *b_SoftDropSubJetInfo_nPFMuon;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_IdxJet;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_nMuHit;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_nTkHit;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_nPixHit;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_nOutHit;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_nTkLwM;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_nPixLwM;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_nMatched;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_chi2;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_chi2Tk;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_isGlobal;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_hist;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_pt;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_eta;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_phi;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_ptrel;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_deltaR;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_ratio;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_ratioRel;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_IP;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_IP2D;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_IPsig;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_IP2Dsig;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_dz;   //!
   TBranch        *b_SoftDropSubJetInfo_PFMuon_GoodQuality;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_FatJetIdx;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nFirstTrkInc;   //!
   TBranch        *b_SoftDropSubJetInfo_Jet_nLastTrkInc;   //!
   TBranch        *b_SoftDropSubJetInfo_nTrkInc;   //!
   TBranch        *b_SoftDropSubJetInfo_TrkInc_pt;   //!
   TBranch        *b_SoftDropSubJetInfo_TrkInc_eta;   //!
   TBranch        *b_SoftDropSubJetInfo_TrkInc_phi;   //!
   TBranch        *b_SoftDropSubJetInfo_TrkInc_ptrel;   //!
   TBranch        *b_SoftDropSubJetInfo_TrkInc_IPsig;   //!
   TBranch        *b_SoftDropSubJetInfo_TrkInc_IP;   //!
   TBranch        *b_FatJetInfo_Jet_nSubJets_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nFirstSJ_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nLastSJ_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nsharedtracks_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nsubjettracks_Pruned;   //!
   TBranch        *b_FatJetInfo_Jet_nsharedsubjettracks_Pruned;   //!
   TBranch        *b_FatJetInfo_nSubJet_Pruned;   //!
   TBranch        *b_FatJetInfo_SubJetIdx_Pruned;   //!
   TBranch        *b_PrunedSubJetInfo_nJet;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_pt;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_genpt;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_residual;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_area;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_jes;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_eta;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_phi;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_mass;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ntracks;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nseltracks;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_flavour;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nbHadrons;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ncHadrons;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Ip2N;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Ip2P;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Ip3N;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Ip3P;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ProbaN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_ProbaP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Proba;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_BprobN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_BprobP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Bprob;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SvxN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_Svx;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SvxNHP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SvxHP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombSvxN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombSvxP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombSvx;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombIVF;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombIVF_P;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_CombIVF_N;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftMuN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftMuP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftMu;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftElN;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftElP;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SoftEl;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_DoubleSV;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_cMVA;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_cMVAv2;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_cMVAv2N;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_cMVAv2P;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_hist1;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_hist2;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_hist3;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_histJet;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_histSvx;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_SV_multi;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nSM;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nSE;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_looseID;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_tightID;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nFirstSE;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nLastSE;   //!
   TBranch        *b_PrunedSubJetInfo_nPFElectron;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_IdxJet;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_pt;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_eta;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_phi;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_ptrel;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_deltaR;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_ratio;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_ratioRel;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_IP;   //!
   TBranch        *b_PrunedSubJetInfo_PFElectron_IP2D;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nFirstSM;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nLastSM;   //!
   TBranch        *b_PrunedSubJetInfo_nPFMuon;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_IdxJet;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_nMuHit;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_nTkHit;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_nPixHit;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_nOutHit;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_nTkLwM;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_nPixLwM;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_nMatched;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_chi2;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_chi2Tk;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_isGlobal;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_hist;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_pt;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_eta;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_phi;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_ptrel;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_deltaR;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_ratio;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_ratioRel;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_IP;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_IP2D;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_IPsig;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_IP2Dsig;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_dz;   //!
   TBranch        *b_PrunedSubJetInfo_PFMuon_GoodQuality;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_FatJetIdx;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nFirstTrkInc;   //!
   TBranch        *b_PrunedSubJetInfo_Jet_nLastTrkInc;   //!
   TBranch        *b_PrunedSubJetInfo_nTrkInc;   //!
   TBranch        *b_PrunedSubJetInfo_TrkInc_pt;   //!
   TBranch        *b_PrunedSubJetInfo_TrkInc_eta;   //!
   TBranch        *b_PrunedSubJetInfo_TrkInc_phi;   //!
   TBranch        *b_PrunedSubJetInfo_TrkInc_ptrel;   //!
   TBranch        *b_PrunedSubJetInfo_TrkInc_IPsig;   //!
   TBranch        *b_PrunedSubJetInfo_TrkInc_IP;   //!

   SubJets(TTree *tree=0);
   virtual ~SubJets();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

