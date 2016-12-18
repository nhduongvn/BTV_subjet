
// This class has been automatically generated on
// Tue Jul 21 11:32:45 2015 by ROOT version 5.34/23
// from TTree ttree/ttree
// found on file: MCTH/QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt50nsRecodebug_MCRUN2_74_V9A/JetTree_Pt-15to30.root
//////////////////////////////////////////////////////////

#ifndef SubJetsAna2D_h
#define SubJetsAna2D_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <SubJets.h>

#define D_nBitTrigger 100
#define D_nPU 100
#define D_ncQuarks 1000
#define D_nbQuarks 1000
#define D_nBHadrons 1000
#define D_nDHadrons 1000
#define D_nDaughters 1000
#define D_nGenlep 1000
#define D_nGenquark 1000
#define D_nGenPruned 1000
#define D_nGenV0 1000
#define D_nJet 1000
#define D_nPFElectron 1000
#define D_nPFMuon 1000
#define D_nTrkInc 1000

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SubJetsAna2D {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nBitTrigger;
   Int_t           BitTrigger[D_nBitTrigger];
   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Float_t         pthat;
   Float_t         mcweight;
   Int_t           BX;
   Int_t           nPV;
   Float_t         PVz;
   Float_t         PVez;
   Float_t         GenPVz;
   Float_t         nPUtrue;
   Int_t           nPU;
   Int_t           PU_bunch[D_nPU];
   Float_t         PU_z[D_nPU];
   Float_t         PU_sumpT_low[D_nPU];
   Float_t         PU_sumpT_high[D_nPU];
   Int_t           PU_ntrks_low[D_nPU];
   Int_t           PU_ntrks_high[D_nPU];
   Int_t           ncQuarks;
   Float_t         cQuark_pT[D_ncQuarks];
   Float_t         cQuark_eta[D_ncQuarks];
   Float_t         cQuark_phi[D_ncQuarks];
   Int_t           cQuark_pdgID[D_ncQuarks];
   Int_t           cQuark_status[D_ncQuarks];
   Int_t           cQuark_fromGSP[D_ncQuarks];
   Int_t           nbQuarks;
   Float_t         bQuark_pT[D_nbQuarks];
   Float_t         bQuark_eta[D_nbQuarks];
   Float_t         bQuark_phi[D_nbQuarks];
   Int_t           bQuark_pdgID[D_nbQuarks];
   Int_t           bQuark_status[D_nbQuarks];
   Int_t           bQuark_fromGSP[D_nbQuarks];
   Int_t           nBHadrons;
   Float_t         BHadron_pT[D_nBHadrons];
   Float_t         BHadron_eta[D_nBHadrons];
   Float_t         BHadron_phi[D_nBHadrons];
   Float_t         BHadron_mass[D_nBHadrons];
   Int_t           BHadron_pdgID[D_nBHadrons];
   Int_t           BHadron_mother[D_nBHadrons];
   Int_t           BHadron_hasBdaughter[D_nBHadrons];
   Float_t         BHadron_SVx[D_nBHadrons];
   Float_t         BHadron_SVy[D_nBHadrons];
   Float_t         BHadron_SVz[D_nBHadrons];
   Int_t           BHadron_nCharged[D_nBHadrons];
   Int_t           BHadron_DHadron1[D_nBHadrons];
   Int_t           BHadron_DHadron2[D_nBHadrons];
   Int_t           nDHadrons;
   Int_t           nDaughters;
   Float_t         DHadron_pT[D_nDHadrons];
   Float_t         DHadron_eta[D_nDHadrons];
   Float_t         DHadron_phi[D_nDHadrons];
   Int_t           DHadron_pdgID[D_nDHadrons];
   Float_t         DHadron_mass[D_nDHadrons];
   Float_t         DHadron_SVx[D_nDHadrons];
   Float_t         DHadron_SVy[D_nDHadrons];
   Float_t         DHadron_SVz[D_nDHadrons];
   Int_t           DHadron_nDaughters[D_nDHadrons];
   Int_t           DHadron_DaughtersPdgID[D_nDaughters];
   Int_t           DHadron_nChargedDaughters[D_nDHadrons];
   Int_t           DHadron_nCharged[D_nDHadrons];
   Int_t           nGenlep;
   Float_t         Genlep_pT[D_nGenlep];
   Float_t         Genlep_eta[D_nGenlep];
   Float_t         Genlep_phi[D_nGenlep];
   Int_t           Genlep_pdgID[D_nGenlep];
   Int_t           Genlep_status[D_nGenlep];
   Int_t           Genlep_mother[D_nGenlep];
   Int_t           nGenquark;
   Float_t         Genquark_pT[D_nGenquark];
   Float_t         Genquark_eta[D_nGenquark];
   Float_t         Genquark_phi[D_nGenquark];
   Int_t           Genquark_pdgID[D_nGenquark];
   Int_t           Genquark_mother[D_nGenquark];
   Int_t           nGenPruned;
   Float_t         GenPruned_pT[D_nGenPruned];
   Float_t         GenPruned_eta[D_nGenPruned];
   Float_t         GenPruned_phi[D_nGenPruned];
   Float_t         GenPruned_mass[D_nGenPruned];
   Int_t           GenPruned_pdgID[D_nGenPruned];
   Int_t           GenPruned_status[D_nGenPruned];
   Int_t           GenPruned_mother[D_nGenPruned];
   Int_t           nGenV0;
   Float_t         GenV0_pT[D_nGenV0];
   Float_t         GenV0_eta[D_nGenV0];
   Float_t         GenV0_phi[D_nGenV0];
   Int_t           GenV0_pdgID[D_nGenV0];
   Float_t         GenV0_SVx[D_nGenV0];
   Float_t         GenV0_SVy[D_nGenV0];
   Float_t         GenV0_SVz[D_nGenV0];
   Int_t           GenV0_nCharged[D_nGenV0];
   Int_t           nJet;
   Float_t         Jet_pt[D_nJet];
   Float_t         Jet_genpt[D_nJet];
   Float_t         Jet_residual[D_nJet];
   Float_t         Jet_area[D_nJet];
   Float_t         Jet_jes[D_nJet];
   Float_t         Jet_eta[D_nJet];
   Float_t         Jet_phi[D_nJet];
   Float_t         Jet_mass[D_nJet];
   Int_t           Jet_ntracks[D_nJet];
   Int_t           Jet_nseltracks[D_nJet];
   Int_t           Jet_flavour[D_nJet];
   Int_t           Jet_nbHadrons[D_nJet];
   Int_t           Jet_ncHadrons[D_nJet];
   Float_t         Jet_Ip2N[D_nJet];
   Float_t         Jet_Ip2P[D_nJet];
   Float_t         Jet_Ip3N[D_nJet];
   Float_t         Jet_Ip3P[D_nJet];
   Float_t         Jet_ProbaN[D_nJet];
   Float_t         Jet_ProbaP[D_nJet];
   Float_t         Jet_Proba[D_nJet];
   Float_t         Jet_BprobN[D_nJet];
   Float_t         Jet_BprobP[D_nJet];
   Float_t         Jet_Bprob[D_nJet];
   Float_t         Jet_SvxN[D_nJet];
   Float_t         Jet_Svx[D_nJet];
   Float_t         Jet_SvxNHP[D_nJet];
   Float_t         Jet_SvxHP[D_nJet];
   Float_t         Jet_CombSvxN[D_nJet];
   Float_t         Jet_CombSvxP[D_nJet];
   Float_t         Jet_CombSvx[D_nJet];
   Float_t         Jet_CombIVF[D_nJet];
   Float_t         Jet_CombIVF_P[D_nJet];
   Float_t         Jet_CombIVF_N[D_nJet];
   Float_t         Jet_SoftMuN[D_nJet];
   Float_t         Jet_SoftMuP[D_nJet];
   Float_t         Jet_SoftMu[D_nJet];
   Float_t         Jet_SoftElN[D_nJet];
   Float_t         Jet_SoftElP[D_nJet];
   Float_t         Jet_SoftEl[D_nJet];
   Float_t         Jet_DoubleSV[D_nJet];
   Float_t         Jet_cMVA[D_nJet];
   Float_t         Jet_cMVAv2[D_nJet];
   Float_t         Jet_cMVAv2N[D_nJet];
   Float_t         Jet_cMVAv2P[D_nJet];
   Int_t           Jet_hist1[D_nJet];
   Int_t           Jet_hist2[D_nJet];
   Int_t           Jet_hist3[D_nJet];
   Int_t           Jet_histJet[D_nJet];
   Int_t           Jet_histSvx[D_nJet];
   Int_t           Jet_SV_multi[D_nJet];
   Int_t           Jet_nSM[D_nJet];
   Int_t           Jet_nSE[D_nJet];
   Int_t           Jet_looseID[D_nJet];
   Int_t           Jet_tightID[D_nJet];
   Int_t           Jet_nFirstSE[D_nJet];
   Int_t           Jet_nLastSE[D_nJet];
   Int_t           nPFElectron;
   Int_t           PFElectron_IdxJet[D_nPFElectron];
   Float_t         PFElectron_pt[D_nPFElectron];
   Float_t         PFElectron_eta[D_nPFElectron];
   Float_t         PFElectron_phi[D_nPFElectron];
   Float_t         PFElectron_ptrel[D_nPFElectron];
   Float_t         PFElectron_deltaR[D_nPFElectron];
   Float_t         PFElectron_ratio[D_nPFElectron];
   Float_t         PFElectron_ratioRel[D_nPFElectron];
   Float_t         PFElectron_IP[D_nPFElectron];
   Float_t         PFElectron_IP2D[D_nPFElectron];
   Int_t           Jet_nFirstSM[D_nJet];
   Int_t           Jet_nLastSM[D_nJet];
   Int_t           nPFMuon;
   Int_t           PFMuon_IdxJet[D_nPFMuon];
   Int_t           PFMuon_nMuHit[D_nPFMuon];
   Int_t           PFMuon_nTkHit[D_nPFMuon];
   Int_t           PFMuon_nPixHit[D_nPFMuon];
   Int_t           PFMuon_nOutHit[D_nPFMuon];
   Int_t           PFMuon_nTkLwM[D_nPFMuon];
   Int_t           PFMuon_nPixLwM[D_nPFMuon];
   Int_t           PFMuon_nMatched[D_nPFMuon];
   Float_t         PFMuon_chi2[D_nPFMuon];
   Float_t         PFMuon_chi2Tk[D_nPFMuon];
   Int_t           PFMuon_isGlobal[D_nPFMuon];
   Int_t           PFMuon_hist[D_nPFMuon];
   Float_t         PFMuon_pt[D_nPFMuon];
   Float_t         PFMuon_eta[D_nPFMuon];
   Float_t         PFMuon_phi[D_nPFMuon];
   Float_t         PFMuon_ptrel[D_nPFMuon];
   Float_t         PFMuon_deltaR[D_nPFMuon];
   Float_t         PFMuon_ratio[D_nPFMuon];
   Float_t         PFMuon_ratioRel[D_nPFMuon];
   Float_t         PFMuon_IP[D_nPFMuon];
   Float_t         PFMuon_IP2D[D_nPFMuon];
   Float_t         PFMuon_IPsig[D_nPFMuon];
   Float_t         PFMuon_IP2Dsig[D_nPFMuon];
   Float_t         PFMuon_dz[D_nPFMuon];
   Int_t           PFMuon_GoodQuality[D_nPFMuon];
   Int_t           Jet_nFirstTrkInc[D_nJet];
   Int_t           Jet_nLastTrkInc[D_nJet];
   Int_t           nTrkInc;
   Float_t         TrkInc_pt[D_nTrkInc];
   Float_t         TrkInc_eta[D_nTrkInc];
   Float_t         TrkInc_phi[D_nTrkInc];
   Float_t         TrkInc_ptrel[D_nTrkInc];
   Float_t         TrkInc_IPsig[D_nTrkInc];
   Float_t         TrkInc_IP[D_nTrkInc];

   // List of branches
   TBranch        *b_nBitTrigger;   //!
   TBranch        *b_BitTrigger;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_mcweight;   //!
   TBranch        *b_BX;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_PVez;   //!
   TBranch        *b_GenPVz;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_PU_bunch;   //!
   TBranch        *b_PU_z;   //!
   TBranch        *b_PU_sumpT_low;   //!
   TBranch        *b_PU_sumpT_high;   //!
   TBranch        *b_PU_ntrks_low;   //!
   TBranch        *b_PU_ntrks_high;   //!
   TBranch        *b_ncQuarks;   //!
   TBranch        *b_cQuark_pT;   //!
   TBranch        *b_cQuark_eta;   //!
   TBranch        *b_cQuark_phi;   //!
   TBranch        *b_cQuark_pdgID;   //!
   TBranch        *b_cQuark_status;   //!
   TBranch        *b_cQuark_fromGSP;   //!
   TBranch        *b_nbQuarks;   //!
   TBranch        *b_bQuark_pT;   //!
   TBranch        *b_bQuark_eta;   //!
   TBranch        *b_bQuark_phi;   //!
   TBranch        *b_bQuark_pdgID;   //!
   TBranch        *b_bQuark_status;   //!
   TBranch        *b_bQuark_fromGSP;   //!
   TBranch        *b_nBHadrons;   //!
   TBranch        *b_BHadron_pT;   //!
   TBranch        *b_BHadron_eta;   //!
   TBranch        *b_BHadron_phi;   //!
   TBranch        *b_BHadron_mass;   //!
   TBranch        *b_BHadron_pdgID;   //!
   TBranch        *b_BHadron_mother;   //!
   TBranch        *b_BHadron_hasBdaughter;   //!
   TBranch        *b_BHadron_SVx;   //!
   TBranch        *b_BHadron_SVy;   //!
   TBranch        *b_BHadron_SVz;   //!
   TBranch        *b_BHadron_nCharged;   //!
   TBranch        *b_BHadron_DHadron1;   //!
   TBranch        *b_BHadron_DHadron2;   //!
   TBranch        *b_nDHadrons;   //!
   TBranch        *b_nDaughters;   //!
   TBranch        *b_DHadron_pT;   //!
   TBranch        *b_DHadron_eta;   //!
   TBranch        *b_DHadron_phi;   //!
   TBranch        *b_DHadron_pdgID;   //!
   TBranch        *b_DHadron_mass;   //!
   TBranch        *b_DHadron_SVx;   //!
   TBranch        *b_DHadron_SVy;   //!
   TBranch        *b_DHadron_SVz;   //!
   TBranch        *b_DHadron_nDaughters;   //!
   TBranch        *b_DHadron_DaughtersPdgID;   //!
   TBranch        *b_DHadron_nChargedDaughters;   //!
   TBranch        *b_DHadron_nCharged;   //!
   TBranch        *b_nGenlep;   //!
   TBranch        *b_Genlep_pT;   //!
   TBranch        *b_Genlep_eta;   //!
   TBranch        *b_Genlep_phi;   //!
   TBranch        *b_Genlep_pdgID;   //!
   TBranch        *b_Genlep_status;   //!
   TBranch        *b_Genlep_mother;   //!
   TBranch        *b_nGenquark;   //!
   TBranch        *b_Genquark_pT;   //!
   TBranch        *b_Genquark_eta;   //!
   TBranch        *b_Genquark_phi;   //!
   TBranch        *b_Genquark_pdgID;   //!
   TBranch        *b_Genquark_mother;   //!
   TBranch        *b_nGenPruned;   //!
   TBranch        *b_GenPruned_pT;   //!
   TBranch        *b_GenPruned_eta;   //!
   TBranch        *b_GenPruned_phi;   //!
   TBranch        *b_GenPruned_mass;   //!
   TBranch        *b_GenPruned_pdgID;   //!
   TBranch        *b_GenPruned_status;   //!
   TBranch        *b_GenPruned_mother;   //!
   TBranch        *b_nGenV0;   //!
   TBranch        *b_GenV0_pT;   //!
   TBranch        *b_GenV0_eta;   //!
   TBranch        *b_GenV0_phi;   //!
   TBranch        *b_GenV0_pdgID;   //!
   TBranch        *b_GenV0_SVx;   //!
   TBranch        *b_GenV0_SVy;   //!
   TBranch        *b_GenV0_SVz;   //!
   TBranch        *b_GenV0_nCharged;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_genpt;   //!
   TBranch        *b_Jet_residual;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_jes;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_ntracks;   //!
   TBranch        *b_Jet_nseltracks;   //!
   TBranch        *b_Jet_flavour;   //!
   TBranch        *b_Jet_nbHadrons;   //!
   TBranch        *b_Jet_ncHadrons;   //!
   TBranch        *b_Jet_Ip2N;   //!
   TBranch        *b_Jet_Ip2P;   //!
   TBranch        *b_Jet_Ip3N;   //!
   TBranch        *b_Jet_Ip3P;   //!
   TBranch        *b_Jet_ProbaN;   //!
   TBranch        *b_Jet_ProbaP;   //!
   TBranch        *b_Jet_Proba;   //!
   TBranch        *b_Jet_BprobN;   //!
   TBranch        *b_Jet_BprobP;   //!
   TBranch        *b_Jet_Bprob;   //!
   TBranch        *b_Jet_SvxN;   //!
   TBranch        *b_Jet_Svx;   //!
   TBranch        *b_Jet_SvxNHP;   //!
   TBranch        *b_Jet_SvxHP;   //!
   TBranch        *b_Jet_CombSvxN;   //!
   TBranch        *b_Jet_CombSvxP;   //!
   TBranch        *b_Jet_CombSvx;   //!
   TBranch        *b_Jet_CombIVF;   //!
   TBranch        *b_Jet_CombIVF_P;   //!
   TBranch        *b_Jet_CombIVF_N;   //!
   TBranch        *b_Jet_SoftMuN;   //!
   TBranch        *b_Jet_SoftMuP;   //!
   TBranch        *b_Jet_SoftMu;   //!
   TBranch        *b_Jet_SoftElN;   //!
   TBranch        *b_Jet_SoftElP;   //!
   TBranch        *b_Jet_SoftEl;   //!
   TBranch        *b_Jet_DoubleSV;   //!
   TBranch        *b_Jet_cMVA;   //!
   TBranch        *b_Jet_cMVAv2;   //!
   TBranch        *b_Jet_cMVAv2N;   //!
   TBranch        *b_Jet_cMVAv2P;   //!
   TBranch        *b_Jet_hist1;   //!
   TBranch        *b_Jet_hist2;   //!
   TBranch        *b_Jet_hist3;   //!
   TBranch        *b_Jet_histJet;   //!
   TBranch        *b_Jet_histSvx;   //!
   TBranch        *b_Jet_SV_multi;   //!
   TBranch        *b_Jet_nSM;   //!
   TBranch        *b_Jet_nSE;   //!
   TBranch        *b_Jet_looseID;   //!
   TBranch        *b_Jet_tightID;   //!
   TBranch        *b_Jet_nFirstSE;   //!
   TBranch        *b_Jet_nLastSE;   //!
   TBranch        *b_nPFElectron;   //!
   TBranch        *b_PFElectron_IdxJet;   //!
   TBranch        *b_PFElectron_pt;   //!
   TBranch        *b_PFElectron_eta;   //!
   TBranch        *b_PFElectron_phi;   //!
   TBranch        *b_PFElectron_ptrel;   //!
   TBranch        *b_PFElectron_deltaR;   //!
   TBranch        *b_PFElectron_ratio;   //!
   TBranch        *b_PFElectron_ratioRel;   //!
   TBranch        *b_PFElectron_IP;   //!
   TBranch        *b_PFElectron_IP2D;   //!
   TBranch        *b_Jet_nFirstSM;   //!
   TBranch        *b_Jet_nLastSM;   //!
   TBranch        *b_nPFMuon;   //!
   TBranch        *b_PFMuon_IdxJet;   //!
   TBranch        *b_PFMuon_nMuHit;   //!
   TBranch        *b_PFMuon_nTkHit;   //!
   TBranch        *b_PFMuon_nPixHit;   //!
   TBranch        *b_PFMuon_nOutHit;   //!
   TBranch        *b_PFMuon_nTkLwM;   //!
   TBranch        *b_PFMuon_nPixLwM;   //!
   TBranch        *b_PFMuon_nMatched;   //!
   TBranch        *b_PFMuon_chi2;   //!
   TBranch        *b_PFMuon_chi2Tk;   //!
   TBranch        *b_PFMuon_isGlobal;   //!
   TBranch        *b_PFMuon_hist;   //!
   TBranch        *b_PFMuon_pt;   //!
   TBranch        *b_PFMuon_eta;   //!
   TBranch        *b_PFMuon_phi;   //!
   TBranch        *b_PFMuon_ptrel;   //!
   TBranch        *b_PFMuon_deltaR;   //!
   TBranch        *b_PFMuon_ratio;   //!
   TBranch        *b_PFMuon_ratioRel;   //!
   TBranch        *b_PFMuon_IP;   //!
   TBranch        *b_PFMuon_IP2D;   //!
   TBranch        *b_PFMuon_IPsig;   //!
   TBranch        *b_PFMuon_IP2Dsig;   //!
   TBranch        *b_PFMuon_dz;   //!
   TBranch        *b_PFMuon_GoodQuality;   //!
   TBranch        *b_Jet_nFirstTrkInc;   //!
   TBranch        *b_Jet_nLastTrkInc;   //!
   TBranch        *b_nTrkInc;   //!
   TBranch        *b_TrkInc_pt;   //!
   TBranch        *b_TrkInc_eta;   //!
   TBranch        *b_TrkInc_phi;   //!
   TBranch        *b_TrkInc_ptrel;   //!
   TBranch        *b_TrkInc_IPsig;   //!
   TBranch        *b_TrkInc_IP;   //!

   SubJetsAna2D(TTree *tree=0);
   virtual ~SubJetsAna2D();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void   Loop(SubJets *FJ,
               TString TagName,
               Float_t  aPtMin = 20.,
               Float_t  aPtMax = 1000.,
               Float_t  aFreeCut = 0.,
               Int_t    aIntCut = 0,
               float    minCutJetPtMax = -1,
               float maxCutJetPtMax = 999999,
               TString  afilename = "output.root",
               TString weightPU_file = "",
               TString weightPthat_file = "",
               TString JSONFile = "",
               bool truePU = true,
               bool WeightTracks = true,
               TString TrigType = "2011", TString period = "All");
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

