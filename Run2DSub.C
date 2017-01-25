/*-------------------------------------------------------------------------------
*
* Parameter "Type" can be
*   -"DATA12A", "DATA12B" ...
*   -"JetTree_Pt-15to30", "JetTree_Pt-30to50"...
*   -"JetTree_Pt-15to30_TP", "JetTree_Pt-30to50_TP"...
*
* Parameter "tTagger" can be "All", "CSV", "CSVJP"...
*   if it is set to all, it will loop on all defined taggers.
*
*-------------------------------------------------------------------------------*/

//#include "SubJetsAna2D.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>
#include "TChain.h"
#include "TFileCollection.h"
#include "TSystem.h"
#include "SubJets.h"
#include "SubJetsAna2D.h"

void Run2DSub(bool isData=false, TString samplePt="All", TString tTrigger="All", TString tTagger="All", bool WeightTracks = false, TString TrigType = "Moriond17", TString period = "9p2invfb", bool runCondor=false)  //period is used to choose PU reweight file
{

 cout << "\n Input setting: " ;
 cout << "\n Data:          " << isData ;
 cout << "\n Sample_pt:     " << samplePt ;
 cout << "\n Trigger:       " << tTrigger ;
 cout << "\n Tagger:        " << tTagger ;
 cout << "\n Weight track:  " << WeightTracks ;
 cout << "\n TrigType:      " << TrigType ;
 cout << "\n Period:        " << period ;
 cout << "\n Run on condor: " << runCondor ;
 

 bool useInputFileList = false ;
 if (samplePt.Contains(".txt")) useInputFileList = true ; 

 gROOT->Reset() ; 

 // Other variables that should be FIX:

  vector<TString> Cuts;
  //Cuts.push_back("JP");
  Cuts.push_back("CSVv2");
  Cuts.push_back("cMVAv2");

// -- General directory where we can find the NTUPLES
  TString ntDir="FileLists_Moriond17"; //in effect when not using file list as input
// -- General directory for the outputs.
  TString oDir = "/uscmst1b_scratch/lpc1/lpctrig/duong/BTagAnalyzer/Output_Moriond17/";
  if (runCondor) oDir = "./" ;

  TString iDir = "MC/";  // First the CMSSW version, then below the full path.
//          iDir       += "QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/";
//  iDir += "QCD_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv1-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/" ;
  iDir += "QCD/" ;
  if (isData) iDir = "Data/JetHT/" ;   

//FIXME
  //TString weightPU_file    = "QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_Data2016B_16June_nPV.weightPU" ; //2.6invfb
  //TString weightPU_file    = "QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_Data2016B_16June_BX_Cont_nPV.weightPU" ;
//  TString weightPthat_file = "QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12.weightPthat";
  TString weightPU_file    = "QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_300_2400_9p2invfb_nPV.weightPU" ; //9.2invfb
  if (period == "7p7invfb") weightPU_file = "QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_300_2400_7p7invfb_nPV.weightPU" ;
  if (period == "7p7to9p2invfb" || period == "Moriond17") weightPU_file = "QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_300_2400_7p7to9p2invfb_nPV.weightPU" ;
  
  TString weightPthat_file = "QCD_TuneCUETP8M1_13TeV_pythia8_RunII_Moriond17.weightPthat" ;
  TString JSONFile = "one";
  
  if (isData) {
    weightPU_file = "one" ;
    weightPthat_file = "one" ;
  }

  vector<TString> PthatRanges;
  if (samplePt == "All" || samplePt == "170to300")  PthatRanges.push_back("/Pt-170to300_all/");
  if (samplePt == "All" || samplePt == "300to470")  PthatRanges.push_back("/Pt-300to470/");
  if (samplePt == "All" || samplePt == "470to600") PthatRanges.push_back("/Pt-470to600/");
  if (samplePt == "All" || samplePt == "600to800")  PthatRanges.push_back("/Pt-600to800_all/");
  if (samplePt == "All" || samplePt == "800to1000") PthatRanges.push_back("/Pt-800to1000_all/");
  if (samplePt == "All" || samplePt == "1000to1400") PthatRanges.push_back("/Pt-1000to1400_all/");
  if (samplePt == "All" || samplePt == "1400to1800") PthatRanges.push_back("/Pt-1400to1800_all/");
  if (samplePt == "All" || samplePt == "1800to2400") PthatRanges.push_back("/Pt-1800to2400_all/");
  if (samplePt == "All" || samplePt == "2400to3200")  PthatRanges.push_back("/Pt-2400to3200_all/");
  //if (samplePt == "All" || samplePt == "3200toInf")  PthatRanges.push_back("/Pt-3200toInf/");

 /////////////for data PthatRanges is used as looping over data set period///////////////
//FIXME  
  if (isData) {
    PthatRanges.clear() ;
    if (samplePt == "All" || samplePt == "Run2016B_23Sep2016-v3") PthatRanges.push_back("/JetHT_Run2016B_23Sep2016-v3/") ;
    if (samplePt == "All" || samplePt == "Run2016C_23Sep2016-v1") PthatRanges.push_back("/JetHT_Run2016C_23Sep2016-v1/") ;
    if (samplePt == "All" || samplePt == "Run2016D_23Sep2016-v1") PthatRanges.push_back("/JetHT_Run2016D_23Sep2016-v1/") ;
    if (samplePt == "All" || samplePt == "Run2016E_23Sep2016-v1") PthatRanges.push_back("/JetHT_Run2016E_23Sep2016-v1/") ;
    if (samplePt == "All" || samplePt == "Run2016F_23Sep2016-v1") PthatRanges.push_back("/JetHT_Run2016F_23Sep2016-v1/") ;
    if (samplePt == "All" || samplePt == "Run2016G_23Sep2016-v1") PthatRanges.push_back("/JetHT_Run2016G_23Sep2016-v1/") ;
    if (samplePt == "All" || samplePt == "Run2016H_PromptReco-v2") PthatRanges.push_back("/JetHT_Run2016H_PromptReco-v2/") ;
    if (samplePt == "All" || samplePt == "Run2016H_PromptReco-v3") PthatRanges.push_back("/JetHT_Run2016H_PromptReco-v3/") ;
  } 

  bool truePU = false;

  float maxCutJetPtMax = 1000000.;
  
  int TrigVal_12[5] =      {  0, 260, 320, 400, 500 }; //500 will not use trigger prescale 
  TString STrigVal_12[5] = { "0", "260", "320", "400", "500"};
  float TrigCut_12[5] =    {  0., 300, 360, 450., 550.}; //not using trigger still use fatJet pT > 400 as minimum cut

  int NTriggers = 5 ;
  int *TrigVal = TrigVal_12;
  TString *STrigVal = STrigVal_12;
  float *TrigCut = TrigCut_12;

//FIXME, no trigger in MC?
  int iTrigMin = 0; // 0 for MC, 1 for DATA
  int iTrigMax = NTriggers;

 // Compile user's analysis class //
 // gROOT->ProcessLine(".L HistogramManager.C+g") ;
 // gROOT->ProcessLine(".L SubJets.C+g") ;
 // gROOT->ProcessLine(".L SubJetsAna2D.C+g") ;
  
  if (gROOT->GetClass("SubJetsAna2D")==0) return;
  if (gROOT->GetClass("SubJets")==0) return;

 
  TChain c("btagana/ttree");
  TChain FJc("btaganaFatJets/ttree");

// //-------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////

// std::vector<int>::iterator it = fifth.begin(); it != fifth.end(); ++it
  //cout << "\n Cuts: " << Cuts.size() ;
  //cout << "\n trig: " << iTrigMin << "  " << iTrigMax ;
  //cout << "\n Pthat: " << PthatRanges.size() ;
  for( int iCut = 0; iCut < Cuts.size(); ++iCut )
  { 
    for( int iTrig = iTrigMin ; iTrig < iTrigMax; ++iTrig )
    {    
      int nPtHat = PthatRanges.size() ;
      if (useInputFileList) nPtHat = 1 ;
      for( int iPthat = 0; iPthat < nPtHat; ++iPthat )
      {
        // Add RootTuple to the chain.
        c.Reset();  // But first reset it otherwise, you accumulate loop after loop.
        FJc.Reset();  // But first reset it otherwise, you accumulate loop after loop.


       //////////////configure input sample//////////////////
        //cout << "\n File is: " << ntDir + iDir + PthatRanges.at(iPthat) + "JetTree_mc_FatJets_Subjets_*.root" ;

        //c.Add( ntDir + iDir + PthatRanges.at(iPthat) + "JetTree_mc_FatJets_Subjets_*.root" );
        //FJc.Add( ntDir + iDir + PthatRanges.at(iPthat) + "JetTree_mc_FatJets_Subjets_*.root" );
        TString fileList = "" ;
        if (!useInputFileList) {
          TString tmp = PthatRanges.at(iPthat) ;
          tmp = tmp.ReplaceAll("/", "") ;
          fileList = ntDir + "/" + iDir + tmp + ".txt" ;
        }
        else { fileList = samplePt ; }
        cout << "\n File list is: " << fileList << endl ;
        TFileCollection fc("fc","list of input root files", fileList) ;
        c.AddFileInfoList((TCollection*)fc.GetList()) ;         
        FJc.AddFileInfoList((TCollection*)fc.GetList()) ;         

        /////////////configure output////////////////////////
        //TString oFileBase = oDir;
        //oFileBase += iDir;  // use the same path as from the input directory
        //oFileBase += "/2DSoftDrop_CutOnEach/TrigType_" + TrigType;  // Add a subdir to differentiate various TrigType
       // oFileBase += PthatRanges.at(iPthat);      // Add subdir for pthat ranges
   

        //gSystem->mkdir(oFileBase,kTRUE);	// Create output Directory if it does not exist

        //cout << "Created dir" << endl;

        //oFileBase += "/JetTree_mc_FatJets_Subjets_";
        //if( weightPU_file.Contains("_Run2015C_") ) oFileBase += "2015C_";
        //if( weightPU_file.Contains("_Run2015D_05Oct2015_v1_SilverNov03_") ) oFileBase += "Run2015D_05Oct2015_v1_SilverNov03_";
        //else if( weightPU_file.Contains("_Run2015D_PromptReco_v4_SilverNov13_") ) oFileBase += "Run2015D_PromptReco_v4_SilverNov13_";
        //else if( weightPU_file.Contains("_Run2015D_PromptReco_v4_Silver_") ) oFileBase += "Run2015D_PromptReco_v4_Silver_";
        //else if( weightPU_file.Contains("_Run2015D_05Oct2015_v1_") ) oFileBase += "Run2015D_05Oct2015_v1_";
        //else if( weightPU_file.Contains("_Run2015D_") ) oFileBase += "2015D_All_";
        //else if( weightPU_file.Contains("_Data2015D_16Dec2015_v1") ) oFileBase += "Data2015D_16Dec2015_v1_";

        TString oFileBase = oDir;
        if (!runCondor) {
          oFileBase += iDir;  // use the same path as from the input directory
          oFileBase += "/2D/TrigType_" + TrigType;  // Add a subdir to differentiate various TrigType
          if (!useInputFileList) oFileBase += PthatRanges.at(iPthat);      // Add subdir for pthat ranges
          else oFileBase += "_" + samplePt ;
          gSystem->mkdir(oFileBase,kTRUE);	// Create output Directory if it does not exist
          cout << "Created dir" << endl;
        }
        if (!isData) {
          oFileBase += "/JetTree_mc_FatJets_Subjets_";
        }
        else {
          oFileBase += "/JetTree_data_FatJets_Subjets_";
        }
        if (runCondor) {
          TString tmp = samplePt ;
          tmp = tmp.ReplaceAll(".txt", "") ;
          oFileBase += tmp + "_" ;
        }
        //FIXME
        if( weightPU_file.Contains("_Run2015C_") ) oFileBase += "2015C_";
        if( weightPU_file.Contains("_Run2015D_05Oct2015_v1_SilverNov03_") ) oFileBase += "Run2015D_05Oct2015_v1_SilverNov03_";
        else if( weightPU_file.Contains("_Run2015D_PromptReco_v4_SilverNov13_") ) oFileBase += "Run2015D_PromptReco_v4_SilverNov13_";
        else if( weightPU_file.Contains("_Run2015D_PromptReco_v4_Silver_") ) oFileBase += "Run2015D_PromptReco_v4_Silver_";
        else if( weightPU_file.Contains("_Run2015D_05Oct2015_v1_") ) oFileBase += "Run2015D_05Oct2015_v1_";
        else if( weightPU_file.Contains("_Run2015D_") ) oFileBase += "2015D_All_";
        else if( weightPU_file.Contains("_Data2015D_16Dec2015_v1") ) oFileBase += "Data2015D_16Dec2015_v1_"; 

        oFileBase += period + "_" ;


        cout << "\n Number of entries of btagana: " << c.GetEntries() ;
        cout << "\n Number of entries of btagana_subjet: " << FJc.GetEntries() ;

        SubJetsAna2D* t = new SubJetsAna2D(&c);
        SubJets* FJt = new SubJets(&FJc);

        TString oFile = oFileBase;
        if( WeightTracks ) oFile += "TW_";	// For TrackWeight
        else oFile += "NW_";			// For NoWeight
        oFile += STrigVal[iTrig] + "_";	// Add the HLTCut in the fileName
        oFile = oFile + Cuts.at(iCut) + ".root";	// Add tagger type
        cout << "Cuts.at(iCut) = " << Cuts.at(iCut) << endl;
        cout << "tTagger = " << tTagger << endl;

        cout << "STrigVal[iTrig] = " << STrigVal[iTrig] << endl;
        cout << "tTrigger = " << tTrigger << endl;
        if( (Cuts.at(iCut) == tTagger || tTagger == "All") &&
            (STrigVal[iTrig] == tTrigger || tTrigger == "All") )
        {
          cout << "Trying to call loop with following parameters " << endl;
          cout << TrigVal[iTrig] << "\t" << TrigCut[iTrig] << endl;
          cout << oFile << endl;
          cout << weightPU_file << "\t" << weightPthat_file << "\t" << truePU << endl;
          //
          //TEMP change to pt cut of 30
          t->Loop(FJt, Cuts.at(iCut),20.,1000.,
             0., TrigVal[iTrig], TrigCut[iTrig], maxCutJetPtMax, oFile,
             weightPU_file, weightPthat_file, JSONFile, truePU, WeightTracks, TrigType, period);
          //
        }
      }
    }
  }
}
