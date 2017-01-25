#ifndef SubJetsAna2D_cxx
#define SubJetsAna2D_cxx
#include "SubJetsAna2D.h"
#include "PS.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "HistogramManager.h"
#include "TLorentzVector.h"
#include <iostream>

#include <fstream>
#include <TString.h>

//TEMP change pt low pT bin to 30
#define minPtBin  20.
#define maxPtBin  1000.
#define NPtBin  ( (int) ((maxPtBin-minPtBin)/10.) )

#define minEtaBin  0.
#define maxEtaBin  2.4
#define NEtaBin  ( (int) ((maxEtaBin-minEtaBin)*10.) )

#define NDelimiters 7

int ReadJSon(unsigned int N_MaxRun, unsigned int N_MaxLumiRangePerRun, unsigned int &RunSize, unsigned int *JS_Run, unsigned int *JS_Run_NLumi, unsigned int *JS_LumiMin, unsigned int *JS_LumiMax, TString fileName)
{
  RunSize = 0;
  unsigned int iLumi = 0;
  int Depth = -1;
  int SubDepth0 = -1;
  int SubDepth2 = -1;
  if( fileName == "one" )
  {
    return 0;
  }
  fileName = "JSonFiles/"+fileName;

  ifstream in;
  in.open(fileName); 
  string line;

  int minI = 1000;

  string Delimiters[NDelimiters] = {"{","}","[","]",":",",","\""};
  int iDelimiterMin = -1; 

  float tmp;

  if( in.is_open() )
  {
    while ( getline(in, line) )
    {
      while(line.size() > 0)
      {
        minI = 1000;
        int tmpI = 1000;

        for( int i = 0; i < NDelimiters; ++i )
        {
          tmpI = line.find_first_of(Delimiters[i]);
          if( tmpI != string::npos && tmpI < minI )
          {
            minI = tmpI;
            iDelimiterMin = i;
          }
        }
        if( Depth == 2 )
        {
          if( iLumi >= N_MaxLumiRangePerRun ) return -2;
          if( Delimiters[iDelimiterMin] == "," )
          {
            ++SubDepth2;
            TString tmpLine = line.substr(0,minI );
            JS_LumiMin[RunSize*N_MaxLumiRangePerRun+iLumi] = tmpLine.Atoi();
          }
          if( Delimiters[iDelimiterMin] == "]" )
          {
            ++SubDepth2;
            TString tmpLine = line.substr(0,minI );
            JS_LumiMax[RunSize*N_MaxLumiRangePerRun+iLumi] = tmpLine.Atoi();
            cout << JS_Run[RunSize] <<  "|||" << JS_LumiMin[RunSize*N_MaxLumiRangePerRun+iLumi] << " ; " << JS_LumiMax[RunSize*N_MaxLumiRangePerRun+iLumi] << endl;
            iLumi++;
          }

        } else
        {
          SubDepth2 = -1;
        }

        if( Delimiters[iDelimiterMin] == "{" ) ++Depth;
        if( Delimiters[iDelimiterMin] == "[" ) ++Depth;
        if( Delimiters[iDelimiterMin] == "}" ) --Depth;
        if( Delimiters[iDelimiterMin] == "]" ) --Depth;


        if( Depth == 0 )
        {
          if( Delimiters[iDelimiterMin] == "\"" ) ++SubDepth0;
          if( Delimiters[iDelimiterMin] == ":" ) ++SubDepth0;

          if( SubDepth0 == 1 ) // New Run
          {
            if( RunSize < N_MaxRun ) JS_Run_NLumi[RunSize] = iLumi;
            //if( RunSize >= 0 && RunSize < N_MaxRun ) JS_Run_NLumi[RunSize] = iLumi;
            RunSize++;
            if( RunSize >= N_MaxRun ) return -3;
            iLumi = 0;
            TString tmpRun = line.substr(0,minI);
            JS_Run[RunSize] = tmpRun.Atoi();
          }
        } else
        {
          SubDepth0 = -1;
        }
        
        line.erase(0,minI+1);

        
      }
    }
  }
  else
  {
    cout << "WARNING: Fail to open file " << fileName << endl;
    RunSize = 0;
    return -1;
  }

  in.close();
  cout << "===============" << endl;
  cout << RunSize << endl;
  return RunSize;
}

  // float WeightNS[18*52];
  // float WeightNS[N_ptbin_nseltrack*N_SelTracks];
void ReadNSelTracksWeight(float *WeightNS, unsigned int N_SelTracks, unsigned int N_ptbin_nseltrack,
                          TString fileName, TString *St_fhseltrackPtMin, TString *St_fhseltrackPtMax)
{
  if( fileName == "one" )
  {
    for( unsigned int i = 0; i < N_ptbin_nseltrack; ++i )
    for( unsigned int j = 0; j < N_SelTracks; ++j )
    {
       WeightNS[i*N_SelTracks+j] = 1;
    }
    return;
  }
  fileName = "Weights/"+fileName;

  for( unsigned int iPt = 0; iPt < N_ptbin_nseltrack; ++iPt )
  {
    ifstream in;
    in.open(fileName+St_fhseltrackPtMin[iPt]+"to"+St_fhseltrackPtMax[iPt]+".txt"); // change name of file

    float tmp;
    unsigned int nlines=0;

    if( in.is_open() )
    {
      while (1) {
        if( nlines >= N_SelTracks )
        {
          in >> tmp; // not used
          if( in.good() )
            cout << "WARNING: file " << fileName << " contains more than " << N_SelTracks << " lines ! " << endl;
          break;
        }
  
        in >> tmp;
        if( !in.good() ) break;
        WeightNS[iPt*N_SelTracks+nlines] = tmp;
        nlines++;
     }
    }
    else cout << "WARNING: Fail to open file " << fileName << endl;

    for( int i = nlines; i < N_SelTracks; ++i )
    {
       cout << "Warning, file " << fileName<< St_fhseltrackPtMin[iPt]<< "to"<< St_fhseltrackPtMax[iPt]<< ".txt" << " does not contain enough lines ! " << endl;
       cout << "Setting remaining weights to 1 " << endl;
       WeightNS[iPt*N_SelTracks+i] = 1;
    }
    in.close();
  }
}

int ReadPUWeight(float *WeightPU, int maxPU, TString fileName = "one")
{
  if( fileName == "one" )
  {
    for( int i = 0; i < maxPU; ++i )
    {
       WeightPU[i] = 1;
    }
    return 1;
  }
  fileName = "Weights/"+fileName;
  ifstream in;
  in.open(fileName); // change name of file

  float tmp;
  int nlines=0;

  if( in.is_open() )
  {
    while (1) {
      if( nlines >= maxPU )
      {
        in >> tmp; // not used
        if( in.good() )
          cout << "WARNING: file " << fileName << " contains more than " << maxPU << " lines ! " << endl;
        break;
      }

      in >> tmp; // not used
      if( !in.good() ) break;
      if( tmp - nlines > 0.1 || tmp - nlines < -0.1 )
      {
        cout << "Warning: Index seems not correct, we have tmp = " << tmp << " while nlines = " << nlines << endl;
      }
  
      in >> tmp; // not used
      if( !in.good() ) break;
      WeightPU[nlines] = tmp;
  
      nlines++;
     }
     return 2;
  }
  else
  {
    cout << "WARNING: Fail to open file " << fileName << endl;

    for( int i = nlines; i < maxPU; ++i )
    {
       WeightPU[i] = 0;
    }
    return 0;
  }
}

void ReadpthatWeight(float *pthatMin, float *pthatMax, float *pthatweight, int maxLines, TString fileName = "one")
{
  if( fileName == "one" )
  {
    for( int i = 0; i < maxLines; ++i )
    {
       pthatMin[i] = 0;
       pthatMax[i] = 999999;
       pthatweight[i] = 1;
    }
    return;
  }

  fileName = "Weights/"+fileName;
  ifstream in;
  in.open(fileName); // change name of file


  float tmp;
  int nlines=0;

  if( in.is_open() )
  {
    while (1) {
      if( nlines >= maxLines )
      {
        in >> tmp; // not used
        if( in.good() )
          cout << "WARNING: file " << fileName << " contains more than " << maxLines << " lines ! " << endl;
        break;
      }

      in >> tmp; // not used
      if( !in.good() ) break;
      pthatMin[nlines] = tmp;
  
      in >> tmp; // not used
      if( !in.good() ) break;
      pthatMax[nlines] = tmp;
  
      in >> tmp; // not used
      if( !in.good() ) break;
      pthatweight[nlines] = tmp;
  
      nlines++;
     }
  }
  else cout << "WARNING: Fail to open file " << fileName << endl;

  for( int i = nlines; i < maxLines; ++i )
  {
     pthatMin[i] = 99999;
     pthatMax[i] = 999999;
     pthatweight[i] = 0;
  }
}

void SetSebStyle()
{
 gStyle->SetTitleFillColor(42);
 gStyle->SetTitleFont(1);
 gStyle->SetStatColor(29);
 gStyle->SetCanvasColor(25);   
 gStyle->SetOptStat(1111111);
 gStyle->SetHistFillColor(5);
}

void SaveInFile(TH1* ahisto, TFile* afile)
{
 if (!ahisto) { std::cout << "!! no histo !!" << std::endl; return ;}
 TDirectory* current = gDirectory ;
 afile->cd();
 ahisto->Write();
 current->cd();
}

void SubJetsAna2D::Loop(SubJets *FJ, TString TagName, float aPtMin, float aPtMax, 
                float aFreeCut, int aIntCut,
                float minCutJetPtMax, float maxCutJetPtMax,
		            TString afilename, TString weightPU_file,
                TString weightPthat_file, TString JSONFile, bool truePU, bool WeightTracks, TString TrigType, TString period)
{
cout << " Re PUFile : " << weightPU_file << endl;
cout << " Re PthatFile : " << weightPthat_file << endl;

TString TagLevel[3] = {"L", "M", "T" };
///////////////////////////////////////////////////////////////////
//                  T I T L E      C A R D S                     //

//$$
  float  TagCut[3]       = {0,0,0};        // tag cut to be determined from TagName and TagLevel
  float  PtMin        = aPtMin;   // pt jet min
  float  PtMax        = aPtMax;   // pt jet max
  float  FreeCut      = aFreeCut; 
  int    IntCut       = aIntCut; 
  TString  filename   = afilename;

// Operating points from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X 
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80X
  if( TagName == "JP" )
  {
    TagCut[0] = .245;
    TagCut[1] = .515;
    TagCut[2] = .760;
  }
  else if( TagName == "CSVv2" )
  {
    //TagCut[0] = .460;
    //TagCut[1] = .800;
    //TagCut[2] = .935;
    TagCut[0] = .5426;
    TagCut[1] = .8484;
    TagCut[2] = .9535;
  }
  else if( TagName == "cMVAv2" )
  {
    //TagCut[0] = -0.715;
    //TagCut[1] =  0.185;
    //TagCut[2] =  0.875;
    TagCut[0] = -0.5884;
    TagCut[1] =  0.4432;
    TagCut[2] =  0.9432;
  }
  else
  {
    cout << "Undefined TagName for " << TagName << endl;
    return;
  }
//$$
  int NtrackMin = 1;         // Taggability
//
//  int NtrackMin = int(FreeCut);
//
//  float VetoPos = 4.;
//  float VetoPos = FreeCut;
  float VetoPos = 1000.;
//$$

  int NpuMin = 0;
//$$
//$$  NpuMin = int(FreeCut);
//$$

//$$
  bool GluonSplit = false;
  bool ReweightNSelTracks = false; // To activate the reweighting according to the nSelTracks.
  //TString weightNS_file = "NSEL_ak5Jets/nselWeight_pt_"; // Set to "one" to set all weights to 1.
  TString weightNS_file = "one"; // Set to "one" to set all weights to 1.

  bool BTagMu  = false;
  bool MuonJet = false;  // need to have MuOp    = false
  bool MuOp    = false;  // need to have MuonJet = false and Away = false
  bool Away    = false;
  bool EtaBin  = false;

  int Year = 2011;
  if( TrigType == "2011" ) Year = 2011;
  else if( TrigType == "2016" || TrigType == "2015" || TrigType == "2012" ||  TrigType == "2015_Cert" || TrigType == "2015_LooserSel") Year = 2012;
  else if(TrigType == "Moriond17") Year = 2016 ;
  //TEMP ignore trigger use jet pT cut
  else if( TrigType == "ignore") Year = 0 ;
  else
  {
    cout << "This TrigType is not foreseen ( " << TrigType << ")" << endl;
  }

  cout << "\n Trigger type is: " << TrigType << "  " << Year << endl ;

  bool AlaLuca = false;
  bool AlaPablo = false;

  bool GSplitMin = false;
  bool GSplitMax = false;

  bool BFragMin = false;
  bool BFragMax = false;

  bool GenTrkMin = false;
  bool GenTrkMax = false;

  bool CFrag = false;
  bool V0weight = false;

  bool dRJetCut = false;

  bool TopAna = false;
//$$

///////////////////////////////////////////////////////////////////
 
//$$
  TH1F::SetDefaultSumw2(kTRUE);
//$$

//**********************************
// Data
//**********************************
  TH1F* hBeforeTrig_MaxJetPt       = new TH1F("hBeforeTrig_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig30_MaxJetPt       = new TH1F("hAfterTrig30_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig60_MaxJetPt       = new TH1F("hAfterTrig60_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig80_MaxJetPt       = new TH1F("hAfterTrig80_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig110_MaxJetPt       = new TH1F("hAfterTrig110_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig150_MaxJetPt       = new TH1F("hAfterTrig150_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig190_MaxJetPt       = new TH1F("hAfterTrig190_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig240_MaxJetPt       = new TH1F("hAfterTrig240_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig260_MaxJetPt       = new TH1F("hAfterTrig260_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig300_MaxJetPt       = new TH1F("hAfterTrig300_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig320_MaxJetPt       = new TH1F("hAfterTrig320_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig370_MaxJetPt       = new TH1F("hAfterTrig370_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig400_MaxJetPt       = new TH1F("hAfterTrig400_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hAfterTrig500_MaxJetPt       = new TH1F("hAfterTrig500_MaxJetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  int NPVMAX = 40 ;
  float NPVMAX_VAL = 40.5 ;
  TH1F* hData_All_nPV         = new TH1F("hData_All_nPV","nb. of PV",NPVMAX,0.5,NPVMAX_VAL);
  TH1F* hData_All_NJets       = new TH1F("hData_All_NJets","nb. of jets",16,-0.5,15.5);
  TH1F* hData_All_tracks      = new TH1F("hData_All_tracks","nb. of tracks",20,0.5,20.5);
  TH1F* hData_All_JetPV       = new TH1F("hData_All_JetPV","#PV",12,0.5,24.5);
  TH2F* hData_All_JetPt       = new TH2F("hData_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hData_All_JetEta      = new TH1F("hData_All_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hData_All_JetRun      = new TH1F("hData_All_JetRun","Run number",400,190700.5,191100.5);
  TH1F* hData_nPV             = new TH1F("hData_nPV","nb. of PV",NPVMAX,0.5,NPVMAX_VAL);
  TH1F* hData_nPV_test             = new TH1F("hData_nPV_test","nb. of PV",NPVMAX,0.5,NPVMAX_VAL);
  TH1F* hData_nPV_noWeight    = new TH1F("hData_nPV_noWeight","nb. of PV",NPVMAX,0.5,NPVMAX_VAL) ;
  TH1F* hData_nPV_noWeight_allEvent = new TH1F("hData_nPV_noWeight_allEvent","nb. of PV",NPVMAX,0.5,NPVMAX_VAL) ;
  TH1F* hData_NJets           = new TH1F("hData_NJets","nb. of jets",16,-0.5,15.5);
  TH1F* hData_NJet30          = new TH1F("hData_NJet30","nb. jets pt>30",16,-0.5,15.5);
  TH1F* hData_JetPV           = new TH1F("hData_JetPV","#PV",12,0.5,24.5);
  TH2F* hData_JetPt           = new TH2F("hData_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hData_JetEta          = new TH1F("hData_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hData_JetRun          = new TH1F("hData_JetRun","Run number",400,190700.5,191100.5);
  TH2F* hData_Tagger           = new TH2F("hData_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hData_JetPt_LE5pv    = new TH1F("hData_JetPt_LE5pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hData_JetEta_LE5pv   = new TH1F("hData_JetEta_LE5pv","|#eta(jet)|",24,0.,2.4);
  TH1F* hData_JetPt_GE6pv    = new TH1F("hData_JetPt_GE6pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hData_JetEta_GE6pv   = new TH1F("hData_JetEta_GE6pv","|#eta(jet)|",24,0.,2.4);
  TH1F* hData_1JetPt           = new TH1F("hData_1JetPt","pt(jet)",30,20.,320.);
  TH1F* hData_2JetPt           = new TH1F("hData_2JetPt","pt(jet)",30,20.,320.);
  TH1F* hData_3JetPt           = new TH1F("hData_3JetPt","pt(jet)",30,20.,320.);
  TH1F* hData_4JetPt           = new TH1F("hData_4JetPt","pt(jet)",30,20.,320.);
  TH1F* hData_PosTag           = new TH1F("hData_PosTag","Tag(+)",50,0.,25.);
  TH1F* hData_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hData_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hData_PosTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hData_PosTag_JetRun[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hData_PosTag_JetPV[i]    = new TH1F("hData_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hData_PosTag_JetPt[i]    = new TH2F("hData_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hData_PosTag_JetEta[i]   = new TH1F("hData_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hData_PosTag_JetRun[i]   = new TH1F("hData_PosTag"+TagLevel[i]+"_JetRun","Run number",70,190400,191100);
  }

  TH1F* hData_NegTag_JetPV[3] = {NULL, NULL, NULL};
  TH2F* hData_NegTag_JetPt[3] = {NULL, NULL, NULL};    
  TH1F* hData_NegTag_JetEta[3] = {NULL, NULL, NULL};   
  TH1F* hData_NegTag_JetRun[3] = {NULL, NULL, NULL};   
  TH1F* hData_NegTag_JetPt_LE5pv[3] = {NULL, NULL, NULL};  
  TH1F* hData_NegTag_JetEta_LE5pv[3] = {NULL, NULL, NULL}; 
  TH1F* hData_NegTag_JetPt_GE6pv[3] = {NULL, NULL, NULL};  
  TH1F* hData_NegTag_JetEta_GE6pv[3] = {NULL, NULL, NULL}; 
  TH1F* hData_NegTag_1JetPt[3] = {NULL, NULL, NULL};    
  TH1F* hData_NegTag_2JetPt[3] = {NULL, NULL, NULL};    
  TH1F* hData_NegTag_3JetPt[3] = {NULL, NULL, NULL};    
  TH1F* hData_NegTag_4JetPt[3] = {NULL, NULL, NULL};    

  for(int i = 0; i < 3; ++i )
  {
    hData_NegTag_JetPV[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hData_NegTag_JetPt[i] = new TH2F("hData_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hData_NegTag_JetEta[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hData_NegTag_JetRun[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_JetRun","Run number",70,190400,191100);

    hData_NegTag_JetPt_LE5pv[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_JetPt_LE5pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
    hData_NegTag_JetEta_LE5pv[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_JetEta_LE5pv","|#eta(jet)|",24,0.,2.4);
    hData_NegTag_JetPt_GE6pv[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_JetPt_GE6pv","pt(jet)",NPtBin, minPtBin, maxPtBin);
    hData_NegTag_JetEta_GE6pv[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_JetEta_GE6pv","|#eta(jet)|",24,0.,2.4);


    hData_NegTag_1JetPt[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_1JetPt","pt(jet)",30,20.,320.);
    hData_NegTag_2JetPt[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_2JetPt","pt(jet)",30,20.,320.);
    hData_NegTag_3JetPt[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_3JetPt","pt(jet)",30,20.,320.);
    hData_NegTag_4JetPt[i] = new TH1F("hData_NegTag"+TagLevel[i]+"_4JetPt","pt(jet)",30,20.,320.);
  }

  TH1F* hData_Trigger     = new TH1F("hData_Trigger","Trigger",100,-0.5,99.5);
  TH1F* hData_Trig1000    = new TH1F("hData_Trig1000","Trig1000",1000,-0.5,999.5);
  TH1F* hData_Trig1_JetPt = new TH1F("hData_Trig1_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hData_Trig2_JetPt = new TH1F("hData_Trig2_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hData_Trig3_JetPt = new TH1F("hData_Trig3_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hData_Trig4_JetPt = new TH1F("hData_Trig4_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hData_Trig5_JetPt = new TH1F("hData_Trig5_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hData_Trig6_JetPt = new TH1F("hData_Trig6_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin);

//**********************************
// All flavours in Monte Carlo
//**********************************
  TH1F* hAllFlav_pthatAll    = new TH1F("hAllFlav_pthatAll","pthat",80,0.,800.);
  TH1F* hAllFlav_pthat       = new TH1F("hAllFlav_pthat","pthat",80,0.,800.);
  TH1F* hAllFlav_pthatTrig   = new TH1F("hAllFlav_pthatTrig","pthat",80,0.,800.);
  TH1F* hAllFlav_All_nPU     = new TH1F("hAllFlav_All_nPU","nb. of PU",60,0.,60.);
  TH1F* hAllFlav_nPU         = new TH1F("hAllFlav_nPU","nb. of PU",60,0.,60.);
  TH1F* hAllFlav_All_Flavour = new TH1F("hAllFlav_All_Flavour","Flavour",22,-0.5,21.5);
  TH1F* hAllFlav_Flavour     = new TH1F("hAllFlav_Flavour","Flavour",22,-0.5,21.5);
  TH1F* hAllFlav_Tagger_Bwd  = new TH1F("hAllFlav_Tagger_Bwd","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_Cwd  = new TH1F("hAllFlav_Tagger_Cwd","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_Tau  = new TH1F("hAllFlav_Tagger_Tau","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_Gam  = new TH1F("hAllFlav_Tagger_Gam","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_K0s  = new TH1F("hAllFlav_Tagger_K0s","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_Lam  = new TH1F("hAllFlav_Tagger_Lam","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_Int  = new TH1F("hAllFlav_Tagger_Int","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_Fak  = new TH1F("hAllFlav_Tagger_Fak","Tagger",100,-25.,25.);
  TH1F* hAllFlav_Tagger_Oth  = new TH1F("hAllFlav_Tagger_Oth","Tagger",100,-25.,25.);

  TH1F* hAllFlav_Gam_JetPV  	   = new TH1F("hAllFlav_Gam_JetPV","#PV",12,0.5,24.5);
  TH2F* hAllFlav_Gam_JetPt	   = new TH2F("hAllFlav_Gam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hAllFlav_Gam_JetEta	   = new TH1F("hAllFlav_Gam_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hAllFlav_K0s_JetPV  	   = new TH1F("hAllFlav_K0s_JetPV","#PV",12,0.5,24.5);
  TH2F* hAllFlav_K0s_JetPt	   = new TH2F("hAllFlav_K0s_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hAllFlav_K0s_JetEta	   = new TH1F("hAllFlav_K0s_JetEta","|#eta(jet)|",24,0.,2.4);

  TH1F* hAllFlav_Gam_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hAllFlav_Gam_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_Gam_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_K0s_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hAllFlav_K0s_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_K0s_NegTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hAllFlav_Gam_NegTag_JetPV[i]   = new TH1F("hAllFlav_Gam_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hAllFlav_Gam_NegTag_JetPt[i]   = new TH2F("hAllFlav_Gam_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Gam_NegTag_JetEta[i]  = new TH1F("hAllFlav_Gam_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hAllFlav_K0s_NegTag_JetPV[i]   = new TH1F("hAllFlav_K0s_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hAllFlav_K0s_NegTag_JetPt[i]   = new TH2F("hAllFlav_K0s_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_K0s_NegTag_JetEta[i]  = new TH1F("hAllFlav_K0s_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }
//
  TH1F* hAllFlav_K0s_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hAllFlav_K0s_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_K0s_PosTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_Lam_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hAllFlav_Lam_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_Lam_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hAllFlav_K0s_PosTag_JetPV[i]   = new TH1F("hAllFlav_K0s_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hAllFlav_K0s_PosTag_JetPt[i]   = new TH2F("hAllFlav_K0s_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_K0s_PosTag_JetEta[i]  = new TH1F("hAllFlav_K0s_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hAllFlav_Lam_PosTag_JetPV[i]   = new TH1F("hAllFlav_Lam_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hAllFlav_Lam_PosTag_JetPt[i]   = new TH2F("hAllFlav_Lam_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Lam_PosTag_JetEta[i]  = new TH1F("hAllFlav_Lam_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }
//
  TH1F* hAllFlav_Lam_JetPV  	   = new TH1F("hAllFlav_Lam_JetPV","#PV",12,0.5,24.5);
  TH2F* hAllFlav_Lam_JetPt	   = new TH2F("hAllFlav_Lam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hAllFlav_Lam_JetEta	   = new TH1F("hAllFlav_Lam_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hAllFlav_Fak_JetPV  	   = new TH1F("hAllFlav_Fak_JetPV","#PV",12,0.5,24.5);
  TH2F* hAllFlav_Fak_JetPt	   = new TH2F("hAllFlav_Fak_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hAllFlav_Fak_JetEta	   = new TH1F("hAllFlav_Fak_JetEta","|#eta(jet)|",24,0.,2.4);

  TH1F* hAllFlav_Lam_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hAllFlav_Lam_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_Lam_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_Fak_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hAllFlav_Fak_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hAllFlav_Fak_NegTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hAllFlav_Lam_NegTag_JetPV[i]   = new TH1F("hAllFlav_Lam_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hAllFlav_Lam_NegTag_JetPt[i]   = new TH2F("hAllFlav_Lam_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Lam_NegTag_JetEta[i]  = new TH1F("hAllFlav_Lam_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hAllFlav_Fak_NegTag_JetPV[i]   = new TH1F("hAllFlav_Fak_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hAllFlav_Fak_NegTag_JetPt[i]   = new TH2F("hAllFlav_Fak_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hAllFlav_Fak_NegTag_JetEta[i]  = new TH1F("hAllFlav_Fak_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }

//**********************************
// udsg-jets
//**********************************
  TH1F* hLightFlav_All_JetPV      = new TH1F("hLightFlav_All_JetPV","#PV",12,0.5,24.5);
  TH2F* hLightFlav_All_JetPt      = new TH2F("hLightFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_All_JetEta     = new TH1F("hLightFlav_All_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hLightFlav_JetPU          = new TH1F("hLightFlav_JetPU","#PU",12,0.5,24.5);
  TH1F* hLightFlav_JetPV          = new TH1F("hLightFlav_JetPV","#PV",12,0.5,24.5);
  TH2F* hLightFlav_JetPt          = new TH2F("hLightFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_JetEta         = new TH1F("hLightFlav_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hLightFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hLightFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hLightFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hLightFlav_NegTag_JetPV[i]   = new TH1F("hLightFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_NegTag_JetPt[i]   = new TH2F("hLightFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_NegTag_JetEta[i]  = new TH1F("hLightFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }

  TH2F* hLightFlav_JetPt_0pu      = new TH2F("hLightFlav_JetPt_0pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2F* hLightFlav_JetPt_GE8pu    = new TH2F("hLightFlav_JetPt_GE8pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2F* hLightFlav_Tagger         = new TH2F("hLightFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_Tagger_Gam     = new TH1F("hLightFlav_Tagger_Gam","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_K0s     = new TH1F("hLightFlav_Tagger_K0s","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_Lam     = new TH1F("hLightFlav_Tagger_Lam","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_Bwd     = new TH1F("hLightFlav_Tagger_Bwd","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_Cwd     = new TH1F("hLightFlav_Tagger_Cwd","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_Tau     = new TH1F("hLightFlav_Tagger_Tau","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_Int     = new TH1F("hLightFlav_Tagger_Int","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_Fak     = new TH1F("hLightFlav_Tagger_Fak","Tagger",100,-25.,25.);
  TH1F* hLightFlav_Tagger_Oth     = new TH1F("hLightFlav_Tagger_Oth","Tagger",100,-25.,25.);
  TH1F* hLightFlav_1JetPt         = new TH1F("hLightFlav_1JetPt","pt(jet)",30,20.,320.);
  TH1F* hLightFlav_2JetPt         = new TH1F("hLightFlav_2JetPt","pt(jet)",30,20.,320.);
  TH1F* hLightFlav_3JetPt         = new TH1F("hLightFlav_3JetPt","pt(jet)",30,20.,320.);
  TH1F* hLightFlav_4JetPt         = new TH1F("hLightFlav_4JetPt","pt(jet)",30,20.,320.);

  TH1F* hLightFlav_PosTag_JetPU[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2F* hLightFlav_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_PosTagger      = new TH1F("hLightFlav_PosTagger","Tagger",1000,0.,25.);
  TH1F* hLightFlav_PosTagger_0pu  = new TH1F("hLightFlav_PosTagger_0pu","Tagger",1000,0.,25.);
  TH1F* hLightFlav_PosTagger_GE8pu= new TH1F("hLightFlav_PosTagger_GE8pu","Tagger",1000,0.,25.);
  TH1F* hLightFlav_PosTag_1JetPt[3] = {NULL, NULL, NULL};  
  TH1F* hLightFlav_PosTag_2JetPt[3] = {NULL, NULL, NULL};  
  TH1F* hLightFlav_PosTag_3JetPt[3] = {NULL, NULL, NULL};  
  TH1F* hLightFlav_PosTag_4JetPt[3] = {NULL, NULL, NULL};  
  TH1F* hLightFlav_BCT_PosTag_JetPV[3] = {NULL, NULL, NULL};   
  TH2F* hLightFlav_BCT_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_BCT_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Gam_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2F* hLightFlav_Gam_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Gam_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_K0s_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2F* hLightFlav_K0s_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_K0s_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Lam_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2F* hLightFlav_Lam_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Lam_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Fak_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2F* hLightFlav_Fak_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Fak_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Fak_PosTag_pthat[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Fak_PosTag_nPU[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Oth_PosTag_JetPV[3] = {NULL, NULL, NULL};
  TH2F* hLightFlav_Oth_PosTag_JetPt[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Oth_PosTag_JetEta[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Oth_PosTag_pthat[3] = {NULL, NULL, NULL};
  TH1F* hLightFlav_Oth_PosTag_nPU[3] = {NULL, NULL, NULL};

  for( int i = 0; i < 3; ++i )
  {
    hLightFlav_PosTag_JetPU[i]   = new TH1F("hLightFlav_PosTag"+TagLevel[i]+"_JetPU","#PU",12,0.5,24.5);
    hLightFlav_PosTag_JetPV[i]   = new TH1F("hLightFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_PosTag_JetPt[i]   = new TH2F("hLightFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_PosTag_JetEta[i]  = new TH1F("hLightFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_PosTag_1JetPt[i]  = new TH1F("hLightFlav_PosTag"+TagLevel[i]+"_1JetPt","pt(jet)",30,20.,320.);
    hLightFlav_PosTag_2JetPt[i]  = new TH1F("hLightFlav_PosTag"+TagLevel[i]+"_2JetPt","pt(jet)",30,20.,320.);
    hLightFlav_PosTag_3JetPt[i]  = new TH1F("hLightFlav_PosTag"+TagLevel[i]+"_3JetPt","pt(jet)",30,20.,320.);
    hLightFlav_PosTag_4JetPt[i]  = new TH1F("hLightFlav_PosTag"+TagLevel[i]+"_4JetPt","pt(jet)",30,20.,320.);
    hLightFlav_BCT_PosTag_JetPV[i]   = new TH1F("hLightFlav_BCT_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_BCT_PosTag_JetPt[i]   = new TH2F("hLightFlav_BCT_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_BCT_PosTag_JetEta[i]  = new TH1F("hLightFlav_BCT_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_Gam_PosTag_JetPV[i]   = new TH1F("hLightFlav_Gam_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_Gam_PosTag_JetPt[i]   = new TH2F("hLightFlav_Gam_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Gam_PosTag_JetEta[i]  = new TH1F("hLightFlav_Gam_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_K0s_PosTag_JetPV[i]   = new TH1F("hLightFlav_K0s_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_K0s_PosTag_JetPt[i]   = new TH2F("hLightFlav_K0s_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_K0s_PosTag_JetEta[i]  = new TH1F("hLightFlav_K0s_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_Lam_PosTag_JetPV[i]   = new TH1F("hLightFlav_Lam_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_Lam_PosTag_JetPt[i]   = new TH2F("hLightFlav_Lam_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Lam_PosTag_JetEta[i]  = new TH1F("hLightFlav_Lam_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_Fak_PosTag_JetPV[i]   = new TH1F("hLightFlav_Fak_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_Fak_PosTag_JetPt[i]   = new TH2F("hLightFlav_Fak_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Fak_PosTag_JetEta[i]  = new TH1F("hLightFlav_Fak_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_Fak_PosTag_pthat[i]  = new TH1F("hLightFlav_Fak_PosTag"+TagLevel[i]+"_pthat","pthat",80,0.,800.);
    hLightFlav_Fak_PosTag_nPU[i]  = new TH1F("hLightFlav_Fak_PosTag"+TagLevel[i]+"_nPU","nb. of PU",60,0.,60.);
    hLightFlav_Oth_PosTag_JetPV[i]   = new TH1F("hLightFlav_Oth_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_Oth_PosTag_JetPt[i]   = new TH2F("hLightFlav_Oth_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Oth_PosTag_JetEta[i]  = new TH1F("hLightFlav_Oth_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_Oth_PosTag_pthat[i]  = new TH1F("hLightFlav_Oth_PosTag"+TagLevel[i]+"_pthat","pthat",80,0.,800.);
    hLightFlav_Oth_PosTag_nPU[i]  = new TH1F("hLightFlav_Oth_PosTag"+TagLevel[i]+"_nPU","nb. of PU",60,0.,60.);
  }


  TH1F* hLightFlav_Gam_JetPV  	     = new TH1F("hLightFlav_Gam_JetPV","#PV",12,0.5,24.5);
  TH2F* hLightFlav_Gam_JetPt	     = new TH2F("hLightFlav_Gam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_Gam_JetEta	     = new TH1F("hLightFlav_Gam_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hLightFlav_K0s_JetPV  	     = new TH1F("hLightFlav_K0s_JetPV","#PV",12,0.5,24.5);
  TH2F* hLightFlav_K0s_JetPt	     = new TH2F("hLightFlav_K0s_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_K0s_JetEta	     = new TH1F("hLightFlav_K0s_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hLightFlav_Lam_JetPV  	     = new TH1F("hLightFlav_Lam_JetPV","#PV",12,0.5,24.5);
  TH2F* hLightFlav_Lam_JetPt	     = new TH2F("hLightFlav_Lam_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_Lam_JetEta	     = new TH1F("hLightFlav_Lam_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hLightFlav_Fak_JetPV  	     = new TH1F("hLightFlav_Fak_JetPV","#PV",12,0.5,24.5);
  TH2F* hLightFlav_Fak_JetPt	     = new TH2F("hLightFlav_Fak_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_Fak_JetEta	     = new TH1F("hLightFlav_Fak_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hLightFlav_Oth_JetPV  	     = new TH1F("hLightFlav_Oth_JetPV","#PV",12,0.5,24.5);
  TH2F* hLightFlav_Oth_JetPt	     = new TH2F("hLightFlav_Oth_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hLightFlav_Oth_JetEta	     = new TH1F("hLightFlav_Oth_JetEta","|#eta(jet)|",24,0.,2.4);

  TH1F* hLightFlav_Fak_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hLightFlav_Fak_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hLightFlav_Fak_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hLightFlav_Oth_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hLightFlav_Oth_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hLightFlav_Oth_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hLightFlav_Oth_NegTag_pthat[3] = {NULL, NULL, NULL };
  TH1F* hLightFlav_Oth_NegTag_nPU[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hLightFlav_Fak_NegTag_JetPV[i]   = new TH1F("hLightFlav_Fak_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_Fak_NegTag_JetPt[i]   = new TH2F("hLightFlav_Fak_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Fak_NegTag_JetEta[i]  = new TH1F("hLightFlav_Fak_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_Oth_NegTag_JetPV[i]   = new TH1F("hLightFlav_Oth_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hLightFlav_Oth_NegTag_JetPt[i]   = new TH2F("hLightFlav_Oth_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hLightFlav_Oth_NegTag_JetEta[i]  = new TH1F("hLightFlav_Oth_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hLightFlav_Oth_NegTag_pthat[i]  = new TH1F("hLightFlav_Oth_NegTag"+TagLevel[i]+"_pthat","pthat",80,0.,800.);
    hLightFlav_Oth_NegTag_nPU[i]  = new TH1F("hLightFlav_Oth_NegTag"+TagLevel[i]+"_nPU","nb. of PU",60,0.,60.);
  }

//**********************************
// gluon-jets
//**********************************
  TH1F* hGluonFlav_All_JetPV      = new TH1F("hGluonFlav_All_JetPV","#PV",12,0.5,24.5);
  TH2F* hGluonFlav_All_JetPt      = new TH2F("hGluonFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hGluonFlav_All_JetEta     = new TH1F("hGluonFlav_All_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hGluonFlav_JetPV          = new TH1F("hGluonFlav_JetPV","#PV",12,0.5,24.5);
  TH2F* hGluonFlav_JetPt          = new TH2F("hGluonFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hGluonFlav_JetEta         = new TH1F("hGluonFlav_JetEta","|#eta(jet)|",24,0.,2.4);
  TH2F* hGluonFlav_Tagger         = new TH2F("hGluonFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hGluonFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hGluonFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hGluonFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hGluonFlav_PosTag = new TH1F("hGluonFlav_PosTag","Tag(+)",50,0.,25.);
  TH1F* hGluonFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hGluonFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hGluonFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hGluonFlav_NegTag_JetPV[i]   = new TH1F("hGluonFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hGluonFlav_NegTag_JetPt[i]   = new TH2F("hGluonFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hGluonFlav_NegTag_JetEta[i]  = new TH1F("hGluonFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);

    hGluonFlav_PosTag_JetPV[i]   = new TH1F("hGluonFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hGluonFlav_PosTag_JetPt[i]   = new TH2F("hGluonFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hGluonFlav_PosTag_JetEta[i]  = new TH1F("hGluonFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }

//**********************************
// uds-jets
//**********************************
  TH1F* hUDSFlav_All_JetPV      = new TH1F("hUDSFlav_All_JetPV","#PV",12,0.5,24.5);
  TH2F* hUDSFlav_All_JetPt      = new TH2F("hUDSFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hUDSFlav_All_JetEta     = new TH1F("hUDSFlav_All_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hUDSFlav_JetPV          = new TH1F("hUDSFlav_JetPV","#PV",12,0.5,24.5);
  TH2F* hUDSFlav_JetPt          = new TH2F("hUDSFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hUDSFlav_JetEta         = new TH1F("hUDSFlav_JetEta","|#eta(jet)|",24,0.,2.4);
  TH2F* hUDSFlav_Tagger         = new TH2F("hUDSFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hUDSFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hUDSFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hUDSFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hUDSFlav_PosTag = new TH1F("hUDSFlav_PosTag","Tag(+)",50,0.,25.);
  TH1F* hUDSFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hUDSFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hUDSFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hUDSFlav_NegTag_JetPV[i]   = new TH1F("hUDSFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hUDSFlav_NegTag_JetPt[i]   = new TH2F("hUDSFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hUDSFlav_NegTag_JetEta[i]  = new TH1F("hUDSFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hUDSFlav_PosTag_JetPV[i]   = new TH1F("hUDSFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hUDSFlav_PosTag_JetPt[i]   = new TH2F("hUDSFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hUDSFlav_PosTag_JetEta[i]  = new TH1F("hUDSFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }

//**********************************
// c-jets
//**********************************
  TH1F* hCFlav_All_JetPV      = new TH1F("hCFlav_All_JetPV","#PV",12,0.5,24.5);
  TH2F* hCFlav_All_JetPt      = new TH2F("hCFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hCFlav_All_JetEta     = new TH1F("hCFlav_All_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hCFlav_JetPV          = new TH1F("hCFlav_JetPV","#PV",12,0.5,24.5);
  TH2F* hCFlav_JetPt          = new TH2F("hCFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hCFlav_JetEta         = new TH1F("hCFlav_JetEta","|#eta(jet)|",24,0.,2.4);
  TH2F* hCFlav_Tagger         = new TH2F("hCFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hCFlav_deltaR         = new TH1F("hCFlav_deltaR","#Delta R gluon split",50,0.,1.);
  TH2F* hCFlav_JetPt_0pu      = new TH2F("hCFlav_JetPt_0pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2F* hCFlav_JetPt_GE8pu    = new TH2F("hCFlav_JetPt_GE8pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hCFlav_PosTagger      = new TH1F("hCFlav_PosTagger","Tagger",1000,0.,25.);
  TH1F* hCFlav_PosTagger_0pu  = new TH1F("hCFlav_PosTagger_0pu","Tagger",1000,0.,25.);
  TH1F* hCFlav_PosTagger_GE8pu= new TH1F("hCFlav_PosTagger_GE8pu","Tagger",1000,0.,25.);
  TH1F* hCFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hCFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hCFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };
  TH1F* hCFlav_PosTag = new TH1F("hCFlav_PosTag","Tag(+)",50,0.,25.);
  TH1F* hCFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hCFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hCFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hCFlav_NegTag_JetPV[i]   = new TH1F("hCFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hCFlav_NegTag_JetPt[i]   = new TH2F("hCFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hCFlav_NegTag_JetEta[i]  = new TH1F("hCFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);

    hCFlav_PosTag_JetPV[i]   = new TH1F("hCFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hCFlav_PosTag_JetPt[i]   = new TH2F("hCFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hCFlav_PosTag_JetEta[i]  = new TH1F("hCFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }

//**********************************
// b-jets
//**********************************
  TH1F* hBFlav_All_JetPV      = new TH1F("hBFlav_All_JetPV","#PV",12,0.5,24.5);
  TH2F* hBFlav_All_JetPt      = new TH2F("hBFlav_All_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hBFlav_All_JetEta     = new TH1F("hBFlav_All_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hBFlav_JetPU          = new TH1F("hBFlav_JetPU","#PU",12,0.5,24.5);
  TH1F* hBFlav_JetPV          = new TH1F("hBFlav_JetPV","#PV",12,0.5,24.5);
  TH2F* hBFlav_JetPt          = new TH2F("hBFlav_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hBFlav_JetPt_etaLT12  = new TH1F("hBFlav_JetPt_etaLT12","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hBFlav_JetPt_etaGT12  = new TH1F("hBFlav_JetPt_etaGT12","pt(jet)",NPtBin, minPtBin, maxPtBin);
  TH1F* hBFlav_JetEta         = new TH1F("hBFlav_JetEta","|#eta(jet)|",24,0.,2.4);
  TH1F* hBFlav_PosTagger      = new TH1F("hBFlav_PosTagger","Tagger",1000,0.,25.);
  TH2F* hBFlav_Tagger         = new TH2F("hBFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hBFlav_deltaR         = new TH1F("hBFlav_deltaR","#Delta R gluon split",50,0.,1.);
  TH2F* hBFlav_JetPt_0pu      = new TH2F("hBFlav_JetPt_0pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH2F* hBFlav_JetPt_GE8pu    = new TH2F("hBFlav_JetPt_GE8pu","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
  TH1F* hBFlav_PosTagger_0pu  = new TH1F("hBFlav_PosTagger_0pu","Tagger",1000,0.,25.);
  TH1F* hBFlav_PosTagger_GE8pu= new TH1F("hBFlav_PosTagger_GE8pu","Tagger",1000,0.,25.);
  TH1F* hBFlav_NegTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hBFlav_NegTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hBFlav_NegTag_JetEta[3] = {NULL, NULL, NULL };

  TH1F* hBFlav_PosTag= new TH1F("hBFlav_PosTag","Tag(+)",50,0.,25.);
  TH1F* hBFlav_PosTag_JetPU[3] = {NULL, NULL, NULL };
  TH1F* hBFlav_PosTag_JetPV[3] = {NULL, NULL, NULL };
  TH2F* hBFlav_PosTag_JetPt[3] = {NULL, NULL, NULL };
  TH1F* hBFlav_PosTag_JetEta[3] = {NULL, NULL, NULL };

  for( int i = 0; i < 3; ++i )
  {
    hBFlav_NegTag_JetPV[i]   = new TH1F("hBFlav_NegTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hBFlav_NegTag_JetPt[i]   = new TH2F("hBFlav_NegTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hBFlav_NegTag_JetEta[i]  = new TH1F("hBFlav_NegTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
    hBFlav_PosTag_JetPU[i]   = new TH1F("hBFlav_PosTag"+TagLevel[i]+"_JetPU","#PU",12,0.5,24.5);
    hBFlav_PosTag_JetPV[i]   = new TH1F("hBFlav_PosTag"+TagLevel[i]+"_JetPV","#PV",12,0.5,24.5);
    hBFlav_PosTag_JetPt[i]   = new TH2F("hBFlav_PosTag"+TagLevel[i]+"_JetPt","pt(jet)",NPtBin, minPtBin, maxPtBin,NEtaBin, minEtaBin, maxEtaBin);
    hBFlav_PosTag_JetEta[i]  = new TH1F("hBFlav_PosTag"+TagLevel[i]+"_JetEta","|#eta(jet)|",24,0.,2.4);
  }

//**********************************
// No flavour
//**********************************
  TH2F* hNoFlav_Tagger        = new TH2F("hNoFlav_Tagger","Tagger",100,-25.,25.,NEtaBin, minEtaBin, maxEtaBin);

//--------- Start use of vector of histos
  //-1° Eta binning
  int NhEta = 3; 
  float fhEtaMax[3] = { 0.8, 1.6, 2.4 };
  float fhEtaMin[3];
  TString St_hEtaMin[3];
  TString St_hEtaMax[3];
  for(int i = 0; i < NhEta; ++i)
  {
    if( i == 0 )
    {
      St_hEtaMin[i] = Form("%1.1f",0.0);
      St_hEtaMin[i].Replace(1,1,"_");
      fhEtaMin[i] = 0.;
    }
    else
    {
      St_hEtaMin[i] = Form("%1.1f",fhEtaMax[i-1]);
      St_hEtaMin[i].Replace(1,1,"_");
      fhEtaMin[i] = fhEtaMax[i-1];
    }
    St_hEtaMax[i] = Form("%1.1f",fhEtaMax[i]);
    St_hEtaMax[i].Replace(1,1,"_");
  } 
  //-2° Pt binning
  /* For december 2013 comparison with ttH
  int NhPt = 6;
  // float fhPtMax[6] = { 30, 40, 60,100,160, 1000}; // For december 2013 comparison with ttH
  float fhPtMin[6];
  */
  // For december 2015 comparison with ttH
  int NhPt = 4;
  float fhPtMax[4] = { 30, 40, 60,1000}; // For december 2013 comparison with ttH
  float fhPtMin[4];
  TString St_hPtMin[4];
  TString St_hPtMax[4];
  for(int i = 0; i < NhPt; ++i)
  {
    if( i == 0 )
    {
      St_hPtMin[i] = Form("%4.0f",0.0);
      St_hPtMin[i].ReplaceAll(" ","0");
      fhPtMin[i] = 20;
    }
    else
    {
      St_hPtMin[i] = Form("%4.0f",fhPtMax[i-1]);
      St_hPtMin[i].ReplaceAll(" ","0");
      fhPtMin[i] = fhPtMax[i-1];
    }
    St_hPtMax[i] = Form("%4.0f",fhPtMax[i]);
    St_hPtMax[i].ReplaceAll(" ","0");
  } 
  //-3° Histograms Pt binning
  TH1F *hData_JetPtBinned_nPV[4];
  //-4° Histograms for Eta and Pt binning
  TH1F *h_AllFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_UDSFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_CFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_BFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_GluonFlavour_All_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_K0s_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_K0s_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Lam_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Lam_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Fak_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Fak_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Gam_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Gam_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Oth_IntTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Oth_IntTag[3][4]; //[NhEta][NhPt];

  TH1F *h_AllFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_UDSFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_CFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_BFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_GluonFlavour_All_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_K0s_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_K0s_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Lam_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Lam_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Fak_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Fak_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Gam_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Gam_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Oth_NegTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Oth_NegTag[3][4]; //[NhEta][NhPt];

  TH1F *h_AllFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_UDSFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_CFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_BFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_GluonFlavour_All_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_K0s_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_K0s_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Lam_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Lam_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Fak_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Fak_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Gam_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Gam_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_AllFlavour_Oth_PosTag[3][4]; //[NhEta][NhPt];
  TH1F *h_LightFlavour_Oth_PosTag[3][4]; //[NhEta][NhPt];
/*
double csvbins[] = { -10.0, 0.0, 0.040, 0.080, 0.120, 0.160, 0.200, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.752, 0.825, 0.898, 0.949, 1.010 };
TH1F *h = new TH1F("h","test",17,csvbins);
*/
  for( int iPt = 0; iPt < NhPt; ++iPt )
  {
    TString hName = "hData_nPV_PtMax_";
    hName += St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
    hData_JetPtBinned_nPV[iPt] = new TH1F(hName,hName,35,0.5,35.5);
  }

  int Ncsvbins = 17;
  double csvbins[18] = { -1.0, 0.0, 0.040, 0.080, 0.120, 0.160, 0.200, 0.244, 0.331, 0.418, 0.505, 0.592, 0.679, 0.752, 0.825, 0.898, 0.949, 1.010 };
  for( int i = 0; i < Ncsvbins+1; ++i ) csvbins[i]*=25.;

  for( int iEta = 0; iEta < NhEta; ++iEta )
  {
    for( int iPt = 0; iPt < NhPt; ++iPt )
    {
      TString hName = "IntTag_AllFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_All_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_All_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_UDSFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_UDSFlavour_All_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_CFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_CFlavour_All_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_BFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_BFlavour_All_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_GluonFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_GluonFlavour_All_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);

      hName = "IntTag_AllFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_K0s_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_K0s_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);

      hName = "IntTag_AllFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Lam_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Lam_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_AllFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Fak_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Fak_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_AllFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Gam_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Gam_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_AllFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Oth_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "IntTag_LightFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Oth_IntTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      //-------
      hName = "NegTag_AllFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_All_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_All_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_UDSFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_UDSFlavour_All_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_CFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_CFlavour_All_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_BFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_BFlavour_All_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_GluonFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_GluonFlavour_All_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);

      hName = "NegTag_AllFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_K0s_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_K0s_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);

      hName = "NegTag_AllFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Lam_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Lam_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_AllFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Fak_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Fak_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_AllFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Gam_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Gam_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_AllFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Oth_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "NegTag_LightFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Oth_NegTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      //-------
      hName = "PosTag_AllFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_All_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_All_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_UDSFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_UDSFlavour_All_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_CFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_CFlavour_All_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_BFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_BFlavour_All_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_GluonFlavour_AllPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_GluonFlavour_All_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);

      hName = "PosTag_AllFlavour_K0sPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_K0s_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_K0sPart_Eta"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_K0s_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);

      hName = "PosTag_AllFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Lam_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_LamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Lam_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_AllFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Fak_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_FakPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Fak_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_AllFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Gam_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_GamPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Gam_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_AllFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_AllFlavour_Oth_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
      hName = "PosTag_LightFlavour_OthPart_Eta_"+St_hEtaMin[iEta]+"to"+St_hEtaMax[iEta]+"_Pt_"+St_hPtMin[iPt]+"to"+St_hPtMax[iPt];
      h_LightFlavour_Oth_PosTag[iEta][iPt] = new TH1F(hName,hName,Ncsvbins,csvbins);
    }
  }

  // Pt Binning for nselTrack reweighting
  unsigned int N_ptbin_nseltrack = 18;
  float fhseltrackPtMin[18] = { 20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1000};
  float fhseltrackPtMax[18];
  TString St_fhseltrackPtMin[18];
  TString St_fhseltrackPtMax[18];
  for(int i = 0; i < N_ptbin_nseltrack; ++i)
  {
    if( i == (N_ptbin_nseltrack-1) )
    {
      fhseltrackPtMax[i] = 10000;
    }
    else
    {
      fhseltrackPtMax[i] = fhseltrackPtMin[i+1];
    }
    St_fhseltrackPtMin[i] = Form("%4.0f",fhseltrackPtMin[i]);
    St_fhseltrackPtMin[i].ReplaceAll(" ","0");

    St_fhseltrackPtMax[i] = Form("%4.0f",fhseltrackPtMax[i]);
    St_fhseltrackPtMax[i].ReplaceAll(" ","0");
  } 
  TH1F *hData_tracks[18]; //[N_ptbin_nseltrack];
  for( int iPt = 0; iPt < N_ptbin_nseltrack; ++iPt )
  {
    TString hName = "hData_tracks_Pt_"+St_fhseltrackPtMin[iPt]+"to"+St_fhseltrackPtMax[iPt];
    hData_tracks[iPt] = new TH1F(hName,hName,51, -0.5, 50.5);
  }

  TH1F* hData_All_nseltracks      = new TH1F("hData_All_nseltracks","nb. of tracks",51,-0.5,50.5); // No reweighting
  TH1F* hData_nseltracks      = new TH1F("hData_nseltracks","nb. of tracks",51,-0.5,50.5);     // PU and pthat reweighting
  TH1F* hData_W_nseltracks      = new TH1F("hData_W_nseltracks","nb. of tracks",51,-0.5,50.5);     // selTracks, PU and pthat reweighting
 //PVH 


 
///////////////////////////////////////////////////////////////////

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  //TEMP
  //nentries = 10000 ;
  std::cout << "Total Entries : " << nentries << std::endl;
  
  Long64_t nbytes = 0, nb = 0;

  bool TagPos[3] = {false, false, false}, TagNeg[3] = {false, false, false}, Veto = false;  
  float varpos, varneg, varPos;
  int numjet = -1, njets = -1;
  int ntagjet = 0;
  int allevents = 0;
  
  int e1 = int(1e+1);
  int e2 = int(1e+2);
  int e3 = int(1e+3);
  int e4 = int(1e+4);
  int e5 = int(1e+5);
  int e6 = int(1e+6);
  int e7 = int(1e+7);
  int e8 = int(1e+8);

//@@
  int itest = 0;
//@@

// Weight in pthat: QCD Fall11 prod
  float ww0 = 1.;  // pthat weight
  float ww = 1.;   // pthat and PU weight and trigger scale factor

// Read the weights related to the trigger prescales for data:
  LTANA::PS ps;
  ps.init(1);
// Get the trigger name from IntCut (to be used ig getPS...)
  TString PS_ToUse = "";
  if ( IntCut ==  40 ) PS_ToUse = "PFJet40";
  else if ( IntCut ==  60 ) PS_ToUse = "PFJet60";
  else if ( IntCut ==  80 ) PS_ToUse = "PFJet80";
  else if ( IntCut ==  140 ) PS_ToUse = "PFJet140";
  else if ( IntCut ==  200 ) PS_ToUse = "PFJet200";
  else if ( IntCut ==  260 ) PS_ToUse = "PFJet260";
  else if ( IntCut ==  320 ) PS_ToUse = "PFJet320";
  else if ( IntCut ==  400 ) PS_ToUse = "PFJet400";
  else 
  {
    cout << "Could not find proper prescale file for IntCut = " << IntCut << "  .Will use trigger prescale of 1" << endl;
  }

// Read pthat weight from file:
  float pthatMin[200], pthatMax[200], pthatweight[200];
  int maxLines = 200;
  ReadpthatWeight(pthatMin, pthatMax, pthatweight, maxLines, weightPthat_file);

// Read PU weight from file:
  float WeightPU[60], WeightPUmin[60], WeightPUmax[60];
  for (int k = 0; k < 60; k++) WeightPU[k] = 1.;

  int maxPU = 60;
  if( weightPU_file != "inJetAna" )
  {
    if ( ReadPUWeight(WeightPU, maxPU, weightPU_file) == 0 ) {
      cout << "\n Fail to read PU weight " << weightPU_file << endl ;
      return;  // If not a success, stop immediately.
    }
  }
  else
  {
    // MC pu reweighting
    if ( Year == 2011 )
    {
      if ( BTagMu )
      {
// WeightPUmin is computed from 2011A alone
// WeightPUmax is computed from 2011A+B but without runs 165088-167913
        WeightPU[0] = 0.42814;
        WeightPU[1] = 0.852011;
        WeightPU[2] = 1.23782;
        WeightPU[3] = 1.50857;
        WeightPU[4] = 1.67743;
        WeightPU[5] = 1.71335;
        WeightPU[6] = 1.63117;
        WeightPU[7] = 1.4899;
        WeightPU[8] = 1.31817;
        WeightPU[9] = 1.13156;
        WeightPU[10] = 0.959846;
        WeightPU[11] = 0.793424;
        WeightPU[12] = 0.658323;
        WeightPU[13] = 0.536642;
        WeightPU[14] = 0.442476;
        WeightPU[15] = 0.364995;
        WeightPU[16] = 0.301626;
        WeightPU[17] = 0.245887;
        WeightPU[18] = 0.199188;
        WeightPU[19] = 0.158801;
        WeightPU[20] = 0.125015;
        WeightPU[21] = 0.0965054;
        WeightPU[22] = 0.0743175;
        WeightPU[23] = 0.0552261;
        WeightPU[24] = 0.0410448;
        WeightPU[25] = 0.0300012;
        WeightPU[26] = 0.0216022;
        WeightPU[27] = 0.01554;
        WeightPU[28] = 0.0109019;
        WeightPU[29] = 0.0077661;
        WeightPU[30] = 0.00527189;
        WeightPU[31] = 0.00367071;
        WeightPU[32] = 0.00247471;
        WeightPU[33] = 0.0016687;
        WeightPU[34] = 0.00223246;
        WeightPUmin[0] = 0.508648;
        WeightPUmin[1] = 1.00906;
        WeightPUmin[2] = 1.45877;
        WeightPUmin[3] = 1.76461;
        WeightPUmin[4] = 1.94019;
        WeightPUmin[5] = 1.9487;
        WeightPUmin[6] = 1.80937;
        WeightPUmin[7] = 1.59264;
        WeightPUmin[8] = 1.33497;
        WeightPUmin[9] = 1.06071;
        WeightPUmin[10] = 0.80788;
        WeightPUmin[11] = 0.577805;
        WeightPUmin[12] = 0.397835;
        WeightPUmin[13] = 0.257752;
        WeightPUmin[14] = 0.162172;
        WeightPUmin[15] = 0.0985342;
        WeightPUmin[16] = 0.0582882;
        WeightPUmin[17] = 0.0332803;
        WeightPUmin[18] = 0.0185821;
        WeightPUmin[19] = 0.0100938;
        WeightPUmin[20] = 0.0053696;
        WeightPUmin[21] = 0.00278412;
        WeightPUmin[22] = 0.00143351;
        WeightPUmin[23] = 0.000709667;
        WeightPUmin[24] = 0.000350299;
        WeightPUmin[25] = 0.00016959;
        WeightPUmin[26] = 8.06722e-05;
        WeightPUmin[27] = 3.82429e-05;
        WeightPUmin[28] = 1.76357e-05;
        WeightPUmin[29] = 8.23732e-06;
        WeightPUmin[30] = 3.65711e-06;
        WeightPUmin[31] = 1.66109e-06;
        WeightPUmin[32] = 7.28646e-07;
        WeightPUmin[33] = 3.18864e-07;
        WeightPUmin[34] = 2.06626e-07;
        WeightPUmax[0] = 0.321623;
        WeightPUmax[1] = 0.633854;
        WeightPUmax[2] = 0.929572;
        WeightPUmax[3] = 1.1604;
        WeightPUmax[4] = 1.34017;
        WeightPUmax[5] = 1.44133;
        WeightPUmax[6] = 1.46356;
        WeightPUmax[7] = 1.44183;
        WeightPUmax[8] = 1.3872;
        WeightPUmax[9] = 1.30012;
        WeightPUmax[10] = 1.20261;
        WeightPUmax[11] = 1.07688;
        WeightPUmax[12] = 0.957174;
        WeightPUmax[13] = 0.824407;
        WeightPUmax[14] = 0.708206;
        WeightPUmax[15] = 0.601179;
        WeightPUmax[16] = 0.506317;
        WeightPUmax[17] = 0.417744;
        WeightPUmax[18] = 0.340911;
        WeightPUmax[19] = 0.272991;
        WeightPUmax[20] = 0.215471;
        WeightPUmax[21] = 0.166586;
        WeightPUmax[22] = 0.128399;
        WeightPUmax[23] = 0.0954633;
        WeightPUmax[24] = 0.0709704;
        WeightPUmax[25] = 0.0518837;
        WeightPUmax[26] = 0.0373623;
        WeightPUmax[27] = 0.0268788;
        WeightPUmax[28] = 0.0188572;
        WeightPUmax[29] = 0.0134334;
        WeightPUmax[30] = 0.00911912;
        WeightPUmax[31] = 0.00634949;
        WeightPUmax[32] = 0.0042807;
        WeightPUmax[33] = 0.00288649;
        WeightPUmax[34] = 0.00107161;
      }
      else if ( IntCut == 300 )
      {
// WeightPUmin is computed from 2011A+B but without runs 178098-180252
// WeightPUmax is computed from 2011A+B but without runs 165088-167913
        WeightPU[0] = 0.0236327;
        WeightPU[1] = 0.186895;
        WeightPU[2] = 0.411244;
        WeightPU[3] = 0.715173;
        WeightPU[4] = 0.976423;
        WeightPU[5] = 1.17936;
        WeightPU[6] = 1.29761;
        WeightPU[7] = 1.36739;
        WeightPU[8] = 1.39617;
        WeightPU[9] = 1.40456;
        WeightPU[10] = 1.41265;
        WeightPU[11] = 1.41307;
        WeightPU[12] = 1.45094;
        WeightPU[13] = 1.46775;
        WeightPU[14] = 1.50175;
        WeightPU[15] = 1.52329;
        WeightPU[16] = 1.5588;
        WeightPU[17] = 1.58469;
        WeightPU[18] = 1.59661;
        WeightPU[19] = 1.57766;
        WeightPU[20] = 1.56339;
        WeightPU[21] = 1.52986;
        WeightPU[22] = 1.50984;
        WeightPU[23] = 1.41909;
        WeightPU[24] = 1.38554;
        WeightPU[25] = 1.34946;
        WeightPU[26] = 1.28489;
        WeightPU[27] = 1.15493;
        WeightPU[28] = 1.03597;
        WeightPU[29] = 0.984955;
        WeightPU[30] = 0.876273;
        WeightPU[31] = 0.772276;
        WeightPU[32] = 0.720553;
        WeightPU[33] = 0.707666;
        WeightPU[34] = 0.378344;
        WeightPUmin[0] = 0.293038;
        WeightPUmin[1] = 0.58781;
        WeightPUmin[2] = 0.86736;
        WeightPUmin[3] = 1.08442;
        WeightPUmin[4] = 1.2539;
        WeightPUmin[5] = 1.35408;
        WeightPUmin[6] = 1.38846;
        WeightPUmin[7] = 1.39185;
        WeightPUmin[8] = 1.37383;
        WeightPUmin[9] = 1.33004;
        WeightPUmin[10] = 1.27531;
        WeightPUmin[11] = 1.18245;
        WeightPUmin[12] = 1.08192;
        WeightPUmin[13] = 0.950148;
        WeightPUmin[14] = 0.822697;
        WeightPUmin[15] = 0.695634;
        WeightPUmin[16] = 0.577288;
        WeightPUmin[17] = 0.465037;
        WeightPUmin[18] = 0.367816;
        WeightPUmin[19] = 0.283842;
        WeightPUmin[20] = 0.214971;
        WeightPUmin[21] = 0.158955;
        WeightPUmin[22] = 0.116885;
        WeightPUmin[23] = 0.0827463;
        WeightPUmin[24] = 0.0584831;
        WeightPUmin[25] = 0.0405946;
        WeightPUmin[26] = 0.0277256;
        WeightPUmin[27] = 0.0188999;
        WeightPUmin[28] = 0.0125536;
        WeightPUmin[29] = 0.00846052;
        WeightPUmin[30] = 0.00543007;
        WeightPUmin[31] = 0.00357258;
        WeightPUmin[32] = 0.00227472;
        WeightPUmin[33] = 0.00144797;
        WeightPUmin[34] = 0.000481535;
        WeightPUmax[0] = 0.138424;
        WeightPUmax[1] = 0.28206;
        WeightPUmax[2] = 0.435051;
        WeightPUmax[3] = 0.582968;
        WeightPUmax[4] = 0.739521;
        WeightPUmax[5] = 0.892853;
        WeightPUmax[6] = 1.03605;
        WeightPUmax[7] = 1.18049;
        WeightPUmax[8] = 1.3209;
        WeightPUmax[9] = 1.43841;
        WeightPUmax[10] = 1.53504;
        WeightPUmax[11] = 1.56652;
        WeightPUmax[12] = 1.56158;
        WeightPUmax[13] = 1.48143;
        WeightPUmax[14] = 1.37646;
        WeightPUmax[15] = 1.24272;
        WeightPUmax[16] = 1.09721;
        WeightPUmax[17] = 0.937959;
        WeightPUmax[18] = 0.785909;
        WeightPUmax[19] = 0.641764;
        WeightPUmax[20] = 0.513977;
        WeightPUmax[21] = 0.401756;
        WeightPUmax[22] = 0.312279;
        WeightPUmax[23] = 0.233716;
        WeightPUmax[24] = 0.174679;
        WeightPUmax[25] = 0.128265;
        WeightPUmax[26] = 0.0927113;
        WeightPUmax[27] = 0.0669142;
        WeightPUmax[28] = 0.0470794;
        WeightPUmax[29] = 0.0336246;
        WeightPUmax[30] = 0.0228792;
        WeightPUmax[31] = 0.0159644;
        WeightPUmax[32] = 0.0107841;
        WeightPUmax[33] = 0.00728486;
        WeightPUmax[34] = 0.00271264;
      }
      else
      {
        WeightPU[0] = 0.0428918;
        WeightPU[1] = 0.335376;
        WeightPU[2] = 0.721533;
        WeightPU[3] = 1.20833;
        WeightPU[4] = 1.55725;
        WeightPU[5] = 1.7342;
        WeightPU[6] = 1.71538;
        WeightPU[7] = 1.58616;
        WeightPU[8] = 1.39386;
        WeightPU[9] = 1.19426;
        WeightPU[10] = 1.02425;
        WeightPU[11] = 0.88516;
        WeightPU[12] = 0.802558;
        WeightPU[13] = 0.735522;
        WeightPU[14] = 0.698884;
        WeightPU[15] = 0.67222;
        WeightPU[16] = 0.662791;
        WeightPU[17] = 0.656711;
        WeightPU[18] = 0.649984;
        WeightPU[19] = 0.634296;
        WeightPU[20] = 0.622947;
        WeightPU[21] = 0.605549;
        WeightPU[22] = 0.594589;
        WeightPU[23] = 0.556591;
        WeightPU[24] = 0.541617;
        WeightPU[25] = 0.526013;
        WeightPU[26] = 0.499605;
        WeightPU[27] = 0.448078;
        WeightPU[28] = 0.40112;
        WeightPU[29] = 0.380669;
        WeightPU[30] = 0.338094;
        WeightPU[31] = 0.297503;
        WeightPU[32] = 0.277176;
        WeightPU[33] = 0.271852;
        WeightPU[34] = 0.145013;
      }
    }
  }

// Read NSelTrack weight from files:

  unsigned int N_SelTracks = 52;
  float WeightNS[18*52];
  ReadNSelTracksWeight(WeightNS, N_SelTracks, N_ptbin_nseltrack, weightNS_file, St_fhseltrackPtMin, St_fhseltrackPtMax);

// Read JSon file:

  unsigned int N_MaxRun = 100;
  unsigned int N_MaxLumiRangePerRun = 1000;
  unsigned int RunSize;
  unsigned int JS_Run[100];
  unsigned int JS_Run_NLumi[100];
  unsigned int JS_LumiMin[100*1000];
  unsigned int JS_LumiMax[100*1000];

  int JSON_RunN = ReadJSon(N_MaxRun, N_MaxLumiRangePerRun, RunSize, JS_Run, JS_Run_NLumi, JS_LumiMin, JS_LumiMax, JSONFile);
  if ( JSON_RunN < 0 )
  {
    cout << "Failed to read JSon file with return value " << JSON_RunN << endl;
    return;
  } 
  
// 

////////////////
// Event loop //
////////////////

  bool testJSONSelect = true;

  int N_repEvts = 0;
  //for (Long64_t jentry=0; jentry<100; jentry++)
  for (Long64_t jentry=0; jentry<nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
   
//----------------------------------------------------------------
//Start selection (remove bad runs, bad PtVspthat...
// Bad run in promptreco data (temporary before reprocessing)

    if (period == "7p7to9p2invfb" && Run <= 276097) {
      //TEMP check period choice
      //cout << "\n =================> Debug: Period == 7p7to9p2invfb and Run <= 276097 skip event " << endl ;
      continue ; 
    }
    
    hData_nPV_noWeight_allEvent->Fill(nPV) ;

    if ( Run > 1 )	// Selections for data files
    {
      // Select using JSon file
      if ( JSON_RunN > 0 )
      {
        bool notInRange = true;
        for( int i = 0; i <= RunSize && notInRange ; ++i)
        {
          for( int j = 0; j < JS_Run_NLumi[i]; ++j )
          {
            if( testJSONSelect )
            {
              cout << Run << " : " << LumiBlock << " =?= " << JS_Run[i] << " : " << JS_LumiMin[i*1000+j] << " success = " << ( Run == JS_Run[i] && ( LumiBlock >= JS_LumiMin[i*1000+j] && LumiBlock <= JS_LumiMax[i*1000+j] ) ) << endl;
              testJSONSelect = false;
            }
            if( Run == JS_Run[i] && ( LumiBlock >= JS_LumiMin[i*1000+j] && LumiBlock <= JS_LumiMax[i*1000+j] ) )
            {
              notInRange = false;
              break;
            }
          }
        }
        if( notInRange ) {
          //cout << "\n===============> Debug: Do JSON check, not in range, event skipped" ;
          continue; 		// Go to next event
        }
      }

      // RunII Certification: 7.3 pb-1 certifiés
      if( TrigType == "2015_Cert" ) 
      {
        if (! ( ( Run == 251244 && ( (LumiBlock >= 85 && LumiBlock <= 86) ||
                                     (LumiBlock >= 88 && LumiBlock <= 93) ||
                                     (LumiBlock >= 96 && LumiBlock <=121) ||
                                     (LumiBlock >=123 && LumiBlock <=156) ||
                                     (LumiBlock >=158 && LumiBlock <=428) ||
                                     (LumiBlock >=430 && LumiBlock <=442) ) ) ||
                ( Run == 251251 && ( (LumiBlock >=  1 && LumiBlock <= 31) ||
                                     (LumiBlock >= 33 && LumiBlock <= 97) ||
                                     (LumiBlock >= 99 && LumiBlock <=167) ) ) ||
                ( Run == 251252 && ( (LumiBlock >=  1 && LumiBlock <=283) ||
                                     (LumiBlock >=285 && LumiBlock <=505) ||
                                     (LumiBlock >=507 && LumiBlock <=554) ) ) ||
                ( Run == 251561 && ( (LumiBlock >=  1 && LumiBlock <= 94) ) )
        ) ) continue;
      }

      // RunII looser run selection corresponding to 40 pb-1
      if( TrigType == "2015_LooserSel" )
      {
        if (! ( ( Run == 251244 && ( (LumiBlock >= 85 && LumiBlock <= 86) ||
                                     (LumiBlock >= 88 && LumiBlock <= 93) ||
                                     (LumiBlock >= 96 && LumiBlock <=121) ||
                                     (LumiBlock >=123 && LumiBlock <=156) ||
                                     (LumiBlock >=158 && LumiBlock <=428) ||
                                     (LumiBlock >=430 && LumiBlock <=442) ) ) ||
                ( Run == 251251 && ( (LumiBlock >=  1 && LumiBlock <= 31) ||
                                     (LumiBlock >= 33 && LumiBlock <= 97) ||
                                     (LumiBlock >= 99 && LumiBlock <=167) ) ) ||
                ( Run == 251252 && ( (LumiBlock >=  1 && LumiBlock <=283) ||
                                     (LumiBlock >=285 && LumiBlock <=505) ||
                                     (LumiBlock >=507 && LumiBlock <=554) ) ) ||
                ( Run == 251561 && ( (LumiBlock >=  1 && LumiBlock <= 94) ) ) ||
                ( Run == 251562 && ( (LumiBlock >=  1 && LumiBlock <=691) ) ) ||
                ( Run == 251636 && ( (LumiBlock >=  1 && LumiBlock <=  3) ||
                                     (LumiBlock >= 45 && LumiBlock <= 45) ) ) ||
                ( Run == 251637 && ( (LumiBlock >=  1 && LumiBlock <=  6) ) ) ||
                ( Run == 251640 && ( (LumiBlock >=  1 && LumiBlock <=  7) ||
                                     (LumiBlock >=123 && LumiBlock <=123) ) ) ||
                ( Run == 251643 && ( (LumiBlock >=  1 && LumiBlock <=606) ) ) ||
                ( Run == 251717 && ( (LumiBlock >= 14 && LumiBlock <= 21) ) ) ||
                ( Run == 251718 && ( (LumiBlock >=  1 && LumiBlock <=288) ) ) ||
                ( Run == 251721 && ( (LumiBlock >=  1 && LumiBlock <= 28) ||
                                     (LumiBlock >=304 && LumiBlock <=304) ) ) ||
                ( Run == 251864 && ( (LumiBlock >=  1 && LumiBlock <= 17) ||
                                     (LumiBlock >=199 && LumiBlock <=199) ) ) ||
                ( Run == 251883 && ( (LumiBlock >= 56 && LumiBlock <= 60) ||
                                     (LumiBlock >= 62 && LumiBlock <=437) ) )
        ) ) continue;
      }
      // Filter laser events which show up with large nJet or large Jet_pt[0]
      if( nJet > 0 )
      {
        if( nJet > 16 || Jet_pt[0] > 1600 ) continue;
      }
    } 
    else if( Run < 0 )
    {
      bool Skip_StrangeEvts = false;
      if( Skip_StrangeEvts )
      {
        if( pthat <  30. ) { if( Jet_pt[0] > 130. ) continue; }		// Go to next event
        else if( pthat <  50. ) { if( Jet_pt[0] > 130. ) continue; }	// Go to next event
        else if( pthat <  80. ) { if( Jet_pt[0] > 180. ) continue; }	// Go to next event
        else if( pthat <  120. ) { if( Jet_pt[0] > 250. ) continue; }	// Go to next event
        else if( pthat <  170. ) { if( Jet_pt[0] > 310. ) continue; }	// Go to next event
        else if( pthat <  300. ) { if( Jet_pt[0] > 460. ) continue; }	// Go to next event

        // Add special trick to remove some kind of repetitive event... Or the result of a bug:
        //   some events with more than 3 jets with values around 382, 267 and 135
        //                    and corresponding eta values around -0.6, -2.12, -1.27
        //                    should be skipped
        if( nJet > 0 && pthat < 200 )
        {
          if ( Jet_pt[0] > 1000 )
          {
            N_repEvts++;
            cout << "Warning: one more strange event type: " << N_repEvts << endl;
            cout << "Event = " << Evt << " \tnJet = " << nJet << endl;
            cout << "Jet0 with pt = " << Jet_pt[0] << " \t and eta = " << Jet_eta[0] << endl;
            continue;		// Go to next event
          }
        }
        if( nJet >= 3 && pthat < 200 )
        {
          bool BadJet0 = false;
          bool BadJet1 = false;
          bool BadJet2 = false;
          Double_t BadJet0_pt, BadJet0_eta;
          Double_t BadJet1_pt, BadJet1_eta;
          Double_t BadJet2_pt, BadJet2_eta;
          for( int iJet = 0; iJet < nJet; ++iJet )
          {
            if( Jet_pt[iJet] > 360 && Jet_pt[iJet] < 400 && Jet_eta[iJet] > -0.70 && Jet_eta[iJet] < -0.60 )
            {
              BadJet0_pt = Jet_pt[iJet];
              BadJet0_eta = Jet_eta[iJet];
              BadJet0 = true;
              cout << "Warning0" << endl;
            }
            if( Jet_pt[iJet] > 250 && Jet_pt[iJet] < 305 && Jet_eta[iJet] > -2.20 && Jet_eta[iJet] < -2.10 )
            {
              BadJet1_pt = Jet_pt[iJet];
              BadJet1_eta = Jet_eta[iJet];
              BadJet1 = true;
              cout << "Warning1" << endl;
            }
            if( Jet_pt[iJet] > 125 && Jet_pt[iJet] < 145 && Jet_eta[iJet] > -1.40 && Jet_eta[iJet] < -1.20 )
            {
              BadJet2_pt = Jet_pt[iJet];
              BadJet2_eta = Jet_eta[iJet];
              BadJet2 = true;
              cout << "Warning2" << endl;
            }
          }
          if( BadJet0 && BadJet1 && BadJet2 )
          {
            N_repEvts++;
            cout << "Warning: one more repetitive event type: " << N_repEvts << endl;
            cout << "Event = " << Evt << " \tnJet = " << nJet << endl;
            cout << "Jet0 with pt = " << BadJet0_pt << " \t and eta = " << BadJet0_eta << endl;
            cout << "Jet1 with pt = " << BadJet1_pt << " \t and eta = " << BadJet1_eta << endl;
            cout << "Jet2 with pt = " << BadJet2_pt << " \t and eta = " << BadJet2_eta << endl;
            continue;		// Go to next event
          }
        }
      }
    } // End of else if( Run < 0 )

    // Filter against monsters in MC
    bool MONSTER = false;
    if ( pthat > 0 && pthat < 120. ) {
      for (int ijet = 0; ijet < nJet; ijet++) {
        if ( Jet_pt[ijet] > pthat * 7. ) MONSTER = true;
      }
    }
    if ( MONSTER ) continue;


//End of bad runs, bad PtVspthat selection 
//----------------------------------------------------------------

// pthat reweighting procedure
//$$
    int npv = nPV, npu = nPU; 
    if ( truePU ) npu = nPUtrue; 
    if ( npv >= 40 ) npv = 40;
    if ( npu >= maxPU-1 ) npu = maxPU-1;


    ww = 0; 
    ww0 = 0.;
    if ( Run >= 0 ) 
    {
      ww0 = 1;    // No reweighting for Data
      ww = 1 ;
      //TEMP turn of trigger scale weight
      //ww = ps.getPSWeight((std::string) PS_ToUse,(int) Run,(int) LumiBlock);    // valid upto PF400 trigger, above this return 1 
    } else {
      bool inPthatBin = false;
      for( int i = 0; i < maxLines; ++i )
      {
        if( pthat > pthatMin[i] && pthat < pthatMax[i] )
        {
          ww0 = pthatweight[i];
          inPthatBin = true;
          break;
        }
      }
      if( !inPthatBin ) continue;
      if( weightPU_file.Contains("_nPV") )      // Then weight according to nPV
      {
        if( nPV > 0 && nPV < 60 ) {
          //TEMP
          ww = ww0 * WeightPU[nPV]; //FIXME change from nPV - 1 to nPV
          //cout << "Weight PU is: " << nPV << "  " << WeightPU[nPV] << endl ;
        }
        else ww = 0;
      } else {					// Weight according to npu
        ww = ww0 * WeightPU[npu]; 
      }
    }
//$$
    hAllFlav_pthatAll->Fill( pthat, 1. );
    hAllFlav_pthat->Fill( pthat , ww0 );
    hData_nPV_test->Fill(nPV, ww) ;
    hData_nPV_noWeight->Fill(npv) ;

//----------------------------------------------------------------
// Trigger selection
//$$

    int fj_getEntry_status = FJ->GetEntry(jentry);
    //cout << "\n Get fat jet tree status: " << fj_getEntry_status ;

    if( FJ->FatJetInfo_nJet > 0 ) hBeforeTrig_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);

// Filter event if maxpt < given cut || > given cut
//TEMP
    //if( FJ->FatJetInfo_nJet > 0 )
    //{
      //if( FJ->FatJetInfo_Jet_pt[0] < minCutJetPtMax ) continue; // Go to next event (Cut associated to the trigger)
      //if( FJ->FatJetInfo_Jet_pt[0] > maxCutJetPtMax ) continue; // Go to next event
    //}
 

    bool Jet30  = false, Jet60  = false, Jet150 = false, Jet190 = false, Jet240 = false;
    bool Jet40  = false, Jet80  = false, Jet140 = false;
    bool Jet200 = false, Jet260 = false, Jet320 = false, Jet400 = false, Jet450 = false, Jet500 = false;
    bool Jet20  = false, Jet70  = false, Jet110 = false, Jet300 = false;

    int njet20 = 0, njet40 = 0, njet60 = 0, njet70 = 0, njet110 = 0, njet300 = 0;
    int njet80 = 0, njet140 = 0, njet200 = 0, njet260 = 0, njet320 = 0;
    int njet400 = 0, njet450 = 0, njet500 = 0;

    int triggerIdx = 0, bitIdx = 0;

    if ( !BTagMu ) {
      if ( Year == 2011 ) {
        triggerIdx = 1;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
        {
          Jet30  = true;
          if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig30_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
        }
   
        triggerIdx = 4;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
        {
          Jet60  = true;
          if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig60_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
        }
   
        triggerIdx = 6;
        bitIdx = int(triggerIdx/32);
        if( Run > 0 )
        {
          if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
          {
            Jet80  = true;
            if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig80_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
          }
        } else {
          if( FJ->FatJetInfo_nJet > 0 && FJ->FatJetInfo_Jet_pt[0] >= 90 ) // Used to simulate the inexisting MC trigger. Also used for Data for compatibility
          {
            Jet80  = true;
            if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig80_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
          }
        }

        triggerIdx =  9;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
        {
          Jet110 = true;
          if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig110_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
        }
   
        triggerIdx = 11;
        bitIdx = int(triggerIdx/32);
        if( Run > 0 )
        {
          if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
          {
            Jet150 = true;
            if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig150_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
          }
        } else {
          if( FJ->FatJetInfo_nJet > 0 && FJ->FatJetInfo_Jet_pt[0] >= 70 ) // Used to simulate the inexisting MC trigger. Also used for Data for compatibility
          {
            Jet150 = true;
            if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig150_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
          }
        }
   
        triggerIdx = 14;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
        {
          Jet190 = true;
          if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig190_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
        }
   
        triggerIdx = 16;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
        {
          Jet240 = true;
          if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig240_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
        }
   
        triggerIdx = 18;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) )
        {
          Jet300 = true;
          if( FJ->FatJetInfo_nJet > 0 ) hAfterTrig300_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0);
        }
      }
      if ( Year == 2012 ) {
        triggerIdx = 2;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;
   
        triggerIdx = 7;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet80  = true;
   
        triggerIdx = 12;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet140 = true;
   
        triggerIdx = 15;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet200 = true;
 
        triggerIdx = 17;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet260 = true;

        triggerIdx = 19;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet320 = true;
   
        triggerIdx = 21;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) {
          Jet400 = true;
          hAfterTrig400_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0) ;
        }

        triggerIdx = 87;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet60 = true;

        triggerIdx = 88;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet450 = true;

        triggerIdx = 89;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) {
          Jet500 = true;
          hAfterTrig500_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0) ;
        }

        for (int ijet = 0; ijet < FJ->FatJetInfo_nJet; ijet++) {
          float residual = 1.;
          float ptjet = FJ->FatJetInfo_Jet_pt[ijet] * residual;
          float etajet = fabs(FJ->FatJetInfo_Jet_eta[ijet]);

          if ( ptjet >  50. ) njet40++;
          if ( ptjet >  70. ) njet60++;
          if ( ptjet > 100. ) njet80++;
          if ( ptjet > 160. ) njet140++;
          if ( ptjet > 220. ) njet200++;
          if ( ptjet > 300. ) njet260++;
          if ( ptjet > 360. ) njet320++;
          //TEMP
          //if ( ptjet > 420. ) njet400++;
          //if ( ptjet > 470. ) njet450++;
          //if ( ptjet > 520. ) njet500++;
          if ( ptjet > 450. ) njet400++;
          if ( ptjet > 500. ) njet450++;
          if ( ptjet > 550. ) njet500++;
        }
        if ( njet40  < 1 ) Jet40  = false;
        if ( njet60  < 1 ) Jet60  = false;
        if ( njet80  < 1 ) Jet80  = false;
        if ( njet140 < 1 ) Jet140 = false;
        if ( njet200 < 1 ) Jet200 = false;
        if ( njet260 < 1 ) Jet260 = false;
        if ( njet320 < 1 ) Jet320 = false;
        if ( njet400 < 1 ) Jet400 = false;
        if ( njet450 < 1 ) Jet450 = false;
        if ( njet500 < 1 ) Jet500 = false;
        
      }

      //////////////////////////////////////////////////////
      
      if ( Year == 2016 ) {

        triggerIdx = 0;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;
   
        triggerIdx = 1;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet60  = true;
        
        triggerIdx = 2;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet80  = true;
   
        triggerIdx = 3;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet140 = true;
   
        triggerIdx = 4;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet200 = true;
 
        triggerIdx = 5;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) {
          Jet260 = true;
          hAfterTrig260_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0) ;
        }

        triggerIdx = 6;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) {
          Jet320 = true;
          hAfterTrig320_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0) ;
        }
   
        triggerIdx = 7;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) {
          Jet400 = true;
          hAfterTrig400_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0) ;
        }

        triggerIdx = 8;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet450 = true;

        triggerIdx = 9;
        bitIdx = int(triggerIdx/32);
        if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) {
          Jet500 = true;
          hAfterTrig500_MaxJetPt->Fill(FJ->FatJetInfo_Jet_pt[0],ww0) ;
        }

        for (int ijet = 0; ijet < FJ->FatJetInfo_nJet; ijet++) {
          float residual = 1.;
          float ptjet = FJ->FatJetInfo_Jet_pt[ijet] * residual;
          float etajet = fabs(FJ->FatJetInfo_Jet_eta[ijet]);
          float maxJetEta = 2.4 ;
          if ( ptjet >  50. && etajet < maxJetEta) njet40++;
          if ( ptjet >  70. && etajet < maxJetEta) njet60++;
          if ( ptjet > 100. && etajet < maxJetEta) njet80++;
          if ( ptjet > 160. && etajet < maxJetEta) njet140++;
          if ( ptjet > 220. && etajet < maxJetEta) njet200++;
          if ( ptjet > 300. && etajet < maxJetEta) njet260++;
          if ( ptjet > 360. && etajet < maxJetEta) njet320++;
          //TEMP
          //if ( ptjet > 420. && etajet < maxJetEta) njet400++;
          //if ( ptjet > 470. && etajet < maxJetEta) njet450++;
          //if ( ptjet > 520. && etajet < maxJetEta) njet500++;
          if ( ptjet > 450. && etajet < maxJetEta) njet400++;
          if ( ptjet > 500. && etajet < maxJetEta) njet450++;
          if ( ptjet > 550. && etajet < maxJetEta) njet500++;
        }

        if ( njet40  < 1 ) Jet40  = false;
        if ( njet60  < 1 ) Jet60  = false;
        if ( njet80  < 1 ) Jet80  = false;
        if ( njet140 < 1 ) Jet140 = false;
        if ( njet200 < 1 ) Jet200 = false;
        if ( njet260 < 1 ) Jet260 = false;
        if ( njet320 < 1 ) Jet320 = false;
        if ( njet400 < 1 ) Jet400 = false;
        if ( njet450 < 1 ) Jet450 = false;
        if ( njet500 < 1 ) Jet500 = false;
        
      }


      //////////////////////////////////////////////////////

      //TEMP bypass trigger, select event based on number of jet passing certain threshold
      if (Year == 0) {
        Jet40  = false;
        Jet60  = false;
        Jet80  = false;
        Jet140 = false;
        Jet200 = false;
        Jet260 = false;
        Jet320 = false;
        Jet400 = false;
        Jet450 = false;
        Jet500 = false;
        for (int ijet = 0; ijet < FJ->FatJetInfo_nJet; ijet++) {
          float residual = 1.;
          float ptjet = FJ->FatJetInfo_Jet_pt[ijet] * residual;
          float etajet = fabs(FJ->FatJetInfo_Jet_eta[ijet]);

          if ( ptjet >  50. ) njet40++;
          if ( ptjet >  70. ) njet60++;
          if ( ptjet > 100. ) njet80++;
          if ( ptjet > 160. ) njet140++;
          if ( ptjet > 220. ) njet200++;
          if ( ptjet > 300. ) njet260++;
          if ( ptjet > 360. ) njet320++;
          //TEMP
          //if ( ptjet > 420. ) njet400++;
          //if ( ptjet > 470. ) njet450++;
          //if ( ptjet > 520. ) njet500++;
          if ( ptjet > 450. ) njet400++;
          if ( ptjet > 500. ) njet450++;
          if ( ptjet > 550. ) njet500++;
        }

        if ( njet40  > 0 ) Jet40  = true;
        if ( njet60  > 0 ) Jet60  = true;
        if ( njet80  > 0 ) Jet80  = true;
        if ( njet140 > 0 ) Jet140 = true;
        if ( njet200 > 0 ) Jet200 = true;
        if ( njet260 > 0 ) Jet260 = true;
        if ( njet320 > 0 ) Jet320 = true;
        if ( njet400 > 0 ) Jet400 = true;
        if ( njet450 > 0 ) Jet450 = true;
        if ( njet500 > 0 ) Jet500 = true;

      }
    }

    else if ( BTagMu ) {

      triggerIdx = 35;
      bitIdx = int(triggerIdx/32);
      if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet20  = true;
   
      triggerIdx = 41;
      bitIdx = int(triggerIdx/32);
      if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet40  = true;
   
      triggerIdx = 44;
      bitIdx = int(triggerIdx/32);
      if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet70  = true;
   
      triggerIdx = 47;
      bitIdx = int(triggerIdx/32);
      if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet110 = true;
   
      triggerIdx = 50;
      bitIdx = int(triggerIdx/32);
      if ( BitTrigger[bitIdx] & ( 1 << (triggerIdx - bitIdx*32) ) ) Jet300 = true;
   
      for (int ijet = 0; ijet < FJ->FatJetInfo_nJet; ijet++) {
        float residual = 1.;
        float ptjet = FJ->FatJetInfo_Jet_pt[ijet] * residual;
        float etajet = fabs(FJ->FatJetInfo_Jet_eta[ijet]);
   
        //hData_Trigger_JetPt->Fill( ptjet , ww );
        //if ( Jet20 )  hData_Trig20_JetPt->Fill(  ptjet , ww );
        //if ( Jet40 )  hData_Trig40_JetPt->Fill(  ptjet , ww );
        //if ( Jet70 )  hData_Trig70_JetPt->Fill(  ptjet , ww );
        //if ( Jet110 ) hData_Trig110_JetPt->Fill( ptjet , ww );
        //if ( Jet300 ) hData_Trig300_JetPt->Fill( ptjet , ww );
   
        if ( ptjet > 20. )  njet20++;
        if ( ptjet > 50. )  njet40++;
        if ( ptjet > 80. )  njet70++;
        if ( ptjet > 120. ) njet110++;
        if ( ptjet > 320. ) njet300++;
      }
      if ( Year == 2012 ) {   
        if ( njet20  < 1 ) Jet20  = false;
        if ( njet40  < 1 ) Jet40  = false;
        if ( njet70  < 1 ) Jet70  = false;
        if ( njet110 < 1 ) Jet110 = false;
        if ( njet300 < 1 ) Jet300 = false;
      }
    }
         
   // Trigger selection
    if ( !BTagMu && !TopAna ) {
      if ( IntCut ==  30 && !Jet30 )  continue;
      if ( IntCut ==  40 && !Jet40 )  continue;
      if ( IntCut ==  60 && !Jet60 )  continue;
      if ( IntCut ==  80 && !Jet80 )  continue;
      //if ( IntCut == 80 && FJ->FatJetInfo_Jet_pt[0] < minCutJetPtMax ) continue;		// Go to next event
      if ( IntCut == 110 && !Jet110 ) continue;
      if ( IntCut == 140 && !Jet140 ) continue;
      if ( IntCut == 150 && !Jet150 ) continue;
      //if ( IntCut == 150 && FJ->FatJetInfo_Jet_pt[0] < minCutJetPtMax ) continue;		// Go to next event
      if ( IntCut == 190 && !Jet190 ) continue;
      if ( IntCut == 200 && !Jet200 ) continue;
      if ( IntCut == 240 && !Jet240 ) continue;
      if ( IntCut == 260 && !Jet260 ) {
        //cout << "\n =========> Debug: Fail 260 trigger, event skipped" ;
        continue;
      }
      if ( IntCut == 300 && !Jet300 ) {
        //cout << "\n =========> Debug: Fail 300 trigger, event skipped" ;
        continue;
      }
      if ( IntCut == 320 && !Jet320 ) {
        //cout << "\n =========> Debug: Fail 320 trigger, event skipped" ;
        continue;
      }
      if ( IntCut == 400 && !Jet400 ) {
        //cout << "\n =========> Debug: Fail 400 trigger, event skipped" ;
        continue;
      }
      if ( IntCut == 450 && !Jet450 ) continue;
      if ( IntCut == 500 && !Jet500 ) {
        //cout << "\n =========> Debug: Fail 500 trigger, event skipped" ;
        continue;
      }
      if ( IntCut == 1 && Year == 2011
           && !Jet30 && !Jet60 && !Jet80 && !Jet110 && !Jet150 
                     && !Jet190 && !Jet240 && !Jet300 ) continue;
      if ( IntCut == 1 && Year == 2012
           && !Jet40 && !Jet80 && !Jet140
                     && !Jet200 && !Jet260 && !Jet320 ) continue;
    }
    else if ( BTagMu ) {
      if ( IntCut == 20  && !Jet20 )  continue;
      if ( IntCut == 40  && !Jet40 )  continue;
      if ( IntCut == 70  && !Jet70 )  continue;
      if ( IntCut == 110 && !Jet110 ) continue;
      if ( IntCut == 300 && !Jet300 ) continue;
      if ( IntCut == 1 && Year == 2011
           && !Jet20 && !Jet40 && !Jet70 && !Jet110 ) continue;
      if ( IntCut == 1 && Year == 2012
           && !Jet20 && !Jet40 && !Jet70 && !Jet110 && !Jet300 ) continue;
    }
         
   
   
   
// End of trigger selection
//----------------------------------------------------------------
//

//----------------------------------------------------------------
// Eta and ptmax selection
    int njet30 = 0;
    float maxJetPt;
    maxJetPt = 0;
    for (int ijet = 0; ijet < FJ->FatJetInfo_nJet; ijet++)
    {
      float residual = 1.;
      if ( Run > 0 ) residual = FJ->FatJetInfo_Jet_residual[ijet];
      float ptjet = FJ->FatJetInfo_Jet_pt[ijet]; // * residual;
      float etajet = fabs(FJ->FatJetInfo_Jet_eta[ijet]);
      // Only used info on jet between EtaMin and EtaMax to make the selection... A bit strange ?
      if ( ptjet > 30. ) njet30++;

      if( ptjet > maxJetPt ) maxJetPt = ptjet;
    }
    //if( maxJetPt < minCutJetPtMax ) continue;				// Go to next event
// End of eta and ptmax selection
//----------------------------------------------------------------


    hAllFlav_pthatTrig->Fill( pthat, 1. );
      
    hData_All_nPV->Fill( npv , ww0 );		// Only pthat reweighted
    for( int iPt = 0; iPt < NhPt; ++iPt )
    {
      if( maxJetPt >= fhPtMin[iPt] && maxJetPt < fhPtMax[iPt] )
      {
        hData_JetPtBinned_nPV[iPt]->Fill( npv , ww0 );
      }
    }

    hData_nPV->Fill( npv , ww );		// pthat and PU reweighted and trigger scale for data
    hAllFlav_All_nPU->Fill( npu , ww0 );
    hAllFlav_nPU->Fill( npu , ww ); //change from nPU to npu
    
//FIXME
    if ( npv >= 40 ) npv = 40; //change from 24 to 40
    if ( npu >= 40 ) npu = 40;

//------------------------------------------------
// Loop on Fat jets
    float tmp_ww = ww;
    njets = FJ->FatJetInfo_nJet;

    for (int ijet = 0; ijet < njets; ijet++) {
      if ( ijet == 0 ) {
        allevents++;
        if ( allevents%10000 == 0 ) std::cout << "events : " << allevents << std::endl ;

        hData_All_NJets->Fill( njets , ww );
        hData_NJets->Fill( numjet , ww );
        hData_NJet30->Fill( njet30 , ww );
   
        numjet = 0;
        ntagjet = 0;
      }
//TEMP no jet mass cut
      if( FJ->FatJetInfo_Jet_massPruned[ijet] > 10000 || FJ->FatJetInfo_Jet_massPruned[ijet] < 50 ) continue; // Skip this jet
      if( fabs(FJ->FatJetInfo_Jet_eta[ijet]) > 2.4 ) continue; // Skip this jet
      if( FJ->FatJetInfo_Jet_pt[ijet] < minCutJetPtMax ) continue; // Go to next jet (Cut associated to the trigger)

      // Loop on "SoftDropSubJet" of this FatJet
      for( int iPSJet = FJ->FatJetInfo_Jet_nFirstSJ_SoftDrop[ijet] ; iPSJet < FJ->FatJetInfo_Jet_nLastSJ_SoftDrop[ijet]; ++ iPSJet )
      {
        int ntracks = FJ->SoftDropSubJetInfo_Jet_ntracks[iPSJet];
        int nseltracks = FJ->SoftDropSubJetInfo_Jet_nseltracks[iPSJet];
        float residual = 1.;
        if ( Run > 0 ) residual = FJ->SoftDropSubJetInfo_Jet_residual[iPSJet];
        float ptjet = FJ->SoftDropSubJetInfo_Jet_pt[iPSJet];// * residual;
        float etajet = fabs(FJ->SoftDropSubJetInfo_Jet_eta[iPSJet]);
        float phijet = FJ->SoftDropSubJetInfo_Jet_phi[iPSJet];
        int flavour = abs( FJ->SoftDropSubJetInfo_Jet_flavour[iPSJet] );
        if( Run > 0 )
        {
          flavour = -1;
        }
        else
        {      
          if ( (flavour != 4) && (flavour != 5)  && (flavour != 21) ) flavour = 1; // This is used if Gluon info is not stored
          // the flavour inforation is not always complete. Anyway, c and b flavour are stored as 4 and 5.
          // if gluon flavor is kept, it is number 21
          // light flavour should contain all the rest (u,d, s)  so the line above should always be correct
        }


        hData_All_nseltracks->Fill(nseltracks, 1.);   // No reweighting
        hData_nseltracks->Fill(nseltracks, tmp_ww);   // pthat and PU reweighting
        hData_W_nseltracks->Fill(nseltracks, ww); // pthat, PU and nseltracks reweighting

//*********************************
// Jet selection

        if ( Run < 0 && FJ->SoftDropSubJetInfo_Jet_genpt[iPSJet] < 8 ) continue;			// Go to next jet
        if ( !(ptjet >= PtMin && ptjet < PtMax) ) continue;		// Go to next jet
        if ( ptjet > maxPtBin) ptjet = maxPtBin -0.01;

        numjet++;

//*********************************

        for( int i = 0; i < 3; ++i )
        {
          TagPos[i] = false;
          TagNeg[i] = false;
        }
        Veto = false;

        varpos = -1000.;
        varneg = -1000.;
        varPos = -1000.;

// Track history
        int cat = 0, catP = 0, catN = 0;
        int cat1 = 0, cat1P = 0, cat1N = 0;
        int cat2 = 0, cat2P = 0, cat2N = 0;

        bool GamP = false, GamN = false, GamPN = false;
        bool K0sP = false, K0sN = false, K0sPN = false;
        bool LamP = false, LamN = false, LamPN = false;
        bool FakP = false, FakN = false, FakPN = false;
        bool OthP = false, OthN = false, OthPN = false;
   
//***************************
// Tagger

// Jet Probability
        if ( TagName == "JP" )
        {
          if ( FJ->SoftDropSubJetInfo_Jet_ProbaP[iPSJet] > 0 ) varPos = 10.*FJ->SoftDropSubJetInfo_Jet_ProbaP[iPSJet];
          if ( FJ->SoftDropSubJetInfo_Jet_ProbaP[iPSJet] > 0 ) varpos = 10.*FJ->SoftDropSubJetInfo_Jet_ProbaP[iPSJet];
          if ( FJ->SoftDropSubJetInfo_Jet_ProbaN[iPSJet] > 0 ) varneg = 10.*FJ->SoftDropSubJetInfo_Jet_ProbaN[iPSJet];
          for( int i = 0; i < 3; ++i )
          {
            if ( FJ->SoftDropSubJetInfo_Jet_ProbaP[iPSJet] > TagCut[i] ) TagPos[i] = true;
            if ( FJ->SoftDropSubJetInfo_Jet_ProbaN[iPSJet] > TagCut[i] ) TagNeg[i] = true;
          }
          cat = FJ->SoftDropSubJetInfo_Jet_histJet[iPSJet];
        }
// Combined CSVv2
        else if ( TagName == "CSVv2" ) {
          if ( FJ->SoftDropSubJetInfo_Jet_CombIVF[iPSJet] > 0. )   varPos = 25.*FJ->SoftDropSubJetInfo_Jet_CombIVF[iPSJet];
          if ( FJ->SoftDropSubJetInfo_Jet_CombIVF_P[iPSJet] > 0. ) varpos = 25.*FJ->SoftDropSubJetInfo_Jet_CombIVF_P[iPSJet];
          if ( FJ->SoftDropSubJetInfo_Jet_CombIVF_N[iPSJet] > 0. ) varneg = 25.*FJ->SoftDropSubJetInfo_Jet_CombIVF_N[iPSJet];
          for( int i = 0; i < 3; ++i )
          {
            if ( FJ->SoftDropSubJetInfo_Jet_CombIVF[iPSJet] > TagCut[i] )   TagPos[i] = true;
            if ( FJ->SoftDropSubJetInfo_Jet_CombIVF_N[iPSJet] > TagCut[i] ) TagNeg[i] = true;
          }
          cat = FJ->SoftDropSubJetInfo_Jet_histJet[iPSJet];
        }
// Combined cMVAv2
        else if ( TagName == "cMVAv2" ) {
          if ( FJ->SoftDropSubJetInfo_Jet_cMVAv2[iPSJet] > -1. )   varPos = 12.5*(1+FJ->SoftDropSubJetInfo_Jet_cMVAv2[iPSJet]);
          if ( FJ->SoftDropSubJetInfo_Jet_cMVAv2P[iPSJet] > -1. ) varpos = 12.5*(1+FJ->SoftDropSubJetInfo_Jet_cMVAv2P[iPSJet]);
          if ( FJ->SoftDropSubJetInfo_Jet_cMVAv2N[iPSJet] > -1. ) varneg = 12.5*(1+FJ->SoftDropSubJetInfo_Jet_cMVAv2N[iPSJet]);
          for( int i = 0; i < 3; ++i )
          {
            if ( FJ->SoftDropSubJetInfo_Jet_cMVAv2[iPSJet] > TagCut[i] )   TagPos[i] = true;
            if ( FJ->SoftDropSubJetInfo_Jet_cMVAv2N[iPSJet] > TagCut[i] ) TagNeg[i] = true;
          }
          cat = FJ->SoftDropSubJetInfo_Jet_histJet[iPSJet];
        }
//$$
        if ( varPos > 24.9 ) varPos = 24.9; 
        if ( varpos > 24.9 ) varpos = 24.9; 
        if ( varneg > 24.9 ) varneg = 24.9; 

        if ( varPos < 0.1 ) varPos = 25*-0.5; 
        if ( varpos < 0.1 ) varpos = 25*-0.5; 
        if ( varneg < 0.1 ) varneg = 25*-0.5; 
//$$

// Track History information
        if      ( cat%e1 == 1		    || cat%e1 == 3)	   catP = 1; // BWeakDecay
        else if ((cat%e2 >= e1 && cat%e2 < 2*e1) || cat%e2 >= 3*e1 ) catP = 2; // CWeakDecay
        else if ((cat%e3 >= e2 && cat%e3 < 2*e2) || cat%e3 >= 3*e2 ) catP = 3; // TauDecay
        else if ((cat%e4 >= e3 && cat%e4 < 2*e3) || cat%e4 >= 3*e3 ) catP = 4; // Conversion
        else if ((cat%e5 >= e4 && cat%e5 < 2*e4) || cat%e5 >= 3*e4 ) catP = 5; // KsDecay
        else if ((cat%e6 >= e5 && cat%e6 < 2*e5) || cat%e6 >= 3*e5 ) catP = 6; // LambdaDecay
        else if ((cat%e7 >= e6 && cat%e7 < 2*e6) || cat%e7 >= 3*e6 ) catP = 7; // Interaction
        else if ((cat%e8 >= e7 && cat%e8 < 2*e7) || cat%e8 >= 3*e7 ) catP = 8; // Fake
        if      ( cat%e1 >= 2 )    catN = 1; // BWeakDecay
        else if ( cat%e2 >= 2*e1 ) catN = 2; // CWeakDecay
        else if ( cat%e3 >= 2*e2 ) catN = 3; // TauDecay
        else if ( cat%e4 >= 2*e3 ) catN = 4; // Conversion
        else if ( cat%e5 >= 2*e4 ) catN = 5; // KsDecay
        else if ( cat%e6 >= 2*e5 ) catN = 6; // LambdaDecay
        else if ( cat%e7 >= 2*e6 ) catN = 7; // Interaction
        else if ( cat%e8 >= 2*e7 ) catN = 8; // Fake

        if ( TagName == "TCHE" || TagName == "TCHP" ) {
          cat1 = FJ->SoftDropSubJetInfo_Jet_hist1[iPSJet];
          if	     ( cat1%e1 == 1			|| cat1%e1 == 3)     cat1P = 1;
          else if ((cat1%e2 >= e1 && cat1%e2 < 2*e1) || cat1%e2 >= 3*e1 ) cat1P = 2;
          else if ((cat1%e3 >= e2 && cat1%e3 < 2*e2) || cat1%e3 >= 3*e2 ) cat1P = 3;
          else if ((cat1%e4 >= e3 && cat1%e4 < 2*e3) || cat1%e4 >= 3*e3 ) cat1P = 4;
          else if ((cat1%e5 >= e4 && cat1%e5 < 2*e4) || cat1%e5 >= 3*e4 ) cat1P = 5;
          else if ((cat1%e6 >= e5 && cat1%e6 < 2*e5) || cat1%e6 >= 3*e5 ) cat1P = 6;
          else if ((cat1%e7 >= e6 && cat1%e7 < 2*e6) || cat1%e7 >= 3*e6 ) cat1P = 7;
          else if ((cat1%e8 >= e7 && cat1%e8 < 2*e7) || cat1%e8 >= 3*e7 ) cat1P = 8;
          else if ((cat1    >= e8 && cat1    < 2*e8) || cat1    >= 3*e8 ) cat1P = 9;
          if	     ( cat1%e1 >= 2 )    cat1N = 1;
          else if ( cat1%e2 >= 2*e1 ) cat1N = 2;
          else if ( cat1%e3 >= 2*e2 ) cat1N = 3;
          else if ( cat1%e4 >= 2*e3 ) cat1N = 4;
          else if ( cat1%e5 >= 2*e4 ) cat1N = 5;
          else if ( cat1%e6 >= 2*e5 ) cat1N = 6;
          else if ( cat1%e7 >= 2*e6 ) cat1N = 7;
          else if ( cat1%e8 >= 2*e7 ) cat1N = 8;
          else if ( cat1    >= 2*e8 ) cat1N = 9;
          if ( cat1P > 0 && ( cat1P < catP || catP == 0 ) ) catP = cat1P; 
          if ( cat1N > 0 && ( cat1N < catN || catN == 0 ) ) catN = cat1N; 
        }

        if ( TagName == "TCHP" || TagName == "SSVHP" ) {
          if ( TagName == "TCHP" ) cat2 = FJ->SoftDropSubJetInfo_Jet_hist2[iPSJet];
          if ( TagName == "SSVHP" ) cat2 = FJ->SoftDropSubJetInfo_Jet_histSvx[iPSJet];
          if	     ( cat2%e1 == 1			|| cat2%e1 == 3)     cat2P = 1;
          else if ((cat2%e2 >= e1 && cat2%e2 < 2*e1) || cat2%e2 >= 3*e1 ) cat2P = 2;
          else if ((cat2%e3 >= e2 && cat2%e3 < 2*e2) || cat2%e3 >= 3*e2 ) cat2P = 3;
          else if ((cat2%e4 >= e3 && cat2%e4 < 2*e3) || cat2%e4 >= 3*e3 ) cat2P = 4;
          else if ((cat2%e5 >= e4 && cat2%e5 < 2*e4) || cat2%e5 >= 3*e4 ) cat2P = 5;
          else if ((cat2%e6 >= e5 && cat2%e6 < 2*e5) || cat2%e6 >= 3*e5 ) cat2P = 6;
          else if ((cat2%e7 >= e6 && cat2%e7 < 2*e6) || cat2%e7 >= 3*e6 ) cat2P = 7;
          else if ((cat2%e8 >= e7 && cat2%e8 < 2*e7) || cat2%e8 >= 3*e7 ) cat2P = 8;
          else if ((cat2    >= e8 && cat2    < 2*e8) || cat2    >= 3*e8 ) cat2P = 9;
          if	     ( cat2%e1 >= 2 )    cat2N = 1;
          else if ( cat2%e2 >= 2*e1 ) cat2N = 2;
          else if ( cat2%e3 >= 2*e2 ) cat2N = 3;
          else if ( cat2%e4 >= 2*e3 ) cat2N = 4;
          else if ( cat2%e5 >= 2*e4 ) cat2N = 5;
          else if ( cat2%e6 >= 2*e5 ) cat2N = 6;
          else if ( cat2%e7 >= 2*e6 ) cat2N = 7;
          else if ( cat2%e8 >= 2*e7 ) cat2N = 8;
          else if ( cat2    >= 2*e8 ) cat2N = 9;
          if ( TagName == "TCHP" ) {
            if ( cat2P > 0 && ( cat2P < catP || catP == 0 ) ) catP = cat2P; 
            if ( cat2N > 0 && ( cat2N < catN || catN == 0 ) ) catN = cat2N;
          } 
          if ( TagName == "SSVHP" ) {
            if ( FJ->SoftDropSubJetInfo_Jet_Svx[iPSJet] > 1.  )  catP = cat2P; 
            if ( FJ->SoftDropSubJetInfo_Jet_SvxN[iPSJet] < -1. ) catN = cat2N;
          } 
        }

        if      ( catP == 4 || catP == 7 ) GamP = true;
        else if ( catP == 5 )              K0sP = true;
        else if ( catP == 6 )              LamP = true;
        else if ( catP == 8 )              FakP = true;
        if      ( catP == 0 || catP == 9 ) OthP = true;
     
        if      ( catN == 4 || catN == 7 ) GamN = true;
        else if ( catN == 5 )              K0sN = true;
        else if ( catN == 6 )              LamN = true;
        else if ( catN == 8 )              FakN = true;
        if      ( catN == 0 || catN == 9 ) OthN = true;
  
        if      ( catP == 4 || catN == 4 || catP == 7 || catN == 7 ) GamPN = true;
        else if ( catP == 5 || catN == 5 )                           K0sPN = true;
        else if ( catP == 6 || catN == 6 )                           LamPN = true;
        else if ( catP == 8 || catN == 8 )                           FakPN = true;
        if      ( catP == 0 || catN == 0 || catP == 9 || catN == 9 ) OthPN = true;
  
//*********************************

        hData_All_tracks->Fill( ntracks , 1. );
        hData_All_JetPV->Fill( npv , 1. );
        hData_All_JetPt->Fill( ptjet, etajet , 1. );
        hData_All_JetEta->Fill( etajet , 1. );
        hData_All_JetRun->Fill( Run , 1. );

        hAllFlav_All_Flavour->Fill( flavour , 1. );

        if ( flavour == 1 || flavour == 21 ) {
          hLightFlav_All_JetPV->Fill( npv , 1. );
          hLightFlav_All_JetPt->Fill( ptjet, etajet , 1. );
          hLightFlav_All_JetEta->Fill( etajet , 1. );
        }

        if (flavour == 21) {				// In principle this will now be empty !
          hGluonFlav_All_JetPV->Fill( npv , 1. );
          hGluonFlav_All_JetPt->Fill( ptjet, etajet , 1. );
          hGluonFlav_All_JetEta->Fill( etajet , 1. );
        }
        else if (flavour == 1) {				// In principle, it will contain all udsg jets !
          hUDSFlav_All_JetPV->Fill( npv , 1. );
          hUDSFlav_All_JetPt->Fill( ptjet, etajet , 1. );
          hUDSFlav_All_JetEta->Fill( etajet , 1. );
        }
        else if (flavour == 4) {				// C jets
          hCFlav_All_JetPV->Fill( npv , 1. );
          hCFlav_All_JetPt->Fill( ptjet, etajet , 1. );
          hCFlav_All_JetEta->Fill( etajet , 1. );
        }
        else if (flavour == 5) {				// B jets
          hBFlav_All_JetPV->Fill( npv , 1. );
          hBFlav_All_JetPt->Fill( ptjet, etajet , 1. );
          hBFlav_All_JetEta->Fill( etajet , 1. );
        }

        //PVH nseltracks
        for( int iPt = 0; iPt < N_ptbin_nseltrack; ++iPt )
        {
          //TString hName = "hData_tracks_Pt_"+St_fhseltrackPtMin[iPt]+"to"+St_fhseltrackPtMax[iPt];
          if( ptjet >= fhseltrackPtMin[iPt] && ptjet < fhseltrackPtMax[iPt] ) hData_tracks[iPt]->Fill(nseltracks,tmp_ww);
        }
   
//*********************************
// Taggability
//$$
        if ( ntracks < NtrackMin ) continue;		// Go to next jet
//$$
        ntagjet++;

        hData_JetPV->Fill( npv , ww ) ;
        hData_JetPt->Fill( ptjet, etajet , ww );
        hData_JetEta->Fill( etajet , ww );
        hData_JetRun->Fill( Run , ww );
        if ( npv <= 5 ) {
          hData_JetPt_LE5pv->Fill( ptjet , ww );
          hData_JetEta_LE5pv->Fill( etajet , ww );
        } 
        else {
          hData_JetPt_GE6pv->Fill( ptjet , ww );
          hData_JetEta_GE6pv->Fill( etajet , ww );
        }
        if ( njet30 == 1 ) hData_1JetPt->Fill( ptjet , ww );
        if ( njet30 == 2 ) hData_2JetPt->Fill( ptjet , ww );
        if ( njet30 == 3 ) hData_3JetPt->Fill( ptjet , ww );
        if ( njet30 >= 4 ) hData_4JetPt->Fill( ptjet , ww );
  
        hAllFlav_Flavour->Fill( flavour , ww );
        if ( GamPN ) {
          hAllFlav_Gam_JetPV->Fill( npv , ww ) ;
          hAllFlav_Gam_JetPt->Fill( ptjet, etajet , ww );
          hAllFlav_Gam_JetEta->Fill( etajet , ww );
        }
        if ( K0sPN ) {
          hAllFlav_K0s_JetPV->Fill( npv , ww ) ;
          hAllFlav_K0s_JetPt->Fill( ptjet, etajet , ww );
          hAllFlav_K0s_JetEta->Fill( etajet , ww );
        }
        if ( LamPN ) {
          hAllFlav_Lam_JetPV->Fill( npv , ww ) ;
          hAllFlav_Lam_JetPt->Fill( ptjet, etajet , ww );
          hAllFlav_Lam_JetEta->Fill( etajet , ww );
        }
        if ( FakPN ) {
          hAllFlav_Fak_JetPV->Fill( npv , ww ) ;
          hAllFlav_Fak_JetPt->Fill( ptjet, etajet , ww );
          hAllFlav_Fak_JetEta->Fill( etajet , ww );
        }
     
        if ( flavour == 1 || flavour == 21 ) {
          hLightFlav_JetPU->Fill( npu , ww0 );
          hLightFlav_JetPV->Fill( npv , ww );
          hLightFlav_JetPt->Fill( ptjet, etajet , ww );
          hLightFlav_JetEta->Fill( etajet , ww );
          if ( nPU < 10 ) hLightFlav_JetPt_0pu->Fill( ptjet , ww );
          if ( nPU >= 20) hLightFlav_JetPt_GE8pu->Fill( ptjet , ww );
          if ( GamPN ) {
            hLightFlav_Gam_JetPV->Fill( npv , ww ) ;
            hLightFlav_Gam_JetPt->Fill( ptjet, etajet , ww );
            hLightFlav_Gam_JetEta->Fill( etajet , ww );
          }
          if ( K0sPN ) {
            hLightFlav_K0s_JetPV->Fill( npv , ww ) ;
            hLightFlav_K0s_JetPt->Fill( ptjet, etajet , ww );
            hLightFlav_K0s_JetEta->Fill( etajet , ww );
          }
          if ( LamPN ) {
            hLightFlav_Lam_JetPV->Fill( npv , ww ) ;
            hLightFlav_Lam_JetPt->Fill( ptjet, etajet , ww );
            hLightFlav_Lam_JetEta->Fill( etajet , ww );
          }
          if ( FakPN ) {
            hLightFlav_Fak_JetPV->Fill( npv , ww ) ;
            hLightFlav_Fak_JetPt->Fill( ptjet, etajet , ww );
            hLightFlav_Fak_JetEta->Fill( etajet , ww );
          }
          if ( OthPN ) {
            hLightFlav_Oth_JetPV->Fill( npv , ww ) ;
            hLightFlav_Oth_JetPt->Fill( ptjet, etajet , ww );
            hLightFlav_Oth_JetEta->Fill( etajet , ww );
          }
          if ( njet30 == 1 ) hLightFlav_1JetPt->Fill( ptjet , ww );
          if ( njet30 == 2 ) hLightFlav_2JetPt->Fill( ptjet , ww );
          if ( njet30 == 3 ) hLightFlav_3JetPt->Fill( ptjet , ww );
          if ( njet30 >= 4 ) hLightFlav_4JetPt->Fill( ptjet , ww );
        }

        if (flavour == 21) {
          hGluonFlav_JetPV->Fill( npv , ww );
          hGluonFlav_JetPt->Fill( ptjet, etajet , ww );
          hGluonFlav_JetEta->Fill( etajet , ww );
        }
        else if (flavour == 1) {
          hUDSFlav_JetPV->Fill( npv , ww );
          hUDSFlav_JetPt->Fill( ptjet, etajet , ww );
          hUDSFlav_JetEta->Fill( etajet , ww );
        }
        else if (flavour == 4) {
          hCFlav_JetPV->Fill( npv , ww );
          hCFlav_JetPt->Fill( ptjet, etajet , ww );
          hCFlav_JetEta->Fill( etajet , ww );
          if ( nPU < 10 ) hCFlav_JetPt_0pu->Fill( ptjet , ww );
          if ( nPU >= 20) hCFlav_JetPt_GE8pu->Fill( ptjet , ww );
        }
        else if (flavour == 5) {
          hBFlav_JetPU->Fill( npu , ww0 );
          hBFlav_JetPV->Fill( npv , ww );
          hBFlav_JetPt->Fill( ptjet, etajet , ww );
          if ( etajet < 1.2 ) hBFlav_JetPt_etaLT12->Fill( ptjet , ww );
          else                hBFlav_JetPt_etaGT12->Fill( ptjet , ww );
          hBFlav_JetEta->Fill( etajet , ww );
          if ( nPU < 10 ) hBFlav_JetPt_0pu->Fill( ptjet , ww );
          if ( nPU >= 20) hBFlav_JetPt_GE8pu->Fill( ptjet , ww );
        }
     
//*********************************
// Tagging

        for( int i = 0; i < 3; ++i )
        {
          if ( TagNeg[i] ) {
            hData_NegTag_JetPV[i]->Fill( npv , ww );
            hData_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hData_NegTag_JetEta[i]->Fill( etajet , ww );
            hData_NegTag_JetRun[i]->Fill( Run , ww );
            if ( npv <= 5 ) {
              hData_NegTag_JetPt_LE5pv[i]->Fill( ptjet , ww );
              hData_NegTag_JetEta_LE5pv[i]->Fill( etajet , ww );
            } 
            else {
              hData_NegTag_JetPt_GE6pv[i]->Fill( ptjet , ww );
              hData_NegTag_JetEta_GE6pv[i]->Fill( etajet , ww );
            }
            if ( njet30 == 1 ) hData_NegTag_1JetPt[i]->Fill( ptjet , ww );
            if ( njet30 == 2 ) hData_NegTag_2JetPt[i]->Fill( ptjet , ww );
            if ( njet30 == 3 ) hData_NegTag_3JetPt[i]->Fill( ptjet , ww );
            if ( njet30 >= 4 ) hData_NegTag_4JetPt[i]->Fill( ptjet , ww );
          } else break;
        }
        for( int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hData_PosTag_JetPV[i]->Fill( npv , ww );
            hData_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hData_PosTag_JetEta[i]->Fill( etajet , ww );
            hData_PosTag_JetRun[i]->Fill( Run , ww );
            if ( K0sPN ) {
              hAllFlav_K0s_PosTag_JetPV[i]->Fill( npv , ww );
              hAllFlav_K0s_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
              hAllFlav_K0s_PosTag_JetEta[i]->Fill( etajet , ww );
            }
            if ( LamPN ) {
              hAllFlav_Lam_PosTag_JetPV[i]->Fill( npv , ww );
              hAllFlav_Lam_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
              hAllFlav_Lam_PosTag_JetEta[i]->Fill( etajet , ww );
            }
          } else break;
        }
        if ( varneg > 0 )
        {
          hData_Tagger->Fill(-varneg , etajet , ww );
        }
        if ( varpos > 0 )
        {
          hData_Tagger->Fill( varpos , etajet , ww );
        }
        if ( varPos > 0 ) hData_PosTag->Fill( varPos , ww );
  
        for( int iEta = 0; iEta < NhEta; ++iEta )
        {
          for( int iPt = 0; iPt < NhPt; ++iPt )
          {
            if( etajet >= fhEtaMin[iEta] && etajet < fhEtaMax[iEta] && ptjet >= fhPtMin[iPt] && ptjet < fhPtMax[iPt] )
            {
              for( int iTagBin = 0; iTagBin <= h_AllFlavour_All_IntTag[iEta][iPt]->GetNbinsX()+1; ++iTagBin )
              {
                Double_t xTagBin = h_AllFlavour_All_IntTag[iEta][iPt]->GetBinCenter(iTagBin);
  
                h_AllFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( K0sN || K0sP ) h_AllFlavour_K0s_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( LamN || LamP ) h_AllFlavour_Lam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( FakN || FakP ) h_AllFlavour_Fak_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( GamN || GamP ) h_AllFlavour_Gam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( OthN || OthP ) h_AllFlavour_Oth_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                if( flavour == 1 || flavour == 21 ) // UDS and Gluon
                {
                  h_LightFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                  if( K0sN || K0sP ) h_LightFlavour_K0s_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                  if( LamN || LamP ) h_LightFlavour_Lam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                  if( FakN || FakP ) h_LightFlavour_Fak_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                  if( GamN || GamP ) h_LightFlavour_Gam_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                  if( OthN || OthP ) h_LightFlavour_Oth_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                }
                if( flavour == 1 ) // UDS only
                {
                  h_UDSFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                }
                if( flavour == 21 ) // Gluon only
                {
                  h_GluonFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                }
                if( flavour == 4 ) // C only
                {
                  h_CFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                }
                if( flavour == 5 ) // B only
                {
                  h_BFlavour_All_IntTag[iEta][iPt]->Fill(xTagBin,ww);
                }
              }
              if( varneg > 0 )
              {
                h_AllFlavour_All_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( K0sN ) h_AllFlavour_K0s_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( LamN ) h_AllFlavour_Lam_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( FakN ) h_AllFlavour_Fak_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( GamN ) h_AllFlavour_Gam_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( OthN ) h_AllFlavour_Oth_NegTag[iEta][iPt]->Fill( varneg, ww );
                if( flavour == 1 || flavour == 21 )
                {
                  h_LightFlavour_All_NegTag[iEta][iPt]->Fill( varneg, ww );
                  if( K0sN ) h_LightFlavour_K0s_NegTag[iEta][iPt]->Fill( varneg, ww );
                  if( LamN ) h_LightFlavour_Lam_NegTag[iEta][iPt]->Fill( varneg, ww );
                  if( FakN ) h_LightFlavour_Fak_NegTag[iEta][iPt]->Fill( varneg, ww );
                  if( GamN ) h_LightFlavour_Gam_NegTag[iEta][iPt]->Fill( varneg, ww );
                  if( OthN ) h_LightFlavour_Oth_NegTag[iEta][iPt]->Fill( varneg, ww );
                }
                if( flavour == 1 ) // UDS only
                {
                  h_UDSFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
                }
                if( flavour == 21 ) // Gluon only
                {
                  h_GluonFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
                }
                if( flavour == 4 ) // C only
                {
                  h_CFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
                }
                if( flavour == 5 ) // B only
                {
                  h_BFlavour_All_NegTag[iEta][iPt]->Fill(varneg,ww);
                }
              }
              if( 1 ) // Used to be: varPos > 0 )
              {
                h_AllFlavour_All_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( K0sP ) h_AllFlavour_K0s_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( LamP ) h_AllFlavour_Lam_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( FakP ) h_AllFlavour_Fak_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( GamP ) h_AllFlavour_Gam_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( OthP ) h_AllFlavour_Oth_PosTag[iEta][iPt]->Fill( varPos, ww );
                if( flavour == 1 || flavour == 21 )
                {
                  h_LightFlavour_All_PosTag[iEta][iPt]->Fill( varPos, ww );
                  if( K0sP ) h_LightFlavour_K0s_PosTag[iEta][iPt]->Fill( varPos, ww );
                  if( LamP ) h_LightFlavour_Lam_PosTag[iEta][iPt]->Fill( varPos, ww );
                  if( FakP ) h_LightFlavour_Fak_PosTag[iEta][iPt]->Fill( varPos, ww );
                  if( GamP ) h_LightFlavour_Gam_PosTag[iEta][iPt]->Fill( varPos, ww );
                  if( OthP ) h_LightFlavour_Oth_PosTag[iEta][iPt]->Fill( varPos, ww );
                }
                if( flavour == 1 ) // UDS only
                {
                  h_UDSFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
                }
                if( flavour == 21 ) // Gluon only
                {
                  h_GluonFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
                }
                if( flavour == 4 ) // C only
                {
                  h_CFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
                }
                if( flavour == 5 ) // B only
                {
                  h_BFlavour_All_PosTag[iEta][iPt]->Fill(varPos,ww);
                }
              }
            }
          }
        }
   
        if ( varneg > 0 ) {
          if      ( catN == 1 ) hAllFlav_Tagger_Bwd->Fill(-varneg , ww );
          else if ( catN == 2 ) hAllFlav_Tagger_Cwd->Fill(-varneg , ww );
          else if ( catN == 3 ) hAllFlav_Tagger_Tau->Fill(-varneg , ww );
          else if ( catN == 4 ) hAllFlav_Tagger_Gam->Fill(-varneg , ww );
          else if ( catN == 5 ) hAllFlav_Tagger_K0s->Fill(-varneg , ww );
          else if ( catN == 6 ) hAllFlav_Tagger_Lam->Fill(-varneg , ww );
          else if ( catN == 7 ) hAllFlav_Tagger_Int->Fill(-varneg , ww );
          else if ( catN == 8 ) hAllFlav_Tagger_Fak->Fill(-varneg , ww );
          else		   hAllFlav_Tagger_Oth->Fill(-varneg , ww );
        }
        if ( varpos > 0 ) {
          if      ( catP == 1 ) hAllFlav_Tagger_Bwd->Fill( varpos , ww );
          else if ( catP == 2 ) hAllFlav_Tagger_Cwd->Fill( varpos , ww );
          else if ( catP == 3 ) hAllFlav_Tagger_Tau->Fill( varpos , ww );
          else if ( catP == 4 ) hAllFlav_Tagger_Gam->Fill( varpos , ww );
          else if ( catP == 5 ) hAllFlav_Tagger_K0s->Fill( varpos , ww );
          else if ( catP == 6 ) hAllFlav_Tagger_Lam->Fill( varpos , ww );
          else if ( catP == 7 ) hAllFlav_Tagger_Int->Fill( varpos , ww );
          else if ( catP == 8 ) hAllFlav_Tagger_Fak->Fill( varpos , ww );
          else		   hAllFlav_Tagger_Oth->Fill( varpos , ww );
        }
   
        if ( flavour == 1 || flavour == 21 ) { // light q+gluon
        for( int i = 0; i < 3; ++i )
        {
          if ( TagNeg[i] ) {
            hLightFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hLightFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hLightFlav_NegTag_JetEta[i]->Fill( etajet , ww );
            if ( FakN ) {
   	      hLightFlav_Fak_NegTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_Fak_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_Fak_NegTag_JetEta[i]->Fill( etajet , ww );
            }
            if ( OthN ) {
   	      hLightFlav_Oth_NegTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_Oth_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_Oth_NegTag_JetEta[i]->Fill( etajet , ww );
   	      hLightFlav_Oth_NegTag_nPU[i]->Fill( npu , ww );
   	      hLightFlav_Oth_NegTag_pthat[i]->Fill( pthat , ww );
            }
          } else break;
        }
        for( int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hLightFlav_PosTag_JetPU[i]->Fill( npu , ww0 );
            hLightFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hLightFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hLightFlav_PosTag_JetEta[i]->Fill( etajet , ww );
            if ( catP >= 1 && catP <= 3 ) {
   	      hLightFlav_BCT_PosTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_BCT_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_BCT_PosTag_JetEta[i]->Fill( etajet , ww );
            }
            if ( GamP ) {
   	      hLightFlav_Gam_PosTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_Gam_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_Gam_PosTag_JetEta[i]->Fill( etajet , ww );
            }
            if ( K0sP ) {
   	      hLightFlav_K0s_PosTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_K0s_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_K0s_PosTag_JetEta[i]->Fill( etajet , ww );
            }
            if ( LamP ) {
   	      hLightFlav_Lam_PosTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_Lam_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_Lam_PosTag_JetEta[i]->Fill( etajet , ww );
            }
            if ( FakP ) {
   	      hLightFlav_Fak_PosTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_Fak_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_Fak_PosTag_JetEta[i]->Fill( etajet , ww );
   	      hLightFlav_Fak_PosTag_nPU[i]->Fill( npu , ww );
   	      hLightFlav_Fak_PosTag_pthat[i]->Fill( pthat , ww );
            }
            if ( OthP ) {
   	      hLightFlav_Oth_PosTag_JetPV[i]->Fill( npv , ww );
   	      hLightFlav_Oth_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
   	      hLightFlav_Oth_PosTag_JetEta[i]->Fill( etajet , ww );
   	      hLightFlav_Oth_PosTag_nPU[i]->Fill( npu , ww );
   	      hLightFlav_Oth_PosTag_pthat[i]->Fill( pthat , ww );
            }
            if ( njet30 == 1 ) hLightFlav_PosTag_1JetPt[i]->Fill( ptjet , ww );
            if ( njet30 == 2 ) hLightFlav_PosTag_2JetPt[i]->Fill( ptjet , ww );
            if ( njet30 == 3 ) hLightFlav_PosTag_3JetPt[i]->Fill( ptjet , ww );
            if ( njet30 >= 4 ) hLightFlav_PosTag_4JetPt[i]->Fill( ptjet , ww );
          } else break;
        }
        if ( varneg > 0 ) hLightFlav_Tagger->Fill(-varneg , etajet , ww );
        if ( varpos > 0 ) hLightFlav_Tagger->Fill( varpos , etajet , ww );
        if ( varPos > 0 ) hLightFlav_PosTagger->Fill( varPos , ww );
        if ( varPos > 0 && nPU < 10  ) hLightFlav_PosTagger_0pu->Fill( varPos , ww );
        if ( varPos > 0 && nPU >= 20 ) hLightFlav_PosTagger_GE8pu->Fill( varPos , ww );
        if ( varneg > 0 ) {
          if      ( catN == 1 ) hLightFlav_Tagger_Bwd->Fill(-varneg , ww );
          else if ( catN == 2 ) hLightFlav_Tagger_Cwd->Fill(-varneg , ww );
          else if ( catN == 3 ) hLightFlav_Tagger_Tau->Fill(-varneg , ww );
          else if ( catN == 4 ) hLightFlav_Tagger_Gam->Fill(-varneg , ww );
          else if ( catN == 5 ) hLightFlav_Tagger_K0s->Fill(-varneg , ww );
          else if ( catN == 6 ) hLightFlav_Tagger_Lam->Fill(-varneg , ww );
          else if ( catN == 7 ) hLightFlav_Tagger_Int->Fill(-varneg , ww );
          else if ( catN == 8 ) hLightFlav_Tagger_Fak->Fill(-varneg , ww );
          else		     hLightFlav_Tagger_Oth->Fill(-varneg , ww );
        }
        if ( varpos > 0 ) {
          if      ( catP == 1 ) hLightFlav_Tagger_Bwd->Fill( varpos , ww );
          else if ( catP == 2 ) hLightFlav_Tagger_Cwd->Fill( varpos , ww );
          else if ( catP == 3 ) hLightFlav_Tagger_Tau->Fill( varpos , ww );
          else if ( catP == 4 ) hLightFlav_Tagger_Gam->Fill( varpos , ww );
          else if ( catP == 5 ) hLightFlav_Tagger_K0s->Fill( varpos , ww );
          else if ( catP == 6 ) hLightFlav_Tagger_Lam->Fill( varpos , ww );
          else if ( catP == 7 ) hLightFlav_Tagger_Int->Fill( varpos , ww );
          else if ( catP == 8 ) hLightFlav_Tagger_Fak->Fill( varpos , ww );
          else		     hLightFlav_Tagger_Oth->Fill( varpos , ww );
        }
      } // light q+gluon

      if ( flavour == 21 ) { // gluon jets
        if ( varpos > 0 ) hGluonFlav_PosTag->Fill( varpos , ww );
        for( int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hGluonFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hGluonFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hGluonFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 ) hGluonFlav_Tagger->Fill(-varneg , etajet , ww );
        if ( varpos > 0 ) hGluonFlav_Tagger->Fill( varpos , etajet , ww );
      } // gluon jets
      else if ( flavour == 1 ) { // uds jets
        if ( varpos > 0 ) hUDSFlav_PosTag->Fill( varpos , ww );
        for( int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hUDSFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hUDSFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hUDSFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 ) hUDSFlav_Tagger->Fill(-varneg , etajet , ww );
        if ( varpos > 0 ) hUDSFlav_Tagger->Fill( varpos , etajet , ww );
      } // uds jets
      else if ( flavour == 4 ) { // c jets
        if ( varpos > 0 ) hCFlav_PosTag->Fill( varpos , ww );
        if ( varPos > 0 ) hCFlav_PosTagger->Fill( varPos , ww );
        if ( varPos > 0 && nPU < 10  ) hCFlav_PosTagger_0pu->Fill( varPos , ww );
        if ( varPos > 0 && nPU >= 20 ) hCFlav_PosTagger_GE8pu->Fill( varPos , ww );
        for( int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hCFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hCFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hCFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 ) hCFlav_Tagger->Fill(-varneg , etajet , ww );
        if ( varpos > 0 ) hCFlav_Tagger->Fill( varpos , etajet , ww );
      } // c jets
      else if ( flavour == 5 ) { // b jets
        if ( varpos > 0 ) hBFlav_PosTag->Fill( varpos , ww );
        if ( varPos > 0 ) hBFlav_PosTagger->Fill( varPos , ww );
        if ( varPos > 0 && nPU < 10  ) hBFlav_PosTagger_0pu->Fill( varPos , ww );
        if ( varPos > 0 && nPU >= 20 ) hBFlav_PosTagger_GE8pu->Fill( varPos , ww );
        for( int i = 0; i < 3; ++i )
        {
          if ( TagPos[i] ) {
            hBFlav_PosTag_JetPU[i]->Fill( npu , ww0 );
            hBFlav_PosTag_JetPV[i]->Fill( npv , ww );
            hBFlav_PosTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hBFlav_PosTag_JetEta[i]->Fill( etajet , ww );
          } else break;
        }
        if ( varneg > 0 ) hBFlav_Tagger->Fill(-varneg , etajet , ww );
        if ( varpos > 0 ) hBFlav_Tagger->Fill( varpos , etajet , ww );
      } // b jets
      else { // no flavour
        if ( varneg > 0 ) hNoFlav_Tagger->Fill(-varneg , etajet , ww );
        if ( varpos > 0 ) hNoFlav_Tagger->Fill( varpos , etajet , ww );
      } // no flavour

//*********************************
// Negative Tag 
//$$
      for( int i = 0; i < 3; ++i )
      {
        if ( TagNeg[i] ) {
          if ( GamN ) {
            hAllFlav_Gam_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_Gam_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_Gam_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( K0sN ) {
            hAllFlav_K0s_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_K0s_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_K0s_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( LamN ) {
            hAllFlav_Lam_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_Lam_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_Lam_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          if ( FakN ) {
            hAllFlav_Fak_NegTag_JetPV[i]->Fill( npv , ww );
            hAllFlav_Fak_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hAllFlav_Fak_NegTag_JetEta[i]->Fill( etajet , ww );
          }
       
          if (flavour == 21) {
            hGluonFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hGluonFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hGluonFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          else if (flavour == 1) {
            hUDSFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hUDSFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hUDSFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          else if (flavour == 4) {
            hCFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hCFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hCFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
          else if (flavour == 5) {
            hBFlav_NegTag_JetPV[i]->Fill( npv , ww );
            hBFlav_NegTag_JetPt[i]->Fill( ptjet, etajet , ww );
            hBFlav_NegTag_JetEta[i]->Fill( etajet , ww );
          }
        } else break;
      }
    } // end loop on SoftDrop Subjet
  } // end loop on FatJets 
} // end loop on Events 

  std::cout << "events : " << allevents << std::endl ;   

//################################################

// Output Postscript

//   TCanvas* c = new TCanvas("c");
  hData_All_NJets -> Draw(); 
//   c->Print("output.ps(");
  
//################################################
  HistogramManager h ;
  
  h.WriteAllHistogramsInFile(filename.Data(),"recreate");
  h.DeleteAllHistograms();	// In order to avoid memory leaks
//################################################
}   



SubJetsAna2D::SubJetsAna2D(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MCTH/QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt50nsRecodebug_MCRUN2_74_V9A/FatJetsAna2D_Pt-15to30.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MCTH/QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt50nsRecodebug_MCRUN2_74_V9A/FatJetsAna2D_Pt-15to30.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("MCTH/QCD_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt50nsRecodebug_MCRUN2_74_V9A/JetTree_Pt-15to30.root:/btagana");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

SubJetsAna2D::~SubJetsAna2D()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SubJetsAna2D::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SubJetsAna2D::LoadTree(Long64_t entry)
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

void SubJetsAna2D::Init(TTree *tree)
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

   fChain->SetBranchAddress("nBitTrigger", &nBitTrigger, &b_nBitTrigger);
   fChain->SetBranchAddress("BitTrigger", BitTrigger, &b_BitTrigger);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("mcweight", &mcweight, &b_mcweight);
   fChain->SetBranchAddress("BX", &BX, &b_BX);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("PVez", &PVez, &b_PVez);
   fChain->SetBranchAddress("GenPVz", &GenPVz, &b_GenPVz);
   fChain->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("PU_bunch", &PU_bunch, &b_PU_bunch);
   fChain->SetBranchAddress("PU_z", &PU_z, &b_PU_z);
   fChain->SetBranchAddress("PU_sumpT_low", &PU_sumpT_low, &b_PU_sumpT_low);
   fChain->SetBranchAddress("PU_sumpT_high", &PU_sumpT_high, &b_PU_sumpT_high);
   fChain->SetBranchAddress("PU_ntrks_low", &PU_ntrks_low, &b_PU_ntrks_low);
   fChain->SetBranchAddress("PU_ntrks_high", &PU_ntrks_high, &b_PU_ntrks_high);
   fChain->SetBranchAddress("ncQuarks", &ncQuarks, &b_ncQuarks);
   fChain->SetBranchAddress("cQuark_pT", &cQuark_pT, &b_cQuark_pT);
   fChain->SetBranchAddress("cQuark_eta", &cQuark_eta, &b_cQuark_eta);
   fChain->SetBranchAddress("cQuark_phi", &cQuark_phi, &b_cQuark_phi);
   fChain->SetBranchAddress("cQuark_pdgID", &cQuark_pdgID, &b_cQuark_pdgID);
   fChain->SetBranchAddress("cQuark_status", &cQuark_status, &b_cQuark_status);
   fChain->SetBranchAddress("cQuark_fromGSP", &cQuark_fromGSP, &b_cQuark_fromGSP);
   fChain->SetBranchAddress("nbQuarks", &nbQuarks, &b_nbQuarks);
   fChain->SetBranchAddress("bQuark_pT", &bQuark_pT, &b_bQuark_pT);
   fChain->SetBranchAddress("bQuark_eta", &bQuark_eta, &b_bQuark_eta);
   fChain->SetBranchAddress("bQuark_phi", &bQuark_phi, &b_bQuark_phi);
   fChain->SetBranchAddress("bQuark_pdgID", &bQuark_pdgID, &b_bQuark_pdgID);
   fChain->SetBranchAddress("bQuark_status", &bQuark_status, &b_bQuark_status);
   fChain->SetBranchAddress("bQuark_fromGSP", &bQuark_fromGSP, &b_bQuark_fromGSP);
   fChain->SetBranchAddress("nBHadrons", &nBHadrons, &b_nBHadrons);
   fChain->SetBranchAddress("BHadron_pT", BHadron_pT, &b_BHadron_pT);
   fChain->SetBranchAddress("BHadron_eta", BHadron_eta, &b_BHadron_eta);
   fChain->SetBranchAddress("BHadron_phi", BHadron_phi, &b_BHadron_phi);
   fChain->SetBranchAddress("BHadron_mass", BHadron_mass, &b_BHadron_mass);
   fChain->SetBranchAddress("BHadron_pdgID", BHadron_pdgID, &b_BHadron_pdgID);
   fChain->SetBranchAddress("BHadron_mother", BHadron_mother, &b_BHadron_mother);
   fChain->SetBranchAddress("BHadron_hasBdaughter", BHadron_hasBdaughter, &b_BHadron_hasBdaughter);
   fChain->SetBranchAddress("BHadron_SVx", BHadron_SVx, &b_BHadron_SVx);
   fChain->SetBranchAddress("BHadron_SVy", BHadron_SVy, &b_BHadron_SVy);
   fChain->SetBranchAddress("BHadron_SVz", BHadron_SVz, &b_BHadron_SVz);
   fChain->SetBranchAddress("BHadron_nCharged", BHadron_nCharged, &b_BHadron_nCharged);
   fChain->SetBranchAddress("BHadron_DHadron1", BHadron_DHadron1, &b_BHadron_DHadron1);
   fChain->SetBranchAddress("BHadron_DHadron2", BHadron_DHadron2, &b_BHadron_DHadron2);
   fChain->SetBranchAddress("nDHadrons", &nDHadrons, &b_nDHadrons);
   fChain->SetBranchAddress("nDaughters", &nDaughters, &b_nDaughters);
   fChain->SetBranchAddress("DHadron_pT", DHadron_pT, &b_DHadron_pT);
   fChain->SetBranchAddress("DHadron_eta", DHadron_eta, &b_DHadron_eta);
   fChain->SetBranchAddress("DHadron_phi", DHadron_phi, &b_DHadron_phi);
   fChain->SetBranchAddress("DHadron_pdgID", DHadron_pdgID, &b_DHadron_pdgID);
   fChain->SetBranchAddress("DHadron_mass", DHadron_mass, &b_DHadron_mass);
   fChain->SetBranchAddress("DHadron_SVx", DHadron_SVx, &b_DHadron_SVx);
   fChain->SetBranchAddress("DHadron_SVy", DHadron_SVy, &b_DHadron_SVy);
   fChain->SetBranchAddress("DHadron_SVz", DHadron_SVz, &b_DHadron_SVz);
   fChain->SetBranchAddress("DHadron_nDaughters", DHadron_nDaughters, &b_DHadron_nDaughters);
   fChain->SetBranchAddress("DHadron_DaughtersPdgID", DHadron_DaughtersPdgID, &b_DHadron_DaughtersPdgID);
   fChain->SetBranchAddress("DHadron_nChargedDaughters", DHadron_nChargedDaughters, &b_DHadron_nChargedDaughters);
   fChain->SetBranchAddress("DHadron_nCharged", DHadron_nCharged, &b_DHadron_nCharged);
   fChain->SetBranchAddress("nGenlep", &nGenlep, &b_nGenlep);
   fChain->SetBranchAddress("Genlep_pT", Genlep_pT, &b_Genlep_pT);
   fChain->SetBranchAddress("Genlep_eta", Genlep_eta, &b_Genlep_eta);
   fChain->SetBranchAddress("Genlep_phi", Genlep_phi, &b_Genlep_phi);
   fChain->SetBranchAddress("Genlep_pdgID", Genlep_pdgID, &b_Genlep_pdgID);
   fChain->SetBranchAddress("Genlep_status", Genlep_status, &b_Genlep_status);
   fChain->SetBranchAddress("Genlep_mother", Genlep_mother, &b_Genlep_mother);
   fChain->SetBranchAddress("nGenquark", &nGenquark, &b_nGenquark);
   fChain->SetBranchAddress("Genquark_pT", Genquark_pT, &b_Genquark_pT);
   fChain->SetBranchAddress("Genquark_eta", Genquark_eta, &b_Genquark_eta);
   fChain->SetBranchAddress("Genquark_phi", Genquark_phi, &b_Genquark_phi);
   fChain->SetBranchAddress("Genquark_pdgID", Genquark_pdgID, &b_Genquark_pdgID);
   fChain->SetBranchAddress("Genquark_mother", Genquark_mother, &b_Genquark_mother);
   fChain->SetBranchAddress("nGenPruned", &nGenPruned, &b_nGenPruned);
   fChain->SetBranchAddress("GenPruned_pT", GenPruned_pT, &b_GenPruned_pT);
   fChain->SetBranchAddress("GenPruned_eta", GenPruned_eta, &b_GenPruned_eta);
   fChain->SetBranchAddress("GenPruned_phi", GenPruned_phi, &b_GenPruned_phi);
   fChain->SetBranchAddress("GenPruned_mass", GenPruned_mass, &b_GenPruned_mass);
   fChain->SetBranchAddress("GenPruned_pdgID", GenPruned_pdgID, &b_GenPruned_pdgID);
   fChain->SetBranchAddress("GenPruned_status", GenPruned_status, &b_GenPruned_status);
   fChain->SetBranchAddress("GenPruned_mother", GenPruned_mother, &b_GenPruned_mother);
   fChain->SetBranchAddress("nGenV0", &nGenV0, &b_nGenV0);
   fChain->SetBranchAddress("GenV0_pT", &GenV0_pT, &b_GenV0_pT);
   fChain->SetBranchAddress("GenV0_eta", &GenV0_eta, &b_GenV0_eta);
   fChain->SetBranchAddress("GenV0_phi", &GenV0_phi, &b_GenV0_phi);
   fChain->SetBranchAddress("GenV0_pdgID", &GenV0_pdgID, &b_GenV0_pdgID);
   fChain->SetBranchAddress("GenV0_SVx", &GenV0_SVx, &b_GenV0_SVx);
   fChain->SetBranchAddress("GenV0_SVy", &GenV0_SVy, &b_GenV0_SVy);
   fChain->SetBranchAddress("GenV0_SVz", &GenV0_SVz, &b_GenV0_SVz);
   fChain->SetBranchAddress("GenV0_nCharged", &GenV0_nCharged, &b_GenV0_nCharged);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_genpt", Jet_genpt, &b_Jet_genpt);
   fChain->SetBranchAddress("Jet_residual", Jet_residual, &b_Jet_residual);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_jes", Jet_jes, &b_Jet_jes);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_ntracks", Jet_ntracks, &b_Jet_ntracks);
   fChain->SetBranchAddress("Jet_nseltracks", Jet_nseltracks, &b_Jet_nseltracks);
   fChain->SetBranchAddress("Jet_flavour", Jet_flavour, &b_Jet_flavour);
   fChain->SetBranchAddress("Jet_nbHadrons", Jet_nbHadrons, &b_Jet_nbHadrons);
   fChain->SetBranchAddress("Jet_ncHadrons", Jet_ncHadrons, &b_Jet_ncHadrons);
   fChain->SetBranchAddress("Jet_Ip2N", Jet_Ip2N, &b_Jet_Ip2N);
   fChain->SetBranchAddress("Jet_Ip2P", Jet_Ip2P, &b_Jet_Ip2P);
   fChain->SetBranchAddress("Jet_Ip3N", Jet_Ip3N, &b_Jet_Ip3N);
   fChain->SetBranchAddress("Jet_Ip3P", Jet_Ip3P, &b_Jet_Ip3P);
   fChain->SetBranchAddress("Jet_ProbaN", Jet_ProbaN, &b_Jet_ProbaN);
   fChain->SetBranchAddress("Jet_ProbaP", Jet_ProbaP, &b_Jet_ProbaP);
   fChain->SetBranchAddress("Jet_Proba", Jet_Proba, &b_Jet_Proba);
   fChain->SetBranchAddress("Jet_BprobN", Jet_BprobN, &b_Jet_BprobN);
   fChain->SetBranchAddress("Jet_BprobP", Jet_BprobP, &b_Jet_BprobP);
   fChain->SetBranchAddress("Jet_Bprob", Jet_Bprob, &b_Jet_Bprob);
   fChain->SetBranchAddress("Jet_SvxN", Jet_SvxN, &b_Jet_SvxN);
   fChain->SetBranchAddress("Jet_Svx", Jet_Svx, &b_Jet_Svx);
   fChain->SetBranchAddress("Jet_SvxNHP", Jet_SvxNHP, &b_Jet_SvxNHP);
   fChain->SetBranchAddress("Jet_SvxHP", Jet_SvxHP, &b_Jet_SvxHP);
   fChain->SetBranchAddress("Jet_CombSvxN", Jet_CombSvxN, &b_Jet_CombSvxN);
   fChain->SetBranchAddress("Jet_CombSvxP", Jet_CombSvxP, &b_Jet_CombSvxP);
   fChain->SetBranchAddress("Jet_CombSvx", Jet_CombSvx, &b_Jet_CombSvx);
   fChain->SetBranchAddress("Jet_CombIVF", Jet_CombIVF, &b_Jet_CombIVF);
   fChain->SetBranchAddress("Jet_CombIVF_P", Jet_CombIVF_P, &b_Jet_CombIVF_P);
   fChain->SetBranchAddress("Jet_CombIVF_N", Jet_CombIVF_N, &b_Jet_CombIVF_N);
   fChain->SetBranchAddress("Jet_SoftMuN", Jet_SoftMuN, &b_Jet_SoftMuN);
   fChain->SetBranchAddress("Jet_SoftMuP", Jet_SoftMuP, &b_Jet_SoftMuP);
   fChain->SetBranchAddress("Jet_SoftMu", Jet_SoftMu, &b_Jet_SoftMu);
   fChain->SetBranchAddress("Jet_SoftElN", Jet_SoftElN, &b_Jet_SoftElN);
   fChain->SetBranchAddress("Jet_SoftElP", Jet_SoftElP, &b_Jet_SoftElP);
   fChain->SetBranchAddress("Jet_SoftEl", Jet_SoftEl, &b_Jet_SoftEl);
   fChain->SetBranchAddress("Jet_DoubleSV", Jet_DoubleSV, &b_Jet_DoubleSV);
   fChain->SetBranchAddress("Jet_cMVA", Jet_cMVA, &b_Jet_cMVA);
   fChain->SetBranchAddress("Jet_cMVAv2", Jet_cMVAv2, &b_Jet_cMVAv2);
   fChain->SetBranchAddress("Jet_cMVAv2N", Jet_cMVAv2N, &b_Jet_cMVAv2N);
   fChain->SetBranchAddress("Jet_cMVAv2P", Jet_cMVAv2P, &b_Jet_cMVAv2P);
   fChain->SetBranchAddress("Jet_hist1", Jet_hist1, &b_Jet_hist1);
   fChain->SetBranchAddress("Jet_hist2", Jet_hist2, &b_Jet_hist2);
   fChain->SetBranchAddress("Jet_hist3", Jet_hist3, &b_Jet_hist3);
   fChain->SetBranchAddress("Jet_histJet", Jet_histJet, &b_Jet_histJet);
   fChain->SetBranchAddress("Jet_histSvx", Jet_histSvx, &b_Jet_histSvx);
   fChain->SetBranchAddress("Jet_SV_multi", Jet_SV_multi, &b_Jet_SV_multi);
   fChain->SetBranchAddress("Jet_nSM", Jet_nSM, &b_Jet_nSM);
   fChain->SetBranchAddress("Jet_nSE", Jet_nSE, &b_Jet_nSE);
   fChain->SetBranchAddress("Jet_looseID", Jet_looseID, &b_Jet_looseID);
   fChain->SetBranchAddress("Jet_tightID", Jet_tightID, &b_Jet_tightID);
   fChain->SetBranchAddress("Jet_nFirstSE", Jet_nFirstSE, &b_Jet_nFirstSE);
   fChain->SetBranchAddress("Jet_nLastSE", Jet_nLastSE, &b_Jet_nLastSE);
   fChain->SetBranchAddress("nPFElectron", &nPFElectron, &b_nPFElectron);
   fChain->SetBranchAddress("PFElectron_IdxJet", PFElectron_IdxJet, &b_PFElectron_IdxJet);
   fChain->SetBranchAddress("PFElectron_pt", PFElectron_pt, &b_PFElectron_pt);
   fChain->SetBranchAddress("PFElectron_eta", PFElectron_eta, &b_PFElectron_eta);
   fChain->SetBranchAddress("PFElectron_phi", PFElectron_phi, &b_PFElectron_phi);
   fChain->SetBranchAddress("PFElectron_ptrel", PFElectron_ptrel, &b_PFElectron_ptrel);
   fChain->SetBranchAddress("PFElectron_deltaR", PFElectron_deltaR, &b_PFElectron_deltaR);
   fChain->SetBranchAddress("PFElectron_ratio", PFElectron_ratio, &b_PFElectron_ratio);
   fChain->SetBranchAddress("PFElectron_ratioRel", PFElectron_ratioRel, &b_PFElectron_ratioRel);
   fChain->SetBranchAddress("PFElectron_IP", PFElectron_IP, &b_PFElectron_IP);
   fChain->SetBranchAddress("PFElectron_IP2D", PFElectron_IP2D, &b_PFElectron_IP2D);
   fChain->SetBranchAddress("Jet_nFirstSM", Jet_nFirstSM, &b_Jet_nFirstSM);
   fChain->SetBranchAddress("Jet_nLastSM", Jet_nLastSM, &b_Jet_nLastSM);
   fChain->SetBranchAddress("nPFMuon", &nPFMuon, &b_nPFMuon);
   fChain->SetBranchAddress("PFMuon_IdxJet", PFMuon_IdxJet, &b_PFMuon_IdxJet);
   fChain->SetBranchAddress("PFMuon_nMuHit", PFMuon_nMuHit, &b_PFMuon_nMuHit);
   fChain->SetBranchAddress("PFMuon_nTkHit", PFMuon_nTkHit, &b_PFMuon_nTkHit);
   fChain->SetBranchAddress("PFMuon_nPixHit", PFMuon_nPixHit, &b_PFMuon_nPixHit);
   fChain->SetBranchAddress("PFMuon_nOutHit", PFMuon_nOutHit, &b_PFMuon_nOutHit);
   fChain->SetBranchAddress("PFMuon_nTkLwM", PFMuon_nTkLwM, &b_PFMuon_nTkLwM);
   fChain->SetBranchAddress("PFMuon_nPixLwM", PFMuon_nPixLwM, &b_PFMuon_nPixLwM);
   fChain->SetBranchAddress("PFMuon_nMatched", PFMuon_nMatched, &b_PFMuon_nMatched);
   fChain->SetBranchAddress("PFMuon_chi2", PFMuon_chi2, &b_PFMuon_chi2);
   fChain->SetBranchAddress("PFMuon_chi2Tk", PFMuon_chi2Tk, &b_PFMuon_chi2Tk);
   fChain->SetBranchAddress("PFMuon_isGlobal", PFMuon_isGlobal, &b_PFMuon_isGlobal);
   fChain->SetBranchAddress("PFMuon_hist", PFMuon_hist, &b_PFMuon_hist);
   fChain->SetBranchAddress("PFMuon_pt", PFMuon_pt, &b_PFMuon_pt);
   fChain->SetBranchAddress("PFMuon_eta", PFMuon_eta, &b_PFMuon_eta);
   fChain->SetBranchAddress("PFMuon_phi", PFMuon_phi, &b_PFMuon_phi);
   fChain->SetBranchAddress("PFMuon_ptrel", PFMuon_ptrel, &b_PFMuon_ptrel);
   fChain->SetBranchAddress("PFMuon_deltaR", PFMuon_deltaR, &b_PFMuon_deltaR);
   fChain->SetBranchAddress("PFMuon_ratio", PFMuon_ratio, &b_PFMuon_ratio);
   fChain->SetBranchAddress("PFMuon_ratioRel", PFMuon_ratioRel, &b_PFMuon_ratioRel);
   fChain->SetBranchAddress("PFMuon_IP", PFMuon_IP, &b_PFMuon_IP);
   fChain->SetBranchAddress("PFMuon_IP2D", PFMuon_IP2D, &b_PFMuon_IP2D);
   fChain->SetBranchAddress("PFMuon_IPsig", PFMuon_IPsig, &b_PFMuon_IPsig);
   fChain->SetBranchAddress("PFMuon_IP2Dsig", PFMuon_IP2Dsig, &b_PFMuon_IP2Dsig);
   fChain->SetBranchAddress("PFMuon_dz", PFMuon_dz, &b_PFMuon_dz);
   fChain->SetBranchAddress("PFMuon_GoodQuality", PFMuon_GoodQuality, &b_PFMuon_GoodQuality);
   fChain->SetBranchAddress("Jet_nFirstTrkInc", Jet_nFirstTrkInc, &b_Jet_nFirstTrkInc);
   fChain->SetBranchAddress("Jet_nLastTrkInc", Jet_nLastTrkInc, &b_Jet_nLastTrkInc);
   fChain->SetBranchAddress("nTrkInc", &nTrkInc, &b_nTrkInc);
   fChain->SetBranchAddress("TrkInc_pt", TrkInc_pt, &b_TrkInc_pt);
   fChain->SetBranchAddress("TrkInc_eta", TrkInc_eta, &b_TrkInc_eta);
   fChain->SetBranchAddress("TrkInc_phi", TrkInc_phi, &b_TrkInc_phi);
   fChain->SetBranchAddress("TrkInc_ptrel", TrkInc_ptrel, &b_TrkInc_ptrel);
   fChain->SetBranchAddress("TrkInc_IPsig", TrkInc_IPsig, &b_TrkInc_IPsig);
   fChain->SetBranchAddress("TrkInc_IP", TrkInc_IP, &b_TrkInc_IP);
   Notify();
}

Bool_t SubJetsAna2D::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SubJetsAna2D::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SubJetsAna2D::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SubJetsAna2D_cxx
