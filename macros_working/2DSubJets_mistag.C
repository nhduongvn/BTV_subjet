#include <fstream>
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <TROOT.h>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TSystem.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Bases/BTagCalibrationStandalone.cc"
#include "Bases/CMS_lumi.C"
#include "Bases/LowStatAndRebin.C"

/* ****************************************************************************************
// Macro used to produce plots like plot 33... from AN 10-382
// 9 december 2015 version
// Use in root:
  tdrTick(0);
  .L 2DSubJets_mistag.C+g
  TH1F retSF, retMistag;
  TCanvas* tc = plot(retSF, retMistag, 3, false, true, 1, 24, "CSVv2", "L", false, 20., 1000., 20., 1000., "Run2015D", 4., false, 0.0, "2DSubJets")
// ou
  LoopPlot(retSF, retMistag, "CSVv2", "2DSubjets",4.)
  
 ****************************************************************************************/

double mygeterr(bool stat, bool ErrorHisto, TH1F* h, int min, int max)
{
   // If the error to combine are statistic (stat=true),
   //    One must use the error propagation: VarTot = Somme(Var)/(N*(N-1))
   //    Also valid for relative errors.
   // If the error is systematic, I propose to use the mean of the error:
   //    Sigma mean = somme(sigma)/mean.
   // I only consider relative errors here.
   // If these errors are directly stored in the histo I just get them by GetBinContent,
   // If on the contrary, the error comes from the histo errors,
   //    I get the relative errors from GetBinError/GetBinContent

   double dentmp = 0;
   double tmp_output = 0;
   double tmp = 0;
   double myNbin = (double) (max - min);
   if( myNbin > 1 )
   {
     for( int i = min; i < max; ++i )
     {
        tmp = h->GetBinContent(i);
        if( !ErrorHisto )
        {
          if( tmp > 0 ) tmp = h->GetBinError(i)/tmp;
          else tmp = 0.;
        }
        if( stat )
        {
          if( tmp > 0 )
          {
            dentmp += 1/(tmp*tmp);
          }
        } else {
          tmp_output += tmp;
        }
     }
     if( stat )
     {
       if (dentmp > 0 ) tmp_output = sqrt( 1/dentmp );
       else tmp_output = 1.;
     }
     else tmp_output /= myNbin;
   }
   else if ( myNbin == 1 )
   {
     tmp = h->GetBinContent(min);
     if( !ErrorHisto )
     {
       if( tmp > 0 ) tmp = h->GetBinError(min)/tmp;
       else tmp = 0.;
     }
     tmp_output = tmp;
   }
   else
   {
     std::cout << "Warning, you requested a mean with less than one bin " << endl;
     tmp_output = 0;
   }
   return tmp_output;
}

double mygetmean(TH1F* h, int min, int max )
{
   double tmp_output = 0;
   double myNbin = (double) (max - min);
   if( myNbin > 0 )
   {
     for( int i = min; i < max; ++i )
     {
        tmp_output += h->GetBinContent(i);
     }
     tmp_output /= myNbin;
   }
   else
   {
     std::cout << "Warning, you requested a mean with maxbin <= minbin " << endl;
     tmp_output = 0;
   }
   return tmp_output;
}

TCanvas* plot(TH1F &retSF, TH1F &retMistag, int nPlot = 3, Bool_t InLoop = false, Bool_t ImprovePlot = false, int MinEta = 0, int MaxEta = 24, TString CutName="TCHE", TString CutStrength="L", bool interactive = false, float ptmin = 20., float ptmax = 1000., float ptminDraw = 20., float ptmaxDraw = 1000., TString Period = "2015C", Double_t xSigma = 4., bool ManualErrors = false, Double_t ManualRelError = 0.0, TString TypeOfJets = "2DSubjets")
{
  cout << "\n Interactive: " << interactive ;
  //TEMP test reco debug
  bool Manual_v0_Error = true ;
  if( ManualErrors )
  {
    cout << "Warning: using manual relative error = " << ManualRelError << endl;
    cout << "If you agree, type \"y\". " << endl;
    char a_char = getchar();
    if( a_char != 'y' ) return NULL;
  }

  int nbins = 33;
  Double_t xbins[33] = {  20.,  30.,  40.,  50.,  60.,
                          80., 100., 120., 140., 160.,
                         180., 200., 220., 240., 260.,
                         280., 320., 360., 400.,
                         440., 480., 520., 560., 600.,
                         650., 700., 750., 800.,
                         850., 900., 950., 990., 1000.};
  setBinVector(nbins, xbins);		// Define rebinning with previous vector.
  vector<Int_t> LowStatBins;		// Used to contain the low statistic bins to set to zero.

  TString EtaRange = "";
  char buffer[50];

  if( MinEta > MaxEta )
  {
    std::cout << "MinEta should be smaller than MaxEta " << endl;
    return NULL;
  }
  TH1F *htemp = new TH1F("htemp","htemp",24,0,2.4); // Temp hist used to reproduce used binning...
  if( MinEta < 1 || MinEta > htemp->GetNbinsX() )
  {
    std::cout << "MinEta out of range. Should be between 1 and " << htemp->GetNbinsX() << " But found to be " << MinEta << endl;
    return NULL;
  }
  if( MaxEta < 1 || MaxEta > htemp->GetNbinsX() )
  {
    std::cout << "MaxEta out of range. Should be between 1 and " << htemp->GetNbinsX() << " But found to be " << MaxEta << endl;
    return NULL;
  }

  float EtaMin = htemp->GetBinLowEdge(MinEta) ;
  float EtaMax = htemp->GetBinLowEdge(MaxEta) + htemp->GetBinWidth(MaxEta) ;

  sprintf(buffer,"%3.1f",EtaMin);
  EtaRange += buffer;
  EtaRange += "_";
  sprintf(buffer,"%3.1f",EtaMax);
  EtaRange += buffer;

  std::cout << "Will Use EtaRange = " << EtaRange << endl;

  sprintf(buffer,"%3.1f",xSigma);
  TString txSigma = "";
  txSigma += buffer;
  std::cout << "When merging pthat file, only use bins with content > " << txSigma << "*error ! " << endl;

  sprintf(buffer,"%4.2f",ManualRelError);
  TString tManualRelError = "";
  tManualRelError += buffer;

  TString oFilePrefix = "xSig_" + txSigma + "_";
  if( ManualErrors ) oFilePrefix += "MErr_" + tManualRelError + "_";
  else oFilePrefix += "AErr_";

  vector<TString> PthatRanges;
  //PthatRanges.push_back("/Pt-15to30/");
  //PthatRanges.push_back("/Pt-30to50/");
  //PthatRanges.push_back("/Pt-50to80/");
  //PthatRanges.push_back("/Pt-80to120/");
  //PthatRanges.push_back("/Pt-120to170/");
  //PthatRanges.push_back("/Pt-170to300/"); //Maybe adding 170 to 300 bins?
  PthatRanges.push_back("/Pt-300to470/");
  PthatRanges.push_back("/Pt-470to600_all/");
  PthatRanges.push_back("/Pt-600to800_all/");
  PthatRanges.push_back("/Pt-800to1000/");
  PthatRanges.push_back("/Pt-1000to1400/");
  PthatRanges.push_back("/Pt-1400to1800/");
  PthatRanges.push_back("/Pt-1800to2400/");
  PthatRanges.push_back("/Pt-2400to3200/");
  //PthatRanges.push_back("/Pt-3200toInf/"); 


  TCanvas *c1 = NULL;
  TString tnPlot = "";
  if( nPlot > 3 || nPlot < 1 ) { std::cout << "Wrong number of plots" << endl; return c1; }
  tnPlot += nPlot;
  TString CutandStrength = " ";
  std::cout << CutandStrength << " " << CutName << " " << CutStrength << endl;
  CutandStrength = CutName + CutStrength;
  if (CutName == "CvsL" || CutName == "CvsB") CutandStrength = "ctag" + CutStrength ;
  std::cout << CutandStrength << " " << CutName << " " << CutStrength << endl;
//--------- Boolean section -----------//
//bool interactive = false;
bool CheckFitMistag = false;
bool CheckFitSF = false;
bool PrintDBvalues = true;
bool PrintDBfunction = true;
bool PrintDBSF_mistag = true;
bool PrintDBSF = true;
bool oldData = false;
bool CorrectForLambdaV0 = true;
//bool CorrectForLambdaV0 = false;
TString oFil = "";
TString SavDir = "";
//--------- End of boolean section ----//

  //TString xtitle = "p_{T} [GeV/c]";
  TString xtitle = "Jet p_{T} [GeV]";

// ****************************************************************************************
// Initialisation of the variables according to the parameters with which "plot" is called
// ****************************************************************************************
  TString DirOfDataFile, DirOfMCFile, DirOfMCFileTH;
  TString RootOfDataFile, RootOfMCFile, JetCut, Lumi, textTitle = "";
  TString JetCuts[10];
  Bool_t includeCorrupted = true;
  int minCorrupted_i[10], maxCorrupted_i[10];
  int NJetCuts = 1 ; //=1 for subjet =8 for normal jet
    //JetCuts[0] = "40"; minCorrupted_i[0] = 0; maxCorrupted_i[0] = ((30+50)-50)/10;
    //JetCuts[1] = "60"; minCorrupted_i[1] = 0; maxCorrupted_i[1] = ((60+50)-50)/10;
    //JetCuts[2] = "80"; minCorrupted_i[2] = (80-50)/10+1; maxCorrupted_i[2] = ((80+50)-50)/10;
    //JetCuts[3] = "140"; minCorrupted_i[3] = (140-50)/10+1; maxCorrupted_i[3] = ((140+50)-50)/10;
    //JetCuts[4] = "200"; minCorrupted_i[4] = (200-50)/10+1; maxCorrupted_i[4] = ((200+50)-50)/10;
    //JetCuts[5] = "260"; minCorrupted_i[5] = (260-50)/10+1; maxCorrupted_i[5] = ((260+80)-50)/10;
    //JetCuts[6] = "320"; minCorrupted_i[6] = (320-50)/10+1; maxCorrupted_i[6] = ((320+80)-50)/10;
    //JetCuts[7] = "400"; minCorrupted_i[7] = (400-50)/10+1; maxCorrupted_i[7] = ((400+80)-50)/10;
    JetCuts[0] = "400"; minCorrupted_i[0] = (400-50)/10+1; maxCorrupted_i[0] = ((400+80)-50)/10; //used for subjet together with changing to =1 at NJetCuts= in line 238
    //JetCuts[0] = "500"; minCorrupted_i[0] = (500-50)/10+1; maxCorrupted_i[0] = ((500+80)-50)/10; //used for subjet together with changing to =1 at NJetCuts= in line 238


// 2015C
if ( Period == "Run2016" )
{
    //DirOfDataFile = "../../outputs/CMSSW_7_6_3/Data/JetHT/Run2015D-16Dec2015-v1/"+TypeOfJets+"/TrigType_2015/";
    //DirOfMCFile = "../../outputs/CMSSW_7_6_3/MC/QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/"+TypeOfJets+"/TrigType_2015/";
    //DirOfMCFileTH = "../../outputs/CMSSW_7_6_3/MC/QCD_TuneCUETP8M1_13TeV_pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12/"+TypeOfJets+"/TrigType_2015/";
    //TEMP
    DirOfDataFile = "../Output_jetMass50_data_7p7/Data/JetHT/JetHT_Run2016BC_7p7invfb/" ;
    DirOfMCFile = "../Output_jetMass50_data_7p7/MC/QCD/" ;
    DirOfMCFileTH = "../Output_jetMass50_data_7p7/MC/QCD/" ;

    SavDir = "Dir_"+TypeOfJets+"_jetMass50_data_7p7/";

    RootOfMCFile   = "JetTree_mc_FatJets_Subjets_sampleList_all_NW_";
    RootOfDataFile = "JetTree_data_FatJets_Subjets_sampleList_all_NW_";

    Lumi = "7.7 fb^{-1} (13 TeV, 25ns)";
}else return NULL;

  gSystem->mkdir(SavDir,kTRUE);
  TString subDir[7] = {"/C", "/hist", "/pdf", "/png", "/root", "/txt", "txtDB"};
  for( int i = 0; i < 7; ++i )
  {
    gSystem->mkdir(SavDir+subDir[i],kTRUE);
  }
  oFil = SavDir + "/txtDB/"+oFilePrefix+CutandStrength+"_"+EtaRange;

// *****************************************************************************

  TString htitle = CutName + CutStrength; // + " tagger " + TypeOfJets; 
  if( CutName == "CSVR" ) htitle = "RCSV" + CutStrength + " tagger";
  if( CutName == "CvsL" || CutName == "CvsB") htitle = "ctag" + CutStrength ;
  if( EtaRange != "0.0_2.4" )
  {
    TString EtaInterval = Form(" ; %1.1f #leq #eta < %1.1f ", EtaMin, EtaMax );
    htitle +=  EtaInterval;
  }
  
  TString bloblo_file = "";
  TString blabla_file = "";
  TString hisFileName = "";

  bloblo_file = SavDir + "/txt/"+oFilePrefix+CutandStrength+ "_" + EtaRange +"_Mistag_stat_scale.txt";  // Name of text output File
  blabla_file = SavDir + "/txt/"+oFilePrefix+CutandStrength+ "_" + EtaRange +"_Mistag.txt";  // Name of text output File
  hisFileName = SavDir + "/hist/"+oFilePrefix+CutandStrength+"_" + EtaRange+".root"; // Name of his output File

 ofstream aout_txt(bloblo_file, ios::out);

// ****************************************************************************************
// End of parameterised initialisation phase
// ****************************************************************************************

 float y;
 float ey;


 string TaggerName = "";
 if ( CutandStrength == "TCHPL" || CutandStrength == "TCHEL" || CutandStrength == "TCHEM" )
    TaggerName = "trackCountingHighEffBJetTags";
 if ( CutandStrength == "TCHPM" || CutandStrength == "TCHPT" || CutandStrength == "TCHET" )
    TaggerName = "trackCountingHighPurBJetTags";
 if ( CutandStrength == "SSVHEM" || CutandStrength == "SSVHPM" )
    TaggerName = "simpleSecondaryVertexHighEffBJetTags";
 if ( CutandStrength == "SSVHET" || CutandStrength == "SSVHPT" )
    TaggerName = "simpleSecondaryVertexHighPurBJetTags";
 if ( CutName == "JP" )
    TaggerName = "jetProbabilityBJetTags";
 if ( CutName == "JBP" )
    TaggerName = "jetBProbabilityBJetTags";
 if ( CutName == "CSV" )
    TaggerName = "combinedSecondaryVertexBJetTags";
 if ( CutName == "CSVv2" )
    TaggerName = CutName;
 if ( CutName == "cMVAv2" )
    TaggerName = CutName;
 if( TaggerName == "" )
 {
   std::cout << "Warning : use Empty TaggerName for some prints... No harm. " << endl;
   if( !InLoop ) getchar();
 }

// Set the cut value for the apropriate working point.
 float wp = 0;
/*
 if ( CutandStrength == "TCHEL" || CutandStrength == "TCHPL" ) wp = 1.70; // TCHE_L
 if ( CutandStrength == "TCHET" || CutandStrength == "TCHPT" ) wp = 3.41; // TCHP_T
 if ( CutandStrength == "TCHEM" ) wp = 3.30; // TCHE_M
 if ( CutandStrength == "TCHPM" ) wp = 1.93; // TCHP_M
 if ( CutandStrength == "JPL" ) wp = .275; // JP_L
 if ( CutandStrength == "JPM" ) wp = .545; // JP_M
 if ( CutandStrength == "JPT" ) wp = .790; // JP_T
 if ( CutandStrength == "JBPL" ) wp = 1.33; // JBP_L
 if ( CutandStrength == "JBPM" ) wp = 2.55; // JBP_M
 if ( CutandStrength == "JBPT" ) wp = 3.74; // JBP_T
 if ( CutandStrength == "CSVL" ) wp = 0.244; // CSV_L
 if ( CutandStrength == "CSVM" ) wp = 0.679; // CSV_M
 if ( CutandStrength == "CSVT" ) wp = 0.898; // CSV_T
 if ( CutandStrength == "CSVv2L" ) wp = 0.605;
 if ( CutandStrength == "CSVv2M" ) wp = 0.890;
 if ( CutandStrength == "CSVv2T" ) wp = 0.970;
 if ( CutandStrength == "SSVHEM" || CutandStrength == "SSVHPM" ) wp = 1.74; // SSVHE_M
 if ( CutandStrength == "SSVHET" || CutandStrength == "SSVHPT" ) wp = 2.00; // SSVHP_T
 if( wp == 0 )
 {
   std::cout << "Warning : use wp = " << wp << " for some prints... No harm. " << endl;
   getchar();
 }
*/

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Initialisation for Systematics

// Flavour fraction
 float Xbb = 0.2;
 float Xcc = 0.2;
 float Xgg = 0.2;

// IP Sign Flip 2011 muon jets 1.1 fb-1
 float Xip = 0;
 if ( CutandStrength == "TCHPL" || CutandStrength == "TCHEL" )
    Xip = 0.041; // loose tc
 if ( CutandStrength == "JPL" )
    Xip = 0.027; // loose jp
 if ( CutandStrength == "JBPL" )
    Xip = 0.023; // loose jbp
 if ( CutandStrength == "CSVL" )
    Xip = 0.045; // loose csv
 if ( CutandStrength == "TCHEM" )
    Xip = 0.014;  // medium tche
 if ( CutandStrength == "TCHPM" )
    Xip = 0.012;  // medium tchp
 if ( CutandStrength == "SSVHEM" || CutandStrength == "SSVHPM" )
    Xip = 0.0067; // medium ssv
 if ( CutandStrength == "JPM" )
    Xip = 0.013;  // medium jp
 if ( CutandStrength == "JBPM" )
    Xip = 0.0040; // medium jbp
 if ( CutandStrength == "CSVM" )
    Xip = 0.016;  // medium csv
 if ( CutandStrength == "TCHPT" || CutandStrength == "TCHET" )
    Xip = 0.0029;  // tight tc
 if ( CutandStrength == "SSVHET" || CutandStrength == "SSVHPT" )
    Xip = 0.00068; // tight ssv
 if ( CutandStrength == "JPT" )
    Xip = 0.0015;  // tight jp
 if ( CutandStrength == "JBPT" )
    Xip = 0.0028;  // tight jbp
 if ( CutandStrength == "CSVT" )
    Xip = 0.0043;  // tight csv

 if( Xip == 0 )
 {
   std::cout << "Warning : use Xip = " << Xip << " for sign flip study " << endl;
   if( !InLoop ) getchar();
 }

// V0 fraction
float SFk0 = 1.3; // Change from 1.4; // to 1.3 As asked by D.Bloch on 8/6/2012.
float SFla = 1.5;
if( !CorrectForLambdaV0 )
{
  SFk0 = 1.;
  SFla = 1.;
}
 float Xk0 = 0.3; // Increase frm 0.15 to 0.4 (0.3)asked by ARC (D.B). Change done on 7(8)/6/2012...//0.15;
 float Xla = 0.50;

// Conversions + secondary interactions
 float Xga = 0.05;

// Fake tracks uncertainty
 float Xfk = 0.5;

// Pile-Up sensitivity in MC
 //float PUsignif = 1.0;

// Event sample systematics
 float ErrSysMis = 0.20;
 float ErrSysSF = 0.10;

 float SF_Trig_syst=0, SF_Run_syst=0;
 float Mistag_Trig_syst=0, Mistag_Run_syst=0;

// SF_Trig_Syst are taken from SF_SystTrigger_Run2015D/Trig*txt
 float syst_lumi = 0.03 ; //add 3% for systematic for change from 7.7 to 9.2 fb
 if( CutandStrength == "CSVv2L" )
 {
   SF_Trig_syst = 0.015; 
   SF_Run_syst = syst_lumi; // change from 0 3%  old comment: From SF_RunSyst/CSVv2L.txt
 }
 if( CutandStrength == "CSVv2M" )
 {
   SF_Trig_syst = 0.032; 
   SF_Run_syst = TMath::Sqrt(0.01409*0.01409 + syst_lumi*syst_lumi) ; // add 3% old comment: From SF_RunSyst/CSVv2M.txt
 }
 if( CutandStrength == "CSVv2T" )
 {
   SF_Trig_syst = 0.12; 
   SF_Run_syst = TMath::Sqrt(0.00068*0.00068 + syst_lumi*syst_lumi) ; // add 3%From SF_RunSyst/CSVv2T.txt
 }
 if( CutandStrength == "cMVAv2L" )
 {
   SF_Trig_syst = 0.010; 
   SF_Run_syst = syst_lumi ; // change to 3% from 0 From SF_RunSyst/CSVv2L.txt
 }
 if( CutandStrength == "cMVAv2M" )
 {
   SF_Trig_syst = 0.041; 
   SF_Run_syst = syst_lumi; // change to 3% from 0 From SF_RunSyst/CSVv2L.txt
 }
 if( CutandStrength == "cMVAv2T" )
 {
   SF_Trig_syst = 0.097; 
   SF_Run_syst = syst_lumi ; // change to 3% from 0 From SF_RunSyst/CSVv2L.txt
 }
 if( CutandStrength == "JPL" )
 {
   SF_Trig_syst = 0.013; 
   SF_Run_syst = 0.00109; // From SF_RunSyst/JPL.txt
 }
 if( CutandStrength == "JPM" )
 {
   SF_Trig_syst = 0.044; 
   SF_Run_syst = 0.00686; // From SF_RunSyst/JPM.txt
 }
 if( CutandStrength == "JPT" )
 {
   SF_Trig_syst = 0.12; 
   SF_Run_syst = 0.00039; // From SF_RunSyst/JPT.txt
 }
 if( CutandStrength == "JBPL" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0444;
 }
 if( CutandStrength == "JBPM" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0562;
 }
 if( CutandStrength == "JBPT" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0523;
 }
 if( CutandStrength == "TCHEL" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0317;
 }
 if( CutandStrength == "TCHEM" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0377;
 }
 if( CutandStrength == "TCHPM" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0348423; //
 }
 if( CutandStrength == "TCHPT" )
 {
   SF_Trig_syst = 0.171207; //From Syst/plot_SF_Syst/Trig_syst.txt
   SF_Run_syst = 0.0348423; // From Syst/plot_SF_Syst/Run_syst.txt

   Mistag_Trig_syst = 0.17382; //From Syst/plot_Mistag_Syst/Trig_syst.txt
   Mistag_Run_syst = 0.0350116; // From Syst/plot_Mistag_Syst/Run_syst.txt
 }
 if( CutandStrength == "SSVHEM" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0201;
 }
 if( CutandStrength == "SSVHPT" )
 {
   SF_Trig_syst = 0.07;
   SF_Run_syst = 0.0355;
 }

 if( SF_Trig_syst == 0 )
 {
   std::cout << "Warning use SF_Trig_syst = " << SF_Trig_syst << " for syst studies " << endl;
   if( !InLoop ) getchar();
 }
 if( SF_Run_syst == 0 )
 {
   std::cout << "Warning use SF_Run_syst = " << SF_Run_syst << " for syst studies " << endl;
   if( !InLoop ) getchar();
 }

 ErrSysSF = sqrt(SF_Trig_syst*SF_Trig_syst + SF_Run_syst*SF_Run_syst);
 ErrSysMis = sqrt(Mistag_Trig_syst*Mistag_Trig_syst + Mistag_Run_syst*Mistag_Run_syst);

 std::cout << CutandStrength << " SF_Trig_syst = " << SF_Trig_syst << " ErrSysSF = " << ErrSysSF << endl;
 std::cout << CutandStrength << " Mistag_Trig_syst = " << Mistag_Trig_syst << " ErrSysMis = " << ErrSysMis << endl;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 Int_t stati=0;
 Bool_t  fit=0;

// Initialise Histograms
//------------------------------------------------------------------------------
// MC with average pile-up and without trigger
 TH1F *h000 = NULL;
 TH1F *h010 = NULL;
 TH1F *h009 = NULL;
 TH1F *h029 = NULL;
 TFile *f2 = NULL;
 for( unsigned long iPthat = 0; iPthat < PthatRanges.size(); ++iPthat )
 {
   TString nameTmp = PthatRanges.at(iPthat) ;
   //nameTmp.ReplaceAll("/","") ;
   //nameTmp.ReplaceAll("-", "_") ;
   //nameTmp = "QCD_" + nameTmp + "_TuneCUETP8M1_13TeV_pythia8/" ;
   f2 = new TFile(DirOfMCFile+nameTmp+RootOfMCFile +"0_"+CutName+".root");
   //TEMP test using only 400 sample
   if( f2 == NULL )
   {
     std::cout << "Error opening " << DirOfMCFile<<PthatRanges.at(iPthat)<<RootOfMCFile <<"0_"<<CutName<<".root" << endl;
     getchar();
     return NULL;
   }// else std::cout << DirOfMCFile<<PthatRanges.at(iPthat)<<RootOfMCFile << "0_" << CutName << ".root opened " << endl;
   f2->cd();


   if( iPthat == 0 )
   {
     GetLowStatBins(f2,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);

     h000 = FindRebinSetLowStatToZero(f2,"hData_JetPt",EtaRange+"_MCavgPUnoTrigger",MinEta,MaxEta,LowStatBins);
     h010 = FindRebinSetLowStatToZero(f2,"hData_NegTag"+CutStrength+"_JetPt", EtaRange+"_MCavgPUnoTrigger",MinEta,MaxEta, LowStatBins );
     h009 = FindRebinSetLowStatToZero(f2,"hLightFlav_JetPt",EtaRange+"_MCavgPUnoTrigger",MinEta,MaxEta,LowStatBins);
     h029 = FindRebinSetLowStatToZero(f2,"hLightFlav_PosTag"+CutStrength+"_JetPt",EtaRange+"_MCavgPUnoTrigger",MinEta,MaxEta,LowStatBins);
   }
   else
   {
     GetLowStatBins(f2,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);
     h000->Add(FindRebinSetLowStatToZero(f2,"hData_JetPt",EtaRange+"_MCavgPUnoTrigger"+PthatRanges.at(iPthat),MinEta,MaxEta,LowStatBins));
     h010->Add(FindRebinSetLowStatToZero(f2,"hData_NegTag"+CutStrength+"_JetPt", EtaRange+"_MCavgPUnoTrigger"+PthatRanges.at(iPthat),MinEta,MaxEta, LowStatBins ));
     h009->Add(FindRebinSetLowStatToZero(f2,"hLightFlav_JetPt", EtaRange+"_MCavgPUnoTrigger"+PthatRanges.at(iPthat),MinEta,MaxEta, LowStatBins ));
     h029->Add(FindRebinSetLowStatToZero(f2,"hLightFlav_PosTag"+CutStrength+"_JetPt", EtaRange+"_MCavgPUnoTrigger"+PthatRanges.at(iPthat),MinEta,MaxEta, LowStatBins ));
   }
   h000->SetDirectory(0);
   h010->SetDirectory(0);
   h009->SetDirectory(0);
   h029->SetDirectory(0);

   f2->Close();
   delete f2;
 }

 if( h000 == NULL ){ std::cout << "h000 NULL " << endl; getchar(); return NULL;}
 if( h010 == NULL ){ std::cout << "h010 NULL " << endl; getchar(); return NULL;}
 if( h009 == NULL ){ std::cout << "h009 NULL " << endl; getchar(); return NULL;}
 if( h029 == NULL ){ std::cout << "h029 NULL " << endl; getchar(); return NULL;}

 //TEMP do not adding data without trigger
 //h000->Reset() ;
 //h010->Reset() ;
 //h009->Reset() ;
 //h029->Reset() ;



 h000->SetTitle("");

  float minPtIntegration = 80;
  float maxPtIntegration = 120;
  float HistPtMin  = h000->GetXaxis()->GetXmin();
  float HistPtMax  = h000->GetXaxis()->GetXmax();
  int   HistPtNbin = h000->GetNbinsX();

 float BinMin[200], BinMax[200]; 
 float Leff[200], LeffMin[200], LeffMax[200], LeffErr[200];
 float Lsf[200], LsfMin[200], LsfMax[200], LsfErr[98+2];
  int iminPtInt = ( minPtIntegration - HistPtMin )/((HistPtMax-HistPtMin)/HistPtNbin) + 1;
  int imaxPtInt = ( maxPtIntegration - HistPtMin )/((HistPtMax-HistPtMin)/HistPtNbin) + 1;

  for( int itmp = 0; itmp <= HistPtNbin+1; ++itmp )
  {
    if( h000->GetBinLowEdge(itmp) >= minPtIntegration )
    {
      iminPtInt = itmp;
      break;
    }
  }
  for( int itmp = 0; itmp <= HistPtNbin+1; ++itmp )
  {
    if( h000->GetBinLowEdge(itmp) >= maxPtIntegration )
    {
      imaxPtInt = itmp;
      break;
    }
  }

  float minPtIntegration2 = 160;
  float maxPtIntegration2 = 320;
  int iminPtInt2 = ( minPtIntegration2 - HistPtMin )/((HistPtMax-HistPtMin)/HistPtNbin) + 1;
  int imaxPtInt2 = ( maxPtIntegration2 - HistPtMin )/((HistPtMax-HistPtMin)/HistPtNbin) + 1;

  for( int itmp = 0; itmp <= HistPtNbin+1; ++itmp )
  {
    if( h000->GetBinLowEdge(itmp) >= minPtIntegration2 )
    {
      iminPtInt2 = itmp;
      break;
    }
  }
  for( int itmp = 0; itmp <= HistPtNbin+1; ++itmp )
  {
    if( h000->GetBinLowEdge(itmp) >= maxPtIntegration2 )
    {
      imaxPtInt2 = itmp;
      break;
    }
  }
  std::cout << "Found Bin extreme edges for 1: " << iminPtInt << "\t" << imaxPtInt << endl;
  std::cout << h000->GetBinLowEdge(iminPtInt) << "\t" << h000->GetBinLowEdge(imaxPtInt) << endl;
  std::cout << "Found Bin extreme edges for 2: " << iminPtInt2 << "\t" << imaxPtInt2 << endl;
  std::cout << h000->GetBinLowEdge(iminPtInt2) << "\t" << h000->GetBinLowEdge(imaxPtInt2) << endl;


  float VminPt[] = { 30,  80, 140, 240, 420, 600};
  float VmaxPt[] = { 80, 140, 240, 420, 600, 1000};
  int   ViminPt[]  = {  0,   0,   0,   0, 0, 0};
  int   VimaxPt[]  = {  0,   0,   0,   0, 0, 0};
  std::cout << "Getting index of histograms corresponding to VminPt and VmaxPt" << endl;
  for( int iV = 0; iV < 6; ++iV )
  {
    for( int itmp = 0; itmp <= HistPtNbin+1; ++itmp )
    {
      if( h000->GetBinLowEdge(itmp) >= VminPt[iV] )
      {
        ViminPt[iV] = itmp;
        break;
      }
    }
    for( int itmp = 0; itmp <= HistPtNbin+1; ++itmp )
    {
      if( h000->GetBinLowEdge(itmp) >= VmaxPt[iV] )
      {
        VimaxPt[iV] = itmp;
        break;
      }
    }
    std::cout << "For VminPt = " << VminPt[iV] << "\t found index = " << ViminPt[iV] << " For which bin low edge is: " << h000->GetBinLowEdge(ViminPt[iV]) << endl;
    std::cout << "For VmaxPt = " << VmaxPt[iV] << "\t found index = " << VimaxPt[iV] << " For which bin low edge is: " << h000->GetBinLowEdge(VimaxPt[iV]) << endl;
    std::cout << endl;
  } 
 

 TH1F* k4= (TH1F*) h000->Clone("k4");
 k4->Reset();
       k4->Divide(h029,h009,1,1,"B"); 

 TH1F* k5= (TH1F*) h000->Clone("k5");
 k5->Reset();
       k5->Divide(h010,h000,1,1,"B");

 TH1F* kRlight= (TH1F*) h000->Clone("kRlight");
 kRlight->Reset();
       kRlight->Divide(k4,k5,1,1);

// stat MC
 TH1F* StatMis= (TH1F*) h000->Clone("StatMis");
 StatMis->Reset();
 TH1F* StatSF= (TH1F*) h000->Clone("StatSF");
 StatSF->Reset();
       for (int i=0; i<HistPtNbin+1; i++) {
	ey = 0.001;
        y = kRlight->GetBinContent(i);
	if (y > 0.) ey = kRlight->GetBinError(i) / y;
	StatMis->SetBinContent(i,ey);
	ey = 0.001;
        y = k5->GetBinContent(i);
	if (y > 0.) ey = k5->GetBinError(i) / y;
	StatSF->SetBinContent(i,ey);
       }
//------------------------------------------------------------------------------
// Loop on MC with average pile-up and with trigger


 TH1F *puMis_Jets[8], *puSF_Jets[8], *bc_Jets[8], *gg_Jets[8];
 TH1F *sf_Jets[8], *fk_Jets[8], *v0_Jets[8], *D00_Jets[8], *D10_Jets[8];

 TH1F *MCNegtag_Jets[8], *cMCNegtag_Jets[8];
 TH1F *MCMistag_Jets[8], *cMCMistag_Jets[8];
 TH1F *MCRlightWeight_Jets[8], *cMCRlightWeight_Jets[8];

 TH1F *SF_Jets[8], *cSF_Jets[8];
 TH1F *DNegtag_Jets[8], *cDNegtag_Jets[8];
 TH1F *DMistag_Jets[8], *cDMistag_Jets[8];

 TH1F* cMCMistag= (TH1F*) h000->Clone("cMCMistag"); cMCMistag->Reset();
 TH1F* cMCNegtag= (TH1F*) h000->Clone("cMCNegtag"); cMCNegtag->Reset();
 TH1F* cDMistag= (TH1F*) h000->Clone("cDMistag"); cDMistag->Reset();
 TH1F* cDNegtag= (TH1F*) h000->Clone("cDNegtag"); cDNegtag->Reset();
 TH1F* cSF= (TH1F*) h000->Clone("cSF"); cSF->Reset();

 TH1F* H00= (TH1F*) h000->Clone("H00"); H00->Reset();
 TH1F* H10= (TH1F*) h000->Clone("H10"); H10->Reset();
 TH1F* H29= (TH1F*) h000->Clone("H29"); H29->Reset();
 TH1F* H09= (TH1F*) h000->Clone("H09"); H09->Reset();
 
 for( int iJetCut = 0; iJetCut < NJetCuts; ++iJetCut)
 {
   puMis_Jets[iJetCut]= (TH1F*) h000->Clone("puMis_"+ JetCuts[iJetCut]); puMis_Jets[iJetCut]->Reset();
   puSF_Jets[iJetCut]= (TH1F*) h000->Clone("puSF_"+JetCuts[iJetCut]); puSF_Jets[iJetCut]->Reset();
   bc_Jets[iJetCut]= (TH1F*) h000->Clone("bc_" + JetCuts[iJetCut]); bc_Jets[iJetCut]->Reset();
   gg_Jets[iJetCut]= (TH1F*) h000->Clone("gg_" + JetCuts[iJetCut]); gg_Jets[iJetCut]->Reset();
   sf_Jets[iJetCut]= (TH1F*) h000->Clone("sf_" + JetCuts[iJetCut]); sf_Jets[iJetCut]->Reset();
   fk_Jets[iJetCut]= (TH1F*) h000->Clone("fk_" + JetCuts[iJetCut]); fk_Jets[iJetCut]->Reset();
   v0_Jets[iJetCut]= (TH1F*) h000->Clone("v0_" + JetCuts[iJetCut]); v0_Jets[iJetCut]->Reset();

   MCNegtag_Jets[iJetCut]= (TH1F*) h000->Clone("MCNegtag_" + JetCuts[iJetCut]); MCNegtag_Jets[iJetCut]->Reset();
   cMCNegtag_Jets[iJetCut]= (TH1F*) h000->Clone("cMCNegtag_" + JetCuts[iJetCut]); cMCNegtag_Jets[iJetCut]->Reset();
   MCMistag_Jets[iJetCut]= (TH1F*) h000->Clone("MCMistag_" + JetCuts[iJetCut]); MCMistag_Jets[iJetCut]->Reset();
   cMCMistag_Jets[iJetCut]= (TH1F*) h000->Clone("cMCMistag_" + JetCuts[iJetCut]); cMCMistag_Jets[iJetCut]->Reset();
   MCRlightWeight_Jets[iJetCut]= (TH1F*) h000->Clone("MCRlightWeight_" + JetCuts[iJetCut]); MCRlightWeight_Jets[iJetCut]->Reset();
   cMCRlightWeight_Jets[iJetCut]= (TH1F*) h000->Clone("cMCRlightWeight_" + JetCuts[iJetCut]); cMCRlightWeight_Jets[iJetCut]->Reset();

   SF_Jets[iJetCut]= (TH1F*) h000->Clone("SF_" + JetCuts[iJetCut]); SF_Jets[iJetCut]->Reset();
   cSF_Jets[iJetCut]= (TH1F*) h000->Clone("cSF_" + JetCuts[iJetCut]); cSF_Jets[iJetCut]->Reset();
   DNegtag_Jets[iJetCut]= (TH1F*) h000->Clone("DNegtag_" + JetCuts[iJetCut]); DNegtag_Jets[iJetCut]->Reset();
   cDNegtag_Jets[iJetCut]= (TH1F*) h000->Clone("cDNegtag_" + JetCuts[iJetCut]); cDNegtag_Jets[iJetCut]->Reset();
   DMistag_Jets[iJetCut]= (TH1F*) h000->Clone("DMistag_" + JetCuts[iJetCut]); DMistag_Jets[iJetCut]->Reset();
   cDMistag_Jets[iJetCut]= (TH1F*) h000->Clone("cDMistag_" + JetCuts[iJetCut]); cDMistag_Jets[iJetCut]->Reset();

   TH1F *h00 = NULL;
   TH1F *h01 = NULL;
   TH1F *h02 = NULL;
   TH1F *h03 = NULL;
   TH1F *h04 = NULL;
   TH1F *h09 = NULL;
 
   TH1F *h10 = NULL;
   TH1F *h11 = NULL;
   TH1F *h12 = NULL;
   TH1F *h13 = NULL;
   TH1F *h14 = NULL;
   TH1F *h19 = NULL;

   TH1F *h20 = NULL;
   TH1F *h21 = NULL;
   TH1F *h24 = NULL;
   TH1F *h29 = NULL;
 
   TFile *fMCTriggered = NULL;
   for( unsigned long iPthat = 0; iPthat < PthatRanges.size(); ++iPthat )
   {
     //delete fMCTriggered;

     TString nameTmp = PthatRanges.at(iPthat) ;
     //nameTmp.ReplaceAll("/","") ;
     //nameTmp.ReplaceAll("-", "_") ;
     //nameTmp = "QCD_" + nameTmp + "_TuneCUETP8M1_13TeV_pythia8/" ;
     fMCTriggered = new TFile(DirOfMCFile+nameTmp+RootOfMCFile+JetCuts[iJetCut] +"_"+CutName+".root");
     //fMCTriggered = new TFile(DirOfMCFile+PthatRanges.at(iPthat)+RootOfMCFile+JetCuts[iJetCut] +"_"+CutName+".root");
     if( fMCTriggered == NULL )
     {
       std::cout << "Error opening " << DirOfMCFile<<PthatRanges.at(iPthat)<<RootOfMCFile << JetCuts[iJetCut] << "_"<<CutName<<".root" << endl;
       getchar();
       return NULL;
     }// else std::cout << DirOfMCFile << JetCuts[iJetCut] <<  "_" << CutName << ".root opened " << endl;
     fMCTriggered->cd();
  
     if( iPthat == 0 )
     {
       TString tNameSufix = EtaRange+"_MCTrig_"+JetCuts[iJetCut];

       GetLowStatBins(fMCTriggered,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);


       h00 = FindRebinSetLowStatToZero(fMCTriggered,"hData_JetPt",tNameSufix,MinEta,MaxEta,LowStatBins);
       h01 = FindRebinSetLowStatToZero(fMCTriggered,"hUDSFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);   // uds-jet
       h02 = FindRebinSetLowStatToZero(fMCTriggered,"hCFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);     // c-jet
       h03 = FindRebinSetLowStatToZero(fMCTriggered,"hBFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);     // b-jet
       h04 = FindRebinSetLowStatToZero(fMCTriggered,"hGluonFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // g-jet
       h09 = FindRebinSetLowStatToZero(fMCTriggered,"hLightFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // udsg-jet

       h10 = FindRebinSetLowStatToZero(fMCTriggered,"hData_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta,LowStatBins);
       h11 = FindRebinSetLowStatToZero(fMCTriggered,"hUDSFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);   // uds neg
       h12 = FindRebinSetLowStatToZero(fMCTriggered,"hCFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);     // c neg
       h13 = FindRebinSetLowStatToZero(fMCTriggered,"hBFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);     // b neg
       h14 = FindRebinSetLowStatToZero(fMCTriggered,"hGluonFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // g neg
       h19 = FindRebinSetLowStatToZero(fMCTriggered,"hLightFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // udsg neg 

       h20 = FindRebinSetLowStatToZero(fMCTriggered,"hData_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);      // all pos
       h21 = FindRebinSetLowStatToZero(fMCTriggered,"hUDSFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);   // uds pos 
       h24 = FindRebinSetLowStatToZero(fMCTriggered,"hGluonFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // g pos
       h29 = FindRebinSetLowStatToZero(fMCTriggered,"hLightFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // udsg pos

     }
     else
     {
       TString tNameSufix = EtaRange+"_"+PthatRanges.at(iPthat)+"_MCTrig_"+JetCuts[iJetCut];

       GetLowStatBins(fMCTriggered,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);


       h00->Add(FindRebinSetLowStatToZero(fMCTriggered,"hData_JetPt",tNameSufix,MinEta,MaxEta,LowStatBins));

       h01->Add(FindRebinSetLowStatToZero(fMCTriggered,"hUDSFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));   // uds-jet
       h02->Add(FindRebinSetLowStatToZero(fMCTriggered,"hCFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));     // c-jet
       h03->Add(FindRebinSetLowStatToZero(fMCTriggered,"hBFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));     // b-jet
       h04->Add(FindRebinSetLowStatToZero(fMCTriggered,"hGluonFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins)); // g-jet
       h09->Add(FindRebinSetLowStatToZero(fMCTriggered,"hLightFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins)); // udsg-jet

       h10->Add(FindRebinSetLowStatToZero(fMCTriggered,"hData_NegTag"+CutStrength+"_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins ));
       h11->Add(FindRebinSetLowStatToZero(fMCTriggered,"hUDSFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       h12->Add(FindRebinSetLowStatToZero(fMCTriggered,"hCFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       h13->Add(FindRebinSetLowStatToZero(fMCTriggered,"hBFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       h14->Add(FindRebinSetLowStatToZero(fMCTriggered,"hGluonFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       h19->Add(FindRebinSetLowStatToZero(fMCTriggered,"hLightFlav_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));

       h20->Add(FindRebinSetLowStatToZero(fMCTriggered,"hData_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       h21->Add(FindRebinSetLowStatToZero(fMCTriggered,"hUDSFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       h24->Add(FindRebinSetLowStatToZero(fMCTriggered,"hGluonFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       h29->Add(FindRebinSetLowStatToZero(fMCTriggered,"hLightFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
     }
     h00->SetDirectory(0);
     h01->SetDirectory(0);
     h02->SetDirectory(0);
     h03->SetDirectory(0);
     h04->SetDirectory(0);
     h09->SetDirectory(0);

 
     h10->SetDirectory(0);
     h11->SetDirectory(0);
     h12->SetDirectory(0);
     h13->SetDirectory(0);
     h14->SetDirectory(0);
     h19->SetDirectory(0);

     h20->SetDirectory(0);
     h21->SetDirectory(0);
     h24->SetDirectory(0);
     h29->SetDirectory(0);
     fMCTriggered->Close();

     delete fMCTriggered;
   }

   if( h00 == NULL ){ std::cout << "h00 NULL " << endl; getchar();}
   if( h01 == NULL ){ std::cout << "h01 NULL " << endl; getchar();}
   if( h02 == NULL ){ std::cout << "h02 NULL " << endl; getchar();}
   if( h03 == NULL ){ std::cout << "h03 NULL " << endl; getchar();}
   if( h04 == NULL ){ std::cout << "h04 NULL " << endl; getchar();}
   if( h09 == NULL ){ std::cout << "h09 NULL " << endl; getchar();}
  
   if( h10 == NULL ){ std::cout << "h10 NULL " << endl; getchar();}
   if( h11 == NULL ){ std::cout << "h11 NULL " << endl; getchar();}
   if( h12 == NULL ){ std::cout << "h12 NULL " << endl; getchar();}
   if( h13 == NULL ){ std::cout << "h13 NULL " << endl; getchar();}
   if( h14 == NULL ){ std::cout << "h14 NULL " << endl; getchar();}
   if( h19 == NULL ){ std::cout << "h19 NULL " << endl; getchar();}
  
   if( h20 == NULL ){ std::cout << "h20 NULL " << endl; getchar();}
   if( h21 == NULL ){ std::cout << "h21 NULL " << endl; getchar();}
   if( h24 == NULL ){ std::cout << "h24 NULL " << endl; getchar();}
   if( h29 == NULL ){ std::cout << "h29 NULL " << endl; getchar();}

   // Triggered MC with Track History
   TH1F *g00 = NULL;
   TH1F *g09 = NULL;
   TH1F *g10 = NULL;
   TH1F *g29 = NULL;
   TH1F *g100 = NULL;
   TH1F *g109 = NULL;
   TH1F *g110 = NULL;
   TH1F *g129 = NULL;

   TH1F *g200 = NULL;
   TH1F *g209 = NULL;
   TH1F *g210 = NULL;
   TH1F *g229 = NULL;

   TH1F *g300 = NULL;
   TH1F *g309 = NULL;
   TH1F *g310 = NULL;
   TH1F *g329 = NULL;

   TH1F *g400 = NULL;
   TH1F *g409 = NULL;
   TH1F *g410 = NULL;
   TH1F *g429 = NULL;

   TH1F *g529 = NULL;
   TH1F *g519 = NULL;
 
   TFile *fMCTriggeredTH = NULL;
   for( unsigned long iPthat = 0; iPthat < PthatRanges.size(); ++iPthat )
   {
     //delete fMCTriggeredTH;
     TString nameTmp = PthatRanges.at(iPthat) ;
     //nameTmp.ReplaceAll("/","") ;
     //nameTmp.ReplaceAll("-", "_") ;
     //nameTmp = "QCD_" + nameTmp + "_TuneCUETP8M1_13TeV_pythia8/" ;
     //fMCTriggeredTH = new TFile(DirOfMCFileTH+PthatRanges.at(iPthat)+RootOfMCFile+JetCuts[iJetCut]+"_"+CutName+".root");
     fMCTriggeredTH = new TFile(DirOfMCFileTH+nameTmp+RootOfMCFile+JetCuts[iJetCut]+"_"+CutName+".root");
     if( fMCTriggeredTH == NULL )
     {
       std::cout << "Error opening " << DirOfMCFileTH<<PthatRanges.at(iPthat)<<RootOfMCFile << JetCuts[iJetCut] << "_"<<CutName<<".root" << endl;
       getchar();
       return NULL;
     }// else std::cout << DirOfMCFileTH << JetCuts[iJetCut] <<  "_" << CutName << ".root opened " << endl;
     fMCTriggeredTH->cd();
  
     if( iPthat == 0 )
     {
       GetLowStatBins(fMCTriggeredTH,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);

       TString tNameSufix = EtaRange+"_MCTHTrig_"+JetCuts[iJetCut];

       g00 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hData_JetPt",tNameSufix,MinEta,MaxEta,LowStatBins);
       g09 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins );
       g10 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hData_NegTag"+CutStrength+"_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins );
       g29 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_PosTag"+CutStrength+"_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins );

       g100 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_K0s_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g109 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_K0s_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g110 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_K0s_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g129 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_K0s_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);

       g200 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Fak_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g209 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Fak_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g210 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Fak_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g229 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Fak_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);

       g300 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Gam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g309 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Gam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g310 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Gam_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g329 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Gam_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);

       g400 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Lam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g409 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Lam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g410 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Lam_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g429 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Lam_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);

       g529 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Oth_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
       g519 = FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Oth_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins);
     } else
     {
       GetLowStatBins(fMCTriggeredTH,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);

       TString tNameSufix = EtaRange+"_"+PthatRanges.at(iPthat)+"_MCTHTrig_"+JetCuts[iJetCut];

       g00->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hData_JetPt",tNameSufix,MinEta,MaxEta,LowStatBins));
       g09->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins ));
       g10->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hData_NegTag"+CutStrength+"_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins ));
       g29->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_PosTag"+CutStrength+"_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins ));

       g100->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_K0s_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g109->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_K0s_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g110->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_K0s_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g129->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_K0s_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));

       g200->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Fak_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g209->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Fak_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g210->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Fak_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g229->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Fak_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));

       g300->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Gam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g309->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Gam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g310->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Gam_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g329->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Gam_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));

       g400->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Lam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g409->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Lam_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g410->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hAllFlav_Lam_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g429->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Lam_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));

       g529->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Oth_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
       g519->Add(FindRebinSetLowStatToZero(fMCTriggeredTH,"hLightFlav_Oth_NegTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins));
     }
     g00->SetDirectory(0);
     g09->SetDirectory(0);
     g10->SetDirectory(0);
     g29->SetDirectory(0);
     g100->SetDirectory(0);
     g109->SetDirectory(0);
     g110->SetDirectory(0);
     g129->SetDirectory(0);

     g200->SetDirectory(0);
     g209->SetDirectory(0);
     g210->SetDirectory(0);
     g229->SetDirectory(0);

     g300->SetDirectory(0);
     g309->SetDirectory(0);
     g310->SetDirectory(0);
     g329->SetDirectory(0);

     g400->SetDirectory(0);
     g409->SetDirectory(0);
     g410->SetDirectory(0);
     g429->SetDirectory(0);

     g529->SetDirectory(0);
     g519->SetDirectory(0);
     fMCTriggeredTH->Close();
     delete fMCTriggeredTH;
   }

   if( g00 == NULL ){ std::cout << "g00 NULL " << endl; getchar();}
   if( g09 == NULL ){ std::cout << "g09 NULL " << endl; getchar();}
   if( g10 == NULL ){ std::cout << "g10 NULL " << endl; getchar();}
   if( g29 == NULL ){ std::cout << "g29 NULL " << endl; getchar();}
   if( g100 == NULL ){ std::cout << "g100 NULL " << endl; getchar();}
   if( g109 == NULL ){ std::cout << "g109 NULL " << endl; getchar();}
   if( g110 == NULL ){ std::cout << "g110 NULL " << endl; getchar();}
   if( g129 == NULL ){ std::cout << "g129 NULL " << endl; getchar();}
   if( g200 == NULL ){ std::cout << "g200 NULL " << endl; getchar();}
   if( g209 == NULL ){ std::cout << "g209 NULL " << endl; getchar();}
   if( g210 == NULL ){ std::cout << "g210 NULL " << endl; getchar();}
   if( g229 == NULL ){ std::cout << "g229 NULL " << endl; getchar();}
   if( g300 == NULL ){ std::cout << "g300 NULL " << endl; getchar();}
   if( g309 == NULL ){ std::cout << "g309 NULL " << endl; getchar();}
   if( g310 == NULL ){ std::cout << "g310 NULL " << endl; getchar();}
   if( g329 == NULL ){ std::cout << "g329 NULL " << endl; getchar();}
   if( g400 == NULL ){ std::cout << "g400 NULL " << endl; getchar();}
   if( g409 == NULL ){ std::cout << "g409 NULL " << endl; getchar();}
   if( g410 == NULL ){ std::cout << "g410 NULL " << endl; getchar();}
   if( g429 == NULL ){ std::cout << "g429 NULL " << endl; getchar();}
   if( g519 == NULL ){ std::cout << "g519 NULL " << endl; getchar();}
   if( g529 == NULL ){ std::cout << "g529 NULL " << endl; getchar();}

   // Triggered MC with average Pile Up
   TH1F* a00 = NULL;
   TH1F* a09 = NULL;
   TH1F* a10 = NULL;
   TH1F* a29 = NULL;
   TFile *fMCTriggered_avPU = NULL;
   for( unsigned long iPthat = 0; iPthat < PthatRanges.size(); ++iPthat )
   {

     TString nameTmp = PthatRanges.at(iPthat) ;
     //nameTmp.ReplaceAll("/","") ;
     //nameTmp.ReplaceAll("-", "_") ;
     //nameTmp = "QCD_" + nameTmp + "_TuneCUETP8M1_13TeV_pythia8/" ;
     fMCTriggered_avPU = new TFile(DirOfMCFile+nameTmp+RootOfMCFile+JetCuts[iJetCut] +"_"+CutName+".root");
     if( fMCTriggered_avPU == NULL )
     {
       std::cout << "Error opening " << DirOfMCFile<<PthatRanges.at(iPthat)<<RootOfMCFile << JetCuts[iJetCut] << "_"<<CutName<<".root" << endl;
       getchar();
       return NULL;
     }// else std::cout << DirOfMCFile << JetCuts[iJetCut] <<  "_" << CutName << ".root opened " << endl;
     fMCTriggered_avPU->cd();
  
     if( iPthat == 0 )
     {
       GetLowStatBins(fMCTriggered_avPU,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);

       TString tNameSufix = EtaRange+"_MCTHTrig_"+JetCuts[iJetCut];

       a00 = FindRebinSetLowStatToZero(fMCTriggered_avPU,"hData_JetPt",tNameSufix,MinEta,MaxEta,LowStatBins);
       a09 = FindRebinSetLowStatToZero(fMCTriggered_avPU,"hLightFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // udsg-jet
       a10 = FindRebinSetLowStatToZero(fMCTriggered_avPU,"hData_NegTag"+CutStrength+"_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins );
       a29 = FindRebinSetLowStatToZero(fMCTriggered_avPU,"hLightFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins); // udsg pos
     } else
     {
       GetLowStatBins(fMCTriggered_avPU,"hData_JetPt", MinEta, MaxEta, xSigma, LowStatBins);

       TString tNameSufix = EtaRange+"_"+PthatRanges.at(iPthat)+"_MCTHTrig_"+JetCuts[iJetCut];

       a00->Add(FindRebinSetLowStatToZero(fMCTriggered_avPU,"hData_JetPt",tNameSufix,MinEta,MaxEta,LowStatBins));
       a09->Add(FindRebinSetLowStatToZero(fMCTriggered_avPU,"hLightFlav_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins)); // udsg-jet
       a10->Add(FindRebinSetLowStatToZero(fMCTriggered_avPU,"hData_NegTag"+CutStrength+"_JetPt", tNameSufix,MinEta,MaxEta, LowStatBins ));
       a29->Add(FindRebinSetLowStatToZero(fMCTriggered_avPU,"hLightFlav_PosTag"+CutStrength+"_JetPt",tNameSufix,MinEta,MaxEta, LowStatBins)); // udsg pos
     }
     a00->SetDirectory(0);
     a09->SetDirectory(0);
     a10->SetDirectory(0);
     a29->SetDirectory(0);
     fMCTriggered_avPU->Close();
     delete fMCTriggered_avPU;
   }

   if( a00 == NULL ){ std::cout << "a00 NULL " << endl; getchar();}
   if( a09 == NULL ){ std::cout << "a09 NULL " << endl; getchar();}
   if( a10 == NULL ){ std::cout << "a10 NULL " << endl; getchar();}
   if( a29 == NULL ){ std::cout << "a29 NULL " << endl; getchar();}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Mistags and Negative Tags

// without V0 reweighting
   TH1F* h4= (TH1F*) h000->Clone("h4"); h4->Reset();
         h4->Divide(h29,h09,1,1,"B"); 
   TH1F* h5= (TH1F*) h000->Clone("h5"); h5->Reset();
         h5->Divide(h10,h00,1,1,"B");
   TH1F* hRlight= (TH1F*) h000->Clone("hRlight"); hRlight->Reset();
         hRlight->Divide(h4,h5,1,1);

// with V0 reweighting

   float nh00 = h00->Integral(); 
   float nh09 = h09->Integral(); 
   float nh10 = h10->Integral(); 
   float nh29 = h29->Integral(); 

   float ng00 = g00->Integral(); 
   float ng09 = g09->Integral(); 
   float ng10 = g10->Integral(); 
   float ng29 = g29->Integral(); 

   float TPscale00 = nh00/ng00;
   float TPscale09 = nh09/ng09;
   float TPscale10 = nh10/ng10;
   float TPscale29 = nh29/ng29;

   float e00 = sqrt(1/nh00 + 1/ng00)*TPscale00;
   float e09 = sqrt(1/nh09 + 1/ng09)*TPscale09;
   float e10 = sqrt(1/nh10 + 1/ng10)*TPscale10;
   float e29 = sqrt(1/nh29 + 1/ng29)*TPscale29;

   float e_Consistence_externe = TPscale00*TPscale00 + TPscale09*TPscale09 + TPscale10*TPscale10 + TPscale29*TPscale29;
   e_Consistence_externe -= (TPscale00+TPscale09+TPscale10+TPscale29)*(TPscale00+TPscale09+TPscale10+TPscale29)/4.;
   e_Consistence_externe = sqrt( e_Consistence_externe/3 );

   std::cout << "rapport 00 = " << TPscale00 << " +- " << e00 << endl;
   std::cout << "rapport 09 = " << TPscale09 << " +- " << e09 << endl;
   std::cout << "rapport 10 = " << TPscale10 << " +- " << e10 << endl;
   std::cout << "rapport 29 = " << TPscale29 << " +- " << e29 << endl;
   std::cout << "Error by external consistency is " << e_Consistence_externe << endl;
   //getchar();

   H00->Sumw2(); 
   H10->Sumw2(); 
   H29->Sumw2(); 
   H09->Sumw2(); 
   H00->Add(h00,g100,1,TPscale09*(SFk0-1)); 
   H10->Add(h10,g110,1,TPscale09*(SFk0-1)); 
   H29->Add(h29,g129,1,TPscale09*(SFk0-1)); 
   H09->Add(h09,g109,1,TPscale09*(SFk0-1)); 
   H00->Add(H00,g400,1,TPscale09*(SFla-1)); 
   H10->Add(H10,g410,1,TPscale09*(SFla-1)); 
   H29->Add(H29,g429,1,TPscale09*(SFla-1)); 
   H09->Add(H09,g409,1,TPscale09*(SFla-1)); 

   MCMistag_Jets[iJetCut]->Divide(H29,H09,1,1,"B");
   MCNegtag_Jets[iJetCut]->Divide(H10,H00,1,1,"B");
   MCRlightWeight_Jets[iJetCut]->Divide(MCMistag_Jets[iJetCut],MCNegtag_Jets[iJetCut],1,1);
   cMCRlightWeight_Jets[iJetCut]->Divide(MCMistag_Jets[iJetCut],MCNegtag_Jets[iJetCut],1,1);

// Correction to MC without trigger 
   for (int i=0; i<HistPtNbin+1; i++)
   {
     float y0 = k4->GetBinContent(i);
     float y1 = h4->GetBinContent(i);
     y = cMCRlightWeight_Jets[iJetCut]->GetBinContent(i);
     ey = cMCRlightWeight_Jets[iJetCut]->GetBinError(i);
     if (y1 > 0.) cMCRlightWeight_Jets[iJetCut]->SetBinContent(i,y * y0 / y1);
     if (y1 > 0.) cMCRlightWeight_Jets[iJetCut]->SetBinError(i,ey * y0 / y1);

     y = MCMistag_Jets[iJetCut]->GetBinContent(i);
     ey = MCMistag_Jets[iJetCut]->GetBinError(i);
     if (y1 > 0.) cMCMistag_Jets[iJetCut]->SetBinContent(i,y * y0 / y1);
     if (y1 > 0.) cMCMistag_Jets[iJetCut]->SetBinError(i,ey * y0 / y1);

     y = MCNegtag_Jets[iJetCut]->GetBinContent(i);
     ey = MCNegtag_Jets[iJetCut]->GetBinError(i);
     y0 = k5->GetBinContent(i);
     y1 = h5->GetBinContent(i);
     if (y1 > 0.) cMCNegtag_Jets[iJetCut]->SetBinContent(i,y * y0 / y1);
     if (y1 > 0.) cMCNegtag_Jets[iJetCut]->SetBinError(i,ey * y0 / y1);
   }

// with average pile-up 
   TH1F* a4= (TH1F*) h000->Clone("a4"); a4->Reset();
         a4->Divide(a29,a09,1,1,"B"); 
   TH1F* a5= (TH1F*) h000->Clone("a5"); a5->Reset();
         a5->Divide(a10,a00,1,1,"B");
   TH1F* Ra= (TH1F*) h000->Clone("Ra"); Ra->Reset();
         Ra->Divide(a4,a5,1,1);

// syst pu
   for (int i=0; i<HistPtNbin+1; i++)
   {
     float dy = 0.;
     if ( hRlight->GetBinContent(i) > 0. )
          dy = TMath::Abs( Ra->GetBinContent(i)/hRlight->GetBinContent(i) - 1. );
     puMis_Jets[iJetCut]->SetBinContent(i,dy);
     dy = 0.;
     if ( h5->GetBinContent(i) > 0. )
 	  dy = TMath::Abs( a5->GetBinContent(i)/h5->GetBinContent(i) - 1. );
     puSF_Jets[iJetCut]->SetBinContent(i,dy);
   }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Flavour fraction

   TH1F* k00= (TH1F*) h000->Clone("k00"); k00->Reset();
   TH1F* k01= (TH1F*) h000->Clone("k01"); k01->Reset();
   TH1F* k02= (TH1F*) h000->Clone("k02"); k02->Reset();
   TH1F* k03= (TH1F*) h000->Clone("k03"); k03->Reset();
   TH1F* k04= (TH1F*) h000->Clone("k04"); k04->Reset();
   TH1F* k11= (TH1F*) h000->Clone("k11"); k11->Reset();
   TH1F* k12= (TH1F*) h000->Clone("k12"); k12->Reset();
   TH1F* k13= (TH1F*) h000->Clone("k13"); k13->Reset();
   TH1F* k14= (TH1F*) h000->Clone("k14"); k14->Reset();
   TH1F* k21= (TH1F*) h000->Clone("k21"); k21->Reset();
   TH1F* k24= (TH1F*) h000->Clone("k24"); k24->Reset();

   k00->Divide(h04,h01,1,1);
   k01->Divide(h01,h00,1,1,"B");
   k02->Divide(h02,h00,1,1,"B");
   k03->Divide(h03,h00,1,1,"B");
   k04->Divide(h04,h00,1,1,"B");
   k11->Divide(h11,h10,1,1,"B");
   k12->Divide(h12,h10,1,1,"B");
   k13->Divide(h13,h10,1,1,"B");
   k14->Divide(h14,h10,1,1,"B");
   k21->Divide(h21,h29,1,1,"B");
   k24->Divide(h24,h29,1,1,"B");

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Systematics

// Compute all relative errors on RLight
// systematics due to flavour fraction
   TH1F* bb= (TH1F*) h000->Clone("bb"); bb->Reset();
   for (int i=0; i<HistPtNbin+1; i++)
   {
     y = k13->GetBinContent(i) - k03->GetBinContent(i);
     y = Xbb * TMath::Abs( y );	//from eq7 of AN-10-382
     ey = Xbb * k13->GetBinError(i);
     bb->SetBinContent(i,y);
     bb->SetBinError(i,ey);
   }

   TH1F* cc= (TH1F*) h000->Clone("cc"); cc->Reset();
   for (int i=0; i<HistPtNbin+1; i++)
   {
     y = k12->GetBinContent(i) - k02->GetBinContent(i);
     y = Xcc * TMath::Abs( y ); //from eq6 of AN-10-382
     ey = Xcc * k12->GetBinError(i);
     cc->SetBinContent(i,y);
     cc->SetBinError(i,ey);
   }

   bc_Jets[iJetCut]->Add(bb,cc);

   for (int i=0; i<HistPtNbin+1; i++)
   {
     float frac = k00->GetBinContent(i);
     float lneg = k11->GetBinContent(i);
     float lpos = k21->GetBinContent(i);
     float gneg = k14->GetBinContent(i);
     float gpos = k24->GetBinContent(i);
     float yl = frac * (lpos - lneg);
     float yg = gpos - gneg;
     y = Xgg * TMath::Abs( yg - yl ); //from eq8 of AN-10-382
     float elneg = k11->GetBinError(i);
     float elpos = k21->GetBinError(i);
     float egneg = k14->GetBinError(i);
     float egpos = k24->GetBinError(i);
     float eyneg = (egneg - frac*elneg) * (egneg - frac*elneg);
     float eypos = (egpos - frac*elpos) * (egpos - frac*elpos);
     ey = Xgg * TMath::Sqrt( eyneg + eypos );
     gg_Jets[iJetCut]->SetBinContent(i,y);
     gg_Jets[iJetCut]->SetBinError(i,ey);
   }

// systematics due to sign flip
   TH1F* h1= (TH1F*) h000->Clone("h1"); h1->Reset();
         h1->Divide(h20,h10,1,1);
   TH1F* h2= (TH1F*) h000->Clone("h2"); h2->Reset();
         h2->Divide(h19,h29,1,1);
         sf_Jets[iJetCut]->Add(h1,h2,Xip,-Xip); //from eq10 of AN-10-382

// systematics due to K0s
   TH1F* h71= (TH1F*) h000->Clone("h71"); h71->Reset();
         h71->Divide(g110,g10,SFk0,1,"B");
   TH1F* h72= (TH1F*) h000->Clone("h72"); h72->Reset();
         h72->Divide(g129,g29,SFk0,1,"B");
   TH1F* h73= (TH1F*) h000->Clone("h73"); h73->Reset();
         h73->Divide(g109,g09,SFk0,1,"B");
   TH1F* h74= (TH1F*) h000->Clone("h74"); h74->Reset();
         h74->Divide(g100,g00,SFk0,1,"B");
   TH1F* ka= (TH1F*) h000->Clone("ka"); ka->Reset();
   for (int i=0; i<h000->GetNbinsX()+1; i++)
   {
     y    = h72->GetBinContent(i) - h71->GetBinContent(i)
          - h73->GetBinContent(i) + h74->GetBinContent(i);
     y = Xk0 * TMath::Abs( y ); //from eq9 of AN-10-382
     ey = h71->GetBinError(i)*h71->GetBinError(i) 
        + h72->GetBinError(i)*h72->GetBinError(i);
     ey = Xk0 * TMath::Sqrt( ey );
     ka->SetBinContent(i,y);
     ka->SetBinError(i,ey);
   }

// systematics due to Lambda
   TH1F* h771= (TH1F*) h000->Clone("h771"); h771->Reset();
         h771->Divide(g410,g10,SFla,1,"B");
   TH1F* h772= (TH1F*) h000->Clone("h772"); h772->Reset();
         h772->Divide(g429,g29,SFla,1,"B");
   TH1F* h773= (TH1F*) h000->Clone("h773"); h773->Reset();
         h773->Divide(g409,g09,SFla,1,"B");
   TH1F* h774= (TH1F*) h000->Clone("h774"); h774->Reset();
         h774->Divide(g400,g00,SFla,1,"B");
   TH1F* la= (TH1F*) h000->Clone("la"); la->Reset();
   for (int i=0; i<h000->GetNbinsX()+1; i++)
   {
     y    = h772->GetBinContent(i) - h771->GetBinContent(i)
          - h773->GetBinContent(i) + h774->GetBinContent(i);
     y = Xla * TMath::Abs( y ); //from eq9 of AN-10-382
     ey = h771->GetBinError(i)*h771->GetBinError(i) 
        + h772->GetBinError(i)*h772->GetBinError(i);
     ey = Xla * TMath::Sqrt( ey );
     la->SetBinContent(i,y);
     la->SetBinError(i,ey);
   }

// systematics due to fakes
/*
   float n8 = g229->Integral(); //GetEntries();
   float n10= g529->Integral() - g519->Integral(); //g529->GetEntries() - g519->GetEntries();
   // !!!!!! ==== En principe il faudrait utiliser les intégrales ci-dessus ====== !!!!!!!!
   // Mais la methode de reweighting pese avec un facteur important des evts tres rares  
   // conduisant ainsi a des erreurs enormes. Du coup, l'utilisation du nombre d'entrees (cf ci dessous)
   // est plus precise !
*/

   float n8 = g229->GetEntries();
   float n10= g529->GetEntries() - g519->GetEntries();
   float SFfk = 1.;
   if ( n8*n10 > 0 ) SFfk = 1. + n10 / n8; 
   std::cout << " SFfk " << SFfk << "\t" << n8 << "\t" << n10 << endl;
   std::cout << endl;
   if( ( ((g529->Integral()) - (g519->Integral())) * (g229->Integral())) > 0 ) std::cout << 1+(g529->Integral() - g519->Integral())/g229->Integral() << endl;
   else std::cout << 1 << endl;
   TH1F* h81= (TH1F*) h000->Clone("h81"); h81->Reset();
         h81->Divide(g210,g10,SFfk,1,"B");

   TH1F* h82= (TH1F*) h000->Clone("h82"); h82->Reset();
         h82->Divide(g229,g29,SFfk,1,"B");
   TH1F* h83= (TH1F*) h000->Clone("h83"); h83->Reset();
         h83->Divide(g209,g09,SFfk,1,"B");
   TH1F* h84= (TH1F*) h000->Clone("h84"); h84->Reset();
         h84->Divide(g200,g00,SFfk,1,"B");
   for (int i=0; i<h000->GetNbinsX()+1; i++)
   {
     y    = h82->GetBinContent(i) - h81->GetBinContent(i)
          - h83->GetBinContent(i) + h84->GetBinContent(i);
     y = Xfk * TMath::Abs( y );
     ey = h81->GetBinError(i)*h81->GetBinError(i) 
        + h82->GetBinError(i)*h82->GetBinError(i);
     ey = Xfk * TMath::Sqrt( ey );
     fk_Jets[iJetCut]->SetBinContent(i,y);
     fk_Jets[iJetCut]->SetBinError(i,ey);
   }

// systematics due to 2ndary interactions
   TH1F* h91= (TH1F*) h000->Clone("h91"); h91->Reset();
         h91->Divide(g310,g10,1,1,"B");
   TH1F* h92= (TH1F*) h000->Clone("h92"); h92->Reset();
         h92->Divide(g329,g29,1,1,"B");
   TH1F* h93= (TH1F*) h000->Clone("h93"); h93->Reset();
         h93->Divide(g309,g09,1,1,"B");
   TH1F* h94= (TH1F*) h000->Clone("h94"); h94->Reset();
         h94->Divide(g300,g00,1,1,"B");
   TH1F* ga= (TH1F*) h000->Clone("ga"); ga->Reset();
   for (int i=0; i<h000->GetNbinsX()+1; i++)
   {
     y    = h92->GetBinContent(i) - h91->GetBinContent(i)
          - h93->GetBinContent(i) + h94->GetBinContent(i);
     y = Xga * TMath::Abs( y );
     ey = h91->GetBinError(i)*h91->GetBinError(i) 
        + h92->GetBinError(i)*h92->GetBinError(i);
     ey = Xga * TMath::Sqrt( ey );
     ga->SetBinContent(i,y);
     ga->SetBinError(i,ey);
   }

   for (int i=0; i<h000->GetNbinsX()+1; i++)
   {
     y = ka->GetBinContent(i)*ka->GetBinContent(i) 
       + la->GetBinContent(i)*la->GetBinContent(i)
       + ga->GetBinContent(i)*ga->GetBinContent(i);
     y = TMath::Sqrt( y );
     ey = ka->GetBinError(i)*ka->GetBinError(i) 
        + la->GetBinError(i)*la->GetBinError(i)
        + ga->GetBinError(i)*ga->GetBinError(i);
     ey = TMath::Sqrt( ey );
     v0_Jets[iJetCut]->SetBinContent(i,y);
     v0_Jets[iJetCut]->SetBinError(i,ey);
   }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Mistag in Data = all neg data * Rlight

// PVH changes up to here !
   TFile *fDATA = NULL;
     //fDATA = new TFile(DirOfDataFile+RootOfDataFile+JetCuts[iJetCut]+"_"+CutName+"_WeightStudy.root");
     fDATA = new TFile(DirOfDataFile+RootOfDataFile+JetCuts[iJetCut]+"_"+CutName+".root");

   if( fDATA == NULL )
   {
     std::cout << "Error opening it" << endl;
     getchar();
   } else
   {
     //std::cout << "Succeded to open " << DirOfDataFile << RootOfDataFile << JetCuts[iJetCut] << "_" << CutName << "_WeightStudy.root" << endl;
   }
   fDATA->cd();


   D00_Jets[iJetCut] = FindRebin(fDATA,"hData_JetPt",EtaRange+"_MCavgPUnoTrigger_"+JetCuts[iJetCut],MinEta,MaxEta);
   D10_Jets[iJetCut] = FindRebin(fDATA,"hData_NegTag"+CutStrength+"_JetPt",EtaRange+"_MCavgPUnoTrigger_"+JetCuts[iJetCut],MinEta,MaxEta);

   D00_Jets[iJetCut]->SetDirectory(0);
   D10_Jets[iJetCut]->SetDirectory(0);

   if( D00_Jets[iJetCut] == NULL ){ std::cout << "D00 NULL " << iJetCut << endl; getchar();}
   if( D10_Jets[iJetCut] == NULL ){ std::cout << "D10 NULL " << iJetCut << endl; getchar();}

   DNegtag_Jets[iJetCut]->Divide(D10_Jets[iJetCut],D00_Jets[iJetCut],1,1,"B");
   for (int i=0; i<h000->GetNbinsX()+1; i++)
   {
     if( DNegtag_Jets[iJetCut]->GetBinContent(i) == 1 ) DNegtag_Jets[iJetCut]->SetBinError(i,1000.);
   }
   // Correction to "unbias" the trigger
   for (int i=0; i<h000->GetNbinsX()+1; i++) {
     float y0 = k5->GetBinContent(i);
     float y1 = h5->GetBinContent(i);
     y = DNegtag_Jets[iJetCut]->GetBinContent(i);
     ey = DNegtag_Jets[iJetCut]->GetBinError(i);
     if (y1 > 0.) cDNegtag_Jets[iJetCut]->SetBinContent(i,y * y0 / y1);
     if (y1 > 0.) cDNegtag_Jets[iJetCut]->SetBinError(i,ey * y0 / y1);
   }
   DMistag_Jets[iJetCut]->Multiply(DNegtag_Jets[iJetCut] ,MCRlightWeight_Jets[iJetCut],1,1); 
   for (int i=0; i<h000->GetNbinsX()+1; i++) {
     ey = DNegtag_Jets[iJetCut] ->GetBinError(i);
     DMistag_Jets[iJetCut]->SetBinError(i,MCRlightWeight_Jets[iJetCut]->GetBinContent(i) * ey);
   }
   cDMistag_Jets[iJetCut]->Multiply(DNegtag_Jets[iJetCut] ,cMCRlightWeight_Jets[iJetCut],1,1); 
   for (int i=0; i<h000->GetNbinsX()+1; i++) {
     ey = DNegtag_Jets[iJetCut] ->GetBinError(i);
     cDMistag_Jets[iJetCut]->SetBinError(i,cMCRlightWeight_Jets[iJetCut]->GetBinContent(i) * ey);
   }

// Data/MC Scale Factor

   SF_Jets[iJetCut]->Divide(DNegtag_Jets[iJetCut],MCNegtag_Jets[iJetCut],1,1); // neg tag Data/MCCorrectedForv0
   cSF_Jets[iJetCut]->Divide(cDNegtag_Jets[iJetCut],cMCNegtag_Jets[iJetCut],1,1); // neg tag Data/MCCorrectedForv0
   TH1F* SFneg= (TH1F*) h000->Clone("SFneg"); SFneg->Reset();
         SFneg->Divide(DNegtag_Jets[iJetCut] ,h5,1,1); 

 }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Fraction of each single jet trigger in Data
 TH1F* D00_total= (TH1F*) h000->Clone("D00_total"); D00_total->Reset();
 D00_total->Sumw2(); 
 TH1F* D10_total= (TH1F*) h000->Clone("D10_total"); D10_total->Reset();
 D10_total->Sumw2(); 
 TH1F* frac_Jets[8]  = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
 TH1F* bc= (TH1F*) h000->Clone("bc"); bc->Reset();
 TH1F* gg= (TH1F*) h000->Clone("gg"); gg->Reset();
 TH1F* sf= (TH1F*) h000->Clone("sf"); sf->Reset();
 TH1F* fk= (TH1F*) h000->Clone("fk"); fk->Reset();
 TH1F* v0= (TH1F*) h000->Clone("v0"); v0->Reset();
 TH1F* puMis= (TH1F*) h000->Clone("puMis"); puMis->Reset();
 TH1F* puSF= (TH1F*) h000->Clone("puSF"); puSF->Reset();

 for( int iJetCut = 0; iJetCut < NJetCuts; iJetCut++)
 {
   D00_total->Add(D00_total,D00_Jets[iJetCut],1,1);
   D10_total->Add(D10_total,D10_Jets[iJetCut],1,1);
 }
 TH1F* DNegtag= (TH1F*) h000->Clone("DNegtag"); DNegtag->Reset();
 DNegtag->Divide(D10_total,D00_total,1,1,"B");
 // Weight the syst errors with frac_Jets
 for( int iJetCut = 0; iJetCut < NJetCuts; iJetCut++)
 {
   frac_Jets[iJetCut]= (TH1F*) h000->Clone("frac_"+JetCuts[iJetCut]); frac_Jets[iJetCut]->Reset();
   frac_Jets[iJetCut]->Divide(D00_Jets[iJetCut],D00_total,1,1,"B");

   bc_Jets[iJetCut]->Multiply(bc_Jets[iJetCut],frac_Jets[iJetCut],1,1); 
   bc->Add(bc,bc_Jets[iJetCut],1,1); 

   gg_Jets[iJetCut]->Multiply(gg_Jets[iJetCut],frac_Jets[iJetCut],1,1);
   gg->Add(gg,gg_Jets[iJetCut],1,1); 

   sf_Jets[iJetCut]->Multiply(sf_Jets[iJetCut],frac_Jets[iJetCut],1,1); 
   sf->Add(sf,sf_Jets[iJetCut],1,1); 

   fk_Jets[iJetCut]->Multiply(fk_Jets[iJetCut],frac_Jets[iJetCut],1,1); 
   fk->Add(fk,fk_Jets[iJetCut],1,1); 

   v0_Jets[iJetCut]->Multiply(v0_Jets[iJetCut],frac_Jets[iJetCut],1,1); 
   v0->Add(v0,v0_Jets[iJetCut],1,1); 

   puMis_Jets[iJetCut]->Multiply(puMis_Jets[iJetCut],frac_Jets[iJetCut],1,1); 
   puMis->Add(puMis,puMis_Jets[iJetCut],1,1); 

   puSF_Jets[iJetCut]->Multiply(puSF_Jets[iJetCut],frac_Jets[iJetCut],1,1); 
   puSF->Add(puSF,puSF_Jets[iJetCut],1,1); 
 }

// Recompute the Systematics (Add statistical error from MC without trigger)

 for (int i=0; i<=h000->GetNbinsX(); i++) {
  float ww = puMis->GetBinContent(i);
  float mcstat = StatMis->GetBinContent(i);	// Relative error on Rlight from MC without trigger
  ww = TMath::Sqrt( ww*ww + mcstat*mcstat );
  puMis->SetBinContent(i,ww);
  puMis->SetBinError(i,0.001);
  ww = puSF->GetBinContent(i);
  mcstat = StatSF->GetBinContent(i);		// Relative error on negative tag from MC without trigger
  ww = TMath::Sqrt( ww*ww + mcstat*mcstat );
  puSF->SetBinContent(i,ww);
  puSF->SetBinError(i,0.001);
 }

 TH1F* totMis= (TH1F*) h000->Clone("totMis"); totMis->Reset();
 TH1F* totSF= (TH1F*) h000->Clone("totSF"); totSF->Reset();
 TH1F* totSFCorr= (TH1F*) h000->Clone("totSFCorr"); totSFCorr->Reset();

 for (int i=0; i<h000->GetNbinsX()+1; i++) {
  float ybc = bc->GetBinContent(i);
  float ygg = gg->GetBinContent(i);
  float ysf = sf->GetBinContent(i);
  float yv0 = v0->GetBinContent(i);
  if( Manual_v0_Error ) 
  {
    //TEMP
    if( CutStrength == "L" ) yv0 = 0.02;
    if( CutStrength == "M" ) yv0 = 0.04;
    if( CutStrength == "T" ) yv0 = 0.05;
    //Update results
    //if( CutStrength == "L" ) yv0 = 0.01;
    //if( CutStrength == "M" ) yv0 = 0.02;
    //if( CutStrength == "T" ) yv0 = 0.02;
    v0->SetBinContent(i,yv0);
  }
  float yfk = fk->GetBinContent(i);
  y = TMath::Sqrt( ybc*ybc + ygg*ygg + ysf*ysf + yv0*yv0 + yfk*yfk );
  float yPU = puMis->GetBinContent(i);
  float yPUSF = puSF->GetBinContent(i);

  float eybc = bc->GetBinError(i);
  float eygg = gg->GetBinError(i);
  float eysf = sf->GetBinError(i);
  float eyv0 = v0->GetBinError(i);
  float eyfk = fk->GetBinError(i);
  ey = TMath::Sqrt( eybc*eybc + eygg*eygg + eysf*eysf + eyv0*eyv0 + eyfk*eyfk );

  totMis->SetBinContent(i,TMath::Sqrt( y*y + ErrSysMis*ErrSysMis + yPU*yPU ) );
  totMis->SetBinError(i,0.001);

  totSF->SetBinContent(i,TMath::Sqrt( y*y + ErrSysSF*ErrSysSF + yPUSF*yPUSF ) );
  totSF->SetBinError(i,0.001);

  totSFCorr->SetBinContent(i,TMath::Sqrt( y*y + yPUSF*yPUSF ) );
  totSFCorr->SetBinError(i,0.001);
 }

 float w_Jets[8]={0., 0., 0., 0., 0., 0., 0., 0.};

 

 float ww = 0, xx = 0;
// Reweight the Data

 for (int i=0; i<=h000->GetNbinsX(); i++)
 {
   xx = 0; ww = 0;
   for( int iJetCut = 0; iJetCut < NJetCuts; iJetCut++)
   {
     if( includeCorrupted || i < minCorrupted_i[iJetCut] || i > maxCorrupted_i[iJetCut] )
     {
     w_Jets[iJetCut] = cDMistag_Jets[iJetCut]->GetBinError(i);
     if (cDMistag_Jets[iJetCut]->GetBinError(i) > 0.) w_Jets[iJetCut]  = 1./(w_Jets[iJetCut]*w_Jets[iJetCut]); 
     ww += w_Jets[iJetCut];
     xx += (cDMistag_Jets[iJetCut]->GetBinContent(i) * w_Jets[iJetCut]);
     }
     else
     {
       std::cout << "Excluded bin " << i << " for Jet " << JetCuts[iJetCut] << endl;
     }
   }
   if ( ww > 0. ) xx = xx / ww;
   if ( ww > 0. ) ww = 1. / TMath::Sqrt( ww );
   if ( xx > 0. ) ww = ww / xx;
   float mcstat = StatMis->GetBinContent(i);
   ww = TMath::Sqrt( ww*ww + mcstat*mcstat ) * xx;
   cDMistag->SetBinContent(i,xx);
   cDMistag->SetBinError(i,ww);

   xx = 0; ww = 0;
   for( int iJetCut = 0; iJetCut < NJetCuts; iJetCut++)
   {
     if( includeCorrupted || i < minCorrupted_i[iJetCut] || i > maxCorrupted_i[iJetCut] )
     {
     w_Jets[iJetCut] = cMCMistag_Jets[iJetCut]->GetBinError(i);
     if (cMCMistag_Jets[iJetCut]->GetBinError(i) > 0.) w_Jets[iJetCut]  = 1./(w_Jets[iJetCut]*w_Jets[iJetCut]); 
     ww += w_Jets[iJetCut];
     xx += (cMCMistag_Jets[iJetCut]->GetBinContent(i) * w_Jets[iJetCut]);
     }
   }
   if ( ww > 0. ) xx = xx / ww;
   if ( ww > 0. ) ww = 1. / TMath::Sqrt( ww );
   cMCMistag->SetBinContent(i,xx);
   cMCMistag->SetBinError(i,ww);

   xx = 0; ww = 0;
   for( int iJetCut = 0; iJetCut < NJetCuts; iJetCut++)
   {
     if( includeCorrupted || i < minCorrupted_i[iJetCut] || i > maxCorrupted_i[iJetCut] )
     {
     w_Jets[iJetCut] = cMCNegtag_Jets[iJetCut]->GetBinError(i);
     if (cMCNegtag_Jets[iJetCut]->GetBinError(i) > 0.) w_Jets[iJetCut]  = 1./(w_Jets[iJetCut]*w_Jets[iJetCut]); 
     ww += w_Jets[iJetCut];
     xx += (cMCNegtag_Jets[iJetCut]->GetBinContent(i) * w_Jets[iJetCut]);
     }
   }
   if ( ww > 0. ) xx = xx / ww;
   if ( ww > 0. ) ww = 1. / TMath::Sqrt( ww );
   cMCNegtag->SetBinContent(i,xx);
   cMCNegtag->SetBinError(i,ww);

   xx = 0; ww = 0;
   for( int iJetCut = 0; iJetCut < NJetCuts; iJetCut++)
   {
     if( includeCorrupted || i < minCorrupted_i[iJetCut] || i > maxCorrupted_i[iJetCut] )
     {
     w_Jets[iJetCut] = cDNegtag_Jets[iJetCut]->GetBinError(i);
     if (cDNegtag_Jets[iJetCut]->GetBinError(i) > 0.) w_Jets[iJetCut]  = 1./(w_Jets[iJetCut]*w_Jets[iJetCut]); 
     ww += w_Jets[iJetCut];
     xx += (cDNegtag_Jets[iJetCut]->GetBinContent(i) * w_Jets[iJetCut]);
     }
   }
   if ( ww > 0. ) xx = xx / ww;
   if ( ww > 0. ) ww = 1. / TMath::Sqrt( ww );
   cDNegtag->SetBinContent(i,xx);
   cDNegtag->SetBinError(i,ww);

 }
 cSF->Divide(cDMistag,cMCMistag,1.,1.);

 //Open file to save histos:
 TFile fHis(hisFileName, "recreate");


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// *****************************************************************************
// Define fit functions --------------
 TF1* Fun=NULL;
 TF1* Gun=NULL;
 TF1* Hun=NULL;

 if( CutandStrength == "TCHEL" || CutandStrength == "TCHPL" || CutandStrength == "TCHPM" || CutandStrength == "SSVHEM" || CutandStrength == "SSVHPM" || CutandStrength == "JBPM" )
      Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",ptmin,ptmax);
 else if ( CutandStrength == "CSVL" )
 {
   Fun = new TF1("Fun","[0]*(1+[1]*x+[2]*x*x)/(1+[3]*x)",ptmin,ptmax);
 }
 else if ( CutandStrength == "CSVT" ) Fun = new TF1("Fun","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)",ptmin,ptmax);
 else Fun = new TF1("Fun","[0]+[1]*x+[2]*x*x",ptmin,ptmax);
      Fun->SetLineColor(1);
      Fun->SetLineStyle(1);
      Fun->SetLineWidth((Width_t)0.5);

 if( CutandStrength == "TCHEL" || CutandStrength == "TCHPL" || CutandStrength == "TCHPM" || CutandStrength == "SSVHEM" || CutandStrength == "SSVHPM" || CutandStrength == "JBPM" )
      Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",ptmin,ptmax);
 else if ( CutandStrength == "CSVL" )
 {
   Gun = new TF1("Gun","[0]*(1+[1]*x+[2]*x*x)/(1+[3]*x)",ptmin,ptmax);
 }
 else if ( CutandStrength == "CSVT" ) Gun = new TF1("Gun","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)",ptmin,ptmax);
 else Gun = new TF1("Gun","[0]+[1]*x+[2]*x*x",ptmin,ptmax);
      Gun->SetLineColor(1);
      Gun->SetLineStyle(2);

 if( CutandStrength == "TCHEL" || CutandStrength == "TCHPL" || CutandStrength == "TCHPM" || CutandStrength == "SSVHEM" || CutandStrength == "SSVHPM" || CutandStrength == "JBPM" )
      Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",ptmin,ptmax);
 else if ( CutandStrength == "CSVL" )
 {
   Hun = new TF1("Hun","[0]*(1+[1]*x+[2]*x*x)/(1+[3]*x)",ptmin,ptmax);
 }
 else if ( CutandStrength == "CSVT" ) Hun = new TF1("Hun","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)",ptmin,ptmax);
 else Hun = new TF1("Hun","[0]+[1]*x+[2]*x*x",ptmin,ptmax);
      Hun->SetLineColor(1);
      Hun->SetLineStyle(2);
//---------
 TF1 *FUn = new TF1("FUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",ptmin+0.01,ptmax);
 TF1 *GUn = new TF1("GUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",ptmin+0.01,ptmax);
 TF1 *HUn = new TF1("HUn","[0]+[1]*x+[2]*x*x+[3]*x*x*x",ptmin+0.01,ptmax);

 if ( CutandStrength == "TCHEL" || CutandStrength == "TCHPL" )
 {
    FUn = new TF1("FUn","[0]+[1]*x+[2]*x*x",ptmin+0.01,ptmax);
    GUn = new TF1("GUn","[0]+[1]*x+[2]*x*x",ptmin+0.01,ptmax);
    HUn = new TF1("HUn","[0]+[1]*x+[2]*x*x",ptmin+0.01,ptmax);
 }
 if( CutName == "TCHE" )
 {
    FUn = new TF1("FUn","[0]*(1+[1]*x+[2]*x*x)+[3]*x*x*x/(1+[4]*x)",ptmin+0.01,ptmax);
    GUn = new TF1("GUn","[0]*(1+[1]*x+[2]*x*x)+[3]*x*x*x/(1+[4]*x)",ptmin+0.01,ptmax);
    HUn = new TF1("HUn","[0]*(1+[1]*x+[2]*x*x)+[3]*x*x*x/(1+[4]*x)",ptmin+0.01,ptmax);
 }
//
 if( CutandStrength == "cMVAv2L"  && EtaMin < 1.5)
 {
    FUn = new TF1("FUn","[0]+[1]*x+[2]*x*x+[3]/x",ptmin+0.01,ptmax);
    GUn = new TF1("GUn","[0]+[1]*x+[2]*x*x+[3]/x",ptmin+0.01,ptmax);
    HUn = new TF1("HUn","[0]+[1]*x+[2]*x*x+[3]/x",ptmin+0.01,ptmax);
    //FUn = new TF1("FUn","[0]+[1]*x+[2]*x*x+([4])/(1+exp(([3]+1)*(x-50)))",ptmin+0.01,ptmax);
    //GUn = new TF1("GUn","[0]+[1]*x+[2]*x*x+([4])/(1+exp(([3]+1)*(x-50)))",ptmin+0.01,ptmax);
    //HUn = new TF1("HUn","[0]+[1]*x+[2]*x*x+([4])/(1+exp(([3]+1)*(x-50)))",ptmin+0.01,ptmax);
 }
//
/*
 if( CutandStrength == "CSVv2M" )
 {
    FUn = new TF1("FUn","[0]",ptmin+0.01,ptmax);
    GUn = new TF1("GUn","[0]",ptmin+0.01,ptmax);
    HUn = new TF1("HUn","[0]",ptmin+0.01,ptmax);
 }
*/
 if( CutandStrength == "CSVv2T" || CutandStrength == "cMVAv2T")
 {
    FUn = new TF1("FUn","[0] + [1]/(x*x)",ptmin+0.01,ptmax);
    GUn = new TF1("GUn","[0] + [1]/(x*x)",ptmin+0.01,ptmax);
    HUn = new TF1("HUn","[0] + [1]/(x*x)",ptmin+0.01,ptmax);
 }
 if( CutandStrength == "JPT" )
 {
    FUn = new TF1("FUn","[0]",ptmin+0.01,ptmax);
    GUn = new TF1("GUn","[0]",ptmin+0.01,ptmax);
    HUn = new TF1("HUn","[0]",ptmin+0.01,ptmax);
 }
      FUn->SetLineColor(1);
      FUn->SetLineStyle(1);
      FUn->SetLineWidth((Width_t)0.5);
      GUn->SetLineColor(1);
      GUn->SetLineStyle(2);
      HUn->SetLineColor(1);
      HUn->SetLineStyle(2);

 TH1F* MisMax= (TH1F*) h000->Clone("MisMax"); MisMax->Reset();
 TH1F* MisMin= (TH1F*) h000->Clone("MisMin"); MisMin->Reset();

 TH1F* SFmax= (TH1F*) h000->Clone("SFmax"); SFmax->Reset();
 TH1F* SFmin= (TH1F*) h000->Clone("SFmin"); SFmin->Reset();

gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(42);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.05);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
gStyle->SetOptStat(stati);
gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);

if (fit) {
gStyle->SetOptFit(111);
gStyle->SetStatW(0.5);
gStyle->SetStatH(0.2);
} else {
gStyle->SetOptFit(0);
gStyle->SetStatW(0.4);
gStyle->SetStatH(0.3);
}
   cout << "Now, start with fits..." << endl;
   //getchar();

 TCanvas *cWork = new TCanvas("cWork", "plots",200,10,700,750);
 cWork->SetFillColor(10);
 cWork->SetFillStyle(4000);
 cWork->SetBorderSize(2);

 TFitResultPtr cdMFitptr = cDMistag->Fit("Fun", "rveeNS");
 TFitResult *tfr = cdMFitptr.Get();
 std::cout << "FitProb for higher chsquare : " << Fun->GetProb() << endl;
 std::cout << "Fit Status: " << tfr->Status() << endl;
 std::cout << "Fit CovMatrixStatus: " << tfr->CovMatrixStatus() << endl;
 cDMistag->Draw();
 Fun->Draw("same");
 cWork->Modified();
 cWork->Update();
 if( CheckFitMistag )
 {
   std::cout << "Fit for cDMistag. Quit if fit is not good enough" << endl;
   getchar();
   if( interactive ) cWork->WaitPrimitive();
 }

 TH1F* Hmistag= (TH1F*) h000->Clone("Hmistag"); Hmistag->Reset();
 Hmistag->Sumw2(); 
  for (int i=0; i<h000->GetNbinsX()+1; i++) {
	y = cDMistag->GetBinContent(i);
	ey = cDMistag->GetBinError(i);
	float syst = y * totMis->GetBinContent(i);
	float statMC = y * StatMis->GetBinContent(i);
	syst = syst*syst - statMC*statMC;
	if (syst > 0.) syst = TMath::Sqrt(syst);
        else
        {
          std::cout << "syst is negative for i = " << i << " syst = " << syst << " statMC = " << statMC << endl;
          syst = 0;
        }
	BinMin[i] = h000->GetBinLowEdge(i); //HistPtMin+(i-1)*binw;
	BinMax[i] = h000->GetBinLowEdge(i+1); //HistPtMin+i*binw;

// rel error on Fun is equal to integralError/integral
// abs error on Fun is equal to integralError/integral * Fun (x) ~ IntegralError/integral * integral/binwidth = IntegralError/Binwidth
        float stat = Fun->IntegralError(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
        float eror = TMath::Sqrt( stat*stat + syst*syst );
        std::cout << "eror = " << eror << " stat = " << stat << " syst = " << syst << " statMC = " << statMC << endl;

        float yfit = Fun->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);

	Hmistag->SetBinContent(i,yfit);
	Hmistag->SetBinError(i,0.);
	MisMax->SetBinContent(i,y + eror);
	MisMax->SetBinError(i,ey);
	MisMin->SetBinContent(i,y - eror);
	MisMin->SetBinError(i,ey);
	if ( y > 0. ) {
	  MisMax->SetBinError(i,ey * (1. + eror/y ));
	  MisMin->SetBinError(i,ey * (1. - eror/y ));
	}
	Leff[i] = Fun->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
       }
       Leff[h000->GetNbinsX()] = cDMistag->GetBinContent(h000->GetNbinsX());
       Leff[0] = cDMistag->GetBinContent(0);


       TFitResultPtr MismaxFitptr = MisMax->Fit("Gun", "rveeS");
       TFitResult *MisMaxtfr = MismaxFitptr.Get();
std::cout << "MisMax FitProb for higher chsquare : " << Gun->GetProb() << endl;
std::cout << "Fit Status: " << MisMaxtfr->Status() << endl;
std::cout << "Fit CovMatrixStatus: " << MisMaxtfr->CovMatrixStatus() << endl;

       cWork->Update();
if( CheckFitMistag )
{
std::cout << "Fit for cDMistag+err. Quit if fit is not good enough" << endl;
getchar();
       if( interactive ) cWork->WaitPrimitive();
}
       for (int i=0; i<h000->GetNbinsX()+1; i++) {
	LeffMax[i] = Gun->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
       }


       TFitResultPtr MisminFitptr = MisMin->Fit("Hun", "rveeS");
       TFitResult *MisMintfr = MismaxFitptr.Get();
std::cout << "MisMin FitProb for higher chsquare : " << Hun->GetProb() << endl;
std::cout << "Fit Status: " << MisMintfr->Status() << endl;
std::cout << "Fit CovMatrixStatus: " << MisMintfr->CovMatrixStatus() << endl;

       cWork->Update();
if( CheckFitMistag )
{
std::cout << "Fit for cDMistag-err. Quit if fit is not good enough" << endl;
getchar();
       if( interactive ) cWork->WaitPrimitive();
}
       for (int i=0; i<h000->GetNbinsX()+1; i++) {
	LeffMin[i] = Hun->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
	LeffErr[i] = (LeffMax[i] - LeffMin[i]) / 2.;
       }
       Leff[h000->GetNbinsX()] = cDMistag->GetBinContent(h000->GetNbinsX());
       LeffErr[h000->GetNbinsX()] = LeffErr[h000->GetNbinsX()-1] * Leff[h000->GetNbinsX()] / Leff[h000->GetNbinsX()-1];


       TFitResultPtr cSFFitptr = cSF->Fit("FUn", "rveeNS", "",ptmin+0.01);
       TFitResult *cSFtfr = cSFFitptr.Get();
std::cout << "cSF FitProb for higher chsquare : " << FUn->GetProb() << endl;
std::cout << "Fit Status: " << cSFtfr->Status() << endl;
std::cout << "Fit CovMatrixStatus: " << cSFtfr->CovMatrixStatus() << endl;

       cWork->Update();
if( CheckFitSF )
{
std::cout << "Fit for cSF. Quit if fit is not good enough" << endl;
getchar();
       if( interactive ) cWork->WaitPrimitive();
}

 TH1F* Hsf= (TH1F*) h000->Clone("Hsf"); Hsf->Reset();
       Hsf->Sumw2(); 
       for (int i=0; i<h000->GetNbinsX()+1; i++) {
	y = cSF->GetBinContent(i);
	ey = cSF->GetBinError(i);
	float syst = y * totSF->GetBinContent(i);
	float statMC = y * StatSF->GetBinContent(i);
	syst = syst*syst - statMC*statMC;
	if (syst > 0.) syst = TMath::Sqrt(syst);
        else if( i != 0 ){ std::cout << "What does it mean, for i= " << i << " syst = " << syst << " statMC= " << statMC << endl; if( interactive ) getchar(); }
        float stat = FUn->IntegralError(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
        float eror = TMath::Sqrt( stat*stat + syst*syst );
        float yfit = FUn->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
        std::cout << "yfit = " << yfit << " for i = " << i << " binw = " << (BinMax[i]-BinMin[i]) << " BinMin = " << BinMin[i] << " BinMax = " << BinMax[i] << endl;
	Hsf->SetBinContent(i,yfit);
	Hsf->SetBinError(i,0.);
	SFmax->SetBinContent(i,yfit + eror);
	SFmax->SetBinError(i,ey);
	SFmin->SetBinContent(i,yfit - eror);
	SFmin->SetBinError(i,ey);
	Lsf[i] = FUn->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
       }

//

       TFitResultPtr SFmaxFitptr = SFmax->Fit("GUn", "rveeNS", "",ptmin+0.01);
       TFitResult *SFmaxtfr = SFmaxFitptr.Get();
std::cout << "SFmax FitProb for higher chsquare : " << GUn->GetProb() << endl;
std::cout << "Fit Status: " << SFmaxtfr->Status() << endl;
std::cout << "Fit CovMatrixStatus: " << SFmaxtfr->CovMatrixStatus() << endl;

       cWork->Update();
if( CheckFitSF )
{
std::cout << "Fit for cSF+err. Quit if fit is not good enough" << endl;
getchar();
       if( interactive ) cWork->WaitPrimitive();
}


       TFitResultPtr SFminFitptr = SFmin->Fit("HUn", "rveeNS", "",ptmin+0.01);
       TFitResult *SFmintfr = SFminFitptr.Get();
std::cout << "SFmin FitProb for higher chsquare : " << HUn->GetProb() << endl;
std::cout << "Fit Status: " << SFmintfr->Status() << endl;
std::cout << "Fit CovMatrixStatus: " << SFmintfr->CovMatrixStatus() << endl;

       cWork->Update();
if( CheckFitSF )
{
std::cout << "Fit for cSF-err. Quit if fit is not good enough" << endl;
getchar();
//
       if( interactive ) cWork->WaitPrimitive();
}
 
 if( nPlot == 1 )
 {
   c1 = new TCanvas("c1", "plots",200,10,700,250);
 } else if ( nPlot == 2 ) {
   c1 = new TCanvas("c1", "plots",200,44,700,500);
 } else {
   c1 = new TCanvas("c1", "plots",200,10,700,750);
 }
 gStyle->SetOptStat(0);
 c1->Range(0,0,1,1);
 c1->SetFillColor(10);
 c1->SetFillStyle(4000);
 c1->SetBorderMode(0);
 c1->SetBorderSize(2);
 c1->SetFrameFillColor(0);
 c1->SetFrameBorderMode(0);

 c1->SetFillColor(10);
 c1->SetFillStyle(4000);
 c1->SetBorderSize(2);

 if( nPlot == 1 )
 {
   c1->Divide(1,1,0.01,0.002);
   c1->cd(1)->SetBottomMargin(0.15);
   c1->cd(1)->SetTopMargin(0.10);
 } else if ( nPlot == 2 ) {
   c1->Divide(1,2,0.01,0.002);
   c1->cd(1)->SetBottomMargin(0.15);
   c1->cd(1)->SetTopMargin(0.10);
   c1->cd(2)->SetBottomMargin(0.15);
   c1->cd(2)->SetTopMargin(0.10);
 } else {
   c1->Divide(1,3,0.01,0.002);
   c1->cd(1)->SetBottomMargin(0.15);
   c1->cd(1)->SetTopMargin(0.10);
   c1->cd(2)->SetBottomMargin(0.15);
   c1->cd(2)->SetTopMargin(0.10);
   c1->cd(3)->SetBottomMargin(0.15);
   c1->cd(3)->SetTopMargin(0.10);
 }

// *****************************************************************************

       totMis->SetLineColor(1);
       totMis->SetMarkerStyle(20);
       totMis->SetMarkerColor(1);
       totMis->SetMarkerSize(1.);
       totMis->SetTickLength(0.03, "YZ");
       totMis->SetTickLength(-0.03,"X");
       totMis->SetLabelOffset(0.023,"X");
       totMis->SetLabelOffset(0.007,"Y");
       totMis->SetLabelSize(0.055, "XYZ");
       totMis->SetLabelFont(42, "XYZ"); 
       totMis->SetTitleFont(42, "XYZ");
       totMis->SetTitleSize(0.07, "XYZ"); 
       totMis->SetTitleOffset(1.0,"X"); 
       totMis->SetTitleOffset(0.6,"Y");
       totMis->GetXaxis()->SetTitle(xtitle);
       totMis->GetXaxis()->SetTitleColor(1);
       totMis->GetXaxis()->SetNdivisions(508);
       totMis->GetYaxis()->SetTitle("relative syst. on Misidentification probability");
       totMis->GetYaxis()->SetTitleColor(1);
       totMis->GetYaxis()->SetNdivisions(509);
       totMis->SetMinimum(0.); 
       if ( CutStrength == "L" ) totMis->SetMaximum(0.24); // loose
       else if ( CutStrength == "M" ) totMis->SetMaximum(0.30); // medium
       else if ( CutStrength == "T" ) totMis->SetMaximum(0.50); // tight
       
       bc->SetLineColor(kRed);
       bc->SetMarkerStyle(21);
       bc->SetMarkerColor(kRed);
       bc->SetMarkerSize(1.);

       gg->SetLineColor(kGreen);
       gg->SetMarkerStyle(22);
       gg->SetMarkerColor(kGreen);
       gg->SetMarkerSize(1.);

       v0->SetLineColor(kBlue);
       v0->SetMarkerStyle(23);
       v0->SetMarkerColor(kBlue);
       v0->SetMarkerSize(1.);

       fk->SetLineColor(1);
       fk->SetMarkerStyle(25);
       fk->SetMarkerColor(1);
       fk->SetMarkerSize(1.);

       sf->SetLineColor(1);
       sf->SetMarkerStyle(24);
       sf->SetMarkerColor(1);
       sf->SetMarkerSize(1.);
 
       puMis->SetLineColor(kOrange);
       puMis->SetMarkerStyle(29);
       puMis->SetMarkerColor(kOrange);
       puMis->SetMarkerSize(1.2);


  std::cout << " You asked an integration from " << minPtIntegration << " to " << maxPtIntegration << endl;
  std::cout << " Binning will actually give you integration from " << bc->GetBinLowEdge(iminPtInt) << " to " << bc->GetBinLowEdge(imaxPtInt) << endl;
 minPtIntegration = bc->GetBinLowEdge(iminPtInt);
 maxPtIntegration = bc->GetBinLowEdge(imaxPtInt);

  std::cout << " You asked a second integration from " << minPtIntegration2 << " to " << maxPtIntegration2 << endl;
  std::cout << " Binning will actually give you integration from " << bc->GetBinLowEdge(iminPtInt2) << " to " << bc->GetBinLowEdge(imaxPtInt2) << endl;
 minPtIntegration2 = bc->GetBinLowEdge(iminPtInt2);
 maxPtIntegration2 = bc->GetBinLowEdge(imaxPtInt2);

     aout_txt << "All values for Pt between " << minPtIntegration << " and " << maxPtIntegration << endl;
     aout_txt << "-----------------" << endl;
     aout_txt << "Syst errors on SF" << endl;
     aout_txt << "-----------------" << endl;
     if( nPlot > 1 )
     {
       c1->cd(2);
     } else {
       c1->cd(1);
     }
       std::cout << "Print i, LowEdge, BinContent for totSF histo " << endl;
       for(int i = 0; i <= totSF->GetNbinsX()+1; ++i )
       {
         std::cout << i << "\t" << totSF->GetBinLowEdge(i) << "\t" << totSF->GetBinContent(i) << endl;
       }
       std::cout << "Print done, type enter to continue " << endl;

       totSF->Draw("E"); 
       totSF->SetLineColor(1);
       totSF->SetMarkerStyle(20);
       totSF->SetMarkerColor(1);
       totSF->SetMarkerSize(1.);
       totSF->SetTickLength(0.03, "YZ");
       totSF->SetTickLength(-0.03,"X");
       totSF->SetLabelOffset(0.023,"X");
       totSF->SetLabelSize(0.055, "XYZ");
       totSF->SetLabelOffset(0.007,"Y");
       totSF->SetLabelSize(0.055, "XYZ");
       totSF->SetLabelFont(42, "XYZ"); 
       totSF->SetTitleFont(42, "XYZ");
       totSF->SetTitleSize(0.07, "XYZ"); 
       totSF->SetTitleOffset(1.0,"X"); 
       totSF->SetTitleOffset(0.6,"Y");
       totSF->GetXaxis()->SetRangeUser(ptminDraw+0.01,ptmaxDraw-0.01);
       totSF->GetXaxis()->SetTitle(xtitle);
       totSF->GetXaxis()->SetTitleColor(1);
       totSF->GetXaxis()->SetNdivisions(508);
       totSF->GetYaxis()->SetTitle("relative syst. on mistag SF");
       totSF->GetYaxis()->SetTitleColor(1);
       totSF->GetYaxis()->SetNdivisions(509);
       totSF->SetMinimum(0.); 
//NOTE set range second plots: relative syst. on mistag SF
       if( CutStrength == "L" ) {
         totSF->SetMaximum(0.13);
         if (CutName == "CvsL" || CutName == "CvsB") totSF->SetMaximum(0.045) ;
       }
       if( CutStrength == "M" ) {
         totSF->SetMaximum(0.3);
         if (CutName == "CvsL" || CutName == "CvsB") totSF->SetMaximum(0.1) ;
       }
       if( CutStrength == "T" ) {
         totSF->SetMaximum(0.4);
         if (CutName == "CvsL" || CutName == "CvsB") totSF->SetMaximum(0.2) ;
       }
       
       bc->Draw("Esame"); 
       aout_txt << mygeterr(false, true, bc, iminPtInt, imaxPtInt ) << "\t:b+c " << endl;
       gg->Draw("Esame"); 
       aout_txt << mygeterr(false, true, gg, iminPtInt, imaxPtInt ) << "\t:gluon" << endl;
       v0->Draw("Esame"); 
       aout_txt << mygeterr(false, true, v0, iminPtInt, imaxPtInt ) << "\t:V0" << endl;
       fk->Draw("Esame"); 
       aout_txt << mygeterr(false, true, fk, iminPtInt, imaxPtInt ) << "\t:fake" << endl;
       sf->Draw("Esame"); 
       aout_txt << mygeterr(false, true, sf, iminPtInt, imaxPtInt ) << "\t:sign flip" << endl;

       puSF->Draw("Esame"); 
       puSF->SetLineColor(kOrange);
       puSF->SetMarkerStyle(29);
       puSF->SetMarkerColor(kOrange);
       puSF->SetMarkerSize(1.2);
       aout_txt << mygeterr(true, true, puSF, iminPtInt, imaxPtInt ) << "\t:MC stat" << endl;

       totSF->Draw("Esame"); 
       aout_txt << mygeterr(false, true, totSF, iminPtInt, imaxPtInt ) << "\t:All + sample" << endl;

  TLegend* leg=NULL;
/*
  if ( CutStrength == "L" )
       leg = new TLegend(0.66,0.38,0.80,0.78); // loose
  else leg = new TLegend(0.66,0.54,0.80,0.94);
*/
    leg = new TLegend(0.66,0.49,0.80,0.89); // loose
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(62);
    leg->AddEntry(totSF,"all+sample","P");
    leg->AddEntry(puSF,"MC stat","P");
    leg->AddEntry(sf," sign flip","P");
    leg->AddEntry(fk,"fake tracks","P");
    leg->AddEntry(v0,"V^{0}+2nd.int.","P");
    leg->AddEntry(bc," b+c jets","P");
    leg->AddEntry(gg," gluon jets","P");
    leg->Draw();

         totSF->Write();
         puSF->Write();

std::cout << endl;
std::cout << htitle << " syst" 
<< "  bc " << (bc->GetBinContent(4)+bc->GetBinContent(5)+bc->GetBinContent(6))/0.03
<< "  g  " << (gg->GetBinContent(4)+gg->GetBinContent(5)+gg->GetBinContent(6))/0.03
<< "  v0 " << (v0->GetBinContent(4)+v0->GetBinContent(5)+v0->GetBinContent(6))/0.03
<< "  fk " << (fk->GetBinContent(4)+fk->GetBinContent(5)+fk->GetBinContent(6))/0.03
<< "  sf " << (sf->GetBinContent(4)+sf->GetBinContent(5)+sf->GetBinContent(6))/0.03
<< "  pu " << (puSF->GetBinContent(4)+puSF->GetBinContent(5)+puSF->GetBinContent(6))/0.03
<< " Tot " << (totSF->GetBinContent(4)+totSF->GetBinContent(5)+totSF->GetBinContent(6))/0.03 << endl;
std::cout << endl;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Mistag in Data = all neg data * Rlight
//NOTE set range first plot
     if( nPlot > 1 )
     {
       c1->cd(1);
     } else {
       c1->cd(1);
       c1->cd(1)->Clear();
     }
       retMistag = *cDMistag;
       retMistag.SetDirectory(0);

       if( CutStrength == "L" ){ 
         cDMistag->SetMinimum(0.); cDMistag->SetMaximum(0.6); 
         if (CutName == "CvsL" || CutName == "CvsB") {
           cDMistag->SetMinimum(0.8); cDMistag->SetMaximum(1.4); 
         }
       }
       if( CutStrength == "M" ){ 
         cDMistag->SetMinimum(0.); cDMistag->SetMaximum(0.1); 
         if (CutName == "CvsL" || CutName == "CvsB") {
           cDMistag->SetMinimum(0.1); cDMistag->SetMaximum(0.4); 
         }
       }
       if( CutStrength == "T" ){ 
         cDMistag->SetMinimum(0.); cDMistag->SetMaximum(0.012); 
         if (CutName == "CvsL" || CutName == "CvsB") {
           cDMistag->SetMinimum(0.0); cDMistag->SetMaximum(0.05); 
         }
       }
       cDMistag->Draw("E"); 
       cDMistag->SetLineColor(2);
       cDMistag->SetMarkerStyle(20);
       cDMistag->SetMarkerColor(2);
       cDMistag->SetMarkerSize(1.);
       cMCMistag->Draw("Esame");
       cMCMistag->SetLineColor(2);
       cMCMistag->SetMarkerStyle(24);
       cMCMistag->SetMarkerColor(2);
       cMCMistag->SetMarkerSize(1.);
       cDMistag->SetTickLength(0.03, "YZ");
       cDMistag->SetTickLength(-0.03,"X");
       cDMistag->SetLabelOffset(0.023,"X");
       cDMistag->SetLabelOffset(0.007,"Y");
       cDMistag->SetLabelSize(0.055, "XYZ");
       cDMistag->SetLabelFont(42, "XYZ"); 
       cDMistag->SetTitleFont(42, "XYZ");
       cDMistag->SetTitleSize(0.07, "XYZ"); 
       cDMistag->SetTitleOffset(1.0,"X"); 
       cDMistag->SetTitleOffset(0.6,"Y");
       cDMistag->GetXaxis()->SetTitle(xtitle);
       cDMistag->GetXaxis()->SetTitleColor(1);
       cDMistag->GetXaxis()->SetRangeUser(ptminDraw+0.01,ptmaxDraw-0.01);
       cDMistag->GetXaxis()->SetNdivisions(508);
       cDMistag->GetYaxis()->SetTitle("Misidentification probability");
       cDMistag->GetYaxis()->SetTitleColor(1);
       cDMistag->GetYaxis()->SetNdivisions(509);
       cDMistag->SetMinimum(0.0); 
/*
       if ( CutandStrength == "TCHEL" || CutandStrength == "TCHPL" ) cDMistag->SetMaximum(0.60); // loose tc
       if ( CutandStrength == "JPL" ) cDMistag->SetMaximum(0.25); // loose jp
       if ( CutandStrength == "JBPL" ) cDMistag->SetMaximum(0.40); // loose jbp csv
       if ( CutandStrength == "CSVL" ) cDMistag->SetMaximum(0.30); // loose jbp csv
       if ( CutandStrength == "TCHEM" ) cDMistag->SetMaximum(0.17); // medium tche
       if ( CutandStrength == "TCHPM" ) cDMistag->SetMaximum(0.25); // medium tchp
       if ( CutandStrength == "SSVHEM" || CutandStrength == "SSVHPM" ) cDMistag->SetMaximum(0.06); // medium ssv
       if ( CutandStrength == "JPM" ) cDMistag->SetMaximum(0.05); // medium jp
       if ( CutandStrength == "JBPM" ) cDMistag->SetMaximum(0.08); // medium jbp
       if ( CutandStrength == "CSVM" ) cDMistag->SetMaximum(0.05); // medium csv
       if ( CutandStrength == "TCHET" || CutandStrength == "TCHPT" ) cDMistag->SetMaximum(0.040); // tight tc
       if ( CutandStrength == "SSVHET" || CutandStrength == "SSVHPT" ) cDMistag->SetMaximum(0.014); // tight ssv
       if ( CutandStrength == "JPT" ) cDMistag->SetMaximum(0.010); // tight jp
       if ( CutandStrength == "JBPT" ) cDMistag->SetMaximum(0.018); // tight jbp
       if ( CutandStrength == "CSVT" ) cDMistag->SetMaximum(0.010); // tight csv
*/
       cDMistag->Draw("Esame"); 
       cMCMistag->Draw("Esame");
       cDMistag->Draw("Esame"); 


    if( CutStrength == "L" ) leg = new TLegend(0.76, 0.175, 0.97, 0.455);
    else leg = new TLegend(0.32, 0.61, 0.53, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(62);
    leg->SetHeader(htitle);
    leg->AddEntry(cDMistag,"Data","P");
    leg->AddEntry(cMCMistag,"MC","P");
    leg->Draw();

  CMS_lumi((TPad*)c1->cd(1),Lumi,10,true);

         Fun->Write();
         Gun->Write();
         Hun->Write();

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Data/MC Scale Factor

     if( nPlot == 2 )
     {
       c1->cd(2);
       c1->cd(2)->Clear();
     } else if ( nPlot == 1 ){
       c1->cd(1);
       c1->cd(1)->Clear();
     } else {
     c1->cd(3);
     }
       retSF = *cSF;
       retSF.SetDirectory(0);

       cSF->Draw("E"); 
       cSF->SetLineColor(2);
       cSF->SetMarkerStyle(20);
       cSF->SetMarkerColor(2);
       cSF->SetMarkerSize(1.);
       cSF->SetTickLength(0.03, "YZ");
       cSF->SetTickLength(-0.03,"X");
       cSF->SetLabelOffset(0.023,"X");
       cSF->SetLabelOffset(0.007,"Y");
       cSF->SetLabelSize(0.055, "XYZ");
       cSF->SetLabelFont(42, "XYZ"); 
       cSF->SetTitleFont(42, "XYZ");
       cSF->SetTitleSize(0.07, "XYZ"); 
       cSF->SetTitleOffset(1.0,"X"); 
       cSF->SetTitleOffset(0.6,"Y");
       cSF->GetXaxis()->SetTitle(xtitle);
       cSF->GetXaxis()->SetTitleColor(1);
       cSF->GetXaxis()->SetNdivisions(508);
       cSF->GetYaxis()->SetTitle("Data / MC misid. SF");
       cSF->GetYaxis()->SetTitleColor(1);
       cSF->GetXaxis()->SetRangeUser(ptminDraw+0.01,ptmaxDraw-0.01);
//NOTE: change range here third plots
       cSF->GetYaxis()->SetRangeUser(0.1,1.65);
       cSF->GetYaxis()->SetNdivisions(508);

       //if( CutStrength == "L" ){ 
         //cSF->SetMinimum(0.6); cSF->SetMaximum(1.4);
       //}
       //if( CutStrength == "M" ){ 
         //cSF->SetMinimum(0.4); cSF->SetMaximum(1.6); 
       //}
       //if( CutStrength == "T" ){ 
         //cSF->SetMinimum(0.4); cSF->SetMaximum(1.6); 
       //}

       FUn->Draw("same"); 

     aout_txt << "-----------------" << endl;
     aout_txt << " Stat error on SF " << endl;
     aout_txt << "-----------------" << endl;
     aout_txt << mygeterr(true, false, cSF, iminPtInt, imaxPtInt ) << "\t: Stat from histo" << endl;

     aout_txt << "-----------------" << endl;
     aout_txt << " unbiased SF " << endl;
     aout_txt << "-----------------" << endl;
     aout_txt << FUn->Integral(minPtIntegration,maxPtIntegration)/(maxPtIntegration-minPtIntegration) << "\t: Mean of interp fct from " << minPtIntegration << " to " << maxPtIntegration << endl;
     aout_txt << "-----------------" << endl;
     aout_txt << "Mistag rate " << endl;
     aout_txt << "-----------------" << endl;
     aout_txt << mygetmean(cDMistag, iminPtInt, imaxPtInt ) << "\t:Mistag rate" << endl;
     aout_txt << mygeterr(true, false, cDMistag, iminPtInt, imaxPtInt ) << "\t:Stat error on Misidentification probability" << endl;

if( ManualErrors )
{
 char buffer[50];
 Double_t manualMax = (1.+ManualRelError);
 sprintf(buffer,"%f",manualMax);
 TString TManualGUn = "";
 TManualGUn += buffer;
 TManualGUn += "* ( " + FUn->GetExpFormula("p") + " )";
 cout << "TManualGUn = " << TManualGUn << endl;
 GUn = new TF1("manualGUn",TManualGUn,ptmin+0.01,ptmax);
 GUn->SetLineColor(1);
 GUn->SetLineStyle(2);

 Double_t manualMin = (1.-ManualRelError);
 sprintf(buffer,"%f",manualMin);
 TString TManualHUn = "";
 TManualHUn += buffer;
 TManualHUn += "* ( " + FUn->GetExpFormula("p") + " )";
 cout << "TManualHUn = " << TManualHUn << endl;
 HUn = new TF1("manualHUn",TManualHUn,ptmin+0.01,ptmax);
 HUn->SetLineColor(1);
 HUn->SetLineStyle(2);
}
       GUn->Draw("same"); 

       for (int i=0; i<h000->GetNbinsX()+1; i++) {
	LsfMax[i] = GUn->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
       }

       HUn->Draw("same"); 
       for (int i=0; i<h000->GetNbinsX()+1; i++) {
	LsfMin[i] = HUn->Integral(BinMin[i],BinMax[i]) / (BinMax[i]-BinMin[i]);
	LsfErr[i] = (LsfMax[i] - LsfMin[i]) / 2.;
       }
       BinMax[h000->GetNbinsX()] = 999.;

   if( CutStrength == "L" ) leg = new TLegend(0.197, 0.171, 0.447, 0.451);
   if( CutStrength == "M" ) leg = new TLegend(0.197, 0.171, 0.447, 0.451);
   if( CutStrength == "T" ) leg = new TLegend(0.197+0.343, 0.171+0.539, 0.447+0.343, 0.451+0.539);
    leg->SetBorderSize(0);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(62);
    leg->SetHeader(htitle);
    leg->AddEntry(cSF,"Data / MC misid. SF","P");
    if( ManualErrors ) leg->AddEntry(GUn,"#pm stat #oplus syst = "+tManualRelError,"L");
    else leg->AddEntry(GUn,"#pm stat #oplus syst","L");
    leg->Draw();
    
         FUn->Write();
         GUn->Write();
         HUn->Write();

  for ( int iJetCut = 0;  iJetCut < NJetCuts; ++iJetCut)
  {
    MCNegtag_Jets[iJetCut]->Write();
    cMCNegtag_Jets[iJetCut]->Write();
    MCMistag_Jets[iJetCut]->Write();
    cMCMistag_Jets[iJetCut]->Write();
    MCRlightWeight_Jets[iJetCut]->Write();
    cMCRlightWeight_Jets[iJetCut]->Write();
   

    SF_Jets[iJetCut]->Write();
    cSF_Jets[iJetCut]->Write();
    DNegtag_Jets[iJetCut]->Write();
    cDNegtag_Jets[iJetCut]->Write();
    DMistag_Jets[iJetCut]->Write();
    cDMistag_Jets[iJetCut]->Write();
  }
  cMCMistag->Write();
  cMCNegtag->Write();
  cDMistag->Write();
  cDNegtag->Write();
  cSF->Write();
  SFmax->Write();
  SFmin->Write();

  fHis.Close();

       float Mis100 = Fun->Integral(50,80)/30.;
       float Mis100Err = (Gun->Integral(50,80)/30. - Hun->Integral(50,80)/30.) / 2;

       float SF100 = FUn->Integral(50,80)/30.;
       float SF100Err = (GUn->Integral(50,80)/30. - HUn->Integral(50,80)/30.) / 2;

std::cout << endl;

float stat1 = 1. / TMath::Sqrt(1./(cDMistag->GetBinError(4)*cDMistag->GetBinError(4))
                              +1./(cDMistag->GetBinError(5)*cDMistag->GetBinError(5))
                              +1./(cDMistag->GetBinError(6)*cDMistag->GetBinError(6)));
float stat2 = 1. / TMath::Sqrt(1./(cSF->GetBinError(4)*cSF->GetBinError(4))
                              +1./(cSF->GetBinError(5)*cSF->GetBinError(5))
                              +1./(cSF->GetBinError(6)*cSF->GetBinError(6)));
float syst1 = TMath::Sqrt(Mis100Err*Mis100Err - stat1*stat1);
float syst2 = TMath::Sqrt(SF100Err*SF100Err - stat2*stat2);

std::cout << htitle << " at 65 GeV: cDMistag = " 
               << Mis100 << " +_ " << stat1 << " +_ " << syst1
 << "   cSF = " << SF100  << " +_ " << stat2 << " +_ " << syst2 << endl;
std::cout << endl;
std::cout << htitle 
<< " syst bc " << (bc->GetBinContent(4)+bc->GetBinContent(5)+bc->GetBinContent(6))/3
<< "  g  " << (gg->GetBinContent(4)+gg->GetBinContent(5)+gg->GetBinContent(6))/3
<< "  v0 " << (v0->GetBinContent(4)+v0->GetBinContent(5)+v0->GetBinContent(6))/3
<< "  fk " << (fk->GetBinContent(4)+fk->GetBinContent(5)+fk->GetBinContent(6))/3
<< "  sf " << (sf->GetBinContent(4)+sf->GetBinContent(5)+sf->GetBinContent(6))/3
<< "  pu " << (puSF->GetBinContent(4)+puSF->GetBinContent(5)+puSF->GetBinContent(6))/3
<< " Tot " << (totSF->GetBinContent(4)+totSF->GetBinContent(5)+totSF->GetBinContent(6))/3 << endl;
std::cout << endl;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 float systL0 = 0., systL1 = 1., systL2 = 0.;
 for (int i=1; i<23; i++) {
  if ( totSF->GetBinContent(i) < systL1 ) systL1 = totSF->GetBinContent(i); 
  if ( totSF->GetBinContent(i) > systL2 ) systL2 = totSF->GetBinContent(i); 
 }
 int nsystL = 0;
 for (int i=1; i<23; i++) {
 if ( totSF->GetBinContent(i) >= systL2 ) continue;
  if ( totSF->GetBinContent(i) > systL0 ) systL0 = totSF->GetBinContent(i); 
 }
 float statL0 = FUn->GetParError(0);
 float sferrL = TMath::Sqrt(statL0*statL0 + systL0*systL0);

std::cout << htitle << " 20<pt<240 cSF = " << FUn->GetParameter(0)  << " +_ " << (GUn->GetParameter(0) - HUn->GetParameter(0)) / 2. << endl;
std::cout << htitle << " 20<pt<240 cSF = " << FUn->GetParameter(0)  << " +_ " << sferrL << endl;
std::cout << nsystL << " " << statL0  << " " << systL0 << " " << systL1 << " " << systL2 << endl;
std::cout << endl;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// Write the ascii file
    
 double mis[10], miserr[10], miser1[10], miser2[10], sfl[10], sflerr[10];
       
 for (int i=0; i<10; i++) {
  mis[i] = 0; miserr[i] = 0; sfl[i] = 0;  sflerr[i] = 0; 
  mis[i]    = Fun->GetParameter(i);
  miserr[i] = (Gun->GetParameter(i) - Hun->GetParameter(i)) / 2.;
  sfl[i]    = FUn->GetParameter(i);
  sflerr[i] = (GUn->GetParameter(i) - HUn->GetParameter(i)) / 2.;
  miser1[i] = Gun->GetParameter(i);
  miser2[i] = Hun->GetParameter(i);
 }
 miser1[0] = miser1[0]/2.;
 miser2[0] = miser2[0]/2.;
 

 ofstream out_txt(blabla_file, ios::out);

 out_txt << TaggerName << endl; // < ---the tagger name
 out_txt << wp << endl; // <---- the WP cut on the discriminator
 out_txt << "PerformancePayloadFromTFormula" << endl;
 out_txt << 4 << endl; // <---- # of results, lighteff and lightefferr
 out_txt << 2 << endl; // <---- # of variables to bin (eta and pt)
 out_txt << "1005 1006 1013 1014" << endl; // results mistag +_ err, SF +_ err
 out_txt << "2 5" << endl; // <----- binning is pt |eta|
if ( CutandStrength != "CSVL" ) {
 out_txt << mis[0] << " + " << mis[1] << "*x"  
		   << " + " << mis[2] << "*x*x"
		   << " + " << mis[3] << "*x*x*x" << endl;
 out_txt << miserr[0] << " + " << miserr[1] << "*x"  
		      << " + " << miserr[2] << "*x*x"
		      << " + " << miserr[3] << "*x*x*x" << endl;
}
if ( CutandStrength == "CSVL" ) {
 out_txt << mis[0] << "*(1+" << mis[1] << "*x"  
		   << "+" << mis[2] << "*x*x)/(1+" << mis[3] << "*x)" << endl;
 out_txt << miser1[0] << "*(1+" << miser1[1] << "*x"  
               << "+" << miser1[2] << "*x*x)/(1+" << miser1[3] << "*x)-"
         << miser2[0] << "*(1+" << miser2[1] << "*x"  
               << "+" << miser2[2] << "*x*x)/(1+" << miser2[3] << "*x)" << endl;
}
 out_txt << sfl[0] << " + " << sfl[1] << "*x"  
		   << " + " << sfl[2] << "*x*x"
		   << " + " << sfl[3] << "*x*x*x" << endl;
 out_txt << sflerr[0] << " + " << sflerr[1] << "*x"  
		      << " + " << sflerr[2] << "*x*x"
		      << " + " << sflerr[3] << "*x*x*x" << endl;
 out_txt << HistPtMin << " " << HistPtMax << " " << EtaMin << " " << EtaMax << endl;
 out_txt << endl;
 out_txt << endl;
       
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  c1->Update();
  std::cout << "Make last changes to plot before saving it as .png file" << endl;
  if( ImprovePlot ) c1->WaitPrimitive();
  std::cout << " passed ImprovePlot" << endl;
  c1->Update();
  std::cout << " Updated" << endl;
  c1->SaveAs(SavDir + "/png/"+oFilePrefix+tnPlot+"_"+CutandStrength+"_"+EtaRange+".png","");
  c1->SaveAs(SavDir + "/pdf/"+oFilePrefix+tnPlot+"_"+CutandStrength+"_"+EtaRange+".pdf","");
  c1->SaveAs(SavDir + "/C/"+oFilePrefix+tnPlot+"_"+CutandStrength+"_"+EtaRange+".C","");
  c1->SaveAs(SavDir + "/root/"+oFilePrefix+tnPlot+"_"+CutandStrength+"_"+EtaRange+".root","");
  std::cout << " Saved" << endl;

  if( PrintDBSF )
  {
     ofstream sfout_txt(oFil+"_SF_Syst.out", ios::out ); // | ios::app );
     sfout_txt << setw(10);
     sfout_txt << CutandStrength << " & $ ";
     sfout_txt << setiosflags(ios::fixed);
     sfout_txt << setprecision(1);

     sfout_txt << 100*mygeterr(false, true, bc, iminPtInt, imaxPtInt ) << "\\\% $ & $ ";	     //"b+c"
     sfout_txt << 100*mygeterr(false, true, gg, iminPtInt, imaxPtInt ) << "\\\% $ & $ ";       //"gluon"
     sfout_txt << 100*mygeterr(false, true, v0, iminPtInt, imaxPtInt ) << "\\\% $ & $ ";       //"V0"
     sfout_txt << 100*mygeterr(false, true, fk, iminPtInt, imaxPtInt ) << "\\\% $ & $ ";       //"mis-meas."
     sfout_txt << 100*mygeterr(false, true, sf, iminPtInt, imaxPtInt ) << "\\\% $ & $ ";       //"sign flip"
     sfout_txt << 100*mygeterr(true, true, puSF, iminPtInt, imaxPtInt ) << "\\\% $ & $ ";     //"Pileup+MC Stat"
     sfout_txt << 100*ErrSysSF << "\\\% $ & $ ";                                          //"Sample"
     sfout_txt << 100*mygeterr(false, true, totSF, iminPtInt, imaxPtInt ) << "\\\% $ \\\\  ";    //"All+Sample"
     sfout_txt << endl;

     ofstream sfout2_txt(oFil+"_SF_Syst2.out", ios::out ); // | ios::app );
     sfout2_txt << setw(10);
     sfout2_txt << CutandStrength << " & $ ";
     sfout2_txt << setiosflags(ios::fixed);
     sfout2_txt << setprecision(1);

     sfout2_txt << 100*mygeterr(false, true, bc, iminPtInt2, imaxPtInt2 ) << "\\\% $ & $ ";	     //"b+c"
     sfout2_txt << 100*mygeterr(false, true, gg, iminPtInt2, imaxPtInt2 ) << "\\\% $ & $ ";       //"gluon"
     sfout2_txt << 100*mygeterr(false, true, v0, iminPtInt2, imaxPtInt2 ) << "\\\% $ & $ ";       //"V0"
     sfout2_txt << 100*mygeterr(false, true, fk, iminPtInt2, imaxPtInt2 ) << "\\\% $ & $ ";       //"mis-meas."
     sfout2_txt << 100*mygeterr(false, true, sf, iminPtInt2, imaxPtInt2 ) << "\\\% $ & $ ";       //"sign flip"
     sfout2_txt << 100*mygeterr(true, true, puSF, iminPtInt2, imaxPtInt2 ) << "\\\% $ & $ ";     //"Pileup+MC Stat"
     sfout2_txt << 100*ErrSysSF << "\\\% $ & $ ";                                          //"Sample"
     sfout2_txt << 100*mygeterr(false, true, totSF, iminPtInt2, imaxPtInt2 ) << "\\\% $ \\\\  ";    //"All+Sample"
     sfout2_txt << endl;
  }
  if( PrintDBSF_mistag )
  {
     ofstream bout_txt(oFil+"_SF_Mistag.out", ios::out ); // | ios::app );
     bout_txt << setw(10);
     bout_txt << CutandStrength << " & $";
     bout_txt << setiosflags(ios::fixed);
     if( CutStrength == "T" ) bout_txt << setprecision(5);
     else bout_txt << setprecision(4);
     bout_txt << mygetmean(cDMistag, iminPtInt, imaxPtInt ) << "\\pm  "; // Mistag rate
     bout_txt << mygetmean(cDMistag, iminPtInt, imaxPtInt )*mygeterr(true, false, cDMistag, iminPtInt, imaxPtInt ) << "$ & $";  // MistagRateStatError

     if( CutStrength == "L" ) bout_txt << setprecision(3);
     else bout_txt << setprecision(2);
     bout_txt << FUn->Integral(minPtIntegration,maxPtIntegration)/(maxPtIntegration-minPtIntegration) << "\\pm"; // Unbiased SF
     bout_txt << FUn->Integral(minPtIntegration,maxPtIntegration)/(maxPtIntegration-minPtIntegration)*mygeterr(true, false, cSF, iminPtInt, imaxPtInt ) << "\\pm ";       // SFStatErrorFromHisto
     bout_txt << FUn->Integral(minPtIntegration,maxPtIntegration)/(maxPtIntegration-minPtIntegration)*mygeterr(false, true, totSF, iminPtInt, imaxPtInt ) << "$ \\\\  ";    //SFSystErr "All+Sample"
     bout_txt << endl;

     ofstream bout2_txt(oFil+"_SF_Mistag2.out", ios::out ); // | ios::app );
     bout2_txt << setw(10);
     bout2_txt << CutandStrength << " & $";
     bout2_txt << setiosflags(ios::fixed);
     if( CutStrength == "T" ) bout2_txt << setprecision(5);
     else bout2_txt << setprecision(4);
     bout2_txt << mygetmean(cDMistag, iminPtInt2, imaxPtInt2 ) << "\\pm  "; // Mistag rate
     bout2_txt << mygetmean(cDMistag, iminPtInt2, imaxPtInt2 )*mygeterr(true, false, cDMistag, iminPtInt2, imaxPtInt2 ) << "$ & $";  // MistagRateStatError

     if( CutStrength == "L" ) bout2_txt << setprecision(3);
     else bout2_txt << setprecision(2);
     bout2_txt << FUn->Integral(minPtIntegration2,maxPtIntegration2)/(maxPtIntegration2-minPtIntegration2) << "\\pm"; // Unbiased SF
     bout2_txt << FUn->Integral(minPtIntegration2,maxPtIntegration2)/(maxPtIntegration2-minPtIntegration2)*mygeterr(true, false, cSF, iminPtInt2, imaxPtInt2 ) << "\\pm ";       // SFStatErrorFromHisto
     bout2_txt << FUn->Integral(minPtIntegration2,maxPtIntegration2)/(maxPtIntegration2-minPtIntegration2)*mygeterr(false, true, totSF, iminPtInt2, imaxPtInt2 ) << "$ \\\\  ";    //SFSystErr "All+Sample"
     bout2_txt << endl;
  }
  if ( PrintDBvalues )
  {
     float digit = 1e6;
        
     ofstream cout_txt(oFil+"_DBValues.out", ios::out ); // | ios::app);

     cout_txt << TaggerName << endl; // < ---the tagger name
     cout_txt << wp << endl; // <---- the WP cut on the discriminator
     cout_txt << "PerformancePayloadFromTable" << endl;
     cout_txt << 4 << endl; // <---- # of results, lighteff and lightefferr
     cout_txt << 2 << endl; // <---- # of variables to bin (eta and pt)
     cout_txt << "1005 1006 ";
     cout_txt << "1013 1014" << endl; // <----- results are light eff +_ error, light SF +_ error
     cout_txt << "5 2" << endl; // <----- binning is |eta| / pt

     for (int i=1; i<h000->GetNbinsX() && BinMin[i] < ptmax; i++) { // <----- eta and Et ranges, light eff +_ error, light SF +_ error
      if ( Lsf[i] >= 1 ) digit = 1e6;
      else digit = 1e5;
      if ((EtaMin == 0 || EtaMin == 1 || EtaMin == 2) && BinMin[i] < 100 && BinMax[i] < 100) {  
        cout_txt << EtaMin    << "     " << EtaMax     << "      " 
                << BinMin[i] << "      " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if ((EtaMin == 0 || EtaMin == 1 || EtaMin == 2) && BinMin[i] < 100) {  
        cout_txt << EtaMin    << "     " << EtaMax     << "      " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if (EtaMin == 0 || EtaMin == 1 || EtaMin == 2) {  
        cout_txt << EtaMin    << "     " << EtaMax     << "     " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if ((TMath::Abs(EtaMax-1.0)<0.01 || TMath::Abs(EtaMax-2.0)<0.01) && BinMin[i] < 100 && BinMax[i] < 100) {  
        cout_txt << EtaMin    << "   " << EtaMax     << "        " 
                << BinMin[i] << "      " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if ((TMath::Abs(EtaMax-1.0)<0.01 || TMath::Abs(EtaMax-2.0)<0.01) && BinMin[i] < 100) {  
       cout_txt << EtaMin    << "   " << EtaMax     << "        " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if (TMath::Abs(EtaMax-1.0)<0.01 || TMath::Abs(EtaMax-2.0)<0.01) {  
        cout_txt << EtaMin    << "   " << EtaMax     << "       " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if ((TMath::Abs(EtaMin-0.8)<0.01 || TMath::Abs(EtaMin-1.6)<0.01) && BinMin[i] < 100 && BinMax[i] < 100) {  
        cout_txt << EtaMin    << "   " << EtaMax     << "      " 
                << BinMin[i] << "      " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if ((TMath::Abs(EtaMin-0.8)<0.01 || TMath::Abs(EtaMin-1.6)<0.01) && BinMin[i] < 100) {  
       cout_txt << EtaMin    << "   " << EtaMax     << "      " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if (TMath::Abs(EtaMin-0.8)<0.01 || TMath::Abs(EtaMin-1.6)<0.01) {  
        cout_txt << EtaMin    << "   " << EtaMax     << "     " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if (BinMin[i] < 100 && BinMax[i] < 100) {  
        cout_txt << EtaMin    << "   " << EtaMax     << "      " 
                << BinMin[i] << "      " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else if (BinMin[i] < 100) {  
        cout_txt << EtaMin    << "   " << EtaMax     << "        " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
      else {  
        cout_txt << EtaMin    << "   " << EtaMax     << "       " 
                << BinMin[i] << "     " << BinMax[i]  << "     " 
         << int(Leff[i]*1e6)/1e6   << "     " << int(LeffErr[i]*1e6)/1e6 << "     " 
         << int(Lsf[i]*digit)/digit    << "     " << int(LsfErr[i]*1e5)/1e5  << endl;
      }
     }
  }
  if ( PrintDBfunction )
  {
     TString TSFUn = FUn->GetExpFormula("p");
     TString TSGUn = GUn->GetExpFormula("p");
     TString TSHUn = HUn->GetExpFormula("p");
     TString sEtaMin = Form("%1.1f",EtaMin);
     TString sEtaMax = Form("%1.1f",EtaMax);

     ofstream dout_txt(oFil+"_DBfunction.out", ios::out ); // | ios::app);
     dout_txt << "if( Atagger == \"" << CutandStrength << "\" && sEtamin == \"" << sEtaMin << "\" && sEtamax == \"" << sEtaMax << "\")" << endl;
     dout_txt << "{" << endl;
     dout_txt << "if( meanminmax == \"mean\" ) tmpSFl = new TF1(\"SFlight\",\""<< TSFUn << "\", 20.,ptmax);" << endl;
     dout_txt << "if( meanminmax == \"min\" ) tmpSFl = new TF1(\"SFlightMin\",\""<< TSHUn << "\", 20.,ptmax);" << endl;
     dout_txt << "if( meanminmax == \"max\" ) tmpSFl = new TF1(\"SFlightMax\",\""<< TSGUn << "\", 20.,ptmax);" << endl;
     dout_txt << "}" << endl;

     //Prepare writing to csv file
     cout << "Prepare writing to csv file" << endl;
     std::string taggr = CutName.Data();
     //if (CutName == "CvsL" || CutName == "CvsB") taggr = "cTag" ;
     BTagCalibration calib(taggr);

     {
       // create entry
       BTagEntry::OperatingPoint op = BTagEntry::OP_LOOSE;
       if( CutStrength == "M" ) op = BTagEntry::OP_MEDIUM;
       if( CutStrength == "T" ) op = BTagEntry::OP_TIGHT;

       BTagEntry::Parameters params(
          op,
          "ToPreciseLater", 
          "central", 
          BTagEntry::FLAV_UDSG, 
          EtaMin, EtaMax,   // eta
          20, ptmax,   // pt TEMP change to 30 from 20
          0, 1  // discr.
       );
        
        
       // two ways to instantiate an entry:
      //BTagEntry e1(params, "2*x");  // construct with function string
      BTagEntry e1(FUn, params);  // or from TF1
      params.sysType = "up";
      BTagEntry eUp(GUn, params);
      params.sysType = "down";
      BTagEntry eDown(HUn, params);

      // add entries to calibration
      calib.addEntry(e1);
      calib.addEntry(eUp);
      calib.addEntry(eDown);

     }
      
     // write to file
     std::ofstream csvOutFile(oFil+"_csv.out", ios::out);
     calib.makeCSV(csvOutFile);
     csvOutFile.close();
     std::cout << "writing to csv file done !" << endl;


     //Beginning of writing to csv file for each error type and for the different pt bins
     // Here with upper and lower values. Only for correlated and uncorr
     cout << "Prepare writing to csv file for each error and pt bins corr and uncorr " << endl;
     //const std::string taggr = CutName.Data(); // Deja fait plus haut.
     BTagCalibration Xcalib(taggr);

     {
       // create entry
       BTagEntry::OperatingPoint op = BTagEntry::OP_LOOSE;
       if( CutStrength == "M" ) op = BTagEntry::OP_MEDIUM;
       if( CutStrength == "T" ) op = BTagEntry::OP_TIGHT;
       //
       BTagEntry::Parameters params(
          op,
          "incl",
          "central",
          BTagEntry::FLAV_UDSG,
          EtaMin, EtaMax,   // eta
          20, ptmax,   // pt
          0, 1  // discr.
       );
  
       for( int iV = 0; iV < 6; ++iV )
       {
         //cout << "Test" << iV << " : " << h000->GetBinLowEdge(ViminPt[iV]) << "\t: " << ViminPt[iV] << endl;
         char buffer[50];
         std::string tmpValue;
         BTagEntry::Parameters Vparams(
          op,
          "incl",
          "central",
          BTagEntry::FLAV_UDSG,
          EtaMin, EtaMax,   // eta
          h000->GetBinLowEdge(ViminPt[iV]), h000->GetBinLowEdge(VimaxPt[iV]),   // pt
          0, 1  // discr.
         );
       
          params.ptMin = h000->GetBinLowEdge(ViminPt[iV]);
         params.ptMax = h000->GetBinLowEdge(VimaxPt[iV]);

         Double_t mean_FUn = FUn->Integral(params.ptMin,params.ptMax)/(params.ptMax-params.ptMin);
         cout << " iV = " << iV << " ptmin = " << params.ptMin << " ptmax = " << params.ptMax << " mean_FUn = " << mean_FUn << endl;
         params.sysType = "central";    //"Sample"
         sprintf(buffer,"%10.8f",mean_FUn);
         tmpValue = buffer;
         BTagEntry central(tmpValue, params);

         params.sysType = "up"; //"All+Sample"
         sprintf(buffer,"%10.8f",mean_FUn*(1.+mygeterr(false, true, totSF, ViminPt[iV], VimaxPt[iV] )));
         tmpValue = buffer;
         BTagEntry up(tmpValue, params);
         params.sysType = "down";       //"All+Sample"
         sprintf(buffer,"%10.8f",mean_FUn*(1.-mygeterr(false, true, totSF, ViminPt[iV], VimaxPt[iV] )));
         tmpValue = buffer;
         BTagEntry down(tmpValue, params);

         params.sysType = "up_correlated";      //"All+Sample"
         sprintf(buffer,"%10.8f",mean_FUn*(1.+mygeterr(false, true, totSFCorr, ViminPt[iV], VimaxPt[iV] )));
         tmpValue = buffer;
         BTagEntry up_correlated(tmpValue, params);
         params.sysType = "down_correlated";    //"All+Sample"
         sprintf(buffer,"%10.8f",mean_FUn*(1.-mygeterr(false, true, totSFCorr, ViminPt[iV], VimaxPt[iV] )));
         tmpValue = buffer;
         BTagEntry down_correlated(tmpValue, params);

         params.sysType = "up_uncorrelated";    //"Sample"
         sprintf(buffer,"%10.8f",mean_FUn*(1.+ErrSysSF));
         tmpValue = buffer;
         BTagEntry up_uncorrelated(tmpValue, params);
         params.sysType = "down_uncorrelated";  //"Sample"
         sprintf(buffer,"%10.8f",mean_FUn*(1.-ErrSysSF));
         tmpValue = buffer;
         BTagEntry down_uncorrelated(tmpValue, params);

         // add entries to calibration
         Xcalib.addEntry(central);
         Xcalib.addEntry(up);
         Xcalib.addEntry(down);
         Xcalib.addEntry(up_correlated);
         Xcalib.addEntry(down_correlated);
         Xcalib.addEntry(up_uncorrelated);
         Xcalib.addEntry(down_uncorrelated);
       }
      
     }

     // write to file
     std::ofstream XcsvOutFile(oFil+"_X_csv.out", ios::out);
     Xcalib.makeCSV(XcsvOutFile);
     XcsvOutFile.close();
     std::cout << "writing to csv file done !" << endl; 


     TString TSFun = Fun->GetExpFormula("p");
     TString TSGun = Gun->GetExpFormula("p");
     TString TSHun = Hun->GetExpFormula("p");

     ofstream dout_txt2(oFil+"_DBmistagfunction.out", ios::out ); // | ios::app);
     dout_txt2 << "if( Atagger == \"" << CutandStrength << "\" && sEtamin == \"" << sEtaMin << "\" && sEtamax == \"" << sEtaMax << "\")" << endl;
     dout_txt2 << "{" << endl;
     dout_txt2 << "if( meanminmax == \"mean\" ) tmpMistag = new TF1(\"Mistag\",\""<< TSFun << "\", 20.,ptmax);" << endl;
     dout_txt2 << "}" << endl;


     
  }
  std::cout << "Just before the return" << endl;

  return c1;
}












void LoopPlot(TH1F &retSF, TH1F &retMistag, TString tagger = "CSVv2", TString TypeOfJets = "2DSubjets", Double_t xSigma = 4.)
{
  Bool_t PrintWarnings = false;
  Bool_t Improveplot = true;
  plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 1,24, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 1, 3, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 4, 6, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 7, 9, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot,10,12, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot,13,15, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot,16,18, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot,19,24, tagger, "L", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
//
  plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 1, 24, tagger, "M", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 1, 8, tagger, "M", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 9, 16, tagger, "M", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
  //plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 17, 24, tagger, "M", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);

  plot(retSF, retMistag, 3, !PrintWarnings, Improveplot, 1, 24, tagger, "T", false, 30., 1000., 30., 1000., "Run2016", xSigma, false, 0.0, TypeOfJets);
}
