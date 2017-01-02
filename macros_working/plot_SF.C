//need to change input and output directory around line 252
{
  gROOT->ProcessLine(".L 2DSubJets_mistag.C+g") ;
  TH1F retSF, retMistag ;
  LoopPlot(retSF, retMistag, "CSVv2", "2DSubjets",4.) ;
  LoopPlot(retSF, retMistag, "cMVAv2", "2DSubjets",4.) ;
}
