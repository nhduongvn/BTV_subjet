{
  //need to change the input and output folder in 2DSubJets_plot_tagger.C arround line 255 also change the lumi 
  gROOT->ProcessLine(".L 2DSubJets_plot_tagger.C++") ;

  plot("CSVv2", "400", 1,  8,  50, "Run2016", 0.1, 0.3) ;
  plot("CSVv2", "400", 9,  16, 50, "Run2016", 0.1, 0.3) ;
  plot("CSVv2", "400", 17, 24, 50, "Run2016", 0.1, 0.3) ;
  plot("CSVv2", "400", 1,  24, 50, "Run2016", 0.1, 0.3) ;

  plot("cMVAv2", "400", 1,  8,  50, "Run2016", 0.1, 0.3) ;
  plot("cMVAv2", "400", 9,  16, 50, "Run2016", 0.1, 0.3) ;
  plot("cMVAv2", "400", 17, 24, 50, "Run2016", 0.1, 0.3) ;
  plot("cMVAv2", "400", 1,  24, 50, "Run2016", 0.1, 0.3) ;


}
