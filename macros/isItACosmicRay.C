/////////////////////// written by Oindree Banerjee ////////////////////////////
/////////////////////// Dec 13 2017 ///////////////////////////////////////////

//#include "FFTtools.h"
//#include "AnitaVersion.h"

void isItACosmicRay(int whichrun = 410, int whichevent = 77213198, int whichanita = 3);

void isItACosmicRay(int whichrun, int whichevent, int whichanita)

{
  AnitaVersion::set(whichanita); 

  //FFTtools::loadWisdom("magicwisdom.dat"); 

  AnitaDataset d(whichrun); 

  FilterStrategy strategy;

  strategy.addOperation(new SimplePassBandFilter(0.2,1.3)); 

  UCorrelator::AnalysisConfig cfg; 

  //cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter;
  //cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg, true);

  //TFile out("test.root","RECREATE");

  //AnitaEventSummary *sum = new AnitaEventSummary; 
  
  d.getEvent(whichevent); 

  //FilteredAnitaEvent fev(d.useful(), &strategy, d.gps(), d.header());

  //analyzer.analyze(&fev,sum,d.truth()); 

  AnitaTemplateMachine atm;
  AnitaTemplateSummary ats;
  cout << "sizeof " << sizeof(AnitaTemplateSummary::coherent) << endl; 
  //atm.loadTemplates(); 

  atm.doTemplateAnalysis(analyzer->getCoherent(AnitaPol::kHorizontal,0),0,0,&ats);
  cout << ats.coherent[0][0].cRay[4] << endl; 
  



  cout << "Goodbye forever." << endl; 


}
