/////////////////////// written by Oindree Banerjee ////////////////////////////
/////////////////////// Dec 13 2017 ///////////////////////////////////////////
/////////////////////// Dec 30 - 31 2017 //////////////////////////////////////

#include "FFTtools.h"
#include "AnitaVersion.h"
#include "AnitaDataset.h"
#include "FilterStrategy.h"
#include "UCFilters.h"
#include "UCUtil.h"
#include "GeomFilter.h"
#include "Correlator.h"
#include "AnalysisConfig.h"
#include "Analyzer.h"
#include "AnitaTemplates.h"

void isItACosmicRay(int whichrun = 364, int whichevent = 64175392, int whichanita = 3);   //test on a payload blast

void isItACosmicRay(int whichrun, int whichevent, int whichanita)

{
  AnitaVersion::set(whichanita); 

  //FFTtools::loadWisdom("magicwisdom.dat"); 

  AnitaDataset d(whichrun); 

  FilterStrategy strategy;

  //strategy.addOperation(new SimplePassBandFilter(0.2,1.3)); 
  strategy.addOperation(0); 
  //strategy.addOperation(geomFilter,true); 

  UCorrelator::AnalysisConfig cfg; 

  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter;
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg, true);

  //TFile out("test.root","RECREATE");

  //AnitaEventSummary *sum = new AnitaEventSummary; 
  
  d.getEvent(whichevent); 

  //FilteredAnitaEvent fev(d.useful(), &strategy, d.gps(), d.header());

  //analyzer.analyze(&fev,sum,d.truth()); 

  AnitaTemplateMachine atm;
  AnitaTemplateSummary ats;
  cout << "sizeof(AnitaTemplateSummary::coherent) is:  " << sizeof(AnitaTemplateSummary::coherent) << endl; 
  atm.loadTemplates(); 

  atm.doTemplateAnalysis(analyzer->getCoherent(AnitaPol::kHorizontal,0),0,0,&ats);
  cout << "abs(ats.coherent[0][0].cRay[4] is: " << abs(ats.coherent[0][0].cRay[4]) << endl; 
  cout << "abs(ats.coherent[0][0].cRay[3] is: " << abs(ats.coherent[0][0].cRay[3]) << endl; 
  cout << "abs(ats.coherent[0][0].cRay[2] is: " << abs(ats.coherent[0][0].cRay[2]) << endl; 
  cout << "abs(ats.coherent[0][0].cRay[1] is: " << abs(ats.coherent[0][0].cRay[1]) << endl; 
  
  cout << "Goodbye forever." << endl; 


}
