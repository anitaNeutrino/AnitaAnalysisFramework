#include "AnitaDataset.h"
#include "NoiseMonitor.h"
#include "FilterStrategy.h"
#include "FilteredAnitaEvent.h"
#include "TFile.h"



void readNoiseTree(int run=352){

  TString fileName = TString::Format("filteredEventNoise%d.root", run);
  NoiseMonitor nm(fileName);

  AnitaDataset d(run);
  const int n = d.N();
  FilterStrategy* fs = new FilterStrategy();

  TH2D* hTest = new TH2D("hTest", "Channel noise vs. entry", 128, 0, n, 128, 0, 100);

  for(int i=0; i < n; i++){
    d.getEntry(i);

    FilteredAnitaEvent fEv(d.useful(), fs, d.gps(),  d.header());
    nm.update(&fEv);

    hTest->Fill(i, nm.getNoise(AnitaPol::kHorizontal, 0));
  }
  hTest->Draw("colz");

}


/** 
 * Use this macro to generate a minimum bias noise tree just the once.
 * The channel RMSs will depend on the filter used.
 * The produced noise tree tries to track this by saving the string of the filter description.
 * 
 * @param run is the run to produce the tree for.
 */

void makeNoiseTree(int run = 352){

  AnitaDataset d(run);

  FilterStrategy* fs = new FilterStrategy();

  double timeScaleSeconds = 10;

  TFile* fNoiseFile = new TFile(TString::Format("filteredEventNoise%d.root", run), "recreate", "", 9);
  
  NoiseMonitor nm(timeScaleSeconds, NoiseMonitor::kUneven, fNoiseFile);
  const int n = d.N();
  const int printEvery = n/100;
  for(int i=0; i < n; i++){
    d.getEntry(i);

    FilteredAnitaEvent fEv(d.useful(), fs, d.gps(),  d.header());
    nm.update(&fEv);

    if((i%printEvery)==0){
      std::cout << "\r" << "done " << i << " of " << n << std::flush;
    }
  }
  std::cout << std::endl << "complete!" << std::endl;

  fNoiseFile->Write();
  fNoiseFile->Close();

  readNoiseTree(run);

}
