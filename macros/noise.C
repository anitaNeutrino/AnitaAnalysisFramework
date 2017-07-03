#include "AnitaDataset.h"
#include "NoiseMonitor.h"
#include "FilterStrategy.h"
#include "FilteredAnitaEvent.h"
#include "TFile.h"

void noise(int run = 352){

  AnitaDataset d(run);

  FilterStrategy* fs = new FilterStrategy();

  double timeScaleSeconds = 1;

  // how compression handles storing all events vs. just min bias?
  // based on the first 10 % of run 352, it seems the full data set is 2.5 times larger  
  // for a factor of 20 more events (although it's basically the same number repeating)
  TFile* fOut1 = new TFile("noiseTestAll.root", "recreate", "", 9);
  TFile* fOut2 = new TFile("noiseTestMB.root", "recreate");
  
  NoiseMonitor nm1(timeScaleSeconds, NoiseMonitor::kUneven, fOut1);
  NoiseMonitor nm2(timeScaleSeconds, NoiseMonitor::kUneven, fOut2);  
  const int n = d.N()/10;
  const int printEvery = n/100;
  for(int i=0; i < n; i++){
    d.getEntry(i);

    FilteredAnitaEvent fEv(d.useful(), fs, d.gps(),  d.header());
    nm1.update(&fEv);

    if(!d.header()->getTriggerBitRF()){
      nm2.update(&fEv);
    }
    
    // if((i%printEvery)==0){
    //   std::cout << "done " << i << " of " << n << std::endl;
    // }
  }

  fOut1->Write();
  fOut1->Close();

  fOut2->Write();
  fOut2->Close();
  
}
