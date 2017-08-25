#include "NoiseMonitor.h"
#include "FilteredAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "FilterStrategy.h"
#include "FilterOperation.h"
#include "AnitaDataset.h"
#include <numeric>

#include "TTree.h"
#include "TFile.h"
#include "TProfile2D.h"
#include <stdlib.h>



UInt_t NoiseMonitor::makeStratHashFromDesc(const FilterStrategy* fs){
  TString hasher;
  for(unsigned i=0; i < fs->nOperations(); i++){
    hasher += fs->getOperation(i)->description();
  }
  UInt_t hash = hasher.Hash();
  return hash;
}


TString NoiseMonitor::getFileName(int run){
  return TString::Format("%s/rms_anita%d_run%d_hash%u.root", fRmsDir, AnitaVersion::get(), run, fHash);
}

void NoiseMonitor::getProfilesFromFile(int run){

  TString theRootPwd = gDirectory->GetPath();

  TString fileName = getFileName(run);
  TFile* fFile = TFile::Open(fileName);

  ProfPair p;
  if(fFile){
    p.set((TProfile2D*) fFile->Get(getHistName(AnitaPol::kHorizontal, run)),
          (TProfile2D*) fFile->Get(getHistName(AnitaPol::kVertical, run)));
  }

  if(!p.get(AnitaPol::kVertical) || !p.get(AnitaPol::kHorizontal)){
    makeProfiles(run);

    fFile = TFile::Open(fileName);
    p.set((TProfile2D*) fFile->Get(getHistName(AnitaPol::kHorizontal, run)),
          (TProfile2D*) fFile->Get(getHistName(AnitaPol::kVertical, run)));
  }

  if(!p.get(AnitaPol::kVertical) || !p.get(AnitaPol::kHorizontal)){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__
              << ", unable to find RMS profile! "
              << "Something bad is about to happen!"
              << std::endl;
  }
  fFiles[run] = fFile;

  fCurrent.set(p);

  gDirectory->cd(theRootPwd);
}


TString NoiseMonitor::getHistName(AnitaPol::AnitaPol_t pol, int run){
  TString name = TString::Format("rms?_%d_%u", run, fHash);
  const Ssiz_t namePolCharPos = name.First('?');
  name[namePolCharPos] = AnitaPol::polAsChar(pol);
  return name;
}


void NoiseMonitor::makeProfiles(int run){

  TString fileName = getFileName(run);

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ", generating RMS profiles for run " << run << " in file " << fileName << std::endl;
  TFile* f = TFile::Open(fileName, "recreate");

  AnitaDataset d(run);
  const int n = d.N();
  RawAnitaHeader* header = NULL;

  d.first();
  header = d.header();
  double startTime = header->realTime;

  d.last();
  header = d.header();
  double endTime = header->realTime + 1; // make sure last event is in bin

  const int nBin = (endTime - startTime)/double(defaultTimeScaleSeconds);

  // temp turn off output, turned back on at the end of the function
  std::vector<bool> original_enable_outputs = fFilterStrat->enable_outputs;
  for(unsigned i=0; i < fFilterStrat->enable_outputs.size(); i++){
    fFilterStrat->enable_outputs.at(i) = false;
  }

  TString name = TString::Format("rms?_%d_%u", run, fHash);
  TString title = TString::Format("RMS for each ?Pol for run %d (FilterStrategy operations descriptions hash = %u);Antenna;realTime", run, fHash);

  const Ssiz_t titlePolCharPos = title.First('?');

  std::vector<TProfile2D*> profs(AnitaPol::kNotAPol, NULL);
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

    TString name = getHistName(pol, run);

    title[titlePolCharPos] = AnitaPol::polAsChar(pol);
    profs.at(pol) = new TProfile2D(name, title, NUM_SEAVEYS, -0.5, NUM_SEAVEYS-0.5, nBin, startTime, endTime);
  }

  const int printEvery = n/100;
  int nextPrint = 0;
  for(int entry=0; entry < n; entry++){
    d.getEntry(entry);
    header = d.header();
    if(header->getTriggerBitRF()==0){

      FilteredAnitaEvent fEv(d.useful(), fFilterStrat, d.gps(), header);
      for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
        AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

        for(int ant=0; ant < NUM_SEAVEYS; ant++){
          const AnalysisWaveform* wf = fEv.getFilteredGraph(ant, pol);
          const TGraphAligned* gr = wf->even();
          double rms = TMath::RMS(gr->GetN(), gr->GetY());
          profs.at(pol)->Fill(ant, header->realTime, rms);
        }
      }
    }
    if(entry == nextPrint){
      std::cerr << "\r" << 100*(entry+1)/n << "% ";
      nextPrint += printEvery;
      nextPrint = nextPrint >= n ? n - 1 : nextPrint;
      if(entry==n-1){
        std::cerr << std::endl;
      }
    }
  }

  // turn the outputs back on, if they were on
  fFilterStrat->enable_outputs = original_enable_outputs;

  f->Write();
  f->Close();

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ", complete!" << std::endl;
}



/** 
 * Default constructor, requires a filter strategy
 * If you pass it a NULL, things probably won't end well
 * 
 * @param fs is the filter strategy we want to use to get the waveform RMSs
 */
NoiseMonitor::NoiseMonitor(FilterStrategy* fs)
    : fFilterStrat(fs)
{

  fRmsDir = getenv("ANITA_RMS_DIR");
  if(!fRmsDir){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__
              << ", can't see env ANITA_RMS_DIR. "
              << "Will look for/make rms profiles in pwd"
              << std::endl;
    fRmsDir = ".";
  }
  fHash = makeStratHashFromDesc(fs);
}




/** 
 * Destructor.
 */
NoiseMonitor::~NoiseMonitor(){

  std::map<int, TFile*>::iterator it;
  for(it = fFiles.begin(); it != fFiles.end(); it++){
    it->second->Close();
  }
}



/** 
 * Set fCurrent if found,
 * 
 * @param realTime is used to find the run, to find the profile in memory
 */
void NoiseMonitor::findProfilesInMemory(UInt_t realTime){

  int run = AnitaDataset::getRunAtTime(double(realTime));

  // first look in the map
  std::map<int, ProfPair>::iterator it = fRunProfiles.find(run);
  if(it!=fRunProfiles.end()){
    // then this the map
    fCurrent.set(it->second);
  }
  else{
    getProfilesFromFile(run);
  }
}


double NoiseMonitor::getRMS(AnitaPol::AnitaPol_t pol, Int_t ant, UInt_t realTime){

  const TProfile2D* p = fCurrent.get(pol);
  
  if(!(p && realTime >= fCurrent.startTime() && realTime < fCurrent.endTime())){
    findProfilesInMemory(realTime); // update fCurrent
    p = fCurrent.get(pol); // should get p again
  }

  // surely find bin should be const?
  Int_t bin = ((TProfile2D*)p)->FindBin(ant, realTime);
  return p->GetBinContent(bin);
}


void NoiseMonitor::ProfPair::set(const TProfile2D* h, const TProfile2D* v){
  H = h;
  V = v;
  // use both H and V so if one is NULL it fucks up
  // presumably the limits are the same.. but I'm not going to check right now
  fStartTime = H->GetYaxis()->GetBinLowEdge(1);
  fEndTime = V->GetYaxis()->GetBinLowEdge(V->GetNbinsY());
}

void NoiseMonitor::ProfPair::set(const ProfPair& other){

  set(other.get(AnitaPol::kHorizontal), other.get(AnitaPol::kVertical));
  
}
