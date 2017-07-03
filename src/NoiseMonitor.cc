#include "NoiseMonitor.h"
#include "FilteredAnitaEvent.h"
#include "RawAnitaHeader.h"
#include <numeric>

#include "TTree.h"
#include "TFile.h"

/** 
 * Default constructor. Chose the time scale to monitor the noise over, exactly which waveform to check.
 * 
 * This class is intended to give you the noise value for an SNR like number.
 *
 * @param timeScaleSeconds is the number of seconds to track minimum bias events over
 * @param waveOption is the internal state of the AnalysisWaveform you want to query for the waveform RMS.
 */
NoiseMonitor::NoiseMonitor(double timeScaleSeconds, WaveOption waveOption, TFile* outFile)
    : fTimeScaleSeconds(timeScaleSeconds), fWaveOption(waveOption), fWriteIndex(-1), fOutFile(outFile),
      fNoiseTree(NULL), fRun(0), fEventNumber(0) {

  // fill noise array will zeros
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      fNoise[polInd][ant] = 0;
    }
  }
  
  prepareOutputNoiseTree();
  
}




/** 
 * Updates the noise is the FilteredAnitaEvent was not an RF trigger. Fills the noise tree if there is one.
 * 
 * @param fEv is the FilteredAnitaEvent to process
 */
void NoiseMonitor::update(const FilteredAnitaEvent* fEv){

  const RawAnitaHeader* header = fEv->getHeader();

  Bool_t isRfTrigger = header->getTriggerBitRF();
  
  if(!isRfTrigger){
    // figure out a sub-second trigger time from header info
    double tt = header->triggerTime;
    double ttNs = header->triggerTimeNs;
    tt += 1e-9*ttNs;

    // this function does a lot of work to remove old events
    // and updates the fWriteIndex, making it safe to use below
    removeOldEvents(tt);
    
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
        const AnalysisWaveform* wf = fEv->getFilteredGraph(ant, pol);

        // pick exact channel configuration based on waveform option even/uneven
        const TGraphAligned* gr = fWaveOption == kEven ? wf->even() : wf->uneven();

        int n = gr->GetN();
        double thisRMS = TMath::RMS(n,  gr->GetY());
        double sumVSquared = thisRMS*n;
        fEventSumVSquared[polInd][ant][fWriteIndex] = sumVSquared;
        fEventNumPoints[polInd][ant][fWriteIndex] = n;

        // sum over (zero meaned) v squared and nun points, to get RMS over recent min bias events
        double newSumVSquared = std::accumulate(fEventSumVSquared[pol][ant].begin(), fEventSumVSquared[pol][ant].end(), 0);
        double newNumPoints = std::accumulate(fEventNumPoints[pol][ant].begin(), fEventNumPoints[pol][ant].end(), 0);
        fNoise[pol][ant] = newNumPoints > 0 ? newSumVSquared / newNumPoints : 0;
      }
    }
  }

  // if we have an output tree, fill it
  if(fNoiseTree){
    fRun = header->run;
    fEventNumber= header->eventNumber;
    fNoiseTree->Fill();
  }
}



/** 
 * Function that replaces all entries in vectors with 0 if the event was too old.
 * Also updates the common vector index (fWriteIndex), to point to an empty index, pushing back vectors if necessary.
 * 
 * @param currentTime is the time of the current event (about to be added in update(FilteredAnitaEvent*))
 * 
 * @return the number of events removed.
 */
int NoiseMonitor::removeOldEvents(double currentTime){

  int numRemoved = 0; 
  double oldestAllowed = currentTime - fTimeScaleSeconds;

  fWriteIndex = -1; // we're going to update the vector index below

  // loop over all events in the vectors, overwriting old events if necessary
  for(unsigned i=0; i < fEventTimes.size(); i++){
   
    if(fEventTimes[i] > 0 && fEventTimes[i] < oldestAllowed){
      fEventTimes[i] = 0; // can use this to decide whether or not to reuse this element of all vectors

      for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
        for(int ant=0; ant < NUM_SEAVEYS; ant++){
          fEventSumVSquared[polInd][ant][i] = 0;
          fEventNumPoints[polInd][ant][i] = 0;
        }
      }
      numRemoved++;
    }

    // the first empty vector entry will be the new write index
    if(fWriteIndex==-1 && fEventTimes[i] == 0){
      fWriteIndex = i;
    }
  }

  // if we didn't find an empty vector entry to write to, then
  // need to push back all vectors and update fWriteIndex
  if(fWriteIndex==-1){
    fEventTimes.push_back(0);
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
        fEventSumVSquared[polInd][ant].push_back(0);
        fEventNumPoints[polInd][ant].push_back(0);
      }
    }
    fWriteIndex = fEventTimes.size() - 1;
  }

  return numRemoved;
}




/** 
 * Try to create the output noise TTree.
 * TODO, be clever and do some checking about whether file was opened in write mode.
 */
void NoiseMonitor::prepareOutputNoiseTree(){

  if(fOutFile && !fNoiseTree){
    
    fOutFile->cd();

    // use the TTree title to provide a short description of the constructor parameters
    TString treeTitle = TString::Format("Tree of channel RMS values averaged over %lf seconds of min bias events using the ", fTimeScaleSeconds);
    treeTitle += fWaveOption == kEven ? "even" : "uneven";
    treeTitle += " waveforms.";

    fNoiseTree = new TTree("minBiasNoiseTree", treeTitle);
    fNoiseTree->Branch("run", &fRun);
    fNoiseTree->Branch("eventNumber", &fEventNumber);
    TString noiseBranchName = TString::Format("noise[%d][%d]", AnitaPol::kNotAPol, NUM_SEAVEYS);
    fNoiseTree->Branch(noiseBranchName, fNoise, noiseBranchName+"/D");
  }
}
