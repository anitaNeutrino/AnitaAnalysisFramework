#ifndef NOISE_MONITOR_H
#define NOISE_MONITOR_H

#include "AnitaConventions.h"
#include "TString.h"
#include "TObject.h"

class FilteredAnitaEvent;
class TFile;
class TTree;


/** 
 * @class NoiseMonitor is designed to track the channel RMS of min bias events.
 * The constructor defines exactly how long to track that value in seconds, 
 * and whether to use the uneven/even (i.e. uninterpolated/interpolated) waveforms.
 * update(FilteredAnitaEvent* fEv) is expected to be called for every event in the analysis.
 * However, noise values will only be updated for non-RF triggers.
 * A non-NULL TFile pointer passed in the constructor means a TTree of RMS values will be saved.
 */
class NoiseMonitor {

 public:


  // Enumarates the possible AnalysisWaveform choice in FilteredAnitaEvent.
  enum WaveOption{
    kUneven,
    kEven
  };

  // The default timescale in seconds
  enum {
    defaultTimeScaleSeconds = 10
  };




  /**
   * @class Wrapper class for simplified ROOT IO, just holds data we'd like to persist
   */
  class FilteredMinBiasEventNoise : public TObject {
   public:
    FilteredMinBiasEventNoise(double timeScaleSeconds = defaultTimeScaleSeconds, WaveOption opt = kUneven)
        : run(0), eventNumber(0), fTimeScaleSeconds(timeScaleSeconds), fWaveOption(opt), fFilterDesc("")
    {
      for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
        for(int ant=0; ant < NUM_SEAVEYS; ant++){
          fNoise[pol][ant] = 0;
        }
      }
    }

    Int_t run; //!< For ROOT output
    UInt_t eventNumber; //< For ROOT output
    double fTimeScaleSeconds; //!< Length of time to track min bias channel RMS
    WaveOption fWaveOption; //< Which AnalysisWaveform in FilteredAnitaEvent?
    TString fFilterDesc; //< The RMS of the waveforms will depend on the filter used, so record it
    double fNoise[AnitaPol::kNotAPol][NUM_SEAVEYS]; //< The min bias channel RMS values
    ClassDef(FilteredMinBiasEventNoise, 1);
  };




  
  

  NoiseMonitor(double timeScaleSeconds = defaultTimeScaleSeconds, WaveOption opt = kUneven, TFile* outFile = NULL); // Default constructor, calculates on the fly
  NoiseMonitor(const char* inFileName); // constructor which tries to read input file.
  virtual ~NoiseMonitor(); // destructor
  
  void update(const FilteredAnitaEvent* fEv);

  /** 
   * Return the RMS noise of the channel in mV^{2} averaged over min bias events in the last fTimeScaleSeconds seconds.
   * 
   * @param pol is the polarisation of the channel
   * @param ant is the antenna
   * 
   * @return RMS noise of the channel in mV^{2}
   */
  double getNoise(AnitaPol::AnitaPol_t pol, int ant) const {return fEventNoise->fNoise[pol][ant];}
  
 protected:

  int fWriteIndex; //!< Common index for all vectors
  int removeOldEvents(double currentTime); //!< Internal function remove events older than fTimeScaleSeconds  

  // vectors to hold per-event noise data
  std::vector<double> fEventSumVSquared[AnitaPol::kNotAPol][NUM_SEAVEYS]; //!< Vector of Sum of voltage squared for all channels (where vector entries are per-event)
  std::vector<int> fEventNumPoints[AnitaPol::kNotAPol][NUM_SEAVEYS];  //!< Vector of number of entries in waveforms for all channels (where vector entries are per-event)
  std::vector<double> fEventTimes; //!< Tracks min-bias event times

  // resulting per-channel noise values
  FilteredMinBiasEventNoise* fEventNoise;
  int fReadMode;

  // Managing ROOT IO
  TFile* fOutFile; //!< Copy of pointer ROOT file given in constructor
  TFile* fInFile; //!< Copy of pointer ROOT file given in constructor
  TTree* fNoiseTree; //!< Created if non-NULL ROOT file is passed in constructor
  void prepareOutputNoiseTree(const char* treeName = getDefaultTreeName(), const char* branchName=getDefaultBranchName()); //!< Initialises fNoiseTree
  Int_t getPrecalculatedNoiseTree(const char* inFileName, const char* treeName = getDefaultTreeName(), const char* branchName=getDefaultBranchName());

  static const char* getDefaultTreeName(){return "filteredNoiseTree";}
  static const char* getDefaultBranchName(){return "filteredNoise";}

};

#endif //NOISE_MONITOR_H
