#include "AnitaConventions.h"

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

  // which waveforms in the FilteredAnitaEvent should be queried?
  enum WaveOption{
    kUneven,
    kEven
  };
  
  
  NoiseMonitor(double timeScaleSeconds = 1, WaveOption opt = kUneven, TFile* outFile = NULL);
  virtual ~NoiseMonitor() {;} //!< Does nothing
  
  void update(const FilteredAnitaEvent* fEv);

  /** 
   * Return the RMS noise of the channel in mV^{2} averaged over min bias events in the last fTimeScaleSeconds seconds.
   * 
   * @param pol is the polarisation of the channel
   * @param ant is the antenna
   * 
   * @return RMS noise of the channel in mV^{2}
   */
  double getNoise(AnitaPol::AnitaPol_t pol, int ant) const {return fNoise[pol][ant];}
  
 protected:

  // timescale and which waveform settings
  double fTimeScaleSeconds;
  WaveOption fWaveOption;

  // vectors to hold per-event noise data
  int fWriteIndex; //!< Common index for all vectors
  std::vector<double> fEventSumVSquared[AnitaPol::kNotAPol][NUM_SEAVEYS]; //!< Vector of Sum of voltage squared for all channels (where vector entries are per-event)
  std::vector<int> fEventNumPoints[AnitaPol::kNotAPol][NUM_SEAVEYS];  //!< Vector of number of entries in waveforms for all channels (where vector entries are per-event)
  std::vector<double> fEventTimes; //!< Tracks min-bias event times

  // resulting per-channel noise values
  double fNoise[AnitaPol::kNotAPol][NUM_SEAVEYS];

  // Optional ROOT output 
  TFile* fOutFile; //!< Copy of pointer ROOT file given in constructor
  TTree* fNoiseTree; //!< Created if non-NULL ROOT file is passed in constructor
  Int_t fRun; //!< For optional ROOT output
  UInt_t fEventNumber; //< For optional ROOT output

  int removeOldEvents(double currentTime); //!< Internal function remove events older than fTimeScaleSeconds
  void prepareOutputNoiseTree(); //!< Initialises fNoiseTree
};
