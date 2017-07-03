#ifndef ANITA_NOISE_SUMMARY
#define ANITA_NOISE_SUMMARY

#include "TObject.h"
#include "Analyzer.h"
#include "AnitaConventions.h"


class AnitaNoiseSummary
{
 public:

  static const int fifoLength = 60; //one minute of noise averaging

  /* Constructor */
  AnitaNoiseSummary();
  
  void zeroInternals();

  
  /* for calculating rms of waveform from a minute average before event capture */
  
  void fillAvgRMSNoise(FilteredAnitaEvent *filtered, double value);
  double getAvgRMSNoise(int phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol);




 private:

  //internals for fifo
  double noiseFIFO[NUM_PHI][NUM_ANTENNA_RINGS][NUM_PHI][fifoLength]; //where the info is saved
  int fifoPosition; //where in the fifo the lst write was to
  bool fifoFillFlag; //whether you've completely filled the fifo yet




#endif
