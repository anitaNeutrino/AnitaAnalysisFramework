#ifndef ANITA_NOISE_SUMMARY
#define ANITA_NOISE_SUMMARY

#include "math.h"
#include <sstream>

#include "TObject.h"
#include "TH2D.h"
#include "TProfile2D.h"

#include "Analyzer.h"
#include "AnitaConventions.h"
#include "FilteredAnitaEvent.h"
#include "AnalysisConfig.h"

/*===================
  A class to store information about the thermal environment */
class AnitaNoiseSummary
{
 public:
  AnitaNoiseSummary();

  virtual ~AnitaNoiseSummary();

  //functions
  void zeroInternals();
  void deleteHists();

  //size
  int fifoLength; //comes from AnitaNoiseMachine
  static const int nPhi = 180; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  static const int nTheta = 100; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic


  //flags
  bool isMinBias; //is it a min bias event?  Can use this to determine where the updates are
  bool mapFifoFillFlag; //has the fifo filled up yet?
  bool rmsFifoFillFlag; //has the fifo filled up yet?

  //data
  double avgRMSNoise[NUM_PHI][NUM_ANTENNA_RINGS][NUM_POLS];

  TProfile2D *avgMapProf[NUM_POLS] = {NULL}; //this ends up taking a HUGE amount of space
  double avgMaps[NUM_POLS][nPhi][nTheta]; //where the array is stored.  should be nPhi*nTheta*NUM_POLS long


 private:

  ClassDefNV(AnitaNoiseSummary,3);
  
};

/*--------------*/



/*=====================
  A class to process and save information about the thermal environment*/
class AnitaNoiseMachine
{
 public:

  static const int fifoLength = 60; //one minute of noise averaging

  //do you want to do the interferometric map?  It might be resource intensive
  bool fillMap = false;

  /* Constructor */
  AnitaNoiseMachine();
  
  void zeroInternals();

  
  /* for calculating rms of waveform from a minute average before event capture */
  void fillAvgRMSNoise(FilteredAnitaEvent *filtered);
  double getAvgRMSNoise(int phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol);

  /* for building up an enormous memory block of histograms from a bunch of events, then making an average */
  void fillAvgMapNoise(UCorrelator::Analyzer *analyzer);

  /* makes a TProfile2D out of the histograms in the fifo */
  TProfile2D *getAvgMapNoiseProfile(AnitaPol::AnitaPol_t pol);

  /* updates the weird map fifo array that I'm trying to get to work */
  void updateAvgMapNoise();

  // basically just moves things into the summary
  void fillNoiseSummary(AnitaNoiseSummary *noiseSummary); //crabcrabcrab

  //crab



 private:

  //internals for time domain waveform rms fifo
  double rmsFifo[NUM_PHI][NUM_ANTENNA_RINGS][NUM_POLS][fifoLength]; //where the info is saved
  int rmsFifoPos; //where in the fifo the NEXT write will be
  bool rmsFifoFillFlag = false; //whether you've completely filled the fifo once

  //internals for interferometric map fifo (probably enormous in memory so maybe make a flag)
  TH2D *mapFifo[NUM_POLS][fifoLength] = {NULL}; //where the info is saved
  int mapFifoPos;  //where in the fifo the NEXT write will be 
  bool mapFifoFillFlag = false; //whether you've completely filled the fifo once.

  //maybe will be more compressed
  static const int nPhi = 180; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  static const int nTheta = 100; //default in UCorrelator::AnalysisConfig, this is hard to make dynamic
  double avgMaps[NUM_POLS][nPhi][nTheta]; //where the array is stored.  should be nPhi*nTheta*NUM_POLS long

  ClassDefNV(AnitaNoiseMachine, 2); 

};
/*------------------*/


#endif
