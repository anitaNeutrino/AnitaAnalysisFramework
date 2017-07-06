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


/*===================
  A class to store information about the thermal environment */
class AnitaNoiseSummary
{
 public:
  AnitaNoiseSummary();

  //functions
  void zeroInternals();
  void deleteHists();

  //size
  int fifoLength;

  //flags
  bool isMinBias; //is it a min bias event?  Can use this to determine where the updates are
  bool mapFifoFillFlag; //has the fifo filled up yet?
  bool rmsFifoFillFlag; //has the fifo filled up yet?

  //data
  double avgRMSNoise[NUM_PHI][NUM_ANTENNA_RINGS][NUM_POLS];

  TProfile2D *avgMap[NUM_POLS] = {NULL}; //this ends up taking a HUGE amount of space



 private:
  ClassDefNV(AnitaNoiseSummary,2);
  
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
  TProfile2D *getAvgMapNoise(AnitaPol::AnitaPol_t pol);


  // basically just moves things into the summary
  void fillNoiseSummary(AnitaNoiseSummary *noiseSummary); //crabcrabcrab
  //crab



 private:

  //internals for time domain waveform rms fifo
  double rmsFifo[NUM_PHI][NUM_ANTENNA_RINGS][NUM_POLS][fifoLength]; //where the info is saved
  int rmsFifoPos; //where in the fifo the NEXT write will be
  bool rmsFifoFillFlag; //whether you've completely filled the fifo once

  //internals for interferometric map fifo (probably enormous in memory so maybe make a flag)
  TH2D *mapFifo[NUM_POLS][fifoLength] = {NULL}; //where the info is saved
  int mapFifoPos;  //where in the fifo the NEXT write will be 
  bool mapFifoFillFlag; //whether you've completely filled the fifo once.

  ClassDefNV(AnitaNoiseMachine, 2); 

};
/*------------------*/


#endif
