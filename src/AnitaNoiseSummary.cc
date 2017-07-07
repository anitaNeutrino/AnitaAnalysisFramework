#include "AnitaNoiseSummary.h"



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaNoiseSummary::AnitaNoiseSummary() {
  for (int poli=0; poli<NUM_POLS; poli++) {
    avgMapProf[poli] = NULL;
  }
  zeroInternals();
}

AnitaNoiseSummary::~AnitaNoiseSummary() {

  deleteHists();
}

void AnitaNoiseSummary::zeroInternals() {

  fifoLength=0;

  memset(avgRMSNoise,0,NUM_PHI*NUM_ANTENNA_RINGS*NUM_POLS*sizeof(double));

  for (int poli=0; poli<NUM_POLS; poli++) {
    if (avgMapProf[poli] != NULL) {
      delete avgMapProf[poli];
      avgMapProf[poli] = NULL;
    }
  }

  return;
}


void AnitaNoiseSummary::deleteHists() {
  
  for (int poli=0; poli<NUM_POLS; poli++) {
    if (avgMapProf[poli] != NULL) {
      delete avgMapProf[poli];
      avgMapProf[poli] = NULL;
    }
  }  

  return;
}

