#include "AnitaEventSummary.h"

ClassImp(AnitaEventSummary)


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaEventSummary::AnitaEventSummary(){
  zeroInternals();
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * Takes care of copying the header and GPS info into the event summary
 */
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* header){

  zeroInternals();

  tagTriggerAsHPolOrVPol(header);
  eventNumber = header->eventNumber;
  run = header->run;
}


void AnitaEventSummary::zeroInternals(){

  run = 0;
  eventNumber = 0;

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t dir=0; dir < maxDirectionsPerPol; dir++){
      peak[polInd][dir].value = 0;
      peak[polInd][dir].theta = 0;
      peak[polInd][dir].value = 0;
      peak[polInd][dir].snr = 0;
      peak[polInd][dir].hwAngle = 0;
      peak[polInd][dir].triggered = 0;
      peak[polInd][dir].masked = 0;
      peak[polInd][dir].latitude = 0;
      peak[polInd][dir].longitude = 0;
      peak[polInd][dir].altitude = 0;

      coherent[polInd][dir].snr = 0; 
      coherent[polInd][dir].peakHilbert = 0; 
      coherent[polInd][dir].peakVal = 0; 
      coherent[polInd][dir].bandwidth = 0;
      coherent[polInd][dir].numAntennasInCoherent = 0;
      
      deconvolved[polInd][dir].snr = 0; 
      deconvolved[polInd][dir].peakHilbert = 0; 
      deconvolved[polInd][dir].peakVal = 0; 
      deconvolved[polInd][dir].bandwidth = 0;
      deconvolved[polInd][dir].numAntennasInCoherent = 0;      
    }
  }

  flags.isGood = 0;
  flags.isRF = 0;
  flags.isSoftwareTrigger = 0;
  flags.isPayloadBlast = 0; 
  flags.nadirFlag = 0; 
  flags.strongCWFlag = 0;
  flags.isHPolTrigger = 0;
  flags.isVPolTrigger = 0;
  flags.pulser = EventFlags::NONE; 
  flags.isVarner = 0; 
  flags.isVarner2 = 0; 
  
  
}


void AnitaEventSummary::tagTriggerAsHPolOrVPol(const RawAnitaHeader* header){

  Int_t numL3sV = 0;
  Int_t numL3sH = 0;  
  for(Int_t phi=0; phi<NUM_PHI; phi++){
    
    if(((header->l3TrigPattern >> phi) & 1) > 0){
      numL3sV++;
    }
    if(((header->l3TrigPatternH >> phi) & 1) > 0){
      numL3sH++;
    }    
  }

  const int numL3TriggersForEvent = 2;
  if(numL3sV >= numL3TriggersForEvent){
    flags.isVPolTrigger = 1;
  }
  if(numL3sH >= numL3TriggersForEvent){
    flags.isHPolTrigger = 1;
  }

}
