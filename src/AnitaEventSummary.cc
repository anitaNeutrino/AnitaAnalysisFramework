#include "AnitaEventSummary.h"

ClassImp(AnitaEventSummary)

const double C_IN_M_NS = 0.299792; 

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
 * Takes care of copying the header info into the event summary
 */
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* header){

  zeroInternals();

  setTriggerInfomation(header);
  eventNumber = header->eventNumber;
  run = header->run;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * Takes care of copying the header and GPS info into the event summary
 */
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* header, UsefulAdu5Pat* pat){

  zeroInternals();

  setTriggerInfomation(header);
  setSourceInformation(pat);
  eventNumber = header->eventNumber;
  run = header->run;
}




void AnitaEventSummary::zeroInternals(){

  run = 0;
  eventNumber = 0;

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    nPeaks[polInd] = 0; 
    for(Int_t dir=0; dir < maxDirectionsPerPol; dir++){
      peak[polInd][dir].value = 0;
      peak[polInd][dir].theta = 0;
      peak[polInd][dir].phi = 0;
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
  flags.isAdu5Trigger = 0;
  flags.isG12Trigger = 0;
  flags.isSoftwareTrigger = 0;
  flags.isMinBiasTrigger = 0;
  flags.isPayloadBlast = 0; 
  flags.nadirFlag = 0; 
  flags.strongCWFlag = 0;
  flags.isHPolTrigger = 0;
  flags.isVPolTrigger = 0;
  flags.pulser = EventFlags::NONE; 
  flags.isVarner = 0; 
  flags.isVarner2 = 0; 

  sun.theta = -999;
  sun.phi = -999;
  sun.distance = -999;
}


void AnitaEventSummary::setTriggerInfomation(const RawAnitaHeader* header){  

  // setting to zero, for testing this function.
  flags.isVPolTrigger = 0;
  flags.isHPolTrigger = 0;
  
  // need two adjacent L3 phi-sectors to trigger ANITA-3
  for(Int_t phi=0; phi<NUM_PHI; phi++){
    Int_t phi2 = (phi+1)%NUM_PHI; // adjacent phi-sector

    Int_t trigBit = (header->l3TrigPattern >> phi) & 1;
    Int_t trigBit2 = (header->l3TrigPattern >> phi2) & 1;

    // std::cout << header->l3TrigPattern << "\t"
    // 	      << phi << "\t" << phi2 << "\t"
    // 	      << trigBit << "\t" << trigBit2 << std::endl;

    if(trigBit > 0 && trigBit2 > 0){
      flags.isVPolTrigger = 1;      
    }

    Int_t trigBitH = (header->l3TrigPatternH >> phi) & 1;
    Int_t trigBitH2 = (header->l3TrigPatternH >> phi2) & 1;

    if(trigBitH > 0 && trigBitH2 > 0){
      flags.isHPolTrigger = 1;      
    }
  }

  flags.isRF = header->getTriggerBitRF();
  flags.isSoftwareTrigger = header->getTriggerBitADU5() | header->getTriggerBitG12() | header->getTriggerBitSoftExt();
  flags.isAdu5Trigger = header->getTriggerBitADU5();
  flags.isG12Trigger = header->getTriggerBitG12();
  flags.isSoftwareTrigger = header->getTriggerBitSoftExt();
  flags.isMinBiasTrigger = flags.isAdu5Trigger || flags.isG12Trigger || header->getTriggerBitG12() || flags.isSoftwareTrigger;

  
}




void AnitaEventSummary::setSourceInformation(UsefulAdu5Pat* pat){


  pat->getSunPosition(sun.phi, sun.theta);
  sun.distance = -999;  // I guess in theory we could compute this! 

  pat->getThetaAndPhiWaveWaisDivide(wais.theta, wais.phi);
  wais.theta *= 180/ TMath::Pi(); 
  wais.phi *= 180/ TMath::Pi(); 
  wais.distance = pat->getWaisDivideTriggerTimeNs() / C_IN_M_NS; 

  pat->getThetaAndPhiWaveLDB(ldb.theta, ldb.phi);
  ldb.distance = pat->getLDBTriggerTimeNs() / C_IN_M_NS; 
  ldb.theta *= 180/ TMath::Pi(); 
  ldb.phi *= 180/ TMath::Pi(); 

  
}
