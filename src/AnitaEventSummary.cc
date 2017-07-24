#include "AnitaEventSummary.h"

#include "RawAnitaHeader.h" 
#include "UsefulAdu5Pat.h"
#include "TruthAnitaEvent.h" 

//are these even necessary anymore? Who knows! 
//ClassImp(AnitaEventSummary)
//ClassImp(AnitaEventSummary::PointingHypothesis)
//ClassImp(AnitaEventSummary::SourceHypothesis)
  

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
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* header)
    : anitaLocation() {

  zeroInternals();

  setTriggerInfomation(header);
  eventNumber = header->eventNumber;
  run = header->run;
  realTime = header->realTime;

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * Takes care of copying the header and GPS info into the event summary
 */
AnitaEventSummary::AnitaEventSummary(const RawAnitaHeader* header, UsefulAdu5Pat* pat, const TruthAnitaEvent * truth) :
    anitaLocation(dynamic_cast<Adu5Pat*>(pat))
{

  zeroInternals();

  setTriggerInfomation(header);
  setSourceInformation(pat,truth);
  eventNumber = header->eventNumber;
  run = header->run;
  realTime = header->realTime;
}




void AnitaEventSummary::zeroInternals(){

  run = 0;
  eventNumber = 0;

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    nPeaks[polInd] = 0; 
    for(Int_t dir=0; dir < maxDirectionsPerPol; dir++)
    {
      memset(&peak[polInd][dir],0,sizeof(PointingHypothesis)); 
      memset(&coherent[polInd][dir],0,sizeof(WaveformInfo)); 
      memset(&deconvolved[polInd][dir],0,sizeof(WaveformInfo)); 
      memset(&coherent_filtered[polInd][dir],0,sizeof(WaveformInfo)); 
      memset(&deconvolved_filtered[polInd][dir],0,sizeof(WaveformInfo)); 
    }
  }
  memset(&flags,0, sizeof(EventFlags)); 
  flags.pulser = EventFlags::NONE; 

  sun.reset();
  wais.reset(); 
  ldb.reset(); 
  mc.reset(); 

}


/** 
 * Utility function to get the polarisation of the higher map peak
 * Useful for TTree::Draw() 
 * 
 * @param peakInd which highest peak. The default = 0, which is the highest (not sure how much sense it makes for peakInd > 0)
 * 
 * @return the polarisation which has the higher interferometric map peak (kNotAPol if peakInd required outside allowed range)
 */
AnitaPol::AnitaPol_t AnitaEventSummary::higherPeakPol(int peakInd){

  AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol;
  if(peakInd >= 0 && peakInd < maxDirectionsPerPol){
    if(peak[AnitaPol::kVertical][peakInd].value >= peak[AnitaPol::kHorizontal][peakInd].value){
      pol = AnitaPol::kVertical;
    }
    else{
      pol = AnitaPol::kHorizontal;
    }
  }
  return pol;
}



/** 
 * Utility function to return a const reference to the higher map peak. 
 * Useful for TTree::Draw() 
 * 
 * @param peakInd peakInd which highest peak. The default = 0, which is the highest (not sure how much sense it makes for peakInd > 0)
 * 
 * @return the 
 */
const AnitaEventSummary::PointingHypothesis& AnitaEventSummary::higherPeak(int peakInd){
  AnitaPol::AnitaPol_t pol = higherPeakPol();
  return peak[pol][peakInd];
}


/** 
 * Utility function to return a const reference to the the unfiltered coherently summed waveform info of the polarisation of the higher map peak
 * Useful for TTree::Draw() 
 * 
 * @param peakInd peakInd which highest peak. The default = 0, which is the highest (not sure how much sense it makes for peakInd > 0)
 * 
 * @return the unfiltered coherently summed waveform info
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::higherCoherent(int peakInd){
  AnitaPol::AnitaPol_t pol = higherPeakPol();
  return coherent[pol][peakInd];
}


/** 
 * Utility function to return a const reference to the the unfiltered and deconvolved coherently summed waveform info of the polarisation of the higher map peak
 * Useful for TTree::Draw() 
 * 
 * @param peakInd peakInd which highest peak. The default = 0, which is the highest (not sure how much sense it makes for peakInd > 0)
 * 
 * @return the unfiltered, deconvolved coherently summed waveform info
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::higherDeconvolved(int peakInd){
  AnitaPol::AnitaPol_t pol = higherPeakPol();
  return deconvolved[pol][peakInd];
}




/** 
 * Utility function to return a const reference to the the filtered coherently summed waveform info of the polarisation of the higher map peak
 * Useful for TTree::Draw() 
 * 
 * @param peakInd peakInd which highest peak. The default = 0, which is the highest (not sure how much sense it makes for peakInd > 0)
 * 
 * @return the filtered coherently summed waveform info
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::higherCoherentFiltered(int peakInd){
  AnitaPol::AnitaPol_t pol = higherPeakPol();
  return coherent_filtered[pol][peakInd];
}


/** 
 * Utility function to return a const reference to the the filtered and deconvolved coherently summed waveform info of the polarisation of the higher map peak
 * Useful for TTree::Draw() 
 * 
 * @param peakInd peakInd which highest peak. The default = 0, which is the highest (not sure how much sense it makes for peakInd > 0)
 * 
 * @return the filtered, deconvolved coherently summed waveform info
 */
const AnitaEventSummary::WaveformInfo& AnitaEventSummary::higherDeconvolvedFiltered(int peakInd){
  AnitaPol::AnitaPol_t pol = higherPeakPol();
  return deconvolved_filtered[pol][peakInd];
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


void AnitaEventSummary::setSourceInformation(UsefulAdu5Pat* pat, const TruthAnitaEvent * truth){


  pat->getSunPosition(sun.phi, sun.theta);
  sun.distance = 150e9;  // I guess in theory we could compute this! 

  pat->getThetaAndPhiWaveWaisDivide(wais.theta, wais.phi);
  wais.theta *= 180/ TMath::Pi(); 
  wais.phi *= 180/ TMath::Pi(); 
  wais.distance = pat->getWaisDivideTriggerTimeNs() * C_IN_M_NS; 

  pat->getThetaAndPhiWaveLDB(ldb.theta, ldb.phi);
  ldb.distance = pat->getLDBTriggerTimeNs() * C_IN_M_NS; 
  ldb.theta *= 180/ TMath::Pi(); 
  ldb.phi *= 180/ TMath::Pi(); 



  if (truth) 
  {
    pat->getThetaAndPhiWave(truth->sourceLon, truth->sourceLat, truth->sourceAlt, mc.theta,mc.phi);
    mc.theta*=TMath::RadToDeg();
    mc.phi*=TMath::RadToDeg();
  }
  
}


void AnitaEventSummary::MCTruth::reset()
{
 theta = -999;
 phi = -999; 
 memset(&wf[0],0,sizeof(WaveformInfo)); 
 memset(&wf[1],0,sizeof(WaveformInfo)); 
}





/** 
 * Utility function to return the linear polarization fraction so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the linear polarization fraction
 */

double AnitaEventSummary::WaveformInfo::linearPolFrac() {

  double value = TMath::Sqrt( pow(Q,2) + pow(U,2) ) / I;

  return value;
}


/** 
 * Utility function to return the linear polarization angle so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the linear polarization angle is degrees
 */

double AnitaEventSummary::WaveformInfo::linearPolAngle() {
  
  double value = (TMath::ATan(U/Q)/2)*TMath::RadToDeg();

  return value;

}



/** 
 * Utility function to return the circular polarization fraction so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the circular polarization fraction
 */

double AnitaEventSummary::WaveformInfo::circPolFrac() {
  
  double value = TMath::Abs(V)/I;

  return value;

}

/** 
 * Utility function to return the total polarization fraction so I can stop typing it
 * Useful for TTree::Draw() 
 * 
 * 
 * @return the circular polarization fraction
 */

double AnitaEventSummary::WaveformInfo::totalPolFrac() {
  
  double value = TMath::Sqrt(pow(Q,2) + pow(U,2) + pow(V,2))/I;

  return value;

}


/** 
 * Make the source hypothesis go back to nonsense thats easy to recognize
 * 
 * 
 * 
 */

void AnitaEventSummary::SourceHypothesis::reset() {

  theta = -999;
  phi = -999;
  distance = -999;

  memset(mapHistoryVal,0,NUM_POLS*sizeof(Double_t));
  memset(mapValue,0,NUM_POLS*sizeof(Double_t));
}
  



/** 
 * Let the constructor do the hard work.
 * 
 * @param pat is ANITA's gps data
 * 
 */
AnitaEventSummary::PayloadLocation::PayloadLocation(const Adu5Pat* pat){
  update(pat);
}


/** 
 * Copy all the values of interest from the Adu5Pat
 * 
 * @param pat is ANITA's gps data
 */
void AnitaEventSummary::PayloadLocation::update(const Adu5Pat* pat){
  if(pat){
    latitude = pat->latitude;
    longitude = pat->longitude;
    altitude = pat->altitude;
    heading = pat->heading;
  }
  else{
    reset();
  }
}
