#include "AnitaNoiseSummary.h"



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaNoiseSummary::AnitaNoiseSummary() {
  zeroInternals();
}




void AnitaNoiseSumamry::zeroInternals() {
  
  memset(&noiseFIFO,0,NUM_PHI*NUM_ANTENNA_RINGS*NUM_PHI*fifoLength*sizeof(double)); 
  fifoPosition = 0;
  fifoFillFlag = false;
}


void AnitaNoiseSummary::fillAvgRMSNoise(FilteredAnitaEvent *filtered, double value) {
  
  fifoPosition++;
  if (fifoPosition >= fifoLength) {
    fifoPosition = 0;
    fifoFillFlag = true;
  }

  for (int phi=0; phi<NUM_PHI; phi++) {
    for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
      AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
      for (int poli=0; (AnitaPol::AnitaPol_t)poli != AnitaPol::kNotAPol; poli++) {
	AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;

	const TGraphAligned *currWave = getRawGraph(phi,ring,pol)->even();
	noiseFIFO[phi][ringi][poli][fifoPosition] = currWave->GetRMS(1);
      }
    }
  }
}


double AnitaNoiseSummary::getAvgRMSNoise(int phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol){

  double value2 = 0;

  for (int pos=0; pos<fifoLength; pos++) {
    value2 += pow(noiseFIFO[phi][(int)ring][(int)pol][pos],2)/fifoLength;
  }
  return TMath::Sqrt(value2);

}
  
	 
