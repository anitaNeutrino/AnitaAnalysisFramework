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


void AnitaNoiseSummary::zeroInternals() {

  fifoLength=0;

  memset(avgRMSNoise,0,NUM_PHI*NUM_ANTENNA_RINGS*NUM_POLS*sizeof(double));

  for (int poli=0; poli<NUM_POLS; poli++) {
    if (avgMap[poli] != NULL) {
      delete avgMap[poli];
      avgMap[poli] = NULL;
    }
  }


  return;
}


void AnitaNoiseSummary::deleteHists() {
  
  for (int poli=0; poli<NUM_POLS; poli++) {
    if (avgMap[poli] != NULL) {
      delete avgMap[poli];
      avgMap[poli] = NULL;
    }
  }
  return;
}
 


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default Constructor
 *
 * Default constructor for ROOT
 */
AnitaNoiseMachine::AnitaNoiseMachine() {
  zeroInternals();
}




void AnitaNoiseMachine::zeroInternals() {
  memset(rmsFifo,0,NUM_PHI*NUM_ANTENNA_RINGS*NUM_POLS*fifoLength*sizeof(double)); 
  rmsFifoPos = 0;
  rmsFifoFillFlag = false;

  for (int poli=0; poli<NUM_POLS; poli++) {
    for (int loc=0; loc<fifoLength; loc++) {
      if (mapFifo[poli][loc] != NULL) {
	delete mapFifo[poli][loc];
	mapFifo[poli][loc] = NULL;
      }
    }
  }
  mapFifoPos = 0;
  rmsFifoFillFlag = false;

}


void AnitaNoiseMachine::fillAvgRMSNoise(FilteredAnitaEvent *filtered) {
  
  rmsFifoPos++;
  if (rmsFifoPos >= fifoLength) {
    rmsFifoPos = 0;
    rmsFifoFillFlag = true;
  }

  for (int phi=0; phi<NUM_PHI; phi++) {
    for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
      AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
      for (int poli=0; (AnitaPol::AnitaPol_t)poli != AnitaPol::kNotAPol; poli++) {
	AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;

	const TGraphAligned *currWave = filtered->getRawGraph(phi,ring,pol)->even();
	double value = currWave->GetRMS(1);
	rmsFifo[phi][ringi][poli][rmsFifoPos] = currWave->GetRMS(1);
      }
    }
  }
}


double AnitaNoiseMachine::getAvgRMSNoise(int phi, AnitaRing::AnitaRing_t ring, AnitaPol::AnitaPol_t pol){

  double value2 = 0;

  int endPoint = fifoLength;
  if (!rmsFifoFillFlag && phi==0 && ring == 0 && pol == 0) {
    std::cout << "WARNING in AnitaNoiseMachine::getAvgRMSNoise(): Fifo hasn't been filled entirely yet, ";
    std::cout << " only gotten " << rmsFifoPos << " out of " << fifoLength << " values" << std::endl;
  }

  for (int pos=0; pos<endPoint; pos++) {
    value2 += pow(rmsFifo[phi][(int)ring][(int)pol][pos],2)/endPoint;
  }
  return sqrt(value2);

}
  
	

void AnitaNoiseMachine::fillAvgMapNoise(UCorrelator::Analyzer *analyzer) {
  
  for (int poli=0; poli<NUM_POLS; poli++) {
    //delete old histogram if it is still there
    if (mapFifo[poli][mapFifoPos] != NULL) {
      delete mapFifo[poli][mapFifoPos];
      mapFifo[poli][mapFifoPos] = NULL;
    }
    
    //syntactically weirdly copy it out of UCorrelator
    mapFifo[poli][mapFifoPos] = (TH2D*)analyzer->getCorrelator()->getHist()->Clone();
  }
    
  
  mapFifoPos++;

  if (mapFifoPos >= fifoLength) {
    mapFifoPos = 0;
    mapFifoFillFlag = true;
  }

  return;

}


TProfile2D* AnitaNoiseMachine::getAvgMapNoise(AnitaPol::AnitaPol_t pol) {

  int nBinX = mapFifo[(int)pol][0]->GetNbinsX();
  int nBinY = mapFifo[(int)pol][0]->GetNbinsY();
  int xMin  = mapFifo[(int)pol][0]->GetXaxis()->GetBinLowEdge(1);
  int yMin  = mapFifo[(int)pol][0]->GetYaxis()->GetBinLowEdge(1);
  int xMax  = mapFifo[(int)pol][0]->GetXaxis()->GetBinUpEdge(nBinX);
  int yMax  = mapFifo[(int)pol][0]->GetYaxis()->GetBinUpEdge(nBinY);

  std::stringstream name;
  name.str("");
  name << "mapAvg_" << pol;
  
  TProfile2D *outProfile  = new TProfile2D(name.str().c_str(),"Average Interferometric Map", nBinX, xMin, xMax, nBinY, yMin, yMax);

  int lengthToDo = fifoLength;
  if (mapFifoFillFlag == false) lengthToDo = mapFifoPos;

  //stolen from Rene
  for (Int_t ix = 1; ix <= nBinX; ix++) {
    Double_t x = outProfile->GetXaxis()->GetBinCenter(ix);
    for (Int_t iy = 1; iy <= nBinY; iy++) {     
      Double_t y = outProfile->GetYaxis()->GetBinCenter(iy);
      for (Int_t i=0; i<lengthToDo; i++) {
	outProfile->Fill(x,y,mapFifo[(int)pol][i]->GetBinContent(ix,iy),1);
      }
    }
  }

  return outProfile;
}


void AnitaNoiseMachine::fillNoiseSummary(AnitaNoiseSummary *noiseSummary) {
  
  noiseSummary->fifoLength = fifoLength;
  noiseSummary->mapFifoFillFlag = mapFifoFillFlag;
  noiseSummary->rmsFifoFillFlag = rmsFifoFillFlag;


  for (int phi=0; phi<NUM_PHI; phi++) {
    for (int ringi=0; (AnitaRing::AnitaRing_t)ringi != AnitaRing::kNotARing; ringi++) {
      AnitaRing::AnitaRing_t ring = (AnitaRing::AnitaRing_t)ringi;
      for (int poli=0; (AnitaPol::AnitaPol_t)poli != AnitaPol::kNotAPol; poli++) {
	AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;

	noiseSummary->avgRMSNoise[phi][ringi][poli] = getAvgRMSNoise(phi,ring,pol);

      }
    }
  }

  if (fillMap) {
    for (int poli=0; (AnitaPol::AnitaPol_t)poli != AnitaPol::kNotAPol; poli++) {
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)poli;
      if (noiseSummary->avgMap[poli] != NULL) {
	delete noiseSummary->avgMap[poli];
	noiseSummary->avgMap[poli] = NULL;
      }
      noiseSummary->avgMap[poli] = getAvgMapNoise(pol);
    }
  }

  return;

}
