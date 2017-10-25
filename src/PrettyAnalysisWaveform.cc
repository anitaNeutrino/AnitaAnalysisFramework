//////////////////////////////////////////////////////////////////////////////
/////  PrettyAnalysisWaveform.cxx        ANITA event reading class                  /////////                                                                    /////
/////  Description:                                                      /////
/////     A simple class for plotting event stuff like waveforms and correlations/////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include "PrettyAnalysisWaveform.h"

ClassImp(PrettyAnalysisWaveform);



//My incredibly dodgy approach to fitting with MINUIT that I'm going to call from within itself
// this is horrible and dangerous and should never bve done, but hey ho there we go.
void CorSumFCNanita4(Int_t& npar, Double_t*gin,
                       Double_t&f, Double_t* par, Int_t flag)
{
  return;
   //par[0] is phiWave
   //par[1] is thetaWave
   //par[2] is numAnts (11 or 19)

  // CorrelationSummaryAnita4* myDodgyCorSumPtr = dynamic_cast<CorrelationSummaryAnita4*>(gMinuit->GetObjectFit());
  // f=myDodgyCorSumPtr->getChiSquared(par[0],par[1],11);

  //  AllCorrelationSummaryAnita4* myAllCorSumPtr = dynamic_cast<AllCorrelationSummaryAnita4*>(gMinuit->GetObjectFit());
  // f=myAllCorSumPtr->getChiSquared(par[0],par[1],11);
  //  std::cout << par[0] << "\t" << par[1] << "\t" << f << std::endl;
}

PrettyAnalysisWaveform::PrettyAnalysisWaveform(UsefulAnitaEvent* uae, FilterStrategy* strategy, Adu5Pat* pat, RawAnitaHeader* header):FilteredAnitaEvent(uae, strategy, pat, header) {
fPassBandFilter=0;
fNotchFilter=0;
fPat = pat;
}

AnalysisWaveform* PrettyAnalysisWaveform::getCorrelation(int chanIndex1, int chanIndex2) 
{
   AnalysisWaveform* wf1 = (AnalysisWaveform*) getFilteredGraph(chanIndex1);
   AnalysisWaveform* wf2 = (AnalysisWaveform*) getFilteredGraph(chanIndex2);
	 //AnalysisWaveform* wf1 = new AnalysisWaveform(gr1->GetN(), gr1->GetX(), gr1->GetY(), 1./2.6);
	 wf1->forceEvenSize(260);
	 
	 ///AnalysisWaveform* wf2 = new AnalysisWaveform(gr2->GetN(), gr2->GetX(), gr2->GetY(), 1./2.6);
	 
	 wf2->forceEvenSize(260);
	
	 AnalysisWaveform* wfCorr = AnalysisWaveform::correlation(wf1, wf2, 0, 1);

   Double_t x1,y1,x2,y2;
   wf1->even()->GetPoint(0,x1,y1);
   
	 wf2->even()->GetPoint(0,x2,y2);
   fWaveOffset=x1-x2;
   
   wf1->even()->GetPoint(1,x2,y2);
   fDeltaT=x2-x1;
   
   //delete gr1;
   //delete gr2;
	 return wfCorr;     
}

AnalysisWaveform *PrettyAnalysisWaveform::getCorrelationInterpolated(int chanIndex1, int chanIndex2, AnitaPol::AnitaPol_t pol,  Int_t npadfreq) 
{
  if(chanIndex1 <0 || chanIndex1>(NUM_DIGITZED_CHANNELS-1))
    std::cerr << "Invalid channel index:\t" << chanIndex1 << "\n";

  if(chanIndex2 <0 || chanIndex2>(NUM_DIGITZED_CHANNELS-1))
    std::cerr << "Invalid channel index:\t" << chanIndex2 << "\n";
		
	 const AnalysisWaveform* wf1 = getFilteredGraph(chanIndex1, pol);
	 const AnalysisWaveform* wf2 = getFilteredGraph(chanIndex2, pol);
	
	 AnalysisWaveform* wf_1 = new AnalysisWaveform(*wf1); 
	 AnalysisWaveform* wf_2 = new AnalysisWaveform(*wf2);
	 if (npadfreq > 0)
	 {
		wf_1->padFreq(npadfreq);
		wf_2->padFreq(npadfreq);
	 }
	
	 wf_1->padEven(1);
	 wf_2->padEven(1);

	 AnalysisWaveform* wfCorr = AnalysisWaveform::correlation(wf_1, wf_2, 0, 1);

   Double_t x1,y1,x2,y2;
   wf_1->even()->GetPoint(0,x1,y1);
   //   std::cout << 1 << "\t" << chanIndex1 << "\t" << x << "\t" << y << "\n";
   wf_2->even()->GetPoint(0,x2,y2);
   fWaveOffset=x1-x2;
   //fprintf(stderr, "x1 = %g, x2 = %g\n", x1, x2);
	 //   std::cout << 2 << "\t" << chanIndex2 << "\t" << x << "\t" << y << "\n";   
   //   grCor->GetPoint(0,x1,y1);
   wf_1->even()->GetPoint(1,x2,y2);
   fDeltaT=x2-x1;
   
   //delete gr1;
   //delete gr2;
	 delete wf_1;
	 delete wf_2;

   return wfCorr;
}


void PrettyAnalysisWaveform::fillSixAntArrays(int ant, int topAnts[3], int bottomAnts[3])
{
   //   std::cerr << "fillSixAntArrays( " << ant << ")\n";
  int top=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftBottom=-1, rightBottom=-1;

  if(ant<16) {
    top=ant;
    bottom=AnitaGeomTool::getAzimuthPartner(top);
  }
  else {
    bottom=ant;
    top=AnitaGeomTool::getAzimuthPartner(bottom);
  } 

  //  std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);
  topAnts[0]=leftTop;
  topAnts[1]=top;
  topAnts[2]=rightTop;
  bottomAnts[0]=leftBottom;
  bottomAnts[1]=bottom;
  bottomAnts[2]=rightBottom;

}


void PrettyAnalysisWaveform::fillNineAntArrays(int ant, int nineAnts[9])
{
  // Top 0-2
  // Middle 3-5
  // Bottom 6-8
  
   //   std::cerr << "fillSixAntArrays( " << ant << ")\n";
  int top=-1,middle=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftMiddle=-1, rightMiddle=-1;
  int leftBottom=-1, rightBottom=-1;

  if(ant<16) {
    top=ant;
    middle=AnitaGeomTool::getAzimuthPartner(top);
    bottom=middle+16;
  } else if (ant<32){
    middle = ant;
    top = AnitaGeomTool::getAzimuthPartner(middle);
    bottom=middle+16;
  } else {
    bottom=ant;
    middle=AnitaGeomTool::getAzimuthPartner(bottom);
    top=AnitaGeomTool::getAzimuthPartner(middle);
  } 

//   std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(middle,leftMiddle,rightMiddle);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);

  nineAnts[0] = leftTop;
  nineAnts[1] = top;
  nineAnts[2] = rightTop;
  nineAnts[3] = leftMiddle;
  nineAnts[4] = middle;
  nineAnts[5] = rightMiddle;
  nineAnts[6] = leftBottom;
  nineAnts[7] = bottom;
  nineAnts[8] = rightBottom;

}


void PrettyAnalysisWaveform::fillNextFourAntArrays(int ant, int nextFourAnts[4])
{

  int top=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftLeftTop=-1, rightRightTop=-1;
  int leftBottom=-1, rightBottom=-1;
  int leftLeftBottom=-1, rightRightBottom=-1;

  if(ant<16) {
    top=ant;
    bottom=AnitaGeomTool::getAzimuthPartner(top);
  }
  else {
    bottom=ant;
    top=AnitaGeomTool::getAzimuthPartner(bottom);
  }
  int crap;
  //  std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);
  AnitaGeomTool::getThetaPartners(leftTop,leftLeftTop,crap);
  AnitaGeomTool::getThetaPartners(rightTop,crap,rightRightTop);
  AnitaGeomTool::getThetaPartners(leftBottom,leftLeftBottom,crap);
  AnitaGeomTool::getThetaPartners(rightBottom,crap,rightRightBottom);
  nextFourAnts[0]=leftLeftTop;
  nextFourAnts[1]=rightRightTop;
  nextFourAnts[2]=leftLeftBottom;
  nextFourAnts[3]=rightRightBottom;

}

void PrettyAnalysisWaveform::fillNextSixAntArrays(int ant, int nextSixAnts[6])
{

  int top=-1,middle=-1,bottom=-1;
  int leftTop=-1, rightTop=-1;
  int leftLeftTop=-1, rightRightTop=-1;
  int leftMiddle=-1, rightMiddle=-1;
  int leftLeftMiddle=-1, rightRightMiddle=-1;
  int leftBottom=-1, rightBottom=-1;
  int leftLeftBottom=-1, rightRightBottom=-1;

  if(ant<16) {
    top=ant;
    middle=AnitaGeomTool::getAzimuthPartner(top);
    bottom=middle+16;
  } else if (ant<32){
    middle = ant;
    top = AnitaGeomTool::getAzimuthPartner(middle);
    bottom=middle+16;
  } else {
    bottom=ant;
    middle=AnitaGeomTool::getAzimuthPartner(bottom);
    top=AnitaGeomTool::getAzimuthPartner(middle);
  }


  int crap;
  //  std::cout << top << "\t" << bottom << std::endl;
  AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
  AnitaGeomTool::getThetaPartners(middle,leftMiddle,rightMiddle);
  AnitaGeomTool::getThetaPartners(bottom,leftBottom,rightBottom);
  AnitaGeomTool::getThetaPartners(leftTop,leftLeftTop,crap);
  AnitaGeomTool::getThetaPartners(rightTop,crap,rightRightTop);
  AnitaGeomTool::getThetaPartners(leftMiddle,leftLeftMiddle,crap);
  AnitaGeomTool::getThetaPartners(rightMiddle,crap,rightRightMiddle);
  AnitaGeomTool::getThetaPartners(leftBottom,leftLeftBottom,crap);
  AnitaGeomTool::getThetaPartners(rightBottom,crap,rightRightBottom);
  nextSixAnts[0]=leftLeftTop;
  nextSixAnts[1]=rightRightTop;
  nextSixAnts[2]=leftLeftMiddle;
  nextSixAnts[3]=rightRightMiddle;
  nextSixAnts[4]=leftLeftBottom;
  nextSixAnts[5]=rightRightBottom;

}


void PrettyAnalysisWaveform::fillNadirArrays(int ant, int nadirAnts[5])
{
  //   std::cerr << "fillSixAntArrays( " << ant << ")\n";
 
  //some logic to convert from top to nadir (on or between antennas)



  int nadirAntNums[NUM_PHI]={32,-1,33,-1,34,-1,35,-1,36,-1,37,-1,38,-1,39,-1};

  int left = -1; int right = -1;
  int leftLeft = -1; int rightRight = -1;
  int antBottom = ant;  

    if(ant<16) {
      antBottom=AnitaGeomTool::getAzimuthPartner(ant);
    }
      
    int nadir =  nadirAntNums[antBottom-16];

    if(nadir == -1){
      
      leftLeft = -1; 
      rightRight = -1;
	left = nadirAntNums[antBottom-17];

	if((antBottom-15)==16){
	right = 32;
      }else{
	right = nadirAntNums[antBottom-15];
      }  

    }else{
      left = -1; 
      right = -1;
      if(nadir==32){
	leftLeft = 39;
      }else{
	leftLeft = nadir - 1;
      }

      if(nadir==39){
	rightRight = 32;
      }else{
	rightRight = nadir + 1;
      }  
    }

    //only need nadir antennas

 // AnitaGeomTool::getThetaPartners(top,leftTop,rightTop);
//   AnitaGeomTool::getThetaPartners(ant,leftBottom,rightBottom);
//   AnitaGeomTool::getThetaPartners(leftTop,leftLeftTop,crap);
//   AnitaGeomTool::getThetaPartners(rightTop,crap,rightRightTop);
//   AnitaGeomTool::getThetaPartners(leftBottom,leftLeftBottom,crap);
//   AnitaGeomTool::getThetaPartners(rightBottom,crap,rightRightBottom);

//   //  std::cout << top << "\t" << bottom << std::endl;
//   AnitaGeomTool::getThetaPartners(nadir,leftTop,rightTop);
	nadirAnts[0]=leftLeft;
	nadirAnts[1]=rightRight;
	nadirAnts[2]=left;
	nadirAnts[3]=right;
	nadirAnts[4]=nadir;
//   nadirAnts[5]=rightRightTop;
//   nadirAnts[6]=ant;
//   nadirAnts[7]=top;
  

}

//Putative Analysis methods
int PrettyAnalysisWaveform::getMaxAntenna(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
  return getMaxAntennaCorrelationRollingAvg(pol, peakPtr);
  //  return getMaxAntennaCorrelation(pol, peakPtr);
  //  getMaxAntennaVSquared(pol,peakPtr);
}   


int PrettyAnalysisWaveform::getMaxAntennaVSquared(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   //Returns the antenna with the maximum power
   //Could consider changng this to make things better
   double maxVal=0;
   int maxAnt=0;
   for(int ant=0;ant<32;ant++) {
      int chanIndex=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
      Double_t rmsChan=TMath::RMS(useful->fNumPoints[chanIndex],useful->fVolts[chanIndex]);
      for(int samp=0;samp<useful->fNumPoints[chanIndex];samp++) {
	 double vSquared=useful->fVolts[chanIndex][samp]*useful->fVolts[chanIndex][samp];
	 vSquared/=rmsChan;
	 if(vSquared>maxVal) {
	    maxVal=vSquared;
	    maxAnt=ant;
	 }
      }
   }
   if(peakPtr) *peakPtr=maxVal;
   return maxAnt;
}

int PrettyAnalysisWaveform::getMaxAntennaCorrelation(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   //Returns the antenna with the lagest peak/rms in the correlation with its azimuth partner antenna
   double maxVal=0;
   int maxAnt=0;
   for(int ant=0;ant<16;ant++) {
      //Loop over the top antennas
      int otherAnt=AnitaGeomTool::getAzimuthPartner(ant);
      int ciTop=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
      int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(otherAnt,pol);

      AnalysisWaveform* wfCor = getCorrelation(ciTop,ciBottom);      
      TGraph* grCor = (TGraph*) wfCor->even();
			Double_t *y = grCor->GetY();
      Double_t peak=TMath::MaxElement(grCor->GetN(),y);
      Double_t rms=TMath::RMS(grCor->GetN(),y);
      //      FFTtools::getPeakRmsSqVal(grCor,peak,rms);
      if((peak/rms)>maxVal) {
	 maxVal=peak/rms;
	 maxAnt=ant;
	 Double_t maxTop=TMath::MaxElement(useful->fNumPoints[ciTop],useful->fVolts[ciTop]);
	 Double_t maxBottom=TMath::MaxElement(useful->fNumPoints[ciBottom],useful->fVolts[ciBottom]);
	 if(maxBottom>maxTop)
	   maxAnt=otherAnt;
      }
      delete grCor;
			delete wfCor;
	 
   }
   if(peakPtr) *peakPtr=maxVal;
   return maxAnt;
}

int PrettyAnalysisWaveform::getMaxAntennaCorrelationRollingAvg(AnitaPol::AnitaPol_t pol, Double_t *peakPtr)
{
   //Returns the antenna at the centre of three phi-secotrs with the largest peak/rms in the correlation with its azimuth partner antenna
   double maxVal=0;
   int maxAnt=0;
   double maxVals[16]={0};
   for(int ant=0;ant<16;ant++) {
      //Loop over the top antennas
      int otherAnt=AnitaGeomTool::getAzimuthPartner(ant);
      int ciTop=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
      int ciMiddle=AnitaGeomTool::getChanIndexFromAntPol(otherAnt,pol);
      int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(otherAnt+16,pol);

      AnalysisWaveform* wfCor1 = getCorrelation(ciTop,ciMiddle);      
      TGraph* grCor1 = (TGraph*) wfCor1->even();
			AnalysisWaveform* wfCor2 = getCorrelation(ciTop,ciBottom);      
      TGraph* grCor2 = (TGraph*) wfCor2->even();
			AnalysisWaveform* wfCor3 = getCorrelation(ciMiddle,ciBottom);      
      TGraph* grCor3 = (TGraph*) wfCor3->even();
			Double_t *y1 = grCor1->GetY();
      Double_t PeakRMS[3];
      PeakRMS[0]=TMath::MaxElement(grCor1->GetN(),y1)/TMath::RMS(grCor1->GetN(),y1);
      Double_t *y2 = grCor2->GetY();
      PeakRMS[1]=TMath::MaxElement(grCor2->GetN(),y2)/TMath::RMS(grCor2->GetN(),y2);

      Double_t *y3 = grCor3->GetY();
      PeakRMS[2]=TMath::MaxElement(grCor3->GetN(),y3)/TMath::RMS(grCor3->GetN(),y3);
      Double_t maxPeakRMS = TMath::MaxElement(3, PeakRMS);

      //      FFTtools::getPeakRmsSqVal(grCor,peak,rms);
      if(maxPeakRMS>maxVals[ant])
	maxVals[ant]=maxPeakRMS;

     delete grCor1;
     delete grCor2;
     delete grCor3;
     delete wfCor1;
     delete wfCor2;
     delete wfCor3;
   }
   
   for(int ant=0;ant<16;ant++) {
     int leftAnt=ant-1;
     if(leftAnt<0) leftAnt=15;
     int rightAnt=ant+1;
     if(rightAnt>15) rightAnt=0;


     int otherAnt=AnitaGeomTool::getAzimuthPartner(ant);
     int ciTop=AnitaGeomTool::getChanIndexFromAntPol(ant,pol);
     int ciMiddle=AnitaGeomTool::getChanIndexFromAntPol(otherAnt,pol);
     int ciBottom=AnitaGeomTool::getChanIndexFromAntPol(otherAnt+16,pol);
     
     Double_t newVal=maxVals[leftAnt]+maxVals[ant]+maxVals[rightAnt];

     if(newVal>maxVal) {
       maxVal=newVal;
       maxAnt=ant;
       Double_t maxTop=TMath::MaxElement(useful->fNumPoints[ciTop],useful->fVolts[ciTop]);
       Double_t maxMiddle=TMath::MaxElement(useful->fNumPoints[ciMiddle],useful->fVolts[ciMiddle]);
       Double_t maxBottom=TMath::MaxElement(useful->fNumPoints[ciBottom],useful->fVolts[ciBottom]);
       if(maxTop>maxBottom && maxTop>maxMiddle) maxAnt = ant;
       else if (maxMiddle>maxBottom) maxAnt = otherAnt;
       else maxAnt=otherAnt+16;
     }
	 
   }
   if(peakPtr) *peakPtr=maxVal;
   return maxAnt;
}

double PrettyAnalysisWaveform::getHighestSnr(int centreAntenna)
{
	while(centreAntenna>16) centreAntenna -= 16;
	double vpptop=TMath::MaxElement(useful->fNumPoints[centreAntenna], useful->fVolts[centreAntenna]) - TMath::MinElement(useful->fNumPoints[centreAntenna], useful->fVolts[centreAntenna]);
	double vppmid=TMath::MaxElement(useful->fNumPoints[centreAntenna+16], useful->fVolts[centreAntenna+16]) - TMath::MinElement(useful->fNumPoints[centreAntenna+16], useful->fVolts[centreAntenna+16]);
	double vppbot=TMath::MaxElement(useful->fNumPoints[centreAntenna+32], useful->fVolts[centreAntenna+32]) - TMath::MinElement(useful->fNumPoints[centreAntenna+32], useful->fVolts[centreAntenna+32]);
	int nPointsRms = useful->fNumPoints[centreAntenna]/5;
	double snrtop=0;
	double snrmid=0;
	double snrbot=0;
	for(int i = 0; i < nPointsRms; i++)
	{
		snrtop += useful->fVolts[centreAntenna][i]*useful->fVolts[centreAntenna][i]/double(nPointsRms);
		snrmid += useful->fVolts[centreAntenna+16][i]*useful->fVolts[centreAntenna+16][i]/double(nPointsRms);
		snrbot += useful->fVolts[centreAntenna+32][i]*useful->fVolts[centreAntenna+32][i]/double(nPointsRms);
	}
	snrtop = vpptop/sqrt(snrtop);
	snrmid = vppmid/sqrt(snrmid);
	snrbot = vppbot/sqrt(snrbot);
	snrtop = (snrtop>snrbot) ? snrtop:snrbot;
	snrtop = (snrtop>snrmid) ? snrtop:snrmid;
	return snrtop;
}
	


CorrelationSummaryAnita4 *PrettyAnalysisWaveform::getCorrelationSummaryAnita4(Int_t centreAnt,AnitaPol::AnitaPol_t pol, Int_t deltaT, Int_t eventNumber)
{
  //Gets the 11 correlations and then takes the max, rms and neighbouring maxima
  if(centreAnt<0)
    centreAnt=getMaxAntenna(pol);
  // Anita 3 has 15 antennas in 5 phi sectors
  // LLT LT CT RT RRT
  // LLM LM CM RM RRM
  // LLB LB CB RB RRB
  // Resulting in 48 correlations

  int nineAnts[9];
  fillNineAntArrays(centreAnt,nineAnts);
  int nextSixAnts[6];
  fillNextSixAntArrays(centreAnt,nextSixAnts);

   CorrelationSummaryAnita4 *theSum = new CorrelationSummaryAnita4(eventNumber, centreAnt, nineAnts,deltaT);
	 theSum->setSnr(getHighestSnr(centreAnt));

   for(int i=0;i<6;i++)
      theSum->nextSixAnts[i]=nextSixAnts[i];

   //Now can make correlations and find max, rms, etc.
   for(int corInd=0;corInd<NUM_ALLCORRELATIONS_ANITA4;corInd++) {
      //      std::cout << corInd << "\t" << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd] << "\n";
      Int_t ci1=AnitaGeomTool::getChanIndexFromAntPol(theSum->firstAnt[corInd],pol);
      Int_t ci2=AnitaGeomTool::getChanIndexFromAntPol(theSum->secondAnt[corInd],pol);
      //      std::cout << corInd << "\t"<< ci1 << " " << ci2 << "  " << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd] <<std::endl;

//       if (ci1*ci2<0) continue; // Linda added this condition
			//AnalysisWaveform* wfCorr=getCorrelationInterpolated(ci1,ci2, deltaT);
			AnalysisWaveform* wfCorr=getCorrelationInterpolated(theSum->firstAnt[corInd],theSum->secondAnt[corInd], pol, deltaT);
			TGraph* grCor = new TGraph(wfCorr->even()->GetN(), wfCorr->even()->GetX(), wfCorr->even()->GetY());

      //      theSum->rmsCorVals[corInd]=grCor->GetRMS(2);

      double *theTimes = grCor->GetX();
      double *theValues = grCor->GetY();
      
      int numPoints=grCor->GetN();
      double rmsVal=TMath::RMS(numPoints,theValues);
      int maxIndex=TMath::LocMax(numPoints,theValues);
//       double maxVal=theValues[maxIndex];

//           Double_t maxVal,rmsVal;
//           Int_t maxIndex;
//           FFTtools::getPeakRmsRectified(grCor,maxVal,rmsVal,&maxIndex);

//           FFTtools::getPeakRmsSqVal(grCor,maxVal,rmsVal,&maxIndex);
//           for(int i=0;i<grCor->GetN();i++) {
// 	    //     	std::cout << i << "\t" << theTimes[i] << "\t" << theValues[i] << "\n";
//      	 if(theValues[i]>maxVal) {
//      	    maxVal=theValues[i];
//      	    maxIndex=i;
//      	 }
//           }
      theSum->rmsCorVals[corInd]=rmsVal;
      theSum->maxCorVals[corInd]=theValues[maxIndex];
      theSum->maxCorTimes[corInd]=theTimes[maxIndex];

      //      std::cout << theSum->firstAnt[corInd] << "\t" << theSum->secondAnt[corInd]
// 		     << "\t" << theSum->maxCorTimes[corInd] 
// 		     << "\t" << theSum->maxCorVals[corInd] << "\t" 
// 		     << "\t" << (theSum->maxCorTimes[corInd]-fWaveOffset)/fDeltaT << "\t"
// 		     << fWaveOffset << "\t" << fDeltaT << std::endl;

      theSum->secondCorVals[corInd][0]=theSum->maxCorVals[corInd];
      theSum->secondCorTimes[corInd][0]=theSum->maxCorTimes[corInd];
      theSum->secondCorVals[corInd][1]=theSum->maxCorVals[corInd];
      theSum->secondCorTimes[corInd][1]=theSum->maxCorTimes[corInd];
      for(int i=maxIndex-1;i>=1;i--) {
	 if(i<1) break;	 
	 if(theValues[i]>=theValues[i-1] && theValues[i]>=theValues[i+1]) {
	    theSum->secondCorVals[corInd][0]=theValues[i];
	    theSum->secondCorTimes[corInd][0]=theTimes[i];
	    break;
	 }	  
      }
      for(int i=maxIndex+1;i<grCor->GetN();i++) {
	 if(i>=grCor->GetN()-1) break;	 
	 if(theValues[i]>=theValues[i-1] && theValues[i]>=theValues[i+1]) {
	    theSum->secondCorVals[corInd][1]=theValues[i];
	    theSum->secondCorTimes[corInd][1]=theTimes[i];
	    break;
	 }	  
			}  
			delete grCor;
			delete wfCorr;
   }
    
   // //Will add a call to
   theSum->fillErrorsAndFit();


// comment out by peng, because this thetaWave phiWave is calculated from a too simple minimization funciton(not a interferometric map). 
//So the result is useless.
   // //Set up MINUIT for the fit
   // static int firstTime=1;
   // if(firstTime) {
   //    gMinuit = new TMinuit(2);
   //    firstTime=0;
   // }
   // gMinuit->SetObjectFit(theSum);  
   // gMinuit->SetFCN(CorSumFCNanita4);
   // double par[2]={theSum->fAntPhi[1][0],0};               // the start values
   // double stepSize[2]={0.01,0.01};          // step sizes 
   // double minVal[2]={0,-1*TMath::PiOver2()};            // minimum bound on parameter 
   // double maxVal[2]={TMath::TwoPi(),TMath::PiOver2()};            // maximum bound on parameter
   // char parName[2][20];
   // sprintf(parName[0],"phiWave");
   // sprintf(parName[1],"thetaWave");
   // for (int i=0; i<2; i++){
   //    gMinuit->DefineParameter(i, parName[i], par[i], stepSize[i], minVal[i], maxVal[i]);
   // }
   
   // Double_t phiWave,thetaWave;
   // Double_t phiWaveErr,thetaWaveErr;
   // //do the fit and get the answers
   // gMinuit->SetPrintLevel(-1);

   // gMinuit->Migrad();       // Minuit's best minimization algorithm   

   // gMinuit->GetParameter(0,phiWave,phiWaveErr);
   // gMinuit->GetParameter(1,thetaWave,thetaWaveErr);

//    Int_t npari,nparx,istat;
//    Double_t fmin,fedm,errdef;
//    gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
// //    std::cout << fmin << "\t" << fedm << "\t" << npari << "\t" << nparx 
// //    	     << "\t" << istat << std::endl;
//    theSum->setFitResults(phiWave,thetaWave,phiWaveErr,thetaWaveErr,fmin);

///comment out ended.

   // The better way is just use adu5 and get expected Wais theta and phi.
  Double_t phiWave,thetaWave;
  UsefulAdu5Pat* usefulPat = new UsefulAdu5Pat(fPat);
  usefulPat->getThetaAndPhiWave(AnitaLocations::LONGITUDE_WAIS_A4, AnitaLocations::LATITUDE_WAIS_A4, AnitaLocations::ALTITUDE_WAIS_A4, thetaWave, phiWave);
  theSum->thetaWave = thetaWave;
  theSum->phiWave = phiWave;

   for(int corInd=0;corInd<NUM_ALLCORRELATIONS_ANITA4;corInd++) {
      //fill the expected Time delay
      theSum->expectedDeltaT[corInd] = theSum->getDeltaTExpected(phiWave, thetaWave, corInd);
    }



   return theSum;
}


AllCorrelationSummaryAnita4 *PrettyAnalysisWaveform::getAllCorrelationSummaryAnita4(Int_t centreAnt,AnitaPol::AnitaPol_t pol, Int_t deltaT, Int_t eventNumber)
{
  if(centreAnt<0) centreAnt=getMaxAntenna(pol);
   AllCorrelationSummaryAnita4 *theSum = new AllCorrelationSummaryAnita4(eventNumber, centreAnt);
   // fill the thetaWave and phiWave
  Double_t phiWave,thetaWave;
  UsefulAdu5Pat* usefulPat = new UsefulAdu5Pat(fPat);
  usefulPat->getThetaAndPhiWave(AnitaLocations::LONGITUDE_WAIS_A4, AnitaLocations::LATITUDE_WAIS_A4, AnitaLocations::ALTITUDE_WAIS_A4, thetaWave, phiWave);
  theSum->thetaWave = thetaWave;
  theSum->phiWave = phiWave;
  int corInd = 0;
   //Now we can make correlations an
  for(int i = 0; i < 14; i++){
    for(int j = i + 1; j < 15; j++){
      int ant1 = theSum->fifteenAnts[i];
      int ant2 = theSum->fifteenAnts[j];
      if((ant2-ant1+16+2)%16>4){  // do not consider ant pair that separate larger than 2 phi sectors
        continue;
      } 
      // start filling the pair
      Int_t ci1=AnitaGeomTool::getChanIndexFromAntPol(ant1,pol);
      Int_t ci2=AnitaGeomTool::getChanIndexFromAntPol(ant2,pol);
      AnalysisWaveform* wfCorr=getCorrelationInterpolated(ant1, ant2, pol, deltaT);
      TGraph* grCor = new TGraph(wfCorr->even()->GetN(), wfCorr->even()->GetX(), wfCorr->even()->GetY());
      double *theTimes = grCor->GetX();
      double *theValues = grCor->GetY();
      int numPoints=grCor->GetN();
      double rmsVal=TMath::RMS(numPoints,theValues);
      int maxIndex=TMath::LocMax(numPoints,theValues);
      theSum->maxCorVals[corInd]=theValues[maxIndex];
      theSum->maxCorTimes[corInd]=theTimes[maxIndex];
      
      theSum->firstAnt[corInd] = ant1;
      theSum->secondAnt[corInd] = ant2;
      corInd++;
      delete grCor;
      delete wfCorr;
    }
  }
  theSum->fillAntPositions();
  for(int corInd=0; corInd<78; corInd++){
    theSum->expectedDeltaT[corInd] = theSum->getDeltaTExpected(phiWave, thetaWave, corInd);
  }
  return theSum;
}
    

 



