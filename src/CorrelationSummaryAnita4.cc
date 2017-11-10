//////////////////////////////////////////////////////////////////////////////
/////  CorrelationSummaryAnita4.cxx        ANITA event reading class     /////
/////  Description:                                                      /////
/////     A simple class for plotting event stuff like waveforms and     /////
/////     correlations adapted to the ANITA-4 geometry                   /////
/////  Author: Linda Cremonesi (l.cremonesi@ucl.ac.uk)                   /////
/////         based on CorrelarationSummary.cxx written by               ///// 
/////         Ryan Nichol (rjn@hep.ucl.ac.uk)                            /////
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//ANITA Includes
#include "CorrelationSummaryAnita4.h"
#include "AnitaGeomTool.h"

//ROOT Includes
#include "TMath.h"
#include "TStyle.h"
#include "TMinuit.h"

#include "AnitaGeomTool.h"


ClassImp(CorrelationSummaryAnita4);


//Global ish things.
AnitaGeomTool *fCSGeomToolAnita4=0;


CorrelationSummaryAnita4::CorrelationSummaryAnita4( )
{

   if(!fCSGeomToolAnita4)
      fCSGeomToolAnita4=AnitaGeomTool::Instance();
}

CorrelationSummaryAnita4::~CorrelationSummaryAnita4( )
{

}
CorrelationSummaryAnita4::CorrelationSummaryAnita4( int teventNumber, int tcentreAnt, int tnineAnts[], int dt, AnitaPol::AnitaPol_t tpol)
   : eventNumber(teventNumber),centreAntenna(tcentreAnt),deltaT(dt), pol(tpol)
{
  for(int i=0;i<9;i++) 
    nineAnts[i]=tnineAnts[i];

   if(!fCSGeomToolAnita4)
      fCSGeomToolAnita4=AnitaGeomTool::Instance();
}


void CorrelationSummaryAnita4::fillErrorsAndFit()
{
   //At some stage this should probably do something, otherwise it won't be much cop.

   //For now we'll go for the zeroth order solution to errors that is will assign them all to be the same
   //and just for teh sake of doing something we'll make it half a time bin (0.5/2.6) in ns.
   for(int i=0;i< NUM_CORRELATIONS_ANITA4;i++) {
      deltaTErr[i]=0.5/2.6; // in ns
   }

  //Now fill the antenna postions (these might become member variables)
   for(int i=0;i< NUM_CORRELATIONS_ANITA4;i++) {
      fAntPhi[i][0]=fCSGeomToolAnita4->getAntPhiPositionRelToAftFore(firstAnt[i], pol);
      fAntPhi[i][1]=fCSGeomToolAnita4->getAntPhiPositionRelToAftFore(secondAnt[i], pol);
      fAntR[i][0]=fCSGeomToolAnita4->getAntR(firstAnt[i], pol);
      fAntR[i][1]=fCSGeomToolAnita4->getAntR(secondAnt[i], pol);
      fAntZ[i][0]=fCSGeomToolAnita4->getAntZ(firstAnt[i], pol);
      fAntZ[i][1]=fCSGeomToolAnita4->getAntZ(secondAnt[i], pol);
   }
   

}

Double_t CorrelationSummaryAnita4::getChiSquared(Double_t tPhiWave, Double_t tThetaWave, Int_t numAnts)
{
   Double_t chiSq=0;
   for(int i=0;i<numAnts;i++) {
      Double_t dtExpect=getDeltaTExpected(tPhiWave,tThetaWave,i);
      chiSq+=(maxCorTimes[i]-dtExpect)*(maxCorTimes[i]-dtExpect)/(deltaTErr[i]*deltaTErr[i]);
      //      std::cout << i << "\t" << chiSq << "\t" << maxCorTimes[i] 
      //		<< "\t" << dtExpect << "\t"
      //		<<  deltaTErr[i] << "\n";
   }
   return chiSq;
}

Double_t CorrelationSummaryAnita4::getDeltaTExpected(Double_t tPhiWave, Double_t tThetaWave, Int_t pairInd)
{

   Double_t tanThetaW=TMath::Tan(tThetaWave);
   Double_t part1=fAntZ[pairInd][0]*tanThetaW - fAntR[pairInd][0] * TMath::Cos(tPhiWave-fAntPhi[pairInd][0]);
   Double_t part2=fAntZ[pairInd][1]*tanThetaW - fAntR[pairInd][1] * TMath::Cos(tPhiWave-fAntPhi[pairInd][1]);
   return  1e9*((TMath::Cos(tThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}

void CorrelationSummaryAnita4::setFitResults(Double_t tPhi, Double_t tTheta, Double_t tPhiErr, Double_t tThetaErr, Double_t tChiSq) {
   phiWave=tPhi;
   thetaWave=tTheta;
   phiWaveErr=tPhiErr;
   thetaWaveErr=tThetaErr;
   chiSq=tChiSq;
}
