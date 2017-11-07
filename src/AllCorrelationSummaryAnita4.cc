//////////////////////////////////////////////////////////////////////////////
/////  allCorrelationSummaryAnita4.cxx        ANITA event reading class     /////
/////  Description:                                                      /////
/////     A simple class for plotting event stuff like waveforms and     /////
/////     correlations adapted to the ANITA-4 geometry                   /////
/////  Author:Peng Cao based on CorrelationSummaryAntia4.cxx          /////
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

//ANITA Includes
#include "AllCorrelationSummaryAnita4.h"

//ROOT Includes
#include "TMath.h"
#include "TStyle.h"
#include "TMinuit.h"

ClassImp(AllCorrelationSummaryAnita4);





AllCorrelationSummaryAnita4::AllCorrelationSummaryAnita4( )
{

   if(!fGeomToolAnita4)
      fGeomToolAnita4=AnitaGeomTool::Instance();
}

AllCorrelationSummaryAnita4::~AllCorrelationSummaryAnita4( )
{

}
AllCorrelationSummaryAnita4::AllCorrelationSummaryAnita4( int teventNumber, int tcentreAnt)
   : eventNumber(teventNumber),centreAntenna(tcentreAnt)
{
   setFifteenAnts(tcentreAnt);
   if(!fGeomToolAnita4)
      fGeomToolAnita4=AnitaGeomTool::Instance();
}


void AllCorrelationSummaryAnita4::fillAntPositions(AnitaPol::AnitaPol_t pol)
{
  //Now fill the antenna postions (these might become member variables)
   int _pol = int(pol);

   for(int i=0;i< NUM_ALLCORRELATIONS_ANITA4;i++) {
      fAntPhi[i][0][_pol]=fGeomToolAnita4->getAntPhiPositionRelToAftFore(firstAnt[i], pol);
      fAntPhi[i][1][_pol]=fGeomToolAnita4->getAntPhiPositionRelToAftFore(secondAnt[i], pol);
      fAntR[i][0][_pol]=fGeomToolAnita4->getAntR(firstAnt[i], pol);
      fAntR[i][1][_pol]=fGeomToolAnita4->getAntR(secondAnt[i], pol);
      fAntZ[i][0][_pol]=fGeomToolAnita4->getAntZ(firstAnt[i], pol);
      fAntZ[i][1][_pol]=fGeomToolAnita4->getAntZ(secondAnt[i], pol);
   }
   

}

Double_t AllCorrelationSummaryAnita4::getDeltaTExpected(Double_t tPhiWave, Double_t tThetaWave, Int_t pairInd, AnitaPol::AnitaPol_t pol)
{
   int _pol = int(pol);

   Double_t tanThetaW=TMath::Tan(tThetaWave);
   Double_t part1=fAntZ[pairInd][0][_pol]*tanThetaW - fAntR[pairInd][0][_pol] * TMath::Cos(tPhiWave-fAntPhi[pairInd][0][_pol]);
   Double_t part2=fAntZ[pairInd][1][_pol]*tanThetaW - fAntR[pairInd][1][_pol] * TMath::Cos(tPhiWave-fAntPhi[pairInd][1][_pol]);
   return  1e9*((TMath::Cos(tThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}

void AllCorrelationSummaryAnita4::setFifteenAnts(int tcentreAnt){
   //top ring
   fifteenAnts[0]=(tcentreAnt+14)%16;
   fifteenAnts[1]=(tcentreAnt+15)%16;
   fifteenAnts[2]=(tcentreAnt+16)%16;
   fifteenAnts[3]=(tcentreAnt+17)%16;
   fifteenAnts[4]=(tcentreAnt+18)%16;
   //mid ring
   fifteenAnts[5]=(tcentreAnt+14)%16+16;
   fifteenAnts[6]=(tcentreAnt+15)%16+16;
   fifteenAnts[7]=(tcentreAnt+16)%16+16;
   fifteenAnts[8]=(tcentreAnt+17)%16+16;
   fifteenAnts[9]=(tcentreAnt+18)%16+16;
   //bot ring
   fifteenAnts[10]=(tcentreAnt+14)%16+32;
   fifteenAnts[11]=(tcentreAnt+15)%16+32;
   fifteenAnts[12]=(tcentreAnt+16)%16+32;
   fifteenAnts[13]=(tcentreAnt+17)%16+32;
   fifteenAnts[14]=(tcentreAnt+18)%16+32;
}
