#ifndef ALLCORRELATIONSUMMARYANITA4_H
#define ALLCORRELATIONSUMMARYANITA4_H

#include "TObject.h"
#include "TGraph.h"
#include "AnitaGeomTool.h"



#define NUM_ALLCORRELATIONS_ANITA4 78 ///< The number of correlations from 15 antennas (5 phi-sectors) in ANITA-4.



//!  This class will have all pairs of antennas correlations for ANITA-4.

class AllCorrelationSummaryAnita4: public TObject
{


 public:
   AllCorrelationSummaryAnita4(); ///< Default constructor
   ~AllCorrelationSummaryAnita4(); ///< Destructor
   
  //! Assignment constructor.
  AllCorrelationSummaryAnita4(int teventNumber, int tcentreAnt);
  void fillAntPositions(AnitaPol::AnitaPol_t pol);
  void setFifteenAnts(int tcentreAnt);
  double getDeltaTExpected(Double_t tPhiWave, Double_t tThetaWave, Int_t pairInd, AnitaPol::AnitaPol_t pol);
  //Simple Event Characteristics
  int eventNumber; ///< The eventNumber.
  int centreAntenna;
  double maxCorVals[NUM_ALLCORRELATIONS_ANITA4];
  double maxCorTimes[NUM_ALLCORRELATIONS_ANITA4];
  double thetaWave;
  double phiWave;
  double expectedDeltaT[NUM_ALLCORRELATIONS_ANITA4];
  int fifteenAnts[15];
  int firstAnt[NUM_ALLCORRELATIONS_ANITA4];
  int secondAnt[NUM_ALLCORRELATIONS_ANITA4];
  double fAntPhi[NUM_ALLCORRELATIONS_ANITA4][2][2];// [#ofCorrelaitons][firstOrSecondAnt][pol]
  double fAntR[NUM_ALLCORRELATIONS_ANITA4][2][2];// [#ofCorrelaitons][firstOrSecondAnt][pol]
  double fAntZ[NUM_ALLCORRELATIONS_ANITA4][2][2];// [#ofCorrelaitons][firstOrSecondAnt][pol]
  ClassDef(AllCorrelationSummaryAnita4,1); ///< One of ROOT's magic macros.

private:
  AnitaGeomTool* fGeomToolAnita4 = 0;//! DOES NOT PERSIST IN ROOT!
};


#endif //ALLCORRELATIONSUMMARYANITA4_H
