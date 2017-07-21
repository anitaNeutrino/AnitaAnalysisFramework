#ifndef CORRELATIONSUMMARYANITA4_H
#define CORRELATIONSUMMARYANITA4_H

#include "TObject.h"
#include "TGraph.h"


#define NUM_CORRELATIONS_ANITA4 52 ///< The number of correlations from 15 antennas (5 phi-sectors) in ANITA-4.



//!  This is a poorly thought out class that was meant to be a summary of peaks of the correlations of some antennas for ANITA-4.
/*!
  This is a poorly thought outclass that was meant to be a summary of peaks of the correlations of some antennas for ANITA-4.
*/
class CorrelationSummaryAnita4: public TObject
{


 public:
   CorrelationSummaryAnita4(); ///< Default constructor
   ~CorrelationSummaryAnita4(); ///< Destructor
   
  //! Assignment constructor.
  /*!
    \param teventNumber The event number.
    \param tcentreAnt The centre antenna used in the correlation.
    \param nineAnts An array of the nine antennas used.
    \param deltaT The sampling period used in the interpolation (or zero)    
  */
  CorrelationSummaryAnita4(int teventNumber, int tcentreAnt, int nineAnts[6], int deltaT=0);
  void fillErrorsAndFit(); ///< The worker function that actually does the fitting.

  //! Tests a given plane wave hypothesis using a either six or ten antennas.
  /*!
    \param tPhiWave The (payload centric) azimuthal angle of the plane wave.
    \param tThetaWave The (payload centric) elevation angle of the plane wave.
    \param numAnts The number of antennas to use in the test.
    \return The chi-squared of this plane wave hypothesis.
  */
  Double_t getChiSquared(Double_t tPhiWave, Double_t tThetaWave, Int_t numAnts);
  
  //! For a given plane wave hypothesis returns the expected time difference between one of the pairs of antennas.
  /*!
    \param tPhiWave The (payload centric) azimuthal angle of the plane wave.
    \param tThetaWave The (payload centric) elevation angle of the plane wave.
    \param pairInd The index of the pair (the numbering system is the 3 top bottom pairs, 4 left-right pairs, 4 diagonal pairs, plus the 4 'neighbour' + the 4 (next phi to neighbour).
    \return The expected deltaT for this plane wave and this pair of antennas.
  */
  Double_t getDeltaTExpected(Double_t tPhiWave, Double_t tThetaWave, Int_t pairInd);
  
  void setFitResults(Double_t tPhi, Double_t tTheta, Double_t tPhiErr, Double_t tThetaErr, Double_t tChiSq); ///< Sets the result of an external fit.

	void setSnr(Double_t s) {snr = s; }
 

  //Simple Event Characteristics
  int eventNumber; ///< The eventNumber.
  int centreAntenna; ///< The number of the centre antenna.
  int nineAnts[9]; ///< The numbers of the nine central antennas.
  int nextSixAnts[6]; ///< The numbers of the six outside antennas.
  int deltaT; ///< The sampling period used.
	double snr;


  //Correlation Thingies
  
  //!  There are 49 correlations formed from a set of 15 antennas (5-phi sectors).
   /*!
     For clarity in the list below the antennas are described using the following names, for cases where there is a nadir antenna in the centre phi sector we have:
    - LLT LT CT RT RRT
    - LLM LM CM RM RRM
    - LLB LB CB RB RRB
    In this scheme the 49 correlations are:
    -  "[0]"  LT - LM
    -  "[1]"  CT - CM
    -  "[2]"  RT - RM
    -  "[3]"  LM - LB
    -  "[4]"  CM - CB
    -  "[5]"  RM - RB
    -  "[6]"  LT - CT
    -  "[7]"  CT - RT
    -  "[8]"  LM - CM
    -  "[9]"  CM - RM
    -  "[10]" LB - CB
    -  "[11]" CB - RB
    -  "[12]" LT - CM
    -  "[13]" RT - CM
    -  "[14]" LB - CM
    -  "[15]" RB - CM
    -  "[16]" LM - CT
    -  "[17]" RM - CT
    -  "[18]" LM - CB
    -  "[19]" RM - CB
    -  "[20]" LLT - CT
    -  "[21]" CT - RRT
    -  "[22]" LLM - CM
    -  "[23]" CM - RRM
    -  "[24]" LLB - CB
    -  "[25]" CB - RRB
    -  "[26]" LLT - LT
    -  "[27]" RT - RRT
    -  "[28]" LLM - LM
    -  "[29]" RM - RRM
    -  "[30]" LLB - LB
    -  "[31]" RB - RRT
    -  "[32]" LLT - LM
    -  "[33]" LLB - LM
    -  "[34]" RM - RRT
    -  "[35]" RM - RRB
    -  "[36]" LLB - LLT
    -  "[37]" LB - LT
    -  "[38]" CB - CT
    -  "[39]" RB - RT
    -  "[40]" RRB - RRT
    -  "[41]" LLT - LB
    -  "[42]" LT - CB
    -  "[43]" CT - RB
    -  "[44]" RB - RRB
    -  "[45]" LT - LLB
    -  "[46]" CT - LB
    -  "[47]" RT - CB
    -  "[48]" RRT - RB
    
   */

  int firstAnt[NUM_CORRELATIONS_ANITA4]; 
  int secondAnt[NUM_CORRELATIONS_ANITA4]; ///< The index of the second antenna in the 49 possible pairs (3 top-middle, 3 middle-bottom, 6 left-right, 6 diagonal, 6 outside-centre, 6 outside-neighbour, 4 diagonal with neighbour, 13 top-bottom combinations).
  double maxCorVals[NUM_CORRELATIONS_ANITA4]; ///< The maximum correlation value for each of the 49 possible correlations (3 top-middle, 3 middle-bottom, 6 left-right, 6 diagonal, 6 outside-centre, 6 outside-neighbour, 4 diagonal with neighbour, 13 top-bottom combinations).
  double maxCorTimes[NUM_CORRELATIONS_ANITA4]; ///< The time of the maximum correlation value for each of the 49 possible correlations (3 top-middle, 3 middle-bottom, 6 left-right, 6 diagonal, 6 outside-centre, 6 outside-neighbour, 4 diagonal with neighbour, 13 top-bottom combinations).
  double rmsCorVals[NUM_CORRELATIONS_ANITA4]; ///< The rms correlation value for each of the 49 possible correlations (3 top-middle, 3 middle-bottom, 6 left-right, 6 diagonal, 6 outside-centre, 6 outside-neighbour, 4 diagonal with neighbour, 13 top-bottom combinations).




  double secondCorVals[NUM_CORRELATIONS_ANITA4][2]; ///< The peak of the next highest correlation values (tore both left and right vals).
  double secondCorTimes[NUM_CORRELATIONS_ANITA4][2]; ///< The time of the next highest correlation values (tore both left and right vals).

  //Time Thingies
  //  double deltaT[NUM_CORRELATIONS_ANITA4]; is maxCorTimes
  double deltaTErr[NUM_CORRELATIONS_ANITA4]; ///< The error on each of the 49 deltaTs (no idea right now what this means).
  
  //Fit Results
  double phiWave; ///< The azimuthal angle of the plane wave (in payload centric coordinates).
  double thetaWave; ///< The elevation angle of the plane wave (in payload centric coordinates).
  double phiWaveErr; ///< The error on the azimuthal angle of the plane wave (in payload centric coordinates).
  double thetaWaveErr; ///< The error on the elevation angle of the plane wave (in payload centric coordinates).
  double chiSq; ///< The chi-squared of the fit.
  int ndf; ///< The number of degrees of freedom -- no frigging idea I just make it up

  //Antenna postion variables for use in fit
  Double_t fAntPhi[NUM_CORRELATIONS_ANITA4][2]; ///< A lookup table for antenna postions.
  Double_t fAntR[NUM_CORRELATIONS_ANITA4][2]; ///< A lookup table for antenna postions.
  Double_t fAntZ[NUM_CORRELATIONS_ANITA4][2]; ///< A lookup table for antenna postions.
  
  


  ClassDef(CorrelationSummaryAnita4,4); ///< One of ROOT's magic macros.
};


#endif //CORRELATIONSUMMARYANITA4_H
