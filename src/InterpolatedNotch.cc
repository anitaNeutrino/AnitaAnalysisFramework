#include "InterpolatedNotch.h"
#include "FFTtools.h"
using namespace std;
void interpolatedNotchFilter::init(){
  notchLower,notchHigher = false;
  lowBin,lowSigma,highBin,highSigma =0;
}
interpolatedNotchFilter::interpolatedNotchFilter(){
  init();
}
TGraph * interpolatedNotchFilter::interpolatedFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq) //interpolated filter. Do interpolation across notch
{
  if(maxFreq>1200){
    maxFreq=1200;
  }
  //TRandom3 random; // for generating random numbers
  //random.SetSeed(fHeadPtr->eventNumber);
  //cout<<"minFreq is "<<minFreq<< " max freq is "<<maxFreq<<"\n";
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();//256
  double frequencyArray[2000];
  double magFFT[2000];
  for(int j0=0;j0<2000;j0++){
    frequencyArray[j0]=0;
    magFFT[j0]=0;
  }
  FFTWComplex *theFFT=FFTtools::doFFT(length,oldY);
  int newLength=(length/2)+1;
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  double tempF=0.;
  double startFreq=0.;
  double endFreq=0.;
  int i_start=0;
  int i_end=0;
  //cout<<"minFreq is "<<minFreq<<" maxFreq is "<<maxFreq<<"\n";
  for(int i=0;i<newLength;i++){
    //cout<<"minFreq is "<<minFreq<<" maxFreq is "<<maxFreq<<" tempF is "<<tempF<<"\n";
    if(tempF >= (minFreq-deltaF) && tempF <=(minFreq) ){
      //  cout<<"setting starts \n";
      startFreq = tempF;
      i_start = i;
    }
    if(tempF >= (maxFreq) && tempF <= (maxFreq+deltaF)){
      //  cout<<"Setting ends \n";
      endFreq = tempF;
      i_end=i;
    }
    tempF+=deltaF;
  }
  tempF=0;
  double slope;
  double intercept;
  double mag_start = sqrt(pow(theFFT[i_start].re,2) + pow(theFFT[i_start].im,2));
  double mag_end = sqrt(pow(theFFT[i_end].re,2) + pow(theFFT[i_end].im,2));
  slope = (mag_end - mag_start)/(endFreq - startFreq);
  intercept = mag_start - (slope*startFreq);
  double mag=0.;
  double phi_old=0.;
  double x=0.;
  double y=0.;
  //cout<<slope<<endl;
  for(int i=0;i<newLength;i++) {
    if(tempF>startFreq && tempF<endFreq ) {
      phi_old = atan2(theFFT[i].im,theFFT[i].re);//old phase
      mag = slope*tempF + intercept;
      x = mag*cos(phi_old);
      y = mag*sin(phi_old);
      theFFT[i].re = x;
      theFFT[i].im = y;
    }//in inside freq range
    tempF+=deltaF;
  }
  double *filteredVals = FFTtools::doInvFFT(length,theFFT);
  TGraph *grFilteredNotch=new TGraph(length,oldX,filteredVals);  //oldX
  // }
  delete [] theFFT;
  delete [] filteredVals;
  return grFilteredNotch;
}//interpolated notch
void interpolatedNotchFilter::cutParams(bool notchLower,Double_t lowBinCenter,Double_t lowBinSigma,bool notchHigher,
                    Double_t highBinCenter,Double_t highBinSigma){
  if(notchLower){
    lowBin = 1000*lowBinCenter;
    lowSigma = 2*1000*lowBinSigma;
    //cout<<"notching "<<lowBin-lowSigma<<" to "<<lowBin+lowSigma<<endl;
  } else{lowBin =0.;lowSigma =0.;}
  if(notchHigher){
    highBin = 1000*highBinCenter;
    highSigma = 2*1000*highBinSigma;
    //cout<<"notching "<<highBin-highSigma<<" to "<<highBin+highSigma<<endl;
  }else{highBin =0.;highSigma =0.;}
}
void interpolatedNotchFilter::process(FilteredAnitaEvent *ev){
  for(int ant = 0;ant<NUM_SEAVEYS;ant++){
    AnalysisWaveform * unfiltV = getWf(ev, ant, AnitaPol::kVertical);
    TGraph *notchLowV = interpolatedFilter((TGraph*)unfiltV->even(), lowBin-lowSigma, lowBin+lowSigma);
    TGraph *notchHighV = interpolatedFilter(notchLowV,highBin-highSigma,highBin+highSigma);
    AnalysisWaveform *notchedV = AnalysisWaveform::makeWf(notchHighV);
    unfiltV->updateEven(notchedV->even());
    AnalysisWaveform * unfiltH = getWf(ev, ant, AnitaPol::kHorizontal);
    TGraph *notchLowH = interpolatedFilter((TGraph*)unfiltH->even(), lowBin-lowSigma, lowBin+lowSigma);
    TGraph *notchHighH = interpolatedFilter(notchLowH,highBin-highSigma,highBin+highSigma);
    AnalysisWaveform *notchedH = AnalysisWaveform::makeWf(notchHighH);
    unfiltH->updateEven(notchedH->even());
    
    //cout<<"notched"<<endl;
  }
}
void interpolatedNotchFilter::processOne(AnalysisWaveform *wf, const RawAnitaHeader *,int, int){
  cout<<"not implemented"<<'\n';
}
