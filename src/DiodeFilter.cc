#include "DiodeFilter.h" 
#include "TMath.h" 


DiodeFilter::DiodeFilter() 
  : response(260)
{

  TGraphAligned * g = response.updateEven(); 

  double A1 = -0.8; 
  double A2 = -0.2; 
  double A3 = 0.00964;  
  double s1 = 2.3; 
  double s2 = 4; 
  double s3 = 7; 
  double t1 = 15; 
  double t2 = 15; 
  double t3 = 18; 


  for (int i = 0; i < g->GetN(); i++) 
  {

    double t = g->GetX()[i]; 
    g->GetY()[i] =  A1 * TMath::Gaus(t,t1,s1) + A2 * TMath::Gaus(t,t2,s2) + A3 * pow(t-t3,2) * TMath::Gaus(t,t3,s3); 
  }

  response.padEven(1); 
}


void DiodeFilter::processOne(AnalysisWaveform * wf, const RawAnitaHeader *, int,int)
{

  //we start with the hilbert envelope of the waveform 

  AnalysisWaveform power(wf->Neven(), wf->hilbertEnvelope()->GetY(), wf->deltaT(), wf->even()->GetX()[0]);  
  TGraphAligned * g = power.updateEven(); 
  for (int i = 0; i < g->GetN(); i++) g->GetY()[i] *= g->GetY()[i]/50;

  if (power.Neven() < response.Neven()) 
  {
    power.forceEvenSize(response.Neven()); 
  }
  else if (wf->Neven() > response.Neven())
  {
    //I guess we'll pad our diode filter... this is probably not the smartest thing in the world but oh well 
    response.forceEvenSize(wf->Neven()); 
  }

  AnalysisWaveform * conv = AnalysisWaveform::convolution(&response, &power); 
  *wf=*conv; 
  delete conv; 
}


