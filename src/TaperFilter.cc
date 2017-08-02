#include "TaperFilter.h" 
#include "TMath.h" 
#include "AnalysisWaveform.h" 

void GaussianTaper::processOne(AnalysisWaveform* g) 
{


  int i = 0; 
  TGraphAligned * gr = g->updateEven(); 
  while(!gr->GetX()[i]){  i++;}
  double tmin = g->even()->GetX()[i]; 


  while (gr->GetX()[i] <= tmin+mean)
  {
    gr->GetY()[i]  *= TMath::Gaus(gr->GetX()[i], tmin+mean, sigma); 
    i++; 
  }

  i = gr->GetN()-1; 
   

  while(!gr->GetX()[i]) { i--; } 

  double tmax = g->even()->GetX()[i]; 

  while (gr->GetX()[i] >= tmax-mean)
  {
    gr->GetY()[i] *= TMath::Gaus(gr->GetX()[i], tmax-mean, sigma); 
    i--; 
  }
}


