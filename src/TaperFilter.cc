#include "TaperFilter.h" 
#include "TMath.h" 
#include "AnalysisWaveform.h" 

void GaussianTaper::processOne(AnalysisWaveform* g) 
{


  TGraphAligned * gr = g->updateEven(); 

  int i = 0; 

  if (filter_beginning) 
  {
    while(!gr->GetX()[i]){  i++;}
    double tmin = g->even()->GetX()[i]; 


    while (gr->GetX()[i] <= tmin+mean)
    {
      gr->GetY()[i]  *= TMath::Gaus(gr->GetX()[i], tmin+mean, sigma); 
      i++; 
    }
  }

  if (filter_end) 
  {
    i = gr->GetN()-1; 
    while(!gr->GetX()[i]) { i--; } 
    double tmax = g->even()->GetX()[i]; 
    while (gr->GetX()[i] >= tmax-mean)
    {
      gr->GetY()[i] *= TMath::Gaus(gr->GetX()[i], tmax-mean, sigma); 
      i--; 
    }
  }
}


