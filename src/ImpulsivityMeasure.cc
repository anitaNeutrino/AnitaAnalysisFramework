#include "ImpulsivityMeasure.h" 
#include "AnalysisWaveform.h" 
#include "TH2.h" 
#include <algorithm> 


double impulsivity::impulsivityMeasure(const AnalysisWaveform * wf, TGraph * distance_cdf,int pt, bool hilbert) 
{

  const TGraphAligned * g = hilbert ? wf->hilbertEnvelope() : wf->even(); 

  if (pt < 0) 
  {
    g->peakVal(&pt,0,-1,true); 
  }

  int N = TMath::Max(pt+1, wf->Neven()-pt+1); 

  if (distance_cdf) distance_cdf->Set(N); 


  double total = g->getSumV2(); 
  double sumv2 = 0; 


  double v= g->GetY()[pt]; 
  sumv2+=v*v; 
  if (distance_cdf)
  {
    distance_cdf->GetX()[0] = 0; 
    distance_cdf->GetY()[0] = sumv2/total; 
  }

  double ysum = sumv2; 

  int i = 0; 
  while(++i < N) 
  {
    if (pt+i < wf->Neven())
    {
      v= g->GetY()[pt+i]; 
      sumv2+=v*v; 
    }

    if (pt -i >= 0) 
    {
      v= g->GetY()[pt-i]; 
      sumv2+=v*v; 
    }

    if (distance_cdf)
    {
      distance_cdf->GetY()[i] = sumv2 / total; 
      distance_cdf->GetX()[i] = i * wf->deltaT(); 
    }

    ysum+= sumv2; 
  }


  return 2 * ysum / (N*total)-1; 
}



TH2* impulsivity::envelopogram(const AnalysisWaveform * wf, TH2 * out, int min, int max, int step, bool hilbert)
{

  /* Sanity check inputs */ 

  //make sure positive 
  if (min < 0) min =-min; 
  if (max < 0) max =-max; 
  if (step < 0) step =-step; 

  // make sure min less than max
  if (min > max) std::swap(min,max); 




  /* Generate the output histogram  */ 
  
  if (out) 
  {

    out->SetBins(wf->Neven(),
        wf->even()->GetX()[0], wf->even()->GetX()[wf->Neven()-1],
        (max-min)/step+1,
        min,max);
  }
  else
  {
    out = new TH2D("envelopogram","Envelopogram",  
        wf->Neven(), wf->even()->GetX()[0], wf->even()->GetX()[wf->Neven()-1],
        (max-min)/step+1, min,max);
       
    out->SetDirectory(0) ; 
  }

  /* loop over waveform bins */
  for (int i = 0; i< wf->Neven(); i++)
  {

    int j = 0; 
    double sumv2 = 0; 
    int ybin = 1; 
    int navg = 0;
    while (abs(j) <= max)
    {
      int ii = i +j; 
      if (ii >=0 && ii < wf->Neven())
      {
        double v = ( hilbert ? wf->hilbertEnvelope() : wf->even() )->GetY()[i+j]; 
        sumv2 += v*v; 
        navg++; 
      }
      
      if (j >= 0) 
      {
        if ((j+1) >= min && ((j+1)-min) % step == 0)
        {
          out->SetBinContent(i,ybin++, sqrt(2*sumv2/navg)); 
        }
        j = -(j+1); 
      }
      else 
      {
        j=-j; 
      }
    }
  }

  return out;

}
