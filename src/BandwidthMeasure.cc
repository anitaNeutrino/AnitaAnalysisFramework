#include "BandwidthMeasure.h" 
#include "AnalysisWaveform.h" 
#include "TGraphAligned.h" 
#include <algorithm> 

double bandwidth::bandwidthMeasure(const AnalysisWaveform * wf, TGraph* testGraph) 
{
  const TGraphAligned* powdb = wf->power();
  std::vector<double> powers;
  double norm = 0;
  for(int i = 0; i < powdb->GetN(); i++)
  {
    if(powdb->GetX()[i] < .180 || powdb->GetX()[i] > 1.1) continue;
    if(.240 <= powdb->GetX()[i] <= .280) continue; //dont look at TUFFs (may make this smarter later)
    if(.355 <= powdb->GetX()[i] <= .395) continue;
    if(.440 <= powdb->GetX()[i] <= .480) continue;
    powers.push_back(powdb->GetY()[i]);
    norm += powdb->GetY()[i];
  }
  int N = powers.size();

  std::sort(powers.begin(), powers.end(), [](double a, double b) { 
      return b < a;
      });

  TGraph* gTemp = new TGraph(N);
  if(testGraph) testGraph->Set(N);
  gTemp->SetPoint(0, 0, powers[0]/norm);
  if(testGraph) testGraph->SetPoint(0,0, powers[0]/norm);
  for(int i = 1; i < N; i++)
  {
    gTemp->SetPoint(i, i, powers[i]/norm + gTemp->GetY()[i-1]);
    if(testGraph) testGraph->SetPoint(i, i, powers[i]/norm + gTemp->GetY()[i-1]);
  }
  double cdf = 0;
  for(int i = 0; i < N; i ++) cdf += gTemp->GetY()[i]/double(N);
  delete gTemp;
  powers.clear();
  cdf = fabs(cdf-1)*2;
  return cdf;
}
