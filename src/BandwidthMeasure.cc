#include "BandwidthMeasure.h" 
#include "AnitaVersion.h"
#include "AnalysisWaveform.h" 
#include "TGraphAligned.h" 
#include <fstream> 

double bandwidth::bandwidthMeasure(const AnalysisWaveform * wf, int timeCheck, TGraph* testGraph) 
{
  const TGraphAligned* powdb = wf->power();
  std::vector<double> powers;
  double norm = 0;
  std::string notchConfig;
  double notch0, notch1, notch2;
  if(AnitaVersion::get() == 4)
  {
    TString dir;
    dir.Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/index.txt", getenv("ANITA_UTIL_INSTALL_DIR"));
    std::ifstream inf(dir.Data());
    std::string tempConfig;
    long notchTime;
    while(inf >> tempConfig >> notchTime)
    {
      if(timeCheck < notchTime) notchConfig = tempConfig;
    }
    std::istringstream s(notchConfig);
    std::getline(s, tempConfig, '_');
    std::getline(s, tempConfig, '_');
    notch0 = std::stoi(tempConfig);
    notch0 /= 1000.;
    std::getline(s, tempConfig, '_');
    notch1 = std::stoi(tempConfig);
    notch1 /= 1000.;
    std::getline(s, tempConfig, '_');
    notch2 = std::stoi(tempConfig);
    notch2 /= 1000.;
    inf.close();
  }
  int skip = 0;
  for(int i = 0; i < powdb->GetN(); i++)
  {
    skip = 0;
    if(powdb->GetX()[i] < .180 || powdb->GetX()[i] > 1.1) skip = 1;
    if(AnitaVersion::get() == 4)
    {
      if(notch0 > 0)
      {
        if(notch0 - .15 <= powdb->GetX()[i] <= notch0 + .15) skip = 1;
      }
      if(notch1 > 0)
      {
        if(notch1 - .15 <= powdb->GetX()[i] <= notch1 + .15) skip = 1;
      }
      if(notch2 > 0)
      {
        if(notch2 - .15 <= powdb->GetX()[i] <= notch2 + .15) skip = 1;
      }
    }
    if(!skip)
    {
      powers.push_back(powdb->GetY()[i]);
      norm += powdb->GetY()[i];
    }
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

double bandwidth::alternateBandwidthMeasure(const AnalysisWaveform * wf, int timeCheck, double powerThreshold, TGraph* testGraph) 
{
  const TGraphAligned* powdb = wf->power();
  std::vector<double> powers;
  double norm = 0;
  std::string notchConfig;
  double notch0, notch1, notch2;
  if(AnitaVersion::get() == 4)
  {
    TString dir;
    dir.Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/index.txt", getenv("ANITA_UTIL_INSTALL_DIR"));
    std::ifstream inf(dir.Data());
    std::string tempConfig;
    long notchTime;
    while(inf >> tempConfig >> notchTime)
    {
      if(timeCheck < notchTime) notchConfig = tempConfig;
    }
    std::istringstream s(notchConfig);
    std::getline(s, tempConfig, '_');
    std::getline(s, tempConfig, '_');
    notch0 = std::stoi(tempConfig);
    notch0 /= 1000.;
    std::getline(s, tempConfig, '_');
    notch1 = std::stoi(tempConfig);
    notch1 /= 1000.;
    std::getline(s, tempConfig, '_');
    notch2 = std::stoi(tempConfig);
    notch2 /= 1000.;
    inf.close();
  }
  int skip = 0;
  for(int i = 0; i < powdb->GetN(); i++)
  {
    skip = 0;
    if(powdb->GetX()[i] < .180 || powdb->GetX()[i] > 1.1) skip = 1;
    if(AnitaVersion::get() == 4)
    {
      if(notch0 > 0)
      {
        if(notch0 - .15 <= powdb->GetX()[i] <= notch0 + .15) skip = 1;
      }
      if(notch1 > 0)
      {
        if(notch1 - .15 <= powdb->GetX()[i] <= notch1 + .15) skip = 1;
      }
      if(notch2 > 0)
      {
        if(notch2 - .15 <= powdb->GetX()[i] <= notch2 + .15) skip = 1;
      }
    }
    if(!skip)
    {
      powers.push_back(powdb->GetY()[i]);
      norm += powdb->GetY()[i];
    }
  }
  int N = powers.size();

  std::sort(powers.begin(), powers.end(), [](double a, double b) { 
      return b < a;
      });

  if(testGraph)
  {
    testGraph->Set(N);
    for(int i = 0; i < N; i++) testGraph->SetPoint(i, i, powers[i]/norm);
  }
  int nbins = -1;
  double integratedPower = 0;
  while(integratedPower < powerThreshold)
  {
    nbins++;
    integratedPower += powers[nbins]/norm;
  }
  powers.clear();
  double val = double(nbins)/double(N) * 1./powerThreshold;
  return val;
}
