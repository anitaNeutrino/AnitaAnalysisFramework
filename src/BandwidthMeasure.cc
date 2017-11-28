#include "BandwidthMeasure.h" 
#include "AnitaVersion.h"
#include "AnalysisWaveform.h" 
#include "TGraphAligned.h" 
#include <fstream> 
#include <cmath> 

double bandwidth::bandwidthMeasure(const AnalysisWaveform * wf, int timeCheck) 
{
  const TGraphAligned* gPow = wf->power();
  std::vector<double> powers;
  double notch0, notch1, notch2;
  bandwidth::checkNotches(timeCheck, notch0, notch1, notch2);

  double norm = bandwidth::fillPowers(gPow, powers, notch0, notch1, notch2);
  std::sort(powers.begin(), powers.end(), [](double a, double b) { 
      return b < a;
      });
  int N = powers.size();

  TGraph* gTemp = new TGraph(N);
  gTemp->SetPoint(0, 0, powers[0]/norm);
  for(int i = 1; i < N; i++)
  {
    gTemp->SetPoint(i, i, powers[i]/norm + gTemp->GetY()[i-1]);
  }
  double cdf = 0;
  for(int i = 0; i < N; i ++) cdf += gTemp->GetY()[i]/double(N);
  delete gTemp;
  powers.clear();
  cdf = fabs(cdf-1)*2;
  return cdf;
}

double bandwidth::giniIndex(const AnalysisWaveform * wf, int timeCheck) 
{
  const TGraphAligned* gPow = wf->power();
  std::vector<double> powers;
  double notch0, notch1, notch2;
  bandwidth::checkNotches(timeCheck, notch0, notch1, notch2);

  double norm = bandwidth::fillPowers(gPow, powers, notch0, notch1, notch2);
  std::sort(powers.begin(), powers.end());

  int N = powers.size();
  
  TGraph* gTemp = new TGraph(N);
  gTemp->SetPoint(0, 0, powers[0]/norm);
  double gini = 1./(double(N)) - gTemp->GetY()[0];
  for(int i = 1; i < N; i++)
  {
    gTemp->SetPoint(i, i, powers[i]/norm + gTemp->GetY()[i-1]);
    gini += (1. * (i+1))/double(N) - gTemp->GetY()[i];
  }
  delete gTemp;
  powers.clear();
  gini /= double(N);
  gini = fabs(gini - 1);
  return gini;
}

double bandwidth::hooverIndex(const AnalysisWaveform * wf, int timeCheck) 
{
  const TGraphAligned* gPow = wf->power();
  std::vector<double> powers;
  double notch0, notch1, notch2;
  bandwidth::checkNotches(timeCheck, notch0, notch1, notch2);

  double norm = bandwidth::fillPowers(gPow, powers, notch0, notch1, notch2);
  std::sort(powers.begin(), powers.end(), [](double a, double b) { 
      return b < a;
      });
  int N = powers.size();
  double meanVal = norm/double(N);

  double hoover = 0;
  for(int i = 0; i < N; i++) hoover += .5 * fabs(powers[i] - meanVal)/norm;
  hoover = fabs(hoover - 1.);
  powers.clear();
  return hoover;
}

double bandwidth::theilIndex(const AnalysisWaveform * wf, int timeCheck) 
{
  const TGraphAligned* gPow = wf->power();
  std::vector<double> powers;
  double notch0, notch1, notch2;
  bandwidth::checkNotches(timeCheck, notch0, notch1, notch2);

  double norm = bandwidth::fillPowers(gPow, powers, notch0, notch1, notch2);
  std::sort(powers.begin(), powers.end(), [](double a, double b) { 
      return b < a;
      });
  int N = powers.size();
  double meanVal = norm/double(N);

  double theil = 0;
  for(int i = 0; i < N; i++) theil += 1./double(N) * powers[i]/meanVal * std::log(powers[i]/meanVal);

  theil = theil/std::log(double(N));
  theil = fabs(theil - 1.);
  powers.clear();
  return theil;
}

double bandwidth::alternateBandwidthMeasure(const AnalysisWaveform * wf, int timeCheck, double powerThreshold) 
{
  const TGraphAligned* gPow = wf->power();
  std::vector<double> powers;
  double notch0, notch1, notch2;
  bandwidth::checkNotches(timeCheck, notch0, notch1, notch2);
  
  double norm = bandwidth::fillPowers(gPow, powers, notch0, notch1, notch2);
  std::sort(powers.begin(), powers.end(), [](double a, double b) { 
      return b < a;
      });

  int N = powers.size();

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

void bandwidth::checkNotches(int timeCheck, double& notch0, double& notch1, double& notch2)
{
  if(AnitaVersion::get() != 4) return;
  TString dir;
  dir.Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/index.txt", getenv("ANITA_UTIL_INSTALL_DIR"));
  std::ifstream inf(dir.Data());
  std::string notchConfig;
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


double bandwidth::fillPowers(const TGraphAligned* powdb, std::vector<double> &powers, double notch0, double notch1, double notch2)
{
  double norm = 0;
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

  return norm;
}
