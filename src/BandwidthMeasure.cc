#include "BandwidthMeasure.h" 
#include "AnitaVersion.h"
#include "AnalysisWaveform.h" 
#include "TGraphAligned.h" 
#include "FFTtools.h"
#include "TFile.h"
#include <fstream> 
#include <cmath> 
#include <iostream>
#include <string>
#include <cstdlib> 

double bandwidth::bandwidthMeasure(const AnalysisWaveform* wf, int timeCheck, TGraph* gTest) 
{
  const TGraphAligned* gPow = wf->power();
  std::vector<double> powers;
  double notch0, notch1, notch2;
  notch0=0; notch1=0; notch2=0;
  //bandwidth::checkNotches(timeCheck, notch0, notch1, notch2);

  double norm = bandwidth::fillPowers(gPow, powers, notch0, notch1, notch2);
  int N = powers.size();
  std::sort(powers.begin(), powers.end());

  TGraph* gTemp = new TGraph(N);
  gTemp->SetPoint(0, 0, powers[N-1]/norm);
  for(int i = N-2; i > -1; i--)
  {
    gTemp->SetPoint(i, i, powers[i]/norm + gTemp->GetY()[i-1]);
  }
  if(gTest)
  {
    for(int i = 0; i < gTemp->GetN(); i++) gTest->SetPoint(i, gTemp->GetX()[i], gTemp->GetY()[i]);
  }
  double cdf = 0;
  for(int i = 0; i < N; i ++) cdf += gTemp->GetY()[i]/double(N);
  delete gTemp;
  powers.clear();
  cdf = fabs(cdf-1)*2;
  return cdf;
}

double bandwidth::differenceFromImpulse(const AnalysisWaveform* wf, int timeCheck, TGraph* gTest) 
{
  const TGraphAligned* gPow = wf->power();
  //TGraph* gImpRe = bandwidth::loadImpulsePower(timeCheck);
  TGraph* gImp = bandwidth::loadImpulsePower(timeCheck);
  TGraph* gImpRe = bandwidth::downsampleImpulse(gImp, gPow);
  bandwidth::normalizePower(gImpRe);
  double retVal = 0;
  double norm = 0;
  for(int i = 0; i < gPow->GetN(); i++)
  {
    if(gPow->GetX()[i] < .17 || gPow->GetX()[i] > 1.1) continue;
    norm+=gPow->GetY()[i];
  }
  for(int i = 0; i < gPow->GetN(); i++)
  {
    if(gTest) gTest->SetPoint(i, gPow->GetX()[i], fabs(gPow->GetY()[i]/norm - gImpRe->GetY()[i]));
    if(gPow->GetX()[i] < .17 || gPow->GetX()[i] > 1.1) continue;
    retVal += fabs(gPow->GetY()[i]/norm - gImpRe->GetY()[i]);
  }
  delete gImp;
  delete gImpRe;

  retVal /= 2;
  retVal = fabs(retVal - 1);
  return retVal;
}

double bandwidth::maxDifferenceFromImpulse(const AnalysisWaveform* wf, int timeCheck, TGraph* gTest) 
{
  const TGraphAligned* gPow = wf->power();
  //TGraph* gImpRe = bandwidth::loadImpulsePower(timeCheck);
  TGraph* gImp = bandwidth::loadImpulsePower(timeCheck);
  TGraph* gImpRe = bandwidth::downsampleImpulse(gImp, gPow);
  bandwidth::normalizePower(gImpRe);
  double retVal = 0;
  double norm = 0;
  for(int i = 0; i < gPow->GetN(); i++)
  {
    if(gPow->GetX()[i] < .17 || gPow->GetX()[i] > 1.1) continue;
    norm+=gPow->GetY()[i];
  }
  double maxDiff = 0;
  for(int i = 0; i < gPow->GetN(); i++)
  {
    double diff = gPow->GetY()[i]/norm - gImpRe->GetY()[i];
    if(gTest) gTest->SetPoint(i, gPow->GetX()[i], fabs(diff));
    if(gPow->GetX()[i] < .17 || gPow->GetX()[i] > 1.1) continue;
    if(diff > maxDiff) maxDiff = diff;
  }
  delete gImp;
  delete gImpRe;

  return maxDiff;
}

double bandwidth::hooverIndex(const AnalysisWaveform * wf, int timeCheck) 
{
  const TGraphAligned* gPow = wf->power();
  std::vector<double> powers;
  double notch0, notch1, notch2;
  bandwidth::checkNotches(timeCheck, notch0, notch1, notch2);

  double norm = bandwidth::fillPowers(gPow, powers, notch0, notch1, notch2);
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
  int N = powers.size();
  double meanVal = norm/double(N);

  double theil = 0;
  for(int i = 0; i < N; i++) theil += 1./double(N) * powers[i]/meanVal * std::log(powers[i]/meanVal);

  theil = theil/std::log(double(N));
  theil = fabs(theil - 1.);
  powers.clear();
  return theil;
}

double bandwidth::powerInBand(const AnalysisWaveform* wf, double minFreq, double maxFreq) 
{
  const TGraphAligned* gPow = wf->power();

  double norm = 0;
  double fraction = 0;
  for(int i = 0; i < gPow->GetN(); i++)
  {
    if(gPow->GetX()[i] > 1.3) break;
    if(gPow->GetX()[i] < minFreq) continue;
    if(gPow->GetX()[i] < 1.3) norm += gPow->GetY()[i];
    if(gPow->GetX()[i] <= maxFreq) fraction += gPow->GetY()[i];
  }
  fraction /= norm;
  return fraction;
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
    if(timeCheck < notchTime)
    {
      notchConfig = tempConfig;
      break;
    }
  }
  std::istringstream s(notchConfig);
  std::getline(s, tempConfig, '_');
  std::getline(s, tempConfig, '_');
  notch0 = atoi(tempConfig.c_str());
  notch0 /= 1000.;
  std::getline(s, tempConfig, '_');
  notch1 = atoi(tempConfig.c_str());
  notch1 /= 1000.;
  std::getline(s, tempConfig, '_');
  notch2 = atoi(tempConfig.c_str());
  notch2 /= 1000.;
  inf.close();
}

TGraph* bandwidth::loadImpulsePower(int timeCheck)
{
  TString dir;
  if(AnitaVersion::get() != 4)
  {
    dir.Form("%s/share/AnitaAnalysisFramework/responses/SingleBRotter/all.imp", getenv("ANITA_UTIL_INSTALL_DIR"));
    TGraph* g = new TGraph(dir.Data());
    TGraph* gpow = FFTtools::makePowerSpectrum(g);
    delete g;
    return gpow;
  }
  else
  {
    dir.Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/index.txt", getenv("ANITA_UTIL_INSTALL_DIR"));
    std::ifstream inf(dir.Data());
    std::string notchConfig;
    std::string tempConfig;
    long notchTime;
    while(inf >> tempConfig >> notchTime)
    {
      if(timeCheck < notchTime)
      {
        notchConfig = tempConfig;
        break;
      }
    }
    dir.Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/averages/%s.root", getenv("ANITA_UTIL_INSTALL_DIR"), notchConfig.c_str());
    TFile f(dir.Data());
    TGraph* gpow =(TGraph*) f.Get("power");
    return gpow;
  }
}

TGraph* bandwidth::downsampleImpulse(TGraph* imp, const TGraphAligned* examp)
{
  TGraph* theReturn = new TGraph(examp->GetN());
  for(int i = 0; i < examp->GetN(); i++)
  {
    theReturn->SetPoint(i, examp->GetX()[i], imp->Eval(examp->GetX()[i]));
  }
  return theReturn;
}

void bandwidth::normalizePower(TGraph* g)
{
  double norm = 0;
  for(int i = 0; i < g->GetN(); i++)
  {
    if(g->GetX()[i] < .17 || g->GetX()[i] > 1.1) continue;
    norm += g->GetY()[i];
  }
  for(int i = 0; i < g->GetN(); i++) g->GetY()[i] = g->GetY()[i]/norm;
  return;
}

double bandwidth::fillPowers(const TGraphAligned* powdb, std::vector<double> &powers, double notch0, double notch1, double notch2)
{
  double norm = 0;
  int skip = 0;
  for(int i = 0; i < powdb->GetN(); i++)
  {
    skip = 0;
    if(powdb->GetX()[i] < .17 || powdb->GetX()[i] > 1.1) skip = 1;
    /*
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
    */
    if(!skip)
    {
      powers.push_back(powdb->GetY()[i]);
      norm += powdb->GetY()[i];
    }
  }

  return norm;
}
