#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TGraph.h"
#include "FFTtools.h"
#include "AnitaTemplates.h"
#include <stdio.h>
#include <sys/stat.h>

void convolveTemplates(char config, const char * numbers)
{
  TString outname;

  outname = Form("data/templates/crTmpltsA4_%s.root", numbers);
  TFile* fOut = new TFile(outname.Data(), "RECREATE");
  
  TFile* inf = new TFile("data/templates/crTmpltsA3.root");

  TGraph* gconv = new TGraph(Form("%s/share/AnitaAnalysisFramework/responses/A4ImpulseTUFFs/averages/%s.imp", getenv("ANITA_UTIL_INSTALL_DIR"), numbers));
  AnalysisWaveform awfc(gconv->GetN(), gconv->GetX(), gconv->GetY(), .1);
  awfc.forceEvenSize(10000);
  (void) awfc.freq();
  double dfc = awfc.deltaF();
  int Nfreqc = awfc.Nfreq();
  FFTWComplex* freqc = awfc.updateFreq();

  for(int disp = 0; disp < 32; disp++)
  {
    TGraph* g = (TGraph*) inf->Get(Form("efield%d", disp));

    AnalysisWaveform awfH(g->GetN(), g->GetX(), g->GetY(), .1);
    awfH.forceEvenSize(10000);
    (void) awfH.freq();
    double dfH = awfH.deltaF();
    int NfreqH = awfH.Nfreq();
    FFTWComplex* freqH = awfH.updateFreq();
    
    //convolve w tuff response
    for(int i = 0; i < NfreqH; i++)
    {
      //if( dfH*i > 1.4) continue;
      freqH[i] *= freqc[i];
    }
    TGraph* gOut = new TGraph(awfH.even()->GetN(), awfH.even()->GetX(), awfH.even()->GetY());
    gOut->SetName(Form("disp%d", disp));
    fOut->cd();
    gOut->Write();

    delete gOut;
    delete g;
  }
  inf->Close();
  fOut->Close();
  delete gconv;
  delete inf;
  delete fOut;
}

void makeConvolvedTemplates()
{
    convolveTemplates('B', "notches_260_375_0");
    convolveTemplates('J', "notches_250_375_0");
    convolveTemplates('C', "notches_260_0_460");
    convolveTemplates('A', "notches_260_0_0");
    convolveTemplates('O', "notches_260_365_0");
    convolveTemplates('G', "notches_260_385_0");
    convolveTemplates('P', "notches_260_375_460");
}
