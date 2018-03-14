#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TGraph.h"
#include "AnalysisWaveform.h"
#include "FFTtools.h"
#include <stdio.h>
#include <sys/stat.h>

void convolveTemplates(char config, const char * numbers)
{
  TString tempname;
  TString outname;
  tempname = Form("%s/share/AnitaAnalysisFramework/templates/crTmpltsA3.root", getenv("ANITA_UTIL_INSTALL_DIR"));
  TFile f(tempname);

  outname = Form("data/templates/crTmpltsA4_%s.root", numbers);
  TFile* fOut = new TFile(outname.Data(), "RECREATE");

  TFile convfile(Form("%s/share/AnitaAnalysisFramework/tuffModels/config%c.root", getenv("ANITA_UTIL_INSTALL_DIR"), config));
  TGraph* gReal = (TGraph*) convfile.Get("gReal");
  TGraph* gImag = (TGraph*) convfile.Get("gImag");
  double phaseShift = TMath::ATan2(gImag->Eval(0), gReal->Eval(0));

  for(int disp = 0; disp < 32; disp++)
  {
    TGraph* g = (TGraph*) f.Get(Form("disp%d", disp));

    AnalysisWaveform awfH(g->GetN(), g->GetY(), g->GetX()[1] - g->GetX()[0], 0);
    (void) awfH.freq();
    double dfH = awfH.deltaF();
    int NfreqH = awfH.Nfreq();
    FFTWComplex* freqH = awfH.updateFreq();

    //convolve w tuff response
    for(int i = 0; i < NfreqH; i++)
    {
      if(dfH*i > 1.5) continue;
      double reTemp = gReal->Eval(dfH*i);
      double imTemp = gImag->Eval(dfH*i);
      FFTWComplex temp(reTemp,imTemp);
      temp.setMagPhase(temp.getAbs(), temp.getPhase() - phaseShift);
      freqH[i] *= temp;
    }
    TGraph* gOut = new TGraph(awfH.even()->GetN(), awfH.even()->GetX(), awfH.even()->GetY());
    gOut->SetName(Form("disp%d", disp));
    fOut->cd();
    gOut->Write();

    delete g;
    delete gOut;
  }

  fOut->Close();
  f.Close();
  convfile.Close();
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
