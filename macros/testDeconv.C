#include "FFTtools.h"
#include "forward/forward.hpp"
#include <vector>

void testDeconv(int which=13, double noise_rms = 0.2)
{ 
  FFTtools::loadWisdom("wisdom.fftw3"); 

  AnitaResponse::AllPassDeconvolution allpass;  
  AnitaResponse::BandLimitedDeconvolution bandlimited(0.18,0.65);  

  TF1 * restoring_beam = new TF1("restoring","[0]*x*x*exp(-x/[1]) * (x > 0)",-100,100); 

  restoring_beam->SetParameters(1,1.43); 
  AnitaResponse::CLEANDeconvolution clean(restoring_beam);

  // ===== START WAVELET DECONVOLUTION

  // do a 4-stage deconvolution
  const unsigned int p{4};

  // initialize a scaling array equal to 0.5 at all levels
  const std::vector<double> scaling(p+1, 0.5);

  // and a shrink parameter of 1. at all levels
  const std::vector<double> rho(p+1, 1.);

  // we want to use the Meyer basis
  const auto type{forward::WaveletType::Meyer};

  // and construct the wavelet deconvolver
  auto wavelet{AnitaResponse::WaveletDeconvolution(p, type, noise_rms, scaling, rho)};

  // ===== END WAVELET DECONVOLUTION

 // clean.enableSaveIntermediate(); 


  TCanvas * crestoring = new TCanvas("restore","restore"); 
  restoring_beam->Draw();
  crestoring->SaveAs("clean_dbg/restoring.png"); 

  AnitaResponse::ResponseManager rm("SingleBRotter",1, &allpass); 

  TFile f("data/templates/crTmpltsA3.root"); 
  TString str; 
  str.Form("efield%d",which); 
  AnalysisWaveform* orig = AnalysisWaveform::makeWf((TGraph*) f.Get(str.Data())); 
  orig->setColor(3); 
  str.Form("disp%d",which); 
  TGraph * conv = (TGraph*) f.Get(str.Data()); 

  FFTtools::ThermalNoise noise(260, 0.18,1.3, noise_rms, 2); 
  for (int i = 0; i < conv->GetN(); i++)
  {
    conv->GetY()[i] += noise.eval(conv->GetX()[i]); 
  }

  orig->updateEven()->setPlottingLimits(1.1,true,100); 
  AnitaResponse::DeconvolutionMethod  *  methods[] = { &allpass, &clean, &bandlimited, &wavelet };
  const char  *  method_names[] = { "allpass", "CLEAN", "bandlimited", "wavelet" };
  int nmethods = sizeof(methods)/sizeof(*methods);


  const AnitaResponse::AbstractResponse * resp = rm.response(0,AnitaPol::kHorizontal); 

  /*
  TCanvas * cresp = new TCanvas("resp","Response"); 
  resp->impulseResponse(0.1,1024)->drawEven("alp"); 
  cresp->SaveAs("clean_dbg/response.png"); 
  */
  


  for (int i = 0; i < nmethods; i++) 
  {
    TCanvas * c = new TCanvas(method_names[i],method_names[i]); 

    c->Divide(2,2); 

    AnalysisWaveform * wf = AnalysisWaveform::makeWf(conv); 
    wf->padEven(1,0); 
    c->cd(1); 
    wf->updateEven()->setPlottingLimits(1.1,true,100); 
    wf->drawEven(); 

    c->cd(2); 

    TStopwatch watch;
    AnalysisWaveform *dwf = resp->deconvolve(wf, methods[i]); 
    watch.Stop(); 
    printf("method: %s: %p\n", method_names[i], dwf); 
    watch.Print(); 

    dwf->setColor(2); 
    dwf->updateEven()->setPlottingLimits(1.1,true,100); 
    dwf->drawEven(); 

    c->cd(3); 
    wf->drawPower(); 
    dwf->drawPower("lsame"); 
    orig->drawPower("lsame"); 
    c->cd(4); 
    orig->drawEven(); 

    str.Form("clean_dbg/%s_%d_%g.png", method_names[i],which,noise_rms); 
    c->SaveAs(str.Data()); 
  }


  TCanvas * cintermediate = new TCanvas("cintermediate","Intermediate",1800,900); 
  clean.getComponents()->updateEven()->setPlottingLimits(1.1,true,100); 
  clean.getComponents()->drawEven(); 
  str.Form("clean_dbg/CLEAN_components_%d_%g.png", which,noise_rms); 
  cintermediate->SaveAs(str.Data()); 
  /*
   *
  unsigned nint = clean.getIntermediateXcorrs()->size(); 
  cintermediate->Clear(); 
  cintermediate->Divide(2,1); 
  for (unsigned i =  0; i < nint; i++) 
  {
    cintermediate->cd(1) ;
    clean.getIntermediateYs()->at(i)->updateEven()->setPlottingLimits(1.1,true,500); 
    clean.getIntermediateYs()->at(i)->drawEven("alp"); 
    elean.getIntermediateSubs()->at(i)->updateEven()->SetLineColor(2); 
    clean.getIntermediateSubs()->at(i)->updateEven()->SetMarkerColor(2); 
    clean.getIntermediateSubs()->at(i)->drawEven("lp same"); 
    cintermediate->cd(2) ;
    clean.getIntermediateXcorrs()->at(i)->updateEven()->setPlottingLimits(1.1,true,500); 
    clean.getIntermediateXcorrs()->at(i)->drawEven("alp"); 
    str.Form("clean_dbg/iter%d.png", i); 
    cintermediate->SaveAs(str.Data()); 
  }
  */


  FFTtools::saveWisdom("wisdom.fftw3"); 
}
