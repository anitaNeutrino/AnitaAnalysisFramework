
#include "FFTtools.h" 
#include "AnalyticSignal.h" 
#include "AnalysisWaveform.h" 

void testWaveform() 
{
  FFTtools::ThermalNoise noise(5000, 0.2/1.3, 1.2/1.3, 10,2); 
  noise.setExtraNoise(0.5); 

  TGraph * g = noise.makeGaussianDtWf(256,1./2.6,0.1); 

  AnalysisWaveform sig(g->GetN(), g->GetX(), g->GetY(), 1./2.6); 


  AnalysisWaveform::PowerCalculationOptions opt; 
  opt.method = AnalysisWaveform::PowerCalculationOptions::WELCH; 
  opt.window_size=128; 

  TCanvas corig("corig","corig",800,800) ; 
  corig.Divide(2,3); 
  corig.cd(1); 
  g->Draw(); 
  corig.cd(2); 
  sig.drawUneven(); 
  corig.cd(3); 
  sig.drawEven(); 
  corig.cd(4); 
  sig.setPowerCalculationOptions(opt);
  sig.drawPowerdB(); 
  corig.cd(5); 
  sig.drawPhase(); 
  corig.cd(6); 
  sig.drawHilbertEnvelope(); 
  

  TCanvas* ccorr = new TCanvas("ccorr","ccorr",800,800) ; 



  TGraph * g2 = noise.makeGaussianDtWf(256, 1./2.6, 0.1); 
  TGraph * g3 = noise.makeGaussianDtWf(256, 1./2.6, 0.1); 


  for (int i = 0; i < 256; i++) 
  {
    g2->GetY()[i] += 10*sin(2 * TMath::Pi() * (g2->GetX()[i] * 0.1 + 5)); 
    g3->GetY()[i] += 10*sin(2 * TMath::Pi() * (g3->GetX()[i] *0.1 )); 
    g2->GetY()[i] += 10*sin(2 * TMath::Pi() * (g2->GetX()[i] * 0.2 + 5)); 
    g3->GetY()[i] += 10*sin(2 * TMath::Pi() * (g3->GetX()[i]* 0.2 )); 
  }

  double scale = g2->GetRMS(2) * g2->GetRMS(2); 
  AnalysisWaveform *sig2 = new AnalysisWaveform(g2->GetN(), g2->GetX(), g2->GetY(), 1./2.6); 
  AnalysisWaveform *sig3 = new AnalysisWaveform(g3->GetN(), g3->GetX(), g3->GetY(), 1./2.6); 

  sig3->forceEvenSize(sig2->even()->GetN()); 

  ccorr->Divide(2,2); 
  TGraph * fftcorr = FFTtools::getCorrelationGraph((TGraph*) sig2->even(), (TGraph*)sig3->even()); 
  sig2->padEven(1); 
  sig3->padEven(1); 
  
  ccorr->cd(1); 
  sig2->drawEven(); 
  ccorr->cd(2); 
  sig3->drawEven(); 

  AnalysisWaveform * corr = AnalysisWaveform::correlation(sig2,sig3,0,1); 
  ccorr->cd(3); 
  corr->drawEven(); 

  ccorr->cd(4); 
  //AnalysisWaveform *copy = new AnalysisWaveform(*corr2); 

  fftcorr->Draw(); 


  TCanvas * cfreq = new TCanvas("cfreq"); 

  TGraph * g10 = noise.makeGaussianDtWf(256,1./2.6,0.1); 
  AnalysisWaveform * sig10 = new AnalysisWaveform(g10->GetN(), g10->GetX(), g10->GetY(), 1./2.6); 


  int nfreq = sig10->Nfreq(); 
  FFTWComplex * fft = sig10->updateFreq(); 
  for (int i = 0; i < nfreq; i++) 
  {

    if (i > nfreq/4) 
    {
      fft[i].re = 0; 
      fft[i].im = 0; 
    }
  }


  cfreq->Divide(2,1); 
  cfreq->cd(1); 
  sig10->drawPower(); 
  cfreq->cd(2); 
  sig10->drawEven(); 


}
