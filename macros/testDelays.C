#include "FFTtools.h" 

/* This is a toy sweep of delays to see what it does to instantaneous stokes parameters */ 


int testDelays (int ndelays = 9, double min = -2, double max = 2,double Hamp = 1, double Vamp = -0.5, double xcorr = 0.5)
{


  FFTtools::ButterworthFilter but(FFTtools::BANDPASS, 2, 0.5/1.3,0.3/1.3); 
  TGraph * imp = but.impulseGraph(65,1./2.6, 20); 

//  imp->Draw(); 

  AnalysisWaveform up( imp->GetN(), imp->GetY(), 1./2.6, 0); 
  up.padFreq(10); 

  AnalysisWaveform * H= new AnalysisWaveform(250,0.1,0); 
  TGraphAligned * g = H->updateEven(); 
  up.evalEven(g->GetN(), g->GetX(), g->GetY()); 
  for (int i = 0; i < g->GetN();i++) g->GetY()[i] *= Hamp; 

  TCanvas * c1 = new TCanvas("wf","wf"); 
  TCanvas * c2 = new TCanvas("stokes","stokes"); 

  int nw = ceil(sqrt(ndelays)); 

  c1->Divide(nw,ceil(ndelays/nw)); 
  c2->Divide(nw,ceil(ndelays/nw)); 

  for (int idelay = 0; idelay < ndelays; idelay++) 
  {
    double delay = min + (max - min) / (ndelays - 1) * idelay; 

    //make a copy so that I can change the title 
    TGraph * plotH = new TGraph(H->even()->GetN(), H->even()->GetX(), H->even()->GetY() ); 
    plotH->SetTitle(TString::Format("Delay = %g ns", delay)); 
    plotH->GetXaxis()->SetTitle("ns"); 
    plotH->GetYaxis()->SetTitle("A"); 

    AnalysisWaveform * V= new AnalysisWaveform(250,0.1,0); 
    g = V->updateEven(); 
    for (int i = 0; i < g->GetN(); i++) 
    {
      g->GetY()[i] = Vamp * up.evalEven(g->GetX()[i] - delay) ;
    }

    g->SetLineColor(3); 
    c1->cd(idelay +1); 
    plotH->Draw("al"); 
    V->drawEven("same"); 
    polarimetry::StokesAnalysis*  stokes = new polarimetry::StokesAnalysis(H,V,xcorr,0); 
    c2->cd(idelay+1); 
    double windowed_I, windowed_U, windowed_Q, windowed_V; 
    stokes->computeWindowedAverage(0.25,&windowed_I,&windowed_Q,&windowed_U,&windowed_V); 
    double angle = 90 / TMath::Pi() * atan2(windowed_U,windowed_Q); 
    double linpolfraction = sqrt(windowed_Q*windowed_Q + windowed_U * windowed_U)/ windowed_I; 
    double cpolfraction = abs(windowed_V / windowed_I); 
    printf ("delay: %g, linpolfrac: %g, cpolfrac: %g, angle: %g, xcorr_time:%g\n", delay, linpolfraction, cpolfraction, angle, stokes->getPeakCorrTime()); 

    stokes->instGraphs().SetTitle(TString::Format("%s, angle=%0.2f deg, linfrac =%g, cfrac = %g", plotH->GetTitle(), angle, linpolfraction, cpolfraction)); 
    stokes->instGraphs().Draw("alp pmc plc"); 
    gPad->BuildLegend(0.6,0.6,0.9,0.9,"","lp" ); 
  }







  return 0; 

}
