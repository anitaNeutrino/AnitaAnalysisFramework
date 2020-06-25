#include "Polarimetry.h" 
#include "AnalysisWaveform.h" 

#include "TMultiGraph.h" 
#include "TAxis.h" 
#include "FFTtools.h" 

polarimetry::StokesAnalysis::StokesAnalysis(const StokesAnalysis & other) 
{

  avgI = other.avgI; 
  avgQ = other.avgQ; 
  avgQ = other.avgQ; 
  avgV = other.avgV; 

  dI = new TGraph(*other.dI); 
  dQ = new TGraph(*other.dQ); 
  dU = new TGraph(*other.dU); 
  dV = new TGraph(*other.dV); 

  cI = new TGraph(*other.cI); 
  cQ = new TGraph(*other.cQ); 
  cU = new TGraph(*other.cU); 
  cV = new TGraph(*other.cV); 

 instantaneous.Add(dI,"lp"); 
 instantaneous.Add(dQ,"lp"); 
 instantaneous.Add(dU,"lp"); 
 instantaneous.Add(dV,"lp"); 
 instantaneous.SetTitle("Instantaneous Stokes; Time (ns); Powerish"); 

 cumulative.Add(cI,"lp"); 
 cumulative.Add(cQ,"lp"); 
 cumulative.Add(cU,"lp"); 
 cumulative.Add(cV,"lp"); 
 cumulative.SetTitle("Cumulative Stokes; Time (ns);Powerish"); 

 peak_corr_val = 0; 
 peak_corr_time = 0; 

}


polarimetry::StokesAnalysis::StokesAnalysis(const AnalysisWaveform * H, const AnalysisWaveform *V, double xcorr, double time_around_peak) 
{
  //figure out if we need to resample both to the same time base
  
  AnalysisWaveform *Hre = 0;
  AnalysisWaveform *Vre = 0;


  bool need_to_resample =V->deltaT() != H->deltaT() ||  H->even()->GetX()[0] != V->even()->GetX()[0] ;

  //we need to resample both onto the same base 
  if ( need_to_resample    || xcorr)
  {
    if (need_to_resample) 
    {
      double dt = TMath::Min(H->deltaT(), V->deltaT()); 
      double t0 = TMath::Max(H->even()->GetX()[0], V->even()->GetX()[0]); 
      double t1 = TMath::Min(H->even()->GetX()[H->Neven()-1], V->even()->GetX()[V->Neven()-1]); 

      int N = (t1-t0)/dt; 

      Hre = new AnalysisWaveform(N, t0, dt); 
      Vre = new AnalysisWaveform(N, t0, dt); 

      TGraphAligned * gHre = Hre->updateEven(); 
      TGraphAligned * gVre = Vre->updateEven(); 
        
      H->evalEven(N, gHre->GetX(), gHre->GetY()); 
      V->evalEven(N, gVre->GetX(), gVre->GetY()); 
    }

    if (xcorr) 
    {

      if (!Hre) Hre = new AnalysisWaveform(*H); 
      if (!Vre) Vre = new AnalysisWaveform(*V); 
      if (Hre->Neven()!= Vre->Neven())
      {
        Vre->forceEvenSize(Hre->Neven()); 
      }

      if (!Hre->checkIfPaddedInTime() || !Vre->checkIfPaddedInTime())
      {
        Hre->padEven(1); 
        Vre->padEven(1); 
      }


      AnalysisWaveform * hcorr = Hre;
      AnalysisWaveform * vcorr = Vre;

      // only correlate the parts around the combined peak of the hilbert envelopes
      if (time_around_peak !=0)
      {
        double maxval = 0; 
        double tmax = 0; 
        int imax = 0; 
        double dt = hcorr->deltaT(); 

        for (int i = 0; i < hcorr->Neven(); i++) 
        {
          double val = hcorr->hilbertEnvelope()->GetY()[i] + vcorr->hilbertEnvelope()->Eval(hcorr->hilbertEnvelope()->GetX()[i]);
          if (val > maxval) 
          {
            maxval = val; 
            tmax = hcorr->hilbertEnvelope()->GetX()[i]; 
            imax = i; 
          }
        }

        //make a copy to correlate just around the peak

        int nhalf = ceil(time_around_peak/dt); 
        int max = TMath::Min(imax + nhalf, hcorr->Neven()-1); 
        int min = TMath::Max(imax - nhalf, 0); 

        hcorr = new AnalysisWaveform( max-min+1, hcorr->even()->GetY()+min, dt,hcorr->even()->GetX()[min]);

        dt = vcorr->deltaT(); 
        imax = round(tmax-vcorr->even()->GetX()[0]) / dt; 
        max = TMath::Min(imax + nhalf, vcorr->Neven()-1); 
        min = TMath::Max(imax - nhalf, 0); 

        vcorr = new AnalysisWaveform( max-min+1, vcorr->even()->GetY()+min, dt,vcorr->even()->GetX()[min]);
        if (vcorr->Neven() != hcorr->Neven()) vcorr->forceEvenSize(hcorr->Neven());

      }

      AnalysisWaveform * xc = AnalysisWaveform::correlation(hcorr,vcorr); 
      double dt = xc->deltaT(); 
      int loc = 0; 
      int half = xc->Neven()/2;
      int min = TMath::Max(0, int(half-ceil(xcorr/dt))); 
      int max = TMath::Min(int(half+ceil(xcorr/dt)), xc->Neven()-1); 

      
      xc->even()->peakVal(&loc,min,max,true); 
      peak_corr_time = xc->even()->GetX()[loc]; 
      peak_corr_val = xc->even()->GetY()[loc]; 
      Vre->updateEven()->shift(half-loc,true); 
      delete xc; 
      if (time_around_peak) 
      {
        delete hcorr; 
        delete vcorr;
      }
    }

     H = Hre; 
     V = Vre; 
  }

  int N = TMath::Min(H->Neven(), V->Neven()); 
  

  dI = new TGraph(N, H->even()->GetX(), H->even()->GetY()); 
  dI->SetTitle("Instantaneous I"); 
  dQ = new TGraph(N, H->even()->GetX(), H->even()->GetY()); 
  dQ->SetTitle("Instantaneous Q"); 
  dU = new TGraph(N, H->even()->GetX(), H->even()->GetY()); 
  dU->SetTitle("Instantaneous U"); 
  dV = new TGraph(N, H->even()->GetX(), H->even()->GetY()); 
  dV->SetTitle("Instantaneous V"); 


  FFTtools::stokesParameters(N, H->even()->GetY(), H->hilbertTransform()->even()->GetY(), 
                                V->even()->GetY(), V->hilbertTransform()->even()->GetY(),
                                &avgI,&avgQ,&avgU,&avgV, 
                                dI->GetY(), dQ->GetY(), dU->GetY(), dV->GetY(), 
                                false); 

  cI = new TGraph(dI->GetN(), dI->GetX(), dI->GetY()); 
  cI->SetTitle("Cumulative I"); 
  cQ = new TGraph(dQ->GetN(), dQ->GetX(), dQ->GetY()); 
  cQ->SetTitle("Cumulative Q"); 
  cU = new TGraph(dU->GetN(), dU->GetX(), dU->GetY()); 
  cU->SetTitle("Cumulative U"); 
  cV = new TGraph(dV->GetN(), dV->GetX(), dV->GetY()); 
  cV->SetTitle("Cumulative V"); 


 for (int i = 1; i < N; i++) 
 {
   cI->GetY()[i] += cI->GetY()[i-1]; 
   cQ->GetY()[i] += cQ->GetY()[i-1]; 
   cU->GetY()[i] += cU->GetY()[i-1]; 
   cV->GetY()[i] += cV->GetY()[i-1]; 
 }
  

 instantaneous.Add(dI,"lp"); 
 instantaneous.Add(dQ,"lp"); 
 instantaneous.Add(dU,"lp"); 
 instantaneous.Add(dV,"lp"); 
 instantaneous.SetTitle("Instantaneous Stokes; Time (ns); Powerish"); 

 cumulative.Add(cI,"lp"); 
 cumulative.Add(cQ,"lp"); 
 cumulative.Add(cU,"lp"); 
 cumulative.Add(cV,"lp"); 
 cumulative.SetTitle("Cumulative Stokes; Time (ns);Powerish"); 

 if (Hre) delete Hre; 
 if (Vre) delete Vre; 

}


int polarimetry::StokesAnalysis::computeWindowedAverage(double Ifrac, double *I, double *Q, double  *U, double *V, double *PoPerr) const
{
  int Imax = TMath::LocMax(dI->GetN(), dI->GetY()); 
  double I0 = dI->GetY()[Imax];

  double thresh = dI->GetY()[Imax] * Ifrac; 

  double Isum = 0; 
  double Qsum = 0; 
  double Usum = 0; 
  double Vsum = 0; 

  double I2sum = 0; 
  double Q2sum = 0; 
  double U2sum = 0; 
  double V2sum = 0; 

  int n = 0; 
  int Iend = 0;
  for (int i = Imax; i < dI->GetN(); i++)
  {
    Iend = i;
    if (dI->GetY()[i] < thresh) break; 
    Isum += dI->GetY()[i]; 
    Qsum += dQ->GetY()[i]; 
    Usum += dU->GetY()[i]; 
    Vsum += dV->GetY()[i]; 
    n++;
  }

  for (int i = Imax-1; i >= 0; i--)
  {
    if (dI->GetY()[i] < thresh) break; 
    Isum += dI->GetY()[i]; 
    Qsum += dQ->GetY()[i]; 
    Usum += dU->GetY()[i]; 
    Vsum += dV->GetY()[i]; 
    n++;
  }
  Isum = Isum/n;
  Qsum = Qsum/n;
  Usum = Usum/n;
  Vsum = Vsum/n;


  double ret_err = 0;
  if(PoPerr)
  {
    for (int i = Imax; i < dI->GetN(); i++)
    {
      if (dI->GetY()[i] < thresh) break; 
      I2sum += (dI->GetY()[i]/I0 - Isum/I0)*(dI->GetY()[i]/I0 - Isum/I0); 
      Q2sum += (dQ->GetY()[i]/I0 - Qsum/I0)*(dQ->GetY()[i]/I0 - Qsum/I0); 
      U2sum += (dU->GetY()[i]/I0 - Usum/I0)*(dU->GetY()[i]/I0 - Usum/I0); 
      V2sum += (dV->GetY()[i]/I0 - Vsum/I0)*(dV->GetY()[i]/I0 - Vsum/I0); 
    }

    for (int i = Imax-1; i >= 0; i--)
    {
      if (dI->GetY()[i] < thresh) break; 
      I2sum += (dI->GetY()[i]/I0 - Isum/I0)*(dI->GetY()[i]/I0 - Isum/I0); 
      Q2sum += (dQ->GetY()[i]/I0 - Qsum/I0)*(dQ->GetY()[i]/I0 - Qsum/I0); 
      U2sum += (dU->GetY()[i]/I0 - Usum/I0)*(dU->GetY()[i]/I0 - Usum/I0); 
      V2sum += (dV->GetY()[i]/I0 - Vsum/I0)*(dV->GetY()[i]/I0 - Vsum/I0); 
    }

    I2sum = I2sum/n;
    Q2sum = Q2sum/n;
    U2sum = U2sum/n;
    V2sum = V2sum/n;

    //printf("max ind = %d\n", Imax);
    double Qc = TMath::Abs(cQ->GetY()[Iend])/cI->GetY()[Iend];
    double Uc = TMath::Abs(cU->GetY()[Iend])/cI->GetY()[Iend];

    ret_err = U2sum * pow(.5 * Qc/((Qc * Qc) + (Uc * Uc)), 2) + Q2sum * pow(.5 * Uc/((Qc * Qc) + (Uc * Uc)), 2);
    //printf("errc = %g\n", sqrt(ret_err) * 180./TMath::Pi());
  }

  if (I) *I = Isum; 
  if (Q) *Q = Qsum; 
  if (U) *U = Usum; 
  if (V) *V = Vsum; 

  if (PoPerr) *PoPerr = sqrt(ret_err) * 180./TMath::Pi(); 

  return n; 
}

