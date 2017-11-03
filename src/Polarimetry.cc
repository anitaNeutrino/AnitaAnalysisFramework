#include "Polarimetry.h" 
#include "AnalysisWaveform.h" 

#include "TMultiGraph.h" 
#include "FFTtools.h" 



polarimetry::StokesAnalysis::StokesAnalysis(const AnalysisWaveform * H, const AnalysisWaveform *V) 
{
  //figure out if we need to resample both to the same time base
  
  AnalysisWaveform *Hre = 0;
  AnalysisWaveform *Vre = 0;


  //we need to resample both onto the same base 
  if (V->deltaT() != H->deltaT()
      ||  H->even()->GetX()[0] != V->even()->GetX()[0] 
      ||  H->Neven() != V->Neven() )
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

    H = Hre; 
    V = Vre; 
  }

  

  dI = new TGraph(H->Neven(), H->even()->GetX(), H->even()->GetY()); 
  dI->SetTitle("Instantaneous I"); 
  dQ = new TGraph(H->Neven(), H->even()->GetX(), H->even()->GetY()); 
  dQ->SetTitle("Instantaneous Q"); 
  dU = new TGraph(H->Neven(), H->even()->GetX(), H->even()->GetY()); 
  dU->SetTitle("Instantaneous U"); 
  dV = new TGraph(H->Neven(), H->even()->GetX(), H->even()->GetY()); 
  dV->SetTitle("Instantaneous V"); 

  int N = H->Neven(); 

  FFTtools::stokesParameters(N, H->even()->GetY(), H->hilbertTransform()->even()->GetY(), 
                                V->even()->GetY(), V->hilbertTransform()->even()->GetY(),
                                0,0,0,0, 
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
  

 instantaneous.Add(dI); 
 instantaneous.Add(dQ); 
 instantaneous.Add(dU); 
 instantaneous.Add(dV); 

 cumulative.Add(cI); 
 cumulative.Add(cQ); 
 cumulative.Add(cU); 
 cumulative.Add(cV); 

 if (Hre) delete Hre; 
 if (Vre) delete Vre; 

}


int polarimetry::StokesAnalysis::computeWindowedAverage(double Ifrac, double *I, double *Q, double  *U, double *V) const
{
  int Imax = TMath::LocMax(dI->GetN(), dI->GetY()); 

  double thresh = Imax * Ifrac; 

  double Isum = 0; 
  double Qsum = 0; 
  double Usum = 0; 
  double Vsum = 0; 

  int n = 0; 
  for (int i = Imax; i < dI->GetN(); i++)
  {
    if (dI->GetY()[i] < thresh) break; 
    Isum += dI->GetY()[i]; 
    Qsum += dI->GetY()[i]; 
    Usum += dI->GetY()[i]; 
    Vsum += dI->GetY()[i]; 
    n++;
  }

  for (int i = Imax-1; i >= 0; i--)
  {
    if (dI->GetY()[i] < thresh) break; 
    Isum += dI->GetY()[i]; 
    Qsum += dI->GetY()[i]; 
    Usum += dI->GetY()[i]; 
    Vsum += dI->GetY()[i]; 
    n++;
  }

  if (I) *I = Isum/n; 
  if (Q) *Q = Qsum/n; 
  if (U) *U = Usum/n; 
  if (V) *V = Vsum/n; 

  return n; 
}





