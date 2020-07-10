#include "FFTtools.h" 





TH2D * hpower =0;
TH2D * hphase =0;
std::vector<TGraph*> vpower; 
std::vector<TGraph*> vphase; 

void setupHists()
{

  TFile f("data/templates/crTmpltsA3.root"); 

  for (int i = 0; i < 27; i++) 
  {
    TGraph * g = FFTtools::cropWave((TGraph*) f.Get(Form("efield%d",i)),0,600); 
    g->SetTitle(Form("Template %d",i)); 

    AnalysisWaveform * wf = new AnalysisWaveform(g->GetN(),g->GetY(), g->GetX()[1]-g->GetX()[0], g->GetX()[0]); 
    wf->setTitle(g->GetTitle()); 

    int roll; 
    wf->hilbertEnvelope()->peakVal(&roll);  
    FFTtools::rotate(wf->updateEven(), -roll); 

    TGraphAligned * gpower = (TGraphAligned*) wf->power(); 
    TGraphAligned * gphase = new TGraphAligned(*wf->phase()); 

    //find the peak of the power spectrum 
    int peak; 
    double val = gpower->peakVal(&peak); 
    FFTtools::unwrap(gphase->GetN(), gphase->GetY(), 2*TMath::Pi()); 
    for (int j = 0; j < gpower->GetN(); j++) 
    {
      if (j < peak) continue; 

      if (gpower->GetY()[j] / val < 1e-3) 
      {
        gphase->Set(j); 
        break; 
      }
    }

    gpower->SetTitle(g->GetTitle());
    gphase->SetTitle(g->GetTitle());
    vpower.push_back(gpower);
    vphase.push_back(gphase); 
  }

  double df = vpower[0]->GetX()[1] - vpower[0]->GetX()[0];

  hpower = new TH2D("hpower","POWER;freq;T", vpower[0]->GetN(), -df/2, vpower[0]->GetX()[vpower[0]->GetN()-1]+df/2,
      27,-13.5,13.5); 
  hphase = new TH2D("hphase","PHASE;freq;T", vpower[0]->GetN(), -df/2, vpower[0]->GetX()[vpower[0]->GetN()-1]+df/2,
      27,-13.5,13.5); 

  hpower->SetDirectory(0); 
  hphase->SetDirectory(0); 

  for (int x = 1; x <= hpower->GetNbinsX(); x++) 
  {
    for (int y = 1; y <= hpower->GetNbinsY(); y++)
    {
       hpower->SetBinContent(x,y, vpower[y-1]->GetY()[x-1]); 
       hphase->SetBinContent(x,y, vphase[y-1]->GetN() < x ? vphase[y-1]->GetY()[vphase[y-1]->GetN()-1] : vphase[y-1]->GetY()[x-1]); 
    }
  }

}

std::complex<double> crfitfn (double f, const double * pars) 
{

  if (!hpower) setupHists(); 

  double T = *pars; 
  if (T < -13 || T > 13) return 0; 
  if (f < 0 || f> hpower->GetXaxis()->GetBinCenter(hpower->GetNbinsX())) return 0; 

//  printf(" %g %g\n",f,T); 
  double pow = hpower->Interpolate(f,T); 
  double phase = hphase->Interpolate(f,T); 
  double amp = sqrt(pow); 

  return std::complex<double>( cos(phase) * amp, sin(phase)*amp); 
}

void makeSeq(bool with_response = false) 
{

  FreqDomainFunction fn(&crfitfn,1); 

  if (with_response) 
  {
    AnitaResponse::ResponseManager * rm = new AnitaResponse::ResponseManager("SingleBRotter",3);
    fn.setResponse(rm->response(0,0)); 
  }


  TF1 * f = fn.makeTF1("cr"); 
  f->SetParameter(0,1); 
  f->SetParameter(1,0); 

  const char * outfile = with_response ? "cr_response.gif" : "cr.gif"; 
  gSystem->Unlink(outfile); 

  TCanvas * c = new TCanvas("c","cr",800,600); 
  for (double T = -5 ; T <=5; T+=0.1)
  {
    f->SetParameter(2,T); 
    f->SetTitle(Form("T=%g",T)); 
    f->Draw(); 

    c->Print(with_response ? "cr_response.gif+" : "cr.gif+"); 
  }
  c->Print(with_response ? "cr_response.gif++" : "cr.gif++"); 
}



void fitTemplates(int event, int anita = 3) 
{

  AnitaDataset d(150); 
  d.getEvent(event); 

  FilteredAnitaEvent fev





}




  
