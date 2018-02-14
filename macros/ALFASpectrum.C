// This must run with AnitaAnalysisFramework loaded
void ALFASpectrum(int eventNumber)
{

  int run = AnitaDataset::getRunContainingEventNumber(eventNumber); 

  AnitaDataset * d = new AnitaDataset(run); 
  d->getEvent(eventNumber); 

  d->useful()->setAlfaFilterFlag(false); 

  TGraph * alfaGraph = d->useful()->getGraph(AnitaRing::kTopRing, 4, AnitaPol::kHorizontal); 

//  printf("%p\n", alfaGraph); 
  AnalysisWaveform *wf_akima = new AnalysisWaveform( alfaGraph->GetN(), alfaGraph->GetX(), alfaGraph->GetY(), 1./2.6, AnalysisWaveform::AKIMA);

//  AnalysisWaveform *wf_sparse = new AnalysisWaveform( alfaGraph->GetN(), alfaGraph->GetX(), alfaGraph->GetY(), 1./2.6, AnalysisWaveform::REGULARIZED_SPARSE_YEN);

  const TGraphAligned * uneven = wf_akima->uneven(); 

  TString fname;
  fname.Form("alfa_%d_uneven.csv",eventNumber); 

  FILE * ff = fopen(fname.Data(),"w"); 
  for (int i = 0; i < uneven->GetN(); i++)
    fprintf(ff,"%g,%g\n", uneven->GetX()[i], uneven->GetY()[i]); 
  fclose(ff); 

  fname.Form("alfa_%d_akima.csv",eventNumber); 

  const TGraphAligned * even = wf_akima->even(); 

  ff = fopen(fname.Data(),"w"); 
  for (int i = 0; i < even->GetN(); i++)
    fprintf(ff,"%g,%g\n", even->GetX()[i], even->GetY()[i]); 
  fclose(ff); 

  /*
  new TCanvas; 
  alfaGraph->Draw("al"); 
  wf_akima->setTitle("akima"); 
  wf_sparse->setTitle("sparse"); 
  wf_akima->setColor(2); 
  wf_sparse->setColor(3); 
//  wf_akima->padFreq(3); 
//  wf_sparse->padFreq(3); 
  wf_akima->drawEven("lsame"); 
  wf_sparse->drawEven("lsame"); 

  */

  TCanvas * c = new TCanvas;
  wf_akima->padEven(10); 
 //wf_sparse->padEven(9); 
//  wf_sparse->setFreqDisplayRange(0.7,0.9);
  wf_akima->drawPower("al"); 
  wf_akima->setFreqDisplayRange(0.7,0.9);
  wf_akima->setTitle(TString::Format("%d",eventNumber)); 

  TGraph *g = new TGraph(2); 
  g->GetX()[0] = 0.792; 
  g->GetX()[1] = 0.792; 
  g->GetY()[1] = 1; 
  g->GetY()[0] = 0; 


//  wf_sparse->drawPower("lsame"); 

  g->Draw("lsame"); 

  c->SaveAs(TString::Format("alfa_%d.png", eventNumber)); 


}

void doSW()
{

  ALFASpectrum(9097052);
  ALFASpectrum(11116678);
  ALFASpectrum(11989351);
  ALFASpectrum(15717143);
  ALFASpectrum(16952278);
  ALFASpectrum(19459892);
  ALFASpectrum(23695250);
  ALFASpectrum(27142556);
  ALFASpectrum(32907849);
  ALFASpectrum(33485018);
  ALFASpectrum(39599214);
  ALFASpectrum(41529159);
  ALFASpectrum(58592888);
  ALFASpectrum(62273737);
  ALFASpectrum(66313840);
  ALFASpectrum(68298850);
  ALFASpectrum(70013895);
  ALFASpectrum(73726756);
  ALFASpectrum(75277753);
  ALFASpectrum(83878009);

}

void doAll() 
{

  ALFASpectrum(9097075); 
  ALFASpectrum(11116669); 
  ALFASpectrum(11989349); 
  ALFASpectrum(15717147); 
  ALFASpectrum(16952229); 
  ALFASpectrum(19459851); 
  ALFASpectrum(23695286); 
  ALFASpectrum(27142546); 
  ALFASpectrum(32907848); 
  ALFASpectrum(33484995); 
  ALFASpectrum(39599205); 
  ALFASpectrum(49599205); 
  ALFASpectrum(41529195); 
  ALFASpectrum(58592863); 
  ALFASpectrum(62273732); 
  ALFASpectrum(66313844); 
  ALFASpectrum(68298837); 
  ALFASpectrum(70013898); 
  ALFASpectrum(73726742); 
  ALFASpectrum(75277769); 
  ALFASpectrum(83877990); 

}
