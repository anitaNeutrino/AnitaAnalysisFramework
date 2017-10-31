{

  AnitaEventFaker faker("IndividualBRotter");
  UsefulAnitaEvent * ev = faker.makePureNoiseEvent(0.1); 
  faker.addSignal(ev, 0,0); 
  TCanvas * c = new TCanvas; 
  c->Divide(8,6); 

  for (int i =0; i < 48; i++) 
  {
    c->cd(i+1); 
    TGraph * gh = ev->getGraph(i, AnitaPol::kHorizontal);
    gh->SetLineColor(3); 
    gh->SetTitle(TString::Format("ant %d\n",i)); 
    gh->Draw("alp"); 
    TGraph * gv = ev->getGraph(i, AnitaPol::kVertical); 
    gv->SetLineColor(4); 
    gv->Draw("lpsame"); 
  }
}
