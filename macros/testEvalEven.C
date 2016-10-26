#include "AnalysisWaveform.h" 
void testEvalEven()
{
  AnalysisWaveform wf; 
  TGraph * g = wf.updateEven();


  for (int i = 0; i < g->GetN(); i++) 
  {
    g->GetY()[i] = sin(g->GetX()[i] * 100-0.1); 
  }

  TGraph *g2  = new TGraph(100); 

  for (int i = 0; i < 100; i++)
  {
    g2->GetX()[i] = 0.9 * i; 
  }
  wf.evalEven(g2->GetN(), g2->GetX(), g2->GetY());

  g2->Draw("alp");
}
