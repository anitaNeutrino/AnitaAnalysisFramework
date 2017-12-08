#include "SensitivityCalculator.h" 
#include "TGraphErrors.h" 
#include "TAxis.h" 

#ifdef USE_ROOSTATS
#include "RooStats/FeldmanCousins.h" 
#include "RooStats/ProfileLikelihoodCalculator.h" 
#endif



SensitivityCalculator::SensitivityCalculator(double CL, bool FC) 
#ifdef USE_ROOSTATS
: w("w")
#endif 
{
  cl = CL;
  use_fc = FC; 

#ifndef USE_ROOSTATS

  fprintf(stderr,"Need to build with USE_ROOSTATS to use Sensitivity Calculator\n"); 


#else

  w.factory("Poisson::obsModel(nobs[0,0,50], sum(prod(s[1,0,50], eff[0.7,0,1]),banthro[1,0,50],bthermal[1,0,50]))"); 
  w.factory("Poisson::thermalModel(nthermal[1,0,10], prod::taubthermal(bthermal, tauthermal[3,0,10]))");
  w.factory("Poisson::anthroModel(nanthro[1,0,10], prod::taubanthro(banthro, tauanthro[2,0,10]))");
  w.factory("Gaussian::effModel(mean_eff[0.7,0,1], eff, sigma_eff[0.05,0,0.2])"); 

  w.factory("PROD::model(obsModel,thermalModel,anthroModel,effModel)"); 

  w.var("tauthermal")->setConstant(); 
  w.var("tauanthro")->setConstant(); 
  w.var("nthermal")->setConstant(); 
  w.var("nanthro")->setConstant(); 
  w.var("mean_eff")->setConstant(); 
  w.var("sigma_eff")->setConstant(); 

  w.defineSet("poi","s"); //parameter of interest
  w.defineSet("nuis","bthermal,banthro,eff"); //nuisance parameters 
  w.defineSet("gobs","nthermal,nanthro,tauthermal,tauanthro,mean_eff,sigma_eff"); //"global observed parameters" 
  w.defineSet("obs","nobs"); //"global observed parameters" 

  setEfficiency(0.7,0.05); 
  setThermalBackground(1,3); 
  setAnthroBackground(0,0.3); 
#endif

}  

void SensitivityCalculator::SensitivityCalculator::setThermalBackground(int nobs, double tau)
{
#ifdef USE_ROOSTATS

  w.var("tauthermal")->setVal(tau); 
  w.var("nthermal")->setVal(nobs); 
#endif
}


void SensitivityCalculator::SensitivityCalculator::setAnthroBackground(int nobs, double tau)
{
#ifdef USE_ROOSTATS
  w.var("tauanthro")->setVal(tau); 
  w.var("nanthro")->setVal(nobs); 
#endif

}

void SensitivityCalculator::SensitivityCalculator::setEfficiency(double mean, double sigma)
{
#ifdef USE_ROOSTATS
  w.var("mean_eff")->setVal(mean); 
  w.var("sigma_eff")->setVal(sigma); 
#endif
}


void SensitivityCalculator::getLimit(int nobs, double * lower, double * upper) 
{
#ifdef USE_ROOSTATS
  w.var("nobs")->setVal(nobs); 
  RooDataSet data("data","",*w.var("nobs")); 
  data.add(*w.var("nobs")); 

  //generate model with what we've observed 
  RooStats::ModelConfig m("limits"); 
  m.SetWorkspace(w); 
  m.SetPdf(*w.pdf("model")); 
  m.SetObservables(*w.var("nobs")); 
  m.SetGlobalObservables(*w.set("gobs")); 
  m.SetParametersOfInterest(*w.set("poi")); 
  m.SetNuisanceParameters(*w.set("nuis")); 

  if (use_fc)
  {
    RooStats::FeldmanCousins fc(data,m); 
    fc.SetConfidenceLevel(cl); 
    fc.FluctuateNumDataEntries(false); // number counting experiment! 
    fc.UseAdaptiveSampling(true); //no idea what this does but it's supposed to make it faster
    fc.SetNBins(501);
    RooStats::PointSetInterval * ci = fc.GetInterval();  

    if (upper) *upper=ci->UpperLimit(*w.var("s")); 
    if (lower) *lower=ci->LowerLimit(*w.var("s")); 
  }
  else
  {
    //use profile likelihood
    RooStats::ProfileLikelihoodCalculator pl(data, m); 
    pl.SetConfidenceLevel(cl); 
    RooStats::LikelihoodInterval * ci = pl.GetInterval(); 
    if (upper) *upper=ci->UpperLimit(*w.var("s")); 
    if (lower) *lower=ci->LowerLimit(*w.var("s")); 
  }
#endif


}


TGraphErrors * SensitivityCalculator::confidenceBands(int start, int stop) 
{
  TGraphErrors * g = new TGraphErrors(stop-start+1); 
  g->SetMarkerStyle(0); 
  g->SetTitle(TString::Format("%g%% bands", cl * 100)); 
  g->GetXaxis()->SetTitle("Nobserved"); 
  g->GetYaxis()->SetTitle("Signal"); 

  int ii = 0;
  for (int i = start; i <= stop; i++)
  {
    double l,u; 
    getLimit(i,&l,&u); 

    g->SetPoint(ii, i, (l+u)/2); 
    g->SetPointError(ii, 0, (u-l)/2); 
    ii++; 
  }

  return g; 

}






