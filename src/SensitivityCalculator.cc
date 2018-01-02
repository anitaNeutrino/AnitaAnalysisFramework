#include "SensitivityCalculator.h" 
#include "TGraphErrors.h" 
#include "TAxis.h" 
#include "TH1.h" 

#ifdef USE_ROOSTATS
#include "RooWorkspace.h" 
#include "RooStats/FeldmanCousins.h" 
#include "RooStats/HypoTestResult.h" 
#include "RooStats/ProfileLikelihoodCalculator.h" 
#endif


SensitivityCalculator::~SensitivityCalculator() 
{
#ifdef USE_ROOSTATS

  delete w; 
#endif
}



SensitivityCalculator::SensitivityCalculator(double CL, bool FC) 
{
  cl = CL;
  use_fc = FC; 

#ifndef USE_ROOSTATS

  fprintf(stderr,"Need to build with USE_ROOSTATS to use Sensitivity Calculator\n"); 


#else
  w = new RooWorkspace("w"); 

  
  /* We observe from a poisson distribution, with a mean of e * s + banthro + bthermal */ 
  w->factory("Poisson::obsModel(nobs[0,0,20], sum(prod(s[1,0,50], eff[0.7,0,1]),banthro[1,0,30],bthermal[0.3,0,30]))"); 

  /*We estimate bthermal from the counts in a sideband that's tauthermal size the signal region, so nthermal (what we observe) is poisson distributed with a mean of tauthermal * bthermal */ 
  w->factory("Poisson::thermalModel(nthermal[1,0,30], prod::taubthermal(bthermal, tauthermal[3,0,10]))");

  /*We estimate banthro from the counts in a sideband that's tauanthro size the signal region, so nanthro (what we observe) is poisson distributed with a mean of tauanthro * banthro */ 
  w->factory("Poisson::anthroModel(nanthro[1,0,30], prod::taubanthro(banthro, tauanthro[1,0,10]))");

  /*Errors on tau anthro */ 
  w->factory("Gaussian::anthroTauErrorModel(mean_tauanthro[1,0,10],tauanthro, sigma_tau_anthro[0.1,0,1])"); 

  /* We assign a systematic error (gaussian) to our analysis efficiency */ 
  w->factory("Gaussian::effModel(mean_eff[0.7,0,1], eff, sigma_eff[0.05,0,0.2])"); 

  /* The likelihood  model is the product of all of the distributions */ 
  w->factory("PROD::model(obsModel,thermalModel,anthroModel,effModel, anthroTauErrorModel)"); 

  w->factory("SUM:total_bg(banthro,bthermal)"); 




  /* Set all auxilliary observables to constants */ 
  w->var("tauthermal")->setConstant(); 
  w->var("mean_tauanthro")->setConstant(); 
  w->var("sigma_tau_anthro")->setConstant(); 
  w->var("nthermal")->setConstant(); 
  w->var("nanthro")->setConstant(); 
  w->var("mean_eff")->setConstant(); 
  w->var("sigma_eff")->setConstant(); 

  w->defineSet("poi","s"); //parameter of interest
  w->defineSet("nuis","bthermal,banthro,eff,tauanthro"); //nuisance parameters 
  w->defineSet("gobs","nthermal,nanthro,tauthermal,mean_tauanthro,sigma_tau_anthro,mean_eff,sigma_eff"); //"global observed parameters" 
  w->defineSet("obs","nobs"); //"main observed parameters" 
  w->defineSet("sens","bthermal,banthro,eff,nobs,tauanthro"); 
  w->defineSet("bkg_anthro","banthro,tauanthro"); 
  w->defineSet("bkg_thermal","bthermal"); 
  w->defineSet("bkg_total","banthro,tauanthro,bthermal"); 

  setEfficiency(0.7,0.05); 
  setThermalBackground(1,3); 
  setAnthroBackground(1,1,0.1); 
#endif

}  

void SensitivityCalculator::SensitivityCalculator::setThermalBackground(int nobs, double tau) 
{
#ifdef USE_ROOSTATS
  w->var("tauthermal")->setVal(tau); 
  w->var("nthermal")->setVal(nobs); 
#endif
}


void SensitivityCalculator::SensitivityCalculator::setAnthroBackground(int nobs, double tau, double sigma_tau)
{
#ifdef USE_ROOSTATS
  w->var("mean_tauanthro")->setVal(tau); 
  w->var("sigma_tau_anthro")->setVal(sigma_tau); 
  w->var("nanthro")->setVal(nobs); 
#endif

}

void SensitivityCalculator::SensitivityCalculator::setEfficiency(double mean, double sigma)
{
#ifdef USE_ROOSTATS
  if (sigma == 0) 
  {
    w->var("eff")->setVal(mean); 
    w->var("eff")->setConstant(true); 

  }
  else
  {
    w->var("eff")->setConstant(false); 
    w->var("mean_eff")->setVal(mean); 
    w->var("sigma_eff")->setVal(sigma); 
  }
#endif
}



double SensitivityCalculator::getLimit(int nobs, double * lower, double * upper) 
{
#ifdef USE_ROOSTATS
  w->var("s")->setConstant(false); 
  w->var("nobs")->setVal(nobs); 
  RooDataSet data("data","",*w->var("nobs")); 
  data.add(*w->var("nobs")); 

  //generate model with what we've observed 
  RooStats::ModelConfig m("mod"); 
  m.SetWorkspace(*w); 
  m.SetPdf(*w->pdf("model")); 
  m.SetObservables(*w->var("nobs")); 
  m.SetGlobalObservables(*w->set("gobs")); 
  m.SetParametersOfInterest(*w->set("poi")); 
  m.SetNuisanceParameters(*w->set("nuis")); 

  if (use_fc)
  {
    RooStats::FeldmanCousins fc(data,m); 
    fc.SetConfidenceLevel(cl); 
    fc.FluctuateNumDataEntries(false); // number counting experiment! 
    fc.UseAdaptiveSampling(true); //no idea what this does but it's supposed to make it faster
    fc.SetNBins(501);
    RooStats::PointSetInterval * ci = fc.GetInterval();  

    if (upper) *upper=ci->UpperLimit(*w->var("s")); 
    if (lower) *lower=ci->LowerLimit(*w->var("s")); 
    return 0;
  }
  else
  {
    //use profile likelihood
    RooStats::ProfileLikelihoodCalculator pl(data, m); 
    pl.SetConfidenceLevel(cl); 
    pl.SetNullParameters(*w->set("poi")); 
    RooStats::LikelihoodInterval * ci = pl.GetInterval(); 
    pl.GetHypoTest()->Print(); 
    if (upper) *upper=ci->UpperLimit(*w->var("s")); 
    if (lower) *lower=ci->LowerLimit(*w->var("s")); 
    return 0; 
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


TH1 * SensitivityCalculator::histBAnthro() 
{
#ifdef USE_ROOSTATS
  RooDataSet * data = w->pdf("model")->generate(*w->set("bkg_anthro"),10000); 
  TH1 * hist =  data->createHistogram("B_{anthro}", *w->var("banthro")); 
  hist->Scale(1./10000); 
  delete data; 
  return hist; 
#else 
  return 0; 
#endif
}

TH1 * SensitivityCalculator::histBThermal() 
{
#ifdef USE_ROOSTATS
  RooDataSet * data = w->pdf("model")->generate(*w->set("bkg_thermal"),10000); 
  TH1 * hist =  data->createHistogram("B_{thermal}", *w->var("bthermal")); 
  hist->Scale(1./10000); 
  delete data; 
  return hist; 
#else 
  return 0; 
#endif 
}

TH1 * SensitivityCalculator::histBTotal() 
{
  TH1 * banthro = histBAnthro(); 
  TH1 * bthermal = histBThermal(); 
  TH1 *  btotal = new TH1D("total", "Total", 100,0,20); 

  for (int i = 0; i < 10000; i++)
  {
    btotal->Fill(banthro->GetRandom() + bthermal->GetRandom()); 
  }
  btotal->Scale(1./btotal->Integral()); 
  delete banthro; 
  delete bthermal; 
  return btotal; 
}




TH1 * SensitivityCalculator::histNObserved(double S)
{
#ifdef USE_ROOSTATS

  w->var("s")->setVal(S); 
  w->var("s")->setConstant(true); 
  RooDataSet * data = w->pdf("model")->generate(*w->set("sens"),10000); 
//  data->Print(); 
  TString str;
  str.Form("observed_s_%g",S); 
  TH1 * hist =  data->createHistogram(str.Data(), *w->var("nobs"),RooFit::Binning(20,0,20)); 
  hist->Scale(1./10000); 
  delete data; 
  w->var("s")->setConstant(false);
  return hist; 
#else
  return 0;
#endif
}




