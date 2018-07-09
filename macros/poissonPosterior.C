//attempt at generating a posterior distribution for the ABCD problem to get a credible interval
// this is purely bayesian with an improper prior, and won't match a frequentist confidence interval. 


const int FLAT_PRIOR = 0; 
const int JEFFREYS_PRIOR = 1; 
const int GAMMA_UNINFORMATIVE_PRIOR = 2; 

const char * priors[] = {"FLAT PRIOR, prior = Gamma(alpha = 1, beta = 0)", "JEFFREYS PRIOR, prior = Gamma(alpha = 0.5, beta = 0)", "GAMMA VAGUE PRIOR, prior = Gamma(alpha = 0, beta = 0)"} ; 

//generates a gamma random variate 
double gammarnd(double alpha, double beta) 
{
  // using the wikipedia algorithm 

  int n = alpha; 
  double d = alpha - n; 

  double sum_integer_part = 0; 

  for (int i = 1; i <=n; i++) sum_integer_part += log(gRandom->Uniform(0,1)); 

  double eta, xi; 
  while (true) 
  {
    double U = gRandom->Uniform(0,1); 
    double V = gRandom->Uniform(0,1); 
    double W = gRandom->Uniform(0,1); 

    if (U < exp(1) / (exp(1) + d)) 
    {
      xi = pow(V,1/d); 
      eta = W * pow(xi,d-1); 
    }
    else
    {
      xi = 1 - log(V); 
      eta = W * exp(-xi); 
    }

    if (eta > pow(xi,d-1) * exp(-xi)) continue; 
    break; 
  }


  return (xi - sum_integer_part) / beta; 
} 




void fillSinglePoisson(int lambda, TH1 * fill,int prior = FLAT_PRIOR , int N = 10000) 
{
  for (int i = 0; i < N; i++) 
  {
    if (prior == FLAT_PRIOR) 
    {
//      fill->Fill(gRandom->Poisson(lambda)); 
      fill->Fill(gammarnd(lambda+1,1)); 
    }
    else if (prior == JEFFREYS_PRIOR) 
    {
      fill->Fill(gammarnd(lambda+0.5,1)); 
    }
    else
    {
      fill->Fill(gammarnd(lambda,1)); 
    }
  }

}
void fillABCD(int nbg, int nsiglike_sideband, int nsideband, TH1 * fill, int prior = FLAT_PRIOR, int N = 10000) 
{
  for (int i = 0; i < N; i++)
  {
    if (prior == FLAT_PRIOR) 
    {
//      fill->Fill(gRandom->Poisson(nbg) * gRandom->Poisson(nsiglike_sideband) / gRandom->Poisson(nsideband)); 
      fill->Fill(gammarnd(nbg+1,1) * gammarnd(nsiglike_sideband+1,1) / gammarnd(nsideband+1,1)); 
    }
    if (prior == JEFFREYS_PRIOR) 
    {
      fill->Fill(gammarnd(nbg+0.5,1) * gammarnd(nsiglike_sideband+0.5,1) / gammarnd(nsideband+0.5,1)); 
    }
    else
    {
      fill->Fill(gammarnd(nbg,1) * gammarnd(nsiglike_sideband,1) / gammarnd(nsideband,1)); 
    }
  }
}


void poissonPosterior(int nbg=1, int nsiglike_sideband = 10, int nsideband =100, int prior = JEFFREYS_PRIOR) 
{
  TH1I fillme("posterior","Posterior test", 100,0,5*nbg*nsiglike_sideband/nsideband); 
  fillABCD(1,10,100,&fillme,prior); 
  fillme.DrawCopy("norml"); 

  const double probsum[]={0.16,0.5,0.84} ; 
  double q[3]; 
  fillme.GetQuantiles(3,q,probsum); 

  printf(" Multiplying posteriors of individual poissons for: \n"); 
  printf("   (nbg=%d) * ( nsiglike_sideband=%d ) / (nsideband =%d)\n", nbg, nsiglike_sideband, nsideband); 
  printf(" %s\n", priors[prior]); 
  printf(" Mean: %g\n",fillme.GetMean()); 
  printf(" RMS: %g\n",fillme.GetRMS()); 
  printf(" Median: %g\n",q[1]); 
  printf(" 16th percentile: %g\n",q[0]); 
  printf(" 84th percentile: %g\n",q[2]); 
  printf(" Median-1 sigma: %g\n",q[1]-q[0]); 
  printf(" Median+1 sigma: %g\n",q[2]-q[1]); 
}