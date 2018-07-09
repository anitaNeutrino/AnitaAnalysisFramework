//attempt at generating a posterior distribution for the ABCD problem to get a credible interval
// this is purely bayesian with an improper prior, and won't match a frequentist confidence interval. 


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


void fillSinglePoisson(int lambda, TH1 * fill, int N = 10000) 
{
  for (int i = 0; i < N; i++) fill->Fill(gammarnd(lambda,1)); 

}
void fillABCD(int nbg, int nsiglike_sideband, int nsideband, TH1 * fill, int N = 100) 
{
  for (int i = 0; i < N; i++) fill->Fill(gammarnd(nbg,1) * gammarnd(nsiglike_sideband,1) / gammarnd(nsideband,1)); 
}


void poissonPosterior(int nbg=1, int nsiglike_sideband = 10, int nsideband =100000) 
{
  TH1I fillme("posterior","Posterior test", 100,0,10); 
  fillABCD(1,10,100,&fillme); 
  fillme.DrawCopy(); 

  const double probsum[]={0.16,0.5,0.84} ; 
  double q[3]; 
  fillme.GetQuantiles(3,q,probsum); 
  printf(" Median: %g\n",q[1]); 
  printf(" -1 sigma: %g\n",q[1]-q[0]); 
  printf(" +1 sigma: %g\n",q[2]-q[1]); 
}
