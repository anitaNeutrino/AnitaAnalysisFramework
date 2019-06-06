#include "FreqDomainFunction.h" 

#include "FFTtools.h" 
#include "FFTWComplex.h" 


FreqDomainFunction::FreqDomainFunction (  std::complex<double> (*fn) (double f, const double * pars), 
        int npars, const AnitaResponse::AbstractResponse * r,  double dt , double tmin, double tmax, bool shift )
: do_shift(shift), t0(tmin), t1(tmax),  the_fn(fn), pars(npars+2), response(r), wf ( 1+(tmax-tmin)/dt, dt, tmin)
{
  wf_init = false; 
  pars[0] = 1; 
  debug = false;
  setParameters(&pars[0]); 
  causal_param = 0; 
}

std::complex<double>* FreqDomainFunction::makeCausal(int N, const std::complex<double> *in, int how, std::complex<double> * out) 
{
  double * real = 0;
  double * imag = 0;
  double * real_prime = 0; 
  double * imag_prime = 0; 

  if (how == 0 || how == 2) 
  {
    real = new double[N]; 
    for (int i = 0; i < N; i++) real[i] = in[i].real(); 
    real_prime = FFTtools::getHilbertTransform(N,real);
  }

  if (how == 1 || how == 2) 
  {
    imag = new double[N]; 
    for (int i = 0; i < N; i++) imag[i] = in[i].imag(); 
    imag_prime = FFTtools::getHilbertTransform(N,imag);
  }

  if (!out) out = new std::complex<double>[N]; 

  for (int i = 0; i < N; i++) 
  {
    if (how == 0) 
    {
      out[i] = std::complex<double>(real[i],real_prime[i]);
    }
    else if (how == 1) 
    {
      out[i] = std::complex<double>(imag_prime[i],imag[i]);
    }

    else if (how == 2) 
    {
      out[i] = std::complex<double>(0.5* (imag_prime[i] + real[i]) , 0.5*(imag[i] + real_prime[i]));
    }
  }

  if (imag) delete[] imag; 
  if (real) delete[] real; 
  if (imag_prime) delete[] imag_prime; 
  if (real_prime) delete[] real_prime; 

  return out; 
}



void FreqDomainFunction::setParameters(const double * p) 
{
  if (p == 0) return; 

  if ( !wf_init || memcmp(p+2, &pars[2], (pars.size()-2) * sizeof(double)))
  {
    if (debug) 
    {
      printf("CHANGING TO: "); 
      for (unsigned i = 0; i < pars.size(); i++) printf(" %g ", p[i]); 
      printf("\n"); 
    }

    memmove(&pars[0],p,pars.size()*sizeof(double)); 
    wf_init = true; 

    int nfreq = wf.Nfreq(); 
    double deltaf = wf.deltaF(); 
    std::complex<double> * Y =  (std::complex<double>*)wf.updateFreq(); 

    for (int i = 0; i < nfreq; i++) 
    {
      double f = i *deltaf; 
      Y[i] = the_fn(f, &pars[2]); 


    }
    //we got to do the kramer's kronig goodness
    if (causal_param) 
    {
        makeCausal(nfreq, Y, causal_param-1, Y); 
    }

    if (response) 
    {
      for (int i = 0; i < nfreq; i++) 
      {
        Y[i] *= std::complex<double> (response->getResponse(i*deltaf)); 
      }
    }

    if (do_shift) FFTtools::rotate(wf.updateEven(), wf.Neven()/2); 
  }
  else
  {
    //set the amp and offset parameters
    pars[0] = p[0]; 
    pars[1] = p[1]; 
  }
}

double FreqDomainFunction::eval(double x, const double * p) 
{
  //see if we need to update wf
  setParameters(p); 
  return pars[0] * wf.evalEven(x-pars[1]) * wf.Nfreq(); 
}

