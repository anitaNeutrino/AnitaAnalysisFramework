#include "SystemResponse.h" 
#include "FFTtools.h" 
#include "TF1.h" 

#include "AnalysisWaveform.h" 


static AnitaResponse::BandLimitedDeconvolution bld(0.18,1.1); 
AnitaResponse::DeconvolutionMethod & AnitaResponse::kDefaultDeconvolution = bld; 

void AnitaResponse::WienerDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  FFTWComplex zero(0,0); 
  for (unsigned i = 0; i < N; i++) 
  {
    double f = i * df; 

    double H2 = response[i].getAbsSq(); 
    double SNR = snr(f,H2,N); 
    printf("SNR(%f) = %f\n", f, SNR); 

    if (SNR <=0)  Y[i] = zero; 
    else
    {
      Y[i] = Y[i] / response[i] * ( H2 / (H2 + 1./SNR)); 
    }
  }
}


AnitaResponse::WienerDeconvolution::WienerDeconvolution(const TGraph *g, const double * sc) 
{
  snr_graph = g; 
  min = g->GetX()[0]; 
  max = g->GetX()[g->GetN()-1]; 
  snr_function = 0; 
  noise_level = 0;
  scale = sc; 
}

AnitaResponse::WienerDeconvolution::WienerDeconvolution(const TF1 *f) 
{
  snr_function = f; 
  noise_level = 0;
  f->GetRange(min,max); 
  snr_graph = 0;
  scale = 0; 
}

AnitaResponse::WienerDeconvolution::WienerDeconvolution(double N) 
{
  snr_graph = 0; 
  scale = 0; 
  snr_function = 0; 
  noise_level = N; 
}



double AnitaResponse::WienerDeconvolution::snr(double f, double R2, int N ) const
{


  if (snr_graph)
  {
    if (f < min || f > max ) return 0; 
    double sc = scale ? *scale : 1; 
    return sc*snr_graph->Eval(f); 
  }
  else if (snr_function) 
  {
//    printf("%f\n",f); 
    if (f < min || f > max ) return 0; 
    return snr_function->Eval(f); 
  }
  else if (noise_level) 
  {
//    printf("%g %g\n",R2,noise_level); 
    return R2/noise_level/N ; 
  }
 
  fprintf(stderr,"Something's wrong...\n"); 

  return -1; 
}

void AnitaResponse::AllPassDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  (void) df; 
  for (unsigned i = 0; i < N; i++) 
  {
    FFTWComplex r = response[i]; 
    r/= r.getAbs(); 
    Y[i]/=r; 
  }
}

void AnitaResponse::ImpulseResponseXCorr::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  (void) df; 
  double total_mag = 0; 
  for (unsigned i = 0; i < N; i++) 
  {
    total_mag += response[i].getAbsSq(); 
  }
  total_mag = sqrt(total_mag); 
  
  for (unsigned i = 0; i < N; i++) 
  {
    FFTWComplex r = response[i].conj(); 
    Y[i]*=r/total_mag; 
  }
}

AnitaResponse::CLEAN::CLEAN() 
{
  fMaxLoops = 1000;
  fLoopGain = 0.02;
}

void AnitaResponse::CLEAN::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{
  
  AnalysisWaveform clean(2*(N-1), Y, df, 0);
  for(int i = 0; i < clean.even()->GetN(); i++)
  {
    clean.updateEven()->GetY()[i] = 0;
  }
  AnalysisWaveform y(2*(N-1), Y, df, 0);
  AnalysisWaveform r(2*(N-1), response, df, 0);
  
  double snr_est = (y.even()->pk2pk())/(y.even()->GetRMS(2)); 
  double thresh = 1./snr_est;
  double loop_gain = fLoopGain;
  int max_loops = fMaxLoops;
  double noise_level = snr_est/3.;
  double E_i = y.even()->getSumV2();
  double E_curr = E_i;

  AnalysisWaveform* xc = 0;
  double max_corr = 1;
  double max_clean = 0;
  int loc = TMath::LocMax(r.Neven(), r.even()->GetY());
  r.updateEven()->shift(r.Neven()/2 - loc, false);
  while(E_curr/E_i > thresh && max_loops != 0)
  {
    max_loops--;
    max_corr = 0;
    double rms1 = y.even()->GetRMS(2);
    double rms2 = r.even()->GetRMS(2);
    xc = AnalysisWaveform::correlation(&y, &r, 0, rms1 * rms2);
    max_corr = TMath::MaxElement(xc->even()->GetN(), xc->even()->GetY());
    double min_corr = TMath::MinElement(xc->even()->GetN(), xc->even()->GetY());
    loc = (max_corr > abs(min_corr)) ? TMath::LocMax(xc->even()->GetN(), xc->even()->GetY()) : TMath::LocMin(xc->even()->GetN(), xc->even()->GetY());
    max_corr = TMath::Max(max_corr, abs(min_corr));
    //printf("adding clean component: %g\n", xc->even()->GetY()[loc] * loop_gain);
    max_clean = TMath::Max(max_clean, max_corr*loop_gain);
    clean.updateEven()->GetY()[loc] += xc->even()->GetY()[loc] * loop_gain;
    r.updateEven()->shift(xc->Neven()/2 - loc, false);
    for(int i = 0; i < y.even()->GetN(); i++)
    {
      y.updateEven()->GetY()[i] -= loop_gain * r.even()->GetY()[i] * xc->even()->GetY()[loc];
    }
    E_curr = y.even()->getSumV2();
    //r.updateEven()->shift(-1*(xc->Neven()/2 - loc), true);
    delete xc;
  }
  //printf("removing clean components below: %g\n", max_clean/noise_level);
  for(int i = 0; i < clean.even()->GetN(); i++)
  {
    clean.updateEven()->GetY()[i] *= (abs(clean.updateEven()->GetY()[i]) > max_clean/noise_level) ? 1 : 0;
    //clean.updateEven()->GetY()[i] += y.even()->GetY()[i];
  }
  for(int i = 0; i < N; i++)
  {
    Y[i] = clean.freq()[i];
  }
  //Y = clean.updateFreq();
}

void AnitaResponse::NaiveDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  (void) df; 
  for (unsigned i = 0; i < N; i++) 
  {
    Y[i]/=response[i]; 
  }
}

void AnitaResponse::BandLimitedDeconvolution::deconvolve(size_t N, double df, FFTWComplex * Y, const FFTWComplex * response) const 
{

  size_t min_i =(size_t) TMath::Max(0,int(min_freq / df)); 
  size_t max_i =(size_t) TMath::Min((int) N-1, int(max_freq / df)); 
  FFTWComplex zero(0,0); 

  for (size_t i = 0; i < min_i; i++) 
  {
    Y[i] = zero; 
  }

  for (size_t i = min_i; i < max_i; i++) 
  {
    Y[i]/=response[i]; 
  }
  for (size_t i = max_i; i < N; i++) 
  {
    Y[i] = zero; 
  }

}



AnitaResponse::Response::Response(int Nfreq, double df, int nangles, const double * the_angles, const FFTWComplex ** the_responses)  
  : Nfreq(Nfreq), df(df)
{
  for (int i = 0; i < nangles; i++) 
  {
    addResponseAtAngle(the_angles[i], the_responses[i]); 
  }

}

void AnitaResponse::Response::addResponseAtAngle(double angle, const FFTWComplex * response) 
{
    FFTWComplex * copy = new FFTWComplex[Nfreq]; 
    memcpy(copy, response, Nfreq * sizeof(FFTWComplex)); 
    responses[angle] = copy; 
//    if (angle == 0) response_at_zero = copy; 
    dirty = true; 
}

AnitaResponse::Response::Response(int Nfreq, double df)
 : Nfreq(Nfreq), df(df) 
{

}

AnitaResponse::Response::Response(const TGraph * timedomain, int npad)
{
  AnalysisWaveform aw(timedomain->GetN(), timedomain->GetX(), timedomain->GetY(), timedomain->GetX()[1] - timedomain->GetX()[0]); 
  aw.padEven(npad); 
  (void)  aw.freq(); 
  df = aw.deltaF(); 
  Nfreq = aw.Nfreq(); 
  addResponseAtAngle(0,aw.freq()); 
}



AnitaResponse::Response::Response(int Nfreq, double df, const FFTWComplex * response)  
  : Nfreq(Nfreq),df(df)
{
  addResponseAtAngle(0, response); 
}


FFTWComplex * AnitaResponse::AbstractResponse::getResponseArray(int N, const double * f, double angle) const
{
  FFTWComplex * answer = new FFTWComplex[N]; 
  for (int i = 0; i < N; i++) answer[i] = getResponse(f[i], angle); 
  return answer; 

}

FFTWComplex * AnitaResponse::AbstractResponse::getResponseArray(int N, double  df, double angle) const
{
  FFTWComplex * answer = new FFTWComplex[N]; 
  for (int i = 0; i < N; i++) answer[i] = getResponse(i*df, angle); 
  return answer; 
}



void AnitaResponse::Response::recompute()  const
{

  int nangles = responses.size(); 
  real.SetBins( Nfreq, 0, df * Nfreq,nangles, -90, 90); 
  imag.SetBins( Nfreq, 0, df * Nfreq, nangles, -90, 90); 
  

 
  if (nangles > 1) 
  {
    double bin_boundaries[nangles+1]; 
    double center_angles[nangles]; 

   //fill in centers of angles
    int i = 0; 
    for (std::map<double, FFTWComplex *>::const_iterator it = responses.begin(); it!=responses.end(); it++)
    {
      center_angles[i++] = it->first; 
    }

    
    bin_boundaries[0] = center_angles[0] - (center_angles[1]-center_angles[0])/2;    
    bin_boundaries[nangles] = center_angles[nangles-1] + (center_angles[nangles-1]-center_angles[nangles-2])/2;    

    for (int i = 1; i < nangles; i++) 
    {
      bin_boundaries[i] = (center_angles[i-1] + center_angles[i])/2; 
    }

    imag.GetYaxis()->Set(nangles, bin_boundaries); 
    real.GetYaxis()->Set(nangles, bin_boundaries); 
  }

  int j = 1; 

  for (std::map<double, FFTWComplex *>::const_iterator it = responses.begin(); it!=responses.end(); it++)
  {
    for (int i = 0; i <  Nfreq; i++) 
    {
      imag.SetBinContent(i+1,j, it->second[i].im); 
      real.SetBinContent(i+1,j, it->second[i].re); 
    }
    j++; 
  }

//  imag.Print(); 

  dirty = false; 
}

FFTWComplex AnitaResponse::CompositeResponse::getResponse(double f, double angle )  const
{

  FFTWComplex answer(1,0); 
  for (size_t i = 0; i < responses.size(); i++) answer*= responses[i]->getResponse(f,angle); 
  return answer; 
}






FFTWComplex AnitaResponse::Response::getResponse(double f, double angle ) const
{
  lock.Lock();
  if (dirty)
  {
    recompute(); 
  }
  lock.UnLock(); 

//  printf("%f %f %f %f %f %f\n", f, angle, real.GetXaxis()->GetXmin(), real.GetXaxis()->GetXmax(), real.GetYaxis()->GetXmin(), real.GetYaxis()->GetXmax()); 
  if ( f > real.GetXaxis()->GetXmax() || f > imag.GetXaxis()->GetXmax()) return FFTWComplex(0,0); 
  double re = real.Interpolate(f,angle); 
  double im = imag.Interpolate(f,angle); 
//  printf("%f %f %f\n", f, re, im); 
  
  return FFTWComplex(re,im); 
}

double AnitaResponse::AbstractResponse::getMagnitude(double f, double angle) const 
{
  return getResponse(f,angle).getAbs(); 
}

double AnitaResponse::AbstractResponse::getPhase(double f, double angle) const 
{
  return getResponse(f,angle).getPhase(); 
}






AnalysisWaveform* AnitaResponse::AbstractResponse::impulseResponse(double dt, int N )  const
{
  AnalysisWaveform * out = new AnalysisWaveform(N, dt); 
  out->updateEven()->GetY()[1] = 1; 
  convolveInPlace(out,0); 
  return out; 
}

void AnitaResponse::AbstractResponse::convolveInPlace(AnalysisWaveform * wf, double angle)  const
{
  int old_size = wf->Neven(); 
//  wf->padEven(2); 
  int nf = wf->Nfreq(); 
  double df = wf->deltaF(); 
  FFTWComplex * fft = wf->updateFreq(); 
  for (int i = 0; i < nf; i++) 
  {
    fft[i] *= getResponse(i*df, angle); 
  }
  wf->updateEven()->Set(old_size); 
}

AnalysisWaveform * AnitaResponse::AbstractResponse::convolve(const AnalysisWaveform * in,  double off_axis_angle ) const
{

  //copy 
  AnalysisWaveform * out = new AnalysisWaveform(*in); 
  convolveInPlace(out, off_axis_angle); 
  return out; 
}

AnalysisWaveform * AnitaResponse::AbstractResponse::deconvolve(const AnalysisWaveform * in,  const DeconvolutionMethod * method, double off_axis_angle) const
{

  //copy 
  AnalysisWaveform * out = new AnalysisWaveform(*in); 
  deconvolveInPlace(out,method,off_axis_angle); 
  return out; 
}


/** caching vars */ 
static __thread const AnitaResponse::AbstractResponse *cache_response = 0; 
static __thread int cache_Nf= 0; 
static __thread double cache_df= 0; 
static __thread double cache_angle= 0; 
static __thread FFTWComplex *cache_V = 0; 


void AnitaResponse::AbstractResponse::deconvolveInPlace(AnalysisWaveform * wf,  const DeconvolutionMethod * method, double off_axis_angle) const
{
//  printf("method: %p\n", method); 
  int old_size = wf->Neven(); 
  wf->padEven(1,0); 
  int nf = wf->Nfreq();
  double df = wf->deltaF(); 
  if (!cache_response || cache_Nf != nf || cache_df != df || !cache_V || cache_angle != off_axis_angle) 
  {
    cache_response = this; 
    cache_Nf = nf; 
    cache_df = df; 
    cache_angle = off_axis_angle; 
    if (cache_V) delete [] cache_V; 
    cache_V = getResponseArray(nf,df,off_axis_angle);
  }
    
  method->deconvolve(nf,df, wf->updateFreq(), cache_V); 
//  wf->updateEven()->Set(old_size); 

}
