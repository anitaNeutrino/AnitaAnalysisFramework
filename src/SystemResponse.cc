#include "SystemResponse.h" 
#include "FFTtools.h" 
#include "TF1.h" 

#include "AnalysisWaveform.h" 


static AnitaResponse::BandLimitedDeconvolution bld(0.18,1.1); 
AnitaResponse::DeconvolutionMethod & AnitaResponse::kDefaultDeconvolution = bld; 


void AnitaResponse::DeconvolutionMethod::deconvolve(int Nf, double df, FFTWComplex * Y, const FFTWComplex * response) const
{
  AnalysisWaveform wf(2*(Nf-1),Y, df,0);
  AnalysisWaveform r(2*(Nf-1),response, df,0);
  deconvolve(&wf,&r); 

}

void AnitaResponse::WienerDeconvolution::deconvolve(AnalysisWaveform *wf, const AnalysisWaveform * response_wf) const 
{

  size_t N = wf->Nfreq(); 
  double df = wf->deltaF(); 
  FFTWComplex * Y = wf->updateFreq(); 

  const FFTWComplex * response = response_wf->freq(); 


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

void AnitaResponse::AllPassDeconvolution::deconvolve(AnalysisWaveform * wf, const AnalysisWaveform * response_wf) const 
{

  size_t N = wf->Nfreq(); 
  FFTWComplex * Y = wf->updateFreq(); 
  const FFTWComplex * response = response_wf->freq();

  for (unsigned i = 0; i < N; i++) 
  {
    FFTWComplex r = response[i]; 
    r/= r.getAbs(); 
    Y[i]/=r; 
  }
}

void AnitaResponse::ImpulseResponseXCorr::deconvolve(AnalysisWaveform * wf, const AnalysisWaveform * response_wf) const 
{

  size_t N = wf->Nfreq(); 
  FFTWComplex * Y = wf->updateFreq(); 
  const FFTWComplex * response = response_wf->freq();
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

AnitaResponse::CLEAN::CLEAN(int max_loops, double loop_gain, double thresh_factor, TString restoring_beam, bool add_residuals, bool only_return_residuals)
{
  fMaxLoops = max_loops;
  fLoopGain = loop_gain;
  fThreshFactor = thresh_factor;
  fRestoringBeam = restoring_beam;
  fAddResiduals = add_residuals;
  fOnlyReturnResiduals = only_return_residuals;
}

void AnitaResponse::CLEAN::deconvolve(AnalysisWaveform * Y, const AnalysisWaveform * resp) const 
{
  //set up analysis waveforms
  AnalysisWaveform clean(*Y);
  for(int i = 0; i < clean.even()->GetN(); i++)
  {
    clean.updateEven()->GetY()[i] = 0;
  }
  AnalysisWaveform y(*Y);
  AnalysisWaveform rbeam(*Y);
  AnalysisWaveform r(*resp);
  AnalysisWaveform* xc = 0;

  //set up hyperparameters
  double noise_est=0;
  double noise_mean=0;
  double n_noise=0;
  for(int i = 0; i < y.Neven(); i++)
  {
    if(y.even()->GetX()[i] < 25 || y.even()->GetX()[i] > 70)
    {
      noise_est += y.even()->GetY()[i] * y.even()->GetY()[i];
      noise_mean +=  y.even()->GetY()[i];
      n_noise++;
    }
  }
  noise_mean /= n_noise;
  noise_est = sqrt(noise_est/n_noise - (noise_mean * noise_mean));
  //double noise_est = TMath::RMS(y.Neven(), y.even()->GetY());
  double snr_est = (1./fThreshFactor) * (y.even()->pk2pk())/(2*noise_est); 
  //printf("snr est = %g, rms = %g, nt1 = %d, vpp = %g\n", snr_est, noise_est, n_noise, y.even()->pk2pk());

  double thresh = 1./snr_est;
  double loop_gain = fLoopGain;
  int max_loops = fMaxLoops;
  double noise_level = snr_est/3.;
  double max_corr = 0;
  double max_clean = 0;
  
  //rescale template to WF
  double max_y = y.even()->peakVal(0, 0, -1, true);
  double max_r = r.even()->peakVal(0, 0, -1, true);

  double scale_to_wf = (TMath::Abs(max_y))/(TMath::Abs(max_r));

  for(int i = 0; i < r.Neven(); i++)
  {
    r.updateEven()->GetY()[i] *= scale_to_wf;
  }

  //rescale
  double normalize = 1./sqrt(y.even()->getSumV2() + r.even()->getSumV2());
  for(int i = 0; i < y.Neven(); i++)
  {
    y.updateEven()->GetY()[i] *= normalize;
  }
  for(int i = 0; i < r.Neven(); i++)
  {
    r.updateEven()->GetY()[i] *= normalize;
  }
  if(fRestoringBeam.Length() > 0)
  {
    TF1 restore("restoring_beam", fRestoringBeam.Data(), y.even()->GetX()[0], y.even()->GetX()[y.Neven()-1]); 
    for(int i = 0; i < rbeam.Neven(); i++)
    {
      rbeam.updateEven()->GetY()[i] = restore.Eval(rbeam.updateEven()->GetX()[i]);
    }

    //normalize to unit sum
    double restore_norm = rbeam.even()->getSumV2();
    for(int i = 0; i < rbeam.Neven(); i++)
    {
      rbeam.updateEven()->GetY()[i] /= restore_norm;
    }
  }

  //shift response to middle for ease
  int loc = TMath::LocMax(r.Neven(), r.even()->GetY());
  r.updateEven()->shift(r.Neven()/2 - loc, false);
  double E_i = y.even()->getSumV2();
  double E_curr = E_i;

  //start CLEANING
  while(((E_curr/E_i) > thresh) && (max_loops != 0))
  {
    max_loops--;
    double rms1 = y.even()->GetRMS(2);
    double rms2 = r.even()->GetRMS(2);
    xc = AnalysisWaveform::correlation(&y, &r, 0, rms1 * rms2);
    max_corr = xc->even()->peakVal(&loc, 0, -1, true);
    clean.updateEven()->GetY()[loc] += xc->even()->GetY()[loc] * loop_gain;
    //shift template to align with peak of correlation
    r.updateEven()->shift(xc->Neven()/2 - loc, false);
    //subtract impulse response from dirty waveform
    for(int i = 0; i < y.even()->GetN(); i++)
    {
      y.updateEven()->GetY()[i] -= loop_gain * r.even()->GetY()[i] * xc->even()->GetY()[loc];
    }
    r.updateEven()->shift(-1*(xc->Neven()/2 - loc), false);
    E_curr = y.even()->getSumV2();
    delete xc;
  }
  max_clean = clean.even()->peakVal(0, 0, -1, true);
  //rescale back up residuals and remove clean components below the noise floor or rescale back up ones that are above
  for(int i = 0; i < y.Neven(); i++)
  {
    y.updateEven()->GetY()[i] *= 1./normalize;
    clean.updateEven()->GetY()[i] *= (fabs(clean.updateEven()->GetY()[i]) > max_clean/noise_level) ? 1./normalize : 0;
  }
  if(fOnlyReturnResiduals)
  {
    for(int i = 0; i < Y->Nfreq(); i++)
    {
      Y->updateFreq()[i] = y.freq()[i];
    }
  }
  else if(fRestoringBeam.Length() > 0)
  {
    //convolve clean components w/ restoring beam
    //xc = AnalysisWaveform::convolution(&rbeam, &clean);
    //xc->even()->peakVal(&loc, 0, -1, true);
    //xc->updateEven()->shift(xc->Neven()/2 - loc, false);
    xc = new AnalysisWaveform(clean.Neven(), clean.freq(), clean.deltaF(), clean.even()->GetX()[0]);
    for(int i = 0; i < xc->Nfreq(); i++)
    {
      xc->updateFreq()[i] *= rbeam.freq()[i];
    }
    double scale_clean = fabs(clean.even()->peakVal(0, 0, -1, 1)/xc->even()->peakVal(0,0,-1,1));
    for(int i = 0; i < xc->Neven(); i++)
    {
      xc->updateEven()->GetY()[i] *= scale_clean;
    }
    //add back residuals
    if(fAddResiduals)
    {
      for(int i = 0; i < xc->Neven(); i++)
      {
        xc->updateEven()->GetY()[i] += y.even()->GetY()[i];
      }
    }
    //return the cleaned stuff
    for(int i = 0; i < Y->Nfreq(); i++)
    {
      Y->updateFreq()[i] = xc->freq()[i];
    }
    delete xc;
  }
  else
  {
    if(fAddResiduals)
    {
      for(int i = 0; i < clean.Neven(); i++)
      {
        clean.updateEven()->GetY()[i] += y.even()->GetY()[i];
      }
    }
    //return the cleaned stuff
    for(int i = 0; i < Y->Nfreq(); i++)
    {
      Y->updateFreq()[i] = clean.freq()[i];
    }
  }
}

void AnitaResponse::NaiveDeconvolution::deconvolve(AnalysisWaveform * wf, const AnalysisWaveform * response_wf) const 
{

  size_t N = wf->Nfreq(); 
  FFTWComplex * Y = wf->updateFreq(); 
  const FFTWComplex * response = response_wf->freq();

  for (unsigned i = 0; i < N; i++) 
  {
    Y[i]/=response[i]; 
  }
}

void AnitaResponse::BandLimitedDeconvolution::deconvolve(AnalysisWaveform * wf, const AnalysisWaveform * response_wf) const 
{

  size_t N = wf->Nfreq(); 
  FFTWComplex * Y = wf->updateFreq(); 
  double df = wf->deltaF(); 
  const FFTWComplex * response = response_wf->freq();

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


FFTWComplex * AnitaResponse::AbstractResponse::getResponseArray(int N, const double * f, double angle, FFTWComplex *answer) const
{
  if (!answer)  answer = new FFTWComplex[N]; 
  for (int i = 0; i < N; i++) answer[i] = getResponse(f[i], angle); 
  return answer; 

}

FFTWComplex * AnitaResponse::AbstractResponse::getResponseArray(int N, double  df, double angle, FFTWComplex *answer) const
{
  if (!answer) answer = new FFTWComplex[N]; 
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
static __thread int cache_Nt= 0; 
static __thread double cache_df= 0; 
static __thread double cache_angle= 0; 
static __thread FFTWComplex *cache_V = 0; 
static __thread AnalysisWaveform * cache_response_wf = 0;


void AnitaResponse::AbstractResponse::deconvolveInPlace(AnalysisWaveform * wf,  const DeconvolutionMethod * method, double off_axis_angle) const
{
//  printf("method: %p\n", method); 
  wf->padEven(1,0); 
  int nt = wf->Neven(); 
  double df = wf->deltaF(); 
  if (!cache_response || cache_Nt != nt || cache_df != df || !cache_V || cache_angle != off_axis_angle) 
  {
    cache_response = this; 
    cache_Nt = nt; 
    cache_df = df; 
    cache_angle = off_axis_angle; 
    if (cache_V) delete [] cache_V; 
    int nf = nt/2+1; 
    cache_V = getResponseArray(nf,df,off_axis_angle);
    if (cache_response_wf) delete cache_response_wf; 
    cache_response_wf = new AnalysisWaveform(nt, cache_V, df,0); 
  }
    
  method->deconvolve(wf, cache_response_wf); 
//  wf->updateEven()->Set(old_size); 

}



/// CLEAN 
//

void AnitaResponse::CLEANDeconvolution::clearIntermediate() const
{

  for (unsigned i = 0; i < save_xcorr.size(); i++) delete save_xcorr[i]; 
  for (unsigned i = 0; i < save_ys.size(); i++) delete save_ys[i]; 
  for (unsigned i = 0; i < save_subtracted.size(); i++) delete save_subtracted[i]; 
  save_xcorr.clear(); 
  save_ys.clear(); 
  save_subtracted.clear(); 
}

void AnitaResponse::CLEANDeconvolution::deconvolve(AnalysisWaveform *y, const AnalysisWaveform * h) const
{


  int Nt = y->Neven(); 
  double dt = y->deltaT(); 
  double t0 = y->even()->GetX()[0]; 
  double T = dt *Nt; 
  int Nf =h->Nfreq(); 
  double df = 1./dt/Nt; 

  //now let's evaluate the restoring beam
  //TODO: cache the last one since usually we'll be doing the same one over and over again 

  if (!cached_restoring || Nt != cached_restoring->Neven() || dt != cached_restoring->deltaT())
  {
    if (cached_restoring) delete cached_restoring;
    cached_restoring = new AnalysisWaveform(Nt,dt,-T/2); 
    TGraphAligned * grestore = cached_restoring->updateEven(); 
    for (int i = 0; i < Nt; i++) 
    {
      grestore->GetY()[i] = r->Eval(i*dt); 
    }
  }


  bool done = false;

  double hrms = h->getRMS(); 
  double yrms = y->getRMS(); 
  double power_ratio = yrms/hrms; 

  cmp = AnalysisWaveform(Nt,dt,t0); 

  int iter = 0; 

  clearIntermediate(); 
  double sqrt_threshold = sqrt(threshold);

  while(!done) 
  {
     AnalysisWaveform * xcorr = AnalysisWaveform::correlation(y,h,0,hrms*yrms);  

     int where; 

     double maxcorr = xcorr->even()->peakVal(&where,0,-1,true); 
     if (xcorr->even()->GetY()[where] < 0) maxcorr *=-1; 

     double lag = xcorr->even()->GetX()[where]-t0;



     //record the component here
     cmp.updateEven()->GetY()[(int) (lag/dt)]+= maxcorr * gain;


     if (debug) 
     {
       save_ys.push_back(new AnalysisWaveform(*y)); 
       save_ys.back()->setTitle(Form("iter%d pk2kpk = %g", iter, save_ys.back()->even()->pk2pk())); 
     }


     //now subtract off impulse response with this offset 
     //We are currently in the frequency domain so let's stay there! 
     const std::complex<double> * H = (const std::complex<double>*) h->freq(); 
     std::complex<double> * Y = (std::complex<double>*) y->updateFreq(); 
//     printf ("\nSTART %d\n",iter); 
//
    AnalysisWaveform * subtracted = debug ? new AnalysisWaveform(*h) : 0; 
    if (debug) save_subtracted.push_back(subtracted); 
    FFTWComplex * subtracted_freq = debug ?  subtracted->updateFreq() : 0; 

     const double factor= (maxcorr*gain*power_ratio); 
     const double expo = -2./Nt*TMath::Pi()*lag/dt; 
     for (int i = 0; i < Nf; i++) 
     {
       //this can be sped up using multiple angle formulas! 
       std::complex<double> subtract = factor * (H[i] * std::exp(std::complex<double>(0,i*expo)));
 //      printf("%g:  %g %g ->", i *df, std::real(Y[i]), std::imag(Y[i]));
       Y[i] -= subtract;
       if (debug) subtracted_freq[i] = subtract; 
  //     printf("%g %g\n", std::real(Y[i]), std::imag(Y[i]));
     }


     if (debug) 
     {
       xcorr->setTitle(Form("lag=%g, max=%g",lag,maxcorr)); 
       save_xcorr.push_back(xcorr) ; 
     }
     else 
     {
       delete xcorr; 
     }

     if (y->getRMS()  / yrms < sqrt_threshold || iter++ > max_iter) break; 
  }


  //now get rid of clean components below threshold

  double cc_max = fabs(cmp.even()->peakVal(0,0,-1,true)); 
  cmp.updateEven()->setBelow(cc_max/noiselevel); 


  if (convolve_residuals) 
  {
    //add clean components to y before convolving 


    FFTWComplex * yf = y->updateFreq(); 
    const FFTWComplex *cf = cmp.freq(); 
    for (int i = 0;  i < Nf; i++) 
    {

      yf[i] += cf[i]; 
    }
 

    AnalysisWaveform * convolved = AnalysisWaveform::convolution(y,cached_restoring, 0, cached_restoring->getRMS() *y->getRMS()); 
    *y = *convolved; 
    delete convolved; 

  }
  else 
  {
    AnalysisWaveform * convolved = AnalysisWaveform::convolution(&cmp,cached_restoring, 0, cached_restoring->getRMS() *cmp.getRMS()); 
    
    //now add to residuals 
    //
    TGraph * g = y->updateEven(); 
    for (int i = 0;  i < Nt; i++) 
    {
      g->GetY()[i] += convolved->evalEven(g->GetX()[i]); 
    }
 

    delete convolved; 
  }
  
}



