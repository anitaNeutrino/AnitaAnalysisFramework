#include "AnalysisWaveform.h" 
#include "FFTtools.h" 
#include "TAxis.h"
#include "RFInterpolate.h" 
#include <stdlib.h>
#include <assert.h>


#define ALIGNMENT 32 

#ifdef ENABLE_VECTORIZE
#include "vectorclass.h" 
#define VEC_T double
#define VEC_N 4
#define VEC Vec4d 
#define IVEC Vec4i 
static VEC INDEX(0,1,2,3); 
#endif

AnalysisWaveform::InterpolationType AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::AKIMA; 
AnalysisWaveform::InterpolationOptions AnalysisWaveform::defaultInterpolationOptions; 


/* Branch prediction help macros. Not sure if these work with clang */ 
#define likely(x)       __builtin_expect((!!x),1)
#define unlikely(x)     __builtin_expect((!!x),0)

const int minPointsForGslInterp = 6;

#ifdef ANITA_ANALYSIS_DEBUG
static bool DEBUGON = false; 
void AnalysisWaveform::enableDebug(bool debug) { DEBUGON = debug; } 
#else
void AnalysisWaveform::enableDebug(bool debug) { fprintf(stderr, "enableDebug(%d) ineffective without ANITA_ANALYSIS_DEBUG defined\n", debug); } 
#endif 


static bool ALLOW_EVEN_TO_UNEVEN = false; 
void AnalysisWaveform::allowEvenToUnevenConversion( bool allow) { ALLOW_EVEN_TO_UNEVEN = allow; }


static __thread bool nag_if_not_zero_padded = true; 


/** give things unique names */ 
volatile static unsigned counter = 0; 


static void fillEven(int N, double * x, double dt, double t0)
{
#ifdef ENABLE_VECTORIZE

  int leftover = N % VEC_N; 
  int nit = N / VEC_N + (leftover ? 1 : 0); 
  VEC vt0(t0); 
  VEC vdt(dt); 

  for (int i = 0; i < nit; i++)
  {
    VEC index(i * VEC_N); 
    index += INDEX; 
    VEC vx = mul_add(index,vdt,vt0); 

    if (i == nit-1 && leftover > 0) 
    {
      vx.store_partial(leftover, &x[i*VEC_N]); 
    }
    else
    {
      vx.store(&x[i * VEC_N]); 
    }
  }
#else
  for (int i = 0; i < N; i++)
  {
    x[i] = i * dt + t0; 
  }
#endif
}


static void setNewSize(TGraphAligned * g, int size )
{
  int old_size = g->GetN(); 
  g->Set(size); 
  double dt = g->GetX()[1] - g->GetX()[0]; 
  if (size > old_size)
  {
    fillEven(size-old_size, g->GetX() + old_size, dt, g->GetX()[old_size-1] + dt); 
  }
}


static int complaints = 0; 

AnalysisWaveform::AnalysisWaveform(int N, const double *x, const double * y, double nominal_dt,
                                   InterpolationType interp_type, InterpolationOptions *opt)
  : g_uneven(N,x,y),  dt(nominal_dt), fft(0), theEvenAkimaInterpolator(0,ROOT::Math::Interpolation::kAKIMA),
    interpolation_type(interp_type), must_update_uneven(false), must_update_freq(true), must_update_even(true), 
    uneven_equals_even(false), hilbert_transform(0), force_even_size(0)
{

  if (opt) interpolation_options = *opt; 
  power_dirty = true; 
  hilbert_dirty = true; 
  even_akima_interpolator_dirty = true; 
  power_db_dirty = true; 
  hilbert_envelope_dirty = true; 
  phase_dirty = true;
  group_delay_dirty = true;
  just_padded = false; 
  fft_len = 0; 

  uid = __sync_fetch_and_add(&counter,1); 
  nameGraphs(); 

  //HACK HACK HACK 
  //check for zero's at the end
  for (int i = 1; i < g_uneven.GetN(); i++)
  {
    if (g_uneven.GetX()[i] < g_uneven.GetX()[i-1])
    {
      if (complaints++ < 100)
      {
        fprintf(stderr,"WARNING: truncating graph to %d points to make monotonic\n",i); 
        if (complaints == 100)
        {
          fprintf(stderr," shutting up about truncating graphs after 100 complaints... I think I've made my point clear!\n"); 
        }
      }
      g_uneven.Set(i); 
      break; 
    }
  }

  //should prevent interpolation error from too few points
  while (g_uneven.GetN() < minPointsForGslInterp) {  
      g_uneven.SetPoint(g_uneven.GetN(), g_uneven.GetX()[g_uneven.GetN()-1] + 1, 0); 
  }
 

}

AnalysisWaveform::AnalysisWaveform(int N, const double * y, double dt, double t0)
  : g_even(N), dt(dt), fft(0), theEvenAkimaInterpolator(0,ROOT::Math::Interpolation::kAKIMA), 
     interpolation_type(AKIMA), must_update_uneven(false), must_update_freq(true), must_update_even(false),
     uneven_equals_even(true), hilbert_transform(0), force_even_size(0)
{
  fillEven(N,g_even.GetX(),dt,t0); 
  memcpy(g_even.GetY(), y, sizeof(double) * N); 
  fft_len = N/2 +1; 
  df = 1./(N * dt); 
  power_dirty = true; 
  hilbert_dirty = true; 
  power_db_dirty = true; 
  hilbert_envelope_dirty = true; 
  even_akima_interpolator_dirty = true; 
  phase_dirty = true;
  group_delay_dirty = true;
  just_padded = false; 
  uid = __sync_fetch_and_add(&counter,1); 
  nameGraphs(); 
  df = 0; 
}


AnalysisWaveform::AnalysisWaveform(int N, double dt, double t0)
  : g_even(N), dt(dt), fft(0), theEvenAkimaInterpolator(0,ROOT::Math::Interpolation::kAKIMA), 
    interpolation_type(AKIMA), must_update_uneven(false), must_update_freq(true),
    must_update_even(false), uneven_equals_even(true), hilbert_transform(0), force_even_size(0)
{
  
  fillEven(N,g_even.GetX(),dt,t0); 
  memset(g_even.GetY(),0, sizeof(double) * N); 
  fft_len = N/2 +1; 
  df = 1./(N * dt); 
  power_dirty = true; 
  hilbert_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true;
  group_delay_dirty = true;
  hilbert_envelope_dirty = true; 
  even_akima_interpolator_dirty = true; 
  uid = __sync_fetch_and_add(&counter,1); 
  nameGraphs(); 
  just_padded = false; 

}

AnalysisWaveform::AnalysisWaveform(int N, const FFTWComplex * f, double df, double t0)
  :  g_even(N), df(df), theEvenAkimaInterpolator(0,ROOT::Math::Interpolation::kAKIMA), interpolation_type(AKIMA), 
     must_update_uneven(false), must_update_freq(false), must_update_even(true), uneven_equals_even(true), 
     hilbert_transform(0), force_even_size(0)
{
  fft_len = N/2 +1; 
  int ret = posix_memalign( (void**) &fft, ALIGNMENT, sizeof(FFTWComplex) * fft_len);
  assert (!ret); 

  memcpy(fft, f, fft_len * sizeof(FFTWComplex)); 
  dt = 1./ (N*df); 

  fillEven(N,g_even.GetX(),dt,t0); 
  memset(g_even.GetY(),0, sizeof(double) * N); 

  even_akima_interpolator_dirty = true; 
  power_dirty = true; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
  group_delay_dirty = true; 
  uid = __sync_fetch_and_add(&counter,1); 
  just_padded = false; 
  even_akima_interpolator_dirty = true; 
}



const TGraphAligned * AnalysisWaveform::uneven() const
{

  if (uneven_equals_even) return even(); 

  return &g_uneven; 
}

const TGraphAligned * AnalysisWaveform::even()  const
{

#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called even()!\n"); 
  if (DEBUGON) printf ("[%p]: \tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n", this, must_update_even, must_update_freq, must_update_uneven); 
  //kaboom 
#endif
  assert (!(must_update_even && must_update_freq &&  must_update_uneven)); 

  if ((must_update_even && must_update_freq)) 
  {
    calculateEvenFromUneven(); 
  }
  else if ((must_update_even))
  {
    calculateEvenFromFreq(); 
  }

  return &g_even; 
}

const FFTWComplex * AnalysisWaveform::freq()  const
{


#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called freq()!\n"); 
  if (DEBUGON) printf ("[%p] \tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n", this, must_update_even, must_update_freq, must_update_uneven); 
#endif

  if (must_update_freq) 
  {
    calculateFreqFromEven(); 
  }

  return fft; 
}




void AnalysisWaveform::updateEven(const TGraph * replace)
{

#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called updateEven(replace)!\n"); 
  if (DEBUGON) printf ("[%p] \tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n", this, must_update_even, must_update_freq, must_update_uneven); 
#endif

  g_even.Set(replace->GetN()); 
  memcpy(g_even.GetX(), replace->GetX(), replace->GetN() * sizeof(double)); 
  memcpy(g_even.GetY(), replace->GetY(), replace->GetN() * sizeof(double)); 
  if (!ALLOW_EVEN_TO_UNEVEN) uneven_equals_even = true; 
  must_update_uneven = !uneven_equals_even; 
  must_update_freq = true; 
  must_update_even = false; 
  even_akima_interpolator_dirty = true; 
  just_padded = false; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 
}


TGraphAligned * AnalysisWaveform::updateEven()
{

#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called updateEven()!\n"); 
  if (DEBUGON) printf ("\tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n", must_update_even, must_update_freq, must_update_uneven); 
#endif 

  TGraphAligned * ev = (TGraphAligned *) even(); 
  if (!ALLOW_EVEN_TO_UNEVEN) uneven_equals_even = true; 
  must_update_uneven = !uneven_equals_even; 
  must_update_freq = true; 
  even_akima_interpolator_dirty = true; 
  just_padded = false; 
  hilbert_envelope_dirty = true; 
  hilbert_dirty = true; 
  return ev; 
}

TGraphAligned * AnalysisWaveform::updateUneven()
{
#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called updateUneven()!\n"); 
  if (DEBUGON) printf ("[%p]: \tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n",this, must_update_even, must_update_freq, must_update_uneven); 
#endif
  TGraphAligned * g = (TGraphAligned*) uneven(); 
  uneven_equals_even = false; 
  must_update_even = true; 
  must_update_freq = true; 
  just_padded = false; 
  even_akima_interpolator_dirty = true; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 

  return g; 
}

FFTWComplex * AnalysisWaveform::updateFreq() 
{
#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called updateFreq()!\n"); 
  if (DEBUGON) printf ("[%p]: \tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n",this, must_update_even, must_update_freq, must_update_uneven); 
#endif
  FFTWComplex * fr = (FFTWComplex*) freq();

  power_dirty = true; 
  power_db_dirty = true; 
  hilbert_dirty = true;
  hilbert_envelope_dirty = true;
  phase_dirty = true; 
  group_delay_dirty = true; 
  just_padded = false; 

  even_akima_interpolator_dirty = true; 
  must_update_even = true; 
  if (!ALLOW_EVEN_TO_UNEVEN) uneven_equals_even = true; 
  must_update_uneven = !uneven_equals_even; 

  return fr; 
}


#define DO_FOR_ALL(THING)   g_uneven.THING; g_hilbert_envelope.THING; g_even.THING; g_power.THING; g_power_db.THING; g_phase.THING; g_group_delay.THING; 
#define DO_FOR_TIME(THING)  g_uneven.THING; g_hilbert_envelope.THING; g_even.THING;
#define DO_FOR_FREQ(THING)  g_power.THING; g_power_db.THING; g_phase.THING; g_group_delay.THING; 

void AnalysisWaveform::setColor(int c) 
{
  DO_FOR_ALL(SetLineColor(c)); 
  DO_FOR_ALL(SetMarkerColor(c)); 
}

void AnalysisWaveform::setTitle(const char * title) 
{
  DO_FOR_ALL(SetTitle(title)); 
}

void AnalysisWaveform::setFreqDisplayRange(double l, double h)
{
  DO_FOR_FREQ(GetXaxis()->SetRangeUser(l,h)); 
}

void AnalysisWaveform::setTimeDisplayRange(double l, double h)
{
  DO_FOR_TIME(GetXaxis()->SetRangeUser(l,h)); 
}

void AnalysisWaveform::setWidth(int w) 
{
  DO_FOR_ALL(SetLineWidth(w)); 

}




void AnalysisWaveform::updateUneven(const TGraph * replace)
{

#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called updateUneven(replace)!\n"); 
  if (DEBUGON) printf ("[%p]: \tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n",this, must_update_even, must_update_freq, must_update_uneven); 
#endif


  g_uneven.Set(replace->GetN()); 
  memcpy(g_uneven.GetX(), replace->GetX(), replace->GetN() * sizeof(double)); 
  memcpy(g_uneven.GetY(), replace->GetY(), replace->GetN() * sizeof(double)); 

  uneven_equals_even = false; 
  just_padded = false; 
  must_update_even = true; 
  must_update_uneven = false; 
  must_update_freq = true; 
  even_akima_interpolator_dirty = true; 
}

void AnalysisWaveform::calculateUnevenFromEven() const
{
  const TGraphAligned * g = even(); 


  if (interpolation_type == AKIMA) 
  {

   
    const ROOT::Math::Interpolator * irp = evenAkimaInterpolator(); 
    for (int i = 0; i < g_uneven.GetN(); i++) 
    {
      if (g_uneven.GetX()[i] < g->GetX()[g->GetN()-1]) 
        g_uneven.GetY()[i] = irp->Eval(g_uneven.GetX()[i]); 
    }
  }
  else 
  {

    FFTtools::getInterpolatedGraphSparseInvert(g, &g_uneven, interpolation_options.max_distance, 0, 0, 
                                                  interpolation_type == SPARSE_YEN ? 0 : interpolation_options.mu, 
                                                  interpolation_options.regularization_order); 


  }

  must_update_uneven = false; 
}

const ROOT::Math::Interpolator * AnalysisWaveform::evenAkimaInterpolator() const
{

  if (even_akima_interpolator_dirty) 
  {
    const TGraphAligned * g = even();
    if(g->GetN() >= minPointsForGslInterp){
      theEvenAkimaInterpolator.SetData(g->GetN(),g->GetX(),g->GetY());
      even_akima_interpolator_dirty= false;
    }
    else{
      // perhaps this needs better handling, but for now...
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ": Insufficient points for Akima interpolation." << std::endl;
    }
  }
  return &theEvenAkimaInterpolator; 
}

void AnalysisWaveform::calculateEvenFromUneven()  const
{


  const TGraphAligned * g = uneven(); 
  int npoints = (g->GetX()[g->GetN()-1] - g->GetX()[0]) / dt; 
  g_even.Set(force_even_size ? force_even_size : npoints); 

  // this has been applied so we can forget about it now 
  force_even_size = 0; 

  double t0 = g->GetX()[0]; 
  for (int i = 0; i < g_even.GetN(); i++) 
  {
    double x = t0 + i * dt; 
    g_even.GetX()[i] = x;   

  }
   
  if (interpolation_type == AKIMA) 
  {
    ROOT::Math::Interpolator irp(g->GetN(), ROOT::Math::Interpolation::kAKIMA);
    if(g->GetN() >= minPointsForGslInterp){
      irp.SetData(g->GetN(),g->GetX(),g->GetY());
    }
    else{
      // perhaps this needs better handling, but for now...
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ": Insufficient points for Akima interpolation." << std::endl;
    }

    // irp.SetData(g->GetN(), g->GetX(), g->GetY());

    for (int i = 0; i < g_even.GetN(); i++) 
    {
      if (g_even.GetX()[i] > g->GetX()[g->GetN()-1])
      {
        memset(g_even.GetY()+i, 0, (g_even.GetN() - i )*sizeof(double)); 
        break ; 
   
//        g_even.GetY()[i] = 0; 
      }

      else
      {
        g_even.GetY()[i] = irp.Eval(g_even.GetX()[i]); 
      }

//      printf("%d %f\n", i, g_even.GetY()[i]); 
    }
  }
  else 
  {
    FFTtools::getInterpolatedGraphSparseInvert(g, &g_even, interpolation_options.max_distance, 0, 0, 
                                                  interpolation_type == SPARSE_YEN ? 0 : interpolation_options.mu, 
                                                  interpolation_options.regularization_order); 

  }
  hilbert_envelope_dirty = true; 

  must_update_even = false; 
}

void AnalysisWaveform::calculateEvenFromFreq() const
{
  FFTtools::doInvFFTNoClobber(g_even.GetN(), fft, g_even.GetY()); 

  double t0 = g_even.GetX()[0]; 

  if (force_even_size) 
  {
    setNewSize(&g_even, force_even_size); 
    force_even_size = 0; //applied so we can forget about it now. 
  }
  if (g_even.GetX()[1] - t0 != dt) 
  {
    fillEven(g_even.GetN(), g_even.GetX(), dt,t0); 
  }
  
  hilbert_envelope_dirty = true; 

  must_update_even = false; 
}

void AnalysisWaveform::calculateFreqFromEven() const
{

  const TGraphAligned * g = even();  // in case it must be converted from uneven

  assert(g->GetN() > 0); 
  dt = g->GetX()[1] - g->GetX()[0]; 
  int old_fft_len = fft_len; 
  fft_len = g->GetN()/2+1;  
  df = 1./ (g->GetN() * dt); 


  if (fft &&  old_fft_len != fft_len) 
  {
    free(fft); 
    fft = 0; 
  }

  if (!fft) 
  {

    int ret = posix_memalign( (void**) &fft, ALIGNMENT, sizeof(FFTWComplex) * fft_len);
    assert(!ret); 
  }

  FFTtools::doFFT(g->GetN(), g->GetY(),fft); 

  must_update_freq = false; 
  power_dirty = true; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
  group_delay_dirty = true; 
}



void AnalysisWaveform::updateFreq(int new_N, const FFTWComplex * new_fft, double new_df )
{

#ifdef ANITA_ANALYSIS_DEBUG
  if (DEBUGON) printf("Called updateFreq(replace)!\n"); 
  if (DEBUGON) printf ("[%p]: \tmust_update_even=%d, must_update_freq=%d, must_update_uneven=%d\n",this, must_update_even, must_update_freq, must_update_uneven); 
#endif

  if (new_N && new_N != g_even.GetN())
  {
    free(fft); 
    fft = 0; 
    fft_len = new_N/2 + 1; 
    int ret = posix_memalign( (void**) &fft, ALIGNMENT, sizeof(FFTWComplex) * fft_len); 
    assert(!ret); 
    g_even.Set(new_N); 
  }

  memcpy(fft, new_fft, fft_len * sizeof(FFTWComplex)); 
  if (new_df)
  {
    df = new_df; 
  }
  dt = 1./ (g_even.GetN() * df); 

  must_update_freq = false; 
  must_update_even = true; 
  if (!ALLOW_EVEN_TO_UNEVEN) uneven_equals_even = true; 
  must_update_uneven = !uneven_equals_even; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
  group_delay_dirty = true; 
  just_padded = false; 
}

const TGraphAligned * AnalysisWaveform::phase() const
{
  const FFTWComplex * the_fft = freq(); //will update if necessary
  if (phase_dirty)
  {

    g_phase.Set(fft_len); 
    for (int i = 0; i < fft_len; i++) 
    {
        g_phase.SetPoint(i, i * df, the_fft[i].getPhase()); 
    }

    phase_dirty = false;
  }

  return &g_phase; 
}

const TGraphAligned * AnalysisWaveform::groupDelay() const
{
  if (group_delay_dirty)
  {

    g_group_delay.adopt(phase()); 
    g_group_delay.GetYaxis()->SetTitle("Group Delay (ns) "); 
    g_group_delay.GetXaxis()->SetTitle("f (Ghz)"); 

    FFTtools::unwrap( g_group_delay.GetN(), g_group_delay.GetY(), 2 * TMath::Pi()); 

    double last = g_group_delay.GetY()[0]; 
    double dw = df * 2 * TMath::Pi(); 
    for (int i = 1; i < fft_len; i++) 
    {
      double current = g_group_delay.GetY()[i] ; 
      g_group_delay.GetY()[i-1] = (last - current) /dw; 
      last = current; 
    }

    g_group_delay.GetY()[g_group_delay.GetN()-1] = 0; 
    group_delay_dirty = false; 
  }

  return &g_group_delay; 
}

const TGraphAligned * AnalysisWaveform::powerdB() const
{
  const TGraphAligned * the_power = power(); //will update if necessary
  if (power_db_dirty)
  {

    g_power_db.adopt(the_power);
    g_power_db.dBize(); 
    g_power_db.GetXaxis()->SetTitle("Frequency"); 
    g_power_db.GetYaxis()->SetTitle("Power (dBm)"); 
    power_db_dirty = false;
  }

  return &g_power_db; 
}


const AnalysisWaveform * AnalysisWaveform::hilbertTransform() const 
{

  (void) freq(); 

  if (hilbert_dirty)
  {
    if (hilbert_transform) delete hilbert_transform; 
    hilbert_transform = new AnalysisWaveform(*this); 

    FFTWComplex * hilbert = hilbert_transform->updateFreq(); 
    for (int i = 0; i < fft_len; i++) 
    {//todo: check if this gets vectorized  
      double temp_im = hilbert[i].im; 
      hilbert[i].im = hilbert[i].re; 
      hilbert[i].re = - temp_im; 
    }

    hilbert_dirty = false; 
  }


  return hilbert_transform; 
}

const TGraphAligned * AnalysisWaveform::hilbertEnvelope() const
{
  const TGraphAligned * the_even = even(); 

  if (hilbert_envelope_dirty)
  {
    g_hilbert_envelope.adopt(the_even); 
    const AnalysisWaveform * hilbert = hilbertTransform(); 

    double * y = g_hilbert_envelope.GetY(); 
    const double * H = hilbert->even()->GetY(); 
    for (int i = 0; i < g_hilbert_envelope.GetN(); i++) 
    {
      y[i] = sqrt(y[i] * y[i] + H[i]*H[i]); 
    }

    hilbert_envelope_dirty = false; 

  }

  return &g_hilbert_envelope; 
}

const TGraphAligned * AnalysisWaveform::power() const
{

  if (power_dirty)
  {

    g_power.GetXaxis()->SetTitle("f (GHz)"); 
    g_power.GetYaxis()->SetTitle("Power (arb)"); 
    if (power_options.method == PowerCalculationOptions::FFT) 
    {
      const FFTWComplex * the_fft = freq(); //will update if necessary
      g_power.Set(fft_len); 
      for (int i = 0; i < fft_len; i++) 
      {
        if (i == 0 || i == fft_len -1)
        {
          g_power.SetPoint(i, i * df, the_fft[i].getAbsSq()/fft_len/50/1000); 
        }
        else
        {
          g_power.SetPoint(i, i * df, the_fft[i].getAbsSq()*2/fft_len/50/1000); 
        }
      }
    }
    else 
    {
      const TGraphAligned * g = even(); 

      if (power_options.method == PowerCalculationOptions::BARTLETT)
      {
        FFTtools::welchPeriodogram(g, power_options.window_size, 0, &FFTtools::RECTANGULAR_WINDOW, true,&g_power); 
      }
      else 
      {
        FFTtools::welchPeriodogram(g, power_options.window_size, 0.5, power_options.type,true,&g_power); 
      }
    }

    power_dirty = false;
  }

  return &g_power; 
}

AnalysisWaveform::~AnalysisWaveform() 
{
  if (fft) 
  {
    free(fft); 
  }

  if (hilbert_transform) 
    delete hilbert_transform; 

}

void AnalysisWaveform::forceEvenSize(int size)
{
  if (must_update_even)  //defer 
  {
    force_even_size = size; 
    return; 
  }
  setNewSize(&g_even, size); 

}

void AnalysisWaveform::evalEven(int N, const double * __restrict t, double * __restrict v, EvenEvaluationType type)  const
{

  if (type == EVAL_LINEAR) 
  {
    const TGraphAligned * g = even(); 
    int Ng = g->GetN(); 
    if (Ng < 1) return; 
    const double *y = g->GetY(); 
  //  __builtin_prefetch(y); 
  //  __builtin_prefetch(t); 
  //  __builtin_prefetch(v,1); 

    double t0 = g->GetX()[0]; 



#ifndef ENABLE_VECTORIZE
    for (int i = 0; i < N; i++) 
    {
      double xval = (t[i]-t0) / dt; 
      int bin_low = int(xval); 

      int branchless_too_low = bin_low < 0; 
      int branchless_too_high = bin_low >= Ng; 
      int branchless_on_edge = bin_low == Ng-1; 

      int bin_high = bin_low+1; 
      bin_low *= (1-branchless_too_high) * (1-branchless_too_low);
      bin_high *= (1-branchless_too_high) * (1-branchless_too_low) * (1-branchless_on_edge);
      double val_low = y[bin_low]; 
      double val_high = y[bin_high]; 
      double frac = xval - bin_low; 
      double val = val_low + frac * (val_high - val_low); 

      v[i] = val * (1-branchless_too_low) * (1-branchless_too_high) * (1-branchless_on_edge) + branchless_on_edge * y[Ng-1]; 
    }
#else

    int leftover = N % VEC_N; 
    int nit = N / VEC_N + (leftover ? 1 : 0);  

    VEC v_t; 
    VEC v_t0(t0); 
    VEC inv_dt(1./dt); 
    IVEC dont_get_out_of_bounds; 

    for (int i = 0; i < VEC_N; i++)
    {
      dont_get_out_of_bounds.insert(i, i >= leftover ? 1 : 0 ); 
    }

    for (int i = 0; i < nit; i++)
    {
      if (i < nit -1 || !leftover)
      {
          v_t.load(t + i * VEC_N); 
      }
      else
      {
        v_t.load_partial(leftover, t+i*VEC_N)  ;
      }

      VEC xval = (v_t - v_t0) * inv_dt; 
      VEC truncated_xval = truncate(xval); 
      VEC frac = xval - truncated_xval; 
      IVEC bin_low = truncate_to_int(truncated_xval); 

      IVEC branchless_too_low = bin_low < IVEC(0); 
      IVEC branchless_too_high = bin_low >= IVEC(Ng); 
      IVEC branchless_on_edge = bin_low== IVEC(Ng-1); 

      int scalar_out_of_bounds =  int(i == nit-1 && leftover > 0); 
      IVEC out_of_bounds = scalar_out_of_bounds * dont_get_out_of_bounds;

      IVEC bin_high = bin_low+1; 

      //Optimistically grab it. Hopefully it's not out of bounds. 
      __builtin_prefetch(y+bin_low[0]); 

      IVEC too_low_or_too_high = (1-branchless_too_low) * (1-branchless_too_high) * (1-out_of_bounds); 
      IVEC kill_it = too_low_or_too_high *  ( 1- branchless_on_edge);  

      bin_low  *= too_low_or_too_high; 
      bin_high *= kill_it; 

      double  lowv[VEC_N]; 
      double  highv[VEC_N]; 
      for (int j = 0; j < VEC_N; j++)
      {
        lowv[j] = y[bin_low[j]]; 
        highv[j] = y[bin_high[j]]; 
      }

      VEC val_low; val_low.load(lowv); 
      VEC val_high; val_high.load(highv); 

      VEC val = mul_add(frac , (val_high - val_low), val_low); 
      
      if (scalar_out_of_bounds)
      {
        val.store_partial(leftover, v + i * VEC_N); 
      }
      else
      {
        val.store(v + i * VEC_N); 
      }
    }

#endif

  }

  else if (type == EVAL_AKIMA) 
  {
    const ROOT::Math::Interpolator * irp = evenAkimaInterpolator(); 

    /*TODO this is probably very slow */ 
    for (int i = 0; i < N; i++) 
    {
      v[i] = irp->Eval(t[i]); 
    }
  }
  else
  {
    fprintf(stderr,"Bad even evaluation type\n"); 
  }

}

double AnalysisWaveform::evalEven(double t, EvenEvaluationType typ)  const
{
  if (typ == EVAL_AKIMA) 
  {
    return evenAkimaInterpolator()->Eval(t); 
  }
  else if (typ == EVAL_LINEAR) 
  {
    const TGraphAligned * g = even(); 
    double t0 = g->GetX()[0]; 
    double xval = (t-t0) / dt; 

    int bin_low = int (xval); 

    if (bin_low < 0) return 0; 
    if (bin_low >= g->GetN()) return 0; 
    if (bin_low ==  g->GetN()-1) return g->GetY()[g->GetN()-1]; 

    int bin_high = bin_low + 1; 
    double val_low = g->GetY()[bin_low]; 
    double val_high = g->GetY()[bin_high]; 

    double frac = xval - bin_low; 

    return val_low + frac*(val_high-val_low); 
  }

  fprintf(stderr, "Bad EvenEvaluationType\n"); 

  return 0; 
}

AnalysisWaveform & AnalysisWaveform::operator=(const AnalysisWaveform & other) 
{
  if (this != &other) 
  {
    g_uneven.adopt(&other.g_uneven); 
    g_even.adopt(&other.g_even); 
    
    dt = other.dt; 
    df = other.df; 
    interpolation_type = other.interpolation_type; 
    interpolation_options = other.interpolation_options; 

    uneven_equals_even = other.uneven_equals_even; 
    must_update_even = other.must_update_even; 
    must_update_freq = other.must_update_freq; 
    must_update_uneven = other.must_update_uneven; 
    just_padded = other.just_padded; 
    force_even_size = other.force_even_size; 
  //  printf("%d\n", force_even_size); 


    //don't bother copying these, they can be regenerated if needed
    power_dirty = true; 
    power_db_dirty = true; 
    phase_dirty = true; 
    group_delay_dirty = true; 
    hilbert_dirty = true; 
    hilbert_envelope_dirty = true; 
    hilbert_transform = 0; 

    
    fft_len = g_even.GetN()/2+1;

    if (!must_update_freq)
    {
      int ret = posix_memalign( (void**) &fft, ALIGNMENT, sizeof(FFTWComplex) * fft_len);
      assert(!ret);
      memcpy(fft, other.fft, fft_len * sizeof(FFTWComplex)); 
    }
    else
    {
      fft = 0; 
    }
  }

  return *this; 
}

AnalysisWaveform::AnalysisWaveform(const AnalysisWaveform & other) 
  :
  g_uneven(other.uneven_equals_even ? 0 : other.g_uneven.GetN()), 
  g_even(other.Neven()) 
{
  

  //first copy over scalars 
  
  dt = other.dt; 
  df = other.df; 

  interpolation_type = other.interpolation_type; 
  interpolation_options = other.interpolation_options; 

  uneven_equals_even = other.uneven_equals_even; 
  must_update_even = other.must_update_even; 
  must_update_freq = other.must_update_freq; 
  must_update_uneven = other.must_update_uneven; 
  just_padded = other.just_padded; 
  force_even_size = other.force_even_size; 
//  printf("%d\n", force_even_size); 


  //don't bother copying these, they can be regenerated if needed
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
  group_delay_dirty = true; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 
  hilbert_transform = 0; 



  // we must copy x from g_uneven if it's not equal to even 
  if (!uneven_equals_even)
  {
    memcpy(g_uneven.GetX(), other.g_uneven.GetX(), g_uneven.GetN() * sizeof(double)); 

    if (!must_update_uneven)
    {
      memcpy(g_uneven.GetY(), other.g_uneven.GetY(), g_uneven.GetN() * sizeof(double)); 
    }
  }
 

  if (!must_update_even)
  {
    memcpy(g_even.GetY(), other.g_even.GetY(), g_even.GetN() * sizeof(double)); 
    memcpy(g_even.GetX(), other.g_even.GetX(), g_even.GetN() * sizeof(double)); 
  }
  
  fft_len = g_even.GetN()/2+1;

  if (!must_update_freq)
  {
    int ret = posix_memalign( (void**) &fft, ALIGNMENT, sizeof(FFTWComplex) * fft_len);
    assert(!ret);
    memcpy(fft, other.fft, fft_len * sizeof(FFTWComplex)); 
  }
  else
  {
    fft = 0; 
  }

}

AnalysisWaveform * AnalysisWaveform::correlation(const AnalysisWaveform *A, const AnalysisWaveform *B, int npad, double scale) 
{
  if (A->Nfreq() != B->Nfreq()) 
  {
    fprintf(stderr,"correlation does not handle the case where A and B are of different lengths (%d vs. %d)!\n", A->Nfreq(), B->Nfreq()); 
    return 0; 
  }


  if(nag_if_not_zero_padded && (!A->checkIfPaddedInTime() || !B->checkIfPaddedInTime()))
  {

    fprintf(stderr,"warning: waveforms don't appear to be padded in time, will be computing circular correlation!\n"); 
  }
  double offset = A->even()->GetX()[0] - B->even()->GetX()[0]; 

  AnalysisWaveform * answer = new AnalysisWaveform(A->Neven(), A->freq(), 
                                                   A->deltaF(), A->even()->GetX()[0]); 

  int N = answer->even()->GetN(); 
  FFTWComplex * update = answer->updateFreq(); 

  const FFTWComplex * Bfreq = B->freq(); 
  double inv = 1./(N*scale); 
  
  //TODO this can be vectorized 
  for (int i = 0; i < B->Nfreq(); i++) 
  {
    FFTWComplex vA = update[i]; 
    FFTWComplex vB = Bfreq[i]; 
    update[i].re =  (vA.re * vB.re + vA.im * vB.im) *inv;
    update[i].im =  (vA.im * vB.re - vA.re * vB.im) *inv; 
//    printf("%f %f\n", update[i].re, update[i].im); 
  }

  answer->padFreq(npad); 

  N = answer->even()->GetN(); 
  TGraphAligned g(N); 

  double dt = answer->deltaT(); 
  memcpy(g.GetY(), answer->even()->GetY() + N/2, N/2 * sizeof(double)); 
  memcpy(g.GetY()+ N/2, answer->even()->GetY(), N/2 * sizeof(double)); 

  for (int i = 0; i < N; i++) 
  {
    g.GetX()[i] =(i - N/2) * dt + offset; 
  }


  answer->updateEven(&g); 


  return answer; 
}

void AnalysisWaveform::padFreq(int npad)
{
  if (npad < 1) return; 



  //new even size
  if (g_even.GetN() ==0) (void) even();  //may need to calculate even if it never has been
  int new_N = g_even.GetN() * (1+npad); 

  FFTWComplex * new_fft =0; 
  int ret = posix_memalign( (void**) &new_fft, ALIGNMENT, sizeof(FFTWComplex) * (new_N/2+1));
  assert(!ret); 

  const FFTWComplex * old_freq = freq(); 
  //printf("%d\n", fft_len); 

  //copy old
  memcpy(new_fft, old_freq, TMath::Min(new_N/2+1,fft_len) * sizeof(FFTWComplex));

  //scale
  double scale = (1 + npad); 
  for (int i =0; i < Nfreq(); i++) { new_fft[i].re *= scale;  new_fft[i].im *= scale; }


  // zero rest 
  memset(new_fft + fft_len, 0, (new_N/2 + 1 - fft_len) * sizeof(FFTWComplex));

  //set new lengths and dt
  fft_len = new_N/2 + 1; 
  g_even.Set(new_N); 
  dt = 1./ (g_even.GetN() * df); 


  //swap in new fft
  free(fft); 
  fft = new_fft; 

  //others need to update now! 
  must_update_even = true; 
  must_update_uneven = !uneven_equals_even; 
  must_update_freq = false; 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
  group_delay_dirty = true; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 
  just_padded = false; 
}



void AnalysisWaveform::padFreqAdd(int npad)
{
  if (npad < 1) return; 



  //new even size
  if (g_even.GetN() ==0) (void) even();  //may need to calculate even if it never has been
  int old_N = g_even.GetN();
  int new_N = g_even.GetN() + npad; 

  FFTWComplex * new_fft =0; 
  int ret = posix_memalign( (void**) &new_fft, ALIGNMENT, sizeof(FFTWComplex) * (new_N/2+1));
  assert(!ret); 

  const FFTWComplex * old_freq = freq(); 
  //printf("%d\n", fft_len); 

  //copy old
  memcpy(new_fft, old_freq, TMath::Min(new_N/2+1,fft_len) * sizeof(FFTWComplex));

  //scale
  double scale = (double(new_N) / old_N); 
  for (int i =0; i < Nfreq(); i++) { 
    new_fft[i].re *= scale;
    new_fft[i].im *= scale; }


  // zero rest 
  memset(new_fft + fft_len, 0, (new_N/2 + 1 - fft_len) * sizeof(FFTWComplex));

  //set new lengths and dt
  fft_len = new_N/2 + 1; 
  g_even.Set(new_N); 
  dt = 1./ (g_even.GetN() * df); 


  //swap in new fft
  free(fft); 
  fft = new_fft; 

  //others need to update now! 
  must_update_even = true; 
  must_update_uneven = !uneven_equals_even; 
  must_update_freq = false; 
  power_dirty = true; 
  power_db_dirty = true; 
  phase_dirty = true; 
  group_delay_dirty = true; 
  hilbert_dirty = true; 
  hilbert_envelope_dirty = true; 
  just_padded = false; 

}


void AnalysisWaveform::setPowerCalculationOptions(PowerCalculationOptions & opt) 
{
  power_options = opt; 
  power_dirty = true; 
  power_db_dirty = true; 
}

void AnalysisWaveform::padEven(int npad, int where)
{
  if (npad < 1) return; 
  TGraphAligned * g = updateEven(); 
  int old_n = g->GetN(); 
  g->Set(g->GetN() *(1+npad)); 
  df /= (1+npad); 



  if (where > 0) 
  {
    for (int i = old_n; i < g->GetN(); i++) 
    {
      g->GetX()[i] = g->GetX()[i-1] + dt; 
      g->GetY()[i] = 0; 
    }
  }
  else  if (where < 0) 
  {
    memcpy(g->GetY() + old_n, g->GetY(), old_n * sizeof(double)); 
    memcpy(g->GetX() + old_n, g->GetX(), old_n * sizeof(double)); 
    for (int i = old_n-1; i >= 0; i--)
    {
      g->GetX()[i] = g->GetX()[i+1] - dt; 
      g->GetY()[i] = 0; 
    }
  }

  else if (where == 0) 
  {
    memmove(g->GetY() + old_n/2, g->GetY(), old_n * sizeof(double)); 
    memmove(g->GetX() + old_n/2, g->GetX(), old_n * sizeof(double)); 

    for (int i = old_n/2-1; i >= 0; i--) 
    {
      g->GetX()[i] = g->GetX()[i+1] - dt; 
      g->GetY()[i] = 0; 
    }
    for (int i = old_n + old_n/2; i < g->GetN(); i++) 
    {
      g->GetX()[i] = g->GetX()[i-1] + dt; 
      g->GetY()[i] = 0; 
    }
  }


  if (npad >= 1) just_padded = true; 
}

bool AnalysisWaveform::checkIfPaddedInTime() const 
{
  if (just_padded) return true; 

  const TGraphAligned *g = even(); 

  for (int i = g->GetN()/2; i <g->GetN(); i++) 
  {
    if (g->GetY()[i] !=0) return false; 
  }
  just_padded = true; 
  return true; 
}

static void doDraw(const TGraphAligned * cg,const char * opt,  int color)
{
  TGraphAligned * g = (TGraphAligned*) cg; 

  if (color >= 0) 
  {
    g->SetLineColor(color); 
    g->SetMarkerColor(color); 
  }
  g->Draw(opt); 
}


void AnalysisWaveform::drawEven(const char * opt, int color) const{ doDraw(even(),opt,color); }
void AnalysisWaveform::drawUneven(const char * opt, int color) const{ doDraw(uneven(),opt,color); }
void AnalysisWaveform::drawPower(const char * opt, int color)const { doDraw(power(),opt,color); }
void AnalysisWaveform::drawPowerdB(const char * opt, int color)const { doDraw(powerdB(),opt,color); }
void AnalysisWaveform::drawPhase(const char * opt, int color) const{ doDraw(phase(),opt,color); }
void AnalysisWaveform::drawHilbertEnvelope(const char * opt, int color) const{ doDraw(hilbertEnvelope(),opt,color); }


int AnalysisWaveform::Nfreq() const
{
  (void) freq(); 
  return fft_len;  
}


void AnalysisWaveform::basisChange(AnalysisWaveform * __restrict x, 
    AnalysisWaveform * __restrict y)  {



  int N = TMath::Min(x->Neven(), y->Neven()); 

  //is this right? who knows
  const TGraphAligned * g_yh = y->hilbertTransform()->even(); 
  const TGraphAligned * g_y = y->even();
  const TGraphAligned * g_xh = x->hilbertTransform()->even();
  const TGraphAligned * g_x = x->even();

  const double one_over_sqrt2 = 1./sqrt(2); 

  double new_x[N] __attribute__((aligned)); 
  double new_y[N] __attribute__((aligned)); 

  for (int i = 0; i < N; i++) 
  {
    new_x[i] = one_over_sqrt2 * ( g_xh->GetY()[i] + g_y->GetY()[i]); 
    new_y[i] = one_over_sqrt2 * ( g_x->GetY()[i] + g_yh->GetY()[i]); 
  }


  memcpy(x->updateEven()->GetY(), new_x, N * sizeof(double)); 
  memcpy(y->updateEven()->GetY(), new_y, N * sizeof(double)); 

  x->forceEvenSize(N); 
  y->forceEvenSize(N); 

}

void AnalysisWaveform::sumDifference( AnalysisWaveform * __restrict x, AnalysisWaveform * __restrict y) 
{
  int N = TMath::Min(x->Neven(), y->Neven()); 

  //is this right? who knows
  const TGraphAligned * g_x = x->even(); 
  const TGraphAligned * g_y = y->even();

  double new_x[N] __attribute__((aligned)); 
  double new_y[N] __attribute__((aligned)); 

  for (int i = 0; i < N; i++) 
  {
    new_x[i] = 0.5 * ( g_x->GetY()[i] + g_y->GetY()[i] ); 
    new_y[i] = 0.5 * ( g_x->GetY()[i] - g_y->GetY()[i] ); 
  }


  memcpy(x->updateEven()->GetY(), new_x, N * sizeof(double)); 
  memcpy(y->updateEven()->GetY(), new_y, N * sizeof(double)); 

  x->forceEvenSize(N); 
  y->forceEvenSize(N); 
}


AnalysisWaveform * AnalysisWaveform::autoCorrelation(int npadtime, int npadfreq, double scale) 
{

  if (npadtime)
  {
    AnalysisWaveform copy(*this); 
    copy.padEven(1); 
    return correlation(&copy,&copy,npadfreq,scale); 
  }

  return correlation(this,this,npadfreq,scale); 
}

void AnalysisWaveform::zeroMeanEven()
{
  TGraphAligned * g = updateEven(); 

  double mean = g->GetMean(2); 
  for (int i = 0; i <g->GetN(); i++) 
  {
    g->GetY()[i]-=mean; 
  }
}

void AnalysisWaveform::setCorrelationNag(bool nag) 
{

  nag_if_not_zero_padded = nag; 
}

void AnalysisWaveform::nameGraphs() 
{
  g_uneven.SetName(TString::Format("wf_uneven_%d", uid)); 
  g_even.SetName(TString::Format("wf_even_%d", uid)); 
  g_hilbert_envelope.SetName(TString::Format("wf_hilbert_env_%d", uid)); 
  g_power.SetName(TString::Format("wf_power_%d", uid)); 
  g_power_db.SetName(TString::Format("wf_power_db_%d", uid)); 
  g_phase.SetName(TString::Format("wf_phase_%d", uid)); 
  g_group_delay.SetName(TString::Format("wf_group_delay_%d", uid)); 
}



