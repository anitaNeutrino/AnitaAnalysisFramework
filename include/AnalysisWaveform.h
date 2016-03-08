#ifndef ANALYSIS_WAVEFORM_HH
#define ANALYSIS_WAVEFORM_HH


/* This class is similar in principle to RFWaveform from FFTtools. 
 *
 * It holds 3 coupled versions of the waveform: 
 *
 *   uneven: An unevenly sampled time-domain waveform 
 *   even: An evenly sampled time-domain waveform
 *   freq An evenly sampled frequency-domain waveform 
 *
 *   Transforming between uneven and even is accomplished via interpolation (which is not strictly invertible) 
 *   Transforming between even and freq is accomplished via FT (which is more or less invertible, up to machine precision)
 *
 *   Because the transformation between even and uneven is not reversible, it usually only makes sense to access
 *   or modify uneven prior to modifying even or freq. If an access to uneven is attempted after modifying even or freq, 
 *
 *
 *   There are two types of accessors, const ones which retrieve const data: 
 *     uneven() 
 *     even() 
 *     freq() 
 *
 *
 *  And ones which allow modification of the data and will then recompute the other two versions the next time they're requested: 
 *
 *   updateUneven(); 
 *   updateEven(); 
 *   updateFreq(); 
 *
 *  Those have two versions, one where you modify the return value and one where you replace the value. 
 *   
 *
 *  In order to make it seem magic, there is a lot of internal state that is complicated to reason about. Hopefully there are no bugs.. .
 *
 *   
 *
 * 
 *
 */  
class FFTWComplex; 

#include "TGraph.h" 

class AnalysisWaveform 
{

  public: 
    enum InterpolationType
    {
      AKIMA, 
      SPARSE_YEN,    
      REGULARIZED_SPARSE_YEN    
    }; 

    /** Interpolation options, relevant only for certain interpolators */ 
    struct InterpolationOptions 
    {

      InterpolationOptions() 
        : max_distance(8), regularization_order(0), mu(1e-3) {; } 

      int max_distance;  
      int regularization_order; 
      double mu; 
    }; 

    /** Constructors */ 
    AnalysisWaveform(int Nt, const double * x, const double * y, double nominal_dt = 1./2.6, InterpolationType type = AKIMA, InterpolationOptions * opt  = 0);  // constructor from unevenly sampled waveform
    AnalysisWaveform(int Nt, const double * y, double dt, double t0);  //constructor from evenly sampled waveform (in which case, even and uneven are initialized to the same) 
    AnalysisWaveform(int Nt, const FFTWComplex * f, double df, double t0);  //constructor from frequency domain waveform (in which case, even and uneven are the initialized to the same) 
    AnalysisWaveform(const AnalysisWaveform & other);  // copy constructor... more subtle than you might think! 

    /* Computes the (circular) correlation (in the frequency domain) of the two waveforms as a new waveform. Note that if you want to
     * correlate two traces, they should be padded first. This does not pad them for you. It is also assumed the two are of the same length.
     *
     * There is no normalization done at all, the frequency values are simply multiplied appropriately
     *
     **/ 
    static AnalysisWaveform * correlation(const AnalysisWaveform * A, const AnalysisWaveform * B); 

    virtual ~AnalysisWaveform(); 

    /** Constant accessors. If you coerce the compiler into allowing modification, coupled waveforms won't be updated */ 
    const TGraph * uneven() const ;
    const TGraph * even()  const; 
    const FFTWComplex * freq() const; 
    const TGraph * power() const; 
    const TGraph * powerdB() const; 
    const TGraph * phase() const; 
    int Nfreq() const { (void) freq(); return fft_len; } 
    int Neven() const { return even()->GetN(); } 
    int Nuneven() const { return uneven()->GetN(); } 
    double deltaT() const { return dt ; } 
    double deltaF() const { return df ; }

    // drawers since drawing is non-const (and we don't care about silly things like axes for constness)
    void drawEven(const char * opt = "") const; 
    void drawUneven(const char * opt = "") const; 
    void drawPower(const char * opt = "") const; 
    void drawPowerdB(const char * opt = "") const; 
    void drawPhase(const char * opt = "") const; 


    /*  Update frequency graph by modifying return value */ 
    FFTWComplex * updateFreq(); 
    /*  Update frequency graph by modifying return value */ 
    void updateFreq(int new_N, const FFTWComplex * new_freq, double new_df = 0) ;
    TGraph * updateEven(); 
    void updateEven(const TGraph * replace_even); 
    TGraph * updateUneven(); 
    void updateUneven(const TGraph * replace_uneven); 


    // pad the even waveform (equivalent to upsampling the frequency) 
    void padEven(int factor); 

    // pad the frequency (equivalent to upsampling the even values)
    void padFreq(int factor); 



  private: 
    void calculateEvenFromUneven() const; 
    void calculateEvenFromFreq() const;  
    void calculateFreqFromEven() const;  
    void calculateUnevenFromEven() const;

    mutable TGraph g_uneven; 
    mutable TGraph g_even; 
    mutable TGraph g_power; 
    mutable TGraph g_power_db; 
    mutable TGraph g_phase; 
    mutable double dt; 
    mutable double df; 
    mutable int fft_len; 
    mutable FFTWComplex * fft; 

    InterpolationType interpolation_type; 
    InterpolationOptions interpolation_options; 
    mutable bool must_update_uneven; 
    mutable bool must_update_freq; 
    mutable bool must_update_even; 
    mutable bool uneven_equals_even; 
    mutable bool power_dirty; 
    mutable bool power_db_dirty; 
    mutable bool phase_dirty; 
}; 


#endif 
